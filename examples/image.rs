use anyhow::{Result, ensure, bail};
use std::fs::File;
use rayon::prelude::*;

use palletize::*;

fn main() -> Result<()> {
    // Arg parsing
    let mut args = std::env::args().skip(1);
    let image_path = args.next().expect("Requires image path arg first");

    let palette: Vec<[u8; 3]> = args.map(|s| parse_hex_color(&s)).collect();
    let palette_cie: Vec<[f32; 3]> = palette.iter().copied().map(px_to_cielab).collect();

    assert!(palette.len() > 0);

    // PNG decode
    let decoder = png::Decoder::new(File::open(image_path)?);
    let (info, mut reader) = decoder.read_info()?;

    ensure!(info.bit_depth == png::BitDepth::Eight, "Only eight-bit images are supported");
    let n_components = match info.color_type {
        png::ColorType::RGB => 3,
        png::ColorType::RGBA => 4,
        ty => bail!("Unsupported color type {:?}", ty),
    };

    let mut buf = vec![0; info.buffer_size()];
    reader.next_frame(&mut buf)?;

    // Palettize
    //buf.par_chunks_exact_mut(n_components).for_each(|px| {
    buf.chunks_exact_mut(n_components).for_each(|px| {
        let rgb = [px[0], px[1], px[2]];
        let lab = px_to_cielab(rgb);
        let mut best_dist = std::f32::MAX;
        let mut best_match = palette[0];
        for (idx, pal) in palette_cie.iter().enumerate() {
            let dist = ciede2000_diff(lab, *pal);
            if dist < best_dist {
                best_match = palette[idx];
                best_dist = dist;
            }
        }
        px[0] = best_match[0];
        px[1] = best_match[1];
        px[2] = best_match[2];
    });

    // Encode and write
    let mut encoder = png::Encoder::new(File::create("out.png")?, info.width, info.height);
    encoder.set_color(info.color_type);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header().unwrap();
    writer.write_image_data(&buf)?;

    Ok(())
}