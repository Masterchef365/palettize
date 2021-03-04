use std::f32::consts::PI;

fn srgb_to_ciexyz(rgb: [f32; 3]) -> [f32; 3] {
    let mut cie = [0.; 3];
    // https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
    const MATRIX: [[f32; 3]; 3] = [
        [0.41239080, 0.35758434, 0.18048079],
        [0.21263901, 0.71516868, 0.07219232],
        [0.01933082, 0.11919478, 0.95053215],
    ];

    for (cie, (rgb, row)) in cie.iter_mut().zip(rgb.iter().zip(MATRIX.iter())) {
        for m in row {
            *cie += *rgb * m;
        }
    }
    cie
}

fn ciexyz_to_cielab([x, y, z]: [f32; 3]) -> [f32; 3] {
    // https://en.wikipedia.org/wiki/CIELAB_color_space#From_CIEXYZ_to_CIELAB[11]
    let [xn, yn, zn]: [f32; 3] = [95.0489, 100., 108.8840];
    let delta: f32 = 6. / 29.;
    let f = |t: f32| {
        if t > delta.powf(3.) {
            t.cbrt()
        } else {
            (t / (3. * delta)) + (4. / 29.)
        }
    };

    let f_y_yn = f(y / yn);
    let l = 116. * f_y_yn - 16.;
    let a = 500. * (f(x / xn) - f_y_yn);
    let b = 200. * (f_y_yn - f(z / zn));
    [l, a, b]
}

fn to_float([r, g, b]: [u8; 3]) -> [f32; 3] {
    let f = |c: u8| c as f32 / u8::MAX as f32;
    [f(r), f(g), f(b)]
}

fn px_to_cielab(rgb: [u8; 3]) -> [f32; 3] {
    ciexyz_to_cielab(srgb_to_ciexyz(to_float(rgb)))
}

fn parse_hex_color(s: &str) -> [u8; 3] {
    let s = s.trim_start_matches('#');
    let v = u32::from_str_radix(s, 16).expect(s);
    [
        ((v & 0xFF0000) >> (2 * 8)) as u8,
        ((v & 0x00FF00) >> (1 * 8)) as u8,
        (v & 0x0000FF) as u8,
    ]
}

fn ciede2000_diff([l1, a1, b1]: [f32; 3], [l2, a2, b2]: [f32; 3]) -> f32 {
    let euclid_dist = |a: f32, b: f32| (a * a + b * b).sqrt();
    let c_ab = (euclid_dist(a1, b1) + euclid_dist(a2, b2)) / 2.;
    let c_ab_7 = c_ab.powf(7.);
    let g = 0.5 * (1. - (c_ab_7 / (c_ab_7 + 25.0f32.powf(7.))));

    let ai = |a: f32| (1. + g) * a;
    let ci = |a: f32, b: f32| euclid_dist(ai(a), b);
    let hi = |a: f32, b: f32| {
        if a == 0. && b == 0. {
            0.
        } else {
            b.atan2(ai(a)).to_degrees() + 360.
        }
    };

    let c1 = ci(a1, b1);
    let c2 = ci(a2, b2);

    let h1 = hi(a1, b1);
    let h2 = hi(a2, b2);
    let hdiff = h2 - h1;

    let dl = l2 - l1;
    let dc = c2 - c1;
    let cc = c2 * c1;

    let dh = if cc == 0. {
        0.
    } else if hdiff.abs() <= 180. {
        hdiff
    } else if hdiff > 180. {
        hdiff - 360.
    } else {
        hdiff + 360.
    };

    let dh_2 = 2. * cc.sqrt() * (dh / 2.).sin();

    let l_bar = (l1 + l2) / 2.;
    let c_bar = (c1 + c2) / 2.;
    let h_bar = if cc == 0. {
        h1 + h2
    } else if hdiff.abs() < 180. {
        (h1 + h2) / 2.
    } else if hdiff.abs() > 180. && (h1 + h2) < 360. {
        (h1 + h2 + 360.) / 2.
    } else {
        (h1 + h2 - 360.) / 2.
    };

    let t = 1. 
        - 0.17 * (h_bar - 30.).to_radians().cos()
        + 0.24 * (2. * h_bar).to_radians().cos()
        + 0.32 * (3. * h_bar + 6.).to_radians().cos()
        - 0.20 * (4. * h_bar - 63.).to_radians().cos();

    dbg!(t);

    let delta_theta = 30. * ((h_bar - 275.) / 25.).exp();
    let rc = 2. * (c_bar.powf(7.) / (c_bar.powf(7.) + 25.0f32.powf(7.))).sqrt();
    let lbsq = (l_bar - 50.).powf(2.);
    let sl = 1. + ((0.0015 * lbsq) / (20. + lbsq).sqrt());

    let sc = 1. + 0.045 * c_bar;
    let sh = 1. + 0.015 * c_bar * t;
    let rt = -(2. * delta_theta).to_radians().sin() * rc;

    (
        (dl / sl).powf(2.)
            + (dc / sc).powf(2.)
            + (dh_2 / sh).powf(2.)
            + (rt * (dc / sc) * (dh / sh))
    ).sqrt()
}

fn main() {
    //let args = std::env::args().skip(1);
    //let palette: Vec<[f32; 3]> = args.map(|s| px_to_cielab(parse_hex_color(&s))).collect();
    let a = [50., 2.6772, -79.7751];
    let b = [50., 0., -82.7485];
    let diff = ciede2000_diff(a, b);
    dbg!(diff);
}

#[test]
fn test_parse_hex_color() {
    assert_eq!(parse_hex_color("00AF00"), [0x00, 0xAF, 0x00]);
    assert_eq!(parse_hex_color("00AFBB"), [0x00, 0xAF, 0xBB]);
    assert_eq!(parse_hex_color("EEAF00"), [0xEE, 0xAF, 0x00]);
}
