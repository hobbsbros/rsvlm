//! Main executable for the Rust Vortex Lattice model.

mod vlm;

use vlm::Vlm;

const N: usize = 250;

fn main() {
    let mut alpha = [0.0; N];
    let mut a = -2.0;
    for i in 0..(N/2) {
        alpha[i] = a;
        alpha[N - i - 1] = a;
        a += 2.0 / ((N/2) as f64);
    }

    let vlm = Vlm::<N>::new(
        10.0,
        1.0,
        alpha,
        [0.0; N],
    );

    let mut a = 0.0;

    while a < 10.0 {
        println!(
            "a = {:.2} | CL = {:.6} | CDi = {:.6} | L/D = {:.6} | e = {:.6}",
            a,
            vlm.cl(a),
            vlm.cd(a),
            vlm.ld(a),
            vlm.spaneff(a),
        );
        a += 1.0;
    }
}