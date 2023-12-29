//! Main executable for the Rust Vortex Lattice model.

mod vlm;

use vlm::Vlm;

const N: usize = 250;

const SPAN: f64 = 10.0;

const ROOT_CHORD: f64 = 1.0;
const TIP_CHORD: f64 = 0.0;

const ROOT_AOA: f64 = 0.0;
const TIP_AOA: f64 = -2.0;

const TIP_DIHEDRAL: f64 = 0.0;
const ROOT_DIHEDRAL: f64 = 0.0;

fn main() {
    let mut alpha = [0.0; N];
    let mut a = TIP_AOA;
    for i in 0..(N/2) {
        alpha[i] = a;
        alpha[N - i - 1] = a;
        a += (ROOT_AOA - TIP_AOA) / ((N/2) as f64);
    }

    let mut beta = [0.0; N];
    let mut b = TIP_DIHEDRAL;
    for i in 0..(N/2) {
        beta[i] = -b;
        beta[N - i - 1] = b;
        b += (ROOT_DIHEDRAL - TIP_DIHEDRAL) / ((N/2) as f64);
    }

    let mut chord = [0.0; N];
    let mut c = TIP_CHORD;
    for i in 0..(N/2) {
        chord[i] = c;
        chord[N - i - 1] = c;
        c += (ROOT_CHORD - TIP_CHORD) / ((N/2) as f64);
    }

    let vlm = Vlm::<N>::new(
        SPAN,
        chord,
        alpha,
        beta,
    );

    let aoa = 4.0;

    println!(
        "AoA = {:.2} | CL = {:.6} | CDi = {:.6} | L/D = {:.6} | AR = {:.6} | e = {:.6}",
        aoa,
        vlm.cl(aoa),
        vlm.cd(aoa),
        vlm.ld(aoa),
        vlm.aspect_ratio(),
        vlm.spaneff(aoa),
    );
}