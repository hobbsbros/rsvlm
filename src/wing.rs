//! Input file parser.

use std::path::Path;

pub struct Wing<const N: usize> {
    /// Name
    name: String,

    /// Reference span (meters)
    b: f64,

    /// Chord distribution (Y: meters, chord: meters)
    c: Vec<(f64, f64)>,

    /// Angle of attack distribution (Y: meters, angle: degrees)
    alpha: Vec<(f64, f64)>,

    /// Dihedral distribution (Y: meters, angle: degrees)
    beta: Vec<(f64, f64)>,
}

impl<const N: usize> Wing<N> {
    /// Constructs a new wing.
    pub fn from(path: &Path) -> Self {
        todo!()
    }

    /// Constructs a chord array using linear interpolation.
    pub fn chord(&self) -> [f64; N] {
        let mut chord = [0.0; N];

        chord
    }

    /// Constructs an angle of attack array using linear interpolation.
    pub fn alpha(&self) -> [f64; N] {
        let mut alpha = [0.0; N];

        alpha
    }

    /// Constructs a dihedral angle array using linear interpolation.
    pub fn beta(&self) -> [f64; N] {
        let mut beta = [0.0; N];

        beta
    }
}