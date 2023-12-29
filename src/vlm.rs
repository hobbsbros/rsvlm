//! Vortex lattice model.

use maria_linalg::{
    Matrix,
    Vector,
};

const PI: f64 = 3.141592653;

/// Converts an angle from degrees to radians.
fn rad(angle: f64) -> f64 {
    angle * PI / 180.0
}

/// Converts an angle from radians to degrees.
fn deg(angle: f64) -> f64 {
    angle * 180.0 / PI
}

pub struct Vlm<const N: usize> {
    /// Span (meters)
    b: f64,

    /// Mean aerodynamic chord (meters)
    c: f64,

    /// Angle of attack (degrees)
    alpha: Vector<N>,

    /// Local dihedral angle (degrees)
    beta: Vector<N>,

    /// Inverse circulation tensor (unitless)
    invq: Matrix<N>,
}

impl<const N: usize> Vlm<N> {
    /// Constructs a new vortex lattice.
    pub fn new(b: f64, c: f64, alpha: [f64; N], beta: [f64; N]) -> Self {
        let mut alpha_rad = [0.0f64; N];
        let mut beta_rad = [0.0f64; N];
        for i in 0..N {
            alpha_rad[i] = rad(alpha[i]);
            beta_rad[i] = rad(beta[i]);
        }

        let mut obj = Self {
            b,
            c,
            alpha: alpha_rad.into(),
            beta: beta_rad.into(),
            invq: Matrix::zero(),
        };
        
        let q = obj.circulation();
        
        obj.invq = q.inverse();

        obj
    }

    /// Constructs the downwash vector, given an angle
    ///     of attack in *degrees*.
    pub fn downwash(&self, aoa: f64) -> Vector<N> {
        let mut k = Matrix::zero();

        let s = self.b / (N as f64);

        for m in 0..N {
            for n in 0..N {
                let i = m as f64;
                let j = n as f64;
                k[(m, n)] = 1.0 / (4.0 * PI * s) * (-1.0 / (j - i + 0.5) + 1.0 / (j - i - 0.5));
            }
        }

        let g = self.solve(aoa);

        k.mult(g)
    }

    /// Constructs the induced angle of attack vector in *degrees*,
    ///     given an angle of attack in *degrees*.
    pub fn induced_aoa(&self, aoa: f64) -> Vector<N> {
        // Compute downwash
        let w = self.downwash(aoa);

        let mut ai = Vector::zero();

        for i in 0..N {
            ai[i] = -deg(w[i].atan2(1.0));
        }

        ai
    }

    /// Constructs the Kronecker delta.
    pub fn kronecker() -> Matrix<N> {
        let mut delta = Matrix::zero();

        for k in 0..N {
            delta[(k, k)] = 1.0;
        }

        delta
    }

    /// Constructs the circulation tensor.
    pub fn circulation(&self) -> Matrix<N> {
        let mut q = Matrix::zero();

        let xi = self.b / (N as f64) / self.c;
        let delta = Self::kronecker();

        for m in 0..N {
            for n in 0..N {
                let i = m as f64;
                let j = n as f64;
                q[(m, n)] += delta[(m, n)] / PI;
                q[(m, n)] += self.beta[m].cos() / (4.0 * PI * xi) * (1.0 / (j - i + 0.5) - 1.0 / (j - i - 0.5));
            }
        }

        q
    }

    /// Computes the distributed circulation strengths,
    ///     given an angle of attack in *degrees*.
    pub fn solve(&self, aoa: f64) -> Vector<N> {
        // Compute the angles of attack
        let mut alpha = self.alpha;
        for i in 0..N {
            alpha[i] += rad(aoa);
        }

        self.invq.mult(alpha)
    }

    /// Computes the lift distribution, given an angle
    ///     of attack in *degrees*.
    pub fn lift_distribution(&self, aoa: f64) -> Vector<N> {
        let mut distr = Vector::zero();
        let g = self.solve(aoa);
        let ai = self.induced_aoa(aoa);

        for i in 0..N {
            distr[i] = g[i] * self.beta[i].cos() * rad(ai[i]).cos() * self.c;
        }

        distr
    }

    /// Computes the coefficient of lift, given an angle
    ///     of attack in *degrees*.
    pub fn cl(&self, aoa: f64) -> f64 {
        let mut cl = 0.0;
        let cldistr = self.lift_distribution(aoa);

        for i in 0..N {
            cl += cldistr[i] * self.b / (N as f64);
        }

        cl / (0.5 * self.b * self.c)
    }
    
    /// Computes the coefficient of drag, given an angle
    ///     of attack in *degrees*.
    pub fn cd(&self, aoa: f64) -> f64 {
        let mut cd = 0.0;
        let g = self.solve(aoa);
        let ai = self.induced_aoa(aoa);

        for i in 0..N {
            cd += 2.0 * g[i] * self.beta[i].cos() * rad(ai[i]).sin() / (N as f64);
        }

        cd
    }

    /// Computes the lift-drag ratio, given an angle of
    ///     attack in *degrees*.
    pub fn ld(&self, aoa: f64) -> f64 {
        self.cl(aoa) / self.cd(aoa)
    }

    /// Computes the span efficiency of the wing,
    ///     given an angle of attack in *degrees*.
    pub fn spaneff(&self, aoa: f64) -> f64 {
        let cl = self.cl(aoa);
        let cd = self.cd(aoa);

        cl.powi(2) / (cd * PI * self.b / self.c)
    }
}