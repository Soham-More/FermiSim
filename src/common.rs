pub type VecD = na::DVector<f64>;

pub mod constants{
    pub const K:f64 = 1.380649e-23;
    pub const Q:f64 = 1.60217663e-19;
    pub const EPSILON_VACCUM:f64 = 8.8542e-12;
    pub const PLANK_CONST:f64 = 6.62607015e-34;
    pub const ELECTRON_MASS:f64 = 9.1093837e-31;

    pub fn thermal_pot(temp:f64) -> f64 { (temp * K) / Q }
    pub fn from_eV(energy:f64) -> f64 { energy * Q }
}

pub mod interp {
    #[derive(Debug)]
    pub enum Types {
        Nearest,
        Linear
    }

    pub use Types::*;

    pub fn nearest1D(x:f64, fi:&Vec<f64>, xi:&Vec<f64>) -> f64
    {
        let (id, _) = xi.iter().enumerate()
        .min_by(|a ,b| (a.1 - x).abs().total_cmp(&f64::abs(b.1 - x)))
        .expect("invalid input: xi is empty!");

        return fi[id];
    }

    // assumes xi is sorted
    // O(n), TODO: make it faster
    // TODO: cover out of range x values 
    pub fn linear1D(x:f64, fi:&Vec<f64>, xi:&Vec<f64>) -> f64
    {
        for (i, &val) in xi.iter().enumerate()
        {
            if val > x
            {
                let slope = (fi[i] - fi[i - 1]) / (xi[i] - xi[i - 1]);
                return fi[i - 1] + slope * (x - xi[i - 1]);
            }
        }

        return *fi.last().expect("invalid input: xi is empty!");
    }
}

pub mod stats {
    use super::constants;

    pub fn fermi_dirac(E:f64, fermi_lvl:f64, temp:f64, degeneracy:f64) -> f64
    {
        let thermal_pot = constants::K * temp;
        let normed_energy = (E - fermi_lvl) / thermal_pot;

        1.0 / ( 1.0 + degeneracy * f64::exp(normed_energy) )
    }
    pub fn fermi_dirac_derivativeE(E:f64, fermi_lvl:f64, temp:f64, degeneracy:f64) -> f64
    {
        let thermal_pot = constants::K * temp;
        let normed_energy = (E - fermi_lvl) / thermal_pot;

        -degeneracy * f64::exp(normed_energy) / ( thermal_pot * (1.0 + degeneracy * f64::exp(normed_energy)).powi(2) )
    }
    pub fn fermi_dirac_derivativeF(E:f64, fermi_lvl:f64, temp:f64, degeneracy:f64) -> f64
    {
        let thermal_pot = constants::K * temp;
        let normed_energy = (E - fermi_lvl) / thermal_pot;

        degeneracy * f64::exp(normed_energy) / ( thermal_pot * (1.0 + degeneracy * f64::exp(normed_energy)).powi(2) )
    }

}

pub mod sim {

    struct SimulationProfile {
        pub temp: f64,
        pub charge_error_tol: f64,
        pub potential_error_reltol: f64,
        pub ss_max_iter: usize,
    }

}

