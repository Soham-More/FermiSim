use core::f64;

use crate::common::*;
use rgsl::fermi_dirac::complete_integrals::*;

#[derive(Debug)]
pub struct CarrrierInfo
{
    pub effectiveMass:f64,
    pub mobility:f64,
}

#[derive(Debug)]
pub struct Bulk
{
    pub electron_affinity:f64,
    pub band_gap:f64,
    pub relative_permitivity:f64,
    pub hole_properties:CarrrierInfo,
    pub electron_properties:CarrrierInfo,
    
    pub Ev:f64,
    pub Ec:f64,
    pub epsilon:f64,
}

impl CarrrierInfo
{
    fn density_of_states(&self, temp:f64) -> f64
    {
        2.0 * (2.0 * f64::consts::PI * self.effectiveMass * constants::K * temp).powf(1.5) / constants::PLANK_CONST.powi(3)
    }
}

impl Bulk {
    pub fn create(electron_affinity:f64, band_gap:f64, relative_permitivity:f64, hole_properties:CarrrierInfo, electron_properties:CarrrierInfo) -> Bulk
    {
        Bulk {
            electron_affinity,
            band_gap,
            relative_permitivity,
            hole_properties,
            electron_properties,
            Ev: -electron_affinity-band_gap,
            Ec: -electron_affinity,
            epsilon: constants::EPSILON_VACCUM * relative_permitivity,
        }
    }

    pub fn create_silicon_300K() -> Bulk
    {
        let hole_properties = CarrrierInfo{
            mobility:0.045,
            effectiveMass:0.48*constants::ELECTRON_MASS,
        };
    
        let electron_properties = CarrrierInfo{
            mobility:0.1,
            effectiveMass:1.08*constants::ELECTRON_MASS,
        };

        Bulk {
            electron_affinity:1.3895213 * constants::Q,
            band_gap:1.14 * constants::Q,
            relative_permitivity:11.68,
            hole_properties,
            electron_properties,
            Ev: (-1.3895213-1.14) * constants::Q,
            Ec: -1.3895213 * constants::Q,
            epsilon: constants::EPSILON_VACCUM * 11.68,
        }
    }

    pub fn create_GaAs_300K() -> Bulk
    {
        // ref: https://www.ioffe.ru/SVA/NSM/Semicond/GaAs/basic.html
        // ref: https://www.ioffe.ru/SVA/NSM/Semicond/GaAs/electric.html
        let hole_properties = CarrrierInfo{
            mobility:0.4,
            effectiveMass:0.51*constants::ELECTRON_MASS,
        };

        let electron_properties = CarrrierInfo{
            mobility:8.5,
            effectiveMass:0.063*constants::ELECTRON_MASS,
        };

        Bulk::create(
            constants::from_eV(4.07), 
            constants::from_eV(1.42),
            12.9,
            hole_properties,
            electron_properties
        )
    }

    pub fn create_AlGaAs_300K(x:f64) -> Bulk
    {
        // ref: https://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/basic.html
        // ref: https://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/ebasic.html
        
        if x < 0.0 && x > 1.0 { panic!("Error: The mole_fraction must be a value between 0 and 1."); }

        if x < 0.45
        {
            let Eg = 1.422 + 1.2475*x;

            let hole_prop = CarrrierInfo{
                mobility:0.37 - 0.97*x + 0.74*x*x,
                effectiveMass:0.64*constants::ELECTRON_MASS,
            };
        
            let elec_prop = CarrrierInfo{
                mobility:0.8 - 2.2*x + x*x,
                effectiveMass:(0.063 + 0.083*x)*constants::ELECTRON_MASS,
            };
        
            Bulk::create(
                constants::from_eV(4.07 - 1.1*x), 
                constants::from_eV(Eg),
                12.9 - 2.84 * x,
                hole_prop,
                elec_prop
            )
        }
        else
        {
            let Eg = 1.9+0.125*x+0.143*x*x;

            let hole_prop = CarrrierInfo{
                mobility:0.37 - 0.97*x + 0.74*x*x,
                effectiveMass:(0.51 + 0.25*x)*constants::ELECTRON_MASS,
            };

            let elec_prop = CarrrierInfo{
                mobility:-0.225 + 1.16*x - 0.72*x*x,
                effectiveMass:(0.85 - 0.14*x)*constants::ELECTRON_MASS,
            };

            Bulk::create(
                constants::from_eV(3.64 - 0.14*x), 
                constants::from_eV(Eg),
                12.9 - 2.84 * x,
                hole_prop,
                elec_prop
            )
        }
    }

    pub fn electron_conc(&self, fermi_lvl:f64, potential:f64, temp:f64) -> f64
    {
        let Ec_potential = self.Ec - constants::Q * potential;
        let normalized_energy  = -(Ec_potential - fermi_lvl) / (constants::K * temp);

        //println!("Electron DOS {:e}, Occupied: {:e}", self.electron_properties.density_of_states(temp), fermi_dirac_half(normalized_energy));

        self.electron_properties.density_of_states(temp) * fermi_dirac_half(normalized_energy)
    }
    pub fn electron_conc_derivative_pot(&self, fermi_lvl:f64, potential:f64, temp:f64) -> f64
    {
        let Ec_potential = self.Ec - constants::Q * potential;
        let normalized_energy  = -(Ec_potential - fermi_lvl) / (constants::K * temp);

        self.electron_properties.density_of_states(temp) * fermi_dirac_mhalf(normalized_energy) * -constants::Q / -(constants::K * temp)
    }

    pub fn electron_charge(&self, fermi_lvl:f64, potential:f64, temp:f64) -> f64
    {
        -constants::Q * self.electron_conc(fermi_lvl, potential, temp)
    }
    pub fn electron_charge_derivative_pot(&self, fermi_lvl:f64, potential:f64, temp:f64) -> f64
    {
        -constants::Q * self.electron_conc_derivative_pot(fermi_lvl, potential, temp)
    }

    pub fn hole_conc(&self, fermi_lvl:f64, potential:f64, temp:f64) -> f64
    {
        let Ev_potential = self.Ev - constants::Q * potential;
        let normalized_energy  = -(fermi_lvl - Ev_potential) / (constants::K * temp);

        self.hole_properties.density_of_states(temp) * fermi_dirac_half(normalized_energy)
    }
    pub fn hole_conc_derivative_pot(&self, fermi_lvl:f64, potential:f64, temp:f64) -> f64
    {
        let Ev_potential = self.Ev - constants::Q * potential;
        let normalized_energy  = -(fermi_lvl - Ev_potential) / (constants::K * temp);

        -self.hole_properties.density_of_states(temp) * fermi_dirac_mhalf(normalized_energy) * -constants::Q / -(constants::K * temp)
    }

    pub fn hole_charge(&self, fermi_lvl:f64, potential:f64, temp:f64) -> f64
    {
        constants::Q * self.hole_conc(fermi_lvl, potential, temp)
    }
    pub fn hole_charge_derivative_pot(&self, fermi_lvl:f64, potential:f64, temp:f64) -> f64
    {
        constants::Q * self.hole_conc_derivative_pot(fermi_lvl, potential, temp)
    }

}


