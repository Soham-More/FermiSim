use crate::common::*;
use constants::Q;

#[derive(Debug)]
pub enum Types
{
    Acceptor,
    Donor
}

#[derive(Debug)]
pub struct Dopant
{
    pub sampled_conc:Vec<f64>,
    pub sampled_x:Vec<f64>,
    pub interp_mode:interp::Types,
    pub doping_type:Types,
    pub dopantE:f64,
    pub degeneracy:f64,
}

impl Dopant
{
    pub fn create(sampled_conc:Vec<f64>, sampled_x:Vec<f64>, interp_mode:interp::Types, doping_type:Types, dopantE:f64, degeneracy:f64) -> Dopant
    {
        Dopant {
            sampled_conc,
            sampled_x,
            interp_mode,
            doping_type,
            dopantE,
            degeneracy
        }
    }

    pub fn create_donor(sampled_conc:Vec<f64>, sampled_x:Vec<f64>, interp_mode:interp::Types, dopantE:f64, degeneracy:f64) -> Dopant
    {
        Dopant {
            sampled_conc,
            sampled_x,
            interp_mode,
            doping_type:Types::Donor,
            dopantE,
            degeneracy
        }
    }
    pub fn create_acceptor(sampled_conc:Vec<f64>, sampled_x:Vec<f64>, interp_mode:interp::Types, dopantE:f64, degeneracy:f64) -> Dopant
    {
        Dopant {
            sampled_conc,
            sampled_x,
            interp_mode,
            doping_type:Types::Acceptor,
            dopantE,
            degeneracy
        }
    }

    // get dopant conc at position x
    pub fn dopant_conc(&self, x:f64) -> f64
    {
        match self.interp_mode {
            interp::Types::Nearest => interp::nearest1D(x, &self.sampled_conc, &self.sampled_x),
            interp::Types::Linear => interp::linear1D(x, &self.sampled_conc, &self.sampled_x),
        }
    }
    // get dopant charge at position x
    pub fn dopant_charge(&self, x:f64) -> f64
    {
        let N = match self.interp_mode {
            interp::Types::Nearest => interp::nearest1D(x, &self.sampled_conc, &self.sampled_x),
            interp::Types::Linear => interp::linear1D(x, &self.sampled_conc, &self.sampled_x),
        };

        match self.doping_type {
            Types::Acceptor => -Q * N,
            Types::Donor => Q * N,
        }
    }

    // calculate the concentration of ionized dopants
    pub fn ionized_conc(&self, x:f64, fermi_lvl:f64, temp:f64) -> f64
    {
        stats::fermi_dirac(self.dopantE, fermi_lvl, temp, self.degeneracy) * self.dopant_conc(x)
    }
    // calculate the derivative of concentration of ionized dopants wrt fermi_lvl
    pub fn ionized_conc_derivative(&self, x:f64, fermi_lvl:f64, temp:f64) -> f64
    {
        stats::fermi_dirac_derivativeF(self.dopantE, fermi_lvl, temp, self.degeneracy) * self.dopant_conc(x)
    }

    // calculate the charge density of ionized dopants
    pub fn ionized_charge(&self, x:f64, fermi_lvl:f64, temp:f64) -> f64
    {
        stats::fermi_dirac(self.dopantE, fermi_lvl, temp, self.degeneracy) * self.dopant_charge(x)
    }
    // calculate the derivative of charge density of ionized dopants wrt fermi_lvl
    pub fn ionized_charge_derivative(&self, x:f64, fermi_lvl:f64, temp:f64) -> f64
    {
        stats::fermi_dirac_derivativeF(self.dopantE, fermi_lvl, temp, self.degeneracy) * self.dopant_charge(x)
    }
}


