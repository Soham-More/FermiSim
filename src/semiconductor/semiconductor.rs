use core::f64;

use crate::common::*;
use crate::fdm1D::Mesh;
use super::bulk::*;
use super::doping::*;

#[derive(Debug)]
pub struct Semiconductor
{
    pub bulk:Bulk,
    pub dopants:Vec<Dopant>,
    begin_pos:f64,
    end_pos:f64,
}

impl Semiconductor {

    pub fn create(bulk:Bulk) -> Semiconductor
    {
        Semiconductor {
            bulk,
            begin_pos:-f64::INFINITY,
            end_pos:f64::INFINITY,
            dopants:Vec::new()
        }
    }

    pub fn set_bulk_range(&mut self, begin_pos:f64, end_pos:f64)
    {
        self.begin_pos = begin_pos;
        self.end_pos = end_pos;
    }

    pub fn is_inside(&self, x:f64) -> bool
    {
        x >= self.begin_pos && x <= self.end_pos
    }

    pub fn electron_conc(&self, x:f64, fermi_lvl:f64, potential:f64, temp:f64) -> f64
    {
        if !self.is_inside(x)
        {
            return 0.0
        }

        self.bulk.electron_conc(fermi_lvl, potential, temp)
    }

    pub fn hole_conc(&self, x:f64, fermi_lvl:f64, potential:f64, temp:f64) -> f64
    {
        if !self.is_inside(x)
        {
            return 0.0
        }

        self.bulk.hole_conc(fermi_lvl, potential, temp)
    }

    pub fn push_dopant(&mut self, dopant:Dopant)
    {
        self.dopants.push(dopant);
    }

    pub fn total_dopant_charge_vec(&self, mesh:&Mesh) -> VecD
    {
        mesh.makeVecFn(|x, _| if self.is_inside(x) {self.dopants.iter().map(|d| d.dopant_charge(x)).sum::<f64>()} else { 0.0 })
    }

    pub fn total_charge(&self, x:f64, fermi_lvl:f64, potential:f64, temp:f64) -> f64
    {
        if !self.is_inside(x)
        {
            return 0.0;
        }

        let mut charge = 0.0;

        for dopant in self.dopants.iter()
        {
            charge += dopant.dopant_charge(x);
        }

        charge + self.bulk.electron_charge(fermi_lvl, potential, temp) + self.bulk.hole_charge(fermi_lvl, potential, temp)
    }

    pub fn total_charge_derivative_pot(&self, x:f64, fermi_lvl:f64, potential:f64, temp:f64) -> f64
    {
        if !self.is_inside(x)
        {
            return 0.0;
        }

        self.bulk.electron_charge_derivative_pot(fermi_lvl, potential, temp) + self.bulk.hole_charge_derivative_pot(fermi_lvl, potential, temp)
    }

    // vectorize this?
    pub fn total_charge_vec(&self, mesh:&Mesh, fermi_lvl:f64, potential:&VecD, temp:f64) -> VecD
    {
        mesh.makeVecFn(|x, i| self.total_charge(x, fermi_lvl, potential[i], temp))
    }

    pub fn total_charge_derivative_pot_vec(&self, mesh:&Mesh, fermi_lvl:f64, potential:&VecD, temp:f64) -> VecD
    {
        mesh.makeVecFn(|x, i| self.total_charge_derivative_pot(x, fermi_lvl, potential[i], temp))
    }

}



