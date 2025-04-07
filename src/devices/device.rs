use crate::common::*;
use crate::fdm1D::tridiag::MatTriDiag;
use crate::semiconductor::*;
use crate::fdm1D::*;
use super::state::*;

pub struct Device {
    bulk_layers:Vec<Semiconductor>,
    poissionProb:PoissionProblem,
    temp:f64,
    last_pos:f64,

    pub mesh:Mesh,
    pub steady_state:State,
    pub full_width:f64,
}

impl Device{

    fn total_charge_vec(&self, fermi_lvl:f64, potential:&VecD) -> VecD
    {
        self.bulk_layers.iter()
            .map(|layer| layer.total_charge_vec(&self.mesh, fermi_lvl, &potential, self.temp))
            .fold(self.mesh.zeroVec(), |acc, charge| acc + charge)
    }

    fn total_charge_derivative_pot_vec(&self, fermi_lvl:f64, potential:&VecD) -> VecD
    {
        self.bulk_layers.iter()
            .map(|layer| layer.total_charge_derivative_pot_vec(&self.mesh, fermi_lvl, &potential, self.temp))
            .fold(self.mesh.zeroVec(), |acc, charge| acc + charge)
    }

    pub fn create(temp:f64) -> Device
    {
        Device {
            bulk_layers:Vec::new(),
            poissionProb:PoissionProblem::default(),
            temp,
            last_pos:0.0,
            mesh:Mesh::new(),
            steady_state:State::default(),
            full_width:0.0,
        }
    }

    pub fn push_bulk_layer(&mut self, layer:Bulk, width:f64, samples:u32)
    {
        self.mesh.extend(
            (0..samples).map(|i| f64::from(i) * (width / f64::from(samples - 1)))
            .collect()
        );
        
        let layer = Semiconductor::create(layer, self.last_pos, self.last_pos + width);
        self.bulk_layers.push(layer);

        self.last_pos += width;
        self.full_width += width;
    }

    pub fn calc_steady_state(&mut self)
    {
        let interface_layer_left = self.bulk_layers.first().expect("No layers initialized! ");
        let interface_layer_right = self.bulk_layers.first().expect("No layers initialized! ");

        // the potential at x = 0 is taken as the reference.
        let fermi_lvl = calcRootBisection(
            interface_layer_left.bulk.Ev, 
            interface_layer_left.bulk.Ec, 
            |mu| interface_layer_left.total_charge(0.0, mu, 0.0, self.temp), 
            1e-5,
            1e2,
            100
        ).expect("No fermi_lvl in max_iter using bisection method");

        let built_in_potential = calcRootBisection(
            -interface_layer_right.bulk.band_gap / constants::Q,
            interface_layer_right.bulk.band_gap / constants::Q,
            |V| interface_layer_right.total_charge(self.full_width, fermi_lvl, V, self.temp), 
            1e-5,
            1e2,
            100
        ).expect("No fermi_lvl in max_iter using bisection method");

        let sample_last_idx = self.mesh.lastIdx();

        let mut potential = self.mesh.zeroVec();
        let mut charge = self.total_charge_vec(fermi_lvl, &potential);
        
        potential[0] = 0.0;
        potential[sample_last_idx] = built_in_potential;
        
        for _ in 0..100 {
            let charge_derivative = self.total_charge_derivative_pot_vec(fermi_lvl, &potential);

            // calculate the residual
            let mut residual = self.poissionProb.residue(&potential, &charge);
            residual[0] = 0.0;
            residual[sample_last_idx] = 0.0;

            let poission_mat = self.poissionProb.operator.clone();
            let mut jacobian:MatTriDiag = (poission_mat.0, &charge_derivative + poission_mat.1, poission_mat.2);

            jacobian.1[0] = 1.0;
            jacobian.2[0] = 0.0;

            jacobian.1[sample_last_idx] = 1.0;
            jacobian.0[sample_last_idx - 1] = 0.0;

            let deltaV = tridiag::solve(&jacobian, &mut self.mesh.zeroVec(), residual);

            let w = 1.5;
            potential -= w * deltaV;

            // reinforce the boundary conditions
            potential[0] = 0.0;
            potential[sample_last_idx] = built_in_potential;

            charge = self.total_charge_vec(fermi_lvl, &potential);
        }
    }

}

pub fn calcRootBisection(mut lower_estimate:f64, mut upper_estimate:f64, f:impl Fn(f64) -> f64, xtol:f64, ftol:f64, max_iter:usize) -> Option<f64>
{
    let mut lower_estimate_value = f(lower_estimate);
    let mut upper_estimate_value = f(upper_estimate);

    for i in 0..max_iter
    {
        let midpoint = (lower_estimate + upper_estimate) / 2.0;
        let deltaX = upper_estimate - lower_estimate;
        let deltaF = upper_estimate_value - lower_estimate_value;
        
        if f64::abs(deltaX/midpoint) < xtol && f64::abs(deltaF) < ftol {
            return Some(midpoint);
        }
        
        let midpoint_value = f(midpoint);

        if midpoint_value*lower_estimate_value > 0.0
        {
            lower_estimate = midpoint;
            lower_estimate_value = midpoint_value;
        }
        else if midpoint_value*upper_estimate_value > 0.0
        {
            upper_estimate = midpoint;
            upper_estimate_value = midpoint_value;
        }
    }

    return None
}


