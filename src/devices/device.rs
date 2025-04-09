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
    epsilon:VecD,
    vacc_Ec:VecD,
    vacc_Ev:VecD,

    pub net_doping:VecD,

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
            epsilon:VecD::default(),
            vacc_Ec:VecD::default(),
            vacc_Ev:VecD::default(),
            net_doping:VecD::default(),
            mesh:Mesh::create(vec![0.0]),
            steady_state:State::default(),
            full_width:0.0,
        }
    }

    pub fn push_bulk_layer(&mut self, mut layer: Semiconductor, width:f64, samples:u32)
    {
        self.mesh.extend(
            (0..samples).map(|i| self.last_pos + f64::from(i + 1) * (width / f64::from(samples)))
            .collect()
        );

        if self.epsilon.is_empty()
        {
            self.epsilon = VecD::from_column_slice(&[layer.bulk.epsilon]);
            self.vacc_Ec = VecD::from_column_slice(&[layer.bulk.Ec]);
            self.vacc_Ev = VecD::from_column_slice(&[layer.bulk.Ev]);
        }

        self.epsilon.extend(
            (0..samples).map(|_| layer.bulk.epsilon)
        );
        self.vacc_Ec.extend(
            (0..samples).map(|_| layer.bulk.Ec)
        );
        self.vacc_Ev.extend(
            (0..samples).map(|_| layer.bulk.Ev)
        );
        
        layer.set_bulk_range(self.last_pos, self.last_pos + width);
        self.bulk_layers.push(layer);

        self.last_pos += width;
        self.full_width += width;
    }

    pub fn calc_steady_state(&mut self, charge_tol:f64, rel_potential_tol:f64, max_iter:usize)
    {
        // prepare the poission problem
        self.poissionProb = PoissionProblem::create(&self.mesh, &self.epsilon);

        let interface_layer_left = self.bulk_layers.first().expect("No layers initialized! ");
        let interface_layer_right = self.bulk_layers.last().expect("No layers initialized! ");

        // the potential at x = 0 is taken as the reference.
        self.steady_state.fermi_lvl = calcRootBisection(
            interface_layer_left.bulk.Ev, 
            interface_layer_left.bulk.Ec, 
            |mu| interface_layer_left.total_charge(0.0, mu, 0.0, self.temp), 
            rel_potential_tol,
            charge_tol,
            100
        ).expect("No fermi_lvl in max_iter using bisection method");

        self.steady_state.built_in_potential = calcRootBisection(
            -interface_layer_right.bulk.band_gap / constants::Q,
            interface_layer_right.bulk.band_gap / constants::Q,
            |V| interface_layer_right.total_charge(self.full_width, self.steady_state.fermi_lvl, V, self.temp), 
            rel_potential_tol,
            charge_tol,
            100
        ).expect("No fermi_lvl in max_iter using bisection method");

        let sample_last_idx = self.mesh.lastIdx();

        let mut potential = self.mesh.zeroVec();
        let mut charge = self.total_charge_vec(self.steady_state.fermi_lvl, &potential);
        
        potential[0] = 0.0;
        potential[sample_last_idx] = self.steady_state.built_in_potential;
        
        for i in 0..max_iter {
            let charge_derivative = self.total_charge_derivative_pot_vec(self.steady_state.fermi_lvl, &potential);

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
            potential -= w * &deltaV;

            // reinforce the boundary conditions
            potential[0] = 0.0;
            potential[sample_last_idx] = self.steady_state.built_in_potential;

            let prev_charge = charge;
            charge = self.total_charge_vec(self.steady_state.fermi_lvl, &potential);

            if (prev_charge - &charge).abs().max() < charge_tol && 
                deltaV.abs().max() / potential.abs().max() < rel_potential_tol &&
                i > 1
            {
                break;
            }

            if i == max_iter - 1
            {
                panic!("Steady state did not converge!");
            }
        }

        self.steady_state.potential = potential.clone();
        self.steady_state.charge = charge.clone();
        self.steady_state.Ec = &self.vacc_Ec - constants::Q * &potential;
        self.steady_state.Ev = &self.vacc_Ev - constants::Q * &potential;
        
        self.steady_state.n = self.mesh.makeVecFn( |x, i| 
                            self.bulk_layers.iter()
                            .map(|y| y.electron_conc(x, self.steady_state.fermi_lvl, self.steady_state.potential[i], self.temp))
                            .sum()
        );

        self.steady_state.p = self.mesh.makeVecFn( |x, i| 
            self.bulk_layers.iter()
            .map(|y| y.hole_conc(x, self.steady_state.fermi_lvl, self.steady_state.potential[i], self.temp))
            .sum()
        );

        self.net_doping = self.bulk_layers.iter()
            .map(|y| y.total_dopant_charge_vec(&self.mesh) / constants::Q)
            .sum();

    }

}

pub fn calcRootBisection(mut lower_estimate:f64, mut upper_estimate:f64, f:impl Fn(f64) -> f64, xtol:f64, ftol:f64, max_iter:usize) -> Option<f64>
{
    let mut lower_estimate_value = f(lower_estimate);
    let mut upper_estimate_value = f(upper_estimate);

    for _ in 0..max_iter
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


