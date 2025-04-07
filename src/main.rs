#![allow(non_snake_case)]

extern crate nalgebra as na;
extern crate nalgebra_lapack as nalpk;
extern crate rgsl;
pub mod common;
pub mod fdm1D;

pub mod semiconductor;
pub mod devices;
pub mod pyvi;

use common::constants::ELECTRON_MASS;
use fdm1D::{tridiag::{self, MatTriDiag}, Mesh, PoissionProblem};
use rgsl::error;
use rgsl::Value;
use semiconductor as sc;
use common::*;
use pyvi::PyVi;

fn test_poission_solver() {
    let mesh = Mesh::create((0..20).map(|i| f64::from(i) * 1.0).collect());
    let mut poissionProb = PoissionProblem::create(&mesh, mesh.makeVec(1.0));

    // charge density
    let rho =  mesh.makeVec(1.0);

    println!("rho: {}", rho);
    
    let soln = poissionProb.solve(&rho, -1.0, 1.0);

    println!("Soln: {}", soln);
    println!("Residual: {}", poissionProb.residue(&soln, &rho));
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

fn error_handling(error_str: &str, file: &str, line: u32, error_value: Value) {
    println!("RGSL [{:?}] '{}:{}': {}", error_value, file, line, error_str);
}

fn main() {
    error::set_error_handler(Some(error_handling));

    let len = 1e-5;
    let temp = 300.0;
    let doping = 1e21;

    let hole_prop = sc::CarrrierInfo{
        mobility:0.045,
        effectiveMass:0.48*ELECTRON_MASS,
    };

    let elec_prop = sc::CarrrierInfo{
        mobility:0.1,
        effectiveMass:1.08*ELECTRON_MASS,
    };

    let silicon = sc::Bulk::create_silicon_300K(hole_prop, elec_prop);

    let boron = sc::Dopant::create_acceptor(
        vec![doping, 0.0], 
        vec![0.0, len], 
        interp::Nearest, 
        silicon.Ev + 0.045*constants::Q, 
        4.0
    );
    let phosphorus = sc::Dopant::create_donor(
        vec![0.0, doping], 
        vec![0.0, len], 
        interp::Nearest, 
        silicon.Ec - 0.045*constants::Q, 
        2.0
    );

    let mut diode = sc::Semiconductor::create(silicon, 0.0, len);

    diode.push_dopant(boron);
    diode.push_dopant(phosphorus);

    let sample_count = 4096u32;
    let sample_last_idx = (sample_count - 1) as usize;
    let mesh = Mesh::create((0..sample_count).map(|i| f64::from(i) * (len / f64::from(sample_count))).collect());

    // the derice description is finished
    // solve the device

    // the potential at x = 0 is taken as the reference.
    let fermi_lvl = calcRootBisection(
        diode.bulk.Ev, 
        diode.bulk.Ec, 
        |mu| diode.total_charge(0.0, mu, 0.0, temp), 
        1e-5,
        1e2,
        100
    ).expect("No fermi_lvl in max_iter using bisection method");

    let built_in_potential = calcRootBisection(
        -diode.bulk.band_gap / constants::Q,
        diode.bulk.band_gap / constants::Q,
        |V| diode.total_charge(len, fermi_lvl, V, temp), 
        1e-5,
        1e2,
        100
    ).expect("No fermi_lvl in max_iter using bisection method");

    println!("Solved fermi level: {0} eV", fermi_lvl / constants::Q);
    println!("Solved built-in potential: {0} V", built_in_potential);

    let mut pyviFile = PyVi::create("data.pyvi");

    pyviFile.create_parameter("x", mesh.asVecD());

    pyviFile.create_section("potential", "x");
    pyviFile.create_section("charge", "x");
    pyviFile.create_section("charge derivative", "x");
    
    let poissionProb = PoissionProblem::create(&mesh, mesh.makeVec(diode.bulk.epsilon));

    let mut potential = mesh.zeroVec();
    let mut charge = diode.total_charge_vec(&mesh, fermi_lvl, &potential, temp);
    
    potential[0] = 0.0;
    potential[sample_last_idx] = built_in_potential;

    let charge_derivative = diode.total_charge_derivative_pot_vec(&mesh, fermi_lvl, &potential, temp);

    pyviFile.push_to_section("potential", potential.clone());
    pyviFile.push_to_section("charge", charge.clone());
    pyviFile.push_to_section("charge derivative", charge_derivative);
    
    for _ in 0..100 {
        let charge_derivative = diode.total_charge_derivative_pot_vec(&mesh, fermi_lvl, &potential, temp);

        // calculate the residual
        let mut residual = poissionProb.residue(&potential, &charge);
        residual[0] = 0.0;
        residual[sample_last_idx] = 0.0;

        let poission_mat = poissionProb.operator.clone();
        let mut jacobian:MatTriDiag = (poission_mat.0, &charge_derivative + poission_mat.1, poission_mat.2);

        jacobian.1[0] = 1.0;
        jacobian.2[0] = 0.0;

        jacobian.1[sample_last_idx] = 1.0;
        jacobian.0[sample_last_idx - 1] = 0.0;

        let deltaV = tridiag::solve(&jacobian, &mut mesh.zeroVec(), residual);

        let w = 1.5;
        potential -= w * deltaV;
        potential[sample_last_idx] = built_in_potential;

        charge = diode.total_charge_vec(&mesh, fermi_lvl, &potential, temp);

        pyviFile.push_to_section("potential", potential.clone());
        pyviFile.push_to_section("charge", charge.clone());
        pyviFile.push_to_section("charge derivative", charge_derivative);
    }
}
