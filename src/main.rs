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
use devices::device::Device;

fn test_poission_solver() {
    let mesh = Mesh::create((0..20).map(|i| f64::from(i) * 1.0).collect());
    let mut poissionProb = PoissionProblem::create(&mesh, &mesh.makeVec(1.0));

    // charge density
    let rho =  mesh.makeVec(1.0);

    println!("rho: {}", rho);
    
    let soln = poissionProb.solve(&rho, -1.0, 1.0);

    println!("Soln: {}", soln);
    println!("Residual: {}", poissionProb.residue(&soln, &rho));
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

    let mut diode = sc::Semiconductor::create(silicon);

    diode.push_dopant(boron);
    diode.push_dopant(phosphorus);

    let sample_count = 4096u32;
    
    let mut device  = Device::create(temp);
    device.push_bulk_layer(diode, len, sample_count);

    device.calc_steady_state(1e2, 1e-5, 100);

    let mut pyviFile = PyVi::create("data.pyvi");

    pyviFile.create_parameter("x", device.mesh.asVecD());

    pyviFile.create_section("potential", "x");
    pyviFile.create_section("charge", "x");
    //pyviFile.create_section("charge derivative", "x");

    pyviFile.push_to_section("potential", device.steady_state.potential.clone());
    pyviFile.push_to_section("charge", device.steady_state.charge.clone());
    //pyviFile.push_to_section("charge derivative", charge_derivative);
    
}
