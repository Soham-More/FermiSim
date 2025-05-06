#![allow(non_snake_case)]

extern crate nalgebra as na;
extern crate nalgebra_lapack as nalpk;
extern crate rgsl;
pub mod common;
pub mod fdm1D;

pub mod semiconductor;
pub mod devices;
pub mod pyvi;

use rgsl::error;
use rgsl::Value;
use semiconductor as sc;
use common::*;
use pyvi::PyVi;
use devices::device::Device;

fn error_handling(error_str: &str, file: &str, line: u32, error_value: Value) {
    println!("RGSL [{:?}] '{}:{}': {}", error_value, file, line, error_str);
}
fn main() {
    error::set_error_handler(Some(error_handling));

    let len = 1e-4;
    let temp = 300.0;
    let doping_a = 1e17;
    let doping_b = 8e23;

    let len_spacer = 0.0;//0.1e-6;

    let GaAs = sc::Bulk::create_GaAs_300K();

    // Al_x Ga_1-x As, x = 0.5
    let AlGaAs = sc::Bulk::create_AlGaAs_300K(1.0);

    let zinc = sc::Dopant::create_acceptor(
        vec![doping_a, doping_a], 
        vec![0.0, len], 
        interp::Nearest, 
        GaAs.Ev + 0.045*constants::Q, 
        4.0
    );
    let silicon = sc::Dopant::create_donor(
        vec![0.0, doping_b, doping_b], 
        vec![len, len + 2.0*len_spacer, 2.0*len], 
        interp::Nearest, 
        AlGaAs.Ec - 0.045*constants::Q, 
        2.0
    );

    let mut bottom_layer = sc::Semiconductor::create(GaAs);
    bottom_layer.push_dopant(zinc);

    let mut top_layer = sc::Semiconductor::create(AlGaAs);
    top_layer.push_dopant(silicon);

    let sample_count = 64*4096u32;
    
    let mut device  = Device::create(temp);
    device.push_bulk_layer(bottom_layer, len, sample_count);
    device.push_bulk_layer(top_layer, len, sample_count);
// */
/* 
    let silicon = sc::Bulk::create_silicon_300K();

    let boron = sc::Dopant::create_acceptor(
        vec![doping_a, doping_a, 0.0], 
        vec![0.0, 0.39*len, 0.40*len], 
        interp::Nearest, 
        silicon.Ev + 0.045*constants::Q, 
        4.0
    );
    let phosphorus = sc::Dopant::create_donor(
        vec![0.0, doping_b, doping_b], 
        vec![0.49*len, 0.5*len, 1.0*len], 
        interp::Nearest, 
        silicon.Ec - 0.045*constants::Q, 
        2.0
    );

    let mut bottom_layer = sc::Semiconductor::create(silicon);
    bottom_layer.push_dopant(boron);
    bottom_layer.push_dopant(phosphorus);
    let sample_count = 4096u32;
    
    let mut device  = Device::create(temp);
    device.push_bulk_layer(bottom_layer, len, sample_count);
    // */

    device.calc_steady_state(1e1, 1e-8, 500);
    println!("built-in potential: {:.4} V", device.steady_state.built_in_potential);
    println!("Steady State fermi-level(relative to Vaccum): {:.4} eV", device.steady_state.fermi_lvl / constants::Q);

    let mut pyviFile = PyVi::create("data.pyvi");

    pyviFile.create_parameter("x", device.mesh.asVecD());

    pyviFile.create_section("potential", "x");
    pyviFile.create_section("charge", "x");
    pyviFile.create_section("Ec", "x");
    pyviFile.create_section("Ev", "x");
    pyviFile.create_section("Fermi-Level", "x");
    pyviFile.create_section("n", "x");
    pyviFile.create_section("p", "x");
    pyviFile.create_section("doping", "x");
    //pyviFile.create_section("charge derivative", "x");

    pyviFile.push_to_section("potential", device.steady_state.potential.clone());
    pyviFile.push_to_section("charge", device.steady_state.charge.clone());
    pyviFile.push_to_section("Ec", device.steady_state.Ec.clone() / constants::Q);
    pyviFile.push_to_section("Ev", device.steady_state.Ev.clone() / constants::Q);
    pyviFile.push_to_section("Fermi-Level", VecD::from_element(device.mesh.len(), device.steady_state.fermi_lvl / constants::Q));
    pyviFile.push_to_section("n", device.steady_state.n);
    pyviFile.push_to_section("p", device.steady_state.p);

    pyviFile.push_to_section("doping", device.net_doping);
    //pyviFile.push_to_section("charge derivative", charge_derivative);
    
}
