// Solve poission equation
#![allow(non_snake_case)]

use crate::common::*;
use super::mesh::Mesh;
use super::tridiag;
use tridiag::MatTriDiag;

#[derive(Default)]
pub struct PoissionProblem
{
    pub epsilon: VecD,
    pub operator: MatTriDiag,
    pub scratch: VecD,
}

impl PoissionProblem
{
    // make a poission problem from a mesh and epsilon
    pub fn create(mesh:&Mesh, epsilon:VecD) -> PoissionProblem
    {
        let mut subdiag = mesh.zeroVec();
        let mut diag = mesh.zeroVec();
        let mut superdiag = mesh.zeroVec();
        // set the operator
        let h = mesh.calcStepVec();

        for i in 1..mesh.lastIdx()
        {
            let h_avg = 0.5 * (h[i] + h[i + 1]);
            let coeff_f = (epsilon[i + 1] + epsilon[i]) / (2.0 * h[i]);
            let coeff_b = (epsilon[i - 1] + epsilon[i]) / (2.0 * h[i - 1]);

            diag[i] = - (coeff_f + coeff_b) / h_avg;
            superdiag[i] = coeff_f / h_avg;
            subdiag[i - 1] = coeff_b / h_avg;
        }

        // set the first and last element
        diag[0] = 1.0;
        diag[mesh.lastIdx()] = 1.0;

        PoissionProblem{
            epsilon,
            operator: (subdiag, diag, superdiag),
            scratch:mesh.zeroVec()
        }
    }

    pub fn solve(&mut self, charge:&VecD, left_bc:f64, right_bc:f64) -> VecD
    {
        let mut load_vector = -charge.clone();
        
        // apply boundary conditions
        load_vector[0] = left_bc;
        load_vector[charge.len() - 1] = right_bc;

        return tridiag::solve(&self.operator, &mut self.scratch, load_vector);
    }

    pub fn residue(&self, potential:&VecD, charge:&VecD) -> VecD
    {
        return tridiag::apply(&self.operator, &potential) + charge;
    }
}



