use crate::common::*;

pub type MatTriDiag = (VecD, VecD, VecD);

// convert the tridiagonal repr to matrix repr
pub fn as_matrix(val:&MatTriDiag) -> na::DMatrix<f64>
{
    let mut mat = na::DMatrix::from_diagonal(&val.1);

    for i in 0..(val.1.len()-1)
    {
        mat[(i, i + 1)] = val.2[i];
    }
    for i in 1..(val.1.len())
    {
        mat[(i, i - 1)] = val.0[i - 1];
    }
    return mat;
}

pub fn apply(A:&MatTriDiag, b:&VecD) -> VecD
{
    let (subdiag, diag, superdiag) = A;
    let N = b.len();
    // checks
    if subdiag.len() != N || diag.len() != N || superdiag.len() != N
    {
        panic!("input dimension mismatch while applying matrix")
    }

    let mut x = VecD::zeros(N);

    // first element
    x[0] = diag[0] * b[0] + superdiag[0] * b[1];
    for i in 1..(N-1)
    {
        x[i] = subdiag[i - 1] * b[i - 1] + diag[i] * b[i] + superdiag[i] * b[i + 1];
    }
    x[N - 1] = subdiag[N - 2] * b[N - 2] + diag[N - 1] * b[N - 1];

    return x;
}

pub fn solve(A:&MatTriDiag, scratch:&mut VecD, mut b:VecD) -> VecD
{
    let (subdiag, diag, superdiag) = A;
    let N = b.len();
    // checks
    if subdiag.len() != N || diag.len() != N || superdiag.len() != N || scratch.len() != N
    {
        panic!("input dimension mismatch in TDMA")
    }

    scratch[0] = superdiag[0] / diag[0];
    b[0] = b[0] / diag[0];

    /* loop from 1 to X - 1 inclusive */
    for ix in 1..=N-1
    {
        if ix < N - 1
        {
            scratch[ix] = superdiag[ix] / (diag[ix] - subdiag[ix - 1] * scratch[ix - 1]);
        }
        b[ix] = (b[ix] - subdiag[ix - 1] * b[ix - 1]) / (diag[ix] - subdiag[ix - 1] * scratch[ix - 1]);
    }

    /* loop from X - 2 to 0 inclusive */
    for ix in (0..=(N-2)).rev()
    {
        b[ix] -= scratch[ix] * b[ix + 1];
    }

    return b;
}

