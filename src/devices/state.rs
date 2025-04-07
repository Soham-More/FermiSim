use crate::common::*;

#[derive(Debug, Default)]
pub struct State
{
    potential:VecD,
    p:VecD,
    n:VecD,
    charge:VecD,
    fermi_lvl:VecD,
}

pub struct TransientFrame
{
    state:State,
    time:f64,
    time_step:f64,
}

