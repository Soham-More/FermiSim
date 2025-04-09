use crate::common::*;

#[derive(Debug, Default)]
pub struct State
{
    pub potential:VecD,
    pub p:VecD,
    pub n:VecD,
    pub charge:VecD,
    pub fermi_lvl:f64,
    pub built_in_potential:f64,
    pub Ec:VecD,
    pub Ev:VecD,
}

pub struct TransientFrame
{
    pub state:State,
    pub time:f64,
    pub time_step:f64,
}

