use crate::common::*;

pub struct Mesh
{
    pub points: Vec<f64>,
}

impl Mesh
{
    pub fn len(&self) -> usize { self.points.len() }
    pub fn lastIdx(&self) -> usize { self.points.len() - 1 }
    // make a vector the size of the meshing
    pub fn zeroVec(&self) -> VecD { VecD::zeros(self.points.len()) }
    // make a vector the size of the meshing, with the value v
    pub fn makeVec(&self, v:f64) -> VecD { VecD::from_element(self.points.len(), v) }
    // make a vector by evaluating a function at every meshing point
    pub fn makeVecFn(&self, f:impl Fn(f64, usize) -> f64) -> VecD { VecD::from_fn(self.points.len(), |i, _| { f(self.points[i], i) }) }

    pub fn asVecD(&self) -> VecD
    {
        self.makeVecFn(|x, _| x)
    }


    pub fn calcStepVec(&self) -> VecD
    {
        let mut step = self.zeroVec();
        for i in 0..self.lastIdx()
        {
            step[i] = self.points[i + 1] - self.points[i];
        }
        return step;
    }

    pub fn new() -> Mesh
    {
        Mesh{
            points:Vec::new(),
        }
    }
    pub fn create(points:Vec<f64>) -> Mesh
    {
        Mesh{
            points,
        }
    }

    pub fn extend(&mut self, points:Vec<f64>)
    {
        self.points.extend(points);
    }
}


