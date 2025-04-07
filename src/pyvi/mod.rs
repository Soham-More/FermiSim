use std::{collections::HashMap, str::FromStr};

use super::common::*;
use std::fs;

struct Section
{
    parameter_name:String,
    value:Vec<VecD>
}

pub struct PyVi
{
    filename:String,
    sections:HashMap<String, Section>,
    parameters:HashMap<String, VecD>
}

fn to_string(name:&str) -> String { String::from_str(name).expect("Error: unexpected error occured!") }

impl PyVi {
    
    pub fn create(filename:&str) -> PyVi
    {
        PyVi { 
            filename:to_string(filename), 
            sections: HashMap::new(), 
            parameters: HashMap::new()
        }
    }

    pub fn create_parameter(&mut self, name:&str, value:VecD)
    {
        self.parameters.insert(to_string(name), value);
    }

    pub fn create_section(&mut self, name:&str, parameter_name:&str)
    {
        self.sections.insert(to_string(name), Section{ parameter_name:to_string(parameter_name), value:Vec::new() });
    }

    pub fn push_to_section(&mut self, name:&str, data:VecD)
    {
        self.sections.get_mut(name).expect("Error: no section named {name}").value.push(data);
    }

    fn vector_as_string(vector:&VecD) -> String
    {
        vector.as_slice().iter()
        .map(|x| format!("{:.17e}", x))
        .collect::<Vec<String>>().join(",")
    }
}

impl Drop for PyVi
{
    fn drop(&mut self) {
        let mut file_contents = String::new();

        file_contents += "[Parameter]\n";

        for (name, value) in &self.parameters
        {
            file_contents += format!("{}:", name).as_str();
            file_contents += Self::vector_as_string(&value).as_str();
            file_contents += "\n";
        }

        file_contents += "[Section]\n";

        for (name, section) in &self.sections
        {
            file_contents += format!("({})->[{}]\n" , name, section.parameter_name).as_str();
            
            for (i, value) in section.value.iter().enumerate()
            {
                file_contents += format!("I[{}]=", i).as_str();
                file_contents += Self::vector_as_string(&value).as_str();
                file_contents += "\n";
            }
            file_contents += "\n";
        }

        fs::write(&self.filename, file_contents).expect("Error: PyVi is unable to write to file");

    }
}

