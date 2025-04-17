/*
  Main Program to Generate Equation of State Using Fermi_Gas_Many Model
  Author: Mohammadreza Zakeri
  Date: 2023-10-01
  Description:
    This program uses the Fermi_Gas_Many model to generate the EOS data
    for neutron star matter over a specified range of baryon number densities.
    The results include energy density, pressure, and particle fractions.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

// Include the Fermi_Gas_Many header
#include <CompactStar/EOS/Fermi_Gas_Many.hpp>

int main()
{
	Zaki::String::Directory dir(__FILE__) ;

    // Create an instance of the Fermi_Gas_Many model
    CompactStar::Fermi_Gas_Many fermi_gas_model;

    fermi_gas_model.FindEOS(200) ;
  	fermi_gas_model.SetWrkDir(dir.ParentDir() +"/results") ;

  	fermi_gas_model.ExportEOS("EOS/Fermi_gas_many.eos") ;

    return 0;
}