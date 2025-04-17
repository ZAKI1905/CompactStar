/*
  Model Abstract class
*/

#include <gsl/gsl_math.h>

#include <Zaki/File/VecSaver.hpp>

#include "CompactStar/EOS/Model.hpp"

//==============================================================
//                        Model class
//==============================================================
// Constructor
CompactStar::Model::Model(const Zaki::Math::Range<double>& in_rho_range) 
: Prog("Model"), valid_rho(in_rho_range)
{ }
//--------------------------------------------------------------
/// Destructor
CompactStar::Model::~Model() { }

//--------------------------------------------------------------
void CompactStar::Model::FindEOS(const size_t& eos_pts) 
{
  eos_results.reserve(eos_pts) ;
  // double log_rho_min = log10(7.225e-11) ;
  // double log_rho_max = log10(7.22e-7) ;
  double log_rho_min = log10(valid_rho.min) ;
  double log_rho_max = log10(valid_rho.max) ;
  double step = (log_rho_max - log_rho_min) / eos_pts;

  double rho_i = 1;
  for (double log_rho_i = log_rho_min;  log_rho_i < log_rho_max; log_rho_i += step)
  {
    rho_i = pow(10, log_rho_i) ;
    // eos_results.push_back({EDens(rho_i), Press(rho_i), rho_i}) ;
    eos_results.push_back(EOSRow(rho_i)) ;
  }
}

//--------------------------------------------------------------
void CompactStar::Model::ExportEOS(const Zaki::String::Directory& f_name) const
{
  if(eos_results.size() == 0)
  {
    Z_LOG_ERROR("EOS hasn't been found yet!") ;
    return ;
  }

  Zaki::String::Directory tmp_dir = f_name ;
  if(set_wrk_dir_flag)
    tmp_dir = wrk_dir  + tmp_dir ;

  Zaki::File::VecSaver saver(tmp_dir) ;
  // saver.SetHeader("g/cm^3,    dyne/cm^2, 1/fm^3") ;
  saver.SetHeader(EOSHeader().c_str()) ;
  saver.Export2D(eos_results, "\t ") ;
}

//--------------------------------------------------------------
std::vector<std::vector<double>> CompactStar::Model::GetEOS() const 
{
  if(eos_results.size() == 0)
  {
    Z_LOG_ERROR("EOS hasn't been found yet!") ;
  }
  return eos_results ;
}
//--------------------------------------------------------------
Zaki::Math::Range<double> CompactStar::Model::GetValidEDensRange()
{
  return {EDens(valid_rho.min), EDens(valid_rho.max) } ;
}

//--------------------------------------------------------------
void CompactStar::Model::SetRhoRange(const Zaki::Math::Range<double>& in_range) 
{
  valid_rho = in_range ;
}

//--------------------------------------------------------------
/// This will generate a row of values for EOS
std::vector<double> CompactStar::Model::EOSRow(const double& rho_i) 
{
  return {EDens(rho_i), Press(rho_i), rho_i} ;
}

//--------------------------------------------------------------
/// This will generate the header for EOS
std::string CompactStar::Model::EOSHeader() const
{
  char tmp_header[300] ;
  snprintf(tmp_header, sizeof(tmp_header), "%-13s\t %-13s\t %-13s", 
            "e(g/cm^3)", "p(dyne/cm^2)", "rho(1/fm^3)" ) ;
  return tmp_header ;
}
//==============================================================
