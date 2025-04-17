/*
  Common Functions
*/

#include <Zaki/Physics/Constants.hpp>
#include <Zaki/File/CSVIterator.hpp> 
#include <Zaki/Util/Logger.hpp>

#include "CompactStar/EOS/Common.hpp"

using namespace Zaki::Physics ;
//--------------------------------------------------------------
namespace CompactStar
{
//--------------------------------------------------------------
extern const double MN = 938.91875434*Zaki::Physics::MEV_2_INV_FM ; // 1/fm 
//==============================================================
/// Input Fermi momentum in 1/fm
/// outputs energy density in 1/fm^4
double FermiEDens(const double& in_ke) 
{
  double mu = sqrt(pow(ELECTRON_M_FM, 2) + in_ke*in_ke) ;
  double out =  mu*in_ke*(mu*mu - pow(ELECTRON_M_FM, 2)/2.);

  out += -pow(ELECTRON_M_FM,4)* log( (mu+in_ke)/ELECTRON_M_FM ) / 2.;

  out *= 1./(4.*M_PI*M_PI) ;

  return out ;
}

//--------------------------------------------------------------
/// Input Fermi momentum in 1/fm
/// outputs pressure in 1/fm^4
double FermiPress(const double& in_ke) 
{
  double mu = sqrt(pow(ELECTRON_M_FM, 2) + in_ke*in_ke) ;
  double out =  mu*in_ke*(mu*mu - 5.*pow(ELECTRON_M_FM, 2)/2.);

  out += 3.*pow(ELECTRON_M_FM, 4)* log( (mu+in_ke)/ELECTRON_M_FM ) / 2.;

  out *= 1./(12.*M_PI*M_PI) ;

  return out ;
}

//--------------------------------------------------------------
// // Extracts two columns of a set of data from a tsv file
// DataSet ColExtractor(const Zaki::String::Directory& f_name,
//   const std::vector<size_t> _idx) 
// {
//   // std::ifstream     file( (wrk_dir + "/" + f_name).Str());
//   std::ifstream     file( f_name.Str().c_str());

//   // Error opening the file
//   if (file.fail()) 
//   {
//     Z_LOG_ERROR("File '"+(f_name).Str() +"' cannot be opened!") ;
//     Z_LOG_ERROR("Importing EOS data failed!") ;
//     exit(EXIT_FAILURE) ;
//     return {};
//   }

//   size_t line_num = 0 ;

//   // Sanity check:
//   for (size_t i = 0; i < _idx.size(); i++)
//   {
//     if( _idx[i] <=0)
//     {
//       Z_LOG_ERROR("Column indices must be positive and non-zero!") ;
//       return {} ;
//     }
//   }
  
//   // if(_x_idx <= 0 || _y_idx <=0)
//   // {
//   //   Z_LOG_ERROR("Column indices must be positive and non-zero!") ;
//   //   return {} ;
//   // }

//   DataSet out_data ;

//   // out_data.x_data.reserve(100) ;
//   // out_data.y_data.reserve(100) ;

//   // Reading the input file
//   for(Zaki::File::CSVIterator loop(file, '\t'); loop != Zaki::File::CSVIterator(); ++loop)
//   {
//     if( (*loop).size() < _x_idx || (*loop).size() < _y_idx)
//     {
//       Z_LOG_ERROR("At least one of the column indices is out of range!") ;
//       break ;
//     }

//     // First line (header)
//     if (line_num == 0)
//     {   
//       out_data.x_label = (*loop)[_x_idx-1].c_str() ;
//       out_data.y_label = (*loop)[_y_idx-1].c_str() ;
//     }
//     else 
//     {
//       // x values
//       out_data.x_data.emplace_back(std::atof((*loop)[_x_idx-1].c_str())) ;

//       // y values
//       out_data.y_data.emplace_back(std::atof((*loop)[_y_idx-1].c_str())) ;

//     }
//     line_num++ ;
//   }

//   Z_LOG_INFO("Columns extracted from: "+ (f_name).Str()+".") ;

//   return out_data;
// }

//--------------------------------------------------------------
} // End of namespace CompactStar
//--------------------------------------------------------------
//==============================================================