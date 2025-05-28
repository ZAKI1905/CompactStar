// -------------------------------------------------------------
//                        CompOSE_EOS Class
// -------------------------------------------------------------
// Created on October 29, 2022
//
// Purpose: 
//  To simplify working with the files for the equation of 
//  states from the CompOSE website:
//        compose.obspm.fr
// -------------------------------------------------------------
// #include <matplotlibcpp.hpp>

#include <Zaki/Physics/Constants.hpp>
// #include <Zaki/String/String_Basic.hpp>
#include <Zaki/File/CSVIterator.hpp> 
#include <Zaki/Vector/Vector_Basic.hpp>

#include "CompactStar/EOS/CompOSE_EOS.hpp"


//==============================================================
//             CompOSE_EOS Class
//==============================================================
// The mapping between CompOSE codes with labels
// Taken from Table 3.3 in CompOSE manual (v3.00) from
//   https://compose.obspm.fr/manual
const std::map<int, CompactStar::CompOSE_EOS::CompOSE_Particle> CompactStar::CompOSE_EOS::Compose_Dict = {
    {0, {"$e^-$", Zaki::Physics::ELECTRON_M_MEV}},
    {1, {"$\\mu^-$", Zaki::Physics::MUON_M_MEV}},
    {10, {"$n$", Zaki::Physics::NEUTRON_M_MEV}}, 
    {11, {"$p^+$", Zaki::Physics::PROTON_M_MEV}}, 
    {20, {"$\\Delta^-$", Zaki::Physics::DELTA_MINUS_M_MEV}}, 
    {21, {"$\\Delta^0$", Zaki::Physics::DELTA_ZERO_M_MEV}}, 
    {22, {"$\\Delta^+$", Zaki::Physics::DELTA_PLUS_M_MEV}}, 
    {23, {"$\\Delta^{++}$", Zaki::Physics::DELTA_PLUSPLUS_M_MEV}},
    {100, {"$\\Lambda^0$", Zaki::Physics::LAMBDA_ZERO_M_MEV}}, 
    {110, {"$\\Sigma^-$", Zaki::Physics::SIGMA_MINUS_M_MEV}}, 
    {111, {"$\\Sigma^0$", Zaki::Physics::SIGMA_ZERO_M_MEV}}, 
    {112, {"$\\Sigma^+$", Zaki::Physics::SIGMA_PLUS_M_MEV}},
    {120, {"$\\Xi^-$", Zaki::Physics::XSI_MINUS_M_MEV}}, 
    {121, {"$\\Xi^0$", Zaki::Physics::XSI_ZERO_M_MEV}}, 
    {500, {"$u$", Zaki::Physics::UP_QUARK_M_MEV}}, 
    {501, {"$d$", Zaki::Physics::DOWN_QUARK_M_MEV}}, 
    {502, {"$s$", Zaki::Physics::STRANGE_QUARK_M_MEV}},
    {600, {"$\\gamma$", 0}}
    } ;
//--------------------------------------------------------------
// Constructor
CompactStar::CompOSE_EOS::CompOSE_EOS()
{
  eos.Reserve(3, 350) ;

  eos[(int)EOS_Idx::e].label = "e(g/cm^3)" ;
  eos[(int)EOS_Idx::p].label = "p(dyne/cm^2)" ;
  eos[(int)EOS_Idx::n].label = "rho(1/fm^3)" ;

}

//--------------------------------------------------------------
/// Destructor
CompactStar::CompOSE_EOS::~CompOSE_EOS() 
{ } 
//--------------------------------------------------------------
// Sets the crust-core transition density
void CompactStar::CompOSE_EOS::SetCrustCoreXDensity(
                              const double& transition_density) 
{
  has_crust = true ;
  crust_core_x_den= transition_density ;
}

//--------------------------------------------------------------
/// Imports the EoS grid data
void CompactStar::CompOSE_EOS::ImportGrid(
                          const Zaki::String::Directory& in_dir, 
                          const CompactStar::Dir_Type in_dir_type ) 
{
  Z_LOG_INFO("Importing the density grid file...") ;

  // .........................................................
  //                  Opening the file
  // .........................................................
  Zaki::String::Directory file_dir = in_dir ;

  if (in_dir_type == Dir_Type::relative)
  {
    file_dir = wrk_dir + "/" + in_dir ;
  }
  
  std::ifstream     file( file_dir.Str() ) ;

  // Error opening the file
  if ( file.fail() ) 
  {
    Z_LOG_ERROR("File '"+ file_dir.Str() +"' cannot be opened!") ;
    Z_LOG_ERROR("Importing file failed!") ;
    exit(EXIT_FAILURE) ;
    return ;
  }

  // .........................................................
  // Clearing the data_set:
  eos[(int)EOS_Idx::n].vals.clear() ;
  // .........................................................
  //                  Reading the file
  // .........................................................
  size_t line_num = 0 ;
  
  bool crust_core_x_found = false ;
  for(Zaki::File::CSVIterator loop(file, ' '); loop != Zaki::File::CSVIterator(); ++loop)
  {
    // ---------------------------------------
    //                First line
    // ---------------------------------------
    if (line_num == 0)
    { 
      // Move on to the next line
      line_num++ ;
      continue ;
    }
    // ---------------------------------------
    //                Second line
    // ---------------------------------------
    if (line_num == 1)
    {
      size_t  col_counter = 0 ;
      for (size_t i = 0; i < (*loop).size(); i++)
      {
        // Ignore empty spaces
        // if ((*loop)[i] == "")
        if( (*loop)[i].empty() )
        {
          continue;
        }

        col_counter++ ;

        if (col_counter == 1)
        {        
          // Reserve the column size
          eos[(int)EOS_Idx::n].Reserve( std::atof((*loop)[i].c_str()) ) ;
          grid_size = std::atof((*loop)[i].c_str())  ;
          break;
        }
      }
      // Move on to the next line
      line_num++ ;
      continue ;
    }
    // ---------------------------------------
    //     Importing the rest of the lines
    // ---------------------------------------
    size_t  col_counter = 0 ;
    for (size_t i = 0; i < (*loop).size(); i++)
    {
      // Ignore empty spaces
      // if ((*loop)[i] == "")
      if( (*loop)[i].empty() )
      {
        continue;
      }

      col_counter++ ;
      if (col_counter == 1)
      {
        double tmp_val = std::atof((*loop)[i].c_str()) ;
        eos[(int)EOS_Idx::n].vals.emplace_back( tmp_val ) ;

        // If there is a crust-core transition, find 
        //  at what index ?
        if (tmp_val >= crust_core_x_den && has_crust && !crust_core_x_found)
        {
          crust_core_x_idx = eos[(int)EOS_Idx::n].vals.size()-1 ;
          // std::cout << "\n\t tmp_val = " << tmp_val << "\n" ;
          // std::cout << "\n\t eos[(int)EOS_Idx::n].vals[xidx] = " << eos[(int)EOS_Idx::n].vals[crust_core_x_idx] << "\n" ;
          // std::cout << "\n\t Index = " << crust_core_x_idx << "\n" ;
          crust_core_x_found = true ;
        }
        break;
      }
    }

    line_num++ ;
    // ---------------------------------------
  }

  // std::cout << "\n\t eos[n].vals.size() = " << eos[(int)EOS_Idx::n].vals.size() << "\n" ;
}

//--------------------------------------------------------------
// Imports the ".thermo" file
void CompactStar::CompOSE_EOS::ImportThermo(
                          const Zaki::String::Directory& in_dir, 
                          const bool& gen_plots, 
                          const CompactStar::Dir_Type in_dir_type ) 
{
  Z_LOG_INFO("Importing the '.thermo' file...") ;

  // .........................................................
  //                  Opening the file
  // .........................................................
  Zaki::String::Directory file_dir = in_dir ;

  if (in_dir_type == Dir_Type::relative)
  {
    file_dir = wrk_dir + "/" + in_dir ;
  }
  
  std::ifstream     file( file_dir.Str() ) ;

  // Error opening the file
  if ( file.fail() ) 
  {
    Z_LOG_ERROR("File '"+ file_dir.Str() +"' cannot be opened!") ;
    Z_LOG_ERROR("Importing file failed!") ;
    exit(EXIT_FAILURE) ;
    return ;
  }

  // .........................................................
  // Clearing the data_set:
  eos[(int)EOS_Idx::e].vals.clear() ;
  eos[(int)EOS_Idx::p].vals.clear() ;
  // .........................................................
  //                  Reading the file
  // .........................................................
  size_t line_num = 0 ;
  
  for(Zaki::File::CSVIterator loop(file, ' '); loop != Zaki::File::CSVIterator(); ++loop)
  { 
    // ---------------------------------------
    //                First line
    // ---------------------------------------
    if (line_num == 0)
    { 
      size_t  first_line_col_counter = 0 ;
      for (size_t i = 0; i < (*loop).size(); i++)
      {
        // if ((*loop)[i] == "")
        if( (*loop)[i].empty() )
        {
          continue;
        }

        first_line_col_counter++ ;

        switch (first_line_col_counter)
        {
        case 1:
          m_n = std::atof((*loop)[i].c_str()) ;
          break;

        case 2:
          m_p = std::atof((*loop)[i].c_str()) ;
          break;

        default:
          break;
        }
      }


      // std::cout << " \t m_n = " << m_n << "\n" ;
      // std::cout << " \t m_p = " << m_p << "\n" ;

      // Move on to the next line
      line_num++ ;
      continue ;
    }
    // ---------------------------------------
    //     Importing the rest of the lines
    // ---------------------------------------
    size_t  col_counter = 0 ;
    for (size_t i = 0; i < (*loop).size(); i++)
    {
      // if ((*loop)[i] == "")
      if( (*loop)[i].empty() )
      {
        continue ;
      }
      col_counter++ ;

      switch (col_counter)
      {
        // Pressure column
        case 4:
        {
          double tmp_p = std::atof((*loop)[i].c_str()) ;
          tmp_p *= eos[(int)EOS_Idx::n][line_num-1] ;
          tmp_p *= Zaki::Physics::MEV_FM3_2_Dyn_CM2 ;

          // std::cout << "\t tmp_p= " << tmp_p << "\n" ;

          eos[(int)EOS_Idx::p].vals.emplace_back( tmp_p ) ;
          break;
        }
        // Energy density column
        case 10:
        {
          // std::cout << "\t col_counter = 10: " << (*loop)[i] << "\n" ;
          double tmp_e  = std::atof((*loop)[i].c_str()) + 1 ;
          tmp_e *= eos[(int)EOS_Idx::n][line_num-1] ;

          // Checking if we are in the crust:
          if (line_num-1 < crust_core_x_idx )
          {
            tmp_e *= 939.5653 ;
          }
          else
          {
            tmp_e *= m_n ; 
          }
          tmp_e *= Zaki::Physics::MEV_FM3_2_G_CM3 ;
          
          eos[(int)EOS_Idx::e].vals.emplace_back( tmp_e ) ;
          break;
        }
        default:
          break;
      }
    }

    line_num++ ;
    // ---------------------------------------
  }

  if(gen_plots)
  {
    Zaki::Vector::DataSet::PlotParam plt_par ;
    plt_par.SetGrid() ;
    plt_par.SetXAxisLabel("$\\varepsilon\\, (\\, g \\cdot {\\rm cm}^{-3}\\, )$") ;
    plt_par.SetYAxisLabel("$p\\, (\\, {\\rm dyne } \\cdot {\\rm cm}^{-2}\\, )$") ;
    eos.SetPlotPars(plt_par) ;
    eos.SetWrkDir(file_dir.ThisFileDir()) ;
    eos.LogLogPlot((int)EOS_Idx::e, (int)EOS_Idx::p, "Pressure_vs_Energy.pdf", name ) ;
  }
}
//--------------------------------------------------------------
// Imports the ".compo" file
void CompactStar::CompOSE_EOS::ImportCompo(
                          const Zaki::String::Directory& in_dir, 
                          const bool& gen_plots,
                          const CompactStar::Dir_Type in_dir_type ) 
{
  Z_LOG_INFO("Importing the '.compo' file...") ;

  // .........................................................
  //                  Opening the file
  // .........................................................
  Zaki::String::Directory file_dir = in_dir ;

  if (in_dir_type == Dir_Type::relative)
  {
    file_dir = wrk_dir + "/" + in_dir ;
  }
  
  std::ifstream     file( file_dir.Str() ) ;

  // Error opening the file
  if ( file.fail() ) 
  {
    Z_LOG_ERROR("File '"+ file_dir.Str() +"' cannot be opened!") ;
    Z_LOG_ERROR("Importing file failed!") ;
    exit(EXIT_FAILURE) ;
    return ;
  }

  // .........................................................
  // Clearing the composition part of the data_set:
  for (size_t i = therm_size; i < comp_size + therm_size ; i++)
  {
    eos[i].vals.clear() ;
  }
  
  // .........................................................
  //                  Reading the file
  // .........................................................
  size_t line_num = 0 ;
  
  for(Zaki::File::CSVIterator loop(file, ' '); loop != Zaki::File::CSVIterator(); ++loop)
  {
    // ---------------------------------------
    //                First line
    // ---------------------------------------
    if (line_num == 0)
    { 
      size_t  col_counter = 0 ;
      for (size_t i = 0; i < (*loop).size(); i++)
      {
        // if ((*loop)[i] == "")
        if( (*loop)[i].empty() )
        {
          continue ;
        }
        col_counter++ ;

        if (col_counter == 5)
        {
          // Setting the number of composition columns
          comp_size = std::atof((*loop)[i].c_str()) ;
        }

        // Setting the labels
        if ( col_counter > 5 && col_counter <= 2*comp_size + 5 && col_counter%2 == 0 )
          eos.AddColumn((*loop)[i]) ;
      }
      // Don't move on to the next line
    }
    // ---------------------------------------
    //     Importing all of the lines
    // ---------------------------------------
    size_t  col_counter = 0 ;
    for (size_t i = 0; i < (*loop).size(); i++)
    {
      // if ((*loop)[i] == "")
      if( (*loop)[i].empty() )
      {
        continue ;
      }
      col_counter++ ;

      if ( col_counter > 5 && col_counter <= 2*comp_size + 5 && col_counter%2 == 1 )
      {
        eos[ therm_size + (col_counter - 5)/2 - 1 ][line_num] = std::atof((*loop)[i].c_str()) ;
      }
    }

    line_num++ ;
    // ---------------------------------------
  }

  eos.SetWrkDir(file_dir.ThisFileDir()) ;
  eos.Export(name + ".eos") ;
  
  if(gen_plots)
  {
    std::vector<std::pair<int, std::string>> tmp_plt_idx ;

    for (int i = therm_size; i < comp_size + therm_size ; i++)
    {
      // Making the plot more legible by removing unimportant contributions
      if ( eos[i].Max() < 1e-4 )
      {
        continue ;
      }

      tmp_plt_idx.emplace_back(i, Compose_Dict.at( std::stoi( eos[i].label.c_str() ) ).name ) ;
    }

    // eos.ResetPlotPars() ;
    Zaki::Vector::DataSet::PlotParam plt_par ;
    plt_par.SetGrid() ;
    plt_par.SetXAxisLabel("$n\\, ( {\\rm fm}^{-3} )$") ;
    plt_par.SetYAxisLabel("$f_i$") ;
    plt_par.SetLegend({"upper left", 0.0, 1.0}) ;
    // plt_par.SetXAxis({1e-3, 2}) ;
    plt_par.SetYAxis({1e-3, 1}) ;

    eos.SetPlotPars(plt_par) ;

    eos.LogLogPlot((int)EOS_Idx::n, tmp_plt_idx, "Composition_vs_Density.pdf", name ) ;
  }
}
//--------------------------------------------------------------
// Imports the ".micro" file
void CompactStar::CompOSE_EOS::ImportMicro(
                          const Zaki::String::Directory& in_dir,
                          const bool& gen_plots, 
                          const CompactStar::Dir_Type in_dir_type ) 
{
  Z_LOG_INFO("Importing the '.micro' file...") ;

  // .........................................................
  //                  Opening the file
  // .........................................................
  Zaki::String::Directory file_dir = in_dir ;

  if (in_dir_type == Dir_Type::relative)
  {
    file_dir = wrk_dir + "/" + in_dir ;
  }
  
  std::ifstream     file( file_dir.Str() ) ;

  // Error opening the file
  if ( file.fail() ) 
  {
    Z_LOG_ERROR("File '"+ file_dir.Str() +"' cannot be opened!") ;
    Z_LOG_ERROR("Importing file failed!") ;
    exit(EXIT_FAILURE) ;
    return ;
  }

  // .........................................................
  // Clearing the m_eff, V_eff, U datasets
  for (size_t i = 0; i < m_eff.Dim().size() ; i++)
  {
    m_eff[i].vals.clear() ;
  }
  
  for (size_t i = 0; i < V_eff.Dim().size() ; i++)
  {
    V_eff[i].vals.clear() ;
  }

  for (size_t i = 0; i < U.Dim().size() ; i++)
  {
    U[i].vals.clear() ;
  }

  // .........................................................
  //                  Reading the file
  // .........................................................
  size_t line_num = 0 ;
  std::vector<size_t> m_eff_col_set ;
  std::vector<size_t> V_eff_col_set ;
  std::vector<size_t> U_col_set ;

  for(Zaki::File::CSVIterator loop(file, ' '); loop != Zaki::File::CSVIterator(); ++loop)
  {
    // size_t  m_eff_counter = 0 ;

    // ---------------------------------------
    //                First line
    // ---------------------------------------
    if (line_num == 0)
    { 
      // Adding the first label (density)
      std::vector<std::string> m_eff_labels = {eos[(int)EOS_Idx::n].label } ;
      std::vector<std::string> V_eff_labels = {eos[(int)EOS_Idx::n].label } ;
      std::vector<std::string> U_labels = {eos[(int)EOS_Idx::n].label } ;

      size_t  col_counter = 0 ;
      for (size_t i = 0; i < (*loop).size(); i++)
      {
        // if ((*loop)[i] == "")
        if( (*loop)[i].empty() )
        {
          continue ;
        }
        col_counter++ ;

        // if (col_counter == 4)
        // {
        //   // Setting the number of micro columns
        //   comp_size = std::atof((*loop)[i].c_str()) ;
        // }

        // ..............................................
        // Setting the labels from odd columns after 4
        // that match the m_eff, V_eff and U labels
        if ( col_counter > 4 && col_counter%2 == 1 )
        {
          // std::cout << "\n\t All: " << std::atoi((*loop)[i].c_str()) << " " ;

          // .....................
          // Code for m_eff:
          // K = I * 1000 + 41
          if ( Zaki::String::EndsWith((*loop)[i], "41") 
                &&
                Compose_Dict.count( (std::atoi((*loop)[i].c_str()) - 41) / 1000 ) != 0 
              )
          {
            m_eff_col_set.emplace_back(col_counter) ;
            m_eff_labels.emplace_back(std::to_string((std::atoi((*loop)[i].c_str()) - 41) / 1000)) ;
          }
          // .....................
          // Code for V_eff:
          // K = I * 1000 + 51
          if ( Zaki::String::EndsWith((*loop)[i], "51") 
                &&
                Compose_Dict.count( (std::atoi((*loop)[i].c_str()) - 51) / 1000 ) != 0 
              )
          {
            V_eff_col_set.emplace_back(col_counter) ;
            V_eff_labels.emplace_back(std::to_string((std::atoi((*loop)[i].c_str()) - 51) / 1000)) ;
          }
          // .....................
          // Code for U:
          // K = I * 1000 + 50
          if ( Zaki::String::EndsWith((*loop)[i], "50") 
                &&
                Compose_Dict.count( (std::atoi((*loop)[i].c_str()) - 50) / 1000 ) != 0 
              )
          {
            U_col_set.emplace_back(col_counter) ;
            U_labels.emplace_back(std::to_string((std::atoi((*loop)[i].c_str()) - 50) / 1000)) ;
          }
          // .....................
        }
      }
      // Don't move on to the next line

      // We now know the number of columns in m_eff
      // (+1) is for the number-density column
      m_eff.Reserve(m_eff_col_set.size()+1, eos[(int)EOS_Idx::n].Size()) ;
      for (size_t i = 0; i < m_eff_col_set.size()+1; i++)
      {
        m_eff[i].label = m_eff_labels[i] ;
      }

      V_eff.Reserve(V_eff_col_set.size()+1, eos[(int)EOS_Idx::n].Size()) ;
      for (size_t i = 0; i < V_eff_col_set.size()+1; i++)
      {
        V_eff[i].label = V_eff_labels[i] ;
      }

      U.Reserve(U_col_set.size()+1, eos[(int)EOS_Idx::n].Size()) ;
      for (size_t i = 0; i < U_col_set.size()+1; i++)
      {
        U[i].label = U_labels[i] ;
      }
    }

    // ---------------------------------------
    //     Importing all of the lines
    // ---------------------------------------
    size_t  col_counter = 0 ;

    // ...........................
    // Temporary row values for m_eff, V_eff and U
    std::vector<double> tmp_m_eff_row ;
    std::vector<double> tmp_V_eff_row ;
    std::vector<double> tmp_U_row ;
    // ...........................

    for (size_t i = 0; i < (*loop).size(); i++)
    {
      // if ((*loop)[i] == "")
      if( (*loop)[i].empty() )
      {
        continue ;
      }
      col_counter++ ;

      // .........................
      // Add the density as the first element
      if (col_counter == 1) // Not an empty line
      {
        tmp_m_eff_row.emplace_back( eos[(int)EOS_Idx::n][line_num] ) ;
        tmp_V_eff_row.emplace_back( eos[(int)EOS_Idx::n][line_num] ) ;
        tmp_U_row.emplace_back( eos[(int)EOS_Idx::n][line_num] ) ;
      }
      // .........................
      // Specific columns that correspond to m_eff
      if ( Zaki::Vector::Exists(col_counter, m_eff_col_set) )
      {
        size_t c_i = i + 1 ;
        while ((*loop)[c_i].empty())
        {
          c_i++ ;
        }
        // The effective mass divided by the particle mass is printed in Compose micro tables
        // so we have to multiply by the particle's mass to get the effective mass
        tmp_m_eff_row.emplace_back( std::atof((*loop)[c_i].c_str()) * 
          Compose_Dict.at( (std::atoi((*loop)[i].c_str()) - 41) / 1000 ).m ) ;
        // std::cout << "\n\t -> " << std::atof((*loop)[i].c_str()) << " " ;
        // m_eff[ (col_counter - 5)/2 ][line_num] = std::atof((*loop)[i].c_str()) ;
      }
      // .........................
      // Specific columns that correspond to V_eff
      // Simpler than finding m_eff because the V_eff is just printed as it is.
      if ( Zaki::Vector::Exists(col_counter, V_eff_col_set) )
      {
        size_t c_i = i + 1 ;
        while ((*loop)[c_i].empty())
        {
          c_i++ ;
        }
        tmp_V_eff_row.emplace_back( std::atof((*loop)[c_i].c_str()) ) ;
      }
      // .........................
      // Specific columns that correspond to U
      if ( Zaki::Vector::Exists(col_counter, U_col_set) )
      {
        // std::cout << "\ncol_counter = " << col_counter 
        //           << ", '" << (*loop)[i+1].c_str() << "' " ;
        size_t c_i = i + 1 ;
        while ((*loop)[c_i].empty())
        {
          // std::cout << "\n empty! " << " " ;
          c_i++ ;
        }
        // std::cout << "\ncol_counter = " << col_counter 
        //           << ", '" << (*loop)[c_i].c_str() << "' " ;
        tmp_U_row.emplace_back( std::atof((*loop)[c_i].c_str()) ) ;
      }
      // .........................
    }

    if (col_counter == 0) // Empty lines
    {
      continue;
    }

    // ...........................
    // Transferring the temporary rows to the datasets
    m_eff.AppendRow(tmp_m_eff_row) ;
    V_eff.AppendRow(tmp_V_eff_row) ;
    U.AppendRow(tmp_U_row) ;
    // ...........................
    // Clearing the temporary rows
    tmp_m_eff_row.clear() ;
    tmp_V_eff_row.clear() ;
    tmp_U_row.clear() ;
    // ...........................

    line_num++ ;
    // ---------------------------------------
  }

  // ...........................
  // Exporting the datasets
  if(m_eff_col_set.size() != 0)
  { 
    //     m_eff
    m_eff.SetWrkDir(file_dir.ThisFileDir()) ;
    m_eff.Export(name + "_m_eff.micro") ;
  }
  if(V_eff_col_set.size() != 0)
  {
    //     V_eff
    V_eff.SetWrkDir(file_dir.ThisFileDir()) ;
    V_eff.Export(name + "_V.micro") ;
  }
  if(U_col_set.size() !=0)
  {
    //     U
    U.SetWrkDir(file_dir.ThisFileDir()) ;
    U.Export(name + "_U.micro") ;
  }
  // ...........................


  if(gen_plots)
  {

    std::vector<std::pair<int, std::string>> tmp_plt_idx ;
    Zaki::Vector::DataSet::PlotParam plt_par ;

    plt_par.SetGrid() ;
    plt_par.SetXAxisLabel("$n\\,\\, ( {\\rm fm}^{-3} )$") ;
    // .................................................
    //                  Plotting M_eff
    // .................................................   
    if(m_eff_col_set.size() != 0)
    { 

      for (int i = 1; i < m_eff.Dim().size() ; i++)
      {
        // Making the plot more legible by removing unimportant contributions
        // Remove the ones with almost constant masses
        if ( (m_eff[i].Max() - m_eff[i].Min() ) / m_eff[i].Max() < 0.05 )
        {
          continue ;
        }

        tmp_plt_idx.emplace_back(i, Compose_Dict.at( std::stoi( m_eff[i].label.c_str() ) ).name ) ;
      }

      plt_par.SetYAxisLabel("$m_B^* \\,\\, ( {\\rm MeV} )$") ;
      plt_par.SetLegend({"lower left", 0.95, 0.0}) ;
      // plt_par.SetXAxis({1e-3, 2}) ;
      // plt_par.SetYAxis({200, 1200}) ;

      m_eff.SetPlotPars(plt_par) ;
      m_eff.Plot(0, tmp_plt_idx, "Meff_vs_Density.pdf", name ) ;
    }
    // ................................................. 


    // .................................................
    //                  Plotting V_eff
    // .................................................   
    if(V_eff_col_set.size() != 0)
    {
      tmp_plt_idx.clear() ;

      for (int i = 1; i < V_eff.Dim().size() ; i++)
      {
        // Making the plot more legible by removing unimportant contributions
        // Remove the ones that are zero (e.g., leptons)
        if ( V_eff[i].Max() == 0 )
        {
          continue ;
        }

        tmp_plt_idx.emplace_back(i, Compose_Dict.at( std::stoi( V_eff[i].label.c_str() ) ).name ) ;
      }

      plt_par.SetYAxisLabel("$\\Sigma^0_B\\,\\, ( {\\rm MeV} )$") ;
      plt_par.SetLegend({"lower right", 1.0, 0.0}) ;
      // plt_par.SetXAxis({1e-3, 2}) ;
      // plt_par.SetYAxis({200, 1200}) ;

      V_eff.SetPlotPars(plt_par) ;
      V_eff.LogLogPlot(0, tmp_plt_idx, "V_vs_Density.pdf", name ) ;
    }
    // ................................................. 

    // .................................................
    //                  Plotting U
    // .................................................   
    if(U_col_set.size() != 0)
    {
      tmp_plt_idx.clear() ;

      for (int i = 1; i < U.Dim().size() ; i++)
      {
        // Making the plot more legible by removing unimportant contributions
        // Remove the ones that are zero (e.g., leptons)
        if ( U[i].Max() == 0 )
        {
          continue ;
        }

        tmp_plt_idx.emplace_back(i, Compose_Dict.at( std::stoi( U[i].label.c_str() ) ).name ) ;
      }

      plt_par.SetYAxisLabel("$U_B\\,\\, ( {\\rm MeV} )$") ;
      plt_par.SetLegend({"lower right", 1.0, 0.0}) ;
      // plt_par.SetXAxis({1e-3, 2}) ;
      // plt_par.SetYAxis({200, 1200}) ;

      U.SetPlotPars(plt_par) ;
      U.Plot(0, tmp_plt_idx, "U_vs_Density.pdf", name ) ;
    }
    // ................................................. 
  }
}
//--------------------------------------------------------------
// Extends the M_eff file to the crust region by
// linearly extrapolating
// The path is by default absolute (not relative to wrk_dir) !
void CompactStar::CompOSE_EOS::ExtendMeffToCrust(
                const Zaki::String::Directory& in_dir, 
                const std::string& f_name,
                const bool& gen_plots, 
                const Dir_Type in_dir_type ) 
{
  // std::cout << "\n dir = " << in_dir + f_name <<  "\n\n" ;

  Zaki::Vector::DataSet m_eff_no_crust(in_dir, f_name) ;

  double n_0 = GetEOS(2).Min() ;
  double n_1 = m_eff_no_crust[0][0] ;

  // Clearing the m_eff
  for (size_t i = 0; i < m_eff.Dim().size() ; i++)
  {
    m_eff[i].vals.clear() ;
  }
  
  // std::cout << "\n m_eff_no_crust.Dim().size() = " << m_eff_no_crust.Dim().size() << "\n \n"  ;

  m_eff.Reserve(m_eff_no_crust.Dim().size(), GetEOS(2).Size()) ;
  
  m_eff[0] = GetEOS(2).GetSubSet([](const double& v){ return v < 0.03 ;}) ;
  m_eff[0].label = "rho(1/fm^3)" ;
  for (size_t i = 1; i < m_eff_no_crust.Dim().size() ; i++)
  {
    m_eff[i].label = m_eff_no_crust[i].label ;
    
    // std::cout << "\n m_eff[i].label = " << m_eff[i].label << " "  ;

    double m_1 = m_eff_no_crust[i][0] ;
    int particle_label = std::stoi( m_eff_no_crust[i].label.c_str() ) ;
    double m_0 = Compose_Dict.at( particle_label ).m ;
    std::string name = Compose_Dict.at( particle_label ).name ;
    double slope = ( m_1 - m_0 ) / ( n_1 - n_0 ) ;
    double interc = m_1 - slope * n_1 ;

    // std::cout << "\n i = " << i << " "  ;
    // std::cout << " name = '"<< name << ", m_0 = " << m_0 
    //         << ", m_1 = " << m_1 << ", slope = " << slope 
    //         << ", y-inter = " << interc ;
    
    for (size_t j = 0; j < m_eff[0].Size(); j++)
    { 
      double m = slope * m_eff[0][j] + interc ;
      m_eff[i].vals.emplace_back(m) ;
    }
  }

  m_eff = m_eff.Append(m_eff_no_crust) ;

  // -------------------------------------------------------------
  // The EoS density doesn't necessarily match with m_eff so we will 
  // interpolate and evaluate at EoS densities just to be safe:
  // Zaki::Vector::DataSet
  for (size_t i = 1; i < m_eff.Dim().size(); i++)
  {
    m_eff.Interpolate(0, i) ;
    m_eff[i] = m_eff.Evaluate(i, GetEOS(2)) ;
    m_eff[i].label = m_eff_no_crust[i].label ;
  }
  m_eff[0] = GetEOS(2) ;
  // -------------------------------------------------------------

  m_eff.SetWrkDir(wrk_dir) ;
  m_eff.Export(name + "_m_eff.micro") ;

  if(gen_plots)
  {
    std::vector<std::pair<int, std::string>> tmp_plt_idx ;
    Zaki::Vector::DataSet::PlotParam plt_par ;

    plt_par.SetGrid() ;
    plt_par.SetXAxisLabel("$n\\,\\, ( {\\rm fm}^{-3} )$") ;
    // .................................................
    //                  Plotting M_eff
    // .................................................   
    for (int i = 1; i < m_eff.Dim().size() ; i++)
    {
      // Making the plot more legible by removing unimportant contributions
      // Remove the ones with almost constant masses
      if ( (m_eff[i].Max() - m_eff[i].Min() ) / m_eff[i].Max() < 0.05 )
      {
        continue ;
      }

      tmp_plt_idx.emplace_back(i, Compose_Dict.at( std::stoi( m_eff[i].label.c_str() ) ).name ) ;
    }

    plt_par.SetYAxisLabel("$m_B^* \\,\\, ( {\\rm MeV} )$") ;
    plt_par.SetLegend({"lower left", 0.95, 0.0}) ;
    // plt_par.SetXAxis({1e-3, 2}) ;
    // plt_par.SetYAxis({200, 1200}) ;

    m_eff.SetPlotPars(plt_par) ;
    m_eff.Plot(0, tmp_plt_idx, "Meff_vs_Density.pdf", name ) ;
    // ................................................. 
  }
}

//--------------------------------------------------------------
// Extends the V_eff file to the crust region by
// linearly extrapolating
// The path is by default absolute (not relative to wrk_dir) !
void CompactStar::CompOSE_EOS::ExtendVeffToCrust(
                const Zaki::String::Directory& in_dir, 
                const std::string& f_name,
                const bool& gen_plots, 
                const Dir_Type in_dir_type ) 
{
  Zaki::Vector::DataSet V_eff_no_crust(in_dir, f_name) ;

  double n_0 = GetEOS(2).Min() ;
  double n_1 = V_eff_no_crust[0][0] ;

  // Clearing the V_eff
  for (size_t i = 0; i < V_eff.Dim().size() ; i++)
  {
    V_eff[i].vals.clear() ;
  }

  V_eff.Reserve(V_eff_no_crust.Dim().size(), GetEOS(2).Size()) ;
  
  V_eff[0] = GetEOS(2).GetSubSet([](const double& v){ return v < 0.03 ;}) ;
  V_eff[0].label = "rho(1/fm^3)" ;

  for (size_t i = 1; i < V_eff_no_crust.Dim().size() ; i++)
  {
    V_eff[i].label = V_eff_no_crust[i].label ;


    double V_1 = V_eff_no_crust[i][0] ;
    double V_0 = 0 ;

    double slope = ( V_1 - V_0 ) / ( n_1 - n_0 ) ;
    double interc = V_1 - slope * n_1 ;

    // std::cout << "\n i = " << i << " "  ;
    // std::cout <<" V_0 = " << V_0 
            // << ", V_1 = " << V_1 << ", slope = " << slope 
            // << ", y-inter = " << interc ;
    
    for (size_t j = 0; j < V_eff[0].Size(); j++)
    { 
      double v = slope * V_eff[0][j] + interc ;
      V_eff[i].vals.emplace_back(v) ;
    }
  }

  V_eff = V_eff.Append(V_eff_no_crust) ;

  // -------------------------------------------------------------
  // The EoS density doesn't necessarily match with m_eff so we will 
  // interpolate and evaluate at EoS densities just to be safe:
  // Zaki::Vector::DataSet
  for (size_t i = 1; i < V_eff.Dim().size(); i++)
  {
    V_eff.Interpolate(0, i) ;
    V_eff[i] = V_eff.Evaluate(i, GetEOS(2)) ;
    V_eff[i].label = V_eff_no_crust[i].label ;
  }
  V_eff[0] = GetEOS(2) ;
  // -------------------------------------------------------------

  V_eff.SetWrkDir(wrk_dir) ;
  V_eff.Export(name + "_V.micro") ;

  if(gen_plots)
  {
    std::vector<std::pair<int, std::string>> tmp_plt_idx ;
    Zaki::Vector::DataSet::PlotParam plt_par ;

    plt_par.SetGrid() ;
    plt_par.SetXAxisLabel("$n\\,\\, ( {\\rm fm}^{-3} )$") ;
    // .................................................
    //                  Plotting V_eff
    // .................................................   
    for (int i = 1; i < V_eff.Dim().size() ; i++)
    {
      // Making the plot more legible by removing unimportant contributions
      // Remove the ones that are zero (e.g., leptons)
      if ( V_eff[i].Max() == 0 )
      {
        continue ;
      }

      tmp_plt_idx.emplace_back(i, Compose_Dict.at( std::stoi( V_eff[i].label.c_str() ) ).name ) ;
    }

    plt_par.SetYAxisLabel("$\\Sigma^0_B\\,\\, ( {\\rm MeV} )$") ;
    plt_par.SetLegend({"lower right", 1.0, 0.0}) ;
    // plt_par.SetXAxis({1e-3, 2}) ;
    // plt_par.SetYAxis({200, 1200}) ;

    V_eff.SetPlotPars(plt_par) ;
    V_eff.Plot(0, tmp_plt_idx, "V_vs_Density.pdf", name ) ;
    // .................................................  
  }
}

//--------------------------------------------------------------
// Extends the U file to the crust region by
// linearly extrapolating
// The path is by default absolute (not relative to wrk_dir) !
void CompactStar::CompOSE_EOS::ExtendUeffToCrust(
                const Zaki::String::Directory& in_dir, 
                const std::string& f_name,
                const bool& gen_plots, 
                const Dir_Type in_dir_type ) 
{
  Zaki::Vector::DataSet U_eff_no_crust(in_dir, f_name) ;

  double n_0 = GetEOS(2).Min() ;
  double n_1 = U_eff_no_crust[0][0] ;

  // Clearing the U_eff
  for (size_t i = 0; i < U.Dim().size() ; i++)
  {
    U[i].vals.clear() ;
  }

  U.Reserve(U_eff_no_crust.Dim().size(), GetEOS(2).Size()) ;
  
  U[0] = GetEOS(2).GetSubSet([](const double& v){ return v < 0.03 ;}) ;
  U[0].label = "rho(1/fm^3)" ;

  for (size_t i = 1; i < U_eff_no_crust.Dim().size() ; i++)
  {
    U[i].label = U_eff_no_crust[i].label ;


    double U_1 = U_eff_no_crust[i][0] ;
    double U_0 = 0 ;

    double slope = ( U_1 - U_0 ) / ( n_1 - n_0 ) ;
    double interc = U_1 - slope * n_1 ;

    // std::cout << "\n i = " << i << " "  ;
    // std::cout <<" V_0 = " << U_0 
    //         << ", V_1 = " << U_1 << ", slope = " << slope 
    //         << ", y-inter = " << interc ;
    
    for (size_t j = 0; j < U[0].Size(); j++)
    { 
      double v = slope * U[0][j] + interc ;
      U[i].vals.emplace_back(v) ;
    }
  }

  U = U.Append(U_eff_no_crust) ;

  // -------------------------------------------------------------
  // The EoS density doesn't necessarily match with m_eff so we will 
  // interpolate and evaluate at EoS densities just to be safe:
  // Zaki::Vector::DataSet
  for (size_t i = 1; i < U.Dim().size(); i++)
  {
    U.Interpolate(0, i) ;
    U[i] = U.Evaluate(i, GetEOS(2)) ;
    U[i].label = U_eff_no_crust[i].label ;
  }
  U[0] = GetEOS(2) ;
  // -------------------------------------------------------------

  U.SetWrkDir(wrk_dir) ;
  U.Export(name + "_U.micro") ;

  if(gen_plots)
  {
    std::vector<std::pair<int, std::string>> tmp_plt_idx ;
    Zaki::Vector::DataSet::PlotParam plt_par ;

    plt_par.SetGrid() ;
    plt_par.SetXAxisLabel("$n\\,\\, ( {\\rm fm}^{-3} )$") ;
   
    // .................................................
    //                  Plotting U
    // .................................................   
    // if(U_col_set.size() != 0)
    // {
      for (int i = 1; i < U.Dim().size() ; i++)
      {
        // Making the plot more legible by removing unimportant contributions
        // Remove the ones that are zero (e.g., leptons)
        if ( U[i].Max() == 0 )
        {
          continue ;
        }

        tmp_plt_idx.emplace_back(i, Compose_Dict.at( std::stoi( U[i].label.c_str() ) ).name ) ;
      }

      plt_par.SetYAxisLabel("$U_B\\,\\, ( {\\rm MeV} )$") ;
      plt_par.SetLegend({"lower right", 1.0, 0.0}) ;
      // plt_par.SetXAxis({1e-3, 2}) ;
      // plt_par.SetYAxis({200, 1200}) ;

      U.SetPlotPars(plt_par) ;
      U.Plot(0, tmp_plt_idx, "U_vs_Density.pdf", name ) ;
    // }
    // ................................................. 
  }
}

//--------------------------------------------------------------
Zaki::Vector::DataSet CompactStar::CompOSE_EOS::GetEOS() const 
{
  return eos ;
}

//--------------------------------------------------------------
Zaki::Vector::DataColumn CompactStar::CompOSE_EOS::GetEOS(const int& idx) const 
{
  return eos[idx] ;
}

//--------------------------------------------------------------
Zaki::Vector::DataColumn CompactStar::CompOSE_EOS::GetEOS(const std::string& label) const 
{
  return eos[label] ;
}

//--------------------------------------------------------------
Zaki::Vector::DataSet CompactStar::CompOSE_EOS::GetMeff() const 
{
  return m_eff ;
}

//-------------------------------------------------------------- 
Zaki::Vector::DataColumn CompactStar::CompOSE_EOS::GetMeff(const int& idx) const 
{
  return m_eff[idx] ;
}

//--------------------------------------------------------------
Zaki::Vector::DataColumn CompactStar::CompOSE_EOS::GetMeff(const std::string& label) const 
{
  return m_eff[label] ;
}

//--------------------------------------------------------------
Zaki::Vector::DataSet CompactStar::CompOSE_EOS::GetVeff() const 
{
  return V_eff ;
}

//--------------------------------------------------------------
Zaki::Vector::DataColumn CompactStar::CompOSE_EOS::GetVeff(const int& idx) const 
{
  return V_eff[idx] ;
}

//--------------------------------------------------------------
Zaki::Vector::DataColumn CompactStar::CompOSE_EOS::GetVeff(const std::string& label) const 
{
  return V_eff[label] ;
}

//--------------------------------------------------------------
/// Returns the single-particle potential dataset
Zaki::Vector::DataSet CompactStar::CompOSE_EOS::GetU() const 
{
  return U ;
}

//--------------------------------------------------------------
/// Returns the single-particle potential from index
Zaki::Vector::DataColumn CompactStar::CompOSE_EOS::GetU(const int& idx) const 
{
  return U[idx] ;
}

//--------------------------------------------------------------
/// Returns the single-particle potential from label
Zaki::Vector::DataColumn CompactStar::CompOSE_EOS::GetU(const std::string& label) const 
{
  return U[label] ;
}
//--------------------------------------------------------------
// Returns the neutron mass constant
double CompactStar::CompOSE_EOS::GetNeutronMass() const 
{
  return m_n ;
}

//--------------------------------------------------------------
// Returns the proton mass constant
double CompactStar::CompOSE_EOS::GetProtonMass() const
{
  return m_p ;
}

//--------------------------------------------------------------
//  Import the EOS in the standard format [ not Compose! ]
// imports everything: e, p, n, composition, m_eff, V_eff, U_eff.
//  in_file: Input file
//  in_dir_type: Relative to wrk_dir or absolute?
void CompactStar::CompOSE_EOS::ImportEOS(
                      const Zaki::String::Directory& in_file, 
                      const bool& gen_plots) 
{
  // std::cout << "\n\t " << wrk_dir + "/" + in_file << "\n" ;
  eos.Import(wrk_dir + "/" + in_file + ".eos") ;
  m_eff.Import(wrk_dir + "/" + in_file + "_m_eff.micro") ;
  V_eff.Import(wrk_dir + "/" + in_file + "_V.micro") ;
  U.Import(wrk_dir + "/" + in_file + "_U.micro") ;
}

//--------------------------------------------------------------
// Plots the Fermi energy of particles as a function of density
// The input directory should by default an absolute path.
void CompactStar::CompOSE_EOS::PlotFermiE(
                const Zaki::String::Directory& in_dir,
                const Dir_Type in_dir_type) const 
{
  // ------------------------------------
  //        Finding Fermi Energy
  // ------------------------------------
  Zaki::Vector::DataColumn fermi_electron =  ( 
    pow(Zaki::Physics::ELECTRON_M_FM, 2) 
    + (3*M_PI*M_PI* eos["0"] * eos[(int)EOS_Idx::n]).pow(2./3.)).sqrt() / Zaki::Physics::MEV_2_INV_FM ;

  Zaki::Vector::DataColumn fermi_muon =  ( 
    pow(Zaki::Physics::MUON_M_FM, 2) 
    + (3*M_PI*M_PI* eos["1"] * eos[(int)EOS_Idx::n]).pow(2./3.)).sqrt() / Zaki::Physics::MEV_2_INV_FM ;

  Zaki::Vector::DataColumn fermi_neutron =  ( 
    m_eff["10"].pow(2) 
    + (3*M_PI*M_PI* eos["10"] * eos[(int)EOS_Idx::n]).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
    ).sqrt() + V_eff["10"]  ;

  Zaki::Vector::DataColumn fermi_lambda =  ( 
    m_eff["100"].pow(2) 
    + (3*M_PI*M_PI* eos["100"] * eos[(int)EOS_Idx::n]).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
    ).sqrt() + V_eff["100"]  ;

  Zaki::Vector::DataSet fermi_ds({eos[(int)EOS_Idx::n], 
                                  fermi_electron, fermi_muon, 
                                  fermi_neutron, fermi_lambda}) ;

  // ------------------------------------
  //              Plotting
  // ------------------------------------
  Zaki::Vector::DataSet::PlotParam plt_par ;

  plt_par.SetGrid() ;
  plt_par.SetXAxisLabel("$n\\,\\, ( {\\rm fm}^{-3} )$") ;
  plt_par.SetYAxisLabel("$E_F\\,\\, ( {\\rm MeV} )$") ;
  plt_par.SetLegend({"upper left", 0.0, 1}) ;
  plt_par.AddAxVLine(6.691204387571e-01, {{"label", "J0348"}}) ; // J0348
  plt_par.AddAxVLine(5.519991401615e-01, {{"ls", "--"}}) ; // J1614
  plt_par.AddAxVLine(3.302221541395e-01, {{"ls", "-."}}) ; // J0737A
  plt_par.AddAxVLine(3.102358492014e-01, {{"ls", ":"}}) ; // J0737B

  fermi_ds.SetPlotPars(plt_par) ;

  Zaki::String::Directory file_dir = in_dir ;
  if (in_dir_type == Dir_Type::relative)
  {
    file_dir = wrk_dir + "/" + in_dir ;
  }

  // std::cout << "\n\t" << file_dir + "/E_Fermi.pdf\n" ;
  // fermi_ds.SemiLogXPlot(0, {{1, "$e^-$"}, {2, "$\\mu^-$"}, {3, "$n$"}, {4, "$\\Lambda$"}}, 
  //                          file_dir + "/" + name + "_E_Fermi.pdf", name) ;
  fermi_ds.SemiLogXPlot(0, {{1, "$e^-$"}}, 
                           file_dir + "/" + name + "_e_E_Fermi.pdf", name) ;
}

//--------------------------------------------------------------
// Exports the Fermi energy of particles as a function of density
// The input directory should by default an absolute path.
void CompactStar::CompOSE_EOS::ExportFermiE(
                const Zaki::String::Directory& in_dir,
                const Dir_Type in_dir_type) const 
{
  // ------------------------------------
  //        Finding Fermi Energy
  // ------------------------------------
  Zaki::Vector::DataColumn fermi_electron =  ( 
    pow(Zaki::Physics::ELECTRON_M_FM, 2) 
    + (3*M_PI*M_PI* eos["0"] * eos[(int)EOS_Idx::n]).pow(2./3.)).sqrt() / Zaki::Physics::MEV_2_INV_FM ;
  
  fermi_electron.label = "FE_e [MeV]" ;

  Zaki::Vector::DataColumn fermi_muon =  ( 
    pow(Zaki::Physics::MUON_M_FM, 2) 
    + (3*M_PI*M_PI* eos["1"] * eos[(int)EOS_Idx::n]).pow(2./3.)).sqrt() / Zaki::Physics::MEV_2_INV_FM ;

  fermi_muon.label = "FE_mu [MeV]" ;

  Zaki::Vector::DataColumn fermi_neutron =  ( 
    m_eff["10"].pow(2) 
    + (3*M_PI*M_PI* eos["10"] * eos[(int)EOS_Idx::n]).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
    ).sqrt() + V_eff["10"]  ;
  fermi_neutron.label = "FE_n [MeV]" ;

  Zaki::Vector::DataColumn fermi_lambda =  ( 
    m_eff["100"].pow(2) 
    + (3*M_PI*M_PI* eos["100"] * eos[(int)EOS_Idx::n]).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
    ).sqrt() + V_eff["100"]  ;
  fermi_lambda.label = "FE_Lam [MeV]" ;

  Zaki::Vector::DataSet fermi_ds({eos[(int)EOS_Idx::n], 
                                  fermi_electron, fermi_muon, 
                                  fermi_neutron, fermi_lambda}) ;

  Zaki::String::Directory file_dir = in_dir ;
  if (in_dir_type == Dir_Type::relative)
  {
    file_dir = wrk_dir + "/" + in_dir ;
  }

  fermi_ds.Export(file_dir + "/" + name + "_E_Fermi.tsv", {},'\t') ;
}

//--------------------------------------------------------------
/// Returns the Fermi energy of particles as a function of density 
Zaki::Vector::DataSet CompactStar::CompOSE_EOS::GetFermiE() const 
{
  // ------------------------------------
  //        Finding Fermi Energy
  // ------------------------------------
  Zaki::Vector::DataColumn fermi_electron =  ( 
    pow(Zaki::Physics::ELECTRON_M_FM, 2) 
    + (3*M_PI*M_PI* eos["0"] * eos[(int)EOS_Idx::n]).pow(2./3.)).sqrt() / Zaki::Physics::MEV_2_INV_FM ;
  
  fermi_electron.label = "FE_e [MeV]" ;

  Zaki::Vector::DataColumn fermi_muon =  ( 
    pow(Zaki::Physics::MUON_M_FM, 2) 
    + (3*M_PI*M_PI* eos["1"] * eos[(int)EOS_Idx::n]).pow(2./3.)).sqrt() / Zaki::Physics::MEV_2_INV_FM ;

  fermi_muon.label = "FE_mu [MeV]" ;

  Zaki::Vector::DataColumn fermi_neutron =  ( 
    m_eff["10"].pow(2) 
    + (3*M_PI*M_PI* eos["10"] * eos[(int)EOS_Idx::n]).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
    ).sqrt() + V_eff["10"]  ;
  fermi_neutron.label = "FE_n [MeV]" ;

  Zaki::Vector::DataColumn fermi_lambda =  ( 
    m_eff["100"].pow(2) 
    + (3*M_PI*M_PI* eos["100"] * eos[(int)EOS_Idx::n]).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
    ).sqrt() + V_eff["100"]  ;
  fermi_lambda.label = "FE_Lam [MeV]" ;

  Zaki::Vector::DataSet fermi_ds({eos[(int)EOS_Idx::n], 
                                  fermi_electron, fermi_muon, 
                                  fermi_neutron, fermi_lambda}) ;

  return fermi_ds ;
}
//--------------------------------------------------------------

//==============================================================
