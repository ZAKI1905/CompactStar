/*
  BNV_Chi class
*/

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include <CompactStar/BNV/BNV_Sequence.hpp>
#include <CompactStar/Core/TOVSolver.hpp>

using namespace CompactStar ;

// ------------------------------------------------------------
// Default Constructor
BNV_Sequence::BNV_Sequence()
{

  SetInitOmega() ;

  // if (Pdot_sign < 0 )
  // {
  //   Pdot_sign_str = "-" ;
  // }
} 

// ------------------------------------------------------------
// Constructor from a sequence file
BNV_Sequence::BNV_Sequence( const Zaki::String::Directory& in_dir, 
                          const std::string& in_f_name)
{
  SetWrkDir(in_dir) ;
  seq.Import(in_dir + "/" + in_f_name) ;
  seq.SetWrkDir(in_dir) ;

  // Making the sequence smoother
  seq.MakeSmooth(35) ;

  // Removing unstable NSs (i.e., after M_max)
  seq.Trim(seq[1].MaxIdx(), seq[1].Size()) ;

  // Interpolating mass, radius, baryon number, and MomI as a func of eps_c
  seq.Interpolate(0, {1, 2, 4, 5}) ;
}

// ------------------------------------------------------------
// Destructor
BNV_Sequence::~BNV_Sequence()
{}
// ------------------------------------------------------------
// // Copy constructor
// BNV_Sequence::BNV_Sequence(const BNV_Sequence& other) 
// {
//   Z_LOG_NOTE("BNV_Sequence copy constructor called: from " + other.PtrStr() + " --> " + PtrStr()) ;
// }

// // ------------------------------------------------------------
// // Assignment operator
// BNV_Sequence& BNV_Sequence::operator=(const BNV_Sequence& other) 
// {
//   Z_LOG_NOTE("BNV_Sequence '=' operator called: from " + other.PtrStr() + " --> " + PtrStr()) ;
// }

// ------------------------------------------------------------
// Set the initial omega set
void BNV_Sequence::SetInitOmega() 
{
  if(!init_omega_is_set)
  {
    // ---------------------------------------------
    //  Set the initial index based on beta dataset
    // ---------------------------------------------
    //      H = 1e8
    // -------------------
    if (mag_field == 1e8)
    {
      init_omega_set = {
                        // {2*M_PI/1e-3, "1 (ms)"}, // Br_idx
                        // {2*M_PI/2e-3, "2 (ms)"}, // Br_idx
                        // {2*M_PI/3e-3, "3 (ms)"},
                        // {2*M_PI/5e-3, "5 (ms)"},
                        // {2*M_PI*1e2, "10 (ms)"}, // Br_idx
                        // {2*M_PI/20e-3, "20 (ms)"},
                        // {2*M_PI/0.1, "100 (ms)"},
                        // {2*M_PI/1, "1 (s)"},
                        // {2*M_PI/5, "5 (s)"},
                        // {2*M_PI/10, "10 (s)"}
                        {2*M_PI*1e3, "1 (ms)"}, // m_omega evol
                        {2*M_PI*1e2, "10 (ms)"}, // m_omega evol
                        {2*M_PI*1e1, "100 (ms)"}, // m_omega evol
                        {2*M_PI*1e0, "1 (s)"}, // m_omega evol
                        // {200, "200 (Hz)"}, 
                        // {180, "180 (Hz)"}, 
                        // {150, "150 (Hz)"}, 
                        // {100, "100 (Hz)"}, 
                        // {50, "50 (Hz)"}, 
                        // {10, "10 (Hz)"} 
                        } ;
    } 
    // -------------------
    //      H = 1e9
    // -------------------
    else if (mag_field == 1e9)
    {
      init_omega_set = {
                        // {2*M_PI/1e-3, "1 (ms)"},
                        // {2*M_PI/2e-3, "2 (ms)"}, 
                        // {2*M_PI/5e-3, "5 (ms)"}, 
                        // {2*M_PI/10e-3, "10 (ms)"}

                        {2*M_PI*1e3, "1 (ms)"}, // m_omega evol
                        {2*M_PI*1e2, "10 (ms)"}, // m_omega evol
                        {2*M_PI*1e1, "100 (ms)"}, // m_omega evol
                        {2*M_PI*1e0, "1 (s)"}, // m_omega evol
                        } ;
    }
    // -------------------
    //      H = 1e10
    // -------------------
    else if (mag_field == 1e10)
    {
      // init_omega_set = {{2*M_PI/1e-1, "0.1 (s)"}, 
      //                   {2*M_PI/5e-1, "0.5 (s)"}, 
      //                   {2*M_PI/1, "1 (s)"}} ;
      // init_omega_set = {{2*M_PI/10e-3, "10 (ms)"}, 
      //                   {2*M_PI/0.2, "0.2 (s)"}, 
      //                   {2*M_PI/1, "1 (s)"}} ;
      init_omega_set = {
                        // {2*M_PI/10e-3, "10 (ms)"},
                        // {2*M_PI/50e-3, "50 (ms)"}, 

                        // {2*M_PI*1e3, "1 (ms)"}, // * schematic 
                        // {2*M_PI/32e-3, "32 (ms)"}, // * schematic
                        // {2*M_PI/200e-3, "200 (ms)"} // * schematic

                        {2*M_PI*1e3, "1 (ms)"}, // m_omega evol
                        {2*M_PI*1e2, "10 (ms)"}, // m_omega evol
                        {2*M_PI*1e1, "100 (ms)"}, // m_omega evol
                        {2*M_PI*1e0, "1 (s)"}, // m_omega evol

                        // {2*M_PI/0.2, "0.2 (s)"}, 
                        // {2*M_PI/0.3, "0.3 (s)"}, 
                        // {2*M_PI/0.4, "0.4 (s)"}, 
                        // {2*M_PI/1, "1 (s)"}
                        } ;
    }
    // -------------------
    //      H = 1e11
    // -------------------
    else if (mag_field == 1e11)
    {
      // init_omega_set = {{2*M_PI/1, "1 (s)"}, 
      //                   {2*M_PI/5, "5 (s)"}} ;
      init_omega_set = {
                        // {2*M_PI/10e-3, "10 (ms)"},
                        // {2*M_PI/50e-3, "50 (ms)"}, 
                        // {2*M_PI/0.2, "0.2 (s)"}, 
                        // {2*M_PI/1, "1 (s)"}
                        
                        {2*M_PI*1e3, "1 (ms)"}, // m_omega evol
                        {2*M_PI*1e2, "10 (ms)"}, // m_omega evol
                        {2*M_PI*1e1, "100 (ms)"}, // m_omega evol
                        {2*M_PI*1e0, "1 (s)"}, // m_omega evol
                        } ;
    }
    // -------------------
    else
    {
      Z_LOG_ERROR("Set init_omega_set!") ;
      init_omega_set = {} ;
    }
    // -----------------------------------------
  }
}

// ------------------------------------------------------------
void BNV_Sequence::SetInitOmegaSet(const std::vector<std::pair<double, std::string>>& init_omegas) 
{
  init_omega_set = init_omegas ;
  init_omega_is_set = true ;
}

// ------------------------------------------------------------
/** Generates a sequence of neutron stars.
*
*
* @param in_dir directory
* @param in_model model number
*/
void BNV_Sequence::GenSequence(const Zaki::String::Directory& in_dir, 
                              const std::string& in_model)
{
  // .......................................................
  //                     Solving TOV
  // .......................................................
  CompactStar::TOVSolver solver ;

  solver.SetWrkDir(in_dir.ParentDir().ParentDir()) ;
  solver.ImportEOS("EOS/CompOSE/" + in_model + "/"+ in_model +".eos") ;

  // solver.Solve( {{1.0e+14, 1.8e+15}, 800, "Linear"}, 
  //               "SpinDown_2023/results/DS(CMF)-1_with_crust_Seq",  
  //               "SD_HiRes") ;      

  // solver.AddNCondition(MassCondition) ;
  solver.SetMaxRadius(40) ;
  // solver.SetMaxRadius(20) ;

  solver.SetRadialRes(90000) ;
  
  // solver.Solve( {{1.6085e+14, 1.8e+15}, 1000, "Linear"}, 
  //               "SpinDown_2023/results/"+in_model+"/Nov_2023",  
  //               "SD_HiRes_Nov_2023") ;   

  // solver.Solve( {{1.8e+15, 2.15e+15}, 214, "Linear"}, 
  //               "SpinDown_2023/results/"+in_model+"/Nov_2023",  
  //               "SD_HiRes_Nov_2023_high") ;  
  
  // solver.Solve( {{1.385e+14, 1.6085e+14}, 30, "Linear"}, 
  //               "SpinDown_2023/results/"+in_model+"/Nov_2023",  
  //               "SD_HiRes_Nov_2023_low") ;  

  // solver.Solve( {{2.986e+14, 2.15e+15}, 2000, "Linear"}, 
  //               "SpinDown_2023/results/"+in_model+"/Nov_2023",  
  //               "SD_HiRes_Nov_2023") ;   

  double tmp_step = 9.257e+11 ;
  // solver.Solve( {{1.541908e+14, 2.986e+14-tmp_step}, 155, "Linear"}, 
  //               "SpinDown_2023/results/"+in_model+"/Nov_2023",  
  //               "SD_HiRes_Nov_2023_med") ;   

  solver.Solve( {{1.393796e+14, 1.541908e+14-tmp_step}, 15, "Linear"}, 
                "SpinDown_2023/results/"+in_model+"/Nov_2023",  
                "SD_HiRes_Nov_2023_low") ;   
  // solver.Solve( {{1.34e+14, 2.15e+15}, 500, "Log"}, 
                // "BNV_2023/results/"+ in_model +"/BNV_Seq",  
                // "BNV_Seq") ;   
  
  // solver.Solve( {{1.337e+15, 1.339e+15}, 10, "Linear"}, 
  // // solver.Solve( {{1.31e+15, 1.36e+15}, 10, "Linear"}, 
  //               "BNV_2023/results/"+ in_model +"/BNV_Seq",  
  //               "BNV_Seq_zoom") ;          

  // solver.Solve( {{1.247e+15, 1.427e+15}, 6, "Linear"}, 
  //               "BNV_2023/results/"+ in_model +"/BNV_Seq",  
  //               "BNV_Seq") ;      
                // 1.346
}

// ------------------------------------------------------------
// Imports the sequence
void BNV_Sequence::Import(const Zaki::String::Directory& in_dir, 
                          const std::string& in_f_name)
{
  SetWrkDir(in_dir) ;
  seq.SetWrkDir(in_dir) ;
  seq.Import(in_f_name) ;

  // seq.RemoveBumps() ;

  // Making the sequence smoother
  seq.MakeSmooth(50) ;

  // Removing unstable NSs (i.e., after M_max)
  seq.Trim(seq[1].MaxIdx(), seq[1].Size()) ;

  // Find_b_M(2.01) ;

  // Interpolating mass, radius, baryon number, and MomI as a func of eps_c
  seq.Interpolate(0, {1, 2, 4, 5}) ;
}

// ------------------------------------------------------------
// Evaluates b_O factor using the central difference method
//  and given the input NS mass
void BNV_Sequence::Find_b_O(const double& in_mass, const int& o_idx)
{
  // Find the closest index 
  int m_idx = seq[1].GetClosestIdx(in_mass) ;

  // The step size
  double h = 2*(seq[0][m_idx+1] - seq[0][m_idx]) ;

  std::cout << "\n    * ================================= * ";
  std::cout << "\n\t Half-step size = " << h/2 << " (g/cm^3) " ;

  // Evaluate the derivative at m_idx 
  double o_der = seq[o_idx][m_idx+1] - seq[o_idx][m_idx-1] ;
        o_der /= h ;

  // Evaluate the 2nd derivative at m_idx 
  double o_2nd_der = seq[o_idx][m_idx+2] + seq[o_idx][m_idx-2] - 2*seq[o_idx][m_idx] ;
        o_2nd_der /= h*h ;

  // Evaluate the 3rd derivative at m_idx 
  double o_3rd_der = seq[o_idx][m_idx+3] - 3*seq[o_idx][m_idx+1] 
                    + 3*seq[o_idx][m_idx-1] - seq[o_idx][m_idx-3] ;
        o_3rd_der /= h*h*h ;

  double o_der_error = - o_3rd_der * h * h / 24.0 ;
  std::cout << "\n    * ================================= * ";
  std::cout << "\n\t O   = " << seq[o_idx][m_idx] << " " ;
  std::cout << "\n\t O'  = " << o_der << " " ;
  std::cout << "\n\t O'' = " << o_2nd_der << " " ;
  std::cout << "\n\t O'''= " << o_3rd_der << " " ;
  std::cout << "\n\t Error ( O'[i] ) = " << o_der_error << " " ;
  // std::cout << "\n    * ================================= * ";

  // Evaluate the B derivative at m_idx 
  double B_der = seq[4][m_idx+1] - seq[4][m_idx-1] ;
        B_der /= h ;

  // Evaluate the 2nd derivative at m_idx 
  double B_2nd_der = seq[4][m_idx+2] + seq[4][m_idx-2] - 2*seq[4][m_idx] ;
        B_2nd_der /= h*h ;

  // Evaluate the 3rd derivative at m_idx 
  double B_3rd_der = seq[4][m_idx+3] - 3*seq[4][m_idx+1] 
                    + 3*seq[4][m_idx-1] - seq[4][m_idx-3] ;
        B_3rd_der /= h*h*h ;

  double B_der_error = - B_3rd_der * h * h / 24.0 ;

  std::cout << "\n    * ================================= * ";
  std::cout << "\n\t B'  = " << B_der << " " ;
  std::cout << "\n\t B'' = " << B_2nd_der << " " ;
  std::cout << "\n\t B'''= " << B_3rd_der << " " ;
  std::cout << "\n\t Error ( B'[i] ) = " << B_der_error << " " ;
  // std::cout << "\n    * ================================= * ";

  double b_O = seq[4][m_idx] * o_der / ( B_der * seq[o_idx][m_idx]  ) ;
  double b_O_error = b_O * sqrt( pow(o_der_error/o_der,2) 
                                  + pow(B_der_error/B_der, 2)) ;
  std::cout << "\n    * ================================= * ";
  std::cout << "\n\t b_O = " << b_O << " +-  " << b_O_error ;
  std::cout << "\n    * ================================= * ";
}

// ------------------------------------------------------------
// Evaluates b_O factor using the central difference method
void BNV_Sequence::Find_b_factors()
{
  // Find the closest index 
  // int m_idx = seq[1].GetClosestIdx(in_mass) ;

  Zaki::Vector::DataSet b_o_ds(9, seq[0].Size()) ;

  for (size_t i = 2; i < seq[0].Size()-2; i++)
  {
    // The step size
    double h = 2*(seq[0][i+1] - seq[0][i]) ;

    // std::cout << "\n    * ================================= * ";
    // std::cout << "\n\t Half-step size = " << h/2 << " (g/cm^3) " ;

    // Evaluate the derivative at i 
    double R_der = seq[2][i+1] - seq[2][i-1] ;
          R_der /= h ;
    double I_der = seq[5][i+1] - seq[5][i-1] ;
          I_der /= h ;

    double M_der = seq[1][i+1] - seq[1][i-1] ;
          M_der /= h ;

    // Evaluate the 2nd derivative at i 
    double I_2nd_der = seq[5][i+2] + seq[5][i-2] - 2*seq[5][i] ;
          I_2nd_der /= h*h ;

    // Evaluate the 3rd derivative at i 
    // double o_3rd_der = seq[o_idx][i+3] - 3*seq[o_idx][i+1] 
    //                   + 3*seq[o_idx][i-1] - seq[o_idx][i-3] ;
    //       o_3rd_der /= h*h*h ;

    // double o_der_error = - o_3rd_der * h * h / 24.0 ;

    // std::cout << "\n    * ================================= * ";
    // std::cout << "\n\t O   = " << seq[o_idx][i] << " " ;
    // std::cout << "\n\t O'  = " << o_der << " " ;
    // std::cout << "\n\t O'' = " << o_2nd_der << " " ;
    // std::cout << "\n\t O'''= " << o_3rd_der << " " ;
    // std::cout << "\n\t Error ( O'[i] ) = " << o_der_error << " " ;
    // std::cout << "\n    * ================================= * ";

    // Evaluate the B derivative at i 
    double B_der = seq[4][i+1] - seq[4][i-1] ;
          B_der /= h ;

    // Evaluate the 2nd derivative at i 
    double B_2nd_der = seq[4][i+2] + seq[4][i-2] - 2*seq[4][i] ;
          B_2nd_der /= h*h ;

    // // Evaluate the 3rd derivative at i 
    // double B_3rd_der = seq[4][i+3] - 3*seq[4][i+1] 
    //                   + 3*seq[4][i-1] - seq[4][i-3] ;
    //       B_3rd_der /= h*h*h ;

    // double B_der_error = - B_3rd_der * h * h / 24.0 ;

    // std::cout << "\n    * ================================= * ";
    // std::cout << "\n\t B'  = " << B_der << " " ;
    // std::cout << "\n\t B'' = " << B_2nd_der << " " ;
    // std::cout << "\n\t B'''= " << B_3rd_der << " " ;
    // std::cout << "\n\t Error ( B'[i] ) = " << B_der_error << " " ;
    // std::cout << "\n    * ================================= * ";

    double b_R = seq[4][i] * R_der / ( B_der * seq[2][i]  ) ;
    double b_I = seq[4][i] * I_der / ( B_der * seq[5][i]  ) ;
    double b_M = seq[4][i] * M_der / ( B_der * seq[1][i]  ) ;

    double b_beta_I = I_2nd_der / (I_der * B_der) 
                       - B_2nd_der / pow(B_der, 2) ;
           b_beta_I *=  seq[4][i] ;

    // double b_O_error = b_O * sqrt( pow(o_der_error/o_der,2) 
    //                                 + pow(B_der_error/B_der, 2)) ;
    // std::cout << "\n    * ================================= * ";
    // std::cout << "\n\t b_O = " << b_O << " +-  " << b_O_error ;
    // std::cout << "\n    * ================================= * ";

    b_o_ds[0].vals.emplace_back(seq[0][i]) ; // 0-eps
    b_o_ds[1].vals.emplace_back(seq[1][i]) ; // 1-M
    b_o_ds[2].vals.emplace_back(seq[2][i]) ; // 2-R
    b_o_ds[3].vals.emplace_back(seq[4][i]) ; // 3-B
    b_o_ds[4].vals.emplace_back(seq[5][i]) ; // 4-I

    b_o_ds[5].vals.emplace_back(b_R) ; // 5-b_R
    b_o_ds[6].vals.emplace_back(b_I) ; // 6-b_I
    b_o_ds[7].vals.emplace_back(b_beta_I) ; // 7-b(beta(I))
    b_o_ds[8].vals.emplace_back(b_M) ; // 8-b_M
  }

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetGrid(true) ;
  plt_par.SetYAxisLabel("$b (O)$") ;
  plt_par.SetYAxis({-100, 3}) ;
  plt_par.SetLegend({"lower right", 1.0, 0.0}) ;
  
  // Adding error bars
  // b_o_ds[3] = b_o_ds[1] + b_o_ds[2] ;
  // b_o_ds[4] = b_o_ds[1] - b_o_ds[2] ;

  b_o_ds.SetPlotPars(plt_par) ;
  b_o_ds.SetWrkDir(GetWrkDir()) ;
  b_o_ds.Plot(1, {{5, "$b (R)$"}, 
                    {6, "$b (I)$"},
                    {7, "$b (\\beta(I))$"}
                    },
                     "b_O_vs_eps_CDM.pdf", "") ;

  b_o_ds[0].label = "ec(g/cm^3)" ;
  b_o_ds[1].label = "M" ;
  b_o_ds[2].label = "R(km)" ;
  b_o_ds[3].label = "B" ;
  b_o_ds[4].label = "I(km^3)" ;
  b_o_ds[5].label = "b(R)" ;
  b_o_ds[6].label = "b(I)" ;
  b_o_ds[7].label = "b(beta(I))" ;
  b_o_ds[8].label = "b(M)" ;

  // b_o_ds.MakeSmooth(10) ;
  b_o_ds.Export("B_Factors.tsv") ;
}

// ------------------------------------------------------------
// Returns the observable with the given index
// Zaki::Vector::DataColumn BNV_Sequence::O(const int& idx) const
// {
//   return seq[idx] ;
// }

// ------------------------------------------------------------
// Returns the central energy density
Zaki::Vector::DataColumn BNV_Sequence::Eps_C() const 
{
  return beta[beta_idx.eps] ;
}

// ------------------------------------------------------------
// Returns Mass
Zaki::Vector::DataColumn BNV_Sequence::M() const
{
  return beta[beta_idx.M] ;
}

// ------------------------------------------------------------
// Returns quantity O as a function of time
double BNV_Sequence::O(const double& time, const int& o_idx) const  
{
  // double B_val = B(time) ;
  // int i = beta[beta_idx.B].GetClosestIdx( B_val ) ;

  // if(B_val > beta[beta_idx.B][i])
  // {
  //   i++ ;
  // }
  
  // double o = (B_val - beta[beta_idx.B][i-1]) * beta[o_idx][i-1] ;
  //       o += (beta[beta_idx.B][i] - B_val) * beta[o_idx][i] ;
  //       o /= beta[beta_idx.B][i] - beta[beta_idx.B][i-1] ;
  
  weighted_idx w_idx = Time_to_Weighted_Idx(time) ;

  double o = w_idx.w_L * beta[o_idx][w_idx.idx-1] 
            + w_idx.w_R * beta[o_idx][w_idx.idx] ;

  return o ;
}

// ------------------------------------------------------------
// Returns mass as a function of time
// [Edited on Sep 29, 2023] : The issue was that 
// closest index generates a step-like behaviour 
// and is not accurate since it can return the same mass
// for two different values of B!
// We now interpolate to avoid this issue.
// We also updated other similar functions. 
double BNV_Sequence::M(const double& time) const  
{
  // return O(time, beta_idx.M);

  // double B_val = B(time) ;
  // int i = beta[beta_idx.B].GetClosestIdx( B_val ) ;
  // return M()[i] ; // pre Sep 29, 2023 version!

  weighted_idx w_idx = Time_to_Weighted_Idx(time) ;

  return  w_idx.w_L * M()[w_idx.idx-1] 
            + w_idx.w_R * M()[w_idx.idx] ;

  // if(B_val > beta[beta_idx.B][i])
  // {
  //   i++ ;
  // }
  
  // std::cout << " [ i = " << i << " ] " ; 

  // double m = (B_val - beta[beta_idx.B][i-1]) * beta[beta_idx.M][i-1] ;
  //       m += (beta[beta_idx.B][i] - B_val) * beta[beta_idx.M][i] ;
  //       m /= beta[beta_idx.B][i] - beta[beta_idx.B][i-1] ;

  // return m ;
}

// ------------------------------------------------------------
// Returns Mass
Zaki::Vector::DataColumn BNV_Sequence::R() const
{
  return beta[beta_idx.R] ;
}

// ------------------------------------------------------------
// Returns radius as a function of time
double BNV_Sequence::R(const double& time) const 
{
  // return O(time, beta_idx.R);

  // int idx = B().GetClosestIdx( B(time) ) ;

  // return R()[Time_to_Idx(time)] ;
  weighted_idx w_idx = Time_to_Weighted_Idx(time) ;

  return  w_idx.w_L * R()[w_idx.idx-1] 
            + w_idx.w_R * R()[w_idx.idx] ;
}

// ------------------------------------------------------------
// Returns Mass
Zaki::Vector::DataColumn BNV_Sequence::I() const
{
  return beta[beta_idx.I] ;
}

// ------------------------------------------------------------
// Returns moment-of-inertia as a function of time
double BNV_Sequence::I(const double& time) const 
{
  // return O(time, beta_idx.I);

  // int idx = B().GetClosestIdx( B(time) ) ;

  // return I()[Time_to_Idx(time)] ;
  weighted_idx w_idx = Time_to_Weighted_Idx(time) ;

  return  w_idx.w_L * I()[w_idx.idx-1] 
            + w_idx.w_R * I()[w_idx.idx] ;
}

// ------------------------------------------------------------
// Returns baryon number
Zaki::Vector::DataColumn BNV_Sequence::B() const 
{
  return beta[beta_idx.B] ;
}

// ------------------------------------------------------------
// Returns baryon number as a function of time
// !!! Assumes that t_0 << time
double BNV_Sequence::B(const double& time) const 
{
  // Adding a cut-off for BNV 
  if ( time > bnv_cutoff_time ) 
    return B()[init_idx] * exp(-gamma_bnv * bnv_cutoff_time) ;

  return B()[init_idx] * exp(-gamma_bnv * time) ;
}

// ------------------------------------------------------------
// Records the time when BNV turns off
void BNV_Sequence::SetBNVCuttoff(const double& time) 
{
  if ( bnv_active_flag && M(time) < bnv_cuttoff_mass )
  {
    bnv_cutoff_time = time ;
    bnv_active_flag = false ;
  }
}

// ------------------------------------------------------------
void BNV_Sequence::Plot_Dimless_O() const
{
  Zaki::Vector::DataSet Scaled_O({seq[0],       // 0-eps
                                  seq[1]/1.4,   // 1-(M/1.4)
                                  seq[2]/12,    // 2-(R/12km)
                                  seq[4]/1e57,  // 3-(B/2e57)
                                  seq[5]/ (70*Zaki::Physics::SUN_M_KM)  
                                    // 4-(I/70 M_Sun * km^2) 
                                  } ) ;
  // Scaled_O[2].label = "$R / 10 km$" ;
  // Scaled_O[3].label = "$B / 2\\times 10^{57} $" ;
  // Scaled_O[4].label = "$I / 120 km^3$" ;

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetGrid(true) ;
  plt_par.SetYAxisLabel("$\\frac{O}{O^*}$") ;
  plt_par.SetYAxis({0, 3}) ;
  plt_par.SetLegend({"lower right", 1.0, 0.0}) ;

  // Scaled_O.MakeSmooth(10) ;
  Scaled_O.SetWrkDir( GetWrkDir() ) ;

  plt_par.SetXAxis({2e+14, 2e+15}) ;
  Scaled_O.SetPlotPars(plt_par) ;
  Scaled_O.SemiLogXPlot(0, {{1, "$M\\,/\\,1.4\\, M_{s}$"}, 
                    {2, "$R\\, /\\, 12\\, km$"},
                    {3, "$B\\, /\\, 10^{57}$"},
                    {4, "$I\\, /\\, 70 \\, M_{s} \\, km^2$"}},
                     "O_vs_ec.pdf", "") ;

  plt_par.SetXAxis({0.1, 2.1}) ;
  Scaled_O.SetPlotPars(plt_par) ;
  Scaled_O[1] *= 1.4 ;
  Scaled_O[1].label = "M" ;
  Scaled_O.Plot(1, {{2, "$R\\, /\\, 12\\, km$"}, 
                    {3, "$B\\, /\\, 10^{57}$"},
                    {4, "$I \\,/\\, 70 \\, M_{s} \\, km^2$"}},
                     "O_vs_M.pdf", "") ;
}

// ------------------------------------------------------------
// Evaluates the beta factor for the given index
void BNV_Sequence::EvalBeta() 
{

  // seq.MakeSmooth(5) ;

  // b_o_ds[0].label = "ec(g/cm^3)" ;
  // b_o_ds[1].label = "M" ;
  // b_o_ds[2].label = "R(km)" ;
  // b_o_ds[3].label = "B" ;
  // b_o_ds[4].label = "I(km^3)" ;
  // b_o_ds[5].label = "b(R)" ;
  // b_o_ds[6].label = "b(I)" ;
  // b_o_ds[7].label = "b(beta(I))" ;

  beta.SetWrkDir(GetWrkDir()) ;
  beta.Import("B_Factors.tsv") ;
  beta.Interpolate(0, beta_idx.B) ;

  // D_eps_B :
  // beta = {seq[0],               // 0-eps
  //         seq[1],               // 1-M
  //         seq[2],               // 2-R
  //         seq[4],               // 3-B
  //         seq[5],               // 4-I
  //         seq.Derivative(1)[1] / seq.Derivative(4)[1],  // 5-beta(M)
  //         seq.Derivative(2)[1] / seq.Derivative(4)[1],  // 6-beta(R)
  //         seq.Derivative(5)[1] / seq.Derivative(4)[1],   // 7-beta(I)
  //         seq[4] / 1.4e57 
  //         } ;
  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetGrid(true) ;
  plt_par.SetLegend({"lower left", 0.0, 0.0}) ;
  plt_par.SetXAxis({0.5, 2.06}) ;
  plt_par.SetYAxis({5e-3, 1.5}) ;
  // plt_par.SetYAxis({-1.5, 1.5}) ;

  Zaki::Vector::DataSet tmp_ds({beta[1], beta[8], beta[5], beta[6], 
                                beta[6]*beta[7], 
                                beta[5] + 6*beta[6], 
                                beta[5] * ( 1 + beta[7] - 6*beta[6])}) ;

  tmp_ds.SetWrkDir(GetWrkDir()) ;
  tmp_ds.SetPlotPars(plt_par) ;
  tmp_ds.SemiLogYPlot(0, {{1, "$b(M)$"}, {2, "$b(R)$"}, {3, "$b(I)$"}, 
                  {4, "$b(I) b(\\beta(I))$"},  {5, "$b(I) + 6 b(R)$"},
                  {6, "$b(I) [1 + b(\\beta(I)) - 6 b(R) ]$"}}, "B_Factors.pdf", "") ;


  // beta.Interpolate(0, beta_idx.B) ;
  // Zaki::Vector::DataSet::PlotParam plt_par ;
  // plt_par.SetGrid(true) ;
  // plt_par.SetYAxisLabel("$\\beta(O)$") ;
  // plt_par.SetLegend({"lower left", 0.0, 0.0}) ;
  // beta.SetPlotPars(plt_par) ;
  // // beta.MakeSmooth(20) ; // Change back to 50 for a smoother result!!
  // beta.SetWrkDir( GetWrkDir() ) ;
  // beta.Plot(1, {{5, "$\\beta(M)$"}, {6, "$\\beta(R)$"}, {7, "$\\beta(I)$"}}, "Beta_O.pdf", "") ;
  // plt_par.SetYAxisLabel("$\\beta(M)$") ;
  // beta.SetPlotPars(plt_par) ;
  // beta.Plot(1, {{5, "$\\beta(M)$"}}, "Beta_M.pdf", "") ;
  // plt_par.SetYAxisLabel("$\\mathcal{O}$") ;
  // beta.SetPlotPars(plt_par) ;
  // beta.Plot(0, {{1, "$M$"}, {8, "$B/1.4\\times 10^{57}$"}}, "M_vs_ec.pdf", "") ;

  // Zaki::Vector::DataSet tmp_ds({seq[0], 5e15 * seq.Derivative(1)[1], 5e15*seq.Derivative(4)[1]/1.4e+57,
  //                                 seq[4] * seq.Derivative(1)[1] / (0.8 * seq.Derivative(4)[1] * seq[1]) }) ;
  // tmp_ds.SetWrkDir(GetWrkDir()) ;
  // plt_par.SetYAxisLabel("$d\\mathcal{O} / d\\varepsilon$") ;
  // tmp_ds.SetPlotPars(plt_par) ;
  // tmp_ds.Plot(0, {{1, "$M'$"}, {2, "$B'/1.4\\times 10^{57}$"}, {3, "$b_M$"}}, "M'_vs_ec.pdf", "") ;

} 
// ------------------------------------------------------------
// Input in yr^-1
void BNV_Sequence::SetBNVRate(const double& in_gamma) 
{
  gamma_bnv = in_gamma / Zaki::Physics::YR_2_SEC ;
}

// ------------------------------------------------------------
// Sets the mass at which BNV shuts down
void BNV_Sequence::SetBNVCuttOffMass(const double& in_M) 
{
  bnv_cuttoff_mass = in_M ;
}

// ------------------------------------------------------------
// Returns Beta_I
// Zaki::Vector::DataColumn BNV_Sequence::Beta_I() const 
// {
//   return beta[beta_idx.beta_I] ;
// }

// ------------------------------------------------------------
// Returns Beta_R
// Zaki::Vector::DataColumn BNV_Sequence::Beta_R() const 
// {
//   return beta[beta_idx.beta_R] ;
// }

// ------------------------------------------------------------
// Returns Beta_M
// Zaki::Vector::DataColumn BNV_Sequence::Beta_M() const 
// {
//   return beta[beta_idx.beta_M] ;
// }

// ------------------------------------------------------------
// Returns b_I
Zaki::Vector::DataColumn BNV_Sequence::b_I() const 
{
  // return beta[beta_idx.beta_I] * beta[beta_idx.B] / beta[beta_idx.I] ;
  return beta[beta_idx.b_I] ;
}

// ------------------------------------------------------------
// Returns b_I as a function of time
double BNV_Sequence::b_I(const double& time) const 
{
  // int idx = B().GetClosestIdx( B(time) ) ;

  // return b_I()[Time_to_Idx(time)] ;

  weighted_idx w_idx = Time_to_Weighted_Idx(time) ;

  return  w_idx.w_L * b_I()[w_idx.idx-1] 
            + w_idx.w_R * b_I()[w_idx.idx] ;

  // double B_val = B(time) ;
  // int i = beta[beta_idx.B].GetClosestIdx( B_val ) ;

  // if(B_val > beta[beta_idx.B][i])
  // {
  //   i++ ;
  // }
  
  // double b_i = (B_val - beta[beta_idx.B][i-1]) * b_I()[i-1] ;
  //       b_i += (beta[beta_idx.B][i] - B_val) * b_I()[i] ;
  //       b_i /= beta[beta_idx.B][i] - beta[beta_idx.B][i-1] ;
  
  // return b_i ;
}

// ------------------------------------------------------------
// Returns b(beta_I)
Zaki::Vector::DataColumn BNV_Sequence::b_beta_I() const 
{
  // return beta[beta_idx.beta_I] * beta[beta_idx.B] / beta[beta_idx.I] ;
  return beta[beta_idx.b_beta_I] ;
}

// ------------------------------------------------------------
// Returns b(beta_I) as a function of time
double BNV_Sequence::b_beta_I(const double& time) const 
{
  // int idx = B().GetClosestIdx( B(time) ) ;

  // return b_beta_I()[Time_to_Idx(time)] ;

  weighted_idx w_idx = Time_to_Weighted_Idx(time) ;

  return  w_idx.w_L * b_beta_I()[w_idx.idx-1] 
            + w_idx.w_R * b_beta_I()[w_idx.idx] ;

  // double B_val = B(time) ;
  // int i = beta[beta_idx.B].GetClosestIdx( B_val ) ;

  // if(B_val > beta[beta_idx.B][i])
  // {
  //   i++ ;
  // }
  
  // double b_i = (B_val - beta[beta_idx.B][i-1]) * b_I()[i-1] ;
  //       b_i += (beta[beta_idx.B][i] - B_val) * b_I()[i] ;
  //       b_i /= beta[beta_idx.B][i] - beta[beta_idx.B][i-1] ;
  
  // return b_i ;
}

// ------------------------------------------------------------
// Returns b_R
Zaki::Vector::DataColumn BNV_Sequence::b_R() const 
{
  return beta[beta_idx.b_R] ;
}

// ------------------------------------------------------------
// Returns b_R as a function of time
double BNV_Sequence::b_R(const double& time) const 
{
  // int idx = B().GetClosestIdx( B(time) ) ;

  // return b_R()[Time_to_Idx(time)] ;

  weighted_idx w_idx = Time_to_Weighted_Idx(time) ;

  return  w_idx.w_L * b_R()[w_idx.idx-1] 
            + w_idx.w_R * b_R()[w_idx.idx] ;

  // double B_val = B(time) ;
  // int i = beta[beta_idx.B].GetClosestIdx( B_val ) ;

  // if(B_val > beta[beta_idx.B][i])
  // {
  //   i++ ;
  // }
  
  // double b_r = (B_val - beta[beta_idx.B][i-1]) * b_R()[i-1] ;
  //       b_r += (beta[beta_idx.B][i] - B_val) * b_R()[i] ;
  //       b_r /= beta[beta_idx.B][i] - beta[beta_idx.B][i-1] ;
  
  // return b_r ;
}

// ------------------------------------------------------------
// Returns b_M
Zaki::Vector::DataColumn BNV_Sequence::b_M() const 
{
  // return beta[beta_idx.beta_M] * beta[beta_idx.B] / beta[beta_idx.M] ;
  return beta[beta_idx.b_M] ;
}

// ------------------------------------------------------------
// Returns Beta(O)
Zaki::Vector::DataColumn BNV_Sequence::Beta(const Zaki::Vector::DataColumn& in_dc) const 
{
  Zaki::Vector::DataSet tmp_ds({beta[0], in_dc}) ;

  tmp_ds.Interpolate(0,1) ;
  tmp_ds.Derivative(1)[1] ;

  return tmp_ds.Derivative(1)[1] / beta.Derivative(beta_idx.B)[1] ;
}

// ------------------------------------------------------------
// Returns b(O)
Zaki::Vector::DataColumn BNV_Sequence::b(const Zaki::Vector::DataColumn& in_dc) const 
{
  Zaki::Vector::DataSet tmp_ds({beta[0], in_dc}) ;
  tmp_ds.Interpolate(0,1) ;

  return (tmp_ds.Derivative(1)[1] / beta.Derivative(beta_idx.B)[1]) * (beta[beta_idx.B] / in_dc) ;
}

// ------------------------------------------------------------
// Returns b(O) as a function of time
double BNV_Sequence::b(const Zaki::Vector::DataColumn& in_dc, const double& time) const 
{
  // int idx = B().GetClosestIdx( B(time) ) ;
  // return b(in_dc)[Time_to_Idx(time)] ;

  weighted_idx w_idx = Time_to_Weighted_Idx(time) ;

  return  w_idx.w_L * b(in_dc)[w_idx.idx-1] 
            + w_idx.w_R * b(in_dc)[w_idx.idx] ;

  // double B_val = B(time) ;
  // int i = beta[beta_idx.B].GetClosestIdx( B_val ) ;

  // if(B_val > beta[beta_idx.B][i])
  // {
  //   i++ ;
  // }
  
  // double b_o = (B_val - beta[beta_idx.B][i-1]) * b(in_dc)[i-1] ;
  //       b_o += (beta[beta_idx.B][i] - B_val) * b(in_dc)[i] ;
  //       b_o /= beta[beta_idx.B][i] - beta[beta_idx.B][i-1] ;
  
  // return b_o ;
}
// ------------------------------------------------------------
// Returns the closest index corrsponding to time
int BNV_Sequence::Time_to_Idx(const double& time) const
{
  return B().GetClosestIdx( B(time) ) ;
}

// ------------------------------------------------------------
// Returns the weighted index corrsponding to time
CompactStar::BNV_Sequence::weighted_idx BNV_Sequence::Time_to_Weighted_Idx(const double& time) const 
{
  weighted_idx w_idx ;

  double B_val = B(time) ;

  w_idx.idx = B().GetClosestIdx( B_val ) ;

  if(B_val > beta[beta_idx.B][w_idx.idx])
  {
    w_idx.idx++ ;
  }
  
  w_idx.w_L  = B_val - beta[beta_idx.B][w_idx.idx-1] ;
  w_idx.w_L /= beta[beta_idx.B][w_idx.idx] - beta[beta_idx.B][w_idx.idx-1] ;

  w_idx.w_R  = beta[beta_idx.B][w_idx.idx] - B_val ;
  w_idx.w_R /= beta[beta_idx.B][w_idx.idx] - beta[beta_idx.B][w_idx.idx-1] ;

  return w_idx ;
  // double b_i = (B_val - beta[beta_idx.B][i-1]) * b_I()[i-1] ;
  //       b_i += (beta[beta_idx.B][i] - B_val) * b_I()[i] ;
  //       b_i /= beta[beta_idx.B][i] - beta[beta_idx.B][i-1] ;
}

// ------------------------------------------------------------
// Returns the critical (death) omega assuming a central dipole H-field (per sec)
double BNV_Sequence::Omega_Death_Central(const double& time)
{
  double omega_D_I = pow(12 /* km */ / R(time), 19./15.) ;
  omega_D_I *= pow(1e8 /* gauss*/ / MagField() , 8./15.) ;
  omega_D_I *= 2*M_PI / 6.1e-3 ;

  return omega_D_I ;
}

// ------------------------------------------------------------
// Returns the critical (death) omega assuming a twisted dipole H-field (per sec)
double BNV_Sequence::Omega_Death_Twisted(const double& time)
{
  double omega_D_II = pow(12 /* km */ / R(time), 17./13.) ;
  omega_D_II *= pow(1e8 /* gauss*/ / MagField() , 7./13.) ;
  omega_D_II *= 2*M_PI / 25.4e-3 ;

  return omega_D_II ;
}

// ------------------------------------------------------------
// Returns the breaking index
Zaki::Vector::DataColumn BNV_Sequence::BrIdx(const double& delta) const 
{
  Zaki::Vector::DataColumn  n =  b_I() * delta / (1 - b_I() * delta )  ;
                           n += 3 ;
    // n += (6 * b_R() * delta - b_I() * pow(delta, 2) * (1+b(Beta_I())) ) / (1 - b_I() * delta).pow(2) ;
    n += (6 * b_R() * delta - b_I() * pow(delta, 2) * (1+b_beta_I()) ) / (1 - b_I() * delta).pow(2) ;

  return n ;
}

// ------------------------------------------------------------
// Returns the breaking index as a function of time and omega
double BNV_Sequence::BrIdx(const double& t, const double& omega) const 
{
  double delta = ( Gamma_BNV(t) / ODE_Coeff(t) ) / pow(omega, 2) ;

  // int idx = Time_to_Idx(t) ;
  weighted_idx w_idx = Time_to_Weighted_Idx(t) ;

  double b_I_val = w_idx.w_L * b_I()[w_idx.idx-1] 
            + w_idx.w_R * b_I()[w_idx.idx] ;
  
  double b_R_val = w_idx.w_L * b_R()[w_idx.idx-1] 
            + w_idx.w_R * b_R()[w_idx.idx] ;

  double b_beta_I_val = w_idx.w_L * b_beta_I()[w_idx.idx-1] 
            + w_idx.w_R * b_beta_I()[w_idx.idx] ;

  double  n = 3 + b_I_val * delta / (1 - b_I_val * delta )  ;
    // n += (6 * b_R_val * delta - b_I_val * pow(delta, 2) * (1+b(Beta_I(), t)) ) / pow(1 - b_I_val * delta, 2) ;
    n += (6 * b_R_val * delta - b_I_val * pow(delta, 2) * (1+b_beta_I_val) ) / pow(1 - b_I_val * delta, 2) ;
    
  return n ;
}

// ------------------------------------------------------------
// Returns the breaking index as a function of mass and omega
double BNV_Sequence::BrIdx_M_Omega(double mass, double omega) 
{
  int idx = M().GetClosestIdx(mass) ;

  double delta = ( gamma_bnv / H2R6_I()[idx] ) / pow(omega, 2) ;
  double b_i = b_I()[idx] ;
  double b_r = b_R()[idx] ;

  double n = 3 ;
         n+= b_i * delta / (1 - b_i * delta ) ; 
        // n += (6 * b_r * delta - b_i * pow(delta, 2) * (1+b(Beta_I())[idx]) ) 
        n += (6 * b_r * delta - b_i * pow(delta, 2) * (1+b_beta_I()[idx]) ) 
              / pow(1 - b_i * delta, 2) ;

  return n ;
}

// ------------------------------------------------------------
int BNV_Sequence::FindInitIdx(const double& in_mass) 
{
  init_idx = M().GetClosestIdx(in_mass) ;
  return init_idx ;
}

// ------------------------------------------------------------
// Returns the strength of magnetic field in gauss
double BNV_Sequence::MagField() const 
{
  return mag_field ;
}

// ------------------------------------------------------------
// Sets the strength of magnetic field in gauss
void BNV_Sequence::SetMagField(const double& in_B) 
{
  mag_field = in_B ;

  SetInitOmega() ;
}

// ------------------------------------------------------------
// Returns the threshold magnetic field at which
// the SD-SU transition occurs in gauss
//        bnv_rate in yr^-1
//        P in seconds
//        mass in solar mass
double BNV_Sequence::ThresholdMagField(const double& bnv_rate, const double& P, const double& mass) const 
{
  int idx = M().GetClosestIdx(mass);

  // Unit = [s^2] / [yr * km^3]
  double H_sqrd = (b_I()[idx] * bnv_rate * P*P * I()[idx]) / (4 * M_PI * M_PI * pow(R()[idx], 6) );

  // Unit = [ s / km^3 ]
  H_sqrd /= Zaki::Physics::YR_2_SEC ;

  // Unit = [ 1 / km^2 ]
  H_sqrd *= Zaki::Physics::LIGHT_C_KM_S ; 

  // Unit = [ 1 / km ]
  double H = sqrt(H_sqrd) ;

  // Unit = [ 1 / cm ]
  H /= 1e5 ;

  // Unit = [gauss]
  H *= 3.479e24 ;

  return H ;
}

// ------------------------------------------------------------
// Returns (H^2 R^6 / I) in seconds
double BNV_Sequence::ODE_Coeff(const double& time) const
{
  double tmp = pow(MagField(), 2) ;

  // Change gauss to km^-1
  tmp *= pow( 1e5 * 1e-24 / 3.479, 2 ) ;

  tmp *= pow(R(time), 6) ;

  tmp /= I(time) ;

  // Change km to seconds
  tmp /= Zaki::Physics::LIGHT_C_KM_S ;

  return tmp ;
}

// ------------------------------------------------------------
// Returns (H^2 R^6 / I) in seconds as a DataColumn 
Zaki::Vector::DataColumn BNV_Sequence::H2R6_I() const 
{
  Zaki::Vector::DataColumn H2R6_I = R().pow(6) ;

  H2R6_I *= pow(MagField(), 2) ;

  // Change gauss to km^-1
  H2R6_I *= pow( 1e5 * 1e-24 / 3.479, 2 ) ;

  H2R6_I /= I() ;

  // Change km to seconds
  H2R6_I /= Zaki::Physics::LIGHT_C_KM_S ;

  return H2R6_I ;
}

// ------------------------------------------------------------
// Returns the BNV rate per seconds
double BNV_Sequence::Gamma_BNV(const double& time) const 
{
  // Adding a cut-off for BNV 
  if ( time > bnv_cutoff_time ) 
    return 0 ;

  return gamma_bnv ;
}

// ------------------------------------------------------------
// Returns the extremum value for angular velocity (per seconds)
double BNV_Sequence::ExtremumOmega(const double& time) const 
{
  // "ODE_Coeff" returns (H^2 R^6 / I) in seconds
  double  omega  = sqrt ( b_I(time) * Gamma_BNV(time) / ODE_Coeff(time) ) ;

  return omega ;
}

// ------------------------------------------------------------
// Returns a dataset with two columns:
// [0]: Mass in solar mass
// [1]: Extremum value for angular velocity (per seconds)
Zaki::Vector::DataSet BNV_Sequence::ExtremumOmega() const 
{
  Zaki::Vector::DataSet  M_omega( { M(), ( b_I() * gamma_bnv / H2R6_I() ).sqrt() } ) ;

  return M_omega ;
}

// ------------------------------------------------------------
// Returns the second derivative of omega as a function of 
//  time and omega 
// [ I cross checked this with my notes from Jul 9, 2023 on Aug 10, 2023. ]
double BNV_Sequence::Omega_2ndDer( const double& t, 
                                  const double& omega) const 
{
  // double b_I_dot = b_I(t) * ( b_I(t) - 1 - b(Beta_I(), t)) * Gamma_BNV(t) ;
  double b_I_dot = b_I(t) * ( b_I(t) - 1 - b_beta_I(t)) * Gamma_BNV(t) ;

  double R_dot_over_R = -b_R(t) * Gamma_BNV(t) ;
  double I_dot_over_I = -b_I(t) * Gamma_BNV(t) ;

  double  omega_ddot  = 6 * R_dot_over_R ;
          omega_ddot += 3 * Omega_1stDer(t, omega) / omega ;
          omega_ddot -= I_dot_over_I ;
          omega_ddot *= - ODE_Coeff(t) * pow(omega, 3) ;
          omega_ddot += b_I_dot * Gamma_BNV(t) * omega ;
          omega_ddot += b_I(t) * Gamma_BNV(t) * Omega_1stDer(t, omega) ;

  return omega_ddot ;
}

// ------------------------------------------------------------
// Returns the second derivative of omega (DataColumn) 
// as a function of omega 
Zaki::Vector::DataColumn BNV_Sequence::Omega_2ndDer(const double& omega) 
  const 
{
  // Zaki::Vector::DataColumn b_I_dot = b_I() * ( b_I() - 1 - b(Beta_I())) * gamma_bnv ;
  Zaki::Vector::DataColumn b_I_dot = b_I() * ( b_I() - 1 - b_beta_I()) * gamma_bnv ;

  Zaki::Vector::DataColumn R_dot_over_R = -b_R() * gamma_bnv ;
  Zaki::Vector::DataColumn I_dot_over_I = -b_I() * gamma_bnv ;

  Zaki::Vector::DataColumn  omega_ddot  = 6 * R_dot_over_R ;
          omega_ddot += 3 * Omega_1stDer(omega) / omega ;
          omega_ddot -= I_dot_over_I ;
          omega_ddot *= - H2R6_I() * pow(omega, 3) ;
          omega_ddot += b_I_dot * gamma_bnv * omega ;
          omega_ddot += b_I() * gamma_bnv * Omega_1stDer(omega) ;

  return omega_ddot ;
}

// ------------------------------------------------------------
// Returns the first derivative of omega as a function of 
//  time and omega
double BNV_Sequence::Omega_1stDer(const double& t, const double& omega) const 
{
  double omega_dot = -ODE_Coeff(t) * pow(omega, 3) ;

  // Adding the BNV term
  omega_dot += b_I(t) * Gamma_BNV(t) * omega ;

  return omega_dot ;
}

// ------------------------------------------------------------
// Returns the first derivative of omega (DataColumn) 
// as a function of omega 
Zaki::Vector::DataColumn BNV_Sequence::Omega_1stDer(const double& omega) const 
{
  Zaki::Vector::DataColumn omega_dot = -H2R6_I() * pow(omega, 3) ;

  // Adding the BNV term
  omega_dot += b_I() * gamma_bnv * omega ;

  return omega_dot ;
}

// ------------------------------------------------------------
// ......................
// Dictionary :
// y[0] = Omega(t)
// f[0] = Omega'(t)
// ......................
int BNV_Sequence::ODE(double t, const double y[], double f[], void *params) 
{
  BNV_Sequence* ns_obj = (BNV_Sequence *)params ; 

  f[0]  = ns_obj->ODE_Coeff(t) ;

  f[0] *= - pow(y[0], 3) ;

  // Adding the BNV term
  f[0] += ns_obj->b_I(t) * ns_obj->Gamma_BNV(t) * y[0] ;

  return GSL_SUCCESS ;
}

// ------------------------------------------------------------
void BNV_Sequence::Solve(const double t_0, const double t_f) 
{
  //---------------------------------------------------
  //          omega_ds Definition
  //---------------------------------------------------
  Zaki::Vector::DataSet omega_ds(1+ 2*init_omega_set.size(), time_res) ;
  omega_ds.SetWrkDir(wrk_dir) ;

  omega_ds[0].label = "t (yr)" ;

  for (size_t i = 1; i < 2*init_omega_set.size(); i+=2)
  {
    omega_ds[i].label = "N_" + init_omega_set[(i-1)/2].second ;
    omega_ds[i+1].label = "A_" + init_omega_set[(i-1)/2].second ;
  }
  
  std::vector<std::pair<int, std::string>> omega_ds_plt_labels ; 
  //---------------------------------------------------
  //---------------------------------------------------
  //          br_idx_ds Definition
  //---------------------------------------------------
  Zaki::Vector::DataSet br_idx_ds(1 + 2*init_omega_set.size(), time_res) ;
  br_idx_ds.SetWrkDir(wrk_dir) ;
  br_idx_ds[0].label = "t (yr)" ;

  std::vector<
              std::pair<int, 
                std::map<std::string, std::string> 
                >>  
    br_idx_ds_plt_labels ;
  //---------------------------------------------------

  //---------------------------------------------------
  //          P_Pdot_ds Definition
  //---------------------------------------------------
  Zaki::Vector::DataSet P_Pdot_ds(2 + 3*init_omega_set.size(), time_res) ;
  P_Pdot_ds.SetWrkDir(wrk_dir) ;
  P_Pdot_ds[0].label = "t (yr)" ;
  P_Pdot_ds[1].label = "M" ;

  std::vector<
              std::pair<int, 
                std::map<std::string, std::string> 
                >>  
    P_Pdot_ds_plt_labels ;
  //---------------------------------------------------

  //---------------------------------------------------
  //          M_Omega_ds Definition
  //---------------------------------------------------
  Zaki::Vector::DataSet M_Omega_ds(5 + init_omega_set.size(), time_res) ;
  M_Omega_ds.SetWrkDir(wrk_dir) ;
  M_Omega_ds[0].label = "t (yr)" ;
  M_Omega_ds[1].label = "M (M_sun)" ;
  M_Omega_ds[2].label = "$\\Omega_e (1/s)$" ;
  M_Omega_ds[3].label = "$\\Omega_{D-I} (1/s)$" ;
  M_Omega_ds[4].label = "$\\Omega_{D-II} (1/s)$" ;

  std::vector<
              std::pair<int, 
                std::map<std::string, std::string> 
                >>  
    M_Omega_ds_plt_labels ;
  
  // This holds the mass values corresponding to specific time intervals
  // it is plotted as vertical lines on the (omega,M) plot, down below.
  std::vector<double> time_mass_stamps ;
  std::vector<double> time_stamps = { 1e-2/gamma_bnv, 
                                     1e-1/gamma_bnv,
                                     5e-1/gamma_bnv, 
                                     1e0/gamma_bnv} ;

  //---------------------------------------------------

  M_Omega_ds_plt_labels.emplace_back(2, 
          (std::map<std::string, std::string>){{"label", "$\\Omega_e$"}, {"ls", "--"}
                                          }) ;
  M_Omega_ds_plt_labels.emplace_back(3, 
          (std::map<std::string, std::string>){{"label", "$\\Omega_{D-I}$"}, {"ls", "--"}
                                          }) ;
  M_Omega_ds_plt_labels.emplace_back(4, 
          (std::map<std::string, std::string>){{"label", "$\\Omega_{D-II}$"}, {"ls", "--"}  
                                          }) ;                                        
  //---------------------------------------------------
  //    Loop over initial omega values
  //---------------------------------------------------
  for (size_t i = 0; i < init_omega_set.size(); i++)
  {

    double omega[1] ;
    omega[0] = init_omega_set[i].first ;
    double in_t = t_0 ;

    //----------------------------------------
    //          GSL ODE SYSTEM SETUP
    //----------------------------------------
    gsl_odeiv2_system ode_sys = {BNV_Sequence::ODE, nullptr, 1, this} ;

    gsl_odeiv2_driver *tmp_driver = gsl_odeiv2_driver_alloc_y_new
        (&ode_sys, gsl_odeiv2_step_rk8pd,
        1e-3, 1e-10, 1e-10);
    //----------------------------------------
    double log_t_0 = log10(t_0) ;
    double log_t_f = log10(t_f) ;

    double step = (log_t_f - log_t_0) / time_res  ;
    double alpha = ODE_Coeff(t_0) ;

    //----------------------------------------
    //            Loop over time
    //----------------------------------------
    for (size_t j = 0; j < time_res; j++ )
    // for (double log_t_i = log_t_0; log_t_i < log_t_f; log_t_i += step )
    {

      double log_t_i = log_t_0 + j * step ;

      double t_i = pow(10, log_t_i) ;

      SetBNVCuttoff(t_i) ;
      
      int status = gsl_odeiv2_driver_apply (tmp_driver, &in_t, t_i, omega);
      
      //----------------------------------------
      //          Recording Data
      //----------------------------------------
      // Recording the time has to be done only once!
      if ( i == 0)
      {
        omega_ds[0].vals.emplace_back(t_i / Zaki::Physics::YR_2_SEC) ;
        br_idx_ds[0].vals.emplace_back(t_i / Zaki::Physics::YR_2_SEC) ;
        P_Pdot_ds[0].vals.emplace_back(t_i / Zaki::Physics::YR_2_SEC) ;
        P_Pdot_ds[1].vals.emplace_back(M(t_i)) ;
        M_Omega_ds[0].vals.emplace_back(t_i / Zaki::Physics::YR_2_SEC) ;
        M_Omega_ds[1].vals.emplace_back(M(t_i)) ;
        // std::cout<< "\t t_i = '" << t_i << "', B = " << B(t_i) << ", " << " M(t_i) = '" << M(t_i) << "'.\n" ;
        M_Omega_ds[2].vals.emplace_back( ExtremumOmega(t_i) ) ;
        M_Omega_ds[3].vals.emplace_back( Omega_Death_Central(t_i) ) ;
        M_Omega_ds[4].vals.emplace_back( Omega_Death_Twisted(t_i) ) ;

        M_Omega_ds[0].label = "t[yr]" ;
        M_Omega_ds[1].label = "M" ;
        M_Omega_ds[2].label = "Omega_ext" ;
        M_Omega_ds[3].label = "Omega_D_I" ;
        M_Omega_ds[4].label = "Omega_D_II" ;

        for (size_t k = 0; k < time_stamps.size(); k++)
        {
          if( abs(log_t_i - log10(time_stamps[k])) < step/2.  )
          {
            time_mass_stamps.emplace_back(M(t_i)) ;
          }
        }
      }

      omega_ds[1 + 2 * i].vals.emplace_back(omega[0] / init_omega_set[i].first) ;

      // .....................
      // Analytical Solution
      // .....................
      omega_ds[2 + 2 * i].vals.emplace_back(
        sqrt(1. / 
              (2*alpha*(t_i - t_0) + pow(init_omega_set[i].first, -2)) 
            ) / init_omega_set[i].first                     
                                    ) ;
      // .....................
  
      // .....................
      //    Braking index
      // .....................
      // br_idx_ds[1 + 2 * i].vals.emplace_back(Omega_2ndDer(t_i, omega[0]) * omega[0]
      //                                  / pow(Omega_1stDer(t_i, omega[0]),2)) ;
      br_idx_ds[2 + 2 * i].vals.emplace_back(BrIdx(t_i, omega[0])) ;
      // br_idx_ds[1 + 2 * i].vals.emplace_back(b_R(t_i)) ;
      // br_idx_ds[1 + 2 * i].vals.emplace_back(ODE_Coeff(t_i)) ;
      // br_idx_ds[2 + 2 * i].vals.emplace_back(b_beta_I(t_i)) ;

      // .....................

      // .....................
      //    P, P_dot
      // .....................
      P_Pdot_ds[2 + 3 * i].vals.emplace_back( 2*M_PI / omega[0]) ;
      P_Pdot_ds[3 + 3 * i].vals.emplace_back( -2*M_PI * Omega_1stDer(t_i, omega[0]) / pow(omega[0], 2)) ;
      // P_Pdot_ds[3 + 3 * i].vals.emplace_back( - Pdot_sign * omega[0] / ( 2 * Zaki::Physics::YR_2_SEC 
                                              // * Omega_1stDer(t_i, omega[0]) )  ) ;
      P_Pdot_ds[4 + 3 * i].vals.emplace_back( - omega[0] / ( 2 * Zaki::Physics::YR_2_SEC 
                                              * Omega_1stDer(t_i, omega[0]) )  ) ;                                         
      P_Pdot_ds[2 + 3 * i].label = "P[" + init_omega_set[i].second + "]";
      P_Pdot_ds[3 + 3 * i].label = "P_dot[" + init_omega_set[i].second + "]";
      P_Pdot_ds[4 + 3 * i].label = "SD_Age[" + init_omega_set[i].second + "]" ;

      // .....................

      // .....................
      //    M_Omega
      // .....................
      M_Omega_ds[5+i].vals.emplace_back(omega[0]) ;
      M_Omega_ds[5+i].label = "Omega[" + init_omega_set[i].second + "]";
      // .....................

      //----------------------------------------

      if (status != GSL_SUCCESS)
      {
        printf("\t-------------------%s-------------------\n", "GSL") ;
        printf ("error, return value=%d\n.", status);
        printf("Omega = %2.2e.\n", omega[0]) ;
        printf("T = %2.2e.\n", t_i) ;
        break;
      }
    }
    //----------------------------------------
    //          End of loop over time
    //----------------------------------------
    printf("* ---------------H = %2.1e gauss---------------- *\n"
           "\tOmega_i = %2.2e.\n", mag_field, init_omega_set[i].first) ;
    printf("\tOmega_f = %2.2e.\n"
           "* ------------------------------- *\n", omega[0]) ;

    gsl_odeiv2_driver_free (tmp_driver) ;

    omega_ds_plt_labels.emplace_back(1 + 2 * i, omega_ds[1 + 2 * i].label) ;
    omega_ds_plt_labels.emplace_back(2 + 2 * i, omega_ds[2 + 2 * i].label) ;
    
    // br_idx_ds_plt_labels.emplace_back(2 + 2 * i, 
    //   (std::map<std::string, std::string>){
    //     {"label", "(33)_" + init_omega_set[i].second}, {"ls", "-"}
    //                                       }) ;
    br_idx_ds_plt_labels.emplace_back(2 + 2 * i, 
      (std::map<std::string, std::string>){
        {"label", init_omega_set[i].second}, {"ls", "--"}
                                          }) ;
    // br_idx_ds_plt_labels.emplace_back(1 + 2 * i, 
    //   (std::map<std::string, std::string>){
    //     {"label", init_omega_set[i].second},
    //                                       }) ;

    P_Pdot_ds_plt_labels.emplace_back(3 + 3 * i, 
      (std::map<std::string, std::string>){
        {"label", init_omega_set[i].second}, {"ls", "--"}
                                          }) ;
    
    M_Omega_ds_plt_labels.emplace_back(5 + i, 
      (std::map<std::string, std::string>){
        {"label", init_omega_set[i].second}
                                          }) ;
  }
  //---------------------------------------------------
  //      End of loop over omega initial values
  //---------------------------------------------------

  //---------------------------------------------------
  //            Plotting Data
  //---------------------------------------------------
  Zaki::Vector::DataSet::PlotParam plt_pars ;
  plt_pars.SetGrid() ;
  plt_pars.SetLegend({"upper left", 0.0, 1.0}) ;
  plt_pars.SetXAxis({(1e-3/gamma_bnv)/ Zaki::Physics::YR_2_SEC, (10/gamma_bnv)/ Zaki::Physics::YR_2_SEC}) ;
  // plt_pars.SetXAxis({5e-4, 1e1}) ;
  // plt_pars.SetYAxis({1e-50, 1e-12}) ;
  // plt_pars.SetYAxis({0.75, 1.25}) ;

  // for (size_t i = 0; i < init_omega_set.size(); i++)
  // {
  //   plt_pars.AddAxHLine(ExtremumOmega(0) / init_omega_set[i].first, {{"ls", "--"}}) ;
  // }

  // std::cout << "\n\n\t ExtremumOmega (t=0) = " << ExtremumOmega(0) << " per second." ;
  // std::cout << "\n\t Extremum Period (t=0) = " << 2*M_PI / ExtremumOmega(0) << " second.\n" ;
  char file_stamps_char[100] ;
  snprintf(file_stamps_char, sizeof(file_stamps_char), "_G=%.1e_H=%.1e_M_i=%.2f_Mc=%.1f", 
          gamma_bnv*Zaki::Physics::YR_2_SEC, mag_field, M()[init_idx], bnv_cuttoff_mass ) ;
  std::string file_stamps_str(file_stamps_char) ;

  // char omega_ds_fname[100] ;
  // snprintf(omega_ds_fname, sizeof(omega_ds_fname), omega_t_evol_f_name + "_G=%.1e_H=%.1e_M_i=%.2f_Mc=%.1f.pdf", 
  //         gamma_bnv*Zaki::Physics::YR_2_SEC, mag_field, M()[init_idx], bnv_cuttoff_mass ) ;

  omega_ds.SetPlotPars(plt_pars) ;
  omega_ds.SemiLogXPlot(0, omega_ds_plt_labels, omega_t_evol_f_name + file_stamps_str + ".pdf") ;

  // char br_idx_ds_fname[100] ;
  // snprintf(br_idx_ds_fname, sizeof(br_idx_ds_fname), br_idx_evol_f_name + "_G=%.1e_H=%.1e_M_i=%.2f_Mc=%.1f.pdf", 
  //         gamma_bnv*Zaki::Physics::YR_2_SEC, mag_field, M()[init_idx], bnv_cuttoff_mass ) ;

  // plt_pars.SetYAxis({0, 6}) ;
  plt_pars.SetYAxis({-10, 100}) ;
  plt_pars.SetXAxis({(1e-3/gamma_bnv)/ Zaki::Physics::YR_2_SEC, (10/gamma_bnv)/ Zaki::Physics::YR_2_SEC}) ;
  plt_pars.SetYAxisLabel("$n$") ;
  // plt_pars.AddAxHLine(3, {{"c", "red"}, {"ls", "-."}}) ;
  plt_pars.SetLegend({"lower left", 0.0, 0.0}) ;
  br_idx_ds.SetPlotPars(plt_pars) ; 

  // br_idx_ds.MakeSmooth(5) ;
  br_idx_ds.SemiLogXPlot(0, br_idx_ds_plt_labels, br_idx_evol_f_name + file_stamps_str + ".pdf") ;
  // br_idx_ds.LogLogPlot(0, br_idx_ds_plt_labels, br_idx_ds_fname) ;

  plt_pars.Reset() ;
  plt_pars.SetGrid() ;
  plt_pars.SetXAxis({(1e-3/gamma_bnv)/ Zaki::Physics::YR_2_SEC, (10/gamma_bnv)/ Zaki::Physics::YR_2_SEC}) ;
  plt_pars.SetYAxis({1e6, 1e11}) ;
  // plt_pars.SetYAxis({4.5e9, 5.1e9}) ;
  // plt_pars.SetYAxis({4e9, 1.2e10}) ;
  plt_pars.SetYAxisLabel("$P / 2\\dot{P}$") ;
  plt_pars.SetXAxisLabel("$t\\, [ yr ]$") ;
  plt_pars.SetLegend({"lower left", 0.0, 0.0}) ;
  P_Pdot_ds.SetPlotPars(plt_pars) ; 

  P_Pdot_ds.MakeSmooth(5) ;

  // char P_Pdot_ds_fname[150] ;
  // snprintf(P_Pdot_ds_fname, sizeof(P_Pdot_ds_fname), age_t_evol_f_name + "_G=%.1e_H=%.1e_M_i=%.2f_Mc=%.1f.pdf", 
  //           gamma_bnv*Zaki::Physics::YR_2_SEC, mag_field, M()[init_idx], bnv_cuttoff_mass ) ;

  P_Pdot_ds.LogLogPlot(0, P_Pdot_ds_plt_labels, age_t_evol_f_name + file_stamps_str + ".pdf") ;

  // snprintf(P_Pdot_ds_fname, sizeof(P_Pdot_ds_fname), age_t_evol_f_name + "_G=%.1e_H=%.1e_M_i=%.2f_Mc=%.1f.tsv", 
  //           gamma_bnv*Zaki::Physics::YR_2_SEC, mag_field, M()[init_idx], bnv_cuttoff_mass ) ;
  P_Pdot_ds.Export(age_t_evol_f_name + file_stamps_str + ".tsv") ;
  // --------------------
  //      M_Omega
  // --------------------
  plt_pars.Reset() ;
  plt_pars.SetGrid() ;
  plt_pars.SetXAxis({0.05, 2.1}) ;
  if(mag_field == 1e8 || mag_field == 1e9)
    plt_pars.SetYAxis({5, 2e4}) ; // M_Omega Evol: H = {1e8, 1e9}
  else
    plt_pars.SetYAxis({1, 2e3}) ; // M_Omega Evol: H = {1e10, 1e11}
  // plt_pars.SetYAxis({10, 1e3}) ; // Schematic
  plt_pars.SetYAxisLabel("$\\Omega\\, [1/s]$") ;
  plt_pars.SetXAxisLabel("$M\\, [ M_{\\rm sun} ]$") ;
  plt_pars.SetLegend({"lower left", 0.0, 0.0}) ;

  for (size_t i = 0; i < time_mass_stamps.size(); i++)
  {
    plt_pars.AddAxVLine(time_mass_stamps[i], {{"ls", "-."}}) ;
  }
  
  // Adding M_TOV guideline
  plt_pars.AddAxVLine(M().Max(), {{"ls", ":"}}) ;

  M_Omega_ds.SetPlotPars(plt_pars) ; 
  // char M_Omega_ds_fname[150] ;
  // snprintf(M_Omega_ds_fname, sizeof(M_Omega_ds_fname), omega_m_evol_f_name + "_G=%.1e_H=%.1e_M_i=%.2f_Mc=%.1f.pdf", 
  //           gamma_bnv*Zaki::Physics::YR_2_SEC, mag_field, M()[init_idx], bnv_cuttoff_mass ) ;

  M_Omega_ds.SemiLogYPlot(1, M_Omega_ds_plt_labels, omega_m_evol_f_name + file_stamps_str + ".pdf") ;
  
  // snprintf(M_Omega_ds_fname, sizeof(M_Omega_ds_fname), omega_m_evol_f_name + "_G=%.1e_H=%.1e_M_i=%.2f_Mc=%.1f.tsv", 
  //           gamma_bnv*Zaki::Physics::YR_2_SEC, mag_field, M()[init_idx], bnv_cuttoff_mass ) ;
  M_Omega_ds.Export(omega_m_evol_f_name + file_stamps_str + ".tsv") ;
  // --------------------

}

// ------------------------------------------------------------
// Sets file names for braking index vs time evolution 
// (see "void Solve(t_i, t_f)" and br_idx_evol_f_name)
void BNV_Sequence::SetBrIdxvsTimeFileName(const std::string& f_name) 
{
  br_idx_evol_f_name = f_name ;
}

// ------------------------------------------------------------
// File names for omega vs time evolution
// (see "void Solve(t_i, t_f)" and omega_t_evol_f_name)
void BNV_Sequence::SetOmegavsTimeFileName(const std::string& f_name) 
{
  omega_t_evol_f_name = f_name ;
}

// ------------------------------------------------------------
// File names for age vs time evolution
// (see "void Solve(t_i, t_f)" and age_t_evol_f_name)
void BNV_Sequence::SetAgevsTimeFileName(const std::string& f_name) 
{
  age_t_evol_f_name = f_name ;
}

// ------------------------------------------------------------
// File names for omega vs mass evolution
// (see "void Solve(t_i, t_f)" and omega_m_evol_f_name)
void BNV_Sequence::SetOmegavsMFileName(const std::string& f_name) 
{
  omega_m_evol_f_name = f_name ;
}

// ------------------------------------------------------------
// Returns the critical period at which the pulsar turns off
//  Taken from Eq. (39) of 
//  "Theory of Pulsars: Polar Gaps, Sparks, and Coherent Microwave Radiation"
//  by Ruderman & Sutherland (1974)
double BNV_Sequence::DeathPeriod(const double& t) const 
{
  double rho_6 = 1 ;
  double chi   = 1. / 15. ;

  double R_6 = R(t) / 10 ;
  double B_12 = MagField() / 1e12 ;

  double p_crit  = 1.7 * pow(B_12, 8./13.) * pow(R_6, 21./13) ;
         p_crit *= pow(rho_6, -4./13.) ;
         p_crit *= pow(15 * chi, -2./3.) ;

  return p_crit ;
}

// ------------------------------------------------------------
// Returns the critical omgea at which the pulsar turns off
//  Taken from Eq. (39) of 
//  "Theory of Pulsars: Polar Gaps, Sparks, and Coherent Microwave Radiation"
//  by Ruderman & Sutherland (1974)
double BNV_Sequence::DeathOmega(const double& t) const 
{
  return 2 * M_PI / DeathPeriod(t) ;
}

// ------------------------------------------------------------
// Returns the second derivative of frequency (at t=0)
//  as a function of the frequency and the first derivative
//  nu and nu_dot must be in Hz and Hz / s
//  output is in [s^-3]
double BNV_Sequence::Freq_ddot(const double& nu, 
                              const double& nu_dot) const 
{
  double nu_ddot = - (b_I(0) + 6 * b_R(0)) * gamma_bnv * nu_dot ;
        nu_ddot += b_I(0) 
                  // * (6*b_R(0) - 1 - b(Beta_I(), 0)) 
                  * (6*b_R(0) - 1 - b_beta_I(), 0) 
                  * pow(gamma_bnv,2) * nu ;
                  
  return nu_ddot ;
}

// ------------------------------------------------------------
// Returns the BNV rate for which a maximum value of 
// the second derivative of frequency (at t=0)
//  is achieved.
//  output is dimensionless
// Has to be multiplied by P_dot / P
double BNV_Sequence::Gamma_Max() const 
{
  double gamma = b_I(0) + 6 * b_R(0) ;
         gamma /= b_I(0) 
                  // * (1 + b(Beta_I(), 0) - 6*b_R(0)) ;
                  * (1 + b_beta_I(0) - 6*b_R(0)) ;
                  
  return gamma ;
}
// ------------------------------------------------------------