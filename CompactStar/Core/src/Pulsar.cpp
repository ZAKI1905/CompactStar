/*
  Particle Class
*/

// // Creating directories
// #include <sys/stat.h>
// #include <filesystem>

// #include <Zaki/Math/IntegralTable.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/Core/Pulsar.hpp"
#include "CompactStar/EOS/CompOSE_EOS.hpp"


// -------------------------------------------------------
std::ostream& operator << (std::ostream & output, 
                const CompactStar::Baryon_Lim& baryon_lim)
{
  char tmp[100] ;
  snprintf( tmp, sizeof(tmp), 
            "%-10s B-Fraction = %10.8e\t  B_dot = %10.8e\t \u0393 = %14.8e  (1/yr)", 
           baryon_lim.name.c_str(), baryon_lim.fraction,
           baryon_lim.b_dot, baryon_lim.limit) ;
  output << tmp ;

  return output ;
}

//--------------------------------------------------------------
namespace CompactStar
{

//==============================================================
//             Pulsar Class
//==============================================================
// Default Constructor
Pulsar::Pulsar()
  : mp({0, 0}), spin_p({0,0}), spin_p_dot({0,0})
{ }
//--------------------------------------------------------------
// Constructor
Pulsar::Pulsar(  const std::string& name, 
                  const Zaki::Math::Quantity& m,
                  const Zaki::Math::Quantity& p,
                  const Zaki::Math::Quantity& pdot)
  : Prog(name, true), mp(m), spin_p(p), spin_p_dot(pdot)
{ }
//--------------------------------------------------------------
// Destructor
Pulsar::~Pulsar() 
{ }

// -------------------------------------------------------
// Sets pulsar mass
void Pulsar::SetMass(const Zaki::Math::Quantity& in_mass) 
{
  mp = in_mass ;
}

// -------------------------------------------------------    
// Sets pulsar spin period in s
void Pulsar::SetSpinP(const Zaki::Math::Quantity& in_spin_p) 
{
  spin_p = in_spin_p ;
}

// -------------------------------------------------------    
// Sets pulsar spin period derivative
void Pulsar::SetSpinPDot(const Zaki::Math::Quantity& in_spin_p_dot) 
{
  spin_p_dot = in_spin_p_dot ;
}

// -------------------------------------------------------
// Returns pulsar spin-period in seconds
Zaki::Math::Quantity Pulsar::GetSpinP() const 
{
  return spin_p ;
}

// -------------------------------------------------------
/// Returns pulsar mass
Zaki::Math::Quantity Pulsar::GetMass() const
{
  return mp ;
}

// -------------------------------------------------------
/// Returns pulsar spin period derivative
Zaki::Math::Quantity Pulsar::GetSpinPDot() const 
{
  return spin_p_dot ;
}
// -------------------------------------------------------
// Pulsar spin period derivative divided by period:
//        (P_dot / P) [1/s]
Zaki::Math::Quantity Pulsar::GetSpinPDot_over_P() const 
{
  return spin_p_dot / spin_p ;
}

// -------------------------------------------------------
/// Sets distance in kpc :
/// between the Solar system Barycentre (SSB)
///  and the position of the isolated pulsar or
///  the binary barycentre in the case of a binary pulsar
void Pulsar::SetDistance(const Zaki::Math::Quantity& in_d) 
{
  d = in_d ;
}

// -------------------------------------------------------
/// Sets proper motion of the pulsar
/// Unit: milli-arcsec per year
void Pulsar::SetProperMotion(const Zaki::Math::Quantity& in_mu) 
{
  mu = in_mu ;
}

// -------------------------------------------------------
/// Returns distance in kpc :
/// between the Solar system Barycentre (SSB)
///  and the position of the isolated pulsar or
///  the binary barycentre in the case of a binary pulsar
Zaki::Math::Quantity Pulsar::GetDistance() const 
{
  return d ;
}

// -------------------------------------------------------
/// returns proper motion of the pulsar
/// Unit: milli-arcsec per year
Zaki::Math::Quantity Pulsar::GetProperMotion() const 
{
  return mu ;
}

// -------------------------------------------------------
int Pulsar::FindProfile(const std::string& model_name, 
                        const Zaki::String::Directory& in_dir)
{
  Z_LOG_INFO("Finding the pulsar profile...") ;

  seq_profile.Import( ( wrk_dir + in_dir ) +
                      model_name + "_Sequence.tsv") ;

  int tmp_idx = seq_profile[1].GetClosestIdx(mp.val) ;
  SeqPoint seq_i( seq_profile.GetDataRows()[tmp_idx].vals ) ;
  Zaki::Vector::DataSet prof_i( (wrk_dir + in_dir) + "/profiles", 
                                model_name + "_" + 
                                std::to_string(tmp_idx) +
                                ".tsv")  ;

  // If the mass value of closest index is less than mp 
  // and is the maximum mass for that EoS, then we should
  // just take the maximum mass profile
  if (seq_profile[1][tmp_idx] <= mp.val 
        && 
      tmp_idx == seq_profile[1].MaxIdx())
  // if(true)
  {
    seq_point =  seq_i ; 
    profile = prof_i ;
  }
  // This is buggy because it uses interpolation 
  // and sometimes it generates densities lower than what
  // is allowed by the EoS
  else
  {
    SeqPoint seq_i_1( seq_profile.GetDataRows()[tmp_idx-1].vals ) ;
    
    // .....................................................
    // Since the pulsar mass and our numerical masses
    // might not match exactly, we can linearly interpolate
    // the profile between i-1 and i indices as
    // x P_{i-1} + (1-x) P_{i}
    
    double m_i_1 = seq_i_1.m ;
    double m_i = seq_i.m ;

    // Finding the proximity of the pulsar mass
    // to the m_i :
    double x = ( m_i - mp.val ) / ( m_i - m_i_1 ) ;

    std::cout << "----------------------------------------------------\n" ;
    std::cout << " The closest pulsar index is '" << tmp_idx << "':\n" ;
    std::cout << " \tM_{"<<tmp_idx-1<<"} = " << m_i_1 << "\n" ;
    std::cout << " \tM_{"<<tmp_idx<<"} = " << m_i << "\n" ;
    std::cout << " \tx = " << x << "\n" ;
    std::cout << "----------------------------------------------------\n" ;

    seq_point =  seq_i_1 * x + seq_i * (1-x) ; 

    std::cout << " Seq ["<<tmp_idx-1<<"] = " << seq_i_1.Str() << "\n" ;
    std::cout << " Seq ["<<tmp_idx<<"] = " << seq_i.Str() << "\n" ;
    std::cout << " Ideal Seq_point = " << seq_point.Str() << "\n" ;
    std::cout << " Ideal Mom_I [km^2 M_sun] = " << seq_point.I / Zaki::Physics::SUN_M_KM << "\n" ;
    
    // Interpolating the profile
    Zaki::Vector::DataSet prof_i_1( (wrk_dir + in_dir) + "/profiles", 
                                model_name + "_" + 
                                std::to_string(tmp_idx -1 ) +
                                ".tsv")  ;

    // std::cout << " R_i = " << prof_i[0][-1] << "\n" ;      
    // std::cout << " R_{i-1} = " << prof_i_1[0][-1] << "\n" ;      

    Zaki::Vector::DataSet short_prof ;
    Zaki::Vector::DataSet long_prof ;
    double x_long, x_short ;
    if( prof_i[0][-1] < prof_i_1[0][-1] )
    {
      short_prof  = prof_i ;
      long_prof   = prof_i_1 ;
      x_long      = x   ;
      x_short     = 1-x ;
    }
    else
    {
      short_prof = prof_i_1 ;
      long_prof  = prof_i ;
      x_long      = 1-x   ;
      x_short     = x ;
    }

    profile.data_set = { short_prof[0] } ;

    // Looping over columns
    for (size_t c = 1; c < short_prof.Dim().size() ; c++)
    {
      // std::cout << " \t c = " << c << "\n" ;
      long_prof.Interpolate(0, c) ;
      Zaki::Vector::DataColumn dc ;
      dc.label = short_prof[c].label ;
      dc.Reserve(short_prof[c].Size()) ;

      // std::cout << "c = " << c << " short_prof[0].Size() = " << short_prof[0].Size() << "\n" ;
      // Looping over rows
      for (size_t i = 0; i < short_prof[0].Size() ; i++)
      {
        dc.vals.emplace_back ( 
                        x_long * long_prof.Evaluate(c, short_prof[0][i] ) +
                        x_short * short_prof[c][i]
                              ) ;
      }
      profile.data_set.emplace_back(dc) ;
    }
  }
  // seq_point =  seq_i_1 + seq_i ; 
  
  // Z_LOG_INFO("Pulsar sequence point is given by:") ;
  // std::cout << "\t " << seq_point.Str() << "\n" ;
  // std::cout << "----------------------------------------------------\n" ;
  // .....................................................

  seq_profile.Interpolate(0, 4) ;
  double dB_over_deps = seq_profile.Derivative(4, seq_point.ec ) ;

  seq_profile.Interpolate(0, 5) ;
  double dI_over_deps = seq_profile.Derivative(5, seq_point.ec ) ;

  eta_I = seq_point.b / dB_over_deps ;
  eta_I *= dI_over_deps / seq_point.I ;  
  // std::cout << "eta_I = " << eta_I << "\n" ;
  // .....................................................

  // Zaki::String::Directory pulsar_dir = wrk_dir + "/" + name ;
  Zaki::String::Directory pulsar_dir = wrk_dir ;
  pulsar_dir.Create() ;
  // ............ Checking/Making the pulsar directory ............
  // if (!std::filesystem::is_directory(pulsar_dir))
  // {
  //   if (mkdir(pulsar_dir.c_str(),
  //               ACCESSPERMS) == -1) 
  //   {
  //     Z_LOG_INFO("Directory '"
  //                 + pulsar_dir
  //                 + "' wasn't created, because: "
  //                 + strerror(errno)+".") ;    
  //   }
  //   else
  //     Z_LOG_INFO("Directory '" + pulsar_dir + "' created.") ;
  // }
  // ..............................................................

  // profile.SetWrkDir(wrk_dir + "/" + name) ;
  profile.SetWrkDir(wrk_dir) ;
  profile.SetPrecision(12) ;
  profile.Export( in_dir + name + ".tsv") ;


  // Setting the heat blanket radius
  // defined by eps ~ 10^{10} [ g / cm^3 ] = 7.4237e-9 km^{-2}
  r_blanket_idx =  profile[4].GetClosestIdx(7.4237e-9) ;
  r_blanket = profile[0][r_blanket_idx] ;

  // The threshold radius for direct Urca: p_f(n) = p_F(p) + p_F(e)
  Zaki::Vector::DataColumn kF_n = (3*M_PI*M_PI*profile[5]*profile["10"]).pow(1.0/3.0) ;
  Zaki::Vector::DataColumn kF_p = (3*M_PI*M_PI*profile[5]*profile["11"]).pow(1.0/3.0) ;
  Zaki::Vector::DataColumn kF_e = (3*M_PI*M_PI*profile[5]*profile["0"]).pow(1.0/3.0) ;

  Zaki::Vector::DataColumn kF_diff = kF_n - kF_p - kF_e ;
  
  int durca_idx = kF_diff.GetClosestIdx(0) ;

  r_durca_thresh = profile[0][durca_idx] ;
  
  // Z_LOG_INFO("Finding the direct Urca threshold.") ;
  std::cout << "\n\n Pulsar::FindPulsar : "
            << " kF_diff = " << kF_diff[durca_idx]
            << ", R_DUrca = " 
            << r_durca_thresh << " [km].\n\n" ;


  r_durca_cond.Resize(profile[0].Size()) ;
  r_durca_cond.Fill(0) ;

  for (size_t i = 0; i < durca_idx; i++)
  {
    r_durca_cond[i] = 1 ;
  }
  
  return tmp_idx ;
}

// -------------------------------------------------------
// Imports the pulsar profile
void Pulsar::ImportProfile(const std::string& model_name, 
                   const Zaki::String::Directory& in_dir) 
{
  seq_profile.Import( ( wrk_dir + in_dir ) +
                      model_name + "_Sequence.tsv") ;

  int tmp_idx = seq_profile[1].GetClosestIdx(mp.val) ;

  seq_point =  seq_profile.GetDataRows()[tmp_idx].vals ;

  seq_profile.Interpolate(0, 4) ;
  double dB_over_deps = seq_profile.Derivative(4, seq_point.ec ) ;

  seq_profile.Interpolate(0, 5) ;
  double dI_over_deps = seq_profile.Derivative(5, seq_point.ec ) ;

  eta_I = seq_point.b / dB_over_deps ;
  eta_I *= dI_over_deps / seq_point.I ;  
  std::cout << "eta_I = " << eta_I << "\n" ;

  profile.Import( (wrk_dir + in_dir) + "/"+name+".tsv")  ;
  profile.SetWrkDir(wrk_dir) ;
}

// -------------------------------------------------------
// Returns the profile of this pulsar
Zaki::Vector::DataSet* Pulsar::GetProfile()
{
  return &profile ;
}

// -------------------------------------------------------
// Returns the metric exponent (nu) in g_tt = exp(2*nu)
Zaki::Vector::DataColumn Pulsar::GetMetricNu() const
{
  return profile[6] ;
}

// -------------------------------------------------------
// Returns the metric exponent (nu) in g_tt = exp(2*nu)
// as a function of radius in km
double Pulsar::GetMetricNu(const double& in_r) const
{ 
  // profile.Interpolate(0, 6) ;

  // Outside the NS
  if (in_r > profile[0].Max())
  {
    return 0.5*log(1 - 2*mp.val * Zaki::Physics::SUN_M_KM/in_r) ;
  }

  int r_idx = profile[0].GetClosestIdx(in_r) ;

  // if (in_r < profile[0].Min())
  // {
  //   Z_LOG_WARNING("The input radius is less than R_min, using R_min instead!") ;
  //   r_val = profile[0].Min() ;
  // }

  return profile[6][r_idx] ;
}

// -------------------------------------------------------
// Returns the sequence point that this pulsar belongs to
SeqPoint Pulsar::GetSeqPoint() const
{
  return seq_point ;
}

// -------------------------------------------------------
/// Returns the sequence profile that this pulsar belongs to
const Zaki::Vector::DataSet* Pulsar::GetSeqProfile() const
{
  return &seq_profile ;
}

// -------------------------------------------------------
// In units of 1 per year
double Pulsar::GetBNVSpinDownLimit() 
{
  // This is in units of 1 per s
  bnv_spin_down_limit = (spin_p_dot/spin_p).val ;

  // convert 1/s to 1/yr
  bnv_spin_down_limit *= 3600*24*365 ;

  bnv_spin_down_limit /= eta_I ;


  return bnv_spin_down_limit ;
}

// -------------------------------------------------------
bool Pulsar::IsBaryon(const std::string& in_label)
{
  if(baryons.find(in_label) != baryons.end())
    return true ; 
  else
    return false ;
}

// -------------------------------------------------------
std::vector<Baryon_Lim> Pulsar::EvalBaryonNumber()
{

  Zaki::Vector::DataColumn r    = profile[0] ;
  Zaki::Vector::DataColumn M_r  = profile[1] * Zaki::Physics::SUN_M_KM ;
  Zaki::Vector::DataColumn nu   = profile[6] ;
  Zaki::Vector::DataColumn rho  = profile[5] ;

  std::vector<Baryon_Lim> B_fraction ;
  B_fraction.reserve(10) ;

  for (size_t i = 7; i < profile.Dim().size() ; i++)
  {
    if ( IsBaryon(profile[i].label) )
    {
      Zaki::Vector::DataColumn b_den = profile[i] * rho ;

      Zaki::Vector::DataColumn integ_fr = 4*M_PI*b_den ;
      integ_fr *= r.pow(2) * (1. - 2*M_r / r ).pow(-0.5) ;

      Zaki::Vector::DataColumn integ_b_dot = integ_fr ;
      integ_b_dot *= exp( nu ) ;

      Zaki::Vector::DataSet integrand({r, integ_fr, integ_b_dot}) ;

      integrand.Interpolate(0, 1) ;
      // integrand.SetIntegrationRelErr(1e-6) ;
      double fr_result = integrand.Integrate(1, {r[0], r[-1]}) ;
      
      integrand.Interpolate(0, 2) ;
      double b_dot_result = integrand.Integrate(2, {r[0], r[-1]}) ;

      //  Correcting for the unit conversion (km/fm)^3 = 1e54 :
      // std::cout << i << ") " << baryons.at(profile[i].label) << " (" 
      //           << profile[i].label <<") = " << fr_result.val * 1e54
      //           << " \u00B1 "<< fr_result.err * 1e54 << "\n" ; 
      B_fraction.emplace_back(  baryons.at(profile[i].label),
                                fr_result * 1e54,
                                b_dot_result * 1e54) ;
    }
  }

  return B_fraction ;
}
//--------------------------------------------------------------
void Pulsar::FindBNVGammaLimits()
{
  std::vector<Baryon_Lim> result ;

  result = EvalBaryonNumber() ;

  std::cout << "\n "<< Zaki::String::Multiply("-", 75) <<"\n" ;
  std::cout << "\n \tBNVSpinDownLimit = "<< GetBNVSpinDownLimit() <<"\n" ;

  std::cout << "\n "<< Zaki::String::Multiply("-", 75) <<"\n" ;
  for (size_t i = 0; i < result.size(); i++)
  {
    result[i].fraction /= seq_point.b ;
    result[i].b_dot    /= seq_point.b ;
    result[i].limit = GetBNVSpinDownLimit() / result[i].b_dot ;

    std::cout << "  " << result[i] << "\n" ;
  }
  std::cout << "\n " << Zaki::String::Multiply("-", 75) << "\n" ;
}
//--------------------------------------------------------------
// Plots the pulsar's relative composition vs radius
void Pulsar::PlotRelativeComposition(const std::string& file, 
                          const Zaki::String::Directory& in_dir)
{
  std::vector<std::pair<int, std::string>> tmp_plt_idx ;
  for (size_t i = 7; i < profile.Dim().size() ; i++)
  {
    // Making the plot more legible by removing unimportant contributions
    if ( profile[i].Max() < 1e-4 )
    {
      continue ;
    }

    tmp_plt_idx.emplace_back(i, CompOSE_EOS::Compose_Dict.at( std::stoi( profile[i].label.c_str() ) ).name ) ;
  }

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetGrid() ;
  // plt_par.SetXAxis({0, profile[0][-1]*1.02}) ;
  plt_par.SetXAxis({0, 14}) ;
  plt_par.SetYAxis({1e-3, 1}) ;
  plt_par.SetXAxisLabel("$R \\, ( \\, km \\, )$") ;
  plt_par.SetYAxisLabel("$f_i$") ;
  plt_par.SetLegend({"upper right", 1.03, 1.0}) ;
  
  profile.SetPlotPars(plt_par) ;
  profile.SetWrkDir(wrk_dir) ;
  profile.SemiLogYPlot(0, tmp_plt_idx, in_dir + file, name ) ;
}

//--------------------------------------------------------------
// Plots the pulsar's composition vs radius
void Pulsar::PlotAbsoluteComposition(const std::string& file, 
                          const Zaki::String::Directory& in_dir) const
{
  Zaki::Vector::DataSet comp_ds({ profile[0] }) ;

  std::vector<std::pair<int, std::string>> tmp_plt_idx ;
  for (size_t i = 7 ; i < profile.Dim().size() ; i++)
  {
    // Making the plot more legible by removing unimportant contributions
    if ( profile[i].Max() < 1e-4 )
    {
      continue ;
    }

    tmp_plt_idx.emplace_back(comp_ds.Dim().size(), CompOSE_EOS::Compose_Dict.at( std::stoi( profile[i].label.c_str() ) ).name ) ;
    comp_ds.data_set.emplace_back( profile[i] * profile[5]) ;
  }

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetGrid() ;
  // plt_par.SetXAxis({0, profile[0][-1]*1.02}) ;
  plt_par.SetXAxis({0, 14}) ;
  plt_par.SetYAxis({1e-3, 1}) ;
  plt_par.SetXAxisLabel("$R \\, [ \\, km \\, ]$") ;
  plt_par.SetYAxisLabel("$n\\, [ \\, {\\rm fm}^{-3} \\, ]$") ;
  plt_par.SetLegend({"upper right", 1.03, 1.0}) ;
  
  comp_ds.SetPlotPars(plt_par) ;
  comp_ds.SetWrkDir(wrk_dir) ;
  comp_ds.SemiLogYPlot(0, tmp_plt_idx, in_dir + file, name ) ;
}

//--------------------------------------------------------------
// Plots the Fermi energies vs radius
void Pulsar::PlotFermiE(const std::string& file, 
                          const Zaki::String::Directory& in_dir) const
{
  // ------------------------------------
  //        Finding Fermi Energy
  // ------------------------------------
  Zaki::Vector::DataColumn fermi_electron =  ( 
    pow(Zaki::Physics::ELECTRON_M_FM, 2) 
    + (3*M_PI*M_PI* profile["0"] * profile[5]).pow(2./3.)).sqrt() / Zaki::Physics::MEV_2_INV_FM ;

  Zaki::Vector::DataColumn fermi_muon =  ( 
    pow(Zaki::Physics::MUON_M_FM, 2) 
    + (3*M_PI*M_PI* profile["1"] * profile[5]).pow(2./3.)).sqrt() / Zaki::Physics::MEV_2_INV_FM ;

  Zaki::Vector::DataSet fermi_ds({profile[5], 
                                  fermi_electron, fermi_muon}) ;

  // ------------------------------------
  //              Plotting
  // ------------------------------------
  Zaki::Vector::DataSet::PlotParam plt_par ;

  plt_par.SetGrid() ;
  plt_par.SetXAxisLabel("$R\\,\\, ( {\\rm km} )$") ;
  plt_par.SetYAxisLabel("$E_F\\,\\, ( {\\rm MeV} )$") ;
  plt_par.SetXAxis({0, 14}) ;
  plt_par.SetLegend({"upper right", 0.0, 1}) ;

  fermi_ds.SetPlotPars(plt_par) ;

  // Zaki::String::Directory file_dir = in_dir ;
  // if (in_dir_type == Dir_Type::relative)
  // {
  //   file_dir = wrk_dir + "/" + in_dir ;
  // }

  // std::cout << "\n\t" << file_dir + "/E_Fermi.pdf\n" ;
  fermi_ds.SetWrkDir(wrk_dir) ;
  fermi_ds.SemiLogXPlot(0, {{1, "$e^-$"}, {2, "$\\mu^-$"}}, 
                           in_dir + file, name) ;


  // // Electron particle code is '0'
  // Zaki::Vector::DataSet fermi_ds({ profile[0], 
  //  ( pow(Zaki::Physics::ELECTRON_M_FM, 2) 
  //   + (3*M_PI*M_PI* profile["0"] * profile[5]).pow(2./3.)).sqrt() / Zaki::Physics::MEV_2_INV_FM 
  //                               }) ;
  // Zaki::Vector::DataSet::PlotParam plt_par ;
  // plt_par.SetGrid() ;
  // // plt_par.SetXAxis({0, profile[0][-1]*1.02}) ;
  // plt_par.SetXAxis({0, 14}) ;
  // // plt_par.SetYAxis({1e-3, 1}) ;
  // plt_par.SetXAxisLabel("$R \\, ( \\, km \\, )$") ;
  // plt_par.SetYAxisLabel("$E_F\\, ( \\, {\\rm MeV} \\, )$") ;
  // // plt_par.SetLegend({"upper right", 1.03, 1.0}) ;
  
  // fermi_ds.SetPlotPars(plt_par) ;
  // fermi_ds.SetWrkDir(wrk_dir) ;
  // fermi_ds.Plot(0, 1, in_dir + file, name ) ;
}

//--------------------------------------------------------------
// Plots the metric function exp[ nu(r) ] vs radius
void Pulsar::PlotExpNu(const std::string& file, 
               const Zaki::String::Directory& in_dir) const
{
  Zaki::Vector::DataSet exp_nu_r({ profile[0], exp(profile[6]) }) ; 
  exp_nu_r.SetWrkDir(wrk_dir) ;

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetGrid() ;
  // plt_par.SetXAxis({0, profile[0][-1]*1.02}) ;
  plt_par.SetXAxis({0, 14}) ;
  plt_par.SetYAxisLabel("$\\exp ( \\nu [r] )$") ;

  exp_nu_r.SetPlotPars(plt_par) ;
  exp_nu_r.Plot(0, 1, in_dir + "/" + file, name) ;
}

//--------------------------------------------------------------
// The surface gravity g = GM e^{−nu(R)} / R^2 
//  in units of 10^{14} cm s^{−2}
double Pulsar::Get_Surface_Gravity() const
{
  // Surface radius
  double R_s = profile[0].Max() ;

  double g_14 = mp.val * Zaki::Physics::SUN_M_KM ;
        g_14 *= exp(-GetMetricNu(R_s) );
        g_14 /= pow(R_s, 2) ;

  // Dividing by 10^{14} cm s^{−2}
  g_14 /= 0.01112650056053619 ;

  return g_14 ;
}

//--------------------------------------------------------------
// Sets pulsar's core temperature in kelvin
// Not red-shifted: as measured in a local frame inside the star
void Pulsar::Set_T_Core(const double& in_T_core) 
{
  T_core = in_T_core ;

  double R_0 = profile[0].Min() ;

  T_blanket = T_core * exp( GetMetricNu(R_0) - GetMetricNu(r_blanket) );

  double g_14 = Get_Surface_Gravity() ;

  // T_surf = 34.9947 * pow(g_14, 0.25) * pow(T_blanket, 50.0/91.0) ;

  // From Eq. (49) of [ arXiv: 0502116 ]
  T_surf = pow(1e24 * g_14 * pow(1.81 * T_blanket / 1e8, 2.42), 0.25) ;
}

//--------------------------------------------------------------
// Sets pulsar's blanket (outer crust) temperature in kelvin
// Not red-shifted: as measured in a local frame inside the star
void Pulsar::Set_T_Blanket(const double& in_T_blanket) 
{
  T_blanket = in_T_blanket ;

  double g_14 = Get_Surface_Gravity() ;

  // T_surf = 34.9947 * pow(g_14, 0.25) * pow(T_blanket, 50.0/91.0) ;

  // From Eq. (49) of [ arXiv: 0502116 ]
  T_surf = pow(1e24 * g_14 * pow(1.81 * T_blanket / 1e8, 2.42), 0.25) ;

  double R_0 = profile[0].Min() ;

  T_core = T_blanket * exp( GetMetricNu(r_blanket) - GetMetricNu(R_0) );
}

//--------------------------------------------------------------
// Sets pulsar's surface temperature in kelvin
// Not red-shifted: as measured in a local frame inside the star
void Pulsar::Set_T_Surf(const double& in_T_surf) 
{
  T_surf = in_T_surf ;

  double Ts6 = T_surf / 1e6 ;
  double g_14 = Get_Surface_Gravity() ;

  T_blanket = 1.288e8 * pow( pow(Ts6, 4) / g_14, 0.455) ;

  double R_0 = profile[0].Min() ;

  T_core = T_blanket * exp( GetMetricNu(r_blanket) - GetMetricNu(R_0) );
}

//--------------------------------------------------------------
// Sets pulsar's apparent (red-shifted) surface temperature (T_s) in kelvin
// as detected by a distant observer
void Pulsar::Set_T_Surf_Apparent(const double& in_T_apparent) 
{
  // Value of nu at the surface of NS
  double nu_R = profile[6][-1] ;

  Set_T_Surf( exp(-nu_R) * in_T_apparent) ;
}

//--------------------------------------------------------------
// Reurns the blanket radius
// defined by eps ~ 10^{10} [ g / cm^3 ]
double Pulsar::Get_R_Blanket() const 
{
  return r_blanket ;
}

//--------------------------------------------------------------
// Reurns the threshold radius for direct Urca
// defined by p_f(n) = p_F(p) + p_F(e)
// This value is set in 'FindProfile'.
double Pulsar::Get_R_Durca() const 
{
  return r_durca_thresh ;
}

//--------------------------------------------------------------
// Returns pulsar's core temperature in kelvin
// Not red-shifted: as measured in a local frame inside the star
double Pulsar::Get_T_Core() const 
{
  return T_core ;
}

//--------------------------------------------------------------
// Returns pulsar's blanket (outer crust) temperature in kelvin
// Not red-shifted: as measured in a local frame inside the star
double Pulsar::Get_T_Blanket() const 
{
  return T_blanket ;

}

//--------------------------------------------------------------
// Returns pulsar's surface temperature in kelvin
// Not red-shifted: as measured in a local frame inside the star
double Pulsar::Get_T_Surf() const 
{
  return T_surf ;
}

//--------------------------------------------------------------
// Returns the apparent (red-shifted) surface temperature (T-infty)
// as detected by a distant observer.
double Pulsar::Get_T_Surf_Apparent() const
{
  // Value of nu at the surface of NS
  double nu_R = profile[6][-1] ;

  return exp(nu_R) * T_surf ;
}

//--------------------------------------------------------------
// Returns the red-shifted temperature (T-infty)
// as detected by a distant observer.
// Constant: assuming that the interior of NS is thermalized
double Pulsar::Get_Redshifted_T() const
{ 
  // Value of nu at r_blanket
  // double nu_r_blanket = GetMetricNu(r_blanket) ; // non-const
  double nu_r_blanket = profile[6][r_blanket_idx] ; // const

  return exp(nu_r_blanket) * T_blanket ;
}

//--------------------------------------------------------------
// Returns the temperature as a function of r(km)
// as measured in a local frame inside the star.
// Assuming that the interior of NS is thermalized
double Pulsar::Get_T(const double& in_r) const
{
  if(in_r > profile[0].Max())
  {
    return 0 ;
  }
  
  if(in_r < profile[0].Min())
  {
    return exp(-profile[6][0]) * Get_Redshifted_T() ;
  }

  int r_idx = profile[0].GetClosestIdx(in_r) ;

  return exp(-profile[6][r_idx]) * Get_Redshifted_T() ;
}

//--------------------------------------------------------------
// Returns the temperature DataColumn as a function of r(km)
// as measured in a local frame inside the star.
// Assuming that the interior of NS is thermalized
Zaki::Vector::DataColumn Pulsar::Get_T() const
{
  return exp(- profile[6]) * Get_Redshifted_T() ;
}

//--------------------------------------------------------------
// Reurns the red-shifted surface emission luminosity L^{infty} 
// (in erg s^{-1})
// due to black body radiation
double Pulsar::Get_Surface_Photon_Lumin()  const
{
  // Stefan-Boltzmann constant = 0.567 erg s^{−1} m^{−2} K^{−4} 
  double sigB = 567037.4419184428 ; // erg s^{−1} km^{−2} K^{−4}

  // in km
  double r_s =  profile[0].Max() ;

  double L_gamma_infty  = 4 * M_PI * pow(r_s, 2) * sigB ;
         L_gamma_infty *=  pow(Get_T_Surf(), 4) ;
         L_gamma_infty *=  exp(2*profile[6][-1]) ;

  return  L_gamma_infty ;
}
//--------------------------------------------------------------
// Reurns the red-shifted surface emission luminosity L^{infty}
// due to black body radiation in units of the red-shifted T^inf: 
// it has to be multiplied by "[ (T^inf_9)^2.42 ]"
// Units: [erg s^-1]
double Pulsar::Get_Surface_Photon_Lumin_T9242() const 
{
  // Stefan-Boltzmann constant = 0.567 erg s^{−1} m^{−2} K^{−4} 
  double sigB = 567037.4419184428 ; // erg s^{−1} km^{−2} K^{−4}

  // in km
  double r_s =  profile[0].Max() ;

  double g_14 = Get_Surface_Gravity() ;

  // surface temperature in units of [ T^infty_9 ]^(2.42)
  // Using Eq. (49) of Ref. [ arXiv: 0502116 ]
  double t_surf_4 = 1e24 ;
        t_surf_4 *= g_14 ;
        t_surf_4 *= pow(18.1 * exp(-profile[6][r_blanket_idx]), 2.42) ;

  double L_gamma_infty  = 4 * M_PI * pow(r_s, 2) * sigB ;
         L_gamma_infty *=  t_surf_4 ;
         L_gamma_infty *=  exp(2*profile[6][-1]) ;

  return  L_gamma_infty ;
}

//--------------------------------------------------------------
// Reurns the neutrino emissivity
// due to DUrca
// Units: [erg cm^-3 s^-1]
double Pulsar::Get_DUrca_Neutrino_Emissivity(const double& in_radius) const
{
  if (in_radius > r_durca_thresh)
  {
    return 0 ;
  }
  
  double T_9 = Get_T(in_radius) / 1e9 ;

  return DUrca_emissivity_prefactor * pow(T_9, 6) ;
}

//--------------------------------------------------------------
// Reurns the neutrino Durca rate defined by
// emissivity divided by [ k.T ]
// Units: [fm^-3 s^-1]
double Pulsar::Get_DUrca_Neutrino_Rate(const double& in_radius) const
{
  // Change [erg cm^-3 s^-1] to [ MeV fm^-3 s^-1]  why??
  // return 624150.64799632 * 1e-39
  //         * Get_DUrca_Neutrino_Emissivity(in_radius)
  //                        / Get_T(in_radius) ;

  // Change [erg K^-1 cm^-3 s^-1] to [ fm^-3 s^-1]
  return 1.380649e-16 * 1e-39
          * Get_DUrca_Neutrino_Emissivity(in_radius)
                         / Get_T(in_radius) ;
}

//--------------------------------------------------------------
// Reurns the neutrino Durca rate defined by
// emissivity divided by [ k.T ]
// Units: [fm^-3 s^-1]
Zaki::Vector::DataColumn Pulsar::Get_DUrca_Neutrino_Rate() const
{
  // // Change [erg cm^-3 s^-1] to [ MeV fm^-3 s^-1]
  // return 624150.64799632 * 1e-39 
  //         * Get_DUrca_Neutrino_Emissivity()
  //                        / Get_T() ;

  // Change [erg K^-1 cm^-3 s^-1] to [ fm^-3 s^-1]
  return 1.380649e-16 * 1e-39 
          * Get_DUrca_Neutrino_Emissivity()
                         / Get_T() ;
}

//--------------------------------------------------------------
// Reurns the neutrino emissivity DataColumn
// due to DUrca
// Units: [erg cm^-3 s^-1]
Zaki::Vector::DataColumn Pulsar::Get_DUrca_Neutrino_Emissivity() const
{ 
  return DUrca_emissivity_prefactor * r_durca_cond * Get_T().pow(6) * 1e-54 ;
}

//--------------------------------------------------------------
// Reurns the neutrino emissivity
// due to DUrca, in units of the red-shifted T^inf: 
// it has to be multiplied by "[ (T^inf_9)^6 ]"
// Units: [erg cm^-3 s^-1]
Zaki::Vector::DataColumn Pulsar::Get_DUrca_Neutrino_Emissivity_T96() const 
{
  return DUrca_emissivity_prefactor * r_durca_cond * exp(-6*profile[6]) ;
}

//--------------------------------------------------------------
// Reurns the neutrino emissivity
// due to MUrca
// Units: [erg cm^-3 s^-1]
double Pulsar::Get_MUrca_Neutrino_Emissivity(const double& in_radius) const
{
  double T_9 = Get_T(in_radius) / 1e9 ;

  return MUrca_emissivity_prefactor * pow(T_9, 8) ;
}

//--------------------------------------------------------------
// Reurns the neutrino emissivity
// due to MUrca
// Units: [erg cm^-3 s^-1]
Zaki::Vector::DataColumn Pulsar::Get_MUrca_Neutrino_Emissivity() const
{
  return MUrca_emissivity_prefactor * Get_T().pow(8) * 1e-72 ;
}

//--------------------------------------------------------------
// Reurns the neutrino emissivity
// due to DUrca, in units of the red-shifted T^inf: 
// it has to be multiplied by "[ (T^inf_9)^8 ]"
// Units: [erg cm^-3 s^-1]
Zaki::Vector::DataColumn Pulsar::Get_MUrca_Neutrino_Emissivity_T98() const 
{
  return MUrca_emissivity_prefactor * exp(-8*profile[6]) ;
}

//--------------------------------------------------------------
// Reurns the neutrino Murca rate defined by
// emissivity divided by [ k.T ]
// Units: [fm^-3 s^-1]
double Pulsar::Get_MUrca_Neutrino_Rate(const double& in_radius) const
{
  // Change [erg cm^-3 s^-1] to [ MeV fm^-3 s^-1]
  // return 624150.64799632 * 1e-39 
  //         * Get_MUrca_Neutrino_Emissivity(in_radius)
  //                        / Get_T(in_radius) ;
                         
  // Change [erg K^-1 cm^-3 s^-1] to [ fm^-3 s^-1]
  return 1.380649e-16 * 1e-39 
          * Get_MUrca_Neutrino_Emissivity(in_radius)
                         / Get_T(in_radius) ;                         
}

//--------------------------------------------------------------
// Reurns the neutrino Murca rate defined by
// emissivity divided by [ k.T ]
// Units: [fm^-3 s^-1]
Zaki::Vector::DataColumn Pulsar::Get_MUrca_Neutrino_Rate() const
{
  // // Change [erg cm^-3 s^-1] to [ MeV fm^-3 s^-1]
  // return 624150.64799632 * 1e-39 
  //         * Get_MUrca_Neutrino_Emissivity()
  //                        / Get_T() ;

  // Change [erg K^-1 cm^-3 s^-1] to [ fm^-3 s^-1]
  return 1.380649e-16 * 1e-39 
          * Get_MUrca_Neutrino_Emissivity()
                         / Get_T() ;        
}

//--------------------------------------------------------------
// Reurns the neutrino emissivity
// due to DUrca and MUrca
// Units: [erg cm^-3 s^-1]
double Pulsar::Get_Neutrino_Emissivity(const double& in_radius) const
{
  return Get_DUrca_Neutrino_Emissivity(in_radius) 
        + Get_MUrca_Neutrino_Emissivity(in_radius) ;
}

//--------------------------------------------------------------
// Reurns the neutrino emissivity
// due to DUrca and MUrca
// Units: [erg cm^-3 s^-1]
Zaki::Vector::DataColumn Pulsar::Get_Neutrino_Emissivity() const
{
  return Get_DUrca_Neutrino_Emissivity() 
        + Get_MUrca_Neutrino_Emissivity() ;
}

//--------------------------------------------------------------
// Reurns the total neutrino rate defined by
// emissivity divided by [ k.T ]
// Units: [fm^-3 s^-1]
double Pulsar::Get_Neutrino_Rate(const double& in_radius) const
{
  // // Change [erg cm^-3 s^-1] to [ MeV fm^-3 s^-1]
  // return 624150.64799632 * 1e-39 
  //         * Get_Neutrino_Emissivity(in_radius)
  //                        / Get_T(in_radius) ;

  // Change [erg K^-1 cm^-3 s^-1] to [ fm^-3 s^-1]
  return 1.380649e-16 * 1e-39 
          * Get_Neutrino_Emissivity(in_radius)
                         / Get_T(in_radius) ;    
}

//--------------------------------------------------------------
// Reurns total neutrino rate defined by
// emissivity divided by [ k.T ]
// Units: [fm^-3 s^-1]
Zaki::Vector::DataColumn Pulsar::Get_Neutrino_Rate() const
{
  // Change [erg cm^-3 s^-1] to [ MeV fm^-3 s^-1]
  // return 624150.64799632 * 1e-39 
  //         * Get_Neutrino_Emissivity()
  //                        / Get_T() ;

  // Change [erg K^-1 cm^-3 s^-1] to [ fm^-3 s^-1]
  return 1.380649e-16 * 1e-39 
          * Get_Neutrino_Emissivity()
                         / Get_T() ;    
}

//--------------------------------------------------------------
// Plots the neutrino emissivity as function of radius
// due to DUrca and MUrca
// Units: [erg cm^-3 s^-1]
void Pulsar::Plot_Neutrino_Emissivity(const std::string& file, 
                          const Zaki::String::Directory& in_dir) const
{
  Zaki::Vector::DataSet nu_emiss_plt({profile[0], Get_DUrca_Neutrino_Emissivity(),
  Get_MUrca_Neutrino_Emissivity(), Get_Neutrino_Emissivity()}) ;

  nu_emiss_plt[1].label = "DUrca" ;
  nu_emiss_plt[2].label = "MUrca" ;
  nu_emiss_plt[3].label = "Tot" ;

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetGrid() ;
  plt_par.SetXAxis({0, 14}) ;
  // plt_par.SetYAxis({1e-3, 1}) ;
  plt_par.SetXAxisLabel("$R \\, ( \\, km \\, )$") ;
  plt_par.SetYAxisLabel("$Q_{\\nu}\\, ( \\, {\\rm erg} {\\rm cm}^{-3} {\\rm s}^{-1}\\, )$") ;
  plt_par.SetLegend({"upper right", 1.03, 1.0}) ;
  
  nu_emiss_plt.SetPlotPars(plt_par) ;
  nu_emiss_plt.SetWrkDir(wrk_dir) ;
  nu_emiss_plt.SemiLogYPlot(0, {{1, {{"ls", "--"}, {"label", "DUrca"}}}, 
                                {2, {{"ls", "--"}, {"label","MUrca"}}}, 
                                {3, {{"ls", "-."}, {"label","Tot"}}}}, in_dir + file, name ) ;
}

//--------------------------------------------------------------
// Plots the neutrino rate as function of radius
// due to DUrca and MUrca
// Units: [ fm^-3 s^-1]
void Pulsar::Plot_Neutrino_Rate(const std::string& file, 
                          const Zaki::String::Directory& in_dir) const
{
  Zaki::Vector::DataSet nu_rate_plt({profile[0], Get_DUrca_Neutrino_Rate(),
  Get_MUrca_Neutrino_Rate(), Get_Neutrino_Rate()}) ;

  nu_rate_plt[1].label = "DUrca" ;
  nu_rate_plt[2].label = "MUrca" ;
  nu_rate_plt[3].label = "Tot" ;

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetGrid() ;
  plt_par.SetXAxis({0, 14}) ;
  // plt_par.SetYAxis({1e-3, 1}) ;
  plt_par.SetXAxisLabel("$R \\, ( \\, km \\, )$") ;
  plt_par.SetYAxisLabel("$\\Gamma_{\\nu}\\, ( \\,  {\\rm s}^{-1} {\\rm fm}^{-3}\\, )$") ;
  plt_par.SetLegend({"upper right", 1.03, 1.0}) ;
  
  nu_rate_plt.SetPlotPars(plt_par) ;
  nu_rate_plt.SetWrkDir(wrk_dir) ;
  nu_rate_plt.SemiLogYPlot(0, {{1, {{"ls", "--"}, {"label", "DUrca"}}}, 
                                {2, {{"ls", "--"}, {"label","MUrca"}}}, 
                                {3, {{"ls", "-."}, {"label","Tot"}}}}, in_dir + file, name ) ;
}

//--------------------------------------------------------------
// Reurns the neutrino luminosity
// due to DUrca and MUrca
// Units: [erg s^-1]
double Pulsar::Get_Neutrino_Lumin() const
{
  Z_LOG_INFO("Calculating neutrino luminosity for '" + name + "'.") ;

  Zaki::Vector::DataSet nu_lumin_ds(profile[0]) ;

  Zaki::Vector::DataColumn r = profile[0]  ;
  Zaki::Vector::DataColumn M_r  = profile[1] * Zaki::Physics::SUN_M_KM ; 
  Zaki::Vector::DataColumn nu  = profile[6] ; 

  // change unit to [ erg km^-3 s^-1 ]
  Zaki::Vector::DataColumn Q_nu = 1e15 * Get_Neutrino_Emissivity() ; 

  // for (size_t i = 0; i < profile[0].Size(); i++)
  // {
  //   // change unit to [ erg km^-3 s^-1 ]
  //   Q_nu.vals.emplace_back( 1e15 * Get_Neutrino_Emissivity(r[i]) ) ;
  // }

  Zaki::Vector::DataColumn nu_lumin_dc = Q_nu ; 

    nu_lumin_dc *= 4 * M_PI * r.pow(2) ;

    // Volume correction factor
    nu_lumin_dc *= (1. - 2* M_r / r).pow(-0.5) ;

    // Red-shift factor
    nu_lumin_dc *= exp(2*nu) ;

  nu_lumin_ds.data_set.emplace_back(nu_lumin_dc) ;

  // nu_lumin.AddColumn("nu_emiss") ;
  // for (size_t i = 0; i < profile[0].Size(); i++)
  // {
  //   double tmp_r = profile[0][i] ;
  //   double M_r  = profile[1][i] * Zaki::Physics::SUN_M_KM ; 

  //   // change unit to [ erg km^-3 s^-1 ]
  //   nu_lumin[1][i] = 1e15 * Get_Neutrino_Emissivity(tmp_r) ;

  //   nu_lumin[1][i] *= 4 * M_PI * pow(tmp_r, 2) ;

  //   // Volume correction factor
  //   nu_lumin[1][i] *= pow(1. - 2* M_r / tmp_r, -0.5) ;

  //   // Red-shift factor
  //   nu_lumin[1][i] *= exp(2*GetMetricNu(tmp_r)) ;
  // }
  
  nu_lumin_ds.Interpolate(0,1) ;

  return nu_lumin_ds.Integrate(1, {r.Min(), r.Max()}) ;
}

//--------------------------------------------------------------
// Reurns the integrated neutrino rate
// due to DUrca in units of the red-shifted T^inf: 
//  it has to be multiplied by "[ (T^inf_9)^5 ]"
// Units: [s^-1]
double Pulsar::Get_DUrca_Neutrino_Rate_T95() const 
{
  Zaki::Vector::DataColumn r = profile[0]  ;
  Zaki::Vector::DataColumn M_r  = profile[1] * Zaki::Physics::SUN_M_KM ; 
  Zaki::Vector::DataColumn nu  = profile[6] ; 

  // durca_ds[0] : Radius (km)
  // durca_ds[1] : 4e27 erg.s^-1.cm^-3
  Zaki::Vector::DataSet durca_ds({r, DUrca_emissivity_prefactor * r_durca_cond }) ;

  // Dividing by [ kB (erg/K) * T(K) ] unit becomes [ s^-1.cm^-3 ]
  durca_ds[1] /= 1e9 * 1.380649e-16 ;

  // Change unit to [ s^-1.km^-3 ]
  durca_ds[1] *= 1e15  ;

  // Changing T_9^5 to (T^inf_9)^5
  durca_ds[1] *= exp(-5*nu)  ;

  // Spherical volume element
  durca_ds[1] *= 4 * M_PI * r.pow(2) ;

  // Volume correction factor
  durca_ds[1] *= (1. - 2* M_r / r).pow(-0.5) ;

  // Red-shift factor
  durca_ds[1] *= exp(nu) ;

  durca_ds.Interpolate(0, 1) ;

  return  durca_ds.Integrate(1, {r.Min(), r.Max()}) ;
}

//--------------------------------------------------------------
// Reurns the integrated neutrino rate
// due to MUrca in units of the red-shifted T^inf: 
//  it has to be multiplied by "[ (T^inf_9)^7 ]"
// Units: [s^-1]
double Pulsar::Get_MUrca_Neutrino_Rate_T97() const 
{
  Zaki::Vector::DataColumn r = profile[0]  ;
  Zaki::Vector::DataColumn M_r  = profile[1] * Zaki::Physics::SUN_M_KM ; 
  Zaki::Vector::DataColumn nu  = profile[6] ; 

  // murca_ds[0] : Radius (km)
  // murca_ds[1] : 3e21 erg.s^-1.cm^-3
  Zaki::Vector::DataSet murca_ds({r}) ;
  murca_ds.AddColumn("MUrca", MUrca_emissivity_prefactor) ;

  // Dividing by [ kB (erg/K) * T(K) ] unit becomes [ s^-1.cm^-3 ]
  murca_ds[1] /= 1e9 * 1.380649e-16 ;

  // Change unit to [ s^-1.km^-3 ]
  murca_ds[1] *= 1e15  ;

  // Changing T_9^7 to (T^inf_9)^7
  murca_ds[1] *= exp(-7*nu)  ;

  // Spherical volume element
  murca_ds[1] *= 4 * M_PI * r.pow(2) ;

  // Volume correction factor
  murca_ds[1] *= (1. - 2* M_r / r).pow(-0.5) ;

  // Red-shift factor
  murca_ds[1] *= exp(nu) ;

  murca_ds.Interpolate(0, 1) ;

  return  murca_ds.Integrate(1, {r.Min(), r.Max()}) ;
}

//--------------------------------------------------------------
} // End of namespace CompactStar
//--------------------------------------------------------------

//==============================================================
