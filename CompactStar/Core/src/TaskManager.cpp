/*
  TaskManager class
*/
#include <cmath>

#include <Confind/ContourFinder.hpp>

#include <Zaki/Math/Math_Core.hpp>

#include "CompactStar/Core/TaskManager.hpp"
#include "CompactStar/MixedStar/DarkCore_Analysis.hpp"
#include "CompactStar/Core/TOVSolver_Thread.hpp"
#include <CompactStar/EOS/Fermi_Gas.hpp>

using namespace CompactStar ;

Zaki::Math::Quantity J0348_p_0432_M(2.01, 0.04) ;
//==============================================================
// Condition for exporting the profiles
bool Mass_Condition(const CompactStar::MixedStar& in_star)
{
  return ((in_star.GetSequence().d.m 
          + in_star.GetSequence().v.m) >=
          (J0348_p_0432_M.val - J0348_p_0432_M.err) ) &&
          (in_star.GetSequence().d.m 
          + in_star.GetSequence().v.m) <=
          (J0348_p_0432_M.val + J0348_p_0432_M.err) ;
}

//==============================================================
bool TrueCondition(const CompactStar::MixedStar& in_star)
{
  return true ;
}

//==============================================================
//                        TaskManager class
//==============================================================
//--------------------------------------------------------------
// Constructor
TaskManager::TaskManager(const size_t& in_num_thrds) 
  : Prog("TaskManager"), 
    num_of_thrds(in_num_thrds), threads(in_num_thrds),
    critical_curve("Critical curve")
{

  t_start = std::chrono::high_resolution_clock::now() ;

  if( in_num_thrds > std::thread::hardware_concurrency() )
  {
    num_of_thrds = std::thread::hardware_concurrency() ;
    Z_LOG_WARNING("The number of threads can't be "
                  "larger than '" 
                  + std::to_string(
                    std::thread::hardware_concurrency()) 
                  +"'.") ;
  }
}
//--------------------------------------------------------------
TaskManager::~TaskManager()
{
  auto t_end = std::chrono::high_resolution_clock::now();

  double elapsed_time_s = std::chrono::duration<double>
                            (t_end-t_start).count();

  std::cout << "Task manager took " << (int) (elapsed_time_s / 60) 
            << " mins, " << std::fmod(elapsed_time_s, 60) 
            << " secs.\n" ;
}

//--------------------------------------------------------------
void TaskManager::SetDarEOSDir(const Zaki::String::Directory& in_dir) 
{
  dar_eos_dir = in_dir ;
}

//--------------------------------------------------------------
void TaskManager::SetVisEOSDir(const Zaki::String::Directory& in_dir) 
{
  vis_eos_dir = in_dir ;
}

//--------------------------------------------------------------
TaskManager* TaskManager::SetExclusionRegion(
  const Zaki::Math::Cond_Polygon& in_c_poly) 
{
  c_poly = in_c_poly ;
  return this ;
}

//--------------------------------------------------------------
TaskManager* TaskManager::SetGrid( const Zaki::Math::Axis& in_v_ax, 
              const Zaki::Math::Axis& in_d_ax) 
{
  v_ax = in_v_ax ;
  d_ax = in_d_ax ;
  return this ;
}

//--------------------------------------------------------------
// TaskManager* TaskManager::SetWrkDir(const Zaki::String::Directory& in_dir) 
// {
//   wrk_dir = in_dir ;
//   return this ;
// }

//--------------------------------------------------------------
void TaskManager::Work() 
{
  Zaki::Util::LogManager::SetLogLevels(Zaki::Util::LogLevel::Warning) ;

  for (size_t i = 0; i < num_of_thrds ; i++)
  {
    threads[i] = std::thread(&TaskManager::Task, this, i+1) ;
  }

  for(auto& t: threads)
  {
    if (t.joinable())
    {
      t.join();
    }
  }

  // Access the sequence info by calling this:
  CompactStar::TOVSolver_Thread::GetSequence() ;

  Zaki::Util::LogManager::SetLogLevels(Zaki::Util::LogLevel::Info) ;
}

//--------------------------------------------------------------
void TaskManager::Task(const int tsk_id) const
{
  CompactStar::TOVSolver_Thread solver(tsk_id) ;
  solver.SetRadialRes(3.0e4) ;
  solver.SetWrkDir(wrk_dir) ;

  // solver.ImportEOS_Dark("EOS/Fermi_Gas_0.8mn.eos") ;
  // solver.ImportEOS_Dark(dar_eos_dir) ;

  // solver.ImportEOS("EOS/CompOSE/DS(CMF)-1_with_crust/DS(CMF)-1_with_crust.eos") ;
  solver.ImportEOS(vis_eos_dir, dar_eos_dir) ;

  // solver.AddMixedCondition(TrueCondition) ;

  const size_t pts_per_thrd = (d_ax.res + 1) / num_of_thrds ; 
  size_t res_i = pts_per_thrd - 1 ;

  if ( tsk_id == num_of_thrds )
    res_i += (d_ax.res + 1) % num_of_thrds ;

  const size_t min_idx_offset = pts_per_thrd*(tsk_id-1) ;
  
  solver.SetMinIdxOffset(min_idx_offset) ;
  solver.SetExclusionRegion(c_poly) ;

  Zaki::Math::Range<double> range_i ;
  if( d_ax.scale == "Log" ){
    double log_step = (log10(d_ax.Max()) - log10(d_ax.Min()) )/ d_ax.res ;

    double log_min = log10(d_ax.Min()) + log_step*pts_per_thrd*(tsk_id-1) ;
    double log_max = log_min + log_step*res_i ;
    range_i = { pow(10,log_min), pow(10,log_max) } ; 
  } 
  else {
    double step = (d_ax.Max() - d_ax.Min())/ d_ax.res ;
    double min  = d_ax.Min() + step*pts_per_thrd*(tsk_id-1) ;
    double max  = min + step*res_i ;
    range_i = { min, max } ; 
  }

  char tmp_1[100] ;
  char tmp_2[100] ;
  snprintf(tmp_1, sizeof(tmp_1), "%.1f_%zux%zu", m_chi, v_ax.res, d_ax.res) ;
  snprintf(tmp_2, sizeof(tmp_2), "%.1f", m_chi) ;
  solver.SetSeqFileName( std::string(tmp_1)) ;
  // std::cout << "\n --> " << tmp_1 << "\n" ;
  solver.Solve_Mixed( v_ax, {range_i, res_i, d_ax.scale},
                      "NStar/Dark_Core/" + std::string(tmp_2)
                      // "testing/No crust/30K_Adaptive"
                      ,  std::string(tmp_1)) ;
}

//--------------------------------------------------------------
void TaskManager::ImportSequence(const Zaki::String::Directory& in_dir) 
{
  sequence_grid.Import(wrk_dir + in_dir) ;
}

//--------------------------------------------------------------
void TaskManager::FindDarkEOS(const double& in_m) 
{
  m_chi = in_m ;

  CompactStar::Fermi_Gas F_gas(in_m*Zaki::Physics::NEUTRON_M_MEV) ;
  F_gas.SetRhoRange({1e-5, 5e+1}) ;

  F_gas.FindEOS(2000) ;
  F_gas.SetWrkDir(wrk_dir) ;

  char tmp[100] ;
  snprintf(tmp, sizeof(tmp), "%.1f", in_m) ;
  dar_eos_dir = "EOS/Fermi_Gas_"+std::string(tmp)+"mn.eos" ;
  F_gas.ExportEOS(dar_eos_dir) ;
}

//--------------------------------------------------------------
void TaskManager::SetChiMass(const double& in_m) 
{
  m_chi = in_m ;

  char tmp_chi[50] ;
  snprintf(tmp_chi, sizeof(tmp_chi), "%.1f", m_chi) ;
  chi_str = tmp_chi ;
}

// //--------------------------------------------------------------
// Zaki::Vector::DataSet* TaskManager::GetSequenceGrid() 
// {
// return &sequence_grid ;
// }

//--------------------------------------------------------------
Zaki::Math::Axis TaskManager::GetVisibleAxis() const
{
  return v_ax ;
}

//--------------------------------------------------------------
Zaki::Math::Axis TaskManager::GetDarkAxis() const
{
  return d_ax ;
}

//--------------------------------------------------------------
void TaskManager::FindCriticalCurve()
{

  // Zaki::Vector::DataSet ds(dir.ParentDir() + "/5_Grid2D", "Full_Mixed_Sequence_50x50.tsv") ;
  // Zaki::Vector::DataSet ds(dir.ParentDir() + "/Round 1", "Mixed_Sequence_200x200.tsv") ;
  // Zaki::Vector::DataSet ds(wrk_dir +"/NStar/Dark_Core/Linear_100x100", "Full_Mixed_Seq_Linear_100x100.tsv") ;
  // Zaki::Vector::DataSet ds(wrk_dir +"/NStar/Dark_Core/Linear_50x50", "Full_Mixed_Seq_Linear_50x50.tsv") ;
  // Zaki::Vector::DataSet ds(wrk_dir +"/NStar/Dark_Core/Linear_256x256", "Full_Mixed_Seq_Linear_256x256.tsv") ;
  // Zaki::Vector::DataSet ds(wrk_dir +"/NStar/Dark_Core/testing/No crust/10000_linear_20km", "Seq_10K_L_20km_nocrust_96x96.tsv") ;
  // Zaki::Vector::DataSet ds(wrk_dir + "/NStar/Dark_Core/testing/No crust/" + rad_res 
  //                                  + "K_Adaptive", rad_res + "K_Adaptive.tsv") ;


  double B_vis_max = sequence_grid[6].Max() ;
  double B_vis_min = sequence_grid[6].Min() ;
  
  double B_dar_max = sequence_grid[12].Max() ;
  double B_tot_max = B_vis_max + B_dar_max ;

  // Zaki::Vector::DataSet B_tot_ds({sequence_grid[0], sequence_grid[1], sequence_grid[6] + sequence_grid[12]}) ;
  // Zaki::Math::GridVals_2D B_tot_grid(B_tot_ds, 0, 1, 2) ;
  Zaki::Math::GridVals_2D B_vis_grid(sequence_grid, 0, v_ax.res, 1, d_ax.res, 6) ;
  Zaki::Math::GridVals_2D B_dar_grid(sequence_grid, 0, v_ax.res, 1, d_ax.res, 12) ;
  Zaki::Vector::DataSet M_tot_ds({sequence_grid[0], sequence_grid[1], sequence_grid[3] + sequence_grid[9]}) ;
  Zaki::Math::GridVals_2D M_tot_grid(M_tot_ds, 0, v_ax.res, 1, d_ax.res, 2) ;
  
  // ...........................
  // Setting up the grid
  // ...........................
  // { {X_min, X_max}, X_Res, X_scale}, {Y_min, Y_max}, Y_Res, Y_scale} }
  // Zaki::Math::Grid2D grid = {{{1e+15, 4e+15}, 99, "Linear"},
                              // {{3e+13, 5e+15}, 99, "Linear"}};
  // Zaki::Math::Grid2D grid = {{{1e+15, 4e+15}, 49, "Linear"},
  //                             {{3e+13, 5e+15}, 49, "Linear"}};
  // Zaki::Math::Grid2D grid = {{{8e+14, 5e+15}, 99, "Log"},
  //                             {{5e+13, 5e+16}, 150, "Log"}};
  // Zaki::Math::Grid2D grid =  {{{1e+15, 4e+15}, 255, "Linear"},
  //                                 {{3e+13, 5e+15}, 255, "Linear"}};
  // Zaki::Math::Grid2D grid = {{{1e+15, 4e+15}, 95, "Linear"},
  //                             {{3e+13, 1.5e+15}, 95, "Linear"}} ;

  //  Zaki::Math::Grid2D grid = {{{1e+15, 2.3e+15}, 159, "Linear"}, 
  //                             {{3e+13, 1.4e+15}, 159, "Linear"}} ;

  CONFIND::ContourFinder con ;
  con.SetGrid({v_ax, d_ax}) ;

  // Setting the work directory
  con.SetWrkDir(wrk_dir) ;
  
  // Setting the contour level values

  std::vector<double> cont_lvl ;
  cont_lvl.reserve(200) ;
  // for (double f = 0.001; f < 0.1; f += 0.01) // B_dar
  // for (double f = 0.67; f < 0.725; f += 0.003) // B_Tot
  // for (double f = 0.89; f < 0.9995; f += 0.01) // B_vis
  for (double f = 0.7; f < 0.96; f += 0.005) // B_vis
  // double step = (0.9995*B_vis_max - B_vis_min*1.8)/20 ;
  // std::cout << "\n\t B_vis_max = "<< B_vis_max << ", B_vis_min = " << B_vis_min << ", step = " << step << "\n" ;
  // for (double f = B_vis_min*1.8 ; f < 0.9995*B_vis_max; f += step) // B_vis
  // for (double m = 2.01-0.04; m < 2.01+0.04; m += 0.005) // M_Tot
  {
    // cont_lvl.emplace_back(m) ;
    // cont_lvl.emplace_back(B_tot_max*f) ;
    cont_lvl.emplace_back(f*B_vis_max) ;
    // cont_lvl.emplace_back(B_dar_max*f) ;
  }

  for (double f = 0.96; f < 0.9997; f += 0.001) // B_vis
  {
    cont_lvl.emplace_back(f*B_vis_max) ;
  }
  // con.SetContVal({B_tot_max*0.7}) ;
  // con.SetContVal({2.01-0.04, 2.01, 2.01+0.04}) ;
  // std::vector<std::string> cont_labs ;
  // for (double f = 2.78; f < 2.79; f += 0.002) // B_Tot
  // {
  //     cont_lvl.emplace_back(1e57*f) ;
  //     char tmp[50] ;
  //     snprintf(tmp, sizeof(tmp), "%.3e", f*1e57) ;
  //     cont_labs.emplace_back(tmp) ;
  // }
  // con.SetContVal({2.697e+57 }) ;
  // con.SetContVal(cont_lvl, cont_labs) ;
  con.SetContVal(cont_lvl) ;
  // ...........................

  // ...........................
  // Generating results
  // ...........................
  con.SetGridVals(&B_vis_grid) ;
  // con.SetGridVals(&B_tot_grid) ;
  // con.SetGridVals(&M_tot_grid) ;

  // cont_lvl.clear() ;
  // for (double f = 0.026; f < 0.275; f += 0.025) // B_dar
  // for (double m = 2.01-0.04; m < 2.01+0.04; m += 0.005) // M_Tot
  // {
    // cont_lvl.emplace_back(B_dar_max*f) ;
    // cont_lvl.emplace_back(m) ;
  // }
  // con.SetContVal(cont_lvl) ;
  // con.SetContVal({2.01-0.04}) ;
  // con.SetGridVals(&B_dar_grid) ;
  // con.SetGridVals(&M_tot_grid) ;
  con.SetPlotConnected() ;

  char tmp[100] ;
  snprintf(tmp, sizeof(tmp), "%.1f", m_chi) ;

  con.Plot("NStar/Dark_Core/" + std::string(tmp) +
                          "/B_vis_conts") ;

  M_tot_grid.Interpolate({v_ax, d_ax}) ;

  std::vector<CONFIND::Cont2D> con_set = con.GetContourSet() ;

  // Zaki::Vector::DataSet critical_pt_set(3, con_set.size()) ;
  // critical_pt_set.SetWrkDir(wrk_dir + 
  //                           "/NStar/Dark_Core/testing/No crust/" +
  //                           rad_res + "K_Adaptive") ;
  // critical_pt_set[0].label = "ec(g/cm^3)" ;
  // critical_pt_set[1].label = "ec_d(g/cm^3)" ;
  // critical_pt_set[2].label = "M_tot" ;

  // std::vector<Zaki::Math::Coord2D> critical_curve ;
  critical_curve.Reserve(con_set.size()) ;

  size_t divisions = 1 ;
  for (size_t c = 0; c < con_set.size(); c++)
  {
    // std::cout << " Contour # : " << c << ", val = " << con_set[c].GetVal() << "\n" ;
    double tmp_max_mass_tot = 0 ;
    
    // critical_pt_set.AppendRow({0,0,0}) ;
    critical_curve.Append({0,0}) ;

    for (size_t p = 0; p < con_set[c].size() - 1; p++)
    {   
      // curv.emplace_back(con_set[c][p], con_set[c][p+1]) ;
      Zaki::Math::Segment seg(con_set[c][p], con_set[c][p+1]) ;

      for (int l = 0; l < divisions ; l++)
      {
        double lam = l * (1.0 / divisions) ;
        double tmp_mass_tot = M_tot_grid.Evaluate(seg.P(lam).x, seg.P(lam).y) ;
        
        if (tmp_mass_tot > tmp_max_mass_tot)
        {
          // std::cout << " seg.P(lam) = " << seg.P(lam) << "\n" ;
          // std::cout << "-> lam = " << lam << "\n" ;
          tmp_max_mass_tot = tmp_mass_tot ;
          
          // critical_pt_set[0][-1] = seg.P(lam).x ;
          // critical_pt_set[1][-1] = seg.P(lam).y ;
          // critical_pt_set[2][-1] = tmp_max_mass_tot ;

          critical_curve[-1] = seg.P(lam) ;
        }
      }
    }
  }

  // critical_pt_set.Plot(0,1, "B_vis_critical_" + rad_res 
  //                           + "_" + std::to_string(divisions) 
  //                           + ".pdf", 
  //                           "Critical values on B_vis contours") ;

  // std::cout << " critical_curve.size() = " << critical_curve.Size() << "\n" ;
  // critical_curve.Plot(wrk_dir + "/NStar/Dark_Core/testing/No crust/" +
  //                           rad_res + "K_Adaptive/Crit_curve.pdf") ;
  // critical_curve = Moving_Ave(critical_curve, 10) ;
  critical_curve.MakeSmooth(10) ;
  
  // Adding the purely visible sector maximum mass point to critical_curve
  size_t m_max_v_idx = sequence_grid[3].MaxIdx() ;
  critical_curve.Append({ sequence_grid[2][m_max_v_idx], sequence_grid[8][m_max_v_idx]}) ;

  critical_curve.Plot(wrk_dir + "/NStar/Dark_Core/" + std::string(tmp) +
                          "/Crit_curve_smooth.pdf") ;
  critical_curve.Export(wrk_dir + "/NStar/Dark_Core/" + std::string(tmp) +
                          "/Crit_curve_smooth.tsv") ;

  // std::cout << " critical_curve.size() = " << critical_curve.Size() << "\n" ;

  
  // for (size_t r = 0; r < critical_pt_set[0].Size()-1 ; r++)
  // {
  //   critical_curve.emplace_back((Zaki::Math::Coord2D){critical_pt_set[0][r], critical_pt_set[1][r]}, 
  //                     (Zaki::Math::Coord2D){critical_pt_set[0][r+1], critical_pt_set[1][r+1]}) ;
  // }

  // FindMtotContour(2.01) ;
  // FindBtotContour() ;

#if 0
  for (size_t c = 0; c < con_set.size(); c++)
  {
    double tmp_max_mass_tot = 0 ;
    critical_pt_set.AppendRow({0,0,0}) ;
    for (size_t p = 0; p < con_set[c].size(); p++)
    {   
      // if (con_set[c][p].x > 3.9e15)
      // {
      //   continue ; 
      // }
      
      double tmp_mass_tot = M_tot_grid.Evaluate(con_set[c][p].x, con_set[c][p].y) ;
      if (tmp_mass_tot > tmp_max_mass_tot)
      {
        tmp_max_mass_tot = tmp_mass_tot ;

        critical_pt_set[0][-1] = con_set[c][p].x ;
        critical_pt_set[1][-1] = con_set[c][p].y ;
        critical_pt_set[2][-1] = tmp_max_mass_tot ;
      }            
    }
  }

  critical_pt_set.Plot(0,1, "B_vis_critical_" + rad_res + ".pdf", 
                            "Critical values on B_vis contours") ;
  // critical_pt_set.Export("B_vis_critical.tsv") ;
#endif

  // con.ExportContour("Dark Core/B_vis_conts/B_Vis_Pts.tsv",  Zaki::File::FileMode::Write) ;
  // con.ExportContour("Dark Core/B_tot_conts/B_Tot_Pts.tsv",  Zaki::File::FileMode::Write) ;

  // con.SetPlotXRange(3e13, 1e15) ;
  // con.SetPlotYRange({3e13, 1.45e15}) ;

  //..............................................................
  // Generating plot using Root
  // File name, plot name, x-axis label, y-axis label
  // con.SetPlotConnected() ;
  
  // If choosing "user" option make sure to set the coordinates!
  // con.MakeLegend(true, "Contours", "user") ;
  // con.GetLegend()->SetX1(0.75) ; con.GetLegend()->SetY1(0.75) ;
  // con.GetLegend()->SetX2(0.90) ; con.GetLegend()->SetY2(0.90) ;
  // con.GetLegend()->SetTextSize(0.025) ;
  // con.Plot("Dark Core/B_tot_conts/B_tot_100", "B_tot (100x100)", "eps_v", "eps_d") ;
  // con.Plot("NStar/Dark_Core/Linear_50x50/B_vis_50", "B_vis (50x50)", "eps_v", "eps_d") ;
  // con.Plot("NStar/Dark_Core/Linear_256x256/B_vis_256", "B_vis (256x256)", "eps_v", "eps_d") ;
  // con.Plot("NStar/Dark_Core/testing/No crust/10000_linear_20km/B_vis_96", "B_vis (96x96)", "eps_v", "eps_d") ;
  // con.Plot("NStar/Dark_Core/testing/No crust/10000_linear_20km/B_tot_50", "B_tot (50x50)", "eps_v", "eps_d") ;
  // con.Plot("NStar/Dark_Core/testing/No crust/" + rad_res 
  //           + "K_Adaptive/B_vis_" + rad_res, ("B_vis (96x96, "+rad_res+")").c_str(), "eps_v", "eps_d") ;
  //..............................................................
}

//--------------------------------------------------------------
void TaskManager::FindMtotContour(const double& in_mass)
{
  Zaki::Vector::DataSet M_tot_ds({sequence_grid[0], sequence_grid[1], sequence_grid[3] + sequence_grid[9]}) ;
  Zaki::Math::GridVals_2D M_tot_grid(M_tot_ds, 0, v_ax.res, 1, d_ax.res, 2) ;   

  CONFIND::ContourFinder con ;
  con.SetGrid({v_ax, d_ax}) ;

  // Setting the work directory
  con.SetWrkDir(wrk_dir) ;
  
  // Setting the contour level values
  con.SetContVal({in_mass}) ;

  // Evaluating the contours
  con.SetGridVals(&M_tot_grid) ;

  std::vector<CONFIND::Cont2D> mass_con_set = con.GetContourSet() ;

  Zaki::Math::Curve2D mass_curve = mass_con_set[0].ConvertToCurve2D() ;

  // char tmp_chi[50] ;
  // snprintf(tmp_chi, sizeof(tmp_chi), "%.1f", m_chi) ;
  // std::string chi_str(tmp_chi) ;

  critical_curve.Import(wrk_dir + "/NStar/Dark_Core/" + chi_str +
                            "/Crit_curve_smooth.tsv") ;

  char tmp[50];
  snprintf(tmp, sizeof(tmp), "%.2f", in_mass) ;

  con.SetPlotConnected() ;
  con.Plot("NStar/Dark_Core/" + chi_str +
            "/M_tot_" + std::string(tmp),
            ("M_tot = " + std::string(tmp)).c_str(),
             "eps_v", "eps_d") ;
  con.ExportContour("NStar/Dark_Core/" + chi_str +
                    "/M_tot_" + std::string(tmp), Zaki::File::FileMode::Write) ;

  // .................................................................
  // Finding the intersection of critical curve and const M contour
  std::vector<Zaki::Math::Coord2D> intersection = critical_curve.Intersection(mass_curve) ;

  if(intersection.size() != 1)
  {
    Z_LOG_ERROR("There should be exactly one intersection!") ;
  } 
  else
  {
    std::vector<Zaki::Math::Coord2D> m_tot_range = { mass_curve[0], intersection[0] } ;

    Zaki::File::VecSaver vec_saver(wrk_dir + "/NStar/Dark_Core/" 
                                  + chi_str + "/Intersection_" +
                                  std::string(tmp) + ".tsv") ;
    vec_saver.Export1D(intersection) ;

    FindBtotContour(in_mass, m_tot_range) ;
  }
  // .................................................................

}

//--------------------------------------------------------------
void TaskManager::FindBtotContour(const double& in_mass, 
              const std::vector<Zaki::Math::Coord2D>& in_m_range)
{
  
  char tmp_m[50] ;
  snprintf(tmp_m, sizeof(tmp_m), "%.2f", in_mass) ;
  std::string m_str(tmp_m) ;

  double B_vis_max = sequence_grid[6].Max() ;
  double B_dar_max = sequence_grid[12].Max() ;
  double B_tot_max = B_vis_max + B_dar_max ;

  Zaki::Vector::DataSet B_tot_ds({sequence_grid[0], sequence_grid[1], sequence_grid[6] + sequence_grid[12]}) ;
  Zaki::Math::GridVals_2D B_tot_grid(B_tot_ds, 0, v_ax.res, 1, d_ax.res, 2) ;
  B_tot_grid.Interpolate({v_ax, d_ax}) ;


  CONFIND::ContourFinder con ;
  con.SetGrid({v_ax, d_ax}) ;

  // Setting the work directory
  con.SetWrkDir(wrk_dir) ;
  
  // Setting the contour level values
  std::vector<double> cont_lvl ;
  
  std::cout << " p_m(0) = " << in_m_range[0] << "\n" ;
  std::cout << " p_m(1) = " << in_m_range[1] << "\n" << std::flush ;

  double cont_lvl_min = B_tot_grid.Evaluate(in_m_range[0].x, in_m_range[0].y) ;
  double cont_lvl_max = B_tot_grid.Evaluate(in_m_range[1].x, in_m_range[1].y) ;

  cont_lvl.reserve(10) ;
  // int division = 5 ;
  // for (double f = 0.67; f < 0.725; f += 0.003) // B_Tot
  for (int i = 0; i < cont_divisions; i++) // B_Tot
  {
    cont_lvl.emplace_back( cont_lvl_min + (cont_lvl_max - cont_lvl_min)*i/(cont_divisions-1) ) ;
  }

  con.SetContVal(cont_lvl) ;
  con.SetGridVals(&B_tot_grid) ;


  std::vector<CONFIND::Cont2D> B_con_set = con.GetContourSet() ;

  std::vector<Zaki::Math::Curve2D> B_curv_stable ;

  // char tmp_chi[50] ;
  // snprintf(tmp_chi, tmp_chi, "%.1f", m_chi) ;
  // std::string chi_str(tmp_chi) ;

  for (size_t i = 0; i < B_con_set.size(); i++)
  {
    Zaki::Math::Curve2D B_curve = B_con_set[i].ConvertToCurve2D() ;
    std::vector<Zaki::Math::Coord2D> intersection = critical_curve.Intersection(B_curve) ;

    if(intersection.size() != 1)
    {
      Z_LOG_ERROR("There should be exactly one intersection!") ;
      continue ;
    }

    std::cout<< " -> " << intersection[0] << "\n" ;
    B_curv_stable.emplace_back(B_curve.Bisect(intersection[0]).first) ;

    char tmp_val[50] ;
    snprintf(tmp_val, sizeof(tmp_val), "%.9e", B_con_set[i].GetVal()) ;
    B_curv_stable[B_curv_stable.size()-1].SetLabel(tmp_val) ;

    B_curv_stable[B_curv_stable.size()-1].Export(wrk_dir + "/NStar/Dark_Core/" +
                            chi_str + "/B_conts/" +
                            m_str + "/B_tot_" + m_str + "_"+std::to_string(i)+".tsv") ;
  }

  Zaki::Math::Curve2D::Plot(B_curv_stable, wrk_dir + "/NStar/Dark_Core/" +
                            chi_str + "/B_conts/" + 
                            m_str + "/B_tot_bisected.pdf") ;


  con.SetPlotConnected() ;
  con.Plot("NStar/Dark_Core/" + chi_str 
            + "/B_conts/" + m_str + "/B_tot_" + m_str,
            "B_tot", "eps_v", "eps_d") ;
}

//--------------------------------------------------------------
void TaskManager::Precision_Task(const double& in_mass) const 
{
  Zaki::Util::LogManager::SetLogLevels(Zaki::Util::LogLevel::Warning) ;

  char tmp_m[50] ;
  snprintf(tmp_m, sizeof(tmp_m), "%.2f", in_mass) ;
  std::string m_str(tmp_m) ;

  for (size_t i = 0; i < cont_divisions; i++)
  {
    CompactStar::Contour B_cont ;
    B_cont.Import(wrk_dir + "/NStar/Dark_Core/" 
                          + chi_str + "/B_conts/" + m_str  
                         + "/B_tot_"+m_str+"_"+std::to_string(i)+".tsv") ;
  
  // std::cout << "\t Val = " << B_cont.val  
  //           << ", Size = " << B_cont.curve.Size() 
  //           << ", label = " << B_cont.curve.GetLabel() 
  //           << "\n" ;

  // return ;

    CompactStar::DarkCore_Analysis darkcore_analysis ;

    CompactStar::TOVSolver solver ;
  // solver.SetWrkDir(dir.ParentDir() +"/results") ;
    solver.SetWrkDir(wrk_dir) ;

    // solver.ImportEOS_Dark("EOS/Fermi_Gas_0.8mn.eos") ;

    solver.ImportEOS(vis_eos_dir, dar_eos_dir) ;
  // solver.AddMixedCondition(TrueCondition) ;
    darkcore_analysis.SetLabel(m_str + "_" + std::to_string(i)) ;
    solver.AddAnalysis(&darkcore_analysis) ;
    solver.Solve_Mixed(B_cont, 
      "NStar/Dark_Core/" + chi_str +
      "/B_conts/" + m_str, 
      m_str + "_" + std::to_string(i)) ;
  }

  Zaki::Util::LogManager::SetLogLevels(Zaki::Util::LogLevel::Info) ;
}
//--------------------------------------------------------------
void TaskManager::FindLimits(const double& in_mass) const 
{
  char tmp_m[50] ;
  snprintf(tmp_m, sizeof(tmp_m), "%.2f", in_mass) ;
  std::string m_str(tmp_m) ;

  CompactStar::DarkCore_Analysis darkcore_analysis ;

  for (size_t i = 0; i < cont_divisions; i++)
  {

    darkcore_analysis.FindLimits(in_mass, 
                wrk_dir + "/NStar/Dark_Core/" 
                + chi_str + "/B_conts/" + 
                m_str + "/BNV_rates_"+m_str+"_"+std::to_string(i)+".tsv") ;
  }

  darkcore_analysis.ExportBNV( wrk_dir + "/NStar/Dark_Core/" 
                + chi_str + "/B_conts/" + m_str ) ;

}

//--------------------------------------------------------------

//==============================================================
