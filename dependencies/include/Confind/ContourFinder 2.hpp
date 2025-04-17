#ifndef ContourFinder_H
#define ContourFinder_H

#if __has_include(<omp.h>)
  #include <omp.h>
  #define Z_OMP 1
#else
  #define Z_OMP 0
#endif

#include <random>

// Root
// #include <TMultiGraph.h>
// #include <TLegend.h>

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <Zaki/Math/Func2D.hpp>
#include <Zaki/Math/MemFuncWrapper.hpp>
#include <Zaki/Math/Math_Core.hpp>
#include <Zaki/Vector/DataSet.hpp>


// Local headers
#include "Confind/Cell.hpp"
#include "Confind/Cont2D.hpp"

//==============================================================

namespace CONFIND
{

//==============================================================
class ContourFinder : public Base
{

  friend class Cell ; 
 //--------------------------------------------------------------
  public:

    //............................................
    /// Constructor
    ContourFinder() ;
    
    /// Destructor
    ~ContourFinder() ;

    /// Copy constructor
    ContourFinder(const ContourFinder &zc2) ;

    /// Assignment operator
    ContourFinder& operator=(const ContourFinder &zc2) ;
    //............................................


    enum Mode
    {
      Normal = 0, Fast, Parallel, Ludicrous, Optimal    
    } ;

    //............................................
    // Setters
    //............................................
    void SetGrid(const Zaki::Math::Grid2D&)  ;
    void SetWidth(const size_t&) ;
    void SetHeight(const size_t&) ;
    void SetDeltas()      ;
    void SetGridVals(const Mode& = Fast)    ;
    void SetGridVals(Zaki::Math::GridVals_2D*)    ;

    void SetFunc(double (*f) (const double, const double) ) ; // Normal funcs 
    void SetMemFunc(Zaki::Math::Func2D*) ;  // Non-static mem-funcs
    void SetContVal(const std::vector<double>&) ;
    void SetContVal(const std::vector<double>&, const std::vector<std::string>& label) ;
    void SetScanMode(const char) ;
    void SetOptimizationTrials(const int&) ;
    void SetThreads(const int&) ;    

    // Plot options
    void SetPlotXLabel(const std::string&) ;
    void SetPlotYLabel(const std::string&) ;
    void SetPlotXRange(const Zaki::Math::Range<double>&) ;
    void SetPlotYRange(const Zaki::Math::Range<double>&) ;
    void SetPlotLabel(const std::string&)  ;
    void SetPlotConnected(const bool=true) ;
    void SetLegendLabels(const std::vector<std::string>& ); 
    void MakeLegend(const bool=true, const char* const=nullptr, const char* const=nullptr) ;
    //............................................

    void Clear() ;

    //............................................
    // Getters
    //............................................
    size_t GetN_X()   const ;
    size_t GetN_Y()   const ;
    std::string GetXScale() const ;
    std::string GetYScale() const ;
    std::pair<double, double> GetDeltas()  const ;
    double GetX_Min() const ;
    double GetX_Max() const ;
    double GetY_Min() const ;
    double GetY_Max() const ;
    char GetScanMode() const;
    
    /// Returns 'cont_set'
    std::vector<Cont2D> GetContourSet() const ;

    std::pair<double, double> ij_2_xy(const size_t i, const size_t j) const ;
    
    void Print() const override;
    void ExportContour(const Zaki::String::Directory& f_name, const Zaki::File::FileMode& mode) ;

    // TMultiGraph* GetGraph() ;
    // TLegend* GetLegend() ;

    void Plot(const Zaki::String::Directory& f_name, const char* const=nullptr, 
              const char* const=nullptr, const char* const=nullptr) ;

 //--------------------------------------------------------------
  private:

    void FindContour(Cont2D& cont) ;
    void FindContourFast(Cont2D& cont, double*) ;
    void FindNextContours(double*) ;
    void FindContourLudicrous() ;
    void ThreadNextTaskLudicrous(Bundle& in_b, const std::vector<double>&) ;
    void ThreadTaskLudicrous(Bundle& in_b, std::vector<double>&) ;
    void FindContourParallel(Cont2D& cont) ;
    void FindOptimalMode() ;
    void ThreadTask(Bundle& in_b) ;

    // flags
    bool set_height_flag      = false ;
    bool set_width_flag       = false ;
    bool set_grid_flag        = false ;
    bool set_grid_vals_flag   = false ;
    bool set_func_flag        = false ;
    bool set_cont_val_flag    = false ;
    bool set_mem_func_flag    = false ;
    bool set_plotX_label_flag = false ;
    bool set_plotY_label_flag = false ;
    bool set_plotX_range_flag = false ;
    bool set_plotY_range_flag = false ;
    bool set_plot_label_flag  = false ;
    bool set_plot_connected_flag = false ;
    bool cpy_cons_called      = false ;
    bool set_scan_mode_flage  = false ;
    bool set_leg_lab_flag     = false ;
    bool make_legend_flag     = false  ;

    char scan_mode = 'X' ; 
    Mode algorithm ;

    //............................................
    /// Pointer to the function
    double (*func) (double, double)= nullptr ;
    /// Unique pointer to the member-function
    std::unique_ptr<Zaki::Math::Func2D> genFuncPtr = nullptr ;
    /// TMultiGraph* is Automatically managed by root
    // TMultiGraph* graph = nullptr ;
    /// TLegend* is Automatically managed by root
    // TLegend* legend = nullptr ;
    //............................................

    std::vector<std::string> legend_label_set ;
    std::string legend_header = "Contours" ;
    bool default_legend_opt = true ; 

    unsigned int optimization_trials = 50 ;
    Mode TimeFunc(std::uniform_real_distribution<double>&, 
                  std::uniform_real_distribution<double>&) ;

    std::default_random_engine random_engine ;

    Zaki::Math::Grid2D grid ;
    unsigned int width = 1000, height = 1000 ;
    double delta_x, delta_y ;
    std::string x_label="X", y_label ="Y", plot_label="" ;
    Zaki::Math::Range<double> plot_x_range, plot_y_range ;
    bool connected_plot = false  ;
    std::vector<Cont2D> cont_set ;

    int req_threads = 1 ;
    unsigned int unfound_contours = 0 ;
};

//==============================================================
/// Interface for the user to set member function pointers
template<typename FuncObj, typename MemFuncPtr >
class MemFuncContWrapper : public Base
{

  private:
    ContourFinder cont_finder;

  public:
    MemFuncContWrapper(const FuncObj& obj, const MemFuncPtr& memFn)
      : Base("MemFuncContWrapper", true)
    {
      Z_LOG_NOTE("Member function wrapper initializing...") ;
      cont_finder.SetMemFunc(new Zaki::Math::MemFuncWrapper<FuncObj, double (FuncObj::*)(double, double)>(obj, memFn)) ;
    }
    ~MemFuncContWrapper(){}

    MemFuncContWrapper(const MemFuncContWrapper& other)
    : Base("MemFuncContWrapper", true), cont_finder(other.cont_finder)
    {
#if CONFIND_BASE_DEBUG_MODE
        Z_OBJ_CCTR(this, (void*)(&other), "MemFuncContWrapper", "MemFuncContWrapper") ;
#endif
    }

    void UpdateMemFunc(const FuncObj& obj, const MemFuncPtr& memFn)
    {
      Z_LOG_NOTE("Updating the member function wrapper...") ;
      cont_finder.SetMemFunc(new Zaki::Math::MemFuncWrapper<FuncObj, double (FuncObj::*)(double, double)>(obj, memFn)) ;
    }

    ContourFinder* operator->() {
      return &cont_finder;
    }
};
//==============================================================

//--------------------------------------------------------------
} // CONFIND namespace
//==============================================================
#endif /*ContourFinder_H*/
