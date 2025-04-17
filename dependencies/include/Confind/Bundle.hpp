#ifndef CONFIND_Bundle_H
#define CONFIND_Bundle_H

// #if __has_include(<omp.h>)
//   #include <omp.h>
//   #define Z_OMP 1
// #else
//   #define Z_OMP 0
// #endif

// #include <random>

// Root
// #include <TMultiGraph.h>
// #include <TLegend.h>


#include <Zaki/Math/Func2D.hpp>
#include <Zaki/Util/ObjObserver.hpp>
#include <Zaki/Math/Math_Core.hpp>
// #include <Zaki/File/VecSaver.hpp>


// Local headers
// #include "Confind/Cell.hpp"
#include "Confind/Cont2D.hpp"
#include "Confind/Base.hpp"

namespace CONFIND
{

//==============================================================
// Bundle Class is used for sending information to each thread
//  in a multi-threaded scenario to avoid racing.
class Bundle : public Base
{
  friend class ContourFinder;
  friend class Cell;

  public:
    // Constructor 0
    Bundle() : Base("Bundle", true), Con(0)
    {
        Z_LOG_NOTE("Bundle constructor: from " + PtrStr() + ".") ;
    }

    // Constructor 1
    Bundle(const Zaki::Math::Grid2D& in_g, const Cont2D& in_c,
    const std::unique_ptr<Zaki::Math::Func2D>& in_mf) 
        : Base("Bundle", true), Grid(in_g), Con(in_c)
    {
        Z_LOG_NOTE("Bundle constructor: from " + PtrStr() + ".") ;
      if(in_mf)
        MemFunc = in_mf->Clone() ;
    }

    // Constructor 2
    Bundle(const Zaki::Math::Grid2D& in_g, const Cont2D& in_c,
    double (*in_f)(const double, const double)) 
        : Base("Bundle", true), Grid(in_g), Con(in_c), Func(in_f)
    {
        Z_LOG_NOTE("Bundle constructor: from " + PtrStr() + ".") ;
    }

    // Constructor 3
    Bundle(const Zaki::Math::Grid2D& in_g, const Cont2D& in_c,
    const std::unique_ptr<Zaki::Math::Func2D>& in_mf,
    double (*in_f)(const double, const double)) 
        : Base("Bundle", true), Grid(in_g), Con(in_c), Func(in_f)
    {
        Z_LOG_NOTE("Bundle constructor: from " + PtrStr() + ".") ;
      if(in_mf)
        MemFunc = in_mf->Clone() ;
    }
    
    // Constructor 4
    Bundle(const Zaki::Math::Grid2D& in_g,
    const std::unique_ptr<Zaki::Math::Func2D>& in_mf,
    double (*in_f)(const double, const double)) 
        : Base("Bundle", true), Grid(in_g), Con(0), Func(in_f)
    {
        Z_LOG_NOTE("Bundle constructor: from " + PtrStr() + ".") ;
      if(in_mf)
        MemFunc = in_mf->Clone() ;
    }

    // Copy constructor
    Bundle(const Bundle& other)
      : Base("Bundle", true), Grid(other.Grid), Con(other.Con), Func(other.Func)
    {
      Z_LOG_NOTE("Bundle copy constructor: from " 
              + other.PtrStr() + " --> " + PtrStr() + ".") ;

#if CONFIND_BASE_DEBUG_MODE
      Z_OBJ_CCTR(this, (void*)(&other), "Bundle", "Bundle") ;
#endif

        if(other.MemFunc)
        MemFunc = other.MemFunc->Clone() ;

    }

    ~Bundle()
    {
        Z_LOG_NOTE("Bundle destructor: from "+ PtrStr() + ".") ;
    }
    
    void AddCont(const Cont2D& in_c) { Con = in_c; }
    void AddGrid(const Zaki::Math::Grid2D& in_g) { Grid = in_g; }

    void AddMemFunc(const std::unique_ptr<Zaki::Math::Func2D>& in_mf) 
    { 
      if(in_mf)
        MemFunc  = in_mf->Clone(); 
      else
        Z_LOG_NOTE("Input member function is nullptr!") ;
    }
    
    void AddFunc(double (*in_f)(const double, const double)) 
    { 
      if(in_f)
        Func = in_f;
      else
        Z_LOG_NOTE("Input function is nullptr!") ;
    }

    Cont2D GetCont() {return Con;}
      
  private:
    Zaki::Math::Grid2D Grid;
    Cont2D Con ;
    std::unique_ptr<Zaki::Math::Func2D> MemFunc = nullptr;
    double (*Func)(const double, const double) = nullptr ;
    // std::unique_ptr<double (*)(const double, const double)> func ;

};

//==============================================================

//--------------------------------------------------------------
} // CONFIND namespace
//==============================================================
#endif /*CONFIND_Bundle_H*/
