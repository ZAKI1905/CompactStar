#ifndef CONFIND_Cont2D_H
#define CONFIND_Cont2D_H

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


#include <Zaki/Physics/Coordinate.hpp>
// #include <Zaki/Math/MemFuncWrapper.hpp>
// #include <Zaki/Math/Math_Core.hpp>
#include <Zaki/File/VecSaver.hpp>

#include "Confind/Base.hpp"

namespace CONFIND
{

//==============================================================
class Cont2D : public Base
{

  friend class ContourFinder  ;
  friend std::ostream  & operator<<(std::ostream &os, const Cont2D& p);
  //--------------------------------------------------------------
  public:
    
    // Constructor
    Cont2D(double) ;

    // Copy Constructor
    Cont2D(const Cont2D&) ;

    void SetColor(const CONFIND::Color&) ;
    void SetColor(const unsigned int&) ;
    void SetFound(const bool=true) ;
    void SetLabel(const std::string&) ;

    CONFIND::Color GetColor() const ;
    double GetVal() const ;
    bool GetFound() const ;

    // temproray public
    void AddPts(const std::vector<Zaki::Physics::Coord3D>&) ;
    void RMDuplicates() ;
    Zaki::Physics::Coord3D operator[](const size_t&) const ;
    Cont2D operator+(const Cont2D& c) const ;
    void operator+=(const Cont2D& c) ;
    size_t size() const ;
    void SortNew() ;

    Zaki::Math::Curve2D ConvertToCurve2D() ;

    /// Converts the point to a DataSet
    Zaki::Vector::DataSet ConvertToDataSet() ;


  private:
    double val;
    std::vector<Zaki::Physics::Coord3D> pts;
    CONFIND::Color color ;
    std::string label ;

    bool set_label_flag = false ;
    bool is_found_flag = false ;


    // The most bottom_left point (used for sorting)
    Zaki::Physics::Coord3D bottom_left ;
    int Orientation(const Zaki::Physics::Coord3D&, const Zaki::Physics::Coord3D&, const Zaki::Physics::Coord3D&) const;


    void Export(const Zaki::String::Directory& f_name, const Zaki::File::FileMode& mode) ;

    void Sort() ;

    bool comp_Orient(const Zaki::Physics::Coord3D &, const Zaki::Physics::Coord3D &) const; 
    bool sort_cw = true ;
    bool already_sorted = false ;
    void Clear() ; 
};

std::ostream& operator << ( std::ostream &output, const Cont2D& c ) ;
//==============================================================

} // CONFIND namespace
//==============================================================
#endif /*CONFIND_Cont2D_H*/