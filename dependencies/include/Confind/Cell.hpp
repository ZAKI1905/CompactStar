#ifndef Cell_H
#define Cell_H

#include <array>

#include <Zaki/Physics/Coordinate.hpp>

#include "Confind/Base.hpp"

//--------------------------------------------------------------
namespace CONFIND
{

//--------------------------------------------------------------
// Forward declaration
class Bundle ;
class ContourFinder ;
//--------------------------------------------------------------

//==============================================================
struct vertex
{
  Zaki::Physics::Coord3D xyz   ;

  // status with respect to a level ( 1: above, 0: on, -1: below)
  int status(const double& lvl) const
  {
    if ( xyz.z > lvl)
    {
      return 1 ;
    }
    else if ( xyz.z < lvl)
    {
      return -1 ;
    }
    else
      return 0 ;
  }

} ;
//==============================================================
struct triangle
{
  // std::vector<vertex> v   ;
  // vertex v[3] ;
  std::array<vertex, 3> v;

  triangle(const vertex& v_1, const vertex& v_2, const vertex& v_3) 
    : v({v_1, v_2, v_3}) {} ;

  int status(const double& cont) const
  {
    std::vector<int> set = { v[0].status(cont), v[1].status(cont), v[2].status(cont)} ;
    int sum    = set[0] + set[1] + set[2] ;
    int mul    = set[0] * set[1] * set[2] ;

    if ( mul == 1 )
    {
      // Two possibilities only:  (-1, -1, +1) or (+1, +1, +1)
      if( sum == -1) // --> (-1, -1, +1)
      // (c-3) Two vertices lie below and one above the contour level: (-1, -1, +1)
        return 3;
      else // --> (+1, +1, +1): 
      // (j-10) All the vertices lie above the contour level.
        return 10;
    }

    else if ( mul == -1 )
    {
      // Two possibilities only:  (+1, +1, -1) or (-1, -1, -1)
      if( sum == 1) // --> (+1, +1, -1)
      // (f-6)  One vertex lies below and two above the contour level: (-1, +1, +1)
        return 6;
      else // --> (-1, -1, -1): 
      // (a-1) All the vertices lie below the contour level.
        return 1;
    }
    
    else if ( sum == 2 ) // : (1, 1, 0)
    // i) One vertex lies on and two above the contour level.
      return 9;
    else if ( sum == -2 ) // : (-1, -1, 0)
    // b) One vertex lies on and two below the contour level.
      return 2;

    else if ( sum == 1 ) // : (+1, 0, 0)
    // h) Two vertices lie on and one above the contour level.
      return 8;
    else if ( sum == -1 ) // : (-1, 0, 0)
    // d) Two vertices lie on and one below the contour level.
      return 4;

    else if (set[0] != 0 || set[1] != 0 || set[2]!= 0 )
    // e) One vertex lies below, one on and one above the contour level.
      return 5;

    // g) Three vertices lie on the contour level.
    else
      return 7 ; 

  }
} ;

//==============================================================
class Cell : public Base
{
  friend class ContourFinder ;
  //------------------------------------------------

  public:
    
    // Constructor
    // Default
    Cell() ;

    // Takes: X_idx, delta_X, Y_idx, delta_Y
    Cell(const size_t& i_x, const double& lx, const size_t& i_y, const double& ly);

    // Takes: X_idx, delta_X, Y_idx, delta_Y, Pointer to Bundle, contour value
    Cell(const size_t& i_x, const double& lx, const size_t& i_y, const double& ly,
          Bundle*, const double& c) ;

    // Copy constructor 
    Cell(const Cell &c2) = delete ;
    
    // Assignment operator
    Cell& operator=(const Cell &c2) = delete ;

    // Destructor
    ~Cell() ;

    //............................................
    // Setters
    //............................................
    void SetIdx(const size_t&, const size_t&)       ;
    void SetSize(const double&, const double&)      ;
    void SetVertex(const size_t&, const double&, const double&, const double&)   ;
    void SetVertex(const size_t&, const Zaki::Physics::Coord3D&)    ;
    void SetVertexZ(const size_t&, const double&)   ;
    void SetTriangles() ;
    void SetBundlePtr(Bundle* contPtrIn) ;
    void SetContourValue(const double&) ;
    //............................................

    void FindVerts() ;
    void EvalCenter() ;

    // evaluates the function values
    double EvalFunc(const double& x, const double& y) ; 

    //............................................
    // Getters
    //............................................
    std::pair<size_t, size_t> GetIdx()   const  ;
    // std::vector<Triangle> GetTriangles()   ;
    vertex operator[](const size_t) const ; // Selecting the i-th vertex
    // Triangle GetDown()  ;   // idx = 0
    // Triangle GetRight() ;   // idx = 1
    // Triangle GetUp()   ;   // idx = 2
    // Triangle GetLeft()  ;   // idx = 3

    // size_t GetXIdx()   const ;
    // size_t GetYIdx()   const ;
    double GetLX()     const ;
    double GetLY()     const ;  
    
    double GetFuncVals(const size_t&) const ;
    double GetContourValue() const ;
    const std::vector<Zaki::Physics::Coord3D>& GetContourCoords() const ;
    int GetStatus() ;
    //............................................


  private:

    // Grid minimum x & y
    double G_x_min, G_y_min ;
    
    // Cell minimum x & y
    double x_min, y_min ;

    // size_t idx_x, idx_y;
    double l_x, l_y ;
    std::pair<size_t, size_t> idx  ;
    double contour_val             ;

    enum CellPts
    {
      Center=0, B_Left, B_Right, T_Right, T_Left
    };

    // 1 : bottom-left, 2: bottom-right, 3: top-right, 4: top-left, 0: center
    std::array<vertex, 5> verts_set = {{ {{0, 0, 0}},
                                      {{0, 0, 0}},
                                      {{0, 0, 0}},
                                      {{0, 0, 0}},
                                      {{0, 0, 0}} }};
    
    std::vector<triangle> triangle_set ;
    std::vector<Zaki::Physics::Coord3D> contour_coords;

    Bundle* BundlePtr = nullptr ;

    // flags
    bool set_idx_flag             = false ;
    std::array<bool,5> set_vertex_flag  = {{false, false, false, false, false}} ;
    std::array<bool,4> set_vertexZ_flag = {{false, false, false, false}} ;
    bool set_bundle_ptr_flag      = false ;
    bool set_size_flag            = false ;
    bool found_verts_flag         = false ;
    bool set_contour_val_flag     = false ;
    bool set_full_constructor     = false ;

    //..........................................................................

    double EvalSimpleFunc(const double& x, const double& y) ;
    double EvalMemFunc(const double& x, const double& y)    ;

    // Case 4: d) One vertex lies below and two on the contour level. (-1, 0, 0)
    // Case 8: h) Two vertices lie on and one above the contour level. (+1, 0, 0)
    void case48(const triangle&)     ; 

    // Case 3: c) Two vertices lie below and one above the contour level.
    // Case 6: f) One vertex lies below and two above the contour level.
    void case36(const triangle&, const int&)  ;  // (-1, -1, +1), (-1, +1, +1)

    // Case 5: e) One vertex lies below, one on and one above the contour level.
    void case5(const triangle&) ;
};

//==============================================================
} // CONFIND namespace
//--------------------------------------------------------------

#endif /*Contour_H*/