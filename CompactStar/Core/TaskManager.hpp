// -*- lsst-c++ -*-
/*
* CompactStar
* See License file at the top of the source tree.
*
* Copyright (c) 2023 Mohammadreza Zakeri
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

/**
 * @file TaskManager.hpp
 *
 * @brief Task manager for multithread parameter scanning.
 *
 * @ingroup Core
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
#ifndef CompactStar_TaskManager_H
#define CompactStar_TaskManager_H

// #include <vector>
#include <thread>
#include <Zaki/Math/Math_Core.hpp>

#include "CompactStar/Core/Prog.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
//==============================================================
namespace CompactStar
{

class MixedStar ;
class MixedSequence ;
// ========================================================
class TaskManager : public Prog
{
  //--------------------------------------------------------------
private:
  std::__1::chrono::steady_clock::time_point t_start ;
  
  Zaki::String::Directory dar_eos_dir = "" ;
  Zaki::String::Directory vis_eos_dir = "" ;
  double m_chi ;
  std::string chi_str = "m_chi" ;
  
  size_t num_of_thrds = 2 ;

  std::vector<std::thread> threads ;
  Zaki::Math::Cond_Polygon c_poly ;
  Zaki::Math::Axis v_ax ;
  Zaki::Math::Axis d_ax ;

  // std::vector<Zaki::Math::Segment> critical_curve ;
  Zaki::Math::Curve2D critical_curve ;

  // MixedSequence sequence_grid ;
  Zaki::Vector::DataSet sequence_grid ;

  void Task(const int tsk_id) const ;
  unsigned short int cont_divisions = 10 ;
  //--------------------------------------------------------------
public:
  TaskManager(const size_t& in_num_thrds);
  ~TaskManager();

  void SetDarEOSDir(const Zaki::String::Directory&) ;
  void SetVisEOSDir(const Zaki::String::Directory&) ;
  void FindDarkEOS(const double& in_m) ;
  void SetChiMass(const double& in_m) ;

  TaskManager* SetExclusionRegion(const Zaki::Math::Cond_Polygon&) ;
  TaskManager* SetGrid( const Zaki::Math::Axis& v_ax, 
                const Zaki::Math::Axis& d_ax) ;
  // TaskManager* SetWrkDir(const Zaki::String::Directory& in_dir) override ;
  Zaki::Math::Axis GetVisibleAxis() const ;
  Zaki::Math::Axis GetDarkAxis() const ;

  void Work() ;
  void ImportSequence(const Zaki::String::Directory& in_dir) ;
  void FindCriticalCurve() ;
  void FindBtotContour(const double& in_mass, const std::vector<Zaki::Math::Coord2D>&) ;
  void FindMtotContour(const double& in_mass) ;
  void Precision_Task(const double& in_mass) const ;
  void FindLimits(const double& in_mass) const ;
};
// ========================================================

//==============================================================

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_TaskManager_H*/
