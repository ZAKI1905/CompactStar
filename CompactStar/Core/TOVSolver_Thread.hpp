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
 * @file TOVSolver_Thread.hpp
 *
 * @brief TOV Solver for multi-threading.
 *
 * @ingroup Core
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Created on May 9, 2022
#ifndef CompactStar_TOVSolver_Thread_H
#define CompactStar_TOVSolver_Thread_H

#include <thread>
// #include <vector>
// #include <gsl/gsl_spline.h>

// #include <Zaki/Math/Math_Core.hpp>

// #include "CompactStar/Core/Prog.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

//==============================================================
namespace CompactStar
{

//==============================================================
class TOVSolver_Thread : public TOVSolver
{
  //--------------------------------------------------------------
  private:
    
    /// The number of TOV object instances
    static inline std::atomic<size_t> tov_counter = 0;

    std::mutex m_mutex ;

    static inline MixedSequence mixed_seq_static ;

    /// Number of threads
    static const short unsigned int num_of_thrds = 2 ;
    
    /// The thread id that created this instance
    const short unsigned int task_id = 0 ;

    /// The index offset due to dividing the job between threads
    size_t min_idx_offset = 0 ;

    Zaki::String::Directory tov_result_dir = "" ;
    std::string seq_file_name ;

  //--------------------------------------------------------------
  public:
    
    TOVSolver_Thread(const short unsigned int& in_task_id) ;
    ~TOVSolver_Thread() ;

    void SetSeqFileName(const std::string& in_file) ;
    
    /// Copy Constructor 
    TOVSolver_Thread(const TOVSolver_Thread &) = delete ;

    /// Assignment operator
    TOVSolver_Thread& operator= (const TOVSolver_Thread&) = delete ;

    /// Sets the work directory for the member objects
    Prog* SetMemWrkDir(const Zaki::String::Directory&) override ;

    /// Exports the mixed sequence
    void ExportMixedSequence(const Zaki::String::Directory&) override ;

    /// Exports the mixed star profile
    void ExportMixedStarProfile(const size_t& v_idx, 
      const size_t& d_idx, const Zaki::String::Directory&) override ;


    void SurfaceIsReached(const size_t& v_idx, 
                          const size_t& d_idx) override ;


    void PrintStatus(const size_t& v_idx,
                     const size_t& d_idx,
                     const size_t& v_res, 
                     const size_t& d_res) override ;

    /// Sets min_idx_offset
    void SetMinIdxOffset(const size_t& in_idx) ;

    /// Returns mixed_seq_static
    static MixedSequence GetSequence() ;
};

//==============================================================

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_TOVSolver_Thread_H*/
