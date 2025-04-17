/*
  Created on May 9, 2022

  TOVSolver_Thread class
*/


#include <Zaki/Util/Instrumentor.hpp>

#include "CompactStar/Core/TOVSolver_Thread.hpp"

#define TOV_SOLVER_THREAD_VERBOSE 0

using namespace CompactStar ;

//==============================================================
//                        TOVSolver_Thread class
//==============================================================
// Constructor
TOVSolver_Thread::TOVSolver_Thread(
  const short unsigned int& in_task_id) : task_id(in_task_id)
{
  SetName("TOVSolver_Thread") ;
  tov_counter++ ;
}

//--------------------------------------------------------------
TOVSolver_Thread::~TOVSolver_Thread()
{
  PROFILE_FUNCTION() ;

  tov_counter-- ;

  {
    std::lock_guard<std::mutex> lock(m_mutex) ;

    mixed_seq_static.Combine(mixed_sequence) ;
  }

  if( tov_counter == 0 )
  {
    mixed_seq_static.Export(tov_result_dir.ThisFileDir() 
                            + "/"+seq_file_name+"_Sequence.tsv") ;
    mixed_seq_static.Clear() ;
  }
}

//--------------------------------------------------------------
void TOVSolver_Thread::SetSeqFileName(const std::string& in_file)
{
  seq_file_name = in_file ;
}

//--------------------------------------------------------------
/// Sets the work directory for the member objects
Prog* TOVSolver_Thread::SetMemWrkDir(const Zaki::String::Directory& in_dir) 
{
  mixed_star.SetWrkDir(in_dir) ;
  mixed_sequence.SetWrkDir(in_dir) ;

  {
    std::lock_guard<std::mutex> lock(m_mutex) ;

    if ( !mixed_seq_static.IsWrkDirSet() )
    {
      mixed_seq_static.SetWrkDir(in_dir) ;
    }
  } 

  return this ;
}

//--------------------------------------------------------------
void TOVSolver_Thread::ExportMixedStarProfile
  ( const size_t& v_idx, const size_t& d_idx, 
    const Zaki::String::Directory& in_dir)
{
  mixed_star.Export(in_dir + "_" + 
                    std::to_string(v_idx) + "_" +
                    std::to_string(d_idx + min_idx_offset)
                    + ".tsv") ;
}

//--------------------------------------------------------------
void TOVSolver_Thread::ExportMixedSequence
  (const Zaki::String::Directory& in_dir)
{
  tov_result_dir = in_dir ;
  // mixed_sequence.Export(in_dir.ThisFileDir() + "/" +
  //                       std::to_string(task_id) +
  //                       "_" + in_dir.ThisFile().Str() ) ;
}

//--------------------------------------------------------------
/// Sets min_idx_offset
void TOVSolver_Thread::SetMinIdxOffset(const size_t& in_idx) 
{
  min_idx_offset = in_idx ;
}

//--------------------------------------------------------------
void TOVSolver_Thread::SurfaceIsReached(const size_t& v_idx, 
                                        const size_t& d_idx)
{
  mixed_star.SurfaceIsReached(v_idx, d_idx+min_idx_offset) ;
}

//--------------------------------------------------------------
void TOVSolver_Thread::PrintStatus(const size_t& in_v_idx, 
  const size_t& in_d_idx, const size_t& in_v_res, const size_t& in_d_res) 
{
  char tmp_term[150] ;
  snprintf(tmp_term, sizeof(tmp_term), "[T-%d] Mixed sequence (%3.lu, %3.lu) out "
                    "of (%3.lu, %3.lu).\r", task_id, in_v_idx+1, 
                    in_d_idx+1, in_v_res+1, in_d_res+1) ;
  std::cout << tmp_term << std::flush ;
}

//--------------------------------------------------------------
/// Returns mixed_seq_static
MixedSequence TOVSolver_Thread::GetSequence() 
{
  return mixed_seq_static ;
}

//--------------------------------------------------------------

//==============================================================
