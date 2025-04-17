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
 * @file Prog.hpp
 *
 * @brief The program's base class.
 *
 * @ingroup Core
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Last edit Aug 31, 2021
#ifndef CompactStar_Prog_H
#define CompactStar_Prog_H

#include "CompactStar/Core/CompactStarConfig.h"

#include <Zaki/Util/Logger.hpp>
// #include <Zaki/Util/Instrumentor.hpp>
// #include <Zaki/Util/Simple_Timer.hpp>
// #include <Zaki/String/Directory.hpp>

#include "stdio.h"
#include <atomic>

#define CompactStar_PROG_DEBUG_MODE 0
//==============================================================
namespace CompactStar
{
//==============================================================
class Prog
{
  //--------------------------------------------------------------
  protected:

    /// Prints the banner of the program
    void ShowBanner();

    /// The number of derived object instances
    static inline std::atomic<size_t> counter = 0;

    /// The work directory
    Zaki::String::Directory wrk_dir = "" ;

    /// The label of the class
    std::string name        = "" ;

    // Flags

    /// If the name is set
    bool set_name_flag = false    ;

    ///  If the work directory is set
    bool set_wrk_dir_flag = false ;
    
    /// Sets the work directory for the member objects
    virtual Prog* SetMemWrkDir(const Zaki::String::Directory&) ;
  //--------------------------------------------------------------
  private:

#if CompactStar_PROG_DEBUG_MODE
    /// Track object creation and destruction
    bool track_objs_flag = true ;
#endif
  //--------------------------------------------------------------
  public:

    /// Default Constructor for faster object creation
    ///  -> No name assignment or object tracking
    Prog() ;
    
    /// Constructor via a name and no tracking options
    /// -> Name is set but no object tracking
    Prog(const std::string&) ;

    /// Constructor via a name and tracking options
    /// This is slower and shouldn't be used for objects 
    /// that are created thousands of time, i.e. particles, etc.
    /// -> Name is set and object tracking is on
    Prog(const std::string&, const bool&) ;

    /// Destructor
    virtual ~Prog() ;

    /// Overloading CLass specific new operator
    static void* operator new(size_t sz) ;

    /// Overloading CLass specific delete operator
    static void operator delete(void* m);

    // Setters

    /// Sets the work directory
    virtual Prog* SetWrkDir(const Zaki::String::Directory&) ;

    /// Sets the name of the class
    virtual void SetName(const std::string&) ;

    // void SetObjTracking(const bool) ;

    // Getters

    /// Returns the work directory
    virtual Zaki::String::Directory GetWrkDir() const ;

    /// Returns set_wrk_dir_flag
    bool IsWrkDirSet() const ;

    /// Returns the name of the class
    virtual std::string GetName() const ;

    /// Print method for the derived classes
    virtual void Print() const ;

    /// String form of the pointer to the object
    std::string PtrStr() const ;

};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_Prog_H*/
