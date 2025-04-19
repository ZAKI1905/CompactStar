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
 * @brief Base program class providing common utilities and directory handling.
 *
 * The Prog class offers a standard interface for application banner display,
 * working directory management, naming, and optional object tracking.
 *
 * @ingroup Core
 * @author Mohammadreza Zakeri
 * @contact M.Zakeri@eku.edu
 */
// Last edit Aug 31, 2021
#ifndef CompactStar_Prog_H
#define CompactStar_Prog_H

#include "CompactStar/Core/CompactStarConfig.h"
#include <Zaki/Util/Logger.hpp>

#include "stdio.h"
#include <atomic>

#define CompactStar_PROG_DEBUG_MODE 0
//==============================================================
namespace CompactStar
{
//==============================================================
//                        Prog Class
//==============================================================
/**
 * @class Prog
 * @brief Base class for CompactStar applications.
 *
 * Provides banner display, working directory and name management,
 * and optional instance tracking for derived classes.
 */
class Prog
{
  //--------------------------------------------------------------
  protected:

    /**
     * @brief Display the program banner to stdout.
     *
     * Prints name and version information if available.
     */
    void ShowBanner();

    /**
     * @brief Count of Prog-derived object instances.
     *
     * Used for optional tracking when debug mode is enabled.
     */
    static inline std::atomic<size_t> counter = 0;

    /**
     * @brief Current working directory for program outputs.
     */
    Zaki::String::Directory wrk_dir = "" ;

    /**
     * @brief Label or name identifier for the class instance.
     */
    std::string name        = "" ;

    /**
     * @brief Flag indicating if the name has been set.
     *
     */
    bool set_name_flag = false    ;

    /**
     * @brief Flag indicating if the working directory has been set.
     */
    bool set_wrk_dir_flag = false ;
    
    /**
     * @brief Internal method to set member objects' working directory.
     *
     * @param dir  Directory to assign to member Prog objects.
     * @return     Pointer to this Prog instance for chaining.
     */
    virtual Prog* SetMemWrkDir(const Zaki::String::Directory&) ;
    
  //--------------------------------------------------------------
  private:

#if CompactStar_PROG_DEBUG_MODE
    /**
     * @brief Toggle for object creation/destruction tracking.
     */
    bool track_objs_flag = true ;
#endif
  //--------------------------------------------------------------
  public:

    /**
     * @brief Default constructor without name or tracking.
     *
     * Fast instantiation: no name assignment or debug tracking.
     */
    Prog() ;
    
    /**
     * @brief Constructor with name only (no tracking).
     *
     * @param prog_name  Name label for this instance.
     */
    Prog(const std::string&) ;

    /**
     * @brief Constructor with name and optional tracking.
     *
     * @param prog_name       Name label for this instance.
     * @param enable_tracking If true, instance creation/destruction are tracked.
     * @note Slower than default constructor.
     * @note Use for debugging or when tracking is needed.
     * @note Not recommended for high-frequency object creation.
     */
    Prog(const std::string&, const bool&) ;

    /**
     * @brief Destructor for cleanup and potential tracking.
     */
    virtual ~Prog() ;

    /**
     * @brief Override global new operator for custom allocations.
     *
     * @param sz  Size in bytes to allocate.
     * @return    Pointer to allocated memory block.
     */
    static void* operator new(size_t sz) ;

    /**
     * @brief Override global delete operator for custom deallocations.
     *
     * @param ptr  Memory block to free.
     */
    static void operator delete(void* m);

    // Setters

    /**
     * @brief Set the working directory for program outputs.
     *
     * @param dir  Desired working directory path.
     * @return     Pointer to this Prog instance for chaining.
     */
    virtual Prog* SetWrkDir(const Zaki::String::Directory&) ;

    /**
     * @brief Assign a name label to this instance.
     *
     * @param prog_name  Descriptive name for object.
     */
    virtual void SetName(const std::string&) ;


    // Getters

    /**
     * @brief Retrieve the current working directory.
     *
     * @return Working directory path.
     */
    virtual Zaki::String::Directory GetWrkDir() const ;

    /**
     * @brief Check if the working directory has been set.
     *
     * @return True if SetWrkDir was called.
     */
    bool IsWrkDirSet() const ;

    /**
     * @brief Get the name label of this instance.
     *
     * @return Name string.
     */
    virtual std::string GetName() const ;

    /**
     * @brief Print object-specific information.
     *
     * Derived classes should override to display status.
     */
    virtual void Print() const ;

    /**
     * @brief Return string form of this object's pointer.
     *
     * Useful for logging and debugging.
     *
     * @return Hexadecimal pointer string.
     */
    std::string PtrStr() const ;

};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_Prog_H*/
