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
 * @file Analysis.hpp
 *
 * @brief Virtual class for analyzing the NS in TOVSolver.
 *
 * @ingroup Core
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@eku.edu
 *
*/

#ifndef CompactStar_Analysis_H
#define CompactStar_Analysis_H

#include "CompactStar/Core/Prog.hpp"

//==============================================================
namespace CompactStar
{

class MixedStar ;
class NStar ;

//==============================================================
/**
 * @class Analysis
 *
 * @brief Abstract base class for implementing neutron star (NS) analysis routines.
 *
 * The `Analysis` class is a virtual interface meant to be extended by specific
 * analysis strategies. It provides hooks for analyzing both pure neutron stars
 * (`NStar`) and dark-visible mixed stars (`MixedStar`) as part of TOV sequence scans.
 * 
 * Results can be exported using the `Export()` method. A label is also associated
 * with the analysis for easy identification.
 */
class Analysis : public Prog
{
  protected:
    /**
     * @brief Label associated with the analysis.
     *
     * This string can be used to tag output files or to identify which type
     * of analysis was performed.
     */
    std::string label ;

  //--------------------------------------------------------------
  public:

    /**
     * @brief Default constructor.
     */
    Analysis();

    /**
     * @brief Virtual destructor.
     */
    ~Analysis() ;

    /**
     * @brief Pure virtual function to analyze a mixed star.
     *
     * This method must be implemented by any derived class.
     *
     * @param in_star Pointer to the MixedStar object to analyze.
     */
    virtual void Analyze(MixedStar* in_star) = 0 ;

    /**
     * @brief Pure virtual function to analyze a neutron star.
     *
     * This method must be implemented by any derived class.
     *
     * @param in_star Pointer to the NStar object to analyze.
     */
    virtual void Analyze(NStar* in_star) = 0 ;

    /**
     * @brief Pure virtual function to export analysis results.
     *
     * @param dir Directory where output should be saved.
     */
    virtual void Export(const Zaki::String::Directory&) = 0 ;

    /**
     * @brief Sets a label for the analysis.
     *
     * @param in_label The label to assign.
     */
    virtual void SetLabel(const std::string& in_label) ;
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_Analysis_H*/