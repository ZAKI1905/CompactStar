// -*- lsst-c++ -*-
/*
* Zaki's Common Library
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
 * @file CSVIterator.hpp
 *
 * @brief Parses input files such as CVS, TSV, etc.
 *
 * @ingroup File
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_File_CSVIterator_H
#define Zaki_File_CSVIterator_H

#include "Zaki/File/CSVRow.hpp"

//--------------------------------------------------------------
namespace Zaki::File
{

//==============================================================
/**
 * Parses input files such as CVS, TSV, etc.
 *
 * This class uses CSVRow to read and parse an input file
 * line-by-line given a delimiter.
 * 
 * @see CSVRow
 */
class CSVIterator
{

  //--------------------------------------------------------------
  public:

    /** Default Constructor.
    */
    CSVIterator() ;

    /** Second Constructor.
    *
    * If no delimiter is specified, comma will be assumed (CSV).
    *
    * @param str input file
    * @param in_delim delimiter
    */
    CSVIterator(std::istream& str, const char& in_delim=',') ;

    /** Sets the delimiter.
    *
    * @param in_delim delimiter
    */
    void SetDelim(const char& in_delim);

    typedef std::input_iterator_tag     iterator_category;
    typedef CSVRow                      value_type;
    typedef std::size_t                 difference_type;
    typedef CSVRow*                     pointer;
    typedef CSVRow&                     reference;
    

    /** Pre Increment operator.
    *
    * @param in_delim delimiter
    * @return incremented CVSIterator
    */
    CSVIterator& operator++() ;
    
    /** Post Increment operator.
    *
    * @param n the integer amount to increment
    * @return incremented CVSIterator
    */
    CSVIterator operator++(int n) ;  

    CSVRow const& operator*()   const ;     
    CSVRow const* operator->()  const ;   

    /** Equality comparison operator.
    *
    * @param rhs the other CSVIterator on the R.H.S
    * @return true if the CSVIterators are equal
    */
    bool operator==(CSVIterator const& rhs) ;
    
    /** Inequality comparison operator.
    *
    * @param rhs the other CSVIterator on the R.H.S
    * @return true if the CSVIterators are not equal
    */
    bool operator!=(CSVIterator const& rhs) ;

  //--------------------------------------------------------------

  private:

    /// @brief  m_str variable
    std::istream*       m_str;

    /// @brief  CSVRow variable
    CSVRow              m_row;

};

//--------------------------------------------------------------
} // End of namespace Zaki::File
//--------------------------------------------------------------

//==============================================================
#endif /*Zaki_File_CSVIterator_H*/
