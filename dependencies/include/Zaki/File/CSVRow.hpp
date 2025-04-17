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
 * @file CSVRow.hpp
 *
 * @brief A row of strings.
 *
 * @ingroup File
 * 
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
#ifndef Zaki_File_CSVRow_H
#define Zaki_File_CSVRow_H

#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

//--------------------------------------------------------------
namespace Zaki::File
{

//==============================================================
/**
 * Contains a list of strings with a delimeter.
 *
 * This class defines rows in an imported file and is used by 
 * CSVIterator for reading and parsing an input file line-by-line
 * given a delimiter.
 * 
 * @see CSVIterator
 */
class CSVRow
{

  //--------------------------------------------------------------
  public:

    CSVRow() ;

    std::string const& operator[](std::size_t index) const ;

    std::size_t size() const ;
      
    void readNextRow(std::istream& str) ;
    void SetDelim(const char& in_delim) ;

  //--------------------------------------------------------------
  private:
    std::vector<std::string>    m_data;
    char delim = ',' ;
};

//--------------------------------------------------------------
} // End of namespace Zaki::File
//--------------------------------------------------------------

std::istream& operator>>(std::istream& str, Zaki::File::CSVRow& data) ;
//==============================================================
#endif /*Zaki_File_CSVRow_H*/
