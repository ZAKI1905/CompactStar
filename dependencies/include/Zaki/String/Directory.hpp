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
 * @file Directory.hpp
 *
 * @brief For operations with directory paths.
 *
 * @ingroup String
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_String_Directory_H
#define Zaki_String_Directory_H

#include "string"
//--------------------------------------------------------------
namespace Zaki::String
{

//==============================================================

class Directory
{

  public:
    Directory(const std::string& in_path) 
      : full_path(in_path)
    {
      // if(*(full_path.end()-1) == '/')
      // {
      //   // strip the extra backslash from the end
      //   full_path = full_path.substr(0, full_path.find_last_of("/")) ;
      // }

      // if(*full_path.begin() == '/')
      // {
      //   // strip the extra backslash from the beginning
      //   full_path = full_path.substr(1) ;
      // }
    };

    // Overloading for const char* input
    Directory(const char* in_path) 
      : full_path(in_path)
    {
      // if(*(full_path.end()-1) == '/')
      //   // strip the extra backslash from the end
      //   full_path = full_path.substr(0, full_path.find_last_of("/")) ;

      // if(*full_path.begin() == '/')
      // {
      //   // strip the extra backslash from the beginning:
      //   full_path = full_path.substr(1) ;
      // }
    };

    // Setters
    void SetPath(const std::string&) ;

    // Getters
    /// Returns the directory to the full path
    static Directory ThisFileDir(const char*) ;

    /// Returns the file name from the full path
    static Directory ThisFile(const char*) ;

    // Non-static versions
    /// Returns the directory to the full path
    Directory ThisFileDir() const ;

    /// Returns the file name from the full path
    Directory ThisFile() const ;

    Directory NoExt() const ;
    Directory ParentDir() const ;
    std::string GetPath() const ;
    std::string Str() const ;
    void Print() const ;

    Directory operator+(const std::string&) const;
    Directory operator+(const char*) const;
    Directory operator+(const Directory&) const;
    // Overloading postfix -- operator
    Directory operator--(int) ;

    /// @brief Creates the directory
    void Create() const ;

  private:
    std::string full_path = "";

};
//------------------------------------------------------------
//==============================================================
//                  ostream << overloading
std::ostream& operator << (std::ostream &output, const Zaki::String::Directory&) ;

//==============================================================
} // End of namespace Zaki::String
//--------------------------------------------------------------

//==============================================================
#endif /*Zaki_String_Directory_H*/