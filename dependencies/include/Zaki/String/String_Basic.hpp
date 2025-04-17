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
 * @file String_Basic.hpp
 *
 * @brief Basic string functions.
 *
 * @ingroup String
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_String_Basic_H
#define Zaki_String_Basic_H

#include "string"
#include "sstream"
#include "vector"

//--------------------------------------------------------------
namespace Zaki::String
{

//==============================================================
template <class T>
std::string ToString(const T& a) 
{
  std::stringstream buff;
  buff << a;
  return buff.str();
}

//==============================================================
// Checks if a string ends in another string
bool EndsWith(std::string const &fullString, std::string const &ending) ;

//==============================================================
// returns a copy of str by stripping character c from it
std::string Strip(const std::string& str, char c) ;

//==============================================================
// Parses the input based on the given delim
std::vector<std::string> Pars(const std::string&, const char*, int=1) ;
// Parses the input based on the given delim list
std::vector<std::string> Pars(const std::string&, const std::vector<std::string>&) ;

//==============================================================
// converts a string list of numbers separated by delim (default = ",")
// to a vector of numbers
void Str2Lst(std::vector<int>&, const std::string&, const char* =",") ;
// Overloading for float numbers
void Str2Lst(std::vector<float>&, const std::string&, const char* =",") ;
// overloading for string lists
void Str2Lst(std::vector<std::string>&, const std::string&, const char* =",") ;

//==============================================================
// Integer range input parser for inputs like:
//  (...=1-10,50-75,300,401)
std::vector<int> RangeParser(const std::string&) ;

//==============================================================
// Copying a single character multiple (n) times
std::string Multiply(const char, size_t) ;

// Copying a character multiple (n) times
std::string Multiply(const char*, size_t) ;

//==============================================================
enum class FGColor
{
  Black, Red, Green, Yellow, Blue, Purple, Cyan, White,
  LBlack, LRed, LGreen, LYellow, LBlue, LPurple, LCyan, LWhite
};

//==============================================================
enum class BGColor
{
  BlackBg, RedBg, GreenBg, YellowBg, BlueBg, PurpleBg, CyanBg, WhiteBg
};

//==============================================================
struct Color
{
  FGColor Foreground;
  BGColor Background;

  // Constructor 1:
  Color(const FGColor& fg, const BGColor& bg) 
    : Foreground(fg), Background(bg) {}

  // Constructor 2:
  Color(const FGColor& fg) 
    : Foreground(fg), Background(BGColor::BlackBg) {}

  // Constructor 3:
  Color() 
    : Foreground(FGColor::White), Background(BGColor::BlackBg) {}
  
  // String form of the color for '<<'
  std::string Str() const ;

  private:
    std::string str_form = "" ;
};

//==============================================================
// std::string ColorDict(Color) ;

//==============================================================
int FindOccurrences(const std::string&, const std::string&) ;

//==============================================================
//                  ostream << overloading
std::ostream& operator << (std::ostream &output, const Zaki::String::Color&) ;

//==============================================================

//--------------------------------------------------------------
} // End of namespace Zaki::String
//--------------------------------------------------------------

//==============================================================
#endif /*Zaki_String_Basic_H*/