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
 * @file Banner.hpp
 *
 * @brief Creates banners showing the software info in terminal.
 *
 * @ingroup String
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_String_Banner_H
#define Zaki_String_Banner_H

#include "string"
#include "vector"

#include "Zaki/String/TextBox.hpp"

//--------------------------------------------------------------
namespace Zaki::String
{

//==============================================================
  //--------------------------------------------------------------
  struct Content
  {
    size_t index = 0 ;
    
    Content(const size_t& idx) : index(idx) {}

    virtual std::string Str() const = 0;

    virtual ~Content() {}
  };

  //--------------------------------------------------------------
  struct Version : public Content
  {
    std::string num;
    std::string date;

    // constructor
    Version(const char* in_num, const char* in_date, const size_t& idx=0) 
      :Content(idx), num(in_num), date(in_date) {}

    std::string Str() const override
    {
      return "Version " + num + "         " + date ;
    }

  };

  //--------------------------------------------------------------
  struct Website : public Content
  {
    std::string label;
    std::string link ;

    // constructor
    Website(const char* in_lab, const char* in_link, size_t idx=0) 
      : Content(idx), label(in_lab), link(in_link){}

    std::string Str() const override 
    {
      return label + ": " + link ;
    }
  };

  //--------------------------------------------------------------
    struct Author : public Content
  {
    std::string first_name;
    std::string last_name ;

    // constructor
    Author(const char* in_first, const char* in_last, size_t idx=0) 
      : Content(idx), first_name(in_first), last_name(in_last) {}

    std::string Str() const override 
    {
      return "Created by " + ToString(first_name[0]) + ". " + last_name ;
    }
  };

  //--------------------------------------------------------------
  struct ProgramName : public Content
  {
    std::string name;

    // constructor
    ProgramName(const char* in_name, size_t idx=0) 
      : Content(idx), name(in_name){}


    std::string Str() const override 
    {
      return "Welcome to " + name ;
    }
  };

  //--------------------------------------------------------------
    struct Misc : public Content
  {
    std::string item;

    // constructor
    Misc(const char* in_item, size_t idx=0) 
      :Content(idx), item(in_item) { }


    std::string Str() const override 
    {
      return item ;
    }
  };

  //------------------------------------------------------------
  class Banner
  {

  public:
    Banner();
    Banner(ProgramName*, Version*, Author*, Website*);

    // void SetName(const char*) ;
    // void SetAuthor(const char*) ;
    // void SetVersion(const Version) ;
    // void SetWebsite(const Website) ;
    void AddContent(Content*) ;
    void SortContents() ;

    TextBox* GetTextBox() ;
    void Show() ;

  private:
    // std::string name ;
    TextBox text_box ;  
    // Version version ;
    // std::string author ;
    // std::vector<std::string> extra_txt ;
    // Website website ;
    std::vector<Content*> contents;
  };
  //------------------------------------------------------------

//--------------------------------------------------------------
} // End of namespace Zaki::String
//--------------------------------------------------------------

//==============================================================
#endif /*Zaki_String_Banner_H*/