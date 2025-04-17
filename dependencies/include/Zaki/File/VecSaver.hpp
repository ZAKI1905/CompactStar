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
 * @file VecSaver.hpp
 *
 * @brief Exports data into files with a specified delimiter.
 *
 * @ingroup File
 *
 * @author Mohammadreza Zakeri 
 * Contact: M.Zakeri@uky.edu
 *
*/
#ifndef Zaki_File_VecSaver_H
#define Zaki_File_VecSaver_H

#include <vector>
#include "zlib.h"

// Zaki::Util
#include "Zaki/Util/Logger.hpp"

// Zaki::String
#include "Zaki/String/Directory.hpp"

//--------------------------------------------------------------
namespace Zaki::File
{

enum class FileMode
{
  //  r  ,  w  , wb         , a     ,  r+       ,  w+        , a+
  Read=0, Write, WriteBinary, Append, ReadUpdate, WriteUpdate, AppendUpdate
};

//==============================================================
/**
 * Exports data into files with a specified delimiter.
 *
 * This class exports data into a file with a given delimiter. 
 * The data types can be float, double, integer, char, string
 * or of a generic class type with a 'Str()' method implemented. 
 *
 */
class VecSaver
{
  private:

    Zaki::String::Directory file_name{"out.dat"} ;
    FileMode fmode = FileMode::Write;
    std::string header = "";
    std::string footer = "";

    bool set_head_flag = false ;
    bool set_foot_flag = false ;

    // Precision for printing double numbers
    int double_precision = 8 ;

    const char* Mode()
    {
      switch (fmode)
      {
        case FileMode::Read:
          return "r+" ;
          break;
          
        case FileMode::Write:
          return "w+" ;
          break;

        case FileMode::WriteBinary:
          return "wb" ;
          break;

        case FileMode::Append:
          return "a+" ;
          break;

        case FileMode::ReadUpdate:
          return "r+" ;
          break;
          
        case FileMode::WriteUpdate:
          return "w+" ;
          break;
        
        case FileMode::AppendUpdate:
          return "a+" ;
          break;

        default:
          return "";
          break;
      }
    }

  public:

    // Constructor 1
    // VecSaver();

    // Constructor 2
    VecSaver(const Zaki::String::Directory& in_fname) 
      : file_name(in_fname) {}

    // Constructor 3
    VecSaver(const Zaki::String::Directory& in_fname, const FileMode& in_mode)
      : file_name(in_fname), fmode(in_mode) {}

    //..............
    // Setters 
    //..............
    VecSaver& SetFileName(const Zaki::String::Directory& in_fname)
    {
      file_name = in_fname ;
      return *this;
    }
    //..............
    VecSaver& SetMode(const FileMode& in_mode)
    {
      fmode = in_mode ;
      return *this;
    }
    //..............
    VecSaver& SetHeader(const std::string& in_head)
    {
      header = in_head ;
      set_head_flag = true ;

      return *this;
    }
    //..............
    // VecSaver& SetHeader(const char* in_head = nullptr)
    // {
    //   if(in_head)
    //   {  header = in_head ;
    //     set_head_flag = true ;
    //   }
    //   return *this;
    // }
    //..............
    VecSaver& SetFooter(const std::string& in_foot)
    {
      footer = in_foot ;
      set_foot_flag = true ;

      return *this;
    }
    //..............
    // VecSaver& SetFooter(const char* in_foot = nullptr)
    // {
    //   if(in_foot)
    //   {
    //     footer = in_foot ;
    //     set_foot_flag = true ;
    //   }

    //   return *this;
    // }
    //..............

    void SetPrecision(const int& digits)
    {
      double_precision = digits ;
    }

    // Exporting 1D vectors
    template <typename T> void Export1D(const std::vector<T>&) ;

    // Exporting 2D vectors
    template <typename T> void Export2D(const std::vector<std::vector<T>>&,
                                        const char* delim=",") ; 
};

//--------------------------------------------------------------
// Exporting 1D vectors
template <typename T> 
void VecSaver::Export1D(const std::vector<T>& list)
{
  if (list.size() ==0 ) return ;

  // Added on [Apr 20, 2023]
  file_name.ThisFileDir().Create() ;

  std::FILE * out_file ;
  out_file = fopen(file_name.Str().c_str(), Mode()) ;

  if (!out_file) 
  { 
    Z_LOG_ERROR(("File: '" + file_name.Str() + "' didn't open!" + "\n").c_str()) ;
    return ;
  }

  std::string new_line = "\n" ;

  // Binary mode 
  if( fmode == FileMode::WriteBinary)
  {
    std::string out_string ;

    //.................................................
    // Preparing the string output
    //.................................................
    // Header
    if(set_head_flag) out_string += header + "\n" ;

    char tmp_char[1000] ;
    // Body
    for(size_t i=0 ; i < list.size() ; ++i)
    { 
      // Avoid adding newline "\n" to the last line
      if (i == list.size()-1)
      {
        new_line = "" ;
      }
      
      if constexpr(std::is_floating_point_v<T>)
        snprintf(tmp_char, sizeof(tmp_char), ("%.*e"+new_line).c_str(), double_precision, list[i]) ;
      else if constexpr(std::is_same_v<T, int> || std::is_same_v<T, const int>)
        snprintf(tmp_char, sizeof(tmp_char), ("%d"+new_line).c_str(), list[i]) ;
      else if constexpr(std::is_same_v<T, std::string> || std::is_same_v<T, const std::string>)
        snprintf(tmp_char, sizeof(tmp_char), ("%s"+new_line).c_str(), list[i].c_str()) ;
      else if constexpr(std::is_member_function_pointer<decltype(&T::Str)>::value)
        snprintf(tmp_char, sizeof(tmp_char), ("%s"+new_line).c_str(), list[i].Str().c_str()) ;
      else 
      {
        Z_LOG_ERROR("Export failed: No 'Str()' method found.");
        break ;
      }
      out_string += tmp_char ;
    }

    // Footer
    if(set_foot_flag) out_string += "\n" + footer ;
    //.................................................

    //.................................................
    //    Compressing the string output
    //.................................................
    const Byte* source = (Byte*)out_string.c_str()  ;
    unsigned long sizeDataOriginal = out_string.size() ;

    unsigned long sizeDataCompressed  = (sizeDataOriginal * 1.1) + 12 ;
    Byte* dataCompressed = (Byte*)malloc( sizeDataCompressed );
    // int zlib_status = 
    compress(dataCompressed, &sizeDataCompressed, 
                                source, sizeDataOriginal) ;

    //.................................................

    //.................................................
    //    Saving the binary compressed data
    //.................................................
    fwrite( dataCompressed, sizeDataCompressed, 1, out_file );
    //.................................................
    // Memory leak resolved on May 6, 2022:
    free(dataCompressed) ;
  }
  else
  {
    // Header
    if(set_head_flag) fprintf(out_file, "%s\n", header.c_str()) ;

    // Body
    for(size_t i=0 ; i < list.size() ; ++i)
    { 
      // Avoid adding newline "\n" to the last line
      if (i == list.size()-1)
      {
        new_line = "" ;
      }
      
      if constexpr(std::is_floating_point_v<T>)
          fprintf(out_file, ("%.*e"+new_line).c_str(), double_precision, list[i]) ;
      else if constexpr(std::is_same_v<T, int> || std::is_same_v<T, const int>)
          fprintf(out_file, ("%d"+new_line).c_str(), list[i]) ;
      else if constexpr(std::is_same_v<T, std::string> || std::is_same_v<T, const std::string>)
          fprintf(out_file, ("%s"+new_line).c_str(), list[i].c_str()) ;
      else if constexpr(std::is_member_function_pointer<decltype(&T::Str)>::value)
          fprintf(out_file, ("%s"+new_line).c_str(), list[i].Str().c_str()) ;
      else 
        Z_LOG_ERROR("Export failed: No 'Str()' method found.");
    }

    // Footer
    if(set_foot_flag) fprintf(out_file, "\n%s", footer.c_str()) ;
  }
  fclose(out_file) ;
}

//--------------------------------------------------------------
// Exporting 2D vectors
template <typename T> 
void VecSaver::Export2D(const std::vector<std::vector<T>>& list, const char* delim)
{
  if (list.size() ==0 ) return ;

  // Added on [Apr 20, 2023]
  file_name.ThisFileDir().Create() ;
  
  std::FILE * out_file ;
  out_file = fopen(file_name.Str().c_str(),  Mode()) ;

  if (!out_file) 
  { 
    Z_LOG_ERROR(("File: '" + file_name.Str() + "' didn't open!" + "\n").c_str()) ;
    return ;
  }
  
  // Header
  if(set_head_flag) fprintf(out_file, "%s\n", header.c_str()) ;

  std::string new_line = "\n" ;
  // Body
  for(size_t i=0; i < list.size(); ++i)
  { 
    // Avoid adding newline "\n" to the last line
    if (i == list.size()-1)
    {
      new_line = "" ;
    }
      
    // We want delimiter after all elements except the
    //  last element of each line:
    for(size_t j=0; j < list[i].size()-1; ++j)
    {
      if constexpr(std::is_floating_point_v<T>)
          fprintf(out_file, "%.*e%s", double_precision, list[i][j], delim) ;
      else if constexpr(std::is_same_v<T, int> || std::is_same_v<T, const int>)
          fprintf(out_file, "%d%s", list[i][j], delim) ;
      else if constexpr(std::is_same_v<T, std::string> || std::is_same_v<T, const std::string>)
          fprintf(out_file, "%s%s", list[i][j].c_str(), delim) ;
      else if constexpr(std::is_member_function_pointer<decltype(&T::Str)>::value)
          fprintf(out_file, "%s%s", list[i][j].Str().c_str(), delim) ;
      else 
        Z_LOG_ERROR("Export failed: No 'Str()' method found.");
    }

    // Last element
    size_t j = list[i].size()-1 ;
    if constexpr(std::is_floating_point_v<T>)
        fprintf(out_file, "%.*e", double_precision, list[i][j]) ;
    else if constexpr(std::is_same_v<T, int> || std::is_same_v<T, const int>)
        fprintf(out_file, "%d", list[i][j]) ;
    else if constexpr(std::is_same_v<T, std::string> || std::is_same_v<T, const std::string>)
        fprintf(out_file, "%s", list[i][j].c_str()) ;
    else if constexpr(std::is_member_function_pointer<decltype(&T::Str)>::value)
        fprintf(out_file, "%s", list[i][j].Str().c_str()) ;
    else
      Z_LOG_ERROR("Export failed: No 'Str()' method found.");

    fprintf(out_file, new_line.c_str());
  }

  // Footer
  if(set_foot_flag) fprintf(out_file, "\n%s", footer.c_str()) ;

  fclose(out_file) ;
}

//--------------------------------------------------------------

//==============================================================
// Records the input quantity (int, float or double) and saves 
// it in a text file.
// template <class T> 
// void SaveVec(const std::vector<T>& list, 
//                       const std::string& out_file_str, const char* mode = "w")
// {
//     if (list.size() ==0 ) return ;

//     std::FILE * out_file ;
//     out_file = fopen((out_file_str + ".dat").c_str(), mode) ;

//     if (!out_file) 
//     { 
//     Z_LOG_ERROR(("File: '" + out_file_str + ".dat' didn't open!" + "\n").c_str()) ;
//     }

//     std::string tmp = "";

//     for(size_t i=0 ; i < list.size() ; ++i)
//     { 
//         if constexpr(std::is_floating_point<T>)
//             fprintf(out_file, "%.4e\n", list[i]) ;
//         else if constexpr(std::is_same_v<T, int>)
//             fprintf(out_file, "%d\n", list[i]) ;
//         else if constexpr(std::is_same_v<T, std::string>)
//             fprintf(out_file, "%s\n", list[i].c_str()) ;
//         else
//             fprintf(out_file, "%s\n", list[i].Str().c_str()) ;
//     }

//   fclose(out_file) ;
// }

//==============================================================
// template <class T> 
// void SaveVec
// ( const std::vector<std::vector<T>>& list, const std::string& out_file_str,
//   const std::string& delim, const char* mode= "w")
// {
//   if (list.size() ==0 ) return ;

//   std::FILE * out_file ;
//   out_file = fopen(out_file_str.c_str(), mode) ;

//   for(size_t i=0; i < list.size(); ++i)
//   { 
//     for(size_t j=0; j < list[i].size(); ++j)
//     {
//       // If we want CSV file, and not the last element of the line
//       if(j != list[i].size()-1 )
//       {
//         if constexpr(std::is_floating_point<T>)
//             fprintf(out_file, "%.4e%s", list[i][j], delim.c_str()) ;
//         else if constexpr(std::is_same_v<T, int>)
//             fprintf(out_file, "%d%s", list[i][j], delim.c_str()) ;
//         else if constexpr(std::is_same_v<T, std::string>)
//             fprintf(out_file, "%s%s", list[i][j].c_str(), delim.c_str()) ;
//         else
//             fprintf(out_file, "%s%s", list[i][j].Str().c_str(), delim.c_str()) ;
//       }
//         // fprintf(out_file,"%.4e%s", list[i][j], delim.c_str() ) ; 
//       else
//       {
//         if constexpr(std::is_floating_point<T>)
//             fprintf(out_file, "%.4e", list[i][j]) ;
//         else if constexpr(std::is_same_v<T, int>)
//             fprintf(out_file, "%d", list[i][j]) ;
//         else if constexpr(std::is_same_v<T, std::string>)
//             fprintf(out_file, "%s", list[i][j].c_str()) ;
//         else
//             fprintf(out_file, "%s", list[i][j].Str().c_str()) ;
//       }
//         // fprintf(out_file,"%.4e", list[i][j]) ; 
//     }
//     fprintf(out_file,"\n");
//   }

//   fclose(out_file) ;
// }

//==============================================================

//--------------------------------------------------------------
} // End of namespace Zaki::File
//--------------------------------------------------------------

//==============================================================
#endif /*Zaki_File_VecSaver_H*/
