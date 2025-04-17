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

#ifndef Zaki_Util_MemoryManager_H
#define Zaki_Util_MemoryManager_H

#include <iostream>
#include <string>
#include <fstream>
#include <mutex>

#include "Zaki/Util/Logger.hpp"
#include "Zaki/String/String_Basic.hpp"

//--------------------------------------------------------------
namespace Zaki::Util
{

// Forward declration
struct MemEntry ;

//==============================================================
// Singleton log manager class
class MemManager 
{

  public:

    // Destructor
    ~MemManager() ;

    // Copy constructor
    MemManager(const MemManager&) = delete ;

    static MemManager& Get();

    static void New(const MemEntry&);
    static void Delete(MemEntry&);


  private:

    // Constructor 1 with the default settings
    MemManager() { }

    void INew(const MemEntry&) ;
    void IDelete(MemEntry&) ;

    std::mutex mtx;
    static inline std::vector<MemEntry> heap_allocated = {} ; 
    long int total_mem_size = 0 ;

};

//==============================================================
struct MemEntry
{
  friend class MemManager ;

  public:
    enum class Mode
    {
      Delete=0, New
    };

  private:

    void Free()
    {
      Z_LOG_NOTE("Freeing heap allocated memory at: "+ Str() + ".");
      free(ptr);
    }

    void FailMessage()
    {
      Z_LOG_NOTE("Failed attempt to free a pointer at '" +
                  Str() + "' that wasn't allocated.");
    }

    std::string Str() const
    {
      std::stringstream ss;
      ss << ptr;
      return ss.str() ;
    }

    void* ptr ;
    size_t size = 0 ;
    // Mode mem_mode;

  public:
    // Constructor 1
    MemEntry(void* in_ptr, const Mode& in_mode)  
     : ptr(in_ptr)
    {
      if(in_mode == Mode::New)
        MemManager::New(*this) ;
      else
        MemManager::Delete(*this) ;
    }

    // Constructor 2
    MemEntry(void* in_ptr, const Mode& in_mode, const size_t& in_size)  
     : ptr(in_ptr), size(in_size)
    {
      if(in_mode == Mode::New)
        MemManager::New(*this) ;
      else
        MemManager::Delete(*this) ;
    }

    // MemEntry(const MemEntry&) = delete ;
    // MemEntry& operator=(const MemEntry&) = delete ;
    bool operator==(const MemEntry& other) const 
    {
      return ptr == other.ptr ;
    }


} ;

//==============================================================
} // End of namespace Zaki::Util
//==============================================================
//                     Memory Manager

#define Z_NEW_HIDDEN(PTR, SIZE, L)  Zaki::Util::MemEntry __CONCAT(memEntry,L)(PTR, Zaki::Util::MemEntry::Mode::New, SIZE)
#define Z_DELETE_HIDDEN(PTR, L)  Zaki::Util::MemEntry __CONCAT(memEntry,L)(PTR, Zaki::Util::MemEntry::Mode::Delete)

#define Z_NEW(PTR, SIZE)  Z_NEW_HIDDEN(PTR, SIZE, __COUNTER__)
#define Z_DELETE(PTR)  Z_DELETE_HIDDEN(PTR, __COUNTER__)

/* #define Z_NEW(PTR) \
//       do { \
//        Zaki::Util::MemEntry __CONCAT(memEntry,__COUNTER__)(PTR, Zaki::Util::MemEntry::Mode::New, SIZE); \
//        return ; \
//        } while (0)

// #define Z_DELETE(PTR) \
//       do { \
//        Zaki::Util::MemEntry __CONCAT(memEntry,__COUNTER__)(PTR, Zaki::Util::MemEntry::Mode::Delete); \
//        return ; \
//        } while (0)
 */
//==============================================================
#endif /*Zaki_Util_MemoryManager_H*/
