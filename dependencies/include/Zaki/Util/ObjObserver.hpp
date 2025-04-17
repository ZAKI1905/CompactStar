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

#ifndef Zaki_Util_ObjObserver_H
#define Zaki_Util_ObjObserver_H

#include <iostream>
#include <string>
#include <chrono>
#include <sstream>
#include <vector>
#include <mutex>
#include <memory> // std::unique_ptr

#include "Zaki/Util/Logger.hpp"
#include "Zaki/String/Directory.hpp"

//--------------------------------------------------------------
namespace Zaki::Util
{

//==============================================================
class ObjEvent
{
  protected:
    virtual ObjEvent* IClone() const = 0 ;
    long long obj_ev_t = 0 ;

  public:
    ObjEvent() {} ;

    virtual ~ObjEvent() {} ;

    std::unique_ptr<ObjEvent> Clone() const 
    { return std::unique_ptr<ObjEvent>(IClone()); }

    void SetTime(const long long& in_t) { obj_ev_t = in_t ;}
    long long GetTime() const { return obj_ev_t;} 

    virtual const std::string Str() const = 0 ;

};

//==============================================================
class ObjCreate : public virtual ObjEvent
{
  protected:
    virtual ObjCreate* IClone() const = 0 ;

  public:
    ObjCreate() {} ;
    virtual ~ObjCreate() {} ;
    // virtual const std::string Str() const = 0 ;
};

//==============================================================
class ObjBiEvent : public virtual ObjEvent
{
  friend class ObjManager ;

  protected:
    virtual ObjBiEvent* IClone() const = 0 ;
    std::unique_ptr<ObjBiEvent> GetReaction() 
    { return std::unique_ptr<ObjBiEvent> (IGetReaction()) ; }

    virtual ObjBiEvent* IGetReaction() const = 0 ;

    void* from = nullptr ;
    void* to = nullptr ;
    std::string name ;


  public:
    ObjBiEvent(void* in_fro, void* in_to, const std::string& in_name)
    : from(in_fro) , to(in_to), name(in_name) {}

    virtual ~ObjBiEvent() {}

    void* GetTo() {return to;}
    void* GetFrom() {return from;}

    // virtual const std::string Str() const = 0 ;
};

//==============================================================
class ObjConstruct : public ObjCreate
{
  protected:
    virtual ObjConstruct* IClone() const override
    {
      return (new ObjConstruct(*this)) ;
    }

  public:
    ObjConstruct() {};
    ~ObjConstruct() {};
    virtual const std::string Str() const override { return "CTR" ;}
};

//==============================================================
// Fixing the forward declration issue with covariance
class IObjCopiedFrom : public ObjBiEvent
{
  private:
    // virtual IObjCopiedTo* IClone() const override 
    // {
    //   return (new IObjCopiedTo(*this)) ;
    // }
    
    // virtual IObjCopiedTo* IGetReaction() const override {} 

    // virtual const std::string Str() const override 
    // { 
    //   return  "";
    // }

  public:
    IObjCopiedFrom(void* in_fr, void* in_to, const std::string& in_name)
      : ObjBiEvent(in_fr, in_to, in_name) {}

};

//==============================================================
class ObjCopyConstruct : public ObjCreate, public ObjBiEvent
{
  protected:
    virtual ObjCopyConstruct* IClone() const override 
    {
      return (new ObjCopyConstruct(*this)) ;
    }
    
    virtual IObjCopiedFrom* IGetReaction() const override ;
    // {
    //   return (new ObjCopiedTo(other_ptr)) ;
    // }

  public:
    ObjCopyConstruct(void* in_fro, void* in_to, const std::string& in_name) 
    : ObjBiEvent(in_fro, in_to, in_name) 
    { 
      // std::cout << "\n From: "<<from << ", to: "<< to << ".\n";
    };
    
    ~ObjCopyConstruct() {};

    virtual const std::string Str() const override 
    { 
      std::stringstream ss ;
      ss << "CFrom, " << name << "_" << from ;
      return  ss.str();
    }
};

//==============================================================
class ObjCopiedFrom : public ObjCreate, public IObjCopiedFrom
{
  protected:
    virtual ObjCopiedFrom* IClone() const override 
    {
      return (new ObjCopiedFrom(*this)) ;
    }

    virtual ObjCopyConstruct* IGetReaction() const override ; 
    // {
    //   return (new ObjCopyConstruct(other_ptr)) ;
    // }

  public:
    ObjCopiedFrom(void* in_fro, void* in_to, const std::string in_to_name) 
    : IObjCopiedFrom(in_fro, in_to, in_to_name) { };
    
    ~ObjCopiedFrom() {};

    virtual const std::string Str() const override 
    { 
      std::stringstream ss ;
      ss << "CTo, " << name << "_" << to ;
      return  ss.str();
    }
};
//==============================================================
class ObjDestruct : public ObjEvent
{
  protected:
    virtual ObjDestruct* IClone() const override 
    {
      return (new ObjDestruct(*this)) ;
    }

  public:
    ObjDestruct() {};
    ~ObjDestruct() {} ;
    virtual const std::string Str() const override { return "DTR" ;};
};

//==============================================================
// Forward declration
class ObjEntry ;
struct ObjEventEntry ;

/// Singleton object manager class for keeping track of the status of objects
class ObjManager 
{

  public:

    // Destructor
    ~ObjManager() ;

    // Copy constructor
    ObjManager(const ObjManager&) = delete ;

    static ObjManager& Get();

    static void Emit(ObjEventEntry&);
    static void SetFile(const Zaki::String::Directory&) ;


  private:

    // Constructor 1 with the default settings
    ObjManager() { }
    std::string CurrentDateTime(const std::string& format) const ;

    void IEmit(ObjEventEntry&) ;
    void SetTimeOffset(const long long&) ;
    std::string GetFileHeader() ;
    // Adding a footer to the file
    void AddFooter() ;

    void Export(const ObjEntry&) ;
    Zaki::String::Directory file_name = "";

    long int time_offset = 0 ;
    bool first_time_writing = true ;
    bool first_time_emitting = true ; 
    bool export_flag = false ;

    std::mutex mtx;
    
    static inline std::vector<ObjEntry> obj_list = {} ;

};

//==============================================================
class ObjEntry
{
  friend class ObjManager ;

  private:
    void* ptr = nullptr ;
    std::string name ;
    std::vector<std::unique_ptr<ObjEvent>> events ;

  public:
    ObjEntry(void* in_ptr, const std::string& in_name, 
      std::unique_ptr<ObjEvent>& creation_ev) 
      : ptr(in_ptr), name(in_name)
      {
        events.emplace_back(creation_ev->Clone()) ;
      }

    void AddEvent(ObjEvent* gen_ev)
    {
      events.emplace_back(gen_ev->Clone()) ;
    }

  std::string Str() const
  {
    std::stringstream ss;
    ss << name << "_" << ptr;
    return ss.str() ;
  }

  // for std::find
  bool operator==(const ObjEntry& other) const 
  {
    return ptr == other.ptr ;
  }

  bool operator==(void* other) const 
  {
    return ptr == other ;
  }
};

//==============================================================
struct ObjEventEntry
{
  friend class ObjManager ;
  friend class ObjEntry   ;

  private:

    // void Stop() ;

    void* ptr ;
    std::string name ;
    std::unique_ptr<ObjEvent> obj_ev_ptr ;
    // long long obj_ev_t ;
    
    // bool stopped = false ;

  public:
    // Constructor 1
    ObjEventEntry(void* in_ptr, const std::string& in_name, ObjEvent* in_mode)  
     : ptr(in_ptr), name(in_name)
    {
      obj_ev_ptr = std::unique_ptr<ObjEvent>(in_mode) ;
      // Z_LOG_INFO("ObjEventEntry constructor!") ;
      // Setting the objEntry creation time
      obj_ev_ptr->SetTime(std::chrono::time_point_cast<std::chrono::microseconds>
                        (std::chrono::high_resolution_clock::now()).time_since_epoch().count() ) ;

      // if(in_mode == Mode::Create)
      ObjManager::Emit(*this) ;
      // else
      //   ObjManager::Destruct(*this) ;
    }

    // destructor
    // ~ObjEntry()
    // {
    //   if(!stopped) Stop() ;
    // }

    // ObjEntry(const ObjEntry&) = delete ;
    // ObjEntry& operator=(const ObjEntry&) = delete ;
    // bool operator==(const ObjEventEntry& other) const 
    // {
    //   return ptr == other.ptr ;
    // }


} ;

//==============================================================
} // End of namespace Zaki::Util
//==============================================================
//                     Object Manager

#define Z_OBJ_CTR_HIDDEN(PTR, N, L)  Zaki::Util::ObjEventEntry __CONCAT(objEvEntry,L)(PTR, N, new Zaki::Util::ObjConstruct)
#define Z_OBJ_CCTR_HIDDEN(PTR, PTR2, N1, N2, L)  Zaki::Util::ObjEventEntry __CONCAT(objEvEntry,L)(PTR, N1, new Zaki::Util::ObjCopyConstruct(PTR2, PTR, N2))
#define Z_OBJ_DTR_HIDDEN(PTR, N, L)  Zaki::Util::ObjEventEntry __CONCAT(objEvEntry,L)(PTR, N, new Zaki::Util::ObjDestruct)

#define Z_OBJ_CTR(PTR, N)          Z_OBJ_CTR_HIDDEN(PTR, N, __COUNTER__)
#define Z_OBJ_CCTR(PTR, PTR2, N1, N2)   Z_OBJ_CCTR_HIDDEN(PTR, PTR2, N1, N2, __COUNTER__)
#define Z_OBJ_DTR(PTR, N)          Z_OBJ_DTR_HIDDEN(PTR, N, __COUNTER__)

/*
// #define Z_OBJ_CTR(PTR, N) \
//   do { \
//     Zaki::Util::ObjEventEntry __CONCAT(objEvEntry,__COUNTER__)(PTR, N, new Zaki::Util::ObjConstruct); \
//     return ; \
//   } while (0)

// #define Z_OBJ_CCTR(PTR, PTR2, N1, N2) \
//   do { \
//     Zaki::Util::ObjEventEntry __CONCAT(objEvEntry,__COUNTER__)(PTR, N1, new Zaki::Util::ObjCopyConstruct(PTR2, PTR, N2)); \
//     return ; \
//   } while (0)

// #define Z_OBJ_DTR(PTR, N) \
//   do { \
//     Zaki::Util::ObjEventEntry __CONCAT(objEvEntry,__COUNTER__)(PTR, N, new Zaki::Util::ObjDestruct); \
//     return ; \
//   } while (0)
*/
//==============================================================
#endif /*Zaki_Util_ObjObserver_H*/
