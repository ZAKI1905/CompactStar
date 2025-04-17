#ifndef Base_H
#define Base_H

#include "Confind/ConfindConfig.h"

#include <Zaki/Util/Instrumentor.hpp>
#include <Zaki/String/Directory.hpp>

#include "Confind/Common.hpp"

#define CONFIND_BASE_DEBUG_MODE 1

namespace CONFIND
{

//==============================================================
class Base {

  //--------------------------------------------------------------
  protected:
    
    /// Prints the banned of the program
    void ShowBanner();

    /// The number of derived object instances
    static std::atomic<size_t> counter;

    /// The work directory
    Zaki::String::Directory wrk_dir{""}    ;
    
    /// The name of the object
    std::string name = ""    ;

    // Flags
    
    /// Returns if the name is set
    bool set_name_flag = false    ;
    
    /// Returns if the work directory is set
    bool set_wrk_dir_flag = false ;
    
    static inline std::string class_name = "Base" ;
    
  //--------------------------------------------------------------
  private:
    
    /// Track object creation and destruction
    bool track_objs_flag = false ;
  //--------------------------------------------------------------
  public:
    
    /// Default Constructor for faster object creation
    ///  -> No name assignment or object tracking
    Base() ;

    /// Constructor via a name and no tracking options
    /// -> Name is set but no object tracking
    Base(const std::string&) ;

    /// Constructor via a name and tracking options
    /// This is slower and shouldn't be used for objects
    /// that are created thousands of time, i.e. particles, etc.
    /// -> Name is set and object tracking is on
    Base(const std::string&, const bool&) ;

    /// Destructor
    virtual ~Base() ;

    /// Overloading CLass specific new operator
    static void* operator new(size_t sz) ;

    /// Overloading CLass specific delete operator
    static void operator delete(void* m);

    /// Sets the work directory
    virtual void SetWrkDir(const Zaki::String::Directory&) ;
    
    /// Sets the name
    virtual void SetName(const std::string&) ;

    /// Returns the work directory
    virtual Zaki::String::Directory GetWrkDir() const ;
    
    /// Returns the name of the class
    virtual std::string GetName() const   ;
    
    /// Print method for the derived classes
    virtual void Print() const            ;
    
    /// Returns the pointer address in  a string form
    std::string PtrStr() const ;
};

//==============================================================
} // Namespace "CONFIND" ends.
//==============================================================
#endif /*Base_H*/
