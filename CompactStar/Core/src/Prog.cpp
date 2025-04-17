/*
  Prog class
*/

// Creating directory
#include <sys/stat.h>

#include <Zaki/String/Banner.hpp>
#include <Zaki/Util/MemoryManager.hpp>
#include <Zaki/Util/ObjObserver.hpp>

#include "CompactStar/Core/Prog.hpp"

using namespace CompactStar ;
//==============================================================
//std::atomic<size_t> Prog::counter = 0 ;
//==============================================================

//--------------------------------------------------------------
/// Default Constructor for faster object creation
///  -> No name assignment or object tracking
Prog::Prog()
{
  if(counter==0)
    ShowBanner() ;

  // if(track_objs_flag)
  //   Z_OBJ_CTR(this, name) ;

  counter++ ;
}

//--------------------------------------------------------------
/// Constructor via a name and no tracking options
/// -> Name is set but no object tracking
Prog::Prog(const std::string& in_name)
{
  name = in_name ;
  set_name_flag = true ;

  if(counter==0)
    ShowBanner() ;

  counter++ ;
}

//--------------------------------------------------------------
/// Constructor via a name and tracking options
/// This is slower and shouldn't be used for objects 
/// that are created thousands of time, i.e. particles, etc.
/// -> Name is set and object tracking is on
/// 
/// NOTE: CompactStar_PROG_DEBUG_MODE macro overrides everything
///  It's useful for the fastest runtime, because then 
///  tracking is not checked in the destructor
Prog::Prog(const std::string& in_name, const bool& tracking)
{
  name = in_name ;
  set_name_flag = true ;
    
#if CompactStar_PROG_DEBUG_MODE
    track_objs_flag = tracking ;
  if(track_objs_flag)
    Z_OBJ_CTR(this, name) ;
#endif

  if(counter==0)
    ShowBanner() ;

  counter++ ;
}

//--------------------------------------------------------------
Prog::~Prog()
{

#if CompactStar_PROG_DEBUG_MODE
std::cout << "\n\t~Prog() called!\n";
  if(track_objs_flag)
    Z_OBJ_DTR(this, name) ;
#endif

}

//--------------------------------------------------------------
void Prog::ShowBanner()
{
  using namespace Zaki::String ;

  Banner banner ;

  ProgramName p_name("CompactStar", 1) ;
  banner.AddContent(&p_name) ;

  Author auth("Mohammadreza", "Zakeri", 4) ;
  banner.AddContent(&auth) ;

  Version ver(CompactStar_VERSION_STR, CompactStar_RELEASE_DATE, 2);
  banner.AddContent(&ver) ;

  Website web("GitHub", "github.com/ZAKI1905/CompactStar", 5) ;
  banner.AddContent(&web) ;

  // Misc misc("\U0001f323  \u2192  \U0001f30D", 3); 
  // banner.AddContent(&misc) ;

  banner.GetTextBox()->SetTextColor({FGColor::LCyan, BGColor::BlackBg}) ;
  banner.GetTextBox()->SetFrameColor({FGColor::LYellow, BGColor::BlackBg}) ;
  banner.GetTextBox()->SetPadColor({FGColor::LCyan, BGColor::BlackBg}) ;

  banner.GetTextBox()->SetAlignment(TextBox::center) ;
  banner.GetTextBox()->SetPadding(5) ;

  banner.GetTextBox()->EnableClearScreen() ;

  banner.Show() ;

}
//--------------------------------------------------------------
Prog* Prog::SetWrkDir(const Zaki::String::Directory& input) 
{
  wrk_dir = input;
  set_wrk_dir_flag = true ;

  input.Create() ;
  // // ............ Creating a directory ............
  // if (mkdir(wrk_dir.Str().c_str(), ACCESSPERMS) == -1) 
  // {
  //   Z_LOG_NOTE("Directory '"+input.Str()+"' wasn't created, because: "+strerror(errno)+".") ;
  // }
  // else
  //   Z_LOG_INFO(("Directory '" + wrk_dir.Str() + "' created.").c_str()); 
  // // .................................................

  Z_LOG_INFO("Work directory set to '"+input.Str()+"'.") ;

  SetMemWrkDir(input) ;

  return this ;
}

//--------------------------------------------------------------
Prog* Prog::SetMemWrkDir(const Zaki::String::Directory& input) 
{
  return this ;
}

//--------------------------------------------------------------
void Prog::SetName(const std::string& input) 
{
  name          = input;
  set_name_flag = true ;
  Z_LOG_INFO("Name set to '"+input+"'.") ;
}

//--------------------------------------------------------------
std::string Prog::GetName() const 
{
  if (!set_name_flag)
    Z_LOG_ERROR("Name not set!") ;
  return name;
}

//--------------------------------------------------------------
Zaki::String::Directory Prog::GetWrkDir() const 
{
  if (!set_wrk_dir_flag)
    Z_LOG_ERROR("Work directory not set!") ;

  return wrk_dir;
}
//--------------------------------------------------------------
/// Returns set_wrk_dir_flag
bool Prog::IsWrkDirSet() const 
{
  return set_wrk_dir_flag ;
}

//--------------------------------------------------------------
void Prog::Print() const
{
  std::cout<<"\n --------------------------------------------------------------\n" ;
  std::cout<<"|         Name: "<<GetName() << "\n";
  std::cout<<" --------------------------------------------------------------\n" ;
}

//--------------------------------------------------------------
// Retrns the pointer address in string form
std::string Prog::PtrStr() const 
{
  std::stringstream ss;
  ss << static_cast<const void*>(this);

  return ss.str() ;
}


//--------------------------------------------------------------
// Overloading Class specific new operator
void* Prog::operator new(size_t sz)
{
  // if > 10 KBytes send a note to the user
  if(sz > 10000)
    Z_LOG_NOTE(("Large size of memory requested: "
                + std::to_string(sz) + " bytes.").c_str()) ;

  void* m = malloc(sz);
  
  Z_NEW(m, sz) ;

  return m;
}

//--------------------------------------------------------------
// Overloading CLass specific delete operator
void Prog::operator delete(void* m)
{
  Z_DELETE(m) ;
}

//--------------------------------------------------------------

//==============================================================
