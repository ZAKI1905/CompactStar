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

#ifndef Zaki_Util_Simple_Timer_H
#define Zaki_Util_Simple_Timer_H

#include <string>
#include <chrono>
// #include <algorithm>
#include <fstream>
#include <mutex>
#include <thread>

#include "Zaki/String/Directory.hpp"

//--------------------------------------------------------------
namespace Zaki::Util
{

//==============================================================
struct Operation
{
    std::string name ;
    double duration ;

    Operation(const std::string& in_name, double in_dur) 
        : name(in_name), duration(in_dur) {} 
};

//==============================================================
//                     Singleton design
class TimeManager
{

    //--------------------------------------------------------------
    private:

    std::string     m_sessionName   = "None";
    std::ofstream   m_OutputStream;
    int             m_ProfileCount = 0;
    std::mutex      m_lock;
    bool            m_activeSession = false;

    //--------------------------------------------------------------
    TimeManager() {}

    //--------------------------------------------------------------
    void IEndSession() ;

    //--------------------------------------------------------------
    void IBeginSession(const std::string&, const Zaki::String::Directory& = "Timer_Results.txt") ;

    //--------------------------------------------------------------
    public:

    //--------------------------------------------------------------
    ~TimeManager()
    {
        EndSession();
    }

    //--------------------------------------------------------------
    static void BeginSession(const std::string&, const Zaki::String::Directory& ="Timer_Results.txt") ;
    //--------------------------------------------------------------
    static void EndSession() ;

    //--------------------------------------------------------------
    void WriteProfile(const Operation& oper);

    //--------------------------------------------------------------
    void WriteHeader() ;

    //--------------------------------------------------------------
    void WriteFooter() ;

    //--------------------------------------------------------------
    static TimeManager& Get() ;
    //--------------------------------------------------------------
};

//==============================================================
class Timer
{
    //--------------------------------------------------------------
    private:

    Operation op ;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTimePt ;
    bool stopped = false ;

    //--------------------------------------------------------------
    public:

    //--------------------------------------------------------------
    Timer(const std::string& in_scope) : op({in_scope, 0})
    {
        startTimePt = std::chrono::high_resolution_clock::now() ;
    }
    //--------------------------------------------------------------
    ~Timer()
    {
        if(!stopped) Stop() ; 
    }

    //--------------------------------------------------------------
    void Stop() ;
};

//==============================================================
//                  ostream << overloading
std::ostream& operator << (std::ostream &output, const Operation& o) ;

//--------------------------------------------------------------
} // End of namespace Zaki::Util
//==============================================================
//                    Profiling Macros
//==============================================================

#if defined(_MSC_VER)
#define Z_PRETTY_FUNCTION __FUNCSIG__
#else
#define Z_PRETTY_FUNCTION __PRETTY_FUNCTION__
#endif
//--------------------------------------------------------------
#define TIMING 1
#if TIMING
    #define Z_TIMER_SCOPE_HIDDEN(name, L) Zaki::Util::Timer __CONCAT(sim_timer,L)(name)
    #define Z_TIMER_SCOPE(name) Z_TIMER_SCOPE_HIDDEN(name, __COUNTER__)

/*
    // #define Z_TIMER_SCOPE(name) \
    // do { \
    //     Zaki::Util::Timer __CONCAT(sim_timer,__COUNTER__)(name); \
    //     return ; \
    // } while (0)
*/
    #define Z_TIMER() Z_TIMER_SCOPE(Z_PRETTY_FUNCTION)
#else
    #define Z_TIMER_SCOPE(name)
    #define Z_TIMER()
#endif

//==============================================================
#endif /*Zaki_Util_Simple_Timer_H*/
