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

#ifndef Zaki_Util_Instrumentor_H
#define Zaki_Util_Instrumentor_H

// Basic instrumentation profiler by Cherno

// Usage: include this header file somewhere in your code (eg. precompiled header), and then use like:
//
// Instrumentor::Get().BeginSession("Session Name");        // Begin session 
// {
//     InstrumentationTimer timer("Profiled Scope Name");   // Place code like this in scopes you'd like to include in profiling
//     // Code
// }
// Instrumentor::Get().EndSession();                        // End Session
//
// You will probably want to macro-fy this, to switch on/off easily and use things like __FUNCSIG__ for the profile name.


// Use "chrome://tracing/" to visualize.

#include <string>
#include <chrono>
#include <fstream>
#include <mutex>
#include <thread>

#include "Zaki/String/Directory.hpp"

//--------------------------------------------------------------
namespace Zaki::Util
{

//==============================================================
struct ProfileResult
{
    std::string Name;
    long long Start, End;
    uint32_t ThreadID;
};

//==============================================================
//                     Singleton design
class Instrumentor
{

    //--------------------------------------------------------------
    private:

    std::string     m_sessionName   = "None";
    std::ofstream   m_OutputStream;
    int             m_ProfileCount = 0;
    std::mutex      m_lock;
    bool            m_activeSession = false;

    Instrumentor() {}

    //--------------------------------------------------------------
    void IBeginSession(const std::string&, const Zaki::String::Directory& = "results.json") ;

    //--------------------------------------------------------------
    void IEndSession() ;

    //--------------------------------------------------------------
    public:

    //--------------------------------------------------------------
    static void BeginSession(const std::string& , const Zaki::String::Directory& = "results.json") ;

    //--------------------------------------------------------------
    ~Instrumentor()
    {
        EndSession();
    }

    //--------------------------------------------------------------
    static void EndSession() ;

    //--------------------------------------------------------------
    void WriteProfile(const ProfileResult&) ;

    //--------------------------------------------------------------
    void WriteHeader() ;

    //--------------------------------------------------------------
    void WriteFooter() ;
    //--------------------------------------------------------------
    static Instrumentor& Get() ;
    //--------------------------------------------------------------
};

//==============================================================
class InstrumentationTimer
{
    //--------------------------------------------------------------
    public:

    InstrumentationTimer(const std::string & name)
            : m_Result({ name, 0, 0, 0 }), m_Stopped(false)
    {
        m_StartTimepoint = std::chrono::high_resolution_clock::now();
    }

    //--------------------------------------------------------------
    ~InstrumentationTimer()
    {
        if (!m_Stopped)
            Stop();
    }

    //--------------------------------------------------------------
    void Stop() ;
    //--------------------------------------------------------------
    private:
        ProfileResult m_Result;
        std::chrono::time_point<std::chrono::high_resolution_clock> m_StartTimepoint;
        bool m_Stopped;
};
//==============================================================

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
#define PROFILING 1
#if PROFILING
    #define PROFILE_SCOPE_HIDDEN(name, L) Zaki::Util::InstrumentationTimer __CONCAT(inst_timer,L)(name)
    #define PROFILE_SCOPE(name) PROFILE_SCOPE_HIDDEN(name, __COUNTER__)

    /* #define PROFILE_SCOPE(name) \
     do { \
        Zaki::Util::InstrumentationTimer __CONCAT(inst_timer, __COUNTER__)(name) ; \
        return ; \
        } while (0)
     */
    #define PROFILE_FUNCTION() PROFILE_SCOPE(Z_PRETTY_FUNCTION)
#else
    #define PROFILE_SCOPE(name)
    #define PROFILE_FUNCTION()
#endif

//==============================================================
#endif /*Zaki_Util_Instrumentor_H*/
