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
 * @file DateTime.hpp
 *
 * @brief Clocks and Date-time definitions and conversions.
 *
 * @ingroup Physics
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/


#ifndef Zaki_Physics_DateTime_H
#define Zaki_Physics_DateTime_H

#include <iostream>

//==============================================================
namespace Zaki::Physics
{

//==============================================================
// Time
struct Clock
{
    int h       ;
    int min     ;
    float sec  ;

    //..........................
    // Constructor 0
    Clock() {} 

    // Constructor 1
    Clock(const int& h_in, const int& min_in, const float& sec_in) : 
    h(h_in), min(min_in), sec(sec_in) {};

    // Constructor 2
    Clock(const double& hour_f) 
    {
      inputHour(hour_f) ;
    }

    //..........................

    void inputHour(const double& hour_f) ;

    double getHours() const
    {  return h + min * 1.0 / 60 + sec * 1.0 / 3600 ; }

    //..................................................
    // Operator overloading
    // +
    Clock operator+(const Clock& c) const
    {
        Clock out_c ;
        out_c.inputHour(getHours() + c.getHours()) ;
        return out_c ;
    }

    // +
    Clock operator+(const double& in_hour) const
    {
        Clock out_c ;
        out_c.inputHour(getHours() + in_hour) ;
        return out_c ;
    }
    //.............
    // -
    Clock operator-(const Clock& c) const
    {
        Clock out_c ;
        out_c.inputHour(getHours() - c.getHours()) ;
        return out_c ;
    }
    // -
    Clock operator-(const double& in_hour) const
    {
        Clock out_c ;
        out_c.inputHour(getHours() - in_hour) ;
        return out_c ;
    }
    //.............
    // ==
    bool operator==(const Clock& c) const
    {
        return (h == c.h && min == c.min && sec == c.sec) ;
    }
    // ==
    bool operator==(const double& in_hour) const
    {
        return in_hour == getHours();
    }
    //.............
    // >
    bool operator>(const Clock& c) const
    { return getHours() > c.getHours(); }
    // >
    bool operator>(const double& in_hour) const
    { return getHours() > in_hour; }
    //.............
    // <
    bool operator<(const Clock& c) const
    { return getHours() < c.getHours(); }
    // <
    bool operator<(const double& in_hour) const
    { return getHours() < in_hour; }
    //.............
    // >=
    bool operator>=(const Clock& c) const
    { return getHours() >= c.getHours(); }
    // >=
    bool operator>=(const double& in_hour) const
    { return getHours() >= in_hour; }
    //.............
    // <=
    bool operator<=(const Clock& c) const
    { return getHours() <= c.getHours(); }
    // <=
    bool operator<=(const double& in_hour) const
    { return getHours() <= in_hour; }
    //..................................................
};

//==============================================================
// Gregorian Date
struct Date
{
    // Negative year would be BCE
    int yr      ;
    int mo    ;
    int day   ;

    Clock cl ;
    
    // Constructor 0
    Date() {}

    // Constructor 1
    Date(const int& yr, const int& mo, const int& day, const Clock& cl) ;

    // Constructor 2
    // Sets clock to 0:0:0
    Date(const int& yr, const int& mo, const int& day) ;

    bool operator>(const Date& other) const
    { return JD() > other.JD(); }

    bool operator>=(const Date& other) const
    { return JD() >= other.JD(); }

    bool operator<(const Date& other) const
    { return JD() < other.JD(); }

    bool operator<=(const Date& other) const
    { return JD() <= other.JD(); }

    bool operator==(const Date& other) const
    { return JD() == other.JD(); }

    double operator+(const Date& other) const
    { return JD() + other.JD(); }

    double operator+(const double& other) const
    { return JD() + other; }

    double operator-(const Date& other) const
    { return JD() - other.JD(); }

    //...........................................
    // Julian dates
    // Ref: https://en.wikipedia.org/wiki/Julian_day
    //...........................................
    // The Julian Day number is the number of days elapsed since
    //  noon (12:00) on January 1, 4713 BC.
    double JD() const ;

    // Reduced JD, Epoch: 12h Nov 16, 1858
    // RedJD = JD − 2400000
    double RedJD() const ;

    // Modified JD, Epoch: 0h Nov 17, 1858
    // Introduced by SAO in 1957
    // ModJD = JD − 2400000.5
    double ModJD() const ;

    // Truncated JD, Epoch: 0h May 24, 1968
    // Introduced by NASA in 1979
    // TrunJD = floor(JD − 2440000.5)
    double TrunJD() const ;

    // UnixTDay is similar to Unix time but it's in days
    // Epoch: 0h Jan 1, 1970
    // UnixTDay = JD − 2440587.5
    double UnixTDay() const ;

    // Unix Time, Epoch: 0h Jan 1, 1970
    // Count of seconds (excluding leap seconds)
    // Astronomical almanac for the year 2001, 2000, p. K2
    // UnixT = (JD − 2440587.5) × 86400
    double UnixT() const ;

    // RootT95 is similar to Unix time
    // but it starts from 1995 which is 
    // for plotting purposes in root.
    // Epoch: 0h Jan 1, 1995
    // RootT95 = JD − {1995, 1, 1, {0,0,0}}.JD()
    double RootT95() const ;
    //...........................................

    std::string StrForm(const std::string& format) const ;

};
//==============================================================
struct DateInterval
{
    Date start, end;

    double duration()
    {
      return end - start ;
    }
};

//==============================================================
//                  ostream << overloading
std::ostream& operator << (std::ostream &output, const Zaki::Physics::Date&) ;
std::ostream& operator << (std::ostream &output, const Zaki::Physics::DateInterval&) ;
std::ostream& operator << (std::ostream &output, const Zaki::Physics::Clock&) ;
//==============================================================

//--------------------------------------------------------------
} // End of namespace Zaki::Physics
//--------------------------------------------------------------


#endif /*Zaki_Physics_DateTime_H*/
