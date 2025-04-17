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
 * @file Vector_Basic.hpp
 *
 * @brief Basic vector operations.
 *
 * @ingroup Vector
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_Vector_Basic_H
#define Zaki_Vector_Basic_H

// Functions for vector manipulations

#include <algorithm>
#include <vector>

//==============================================================
namespace Zaki::Vector
{

//==============================================================
// Checks if an element exists in the list.
template <class T> 
bool Exists(const T& Element, const std::vector<T>& Vec) 
{
    if (std::find(Vec.begin(), Vec.end(), Element) != Vec.end())
        return true ;

    return false ;
}

//==============================================================
/// Checks if an element exists in the list, if yes
/// it will return its index, if not it will return '-1'
template <class T, class U> 
int GetIdx(const T& Element, const std::vector<U>& Vec)
{
    auto it = std::find(Vec.begin(), Vec.end(), Element) ;
 
    // If element was found
    if (it != Vec.end())
    {
        return it - Vec.begin() ;
    }
    else 
    {
        // If the element is not present in the vector
        return -1 ;
    }
}

//==============================================================
// Adds element to a list if it doesn't already exists.
template <class T>
void Add(std::vector<T>& lst, const T &ele)
// Add an element to a vector, 
// only if it doesn't already exist.
{
    if(!Exists(ele ,lst))  lst.push_back(ele) ;
}

//==============================================================
// Removes a list of elements from a T type list
template <class T> 
void Remove(std::vector<T>& list, const std::vector<int>& rm_index_list)
// Removing certain members of a vector, 
// given another vector containing an index set.
{
        int    j = 0 ;
    for(size_t i = 0 ; i < list.size() ; ++i)
    {
        if( Exists((int) i+j, rm_index_list))
        {
            list.erase(list.begin()+i) ;
            i--                        ;
            j++                        ;
        }
    }
}

//==============================================================

//--------------------------------------------------------------
} // End of namespace Zaki::Vector
//==============================================================

#endif /*Zaki_Vector_Basic_H*/
