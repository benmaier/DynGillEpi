/* 
 * The MIT License (MIT)
 * Copyright (c) 2016, Benjamin Maier
 *
 * Permission is hereby granted, free of charge, to any person 
 * obtaining a copy of this software and associated documentation 
 * files (the "Software"), to deal in the Software without 
 * restriction, including without limitation the rights to use, 
 * copy, modify, merge, publish, distribute, sublicense, and/or 
 * sell copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall 
 * be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-
 * INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
 * AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF 
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
 * IN THE SOFTWARE.
 */
#ifndef __UTILS__
#define __UTILS__

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <iterator>
#include <set>
#include <utility>
#include <cmath>
#include <numeric>
#include <random>
#include <ctime>
#include <cstdlib>
#include <tuple>

using namespace std;

struct SI_Result {
    vector < vector < size_t > > I;
    vector < vector < size_t > > SI;
    vector < size_t > hist;
};

//======================================================================
// Typedef
//======================================================================
// Variable types:
typedef unsigned int COUNTER;
typedef unsigned int NODE;
typedef vector<NODE> NODES; // list of nodes
typedef vector<bool> BOOLS;
typedef CONTACT pair<NODE>; // contact (i,j)
typedef vector<CONTACT> CONTACTS; // contacts in a single time-frame
typedef vector<CONTACTS> CONTACTS_LIST; // list of contact lists
// Random number generators:
typedef mt19937 ENG; // use Mersenne Twister 19937 as PRNG engine
typedef uniform_int_distribution<size_t> DIST_INT; // define uniform distribution of integers
typedef uniform_real_distribution<double> DIST_REAL; // define uniform distribution of reals on [0,1)
typedef exponential_distribution<double> DIST_EXP; // define exponential distribution


vector<size_t>::iterator choose_random_unique(
        vector<size_t>::iterator begin, 
        vector<size_t>::iterator end, 
        size_t num_random,
        default_random_engine & generator,
        uniform_real_distribution<double> & distribution
    );

#endif
