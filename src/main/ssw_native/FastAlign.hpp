/* The MIT License
   Copyright (c) 2012-2015 Boston College.
   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:
   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* The 2-clause BSD License
   Copyright 2006 Michael Farrar.
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 *  ssw.c
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 1.2.4
 *	Last revision by Mengyao Zhao on 2022-Apr-17.
 *
 *  The lazy-F loop implementation was derived from SWPS3, which is
 *  MIT licensed under ETH Zürich, Institute of Computational Science.
 *
 *  The core SW loop referenced the swsse2 implementation, which is
 *  BSD licensed under Micharl Farrar.
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "ssw.h"
#include "avx2_functions.h"


#include <iostream>

#include <string>

#include <stack>



#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

using namespace std;

class FastAlign{

public:

        FastAlign(int gap_open, int gap_extend, char scoring_matrix, ){
        }
        void initialize_NUCC_matrix() {

            match = new int[256][256];
            match['A']['A'] = 5;
            match['A']['T'] = -4;
            match['A']['G'] = -4;
            match['A']['C'] = -4;
            match['A']['S'] = -4;
            match['A']['W'] = 1;
            match['A']['R'] = 1;
            match['A']['Y'] = -4;
            match['A']['K'] = -4;
            match['A']['M'] = 1;
            match['A']['B'] = -4;
            match['A']['V'] = -1;
            match['A']['H'] = -1;
            match['A']['D'] = -1;
            match['A']['N'] = -2;
            match['T']['A'] = -4;
            match['T']['T'] = 5;
            match['T']['G'] = -4;
            match['T']['C'] = -4;
            match['T']['S'] = -4;
            match['T']['W'] = 1;
            match['T']['R'] = -4;
            match['T']['Y'] = 1;
            match['T']['K'] = 1;
            match['T']['M'] = -4;
            match['T']['B'] = -1;
            match['T']['V'] = -4;
            match['T']['H'] = -1;
            match['T']['D'] = -1;
            match['T']['N'] = -2;
            match['G']['A'] = -4;
            match['G']['T'] = -4;
            match['G']['G'] = 5;
            match['G']['C'] = -4;
            match['G']['S'] = 1;
            match['G']['W'] = -4;
            match['G']['R'] = 1;
            match['G']['Y'] = -4;
            match['G']['K'] = 1;
            match['G']['M'] = -4;
            match['G']['B'] = -1;
            match['G']['V'] = -1;
            match['G']['H'] = -4;
            match['G']['D'] = -1;
            match['G']['N'] = -2;
            match['C']['A'] = -4;
            match['C']['T'] = -4;
            match['C']['G'] = -4;
            match['C']['C'] = 5;
            match['C']['S'] = 1;
            match['C']['W'] = -4;
            match['C']['R'] = -4;
            match['C']['Y'] = 1;
            match['C']['K'] = -4;
            match['C']['M'] = 1;
            match['C']['B'] = -1;
            match['C']['V'] = -1;
            match['C']['H'] = -1;
            match['C']['D'] = -4;
            match['C']['N'] = -2;
            match['S']['A'] = -4;
            match['S']['T'] = -4;
            match['S']['G'] = 1;
            match['S']['C'] = 1;
            match['S']['S'] = -1;
            match['S']['W'] = -4;
            match['S']['R'] = -2;
            match['S']['Y'] = -2;
            match['S']['K'] = -2;
            match['S']['M'] = -2;
            match['S']['B'] = -1;
            match['S']['V'] = -1;
            match['S']['H'] = -3;
            match['S']['D'] = -3;
            match['S']['N'] = -1;
            match['W']['A'] = 1;
            match['W']['T'] = 1;
            match['W']['G'] = -4;
            match['W']['C'] = -4;
            match['W']['S'] = -4;
            match['W']['W'] = -1;
            match['W']['R'] = -2;
            match['W']['Y'] = -2;
            match['W']['K'] = -2;
            match['W']['M'] = -2;
            match['W']['B'] = -3;
            match['W']['V'] = -3;
            match['W']['H'] = -1;
            match['W']['D'] = -1;
            match['W']['N'] = -1;
            match['R']['A'] = 1;
            match['R']['T'] = -4;
            match['R']['G'] = 1;
            match['R']['C'] = -4;
            match['R']['S'] = -2;
            match['R']['W'] = -2;
            match['R']['R'] = -1;
            match['R']['Y'] = -4;
            match['R']['K'] = -2;
            match['R']['M'] = -2;
            match['R']['B'] = -3;
            match['R']['V'] = -1;
            match['R']['H'] = -3;
            match['R']['D'] = -1;
            match['R']['N'] = -1;
            match['Y']['A'] = -4;
            match['Y']['T'] = 1;
            match['Y']['G'] = -4;
            match['Y']['C'] = 1;
            match['Y']['S'] = -2;
            match['Y']['W'] = -2;
            match['Y']['R'] = -4;
            match['Y']['Y'] = -1;
            match['Y']['K'] = -2;
            match['Y']['M'] = -2;
            match['Y']['B'] = -1;
            match['Y']['V'] = -3;
            match['Y']['H'] = -1;
            match['Y']['D'] = -3;
            match['Y']['N'] = -1;
            match['K']['A'] = -4;
            match['K']['T'] = 1;
            match['K']['G'] = 1;
            match['K']['C'] = -4;
            match['K']['S'] = -2;
            match['K']['W'] = -2;
            match['K']['R'] = -2;
            match['K']['Y'] = -2;
            match['K']['K'] = -1;
            match['K']['M'] = -4;
            match['K']['B'] = -1;
            match['K']['V'] = -3;
            match['K']['H'] = -3;
            match['K']['D'] = -1;
            match['K']['N'] = -1;
            match['M']['A'] = 1;
            match['M']['T'] = -4;
            match['M']['G'] = -4;
            match['M']['C'] = 1;
            match['M']['S'] = -2;
            match['M']['W'] = -2;
            match['M']['R'] = -2;
            match['M']['Y'] = -2;
            match['M']['K'] = -4;
            match['M']['M'] = -1;
            match['M']['B'] = -3;
            match['M']['V'] = -1;
            match['M']['H'] = -1;
            match['M']['D'] = -3;
            match['M']['N'] = -1;
            match['B']['A'] = -4;
            match['B']['T'] = -1;
            match['B']['G'] = -1;
            match['B']['C'] = -1;
            match['B']['S'] = -1;
            match['B']['W'] = -3;
            match['B']['R'] = -3;
            match['B']['Y'] = -1;
            match['B']['K'] = -1;
            match['B']['M'] = -3;
            match['B']['B'] = -1;
            match['B']['V'] = -2;
            match['B']['H'] = -2;
            match['B']['D'] = -2;
            match['B']['N'] = -1;
            match['V']['A'] = -1;
            match['V']['T'] = -4;
            match['V']['G'] = -1;
            match['V']['C'] = -1;
            match['V']['S'] = -1;
            match['V']['W'] = -3;
            match['V']['R'] = -1;
            match['V']['Y'] = -3;
            match['V']['K'] = -3;
            match['V']['M'] = -1;
            match['V']['B'] = -2;
            match['V']['V'] = -1;
            match['V']['H'] = -2;
            match['V']['D'] = -2;
            match['V']['N'] = -1;
            match['H']['A'] = -1;
            match['H']['T'] = -1;
            match['H']['G'] = -4;
            match['H']['C'] = -1;
            match['H']['S'] = -3;
            match['H']['W'] = -1;
            match['H']['R'] = -3;
            match['H']['Y'] = -1;
            match['H']['K'] = -3;
            match['H']['M'] = -1;
            match['H']['B'] = -2;
            match['H']['V'] = -2;
            match['H']['H'] = -1;
            match['H']['D'] = -2;
            match['H']['N'] = -1;
            match['D']['A'] = -1;
            match['D']['T'] = -1;
            match['D']['G'] = -1;
            match['D']['C'] = -4;
            match['D']['S'] = -3;
            match['D']['W'] = -1;
            match['D']['R'] = -1;
            match['D']['Y'] = -3;
            match['D']['K'] = -1;
            match['D']['M'] = -3;
            match['D']['B'] = -2;
            match['D']['V'] = -2;
            match['D']['H'] = -2;
            match['D']['D'] = -1;
            match['D']['N'] = -1;
            match['N']['A'] = -2;
            match['N']['T'] = -2;
            match['N']['G'] = -2;
            match['N']['C'] = -2;
            match['N']['S'] = -1;
            match['N']['W'] = -1;
            match['N']['R'] = -1;
            match['N']['Y'] = -1;
            match['N']['K'] = -1;
            match['N']['M'] = -1;
            match['N']['B'] = -1;
            match['N']['V'] = -1;
            match['N']['H'] = -1;
            match['N']['D'] = -1;
            match['N']['N'] = -1;

    }

    /**
     * Initializes the BLOSUM62 scoring matrix.
     */
     void initialize_BLOSUM_matrix() {


        match['A']['A'] = 4;
        match['A']['R'] = -1;
        match['A']['N'] = -2;
        match['A']['D'] = -2;
        match['A']['C'] = 0;
        match['A']['Q'] = -1;
        match['A']['E'] = -1;
        match['A']['G'] = 0;
        match['A']['H'] = -2;
        match['A']['I'] = -1;
        match['A']['L'] = -1;
        match['A']['K'] = -1;
        match['A']['M'] = -1;
        match['A']['F'] = -2;
        match['A']['P'] = -1;
        match['A']['S'] = 1;
        match['A']['T'] = 0;
        match['A']['W'] = -3;
        match['A']['Y'] = -2;
        match['A']['V'] = 0;
        match['A']['B'] = -2;
        match['A']['Z'] = -1;
        match['A']['X'] = 0;
        match['A']['*'] = -4;


        match['R']['A'] = -1;
        match['R']['R'] = 5;
        match['R']['N'] = 0;
        match['R']['D'] = -2;
        match['R']['C'] = -3;
        match['R']['Q'] = 1;
        match['R']['E'] = 0;
        match['R']['G'] = -2;
        match['R']['H'] = 0;
        match['R']['I'] = -3;
        match['R']['L'] = -2;
        match['R']['K'] = 2;
        match['R']['M'] = -1;
        match['R']['F'] = -3;
        match['R']['P'] = -2;
        match['R']['S'] = -1;
        match['R']['T'] = -1;
        match['R']['W'] = -3;
        match['R']['Y'] = -2;
        match['R']['V'] = -3;
        match['R']['B'] = -1;
        match['R']['Z'] = 0;
        match['R']['X'] = -1;
        match['R']['*'] = -4;


        match['N']['A'] = -2;
        match['N']['R'] = 0;
        match['N']['N'] = 6;
        match['N']['D'] = 1;
        match['N']['C'] = -3;
        match['N']['Q'] = 0;
        match['N']['E'] = 0;
        match['N']['G'] = 0;
        match['N']['H'] = 1;
        match['N']['I'] = -3;
        match['N']['L'] = -3;
        match['N']['K'] = 0;
        match['N']['M'] = -2;
        match['N']['F'] = -3;
        match['N']['P'] = -2;
        match['N']['S'] = 1;
        match['N']['T'] = 0;
        match['N']['W'] = -4;
        match['N']['Y'] = -2;
        match['N']['V'] = -3;
        match['N']['B'] = 3;
        match['N']['Z'] = 0;
        match['N']['X'] = -1;
        match['N']['*'] = -4;


        match['D']['A'] = -2;
        match['D']['R'] = -2;
        match['D']['N'] = 1;
        match['D']['D'] = 6;
        match['D']['C'] = -3;
        match['D']['Q'] = 0;
        match['D']['E'] = 2;
        match['D']['G'] = -1;
        match['D']['H'] = -1;
        match['D']['I'] = -3;
        match['D']['L'] = -4;
        match['D']['K'] = -1;
        match['D']['M'] = -3;
        match['D']['F'] = -3;
        match['D']['P'] = -1;
        match['D']['S'] = 0;
        match['D']['T'] = -1;
        match['D']['W'] = -4;
        match['D']['Y'] = -3;
        match['D']['V'] = -3;
        match['D']['B'] = 4;
        match['D']['Z'] = 1;
        match['D']['X'] = -1;
        match['D']['*'] = -4;


        match['C']['A'] = 0;
        match['C']['R'] = -3;
        match['C']['N'] = -3;
        match['C']['D'] = -3;
        match['C']['C'] = 9;
        match['C']['Q'] = -3;
        match['C']['E'] = -4;
        match['C']['G'] = -3;
        match['C']['H'] = -3;
        match['C']['I'] = -1;
        match['C']['L'] = -1;
        match['C']['K'] = -3;
        match['C']['M'] = -1;
        match['C']['F'] = -2;
        match['C']['P'] = -3;
        match['C']['S'] = -1;
        match['C']['T'] = -1;
        match['C']['W'] = -2;
        match['C']['Y'] = -2;
        match['C']['V'] = -1;
        match['C']['B'] = -3;
        match['C']['Z'] = -3;
        match['C']['X'] = -2;
        match['C']['*'] = -4;


        match['Q']['A'] = -1;
        match['Q']['R'] = 1;
        match['Q']['N'] = 0;
        match['Q']['D'] = 0;
        match['Q']['C'] = -3;
        match['Q']['Q'] = 5;
        match['Q']['E'] = 2;
        match['Q']['G'] = -2;
        match['Q']['H'] = 0;
        match['Q']['I'] = -3;
        match['Q']['L'] = -2;
        match['Q']['K'] = 1;
        match['Q']['M'] = 0;
        match['Q']['F'] = -3;
        match['Q']['P'] = -1;
        match['Q']['S'] = 0;
        match['Q']['T'] = -1;
        match['Q']['W'] = -2;
        match['Q']['Y'] = -1;
        match['Q']['V'] = -2;
        match['Q']['B'] = 0;
        match['Q']['Z'] = 3;
        match['Q']['X'] = -1;
        match['Q']['*'] = -4;


        match['E']['A'] = -1;
        match['E']['R'] = 0;
        match['E']['N'] = 0;
        match['E']['D'] = 2;
        match['E']['C'] = -4;
        match['E']['Q'] = 2;
        match['E']['E'] = 5;
        match['E']['G'] = -2;
        match['E']['H'] = 0;
        match['E']['I'] = -3;
        match['E']['L'] = -3;
        match['E']['K'] = 1;
        match['E']['M'] = -2;
        match['E']['F'] = -3;
        match['E']['P'] = -1;
        match['E']['S'] = 0;
        match['E']['T'] = -1;
        match['E']['W'] = -3;
        match['E']['Y'] = -2;
        match['E']['V'] = -2;
        match['E']['B'] = 1;
        match['E']['Z'] = 4;
        match['E']['X'] = -1;
        match['E']['*'] = -4;


        match['G']['A'] = 0;
        match['G']['R'] = -2;
        match['G']['N'] = 0;
        match['G']['D'] = -1;
        match['G']['C'] = -3;
        match['G']['Q'] = -2;
        match['G']['E'] = -2;
        match['G']['G'] = 6;
        match['G']['H'] = -2;
        match['G']['I'] = -4;
        match['G']['L'] = -4;
        match['G']['K'] = -2;
        match['G']['M'] = -3;
        match['G']['F'] = -3;
        match['G']['P'] = -2;
        match['G']['S'] = 0;
        match['G']['T'] = -2;
        match['G']['W'] = -2;
        match['G']['Y'] = -3;
        match['G']['V'] = -3;
        match['G']['B'] = -1;
        match['G']['Z'] = -2;
        match['G']['X'] = -1;
        match['G']['*'] = -4;


        match['H']['A'] = -2;
        match['H']['R'] = 0;
        match['H']['N'] = 1;
        match['H']['D'] = -1;
        match['H']['C'] = -3;
        match['H']['Q'] = 0;
        match['H']['E'] = 0;
        match['H']['G'] = -2;
        match['H']['H'] = 8;
        match['H']['I'] = -3;
        match['H']['L'] = -3;
        match['H']['K'] = -1;
        match['H']['M'] = -2;
        match['H']['F'] = -1;
        match['H']['P'] = -2;
        match['H']['S'] = -1;
        match['H']['T'] = -2;
        match['H']['W'] = -2;
        match['H']['Y'] = 2;
        match['H']['V'] = -3;
        match['H']['B'] = 0;
        match['H']['Z'] = 0;
        match['H']['X'] = -1;
        match['H']['*'] = -4;


        match['I']['A'] = -1;
        match['I']['R'] = -3;
        match['I']['N'] = -3;
        match['I']['D'] = -3;
        match['I']['C'] = -1;
        match['I']['Q'] = -3;
        match['I']['E'] = -3;
        match['I']['G'] = -4;
        match['I']['H'] = -3;
        match['I']['I'] = 4;
        match['I']['L'] = 2;
        match['I']['K'] = -3;
        match['I']['M'] = 1;
        match['I']['F'] = 0;
        match['I']['P'] = -3;
        match['I']['S'] = -2;
        match['I']['T'] = -1;
        match['I']['W'] = -3;
        match['I']['Y'] = -1;
        match['I']['V'] = 3;
        match['I']['B'] = -3;
        match['I']['Z'] = -3;
        match['I']['X'] = -1;
        match['I']['*'] = -4;


        match['L']['A'] = -1;
        match['L']['R'] = -2;
        match['L']['N'] = -3;
        match['L']['D'] = -4;
        match['L']['C'] = -1;
        match['L']['Q'] = -2;
        match['L']['E'] = -3;
        match['L']['G'] = -4;
        match['L']['H'] = -3;
        match['L']['I'] = 2;
        match['L']['L'] = 4;
        match['L']['K'] = -2;
        match['L']['M'] = 2;
        match['L']['F'] = 0;
        match['L']['P'] = -3;
        match['L']['S'] = -2;
        match['L']['T'] = -1;
        match['L']['W'] = -2;
        match['L']['Y'] = -1;
        match['L']['V'] = 1;
        match['L']['B'] = -4;
        match['L']['Z'] = -3;
        match['L']['X'] = -1;
        match['L']['*'] = -4;


        match['K']['A'] = -1;
        match['K']['R'] = 2;
        match['K']['N'] = 0;
        match['K']['D'] = -1;
        match['K']['C'] = -3;
        match['K']['Q'] = 1;
        match['K']['E'] = 1;
        match['K']['G'] = -2;
        match['K']['H'] = -1;
        match['K']['I'] = -3;
        match['K']['L'] = -2;
        match['K']['K'] = 5;
        match['K']['M'] = -1;
        match['K']['F'] = -3;
        match['K']['P'] = -1;
        match['K']['S'] = 0;
        match['K']['T'] = -1;
        match['K']['W'] = -3;
        match['K']['Y'] = -2;
        match['K']['V'] = -2;
        match['K']['B'] = 0;
        match['K']['Z'] = 1;
        match['K']['X'] = -1;
        match['K']['*'] = -4;


        match['M']['A'] = -1;
        match['M']['R'] = -1;
        match['M']['N'] = -2;
        match['M']['D'] = -3;
        match['M']['C'] = -1;
        match['M']['Q'] = 0;
        match['M']['E'] = -2;
        match['M']['G'] = -3;
        match['M']['H'] = -2;
        match['M']['I'] = 1;
        match['M']['L'] = 2;
        match['M']['K'] = -1;
        match['M']['M'] = 5;
        match['M']['F'] = 0;
        match['M']['P'] = -2;
        match['M']['S'] = -1;
        match['M']['T'] = -1;
        match['M']['W'] = -1;
        match['M']['Y'] = -1;
        match['M']['V'] = 1;
        match['M']['B'] = -3;
        match['M']['Z'] = -1;
        match['M']['X'] = -1;
        match['M']['*'] = -4;


        match['F']['A'] = -2;
        match['F']['R'] = -3;
        match['F']['N'] = -3;
        match['F']['D'] = -3;
        match['F']['C'] = -2;
        match['F']['Q'] = -3;
        match['F']['E'] = -3;
        match['F']['G'] = -3;
        match['F']['H'] = -1;
        match['F']['I'] = 0;
        match['F']['L'] = 0;
        match['F']['K'] = -3;
        match['F']['M'] = 0;
        match['F']['F'] = 6;
        match['F']['P'] = -4;
        match['F']['S'] = -2;
        match['F']['T'] = -2;
        match['F']['W'] = 1;
        match['F']['Y'] = 3;
        match['F']['V'] = -1;
        match['F']['B'] = -3;
        match['F']['Z'] = -3;
        match['F']['X'] = -1;
        match['F']['*'] = -4;


        match['P']['A'] = -1;
        match['P']['R'] = -2;
        match['P']['N'] = -2;
        match['P']['D'] = -1;
        match['P']['C'] = -3;
        match['P']['Q'] = -1;
        match['P']['E'] = -1;
        match['P']['G'] = -2;
        match['P']['H'] = -2;
        match['P']['I'] = -3;
        match['P']['L'] = -3;
        match['P']['K'] = -1;
        match['P']['M'] = -2;
        match['P']['F'] = -4;
        match['P']['P'] = 7;
        match['P']['S'] = -1;
        match['P']['T'] = -1;
        match['P']['W'] = -4;
        match['P']['Y'] = -3;
        match['P']['V'] = -2;
        match['P']['B'] = -2;
        match['P']['Z'] = -1;
        match['P']['X'] = -2;
        match['P']['*'] = -4;


        match['S']['A'] = 1;
        match['S']['R'] = -1;
        match['S']['N'] = 1;
        match['S']['D'] = 0;
        match['S']['C'] = -1;
        match['S']['Q'] = 0;
        match['S']['E'] = 0;
        match['S']['G'] = 0;
        match['S']['H'] = -1;
        match['S']['I'] = -2;
        match['S']['L'] = -2;
        match['S']['K'] = 0;
        match['S']['M'] = -1;
        match['S']['F'] = -2;
        match['S']['P'] = -1;
        match['S']['S'] = 4;
        match['S']['T'] = 1;
        match['S']['W'] = -3;
        match['S']['Y'] = -2;
        match['S']['V'] = -2;
        match['S']['B'] = 0;
        match['S']['Z'] = 0;
        match['S']['X'] = 0;
        match['S']['*'] = -4;


        match['T']['A'] = 0;
        match['T']['R'] = -1;
        match['T']['N'] = 0;
        match['T']['D'] = -1;
        match['T']['C'] = -1;
        match['T']['Q'] = -1;
        match['T']['E'] = -1;
        match['T']['G'] = -2;
        match['T']['H'] = -2;
        match['T']['I'] = -1;
        match['T']['L'] = -1;
        match['T']['K'] = -1;
        match['T']['M'] = -1;
        match['T']['F'] = -2;
        match['T']['P'] = -1;
        match['T']['S'] = 1;
        match['T']['T'] = 5;
        match['T']['W'] = -2;
        match['T']['Y'] = -2;
        match['T']['V'] = 0;
        match['T']['B'] = -1;
        match['T']['Z'] = -1;
        match['T']['X'] = 0;
        match['T']['*'] = -4;


        match['W']['A'] = -3;
        match['W']['R'] = -3;
        match['W']['N'] = -4;
        match['W']['D'] = -4;
        match['W']['C'] = -2;
        match['W']['Q'] = -2;
        match['W']['E'] = -3;
        match['W']['G'] = -2;
        match['W']['H'] = -2;
        match['W']['I'] = -3;
        match['W']['L'] = -2;
        match['W']['K'] = -3;
        match['W']['M'] = -1;
        match['W']['F'] = 1;
        match['W']['P'] = -4;
        match['W']['S'] = -3;
        match['W']['T'] = -2;
        match['W']['W'] = 11;
        match['W']['Y'] = 2;
        match['W']['V'] = -3;
        match['W']['B'] = -4;
        match['W']['Z'] = -3;
        match['W']['X'] = -2;
        match['W']['*'] = -4;


        match['Y']['A'] = -2;
        match['Y']['R'] = -2;
        match['Y']['N'] = -2;
        match['Y']['D'] = -3;
        match['Y']['C'] = -2;
        match['Y']['Q'] = -1;
        match['Y']['E'] = -2;
        match['Y']['G'] = -3;
        match['Y']['H'] = 2;
        match['Y']['I'] = -1;
        match['Y']['L'] = -1;
        match['Y']['K'] = -2;
        match['Y']['M'] = -1;
        match['Y']['F'] = 3;
        match['Y']['P'] = -3;
        match['Y']['S'] = -2;
        match['Y']['T'] = -2;
        match['Y']['W'] = 2;
        match['Y']['Y'] = 7;
        match['Y']['V'] = -1;
        match['Y']['B'] = -3;
        match['Y']['Z'] = -2;
        match['Y']['X'] = -1;
        match['Y']['*'] = -4;


        match['V']['A'] = 0;
        match['V']['R'] = -3;
        match['V']['N'] = -3;
        match['V']['D'] = -3;
        match['V']['C'] = -1;
        match['V']['Q'] = -2;
        match['V']['E'] = -2;
        match['V']['G'] = -3;
        match['V']['H'] = -3;
        match['V']['I'] = 3;
        match['V']['L'] = 1;
        match['V']['K'] = -2;
        match['V']['M'] = 1;
        match['V']['F'] = -1;
        match['V']['P'] = -2;
        match['V']['S'] = -2;
        match['V']['T'] = 0;
        match['V']['W'] = -3;
        match['V']['Y'] = -1;
        match['V']['V'] = 4;
        match['V']['B'] = -3;
        match['V']['Z'] = -2;
        match['V']['X'] = -1;
        match['V']['*'] = -4;


        match['B']['A'] = -2;
        match['B']['R'] = -1;
        match['B']['N'] = 3;
        match['B']['D'] = 4;
        match['B']['C'] = -3;
        match['B']['Q'] = 0;
        match['B']['E'] = 1;
        match['B']['G'] = -1;
        match['B']['H'] = 0;
        match['B']['I'] = -3;
        match['B']['L'] = -4;
        match['B']['K'] = 0;
        match['B']['M'] = -3;
        match['B']['F'] = -3;
        match['B']['P'] = -2;
        match['B']['S'] = 0;
        match['B']['T'] = -1;
        match['B']['W'] = -4;
        match['B']['Y'] = -3;
        match['B']['V'] = -3;
        match['B']['B'] = 4;
        match['B']['Z'] = 1;
        match['B']['X'] = -1;
        match['B']['*'] = -4;


        match['Z']['A'] = -1;
        match['Z']['R'] = 0;
        match['Z']['N'] = 0;
        match['Z']['D'] = 1;
        match['Z']['C'] = -3;
        match['Z']['Q'] = 3;
        match['Z']['E'] = 4;
        match['Z']['G'] = -2;
        match['Z']['H'] = 0;
        match['Z']['I'] = -3;
        match['Z']['L'] = -3;
        match['Z']['K'] = 1;
        match['Z']['M'] = -1;
        match['Z']['F'] = -3;
        match['Z']['P'] = -1;
        match['Z']['S'] = 0;
        match['Z']['T'] = -1;
        match['Z']['W'] = -3;
        match['Z']['Y'] = -2;
        match['Z']['V'] = -2;
        match['Z']['B'] = 1;
        match['Z']['Z'] = 4;
        match['Z']['X'] = -1;
        match['Z']['*'] = -4;


        match['X']['A'] = 0;
        match['X']['R'] = -1;
        match['X']['N'] = -1;
        match['X']['D'] = -1;
        match['X']['C'] = -2;
        match['X']['Q'] = -1;
        match['X']['E'] = -1;
        match['X']['G'] = -1;
        match['X']['H'] = -1;
        match['X']['I'] = -1;
        match['X']['L'] = -1;
        match['X']['K'] = -1;
        match['X']['M'] = -1;
        match['X']['F'] = -1;
        match['X']['P'] = -2;
        match['X']['S'] = 0;
        match['X']['T'] = 0;
        match['X']['W'] = -2;
        match['X']['Y'] = -1;
        match['X']['V'] = -1;
        match['X']['B'] = -1;
        match['X']['Z'] = -1;
        match['X']['X'] = -1;
        match['X']['*'] = -4;


        match['*']['A'] = -4;
        match['*']['R'] = -4;
        match['*']['N'] = -4;
        match['*']['D'] = -4;
        match['*']['C'] = -4;
        match['*']['Q'] = -4;
        match['*']['E'] = -4;
        match['*']['G'] = -4;
        match['*']['H'] = -4;
        match['*']['I'] = -4;
        match['*']['L'] = -4;
        match['*']['K'] = -4;
        match['*']['M'] = -4;
        match['*']['F'] = -4;
        match['*']['P'] = -4;
        match['*']['S'] = -4;
        match['*']['T'] = -4;
        match['*']['W'] = -4;
        match['*']['Y'] = -4;
        match['*']['V'] = -4;
        match['*']['B'] = -4;
        match['*']['Z'] = -4;
        match['*']['X'] = -4;
        match['*']['*'] = 1;
    }




        static VEC_TYPE* qP_word (const int8_t* read_num,
        const int8_t* mat,
        const int32_t readLen,
        const int32_t n) {

            int32_t segLen = (readLen + ((sizeof(VEC_TYPE) / sizeof(int16_t)) - 1)) /
                             (sizeof(VEC_TYPE) / sizeof(int16_t)); /* number of segment */
            VEC_TYPE *vProfile = (VEC_TYPE *) aligned_alloc(sizeof(VEC_TYPE), n * segLen * sizeof(VEC_TYPE));
            int16_t *t = (int16_t *) vProfile;
            int32_t nt, i, j;
            int32_t segNum;

            /* Generate query profile rearrange query sequence & calculate the weight of match/mismatch */
            for (nt = 0; LIKELY(nt < n); nt++) {
                for (i = 0; i < segLen; i++) {
                    j = i;
                    for (segNum = 0; LIKELY(segNum < (sizeof(VEC_TYPE) / sizeof(int16_t))); segNum++) {
                        *t++ = j >= readLen ? 0 : mat[nt * n + read_num[j]];
                        j += segLen;
                    }
                }
            }
            return vProfile;
        }

        static alignment_end* ssw_int16 (const int8_t* ref,
        int32_t refLen,
        int32_t readLen,
        const uint8_t weight_gapO, /* will be used as - */
        const uint8_t weight_gapE, /* will be used as - */
        const VEC_TYPE* vProfile) {


            int16_t max = SHRT_MIN;                             /* the max alignment score */
            int32_t end_read = readLen - 1;
            int32_t end_ref = 0; /* 1_based best alignment ending point; Initialized as isn't aligned - 0. */
            int32_t segLen = (readLen + ((sizeof(VEC_TYPE) / sizeof(int16_t)) - 1)) /
                             (sizeof(VEC_TYPE) / sizeof(int16_t)); /* number of segment */

            alignment_end *best = (alignment_end *) calloc(1, sizeof(alignment_end));

            int32_t segLen_per_btvec = segLen / 8;

            int32_t seg_left = segLen - (8 * segLen_per_btvec);




            /* Define 16 byte 0 vector. */
            VEC_TYPE vNeg = VEC_SET1(SHRT_MIN);
            VEC_TYPE vZero = VEC_SET1(0);

            VEC_TYPE *pvHStore = (VEC_TYPE *) aligned_alloc(sizeof(VEC_TYPE), segLen * sizeof(VEC_TYPE));
            //memset(pvHStore, 0, segLen*sizeof(VEC_TYPE));
            VEC_TYPE *pvHLoad = (VEC_TYPE *) aligned_alloc(sizeof(VEC_TYPE), segLen * sizeof(VEC_TYPE));
            //memset(pvHLoad, 0, segLen*sizeof(VEC_TYPE));
            VEC_TYPE *pvE = (VEC_TYPE *) aligned_alloc(sizeof(VEC_TYPE), segLen * sizeof(VEC_TYPE));
            //int itr;
            /*for (itr = 0; itr < segLen; itr++){
                pvE[itr] = vNeg;
            }*/
            VEC_TYPE *pvHmax = (VEC_TYPE *) aligned_alloc(sizeof(VEC_TYPE), segLen * sizeof(VEC_TYPE));

            int32_t i, j, k, l;
            /* 16 byte insertion begin vector */
            VEC_TYPE vGapOE = VEC_SET1(weight_gapO + weight_gapE);

            /* 16 byte insertion extension vector */
            VEC_TYPE vGapE = VEC_SET1(weight_gapE);

            VEC_TYPE vMaxScore = vNeg; /* Trace the highest score of the whole SW matrix. */
            VEC_TYPE vMaxMark = vNeg; /* Trace the highest score till the previous column. */

            VEC_TYPE vTemp;
            int32_t begin = 0, end = refLen, step = 1;

            VEC_TYPE bt_vec_init = VEC_SET1(0);
            VEC_TYPE d_dir_init = VEC_SET1(0x0001);
            VEC_TYPE i_dir_init = VEC_SET1(0x0003);


            VEC_TYPE *dir_mat_align = (VEC_TYPE *) aligned_alloc(sizeof(VEC_TYPE), segLen * refLen * sizeof(VEC_TYPE));

            //memset(dir_vec_align, 0, segLen*refLen*sizeof(VEC_TYPE));



            VEC_TYPE bt_vec_cache[4];

            int first_time = 0;


            /* outer loop to process the reference sequence */

            for (i = begin; LIKELY(i != end); i += step) {
                uint64_t cmp;
                VEC_TYPE e, vF = vNeg; /* Initialize F value to 0.
							   Any errors to vH values will be corrected in the Lazy_F loop.
							 */
                VEC_TYPE vH = pvHStore[segLen - 1];


                vH = VEC_SLL(vH, 2); /* Shift the 128-bit value in vH left by 2 byte. */

                /* Swap the 2 H buffers. */
                VEC_TYPE *pv = pvHLoad;

                VEC_TYPE vMaxColumn = vNeg; /* vMaxColumn is used to record the max values of column i. */

                const VEC_TYPE *vP = vProfile + ref[i] * segLen; /* Right part of the vProfile */
                pvHLoad = pvHStore;
                pvHStore = pv;


                /* inner loop to process the query sequence */
                for (l = 0; l < segLen_per_btvec; l++) {
                    VEC_TYPE bt_vec = bt_vec_init;
                    VEC_TYPE i_dir = i_dir_init;
                    VEC_TYPE d_dir = d_dir_init;
                    for (j = 0; j < 8; j++) {
                        e = VEC_LOAD(pvE + (l * 8) + j);

                        if (first_time) {
                            vH = vZero;
                            e = vNeg;
                        }
                        vH = VEC_ADDS(vH, VEC_LOAD(vP + (l * 8) + j));

                        /* Get max from vH, vE and vF. */




                        vH = VEC_MAX(vH, e);
                        vH = VEC_MAX(vH, vF);

                        VEC_TYPE mask = VEC_CMPEQ(vH, vF);

                        bt_vec = VEC_OR(bt_vec, VEC_AND(mask, i_dir));

                        mask = VEC_CMPEQ(vH, e);

                        bt_vec = VEC_OR(bt_vec, VEC_AND(mask, d_dir));

                        i_dir = VEC_SLL_ELE(i_dir, 2);
                        d_dir = VEC_SLL_ELE(d_dir, 2);

                        vMaxColumn = VEC_MAX(vMaxColumn, vH);

                        /* Save vH values. */

                        VEC_STORE(pvHStore + (l * 8) + j, vH);

                        /* Update vE value. */
                        //vH = _mm_subs_epu16(vH, vGapOE); /* saturation arithmetic, result >= 0 */
                        vH = VEC_SUBS(vH, vGapOE); /* saturation arithmetic, result >= 0 */
                        //e = _mm_subs_epu16(e, vGapE);
                        e = VEC_SUBS(e, vGapE);
                        e = VEC_MAX(e, vH);
                        VEC_STORE(pvE + (l * 8) + j, e);

                        /* Update vF value. */
                        //vF = _mm_subs_epu16(vF, vGapE);
                        vF = VEC_SUBS(vF, vGapE);
                        vF = VEC_MAX(vF, vH);

                        /* Load the next vH. */
                        vH = VEC_LOAD(pvHLoad + (l * 8) + j);
                    }

                    //VEC_STORE(bt_vec_cache_align + (l*sizeof(VEC_TYPE)), bt_vec);
                    bt_vec_cache[l] = bt_vec;
                }

                VEC_TYPE bt_vec = bt_vec_init;
                VEC_TYPE i_dir = i_dir_init;
                VEC_TYPE d_dir = d_dir_init;

                for (j = segLen_per_btvec * 8; j < segLen; j++) {
                    e = VEC_LOAD(pvE + (l * 8) + j);

                    if (first_time) {
                        vH = vZero;
                        e = vNeg;
                    }
                    vH = VEC_ADDS(vH, VEC_LOAD(vP + j));

                    /* Get max from vH, vE and vF. */
                    e = VEC_LOAD(pvE + j);

                    vH = VEC_MAX(vH, e);
                    vH = VEC_MAX(vH, vF);

                    VEC_TYPE mask = VEC_CMPEQ(vH, vF);

                    bt_vec = VEC_OR(bt_vec, VEC_AND(mask, i_dir));

                    mask = VEC_CMPEQ(vH, e);

                    bt_vec = VEC_OR(bt_vec, VEC_AND(mask, d_dir));

                    i_dir = VEC_SLL_ELE(i_dir, 2);
                    d_dir = VEC_SLL_ELE(d_dir, 2);

                    vMaxColumn = VEC_MAX(vMaxColumn, vH);

                    /* Save vH values. */

                    VEC_STORE(pvHStore + j, vH);

                    /* Update vE value. */
                    //vH = _mm_subs_epu16(vH, vGapOE); /* saturation arithmetic, result >= 0 */
                    vH = VEC_SUBS(vH, vGapOE); /* saturation arithmetic, result >= 0 */
                    //e = _mm_subs_epu16(e, vGapE);
                    e = VEC_SUBS(e, vGapE);
                    e = VEC_MAX(e, vH);
                    VEC_STORE(pvE + j, e);

                    /* Update vF value. */
                    //vF = _mm_subs_epu16(vF, vGapE);
                    vF = VEC_SUBS(vF, vGapE);
                    vF = VEC_MAX(vF, vH);

                    /* Load the next vH. */
                    vH = VEC_LOAD(pvHLoad + j);
                }



                //VEC_STORE(bt_vec_cache_align + (l*sizeof(VEC_TYPE)), bt_vec);
                bt_vec_cache[l] = bt_vec;


                /* Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
                for (k = 0; LIKELY(k < 8); ++k) {
                    vF = VEC_SLL(vF, 2);
                    for (l = 0; l < segLen_per_btvec; l++) {
                        VEC_TYPE i_dir = i_dir_init;
                        //VEC_TYPE bt_vec = VEC_LOAD(bt_vec_cache_align + (l*sizeof(VEC_TYPE)));
                        VEC_TYPE bt_vec = bt_vec_cache[l];
                        for (j = 0; j < 8; ++j) {
                            vH = VEC_LOAD(pvHStore + (l * 8) + j);
                            VEC_TYPE mask = VEC_CMPGT(vF, vH);
                            bt_vec = VEC_OR(bt_vec, VEC_AND(mask, i_dir));
                            i_dir = VEC_SLL_ELE(i_dir, 2);
                            vH = VEC_MAX(vH, vF);
                            vMaxColumn = VEC_MAX(vMaxColumn, vH); //newly added line
                            VEC_STORE(pvHStore + (l * 8) + j, vH);
                            //vH = _mm_subs_epu16(vH, vGapOE);
                            vH = VEC_SUBS(vH, vGapOE);
                            //vF = _mm_subs_epu16(vF, vGapE);
                            vF = VEC_SUBS(vF, vGapE);
                            if (UNLIKELY(!VEC_MOVEMASK(VEC_CMPGT(vF, vH)))) goto end;
                        }
                        bt_vec_cache[l] = bt_vec;
                        //VEC_STORE((segLen*(sizeof(VEC_TYPE)/sizeof(int16_t))*sizeof(int16_t)*i) + (l*sizeof(VEC_TYPE)), bt_vec);

                    }
                    VEC_TYPE i_dir = i_dir_init;
                    //VEC_TYPE bt_vec = VEC_LOAD(bt_vec_cache_align + (l*sizeof(VEC_TYPE)));
                    VEC_TYPE bt_vec = bt_vec_cache[l];
                    for (j = segLen_per_btvec * 8; j < segLen; j++) {
                        vH = VEC_LOAD(pvHStore + j);
                        VEC_TYPE mask = VEC_CMPGT(vF, vH);
                        bt_vec = VEC_OR(bt_vec, VEC_AND(mask, i_dir));
                        i_dir = VEC_SLL_ELE(i_dir, 2);
                        vH = VEC_MAX(vH, vF);
                        vMaxColumn = VEC_MAX(vMaxColumn, vH); //newly added line
                        VEC_STORE(pvHStore + j, vH);
                        //vH = _mm_subs_epu16(vH, vGapOE);
                        vH = VEC_SUBS(vH, vGapOE);
                        //vF = _mm_subs_epu16(vF, vGapE);
                        vF = VEC_SUBS(vF, vGapE);
                        if (UNLIKELY(!VEC_MOVEMASK(VEC_CMPGT(vF, vH)))) goto end;
                    }
                    bt_vec_cache[l] = bt_vec;
                }


                end:
                bt_vec_cache[l] = bt_vec;
                for (l = 0; l < segLen_per_btvec + (seg_left > 0 ? 1 : 0); l++) {
                    VEC_STORE(dir_mat_align + l + (i * segLen), bt_vec_cache[l]);
                }
                vMaxScore = VEC_MAX(vMaxScore, vMaxColumn);
                vTemp = VEC_CMPEQ(vMaxMark, vMaxScore);
                cmp = VEC_MOVEMASK(vTemp);
                if (cmp != (uint64_t) VEC_ALL_ONES) {
                    int16_t temp;
                    vMaxMark = vMaxScore;
                    max_of_vec(temp, vMaxScore);
                    vMaxScore = vMaxMark;
                    if (LIKELY(temp > max)) {
                        max = temp;
                        end_ref = i;
                        for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
                    }
                }
                first_time = 0;
            }

            /* Trace the alignment ending position on read. */
            int16_t *t = (int16_t *) pvHmax;
            int32_t column_len = segLen * (sizeof(VEC_TYPE) / sizeof(int16_t));
            for (i = 0; LIKELY(i < column_len); ++i, ++t) {
                int32_t temp;
                if (*t == max) {
                    temp = i / (sizeof(VEC_TYPE) / sizeof(int16_t)) + i % (sizeof(VEC_TYPE) / sizeof(int16_t)) * segLen;
                    if (temp < end_read) end_read = temp;
                }
            }

            free(pvHmax);
            free(pvE);
            free(pvHLoad);
            free(pvHStore);


            best->score = max;
            best->ref = end_ref + 1;
            best->read = end_read + 1;
            best->dir_mat = (uint8_t *) dir_mat_align;

            return best;
        }


        s_profile* ssw_init (const int8_t* read, const int32_t readLen, const int8_t* mat, const int32_t n) {
            s_profile *p = (s_profile *) calloc(1, sizeof(struct _profile));
            p->profile_word = 0;
            p->profile_word = qP_word(read, mat, readLen, n);
            p->read = read;
            p->mat = mat;
            p->readLen = readLen;
            p->n = n;
            return p;
        }

        void init_destroy (s_profile* p) {
            free(p->profile_word);
            free(p);
        }

        s_align* ssw_align (const s_profile* prof,
        const int8_t* ref,
        int32_t refLen,
        const int8_t weight_gapO,
        const int8_t weight_gapE) {

            alignment_end *best = 0;
            int32_t readLen = prof->readLen;
            s_align *r = (s_align *) calloc(1, sizeof(s_align));
            r->ref_begin = -1;
            r->read_begin = -1;
            // Find the alignment scores and ending positions
            best = ssw_int16(ref, refLen, readLen, weight_gapO, weight_gapE, prof->profile_word);


            r->score = best->score;
            r->ref_end = best->ref; // 0_based, always count from the input seq begin
            r->read_end = best->read;   // 0_based, count from the alignment begin (aligned length of the read)

            free(best->dir_mat);
            free(best);

            return r;
        }

        void align_destroy (s_align* a) {
            free(a);
        }

private:

        int16_t score;
        int32_t ref_end;     //0-based position
        int32_t read_end;    //alignment ending position on read, 0-based
        uint8_t *dir_mat;
        VEC_TYPE *profile_word;    // 0: none
        const int8_t *read;
        const int8_t *mat;
        int32_t readLen;
        int32_t n;
        int match[][];


        string cigar;
        stack<char> operation_stack;
        stack<int> count_stack;
        int MAX_LENGTH;
        int GAP_OPEN;
        int GAP_EXT;
        int max_i;
        int max_j;
        int deletions;
        int insertions;
        int similarity;
        double identity;
        char TYPE;
        int offset;
        int range_len;
        int mismatch_penalty, insertion_penalty;
        int CLIPPING_STRIGENCY;
        bool DEBUG = false;

}

