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

#ifdef __ARM_NEON // (M1)
#include "sse2neon.h"
#else // x86 (Intel)
#include <emmintrin.h>
#endif


#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/* Convert the coordinate in the scoring matrix into the coordinate in one line of the band. */
#define set_u(u, w, i, j) { int x=(i)-(w); x=x>0?x:0; (u)=(j)-x+1; }

/* Convert the coordinate in the direction matrix into the coordinate in one line of the band. */
#define set_d(u, w, i, j, p) { int x=(i)-(w); x=x>0?x:0; x=(j)-x; (u)=x*3+p; }

/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

typedef struct {
	int16_t score;
	int32_t ref;	 //0-based position
	int32_t read;    //alignment ending position on read, 0-based
} alignment_end;

typedef struct {
	uint32_t* seq;
	int32_t length;
} cigar;

struct _profile{
	__m128i* profile_byte;	// 0: none
	__m128i* profile_word;	// 0: none
	const int8_t* read;
	const int8_t* mat;
	int32_t readLen;
	int32_t n;
};

/* array index is an ASCII character value from a CIGAR,
   element value is the corresponding integer opcode between 0 and 8 */
const uint8_t encoded_ops[] = {
	0,         0,         0,         0,
	0,         0,         0,         0,
	0,         0,         0,         0,
	0,         0,         0,         0,
	0,         0,         0,         0,
	0,         0,         0,         0,
	0,         0,         0,         0,
	0,         0,         0,         0,
	0 /*   */, 0 /* ! */, 0 /* " */, 0 /* # */,
	0 /* $ */, 0 /* % */, 0 /* & */, 0 /* ' */,
	0 /* ( */, 0 /* ) */, 0 /* * */, 0 /* + */,
	0 /* , */, 0 /* - */, 0 /* . */, 0 /* / */,
	0 /* 0 */, 0 /* 1 */, 0 /* 2 */, 0 /* 3 */,
	0 /* 4 */, 0 /* 5 */, 0 /* 6 */, 0 /* 7 */,
	0 /* 8 */, 0 /* 9 */, 0 /* : */, 0 /* ; */,
	0 /* < */, 7 /* = */, 0 /* > */, 0 /* ? */,
	0 /* @ */, 0 /* A */, 0 /* B */, 0 /* C */,
	2 /* D */, 0 /* E */, 0 /* F */, 0 /* G */,
	5 /* H */, 1 /* I */, 0 /* J */, 0 /* K */,
	0 /* L */, 0 /* M */, 3 /* N */, 0 /* O */,
	6 /* P */, 0 /* Q */, 0 /* R */, 4 /* S */,
	0 /* T */, 0 /* U */, 0 /* V */, 0 /* W */,
	8 /* X */, 0 /* Y */, 0 /* Z */, 0 /* [ */,
	0 /* \ */, 0 /* ] */, 0 /* ^ */, 0 /* _ */,
	0 /* ` */, 0 /* a */, 0 /* b */, 0 /* c */,
	0 /* d */, 0 /* e */, 0 /* f */, 0 /* g */,
	0 /* h */, 0 /* i */, 0 /* j */, 0 /* k */,
	0 /* l */, 0 /* m */, 0 /* n */, 0 /* o */,
	0 /* p */, 0 /* q */, 0 /* r */, 0 /* s */,
	0 /* t */, 0 /* u */, 0 /* v */, 0 /* w */,
	0 /* x */, 0 /* y */, 0 /* z */, 0 /* { */,
	0 /* | */, 0 /* } */, 0 /* ~ */, 0 /*  */
};



static __m128i* qP_word (const int8_t* read_num,
				  const int8_t* mat,
				  const int32_t readLen,
				  const int32_t n) {

	int32_t segLen = (readLen + 7) / 8;
	__m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
	int16_t* t = (int16_t*)vProfile;
	int32_t nt, i, j;
	int32_t segNum;

	/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch */
	for (nt = 0; LIKELY(nt < n); nt ++) {
		for (i = 0; i < segLen; i ++) {
			j = i;
			for (segNum = 0; LIKELY(segNum < 8) ; segNum ++) {
				*t++ = j>= readLen ? 0 : mat[nt * n + read_num[j]];
				j += segLen;
			}
		}
	}
	return vProfile;
}

static alignment_end* sw_sse2_word (const int8_t* ref,
							 int32_t refLen,
							 int32_t readLen,
							 const uint8_t weight_gapO, /* will be used as - */
							 const uint8_t weight_gapE, /* will be used as - */
							 const __m128i* vProfile) {


#define max8(m, vm) (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 8)); \
					(vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 4)); \
					(vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 2)); \
					(m) = _mm_extract_epi16((vm), 0)

	int16_t max = SHRT_MIN;		                     /* the max alignment score */
	int32_t end_read = readLen - 1;
	int32_t end_ref = 0; /* 1_based best alignment ending point; Initialized as isn't aligned - 0. */
	int32_t segLen = (readLen + 7) / 8; /* number of segment */

	/* array to record the largest score of each reference position */
	int16_t* maxColumn = (int16_t*) calloc(refLen, 2);


	/* Define 16 byte 0 vector. */
	__m128i vNeg = _mm_set1_epi16(SHRT_MIN);

	__m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
	__m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
	__m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
	int itr;
	for (itr = 0; itr < segLen; itr++){
	    pvE[itr] = vNeg;
	}
	__m128i* pvHmax = (__m128i*) calloc(segLen, sizeof(__m128i));

	int32_t i, j, k;
	/* 16 byte insertion begin vector */
	__m128i vGapOE = _mm_set1_epi16(weight_gapO+weight_gapE);

	/* 16 byte insertion extension vector */
	__m128i vGapE = _mm_set1_epi16(weight_gapE);

	__m128i vMaxScore = vNeg; /* Trace the highest score of the whole SW matrix. */
	__m128i vMaxMark = vNeg; /* Trace the highest score till the previous column. */
	__m128i vTemp;
	int32_t begin = 0, end = refLen, step = 1;

	/* outer loop to process the reference sequence */

	for (i = begin; LIKELY(i != end); i += step) {
		int32_t cmp;
		__m128i e, vF = vNeg; /* Initialize F value to 0.
							   Any errors to vH values will be corrected in the Lazy_F loop.*/
		__m128i vH = pvHStore[segLen - 1];
		vH = _mm_slli_si128 (vH, 2); /* Shift the 128-bit value in vH left by 2 byte. */

		/* Swap the 2 H buffers. */
		__m128i* pv = pvHLoad;

		__m128i vMaxColumn = vNeg; /* vMaxColumn is used to record the max values of column i. */

		const __m128i* vP = vProfile + ref[i] * segLen; /* Right part of the vProfile */
		pvHLoad = pvHStore;
		pvHStore = pv;

		/* inner loop to process the query sequence */
		for (j = 0; LIKELY(j < segLen); j ++) {
			//vH = _mm_adds_epi16(vH, _mm_load_si128(vP + j));
			vH = _mm_adds_epi16(vH, _mm_load_si128(vP + j));

			/* Get max from vH, vE and vF. */
			e = _mm_load_si128(pvE + j);
			vH = _mm_max_epi16(vH, e);
			vH = _mm_max_epi16(vH, vF);
			vMaxColumn = _mm_max_epi16(vMaxColumn, vH);

			/* Save vH values. */
			_mm_store_si128(pvHStore + j, vH);

			/* Update vE value. */
			//vH = _mm_subs_epu16(vH, vGapOE); /* saturation arithmetic, result >= 0 */
			vH = _mm_subs_epi16(vH, vGapOE); /* saturation arithmetic, result >= 0 */
			//e = _mm_subs_epu16(e, vGapE);
			e = _mm_subs_epi16(e, vGapE);
			e = _mm_max_epi16(e, vH);
			_mm_store_si128(pvE + j, e);

			/* Update vF value. */
			//vF = _mm_subs_epu16(vF, vGapE);
			vF = _mm_subs_epi16(vF, vGapE);
			vF = _mm_max_epi16(vF, vH);

			/* Load the next vH. */
			vH = _mm_load_si128(pvHLoad + j);
		}

		/* Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
		for (k = 0; LIKELY(k < 8); ++k) {
			vF = _mm_slli_si128 (vF, 2);
			for (j = 0; LIKELY(j < segLen); ++j) {
				vH = _mm_load_si128(pvHStore + j);
				vH = _mm_max_epi16(vH, vF);
				vMaxColumn = _mm_max_epi16(vMaxColumn, vH); //newly added line
				_mm_store_si128(pvHStore + j, vH);
				//vH = _mm_subs_epu16(vH, vGapOE);
				vH = _mm_subs_epi16(vH, vGapOE);
				//vF = _mm_subs_epu16(vF, vGapE);
				vF = _mm_subs_epi16(vF, vGapE);
				if (UNLIKELY(! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH)))) goto end;
			}
		}

end:
		vMaxScore = _mm_max_epi16(vMaxScore, vMaxColumn);
		vTemp = _mm_cmpeq_epi16(vMaxMark, vMaxScore);
		cmp = _mm_movemask_epi8(vTemp);
		if (cmp != 0xffff) {
			int16_t temp;
			vMaxMark = vMaxScore;
			max8(temp, vMaxScore);
			vMaxScore = vMaxMark;

			if (LIKELY(temp > max)) {
				max = temp;
				end_ref = i;
				for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
			}
		}

		/* Record the max score of current column. */
		max8(maxColumn[i], vMaxColumn);
	}

	/* Trace the alignment ending position on read. */
	int16_t *t = (int16_t*)pvHmax;
	int32_t column_len = segLen * 8;
	for (i = 0; LIKELY(i < column_len); ++i, ++t) {
		int32_t temp;
		if (*t == max) {
			temp = i / 8 + i % 8 * segLen;
			if (temp < end_read) end_read = temp;
		}
	}

	free(pvHmax);
	free(pvE);
	free(pvHLoad);
	free(pvHStore);

	alignment_end* best = (alignment_end*) calloc(1, sizeof(alignment_end));
	best->score = max;
	best->ref = end_ref+1;
	best->read = end_read+1;

	free(maxColumn);
	return best;
}



s_profile* ssw_init (const int8_t* read, const int32_t readLen, const int8_t* mat, const int32_t n) {
	s_profile* p = (s_profile*)calloc(1, sizeof(struct _profile));
	p->profile_word = 0;
    p->profile_word = qP_word (read, mat, readLen, n);
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

	alignment_end* best = 0;
	int32_t readLen = prof->readLen;
	s_align* r = (s_align*)calloc(1, sizeof(s_align));
	r->ref_begin = -1;
	r->read_begin = -1;
	// Find the alignment scores and ending positions
	best = sw_sse2_word(ref, refLen, readLen, weight_gapO, weight_gapE, prof->profile_word);

	r->score = best->score;
	r->ref_end = best->ref; // 0_based, always count from the input seq begin
	r->read_end = best->read;   // 0_based, count from the alignment begin (aligned length of the read)

	free(best);

	return r;
}

void align_destroy (s_align* a) {
	free(a);
}

