/*
 *  ssw.h
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 1.2.3
 *	Last revision by Mengyao Zhao on 2022-Apr-15.
 *
 */

#ifndef SSW_H
#define SSW_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#ifdef __ARM_NEON // (M1)
#include "sse2neon.h"
#else // x86 (Intel)
#include <emmintrin.h>
#endif


#ifdef __cplusplus
extern "C" {
#endif	// __cplusplus

#define MAPSTR "MIDNSHP=X"
#ifndef BAM_CIGAR_SHIFT
#define BAM_CIGAR_SHIFT 4u
#endif

extern const uint8_t encoded_ops[];

/*!	@typedef	structure of the query profile	*/
struct _profile;
typedef struct _profile s_profile;

/*!	@typedef	structure of the alignment result
    @field	score	the best alignment score
    @field	ref_begin	0-based best alignment beginning position on reference;	ref_begin1 = -1 when the best alignment beginning
						position is not available
    @field	ref_end 	0-based best alignment ending position on reference
    @field	read_begin	0-based best alignment beginning position on read; read_begin1 = -1 when the best alignment beginning
						position is not available
    @field	read_end	0-based best alignment ending position on read
*/
typedef struct {
	int16_t score;
	int32_t ref_begin;
	int32_t ref_end;
	int32_t	read_begin;
	int32_t read_end;
} s_align;


/*!	@function	Create the query profile using the query sequence.
	@param	read	pointer to the query sequence; the query sequence needs to be numbers
	@param	readLen	length of the query sequence
	@param	mat	pointer to the substitution matrix; mat needs to be corresponding to the read sequence
	@param	n	the square root of the number of elements in mat (mat has n*n elements)
	@param	score_size	estimated Smith-Waterman score; if your estimated best alignment score is surely < 255 please set 0; if
						your estimated best alignment score >= 255, please set 1; if you don't know, please set 2
	@return	pointer to the query profile structure
	@note	example for parameter read and mat:
			If the query sequence is: ACGTATC, the sequence that read points to can be: 1234142
			Then if the penalty for match is 2 and for mismatch is -2, the substitution matrix of parameter mat will be:
			//A  C  G  T
			  2 -2 -2 -2 //A
			 -2  2 -2 -2 //C
			 -2 -2  2 -2 //G
			 -2 -2 -2  2 //T
			mat is the pointer to the array {2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2}
*/
s_profile* ssw_init (const int8_t* read, const int32_t readLen, const int8_t* mat, const int32_t n);

/*!	@function	Release the memory allocated by function ssw_init.
	@param	p	pointer to the query profile structure
*/
void init_destroy (s_profile* p);

// @function	ssw alignment.
/*!	@function	Do Striped Smith-Waterman alignment.
	@param	prof	pointer to the query profile structure
	@param	ref	pointer to the target sequence; the target sequence needs to be numbers and corresponding to the mat parameter of
				function ssw_init
	@param	refLen	length of the target sequence
	@param	weight_gapO	the absolute value of gap open penalty
	@param	weight_gapE	the absolute value of gap extension penalty
	@return	pointer to the alignment result structure
*/
s_align* ssw_align (const s_profile* prof,
					const int8_t* ref,
					int32_t refLen,
					const int8_t weight_gapO,
					const int8_t weight_gapE);

/*!	@function	Release the memory allocated by function ssw_align.
	@param	a	pointer to the alignment result structure
*/
void align_destroy (s_align* a);



#ifdef __cplusplus
}
#endif	// __cplusplus

#endif	// SSW_H