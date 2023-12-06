#ifndef SSE2_FUNCTIONS_H
#define SSE2_FUNCTIONS_H

#ifdef __ARM_NEON // (M1)
#include "sse2neon.h"
#else // x86 (Intel)
#include <emmintrin.h>
#endif


#define VEC_TYPE __m128i

#define VEC_SIZE_IN_BITS 4

#define VEC_ALL_ONES 0xFFFF

#define VEC_SET1(val)  _mm_set1_epi16(val)

#define VEC_SLL(v, count) _mm_slli_si128(v, count)

#define VEC_SLL_ELE(v, count) _mm_slli_epi16(v, count)

#define VEC_LOAD(addr) _mm_load_si128((__m128i*)addr)

#define VEC_STORE(addr, v) _mm_store_si128((__m128i*)addr, v)

#define VEC_ADDS(v1, v2)  _mm_adds_epi16(v1, v2)

#define VEC_SUBS(v1, v2) _mm_subs_epi16(v1, v2)  // v1 - v2


#define VEC_MAX(v1, v2) _mm_max_epi16 (v1, v2)


#define VEC_CMPEQ(v1, v2) _mm_cmpeq_epi16(v1, v2)

#define VEC_CMPGT(v1, v2) _mm_cmpgt_epi16(v1, v2)

#define VEC_AND(v1, v2) _mm_and_si128 (v1, v2)

#define VEC_OR(v1, v2) _mm_or_si128 (v1, v2)


#define VEC_MOVEMASK(v) _mm_movemask_epi8(v)

#define max_of_vec(m, vm) (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 8)); \
					(vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 4)); \
					(vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 2)); \
					(m) = _mm_extract_epi16((vm), 0)


#endif