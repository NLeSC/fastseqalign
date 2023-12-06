#ifndef AVX2_FUNCTIONS_H
#define AVX2_FUNCTIONS_H

#include <immintrin.h>

#define VEC_TYPE __m256i

#define VEC_SIZE_IN_BITS 8

#define VEC_ALL_ONES 0xFFFFFFFF

#define VEC_SET1(val)  _mm256_set1_epi16(val)

#define VEC_SLL(v, count) _mm256_alignr_epi8(v, _mm256_permute2x128_si256(v, v, _MM_SHUFFLE(0,0,3,0) ),16-count);

#define VEC_SLL_ELE(v, count) _mm256_slli_epi16(v, count)

#define VEC_LOAD(addr) _mm256_load_si256((__m256i*)addr)

#define VEC_STORE(addr, v) _mm256_store_si256((__m256i*)addr, v)

#define VEC_ADDS(v1, v2)  _mm256_adds_epi16(v1, v2)

#define VEC_SUBS(v1, v2) _mm256_subs_epi16(v1, v2)  // v1 - v2


#define VEC_MAX(v1, v2) _mm256_max_epi16(v1, v2)


#define VEC_CMPEQ(v1, v2) _mm256_cmpeq_epi16(v1, v2)

#define VEC_CMPGT(v1, v2) _mm256_cmpgt_epi16(v1, v2)

#define VEC_AND(v1, v2) _mm256_and_si256 (v1, v2)

#define VEC_OR(v1, v2) _mm256_or_si256 (v1, v2)


#define VEC_MOVEMASK(v) _mm256_movemask_epi8(v)

#define max_of_vec(m, vm) (vm) = _mm256_max_epi16((vm), _mm256_permute2x128_si256((vm), (vm), 0x81));\
                    (vm) = _mm256_max_epi16((vm), _mm256_srli_si256((vm), 8)); \
					(vm) = _mm256_max_epi16((vm), _mm256_srli_si256((vm), 4)); \
					(vm) = _mm256_max_epi16((vm), _mm256_srli_si256((vm), 2)); \
					(m) = _mm256_extract_epi16((vm), 0)


#endif