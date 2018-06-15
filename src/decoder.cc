/***********************************************************************

Marlin: A Fast Entropy Codec

MIT License

Copyright (c) 2018 Manuel Martinez Torres

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

***********************************************************************/

//#include "decode.h"
#include "dictionary.hpp"



/*		template<typename T, size_t CO>
		__attribute__ ((target ("bmi2")))
		size_t decode8(const uint8_t * const i8start, const uint8_t * const i8end, uint8_t * const o8start, uint8_t * const o8end) const {
			
			      uint8_t *o8 = o8start;
			const uint8_t *i8 = i8start;
			
			// Special case, same size! this means the block is uncompressed.
			if (i8end-i8start == o8end-o8start) {
				memcpy(o8start,i8start,i8end-i8start);
				return o8end-o8start;
			}

			// Special case, size 1! this means the block consists all of just one symbol.
			if (i8end-i8start == 1) {
				memset(o8start,*i8start,o8end-o8start);
				return o8end-o8start;
			}

			memset(o8start,mostCommonSourceSymbol,o8end-o8start);
//			return o8end-o8start;
			
			// Decode the Marlin Section
			{

				const uint8_t *endMarlin = i8end - (o8end-o8start)*shift/8;

//				const uint32_t overlappingMask = (1<<(8+O))-1;
				constexpr const uint32_t overlappingMask = (1<<(8+CO))-1;
//				constexpr const T clearSizeMask = T(-1)>>8;
				constexpr const T clearSizeMask = 0;
				uint64_t value = 0;

				while (i8<endMarlin-9) {
					
					uint32_t v32 = (*(const uint32_t *)i8);
					if (((v32 - 0x01010101UL) & ~v32 & 0x80808080UL)) { // Fast test for zero

						uint8_t in = *i8++;
						if (in==0) {
							*o8++ = *i8++;
							value = (value<<8) + 0;
						} else {
							value = (value<<8) + in;
							T v = ((const T *)D)[value & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}
						
						in = *i8++;
						if (in==0) {
							*o8++ = *i8++;
							value = (value<<8) + 0;
						} else {
							value = (value<<8) + in;
							T v = ((const T *)D)[value & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}
						
						in = *i8++;
						if (in==0) {
							*o8++ = *i8++;
							value = (value<<8) + 0;
						} else {
							value = (value<<8) + in;
							T v = ((const T *)D)[value & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}
						
						in = *i8++;
						if (in==0) {
							*o8++ = *i8++;
							value = (value<<8) + 0;
						} else {
							value = (value<<8) + in;
							T v = ((const T *)D)[value & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}
						
					} else { // Has no zeroes! hurray!
						i8+=4;
						//clearSizeMask = 0;
						value = (value<<32) +  v32; //__builtin_bswap32(v32);
						{
							T v = ((const T *)D)[(value>>24) & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
							
						}

						{
							T v = ((const T *)D)[(value>>16) & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}

						{
							T v = ((const T *)D)[(value>>8) & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}

						{
							T v = ((const T *)D)[value & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}
					//}
				}
				
				while (i8<endMarlin) {
					uint8_t in = *i8++;
					if (in==0) {
						*o8++ = *i8++;
					} else {
						value = (value<<8) + in;
						const T *v = &((const T *)D)[value & overlappingMask];
						memcpy(o8, v, std::min(sizeof(T)-1,size_t(*v >> ((sizeof(T)-1)*8))));
						o8 += *v >> ((sizeof(T)-1)*8);
					}
				}				
				//if (endMarlin-i8 != 0) std::cerr << " {" << endMarlin-i8 << "} "; // SOLVED! PROBLEM IN THE CODE
				//if (o8end-o8 != 0) std::cerr << " [" << o8end-o8 << "] "; // SOLVED! PROBLEM IN THE CODE
			}

			// Decode residuals
			if (shift) {
				uint64_t mask=0;
				for (size_t i=0; i<8; i++)
					mask |= ((1ULL<<shift)-1)<<(8ULL*i);
				
				uint64_t *o64    = (uint64_t *)o8start;
				uint64_t *o64end = (uint64_t *)o8end;

				while (o64 != o64end) {
					*o64++ += _pdep_u64(*(const uint64_t *)i8, mask);
					i8 += shift;
				}
			}
			return o8end-o8start;
		}*/
template<typename T>
static ssize_t decodeMarlin8(const uint8_t* const D, const uint8_t O, uint8_t shift, uint8_t* dst, const size_t dstSize, const uint8_t* src, const size_t srcSize) {
	
		  uint8_t *o8 = dst;
	const uint8_t *i8 = src;
	const uint8_t *i8end = src + srcSize;
	
	// Decode the Marlin Section
	{

		const uint8_t *endMarlin = src + dstSize - dstSize * shift / 8;

		const uint32_t overlappingMask = (1<<(8+O))-1;
		constexpr const T clearSizeMask = T(-1)>>8; // TODO: Check this value
		uint64_t value = 0;

		while (i8<endMarlin-9) {
			
			uint32_t v32 = (*(const uint32_t *)i8);
			i8+=4;
			value = (value<<32) +  v32; //__builtin_bswap32(v32);
			{
				T v = ((const T *)D)[(value>>24) & overlappingMask];
				*((T *)o8) = v & clearSizeMask;
				o8 += v >> ((sizeof(T)-1)*8);
				
			}

			{
				T v = ((const T *)D)[(value>>16) & overlappingMask];
				*((T *)o8) = v & clearSizeMask;
				o8 += v >> ((sizeof(T)-1)*8);
			}

			{
				T v = ((const T *)D)[(value>>8) & overlappingMask];
				*((T *)o8) = v & clearSizeMask;
				o8 += v >> ((sizeof(T)-1)*8);
			}

			{
				T v = ((const T *)D)[value & overlappingMask];
				*((T *)o8) = v & clearSizeMask;
				o8 += v >> ((sizeof(T)-1)*8);
			}
		}
		
		while (i8<endMarlin) {
			uint8_t in = *i8++;
			if (in==0) {
				*o8++ = *i8++;
			} else {
				value = (value<<8) + in;
				const T *v = &((const T *)D)[value & overlappingMask];
				memcpy(o8, v, std::min(sizeof(T)-1,size_t(*v >> ((sizeof(T)-1)*8))));
				o8 += *v >> ((sizeof(T)-1)*8);
			}
		}				
		//if (endMarlin-i8 != 0) std::cerr << " {" << endMarlin-i8 << "} "; // SOLVED! PROBLEM IN THE CODE
		//if (o8end-o8 != 0) std::cerr << " [" << o8end-o8 << "] "; // SOLVED! PROBLEM IN THE CODE
	}

	return o8-dst;
}

__attribute__ ((target ("bmi2")))
static ssize_t decodeShift(uint8_t shift, uint8_t* dst, const size_t dstSize, const uint8_t* src) {
	
	// Decode residuals
	uint64_t mask=0;
	for (size_t i=0; i<8; i++)
		mask |= ((1ULL<<shift)-1)<<(8ULL*i);
	
	uint64_t *o64    = (uint64_t *)dst;
	uint64_t *o64end = (uint64_t *)(dst+dstSize);

	while (o64 != o64end) {
		*o64++ += _pdep_u64(*(const uint64_t *)src, mask);
		src += shift;
	}
	
	return dstSize;
}

ssize_t Marlin_decompress(const MarlinDictionary *dict, uint8_t* dst, size_t dstSize, const uint8_t* src, size_t srcSize) {
	
	// Special case, same size! this means the block is uncompressed.
	if (dstSize == srcSize) {
		memcpy(dst, src, dstSize);
		return dstSize;
	}

	// Special case: the entire block is made of one symbol!
	if (srcSize == 1) {
		memset(dst, *src, dstSize);
		return dstSize;
	}

	// Special case: if dstSize is not multiple of 8, we force it to be.
	while (dstSize % 8 != 0) {
		*dst++ = *src++;
		dstSize--; srcSize--;
	}
	
	// Initialization, which might be optional
	memset(dst, dict->mostCommonSourceSymbol, dstSize);
	
	if (dict->K==8 and dict->maxWordSize==3) {
		decodeMarlin8<uint32_t>(dict->decoderTablePointer, dict->O, dict->shift, dst, dstSize, src, srcSize);
	} else if (dict->K==8 and dict->maxWordSize==7) {
		decodeMarlin8<uint64_t>(dict->decoderTablePointer, dict->O, dict->shift, dst, dstSize, src, srcSize);
	} else {
		return -1;
		//decodeGeneric(dict, dst, dstSize, src, srcSize);
	}
	
	return decodeShift(dict->shift, dst, dstSize, src + srcSize - dstSize*dict->shift / 8);
}


