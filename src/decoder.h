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

#ifndef MARLIN_DECODER_H
#define MARLIN_DECODER_H

#include "dictionary.h"

namespace {


	class Decoder {

		using SourceSymbol = Dictionary::SourceSymbol;
		using MarlinSymbol = Dictionary::MarlinSymbol;
		
		const size_t shift;
		const size_t O;
		const size_t maxWordSize;
		
		std::vector<SourceSymbol> decoderTable;
		const SourceSymbol * const D;
		const SourceSymbol mostCommonSourceSymbol;

		
		template<typename T, size_t CO>
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
/*					if (((v32 - 0x01010101UL) & ~v32 & 0x80808080UL)) { // Fast test for zero

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
						
					} else { // Has no zeroes! hurray!*/
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
		}
		
	public:
		
		
		Decoder(const Dictionary &dict, const Configuration &) :
			shift(dict.alphabet.shift),
			O(dict.O),
			maxWordSize(dict.maxWordSize),
			decoderTable(dict.words.size()*(maxWordSize+1)),
			D(decoderTable.data()),
			mostCommonSourceSymbol(dict.alphabet.marlinSymbols.front().sourceSymbol) {
				
			
			SourceSymbol *d = &decoderTable.front();
			for (size_t i=0; i<dict.words.size(); i++) {
				for (size_t j=0; j<maxWordSize; j++)
					*d++ = (dict.words[i].size()>j ? dict.words[i][j] : SourceSymbol(0));
				*d++ = dict.words[i].size();
			}
		}
		
		size_t operator()(const uint8_t * const i8start, const uint8_t * const i8end, uint8_t * const o8start, uint8_t * const o8end) const {
			
			if (maxWordSize==3) {
				switch (O) {
					case   0: return decode8<uint32_t,0>(i8start, i8end, o8start, o8end);
					case   1: return decode8<uint32_t,1>(i8start, i8end, o8start, o8end);
					case   2: return decode8<uint32_t,2>(i8start, i8end, o8start, o8end);
					case   3: return decode8<uint32_t,3>(i8start, i8end, o8start, o8end);
					case   4: return decode8<uint32_t,4>(i8start, i8end, o8start, o8end);
				}
			}

			if (maxWordSize==7) {
				switch (O) {
					case   0: return decode8<uint64_t,0>(i8start, i8end, o8start, o8end);
					case   1: return decode8<uint64_t,1>(i8start, i8end, o8start, o8end);
					case   2: return decode8<uint64_t,2>(i8start, i8end, o8start, o8end);
					case   3: return decode8<uint64_t,3>(i8start, i8end, o8start, o8end);
					case   4: return decode8<uint64_t,4>(i8start, i8end, o8start, o8end);
				}
			}
			throw std::runtime_error ("unsupported maxWordSize");	
		}
	};
	
}

#endif //MARLIN_DECODER_H
