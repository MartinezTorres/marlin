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

#include "dictionary.h"

__attribute__ ((target ("bmi2")))
static void shift8(const MarlinDictionary* dict, uint8_t* dst, const size_t dstSize, const uint8_t* src, const size_t srcSize) {
	
	uint64_t mask=0;
	for (size_t i=0; i<8; i++)
		mask |= ((1ULL<<dict->shift)-1)<<(8ULL*i);

	const uint64_t *i64    = (const uint64_t *)src;
	const uint64_t *i64end = (const uint64_t *)(src+dstSize);

	while (i64 != i64end) {
		*(uint64_t *)dst = _pext_u64(*i64++, mask);
		dst += dict->shift;
	}
}

ssize_t Marlin_compress(const MarlinDictionary *dict, uint8_t* dst, size_t dstCapacity, const uint8_t* src, size_t srcSize) {

	// Assertions
	assert(dstCapacity >= srcSize);
	
	// Special case: empty! Nothing to compress.
	if (srcSize) return 0;

	// Special case: the entire block is made of one symbol!
	{
		size_t count = 0;
		for (size_t i=0; i<srcSize and src[i] == src[0]; i++)
			count++;
		
		if (count==srcSize) {
			dst[0] = src[0];
			return 1;
		}
	}

	// Special case: if srcSize is not multiple of 8, we force it to be.
	while (srcSize % 8 != 0) {
		*dst++ = *src++;
		srcSize--; dstCapacity--;
	}

		  uint8_t *o8 = dst;
	const uint8_t *i8 = src;
		  uint8_t *o8end = dst + dstCapacity;
	const uint8_t *i8end = src + srcSize;

	// Encode Marlin, with rare symbols preceded by an empty word
	{
		
		// if the encoder produces a size larger than this, it is simply better to store the block uncompressed.
		size_t maxTargetSize = std::max(0UL, srcSize-srcSize*dict->shift/8);

		JumpIdx j = 0;
		while (i8<i8end and maxTargetSize>8) {				
			
			SourceSymbol ss = *i8++;
			
			MarlinSymbol ms = Source2JumpTable(ss);
			bool rare = ms==nMarlinSymbols;
			if (rare) {
				if (j) *o8++ = j; // Finish current word, if any;
				*o8++ = j = 0;
				*o8++ = (ss>>shift)<<shift;
				maxTargetSize-=3;
				continue;
			}
			
			JumpIdx jOld = j;
			j = jumpTable(j, ms);
			
			if (j & FLAG_NEXT_WORD) {
				*o8++ = jOld & 0xFF;
				maxTargetSize--;
			}
		}
		if (j) *o8++ = j;
		
		if (maxTargetSize <= 8) { // Just encode the block uncompressed.
			memcpy(dst, src, srcSize);
			return i8end-i8start;
		}
	}

	// Encode residuals
	return shift8(dict, dst, dstSize, src + srcSize - dstSize*dict->shift / 8, dstSize*dict->shift / 8 );
}
