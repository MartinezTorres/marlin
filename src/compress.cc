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

#include "dictionary.hpp"

namespace {

	__attribute__ ((target ("bmi2")))
	void shift8(const MarlinDictionary &dict, uint8_t* dst, const uint8_t* src, const size_t srcSize) {
		
		uint64_t mask=0;
		for (size_t i=0; i<8; i++)
			mask |= ((1ULL<<dict.shift)-1)<<(8ULL*i);

		const uint64_t *i64    = (const uint64_t *)src;
		const uint64_t *i64end = (const uint64_t *)(src+srcSize);

		while (i64 != i64end) {
			*(uint64_t *)dst = _pext_u64(*i64++, mask);
			dst += dict.shift;
		}
	}
	
	
}


//std::shared_ptr<std::vector<MarlinDictionary::WordIdx>> MarlinDictionary::buildCompressorTable() const

ssize_t MarlinDictionary::compress(uint8_t* dst, size_t dstCapacity, const uint8_t* src, size_t srcSize) const {

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


		  uint8_t *o8    = dst;
	const uint8_t *i8    = src;
	const uint8_t *i8end = src + srcSize;

	// Special case: if srcSize is not multiple of 8, we force it to be.
	while ( (i8end-i8) % 8 != 0) {
		*o8++ = *i8++;
	}


	// Encode Marlin, with rare symbols preceded by an empty word
	{
		
		// if the encoder produces a size larger than this, it is simply better to store the block uncompressed.
		size_t maxTargetSize = std::max(size_t(8), srcSize-srcSize*shift/8) - 8;

		WordIdx j = 0; 
		while (i8<i8end) {				
			
			SourceSymbol ss = *i8++;
			
			MarlinSymbol ms = Source2JumpTableShifted[ss>>shift];
			bool isRareSymbol = ms==nMarlinSymbols;
			if (isRareSymbol) {
				if (j) *o8++ = j; // Finish current word, if any;
				*o8++ = j = 0;
				*o8++ = ss & ((1<<shift)-1);
				maxTargetSize-=3;
				if (i8end-i8 > maxTargetSize) { // Just encode the block uncompressed.
					memcpy(dst, src, srcSize);
					return srcSize;
				}
				continue;
			}
			
			WordIdx jOld = j;
			j = jump(j, ms);
			
			if (j & FLAG_NEXT_WORD) 
				*o8++ = jOld & 0xFF;
		}
		if (j) *o8++ = j;
	}

	// Encode residuals
	shift8(dict, o8, src, srcSize);
	
	return o8 - dst; 
}
