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

#include <marlin.h>

#include <cstring>
#include <algorithm>
#include <cassert>


using namespace marlin;

namespace {

template<typename TSource, typename MarlinIdx>
struct TDecompress : TMarlin<TSource,MarlinIdx> {
	
	using typename TMarlin<TSource, MarlinIdx>::Word;
	using typename TMarlin<TSource, MarlinIdx>::CompressorTableIdx;
	using typename TMarlin<TSource, MarlinIdx>::MarlinSymbol;

	using TMarlin<TSource, MarlinIdx>::K;
	using TMarlin<TSource, MarlinIdx>::O;
	using TMarlin<TSource, MarlinIdx>::shift;
	using TMarlin<TSource, MarlinIdx>::maxWordSize;
	using TMarlin<TSource, MarlinIdx>::conf;
	
	using TMarlin<TSource, MarlinIdx>::words;
	using TMarlin<TSource, MarlinIdx>::sourceAlphabet;
	using TMarlin<TSource, MarlinIdx>::marlinAlphabet;
	using TMarlin<TSource, MarlinIdx>::decompressorTablePointer;
	
	
	__attribute__ ((target ("bmi2")))
	ssize_t shift8(View<const uint8_t> src, View<TSource> dst) const {
		
		// Decode residuals
		uint64_t mask=0;
		for (size_t i=0; i<8; i++)
			mask |= ((1ULL<<shift)-1)<<(8ULL*i);
		
		const uint8_t *i8 = reinterpret_cast<const uint8_t *>(src.start);
		uint64_t *o64    = reinterpret_cast<uint64_t *>(dst.start);
		uint64_t *o64end = reinterpret_cast<uint64_t *>(dst.end);

		while (o64 != o64end) {
			*o64++ += _pdep_u64(*reinterpret_cast<const uint64_t *>(i8), mask);
			i8 += shift;
		}
		
		return reinterpret_cast<TSource *>(o64) - dst.start;
	}

		
	template<typename T>
	size_t decompress8(View<const uint8_t> src, View<TSource> dst) const {
		
		const uint8_t *i8    = src.start;
			  TSource *o8    = dst.start;

					const uint32_t overlappingMask = (1<<(8+O))-1;
	//	constexpr const uint32_t overlappingMask = (1<<(8+CO))-1;
					constexpr const T clearSizeMask = T(-1)>>8;
	//	constexpr const T clearSizeMask = 0;
		uint64_t value = 0;

		auto D = decompressorTablePointer;

		while (i8<src.end-9) {
			
			uint32_t v32 = (*(const uint32_t *)i8);
			/*if (((v32 - 0x01010101UL) & ~v32 & 0x80808080UL)) { // Fast test for zero

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
				
			} else { // Has no zeroes! hurray! */
				i8+=4;
				//clearSizeMask = 0;
				value = (value<<32) +  __builtin_bswap32(v32);
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
		
		while (i8<src.end) {
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
		// if (endMarlin-i8 != 0) std::cerr << " {" << endMarlin-i8 << "} "; // SOLVED! PROBLEM IN THE CODE
		// if (o8end-o8 != 0) std::cerr << " [" << o8end-o8 << "] "; // SOLVED! PROBLEM IN THE CODE

		return dst.nElements();
	}



	size_t decompressSlow(View<const uint8_t> src, View<TSource> dst) const {
		
		const uint8_t *i8 = src.start;
			  TSource *o8 = dst.start;
		
		uint64_t value = 0;
		uint64_t valueBits = O;
		while (i8 < src.end) {

			while (valueBits < K+O) {
				value += uint64_t(*i8++) << (64-8-valueBits);
				valueBits += 8;
			}
			
			size_t wordIdx = value >> (64-(K+O));
			if (not wordIdx) { printf("%d %d\n", i8-src.start, wordIdx); exit(-1); }
			
			value = value << K;
			valueBits -= K;
			
			for (auto &&c : words[wordIdx])
				*o8++ = marlinAlphabet[c].sourceSymbol;
		}
		
		return dst.nElements();
	}



	std::unique_ptr<std::vector<TSource>> buildDecompressorTable() const {
		
		auto ret = std::make_unique<std::vector<TSource>>(words.size()*(maxWordSize+1));
		
		TSource *d = &ret->front();
		for (size_t i=0; i<words.size(); i++) {
			for (size_t j=0; j<maxWordSize; j++)
				*d++ = (words[i].size()>j ? marlinAlphabet[words[i][j]].sourceSymbol : TSource(0));
			*d++ = words[i].size();
		}
		return ret;
	}


	ssize_t decompress(View<const uint8_t> src, View<TSource> dst) const {

		// Special case: empty block!
		if (dst.nBytes() == 0 or src.nBytes() == 0) {
			if (src.nBytes() or dst.nBytes()) return -1; // TODO: Error code
			return 0;
		}

		// Special case: the entire block is smaller than the size of a single symbol!
		if (src.nBytes() < sizeof(TSource)) return -1;
		
		// Special case: the entire block is made of one symbol!
		if (src.nBytes() == sizeof(TSource)) {
			TSource s = *reinterpret_cast<const TSource *&>(src.start);
			for (size_t i=0; i<dst.nElements(); i++)
				dst.start[i] = s;
			return dst.nElements();
		}

		// Special case, same size! this means the block is uncompressed.
		if (dst.nBytes() == src.nBytes()) {
			memcpy(dst.start, src.start, src.nBytes());
			return dst.nElements();
		}

		// Special case: if dstSize is not multiple of 8, we force it to be.
		size_t padding = 0;
		while ( dst.nBytes() % 8 != 0) {
			
			*dst.start++ = *reinterpret_cast<const TSource *&>(src.start)++;			
			padding += sizeof(TSource);
		}
		
		// Initialization, which might be optional
		for (size_t i=0; i<dst.nElements(); i++)
			dst.start[i] = marlinAlphabet.front().sourceSymbol;

		if (false) {
			decompressSlow(src, dst);
		} else if (K==8 and maxWordSize==3) {
			decompress8<uint32_t>(src, dst);
		} else if (K==8 and maxWordSize==7) {
			decompress8<uint64_t>(src, dst);
		} else {
			decompressSlow(src, dst);
		}
				
		return padding + shift8(src, dst);
	}
};
}

template<typename TSource, typename MarlinIdx>
std::unique_ptr<std::vector<TSource>> TMarlin<TSource,MarlinIdx>::buildDecompressorTable() const {
	return static_cast<const TDecompress<TSource,MarlinIdx> *>(this)->buildDecompressorTable();
}

template<typename TSource, typename MarlinIdx>
ssize_t TMarlin<TSource,MarlinIdx>::decompress(View<const uint8_t> src, View<TSource> dst) const {
	return static_cast<const TDecompress<TSource,MarlinIdx> *>(this)->decompress(src,dst);
}

////////////////////////////////////////////////////////////////////////
//
// Explicit Instantiations
#include "instantiations.h"
INSTANTIATE_MEMBER(buildDecompressorTable() const -> std::unique_ptr<std::vector<typename TMarlin::TSource_Type>>)	
INSTANTIATE_MEMBER(decompress(View<const uint8_t> src, View<typename TMarlin::TSource_Type> dst) const -> ssize_t)	

