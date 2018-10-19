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
__attribute__ ((target ("bmi2")))
ssize_t shift8(const TMarlinDecompress<TSource,MarlinIdx> &decompressor, View<const uint8_t> src, View<TSource> dst) {
	
	// Decode residuals
	uint64_t mask=0;
	for (size_t i=0; i<8; i++)
		mask |= ((1ULL<<decompressor.shift)-1)<<(8ULL*i);
	
	const uint8_t *i8 = reinterpret_cast<const uint8_t *>(src.start);
	uint64_t *o64    = reinterpret_cast<uint64_t *>(dst.start);
	uint64_t *o64end = reinterpret_cast<uint64_t *>(dst.end);

	while (o64 != o64end) {
		*o64++ |= _pdep_u64(*reinterpret_cast<const uint64_t *>(i8), mask);
		i8 += decompressor.shift;
	}
	
	return reinterpret_cast<TSource *>(o64) - dst.start;
}

template<typename T, typename TSource, typename MarlinIdx>
size_t decompress8_skip(
	const TMarlinDecompress<TSource,MarlinIdx> &decompressor, 
	View<const uint8_t> src, 
	View<TSource> dst) {
	
	const uint8_t *i8    = src.start;
		  TSource *o8    = dst.start;

	const uint32_t overlappingMask = (1<<(8+decompressor.O))-1;
	uint64_t value = 0;

	auto D = decompressor.decompressorTablePointer;

	while (i8<src.end-10) {
		
		uint32_t v32 = (*(const uint32_t *)i8);
		i8 += sizeof(uint32_t);
		value = (value<<32) +  __builtin_bswap32(v32);
		{
			T v = ((const T *)D)[(value>>24) & overlappingMask];
			*((T *)o8) = v;
			o8 += v >> ((sizeof(T)-1)*8);
			
		}

		{
			T v = ((const T *)D)[(value>>16) & overlappingMask];
			*((T *)o8) = v;
			o8 += v >> ((sizeof(T)-1)*8);
		}

		{
			T v = ((const T *)D)[(value>>8) & overlappingMask];
			*((T *)o8) = v;
			o8 += v >> ((sizeof(T)-1)*8);
		}

		{
			T v = ((const T *)D)[value & overlappingMask];
			*((T *)o8) = v;
			o8 += v >> ((sizeof(T)-1)*8);
		}
	}
	
	while (i8<src.end) {
		value = (value<<8) + *i8++;
		{
			T v = ((const T *)D)[value & overlappingMask];
			size_t sz = v >> ((sizeof(T)-1)*8);
			memcpy(o8, &v, sz);
//			*((T *)o8) = v;
			
			o8 += sz;
		}
	}				
	// if (endMarlin-i8 != 0) std::cerr << " {" << endMarlin-i8 << "} "; // SOLVED! PROBLEM IN THE CODE
	// if (o8end-o8 != 0) std::cerr << " [" << o8end-o8 << "] "; // SOLVED! PROBLEM IN THE CODE

	return dst.nElements();
}

template<typename T, size_t KK, typename TSource, typename MarlinIdx>
size_t decompressKK(
	const TMarlinDecompress<TSource,MarlinIdx> &decompressor, 
	View<const uint8_t> src, 
	View<TSource> dst) {
	
	const uint8_t *i8    = src.start;
		  TSource *o8    = dst.start;

	const uint64_t overlappingMask = (1<<(decompressor.K+decompressor.O))-1;
	const T clearSizeMask = T(-1)>>8;
	const T clearSizeOverlay = T(decompressor.marlinMostCommonSymbol) << ((sizeof(T)-1)*8);
	uint64_t value = 0;

	auto D = decompressor.decompressorTablePointer;

	constexpr size_t INCREMENT = KK<8?KK:KK/2;
	constexpr size_t INCREMENTSHIFT = INCREMENT*8;

	while (i8<src.end-INCREMENT-20) {

		uint64_t vRead = 
			(INCREMENT<=4?
				__builtin_bswap32(*(const uint32_t *)i8):
				__builtin_bswap64(*(const uint64_t *)i8));
		i8 += INCREMENT;
		value = (value<<INCREMENTSHIFT) +  (vRead>>((INCREMENT<=4?32:64)-INCREMENTSHIFT));

		if (KK<8) {
			{
				T v = ((const T *)D)[(value>>(7*(KK%8))) & overlappingMask];
				*((T *)o8) = (v & clearSizeMask) + clearSizeOverlay;
				o8 += v >> ((sizeof(T)-1)*8);
				
			}

			{
				T v = ((const T *)D)[(value>>(6*(KK%8))) & overlappingMask];
				*((T *)o8) = (v & clearSizeMask) + clearSizeOverlay;
				o8 += v >> ((sizeof(T)-1)*8);
			}

			{
				T v = ((const T *)D)[(value>>(5*(KK%8))) & overlappingMask];
				*((T *)o8) = (v & clearSizeMask) + clearSizeOverlay;
				o8 += v >> ((sizeof(T)-1)*8);
			}

			{
				T v = ((const T *)D)[(value>>(4*(KK%8))) & overlappingMask];
				*((T *)o8) = (v & clearSizeMask) + clearSizeOverlay;
				o8 += v >> ((sizeof(T)-1)*8);
			}
		}

		{
			{
				T v = ((const T *)D)[(value>>(3*KK)) & overlappingMask];
				*((T *)o8) = (v & clearSizeMask) + clearSizeOverlay;
				o8 += v >> ((sizeof(T)-1)*8);
				
			}

			{
				T v = ((const T *)D)[(value>>(2*KK)) & overlappingMask];
				*((T *)o8) = (v & clearSizeMask) + clearSizeOverlay;
				o8 += v >> ((sizeof(T)-1)*8);
			}

			{
				T v = ((const T *)D)[(value>>(1*KK)) & overlappingMask];
				*((T *)o8) = (v & clearSizeMask) + clearSizeOverlay;
				o8 += v >> ((sizeof(T)-1)*8);
			}

			{
				T v = ((const T *)D)[(value>>(0*KK)) & overlappingMask];
				*((T *)o8) = (v & clearSizeMask) + clearSizeOverlay;
				o8 += v >> ((sizeof(T)-1)*8);
			}
		}
	}
	
	uint64_t valueBits = decompressor.O;
	while (i8 < src.end or valueBits>=decompressor.K+decompressor.O) {
		
		while (valueBits < decompressor.K+decompressor.O) {
			value = (value<<8) + uint64_t(*i8++);
			valueBits += 8;
		}
		
		size_t wordIdx = (value >> (valueBits-(decompressor.K+decompressor.O))) & overlappingMask;
		
		valueBits -= decompressor.K;
					
		{
			T v = ((const T *)D)[wordIdx];
//			*((T *)o8) = (v & clearSizeMask) + clearSizeOverlay;
//			o8 += v >> ((sizeof(T)-1)*8);
			
			
			size_t sz = v >> ((sizeof(T)-1)*8);
			
			T vv = (v & clearSizeMask) + clearSizeOverlay;
			memcpy(o8, &vv, std::min(sz,sizeof(T)-1));
//			*((T *)o8) = v;
			
			o8 += sz;
			
			
		}
	}
	// if (endMarlin-i8 != 0) std::cerr << " {" << endMarlin-i8 << "} "; // SOLVED! PROBLEM IN THE CODE
	// if (o8end-o8 != 0) std::cerr << " [" << o8end-o8 << "] "; // SOLVED! PROBLEM IN THE CODE

	return dst.nElements();
}

template<typename T, typename TSource, typename MarlinIdx>
size_t decompressFast(
	const TMarlinDecompress<TSource,MarlinIdx> &decompressor, 
	View<const uint8_t> src, View<TSource> dst) {
	
	auto K = decompressor.K;
	auto O = decompressor.O;
	
	if (K==8 and decompressor.isSkip) return decompress8_skip<T>(decompressor,src,dst);
	if (K==8) return decompressKK<T,8>(decompressor,src,dst);

	if (K==7) return decompressKK<T,7>(decompressor,src,dst);
	if (K==6) return decompressKK<T,6>(decompressor,src,dst);
	if (K==5) return decompressKK<T,5>(decompressor,src,dst);
	if (K==4) return decompressKK<T,4>(decompressor,src,dst);

	if (K==10) return decompressKK<T,10>(decompressor,src,dst);
	if (K==12) return decompressKK<T,12>(decompressor,src,dst);
	if (K==14) return decompressKK<T,14>(decompressor,src,dst);
	
	const uint8_t *i8    = src.start;
		  TSource *o8    = dst.start;

	const T clearSizeMask = T(-1)>>8;
	const T clearSizeOverlay = T(decompressor.marlinMostCommonSymbol) << ((sizeof(T)-1)*8);

	auto D = decompressor.decompressorTablePointer;

	uint64_t value = 0;
	uint64_t valueBits = O;
	while (i8 < src.end or valueBits>=K+O) {
		while (valueBits < K+O) {
			value += uint64_t(*i8++) << (64-8-valueBits);
			valueBits += 8;
		}
		
		size_t wordIdx = value >> (64-(K+O));
		
		value = value << K;
		valueBits -= K;

		{
			T v = ((const T *)D)[wordIdx];
//			*((T *)o8) = (v & clearSizeMask) + clearSizeOverlay;
//			o8 += v >> ((sizeof(T)-1)*8);

			size_t sz = v >> ((sizeof(T)-1)*8);
			
			T vv = (v & clearSizeMask) + clearSizeOverlay;
			memcpy(o8, &vv, std::min(sz,sizeof(T)-1));
//			*((T *)o8) = v;
			
			o8 += sz;
		}
	}
	
	return dst.nElements();
}

template<typename TSource, typename MarlinIdx>
size_t decompressSlow(
	const TMarlinDecompress<TSource,MarlinIdx> &decompressor, 
	View<const uint8_t> src, View<TSource> dst) {
	
	auto K = decompressor.K;
	auto O = decompressor.O;
	auto maxWordSize = decompressor.maxWordSize;
	
	const uint8_t *i8    = src.start;
		  TSource *o8    = dst.start;

	auto D = decompressor.decompressorTablePointer;
	
	size_t f = 1;
	while ((1<<f)<maxWordSize) f++;

	uint64_t value = 0;
	uint64_t valueBits = O;
	while (i8 < src.end or valueBits>=K+O) {
		while (valueBits < K+O) {
			value += uint64_t(*i8++) << (64-8-valueBits);
			valueBits += 8;
		}
		
		size_t wordIdx = value >> (64-(K+O));
		
		value = value << K;
		valueBits -= K;

		{
			memcpy(o8,&D[wordIdx<<f],maxWordSize); //fails for uint16_t
			o8 += D[(wordIdx<<f)+maxWordSize];
		}
	}
	
	return dst.nElements();
}
/*
template<typename TSource, typename MarlinIdx>
size_t decompressSlow(const TMarlinCompress<TSource,MarlinIdx> &decompressor, View<const uint8_t> src, View<TSource> dst) {
	
	const uint8_t *i8 = src.start;
		  TSource *o8 = dst.start;
	
	uint64_t value = 0;
	uint64_t valueBits = O;
	while (i8 < src.end or valueBits>=K+O) {
		while (valueBits < K+O) {
			value += uint64_t(*i8++) << (64-8-valueBits);
			valueBits += 8;
		}
		
		size_t wordIdx = value >> (64-(K+O));
		
		value = value << K;
		valueBits -= K;
					
		for (auto &&c : words[wordIdx])
			*o8++ = marlinAlphabet[c].sourceSymbol;
	}
	
	return dst.nElements();
}*/


}

template<typename TSource, typename MarlinIdx>
std::unique_ptr<std::vector<TSource>> TMarlinDecompress<TSource,MarlinIdx>::buildDecompressorTable(
	const TMarlinDictionary<TSource,MarlinIdx> &dictionary) const {
	
	auto &&marlinAlphabet = dictionary.marlinAlphabet;
	auto &&words = dictionary.words;
	
	auto ret = std::make_unique<std::vector<TSource>>(words.size()*(maxWordSize+1));
	
	TSource *d = &ret->front();
	for (size_t i=0; i<words.size(); i++) {
		for (size_t j=0; j<maxWordSize; j++)
			*d++ = (words[i].size()>j ? marlinAlphabet[words[i][j]].sourceSymbol : TSource(0));
		*d++ = words[i].size();
	}
	return ret;
}

template<typename TSource, typename MarlinIdx>
ssize_t TMarlinDecompress<TSource,MarlinIdx>::decompress(View<const uint8_t> src, View<TSource> dst) const {

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
	
	ssize_t nUnrepresentedSymbols = *src.start++;
	
	ssize_t unrepresentedSize = nUnrepresentedSymbols * ( sizeof(TSource) + (
		dst.nElements() < 0x100 ? sizeof(uint8_t) :
		dst.nElements() < 0x10000 ? sizeof(uint16_t) :
		dst.nElements() < 0x100000000ULL ? sizeof(uint32_t) :sizeof(uint64_t)
		));
	
	ssize_t residualSize = dst.nElements()*shift/8;

	ssize_t marlinSize = src.end-src.start-unrepresentedSize-residualSize;

	
	// Initialization, which might be optional
	if (not isSkip) {
		TSource s = marlinMostCommonSymbol;
		for (size_t i=0; i<dst.nElements(); i++)
			dst.start[i] = s;
	}

	View<const uint8_t> marlinSrc =
		marlin::make_view(src.start,src.start+marlinSize);
	View<const uint8_t> unrepresentedSrc = 
		marlin::make_view(marlinSrc.end,marlinSrc.end+unrepresentedSize);
	View<const uint8_t> shiftSrc  = 
		marlin::make_view(unrepresentedSrc.end,unrepresentedSrc.end+residualSize);

	if (false) {
		//decompressSlow(marlinSrc, dst);
	} else if (maxWordSize==3) {
		decompressFast<uint32_t>(*this,marlinSrc, dst);
	} else if (maxWordSize==7) {
		decompressFast<uint64_t>(*this,marlinSrc, dst);
	} else {
		//printf("Slow because: %lu %lu\n",K, maxWordSize);
		//decompressSlow(marlinSrc, dst);
	}
	
	//if (nUnrepresentedSymbols) printf("%u %u %u\n",  nUnrepresentedSymbols, unrepresentedSize, dst.nElements());
	// Place unrepresented symbols
	while (unrepresentedSrc.start < unrepresentedSrc.end) {

		size_t idx;
		if (dst.nElements() < 0x100) {
			idx = *reinterpret_cast<const uint8_t *&>(unrepresentedSrc.start)++;	
		} else if (dst.nElements() < 0x10000) {
			idx = *reinterpret_cast<const uint16_t *&>(unrepresentedSrc.start)++;	
		} else if (dst.nElements() < 0x100000000ULL) {
			idx = *reinterpret_cast<const uint32_t *&>(unrepresentedSrc.start)++;	
		} else {
			idx = *reinterpret_cast<const uint64_t *&>(unrepresentedSrc.start)++;	
		}
//			printf("%d \n", idx, *reinterpret_cast<const TSource *&>(unrepresentedSrc.start));
		dst.start[idx] = *reinterpret_cast<const TSource *&>(unrepresentedSrc.start)++;
//			printf("%llu %llu\n", unrepresentedSrc.start, unrepresentedSrc.end);
	}
			
	return padding + shift8(*this, shiftSrc, dst);
}

////////////////////////////////////////////////////////////////////////
//
// Explicit Instantiations
#include "instantiations.h"
INSTANTIATE(TMarlinDecompress)	

//INSTANTIATE_MEMBER(TMarlinDecompress,buildDecompressorTable(const TMarlinDictionary<TSource,MarlinIdx> &dictionary) const -> std::unique_ptr<std::vector<typename TMarlin::TSource_Type>>)	
//INSTANTIATE_MEMBER(TMarlinDecompress,decompress(View<const uint8_t> src, View<typename TMarlin::TSource_Type> dst) const -> ssize_t)	

