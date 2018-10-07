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

using namespace marlin;

namespace {

template<typename TSource, typename MarlinIdx>
struct TCompress : TMarlin<TSource,MarlinIdx> {
	
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
	using TMarlin<TSource, MarlinIdx>::compressorTablePointer;
	


	constexpr static const size_t FLAG_NEXT_WORD = 1UL<<(8*sizeof(CompressorTableIdx)-1);

	__attribute__ ((target ("bmi2")))
	ssize_t shift8(View<const TSource> src, View<uint8_t> dst) const {
		
		uint64_t mask=0;
		for (size_t i=0; i<8; i++)
			mask |= ((1ULL<<shift)-1)<<(8ULL*i);

		const uint64_t *i64    = reinterpret_cast<const uint64_t *>(src.start);
		const uint64_t *i64end = reinterpret_cast<const uint64_t *>(src.end);

		uint8_t *o8 = dst.start;
		
		while (i64 != i64end) {
			*reinterpret_cast<uint64_t *>(o8) = _pext_u64(*i64++, mask);
			o8 += shift;
		}
		return o8 - dst.start;
	}
	
	class JumpTable {

		constexpr static const size_t unalignment = 8; // Too much aligned reads break cache
		const size_t alphaStride;  // Bit stride of the jump table corresponding to the word dimension
		const size_t wordStride;  // Bit stride of the jump table corresponding to the word dimension
	public:
		
		JumpTable(size_t keySize, size_t overlap, size_t nAlpha) :
			alphaStride(std::ceil(std::log2(nAlpha))),
			wordStride(keySize+overlap) {}
			
		template<typename T>
		void initTable(T &table) {
			table = T(((1<<wordStride)+unalignment)*(1<<alphaStride),CompressorTableIdx(-1));
		}
		
		template<typename T, typename T0, typename T1>
		T &operator()(T *table, const T0 &word, const T1 &nextLetter) const { 
			return table[(word&((1<<wordStride)-1))+(nextLetter*((1<<wordStride)+unalignment))];
		}
	};

	ssize_t compressMarlin8(View<const TSource> src, View<uint8_t> dst) const {
		
		MarlinIdx unrepresentedSymbolToken = marlinAlphabet.size();
		JumpTable jump(K, O, unrepresentedSymbolToken+1);	
		
		std::array<MarlinIdx, 1U<<(sizeof(TSource)*8)> source2marlin;
		source2marlin.fill(unrepresentedSymbolToken);
		for (size_t i=0; i<marlinAlphabet.size(); i++)
			source2marlin[marlinAlphabet[i].sourceSymbol>>shift] = i;

			  uint8_t *out   = dst.start;
		const TSource *in    = src.start;

		CompressorTableIdx j = 0; 
		while (in<src.end) {
			
			if (dst.end-out<16) return -1;	// TODO: find the exact value
			
			TSource ss = *in++;
			
			MarlinIdx ms = source2marlin[ss>>shift];
			bool isUnrepresented = ms==unrepresentedSymbolToken;
			if (isUnrepresented) {
				//printf("unrepresented\n");
				if (j) *out++ = j; // Finish current word, if any;
				*out++ = j = 0;
				*reinterpret_cast<TSource *&>(out)++ = ss;
				//*out++ = ss & ((1<<shift)-1);
				continue;
			}
			
			CompressorTableIdx jOld = j;
			j = jump(compressorTablePointer, j, ms);
			
			if (j & FLAG_NEXT_WORD) 
				*out++ = jOld & 0xFF;
		}
		if (j) *out++ = j;
		
		return out - dst.start;
	}
	
	ssize_t compressMarlinReference(View<const TSource> src, View<uint8_t> dst) const {
		
		MarlinIdx unrepresentedSymbolToken = marlinAlphabet.size();
		JumpTable jump(K, O, unrepresentedSymbolToken+1);
		
		std::array<MarlinIdx, 1U<<(sizeof(TSource)*8)> source2marlin;
		source2marlin.fill(unrepresentedSymbolToken);
		for (size_t i=0; i<marlinAlphabet.size(); i++)
			source2marlin[marlinAlphabet[i].sourceSymbol>>shift] = i;

			  uint8_t *out   = dst.start;
		const TSource *in    = src.start;


		std::map<Word, size_t> wordMap;
		for (size_t i=0; i<words.size(); i++) 
			wordMap[words[i]] = i;
		
		uint32_t value = 0;
		 int32_t valueBits = 0;
		Word word;
		while (in<src.end) {				
			
			TSource ss = *in++;
			MarlinIdx ms = source2marlin[ss>>shift];
			
			word.push_back(ms);
			
			if (wordMap.count(word) == 0) {
				
				word.pop_back();

				value += wordMap[word] << (32 - K - valueBits);
				valueBits += K;
				
				word.clear();
				word.push_back(ms);
			}

			while (valueBits>8) {
				*out++ = value >> 24;
				value = value << 8;
				valueBits -= 8;
			}		
			
		}
		if (word.size()) {
			value += wordMap[word] << (32 - K - valueBits);
			valueBits += K;
		}
		
		while (valueBits>0) {
			*out++ = value >> 24;
			value = value << 8;
			valueBits -= 8;
		}
		
		return out - dst.start;
	}

	std::unique_ptr<std::vector<CompressorTableIdx>> buildCompressorTable() const {

		MarlinIdx unrepresentedSymbolToken = marlinAlphabet.size();

		auto ret = std::make_unique<std::vector<CompressorTableIdx>>();
		JumpTable jump(K, O, unrepresentedSymbolToken+1);
		jump.initTable(*ret);
		
		std::array<MarlinIdx, 1U<<(sizeof(MarlinIdx)*8)> source2marlin;
		source2marlin.fill(unrepresentedSymbolToken);
		for (size_t i=0; i<marlinAlphabet.size(); i++)
			source2marlin[marlinAlphabet[i].sourceSymbol>>shift] = i;
		
		const size_t NumChapters = 1<<O;
		const size_t ChapterSize = 1<<K;
		std::vector<std::map<Word, size_t>> positions(NumChapters);

		// Init the mapping (to know where each word goes)
		for (size_t k=0; k<NumChapters; k++)
			for (size_t i=k*ChapterSize; i<(k+1)*ChapterSize; i++)
				positions[k][words[i]] = i;
				
		// Link each possible word to its continuation
		for (size_t k=0; k<NumChapters; k++) {
			for (size_t i=k*ChapterSize; i<(k+1)*ChapterSize; i++) {
				auto word = words[i];
				size_t wordIdx = i;
				while (not word.empty()) {
					TSource lastSymbol = word.back();						
					word.pop_back();
					if (not positions[k].count(word)) throw(std::runtime_error("This word has no parent. SHOULD NEVER HAPPEN!!!"));
					size_t parentIdx = positions[k][word];
					jump(&ret->front(), parentIdx, lastSymbol) = wordIdx;
					wordIdx = parentIdx;
				}
			}
		}
					
		//Link between inner dictionaries
		for (size_t k=0; k<NumChapters; k++)
			for (size_t i=k*ChapterSize; i<(k+1)*ChapterSize; i++)
				for (size_t j=0; j<marlinAlphabet.size(); j++)
					if (jump(&ret->front(),i,j) == CompressorTableIdx(-1)) // words that are not parent of anyone else.
						jump(&ret->front(),i,j) = positions[i%NumChapters][Word(1,j)] + FLAG_NEXT_WORD;
											
		return ret;
	}

	ssize_t compress(View<const TSource> src, View<uint8_t> dst) const {

		// Assertions
		if (dst.nBytes() < src.nBytes()) return -1; //TODO: Real error codes
		
		// Special case: empty! Nothing to compress.
		if (src.nElements()==0) return 0;

		// Special case: the entire block is made of one symbol!
		{
			size_t count = 0;
			for (size_t i=0; i<src.nElements() and src.start[i] == src.start[0]; i++)
				count++;
			
			if (count==src.nElements()) {
				reinterpret_cast<TSource *>(dst.start)[0] = src.start[0];
				return sizeof(TSource);
			}
		}

		// Special case: if srcSize is not multiple of 8, we force it to be.
		size_t padding = 0;
		while (src.nBytes() % 8 != 0) {
			
			*reinterpret_cast<TSource *&>(dst.start)++ = *src.start++;			
			padding += sizeof(TSource);
		}

		// Encode Marlin, with rare symbols preceded by an empty word
		View<uint8_t> marlinDst = marlin::make_view(dst.start,dst.end-(src.nElements()*shift/8));
		ssize_t marlinSize;
		if (K==8) {
			marlinSize = compressMarlin8(src, marlinDst);
		} else {
			marlinSize = compressMarlinReference(src, marlinDst);
		}
		
		// If the encoded size is negative means that Marlin could not provide any meaningful compression, and the whole stream will be copied.
		if (marlinSize < 0) {
			memcpy(dst.start,src.start,src.nBytes());
			return padding + src.nBytes();
		}
		
		// Encode residuals
		dst.start += marlinSize;
		size_t residualSize = shift8(src, dst);
		
		return padding + marlinSize + residualSize; 
	}

};
}

template<typename TSource, typename MarlinIdx>
auto TMarlin<TSource,MarlinIdx>::buildCompressorTable() const -> std::unique_ptr<std::vector<CompressorTableIdx>> {
	return reinterpret_cast<const TCompress<TSource,MarlinIdx> *>(this)->buildCompressorTable();
}

template<typename TSource, typename MarlinIdx>
ssize_t TMarlin<TSource,MarlinIdx>::compress(View<const TSource> src, View<uint8_t> dst) const {
	return reinterpret_cast<const TCompress<TSource,MarlinIdx> *>(this)->compress(src,dst);
}

////////////////////////////////////////////////////////////////////////
//
// Explicit Instantiations
#include "instantiations.h"
INSTANTIATE_MEMBER(buildCompressorTable() const -> std::unique_ptr<std::vector<CompressorTableIdx>>)	
INSTANTIATE_MEMBER(compress(View<const typename TMarlin::TSource_Type> src, View<uint8_t> dst) const -> ssize_t)	

