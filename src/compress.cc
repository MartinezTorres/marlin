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

namespace {

	constexpr static const size_t FLAG_NEXT_WORD = 1UL<<(8*sizeof(MarlinDictionary::CompressorTableIdx)-1);

	__attribute__ ((target ("bmi2")))
	void shift8(uint8_t shift, uint8_t* dst, const uint8_t* src, const size_t srcSize) {
		
		uint64_t mask=0;
		for (size_t i=0; i<8; i++)
			mask |= ((1ULL<<shift)-1)<<(8ULL*i);

		const uint64_t *i64    = (const uint64_t *)src;
		const uint64_t *i64end = (const uint64_t *)(src+srcSize);

		while (i64 != i64end) {
			*(uint64_t *)dst = _pext_u64(*i64++, mask);
			dst += shift;
		}
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
			table = T(((1<<wordStride)+unalignment)*(1<<alphaStride),MarlinDictionary::CompressorTableIdx(-1));
		}
		
		template<typename T, typename T0, typename T1>
		T &operator()(T *table, const T0 &word, const T1 &nextLetter) const { 
			return table[(word&((1<<wordStride)-1))+(nextLetter*((1<<wordStride)+unalignment))];
		}
	};
}


std::unique_ptr<std::vector<MarlinDictionary::CompressorTableIdx>> MarlinDictionary::buildCompressorTable() const {
	
	using MarlinIdx = SourceSymbol;

	MarlinIdx unrepresentedSymbolToken = marlinAlphabet.size();

	auto ret = std::make_unique<std::vector<CompressorTableIdx>>();
	JumpTable jump(K, O, unrepresentedSymbolToken);
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
				SourceSymbol lastSymbol = word.back();						
				word.pop_back();
				if (not positions[k].count(word)) throw(std::runtime_error("This word has no parent. SHOULD NEVER HAPPEN!!!"));
				size_t parentIdx = positions[k][word];
				jump(&ret->front(), parentIdx, source2marlin[lastSymbol]) = wordIdx;
				wordIdx = parentIdx;
			}
		}
	}
				
	//Link between inner dictionaries
	for (size_t k=0; k<NumChapters; k++)
		for (size_t i=k*ChapterSize; i<(k+1)*ChapterSize; i++)
			for (size_t j=0; j<marlinAlphabet.size(); j++)
				if (jump(&ret->front(),i,j)==CompressorTableIdx(-1)) // words that are not parent of anyone else.
					jump(&ret->front(),i,j) = positions[i%NumChapters][Word(1,marlinAlphabet[j].sourceSymbol)] + FLAG_NEXT_WORD;
						
	return ret;
}

ssize_t MarlinDictionary::compress(uint8_t* dst, size_t dstCapacity, const uint8_t* src, size_t srcSize) const {

	// Assertions
	assert(dstCapacity >= srcSize);
	
	// Special case: empty! Nothing to compress.
	if (srcSize==0) return 0;

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

	using MarlinIdx = SourceSymbol;

	MarlinIdx unrepresentedSymbolToken = marlinAlphabet.size();
	std::array<MarlinIdx, 1U<<(sizeof(MarlinIdx)*8)> source2marlin;
	source2marlin.fill(unrepresentedSymbolToken);
	for (size_t i=0; i<marlinAlphabet.size(); i++)
		source2marlin[marlinAlphabet[i].sourceSymbol>>shift] = i;
		
	JumpTable jump(K, O, unrepresentedSymbolToken);

	// Encode Marlin, with rare symbols preceded by an empty word
	{
		
		// if the encoder produces a size larger than this, it is simply better to store the block uncompressed.
		size_t maxTargetSize = std::max(size_t(8), srcSize-srcSize*shift/8) - 8;

		CompressorTableIdx j = 0; 
		while (i8<i8end) {				
			
			SourceSymbol ss = *i8++;
			
			MarlinIdx ms = source2marlin[ss>>shift];
			bool isUnrepresented = ms==unrepresentedSymbolToken;
			if (isUnrepresented) {
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
			
			CompressorTableIdx jOld = j;
			j = jump(decompressorTablePointer, j, ms);
			
			if (j & FLAG_NEXT_WORD) 
				*o8++ = jOld & 0xFF;
		}
		if (j) *o8++ = j;
	}

	// Encode residuals
	shift8(shift, o8, src, srcSize);
	
	return o8 - dst; 
}
