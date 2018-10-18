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

#define   LIKELY(condition) (__builtin_expect(static_cast<bool>(condition), 1))
#define UNLIKELY(condition) (__builtin_expect(static_cast<bool>(condition), 0))

using namespace marlin;

namespace {


template<typename TSource, typename MarlinIdx>
__attribute__ ((target ("bmi2")))
ssize_t shift8(const TMarlinCompress<TSource,MarlinIdx> &compressor, View<const TSource> src, View<uint8_t> dst) {
	
	uint64_t mask=0;
	for (size_t i=0; i<8; i++)
		mask |= ((1ULL<<compressor.shift)-1)<<(8ULL*i);

	const uint64_t *i64    = reinterpret_cast<const uint64_t *>(src.start);
	const uint64_t *i64end = reinterpret_cast<const uint64_t *>(src.end);

	uint8_t *o8 = dst.start;
	
	while (i64 != i64end) {
		*reinterpret_cast<uint64_t *>(o8) = _pext_u64(*i64++, mask);
		o8 += compressor.shift;
	}
	return o8 - dst.start;
}



class JumpTable {

	const size_t alphaStride;  // Bit stride of the jump table corresponding to the word dimension
	const size_t wordStride;  // Bit stride of the jump table corresponding to the word dimension
public:
	
	JumpTable(size_t keySize, size_t overlap, size_t nAlpha) :
		alphaStride(std::ceil(std::log2(nAlpha))),
		wordStride(keySize+overlap) {}
		
	template<typename T>
	void initTable(std::vector<T> &table) {
		table = std::vector<T>(((1<<wordStride))*(1<<alphaStride),T(-1));
	}
	
	template<typename T, typename T0, typename T1>
	inline T &operator()(T *table, const T0 &word, const T1 &nextLetter) const { 
		auto v = (word&((1<<wordStride)-1))+(nextLetter<<wordStride);
//			auto v = ((word&((1<<wordStride)-1))<<alphaStride)+nextLetter;
		return table[v];
	}
};


template<typename TSource, typename MarlinIdx>
ssize_t compressMarlin8 (
	const TMarlinCompress<TSource,MarlinIdx> &compressor,
	View<const TSource> src, 
	View<uint8_t> dst, 
	std::vector<size_t> &unrepresentedSymbols)
{
	
	JumpTable jump(compressor.K, compressor.O, compressor.unrepresentedSymbolToken+1);	
	
		  uint8_t *out   = dst.start;
	const TSource *in    = src.start;

	typename TMarlinCompress<TSource,MarlinIdx>::CompressorTableIdx j = 0; 

	
	//We look for the word that sets up the machine state.
	{		
		TSource ss = *in++;
		
		MarlinIdx ms = compressor.source2marlin[ss>>compressor.shift];
		if (ms==compressor.unrepresentedSymbolToken) {
			unrepresentedSymbols.push_back(in-src.start-1);
			ms = 0; // 0 must be always the most probable symbol;
			//printf("%04x %02x\n", in-src.start-1, ss);
		}
		
		j = compressor.compressorTableInitPointer[ms];
	}
	
	while (in<src.end) {
		
		MarlinIdx ms = compressor.source2marlin[(*in++)>>compressor.shift];
		if (ms==compressor.unrepresentedSymbolToken) {
			unrepresentedSymbols.push_back(in-src.start-1);
			ms = 0; // 0 must be always the most probable symbol;
			//printf("%04x %02x\n", in-src.start-1, ss);
		}
		
		*out = j & 0xFF;
		j = jump(compressor.compressorTablePointer, j, ms);
		
		if (j & compressor.FLAG_NEXT_WORD) {
			out++;
		}

		if (dst.end-out<16) return -1;	// TODO: find the exact value
	}


	//if (not (j & FLAG_NEXT_WORD)) 
	*out++ = j & 0xFF;
	
	return out - dst.start;
}

template<typename TSource, typename MarlinIdx>
ssize_t compressMarlinFast(
	const TMarlinCompress<TSource,MarlinIdx> &compressor,
	View<const TSource> src, 
	View<uint8_t> dst, 
	std::vector<size_t> &unrepresentedSymbols) 
{
	
	JumpTable jump(compressor.K, compressor.O, compressor.unrepresentedSymbolToken+1);	
	
		  uint8_t *out   = dst.start;
	const TSource *in    = src.start;

	typename TMarlinCompress<TSource,MarlinIdx>::CompressorTableIdx j = 0; 

	
	//We look for the word that sets up the machine state.
	{		
		TSource ss = *in++;
		
		MarlinIdx ms = compressor.source2marlin[ss>>compressor.shift];
		if (ms==compressor.unrepresentedSymbolToken) {
			unrepresentedSymbols.push_back(in-src.start-1);
			ms = 0; // 0 must be always the most probable symbol;
			//printf("%04x %02x\n", in-src.start-1, ss);
		}
		
		j = compressor.compressorTableInitPointer[ms];
	}

	uint32_t value = 0;
	 int32_t valueBits = 0;
	
	while (in<src.end) {
		
		if (dst.end-out<16) return -1;	// TODO: find the exact value
		
		TSource ss = *in++;
		
		MarlinIdx ms = compressor.source2marlin[ss>>compressor.shift];
		if (ms==compressor.unrepresentedSymbolToken) {
			unrepresentedSymbols.push_back(in-src.start-1);
			ms = 0; // 0 must be always the most probable symbol;
			//printf("%04x %02x\n", in-src.start-1, ss);
		}
		
		auto jOld = j;
		j = jump(compressor.compressorTablePointer, j, ms);
		
		if (j & compressor.FLAG_NEXT_WORD) {
			
			value |= ((jOld | compressor.FLAG_NEXT_WORD) ^ compressor.FLAG_NEXT_WORD) << (32 - compressor.K - valueBits);
			valueBits += compressor.K;
		}
		
		while (valueBits>8) {
			*out++ = value >> 24;
			value = value << 8;
			valueBits -= 8;
		}
	}

	value |= ((j | compressor.FLAG_NEXT_WORD) ^ compressor.FLAG_NEXT_WORD) << (32 - compressor.K - valueBits);
	valueBits += compressor.K;
	
	while (valueBits>0) {
		*out++ = value >> 24;
		value = value << 8;
		valueBits -= 8;
	}
	
	return out - dst.start;
}	

}

template<typename TSource, typename MarlinIdx>
auto TMarlinCompress<TSource,MarlinIdx>::buildCompressorTableInit(const TMarlinDictionary<TSource,MarlinIdx> &dictionary) const -> std::unique_ptr<std::vector<CompressorTableIdx>> {

	auto ret = std::make_unique<std::vector<CompressorTableIdx>>();
		
	for (size_t ms=0; ms<dictionary.marlinAlphabet.size(); ms++) {		

		for (size_t i=0; i<size_t(1<<K); i++) {
			if (dictionary.words[i].size()==1 and dictionary.words[i].front() == ms) {
				ret->push_back(i);
				break;
			}
		}
	}
	return ret;
}

template<typename TSource, typename MarlinIdx>
auto TMarlinCompress<TSource,MarlinIdx>::buildCompressorTable(const TMarlinDictionary<TSource,MarlinIdx> &dictionary) const -> std::unique_ptr<std::vector<CompressorTableIdx>> {

	auto ret = std::make_unique<std::vector<CompressorTableIdx>>();
	JumpTable jump(K, O, unrepresentedSymbolToken+1);
	jump.initTable(*ret);
	
	const size_t NumChapters = 1<<O;
	const size_t ChapterSize = 1<<K;
	std::vector<std::map<Word, size_t>> positions(NumChapters);

	// Init the mapping (to know where each word goes)
	for (size_t k=0; k<NumChapters; k++)
		for (size_t i=k*ChapterSize; i<(k+1)*ChapterSize; i++)
			positions[k][dictionary.words[i]] = i;
			
	// Link each possible word to its continuation
	for (size_t k=0; k<NumChapters; k++) {
		for (size_t i=k*ChapterSize; i<(k+1)*ChapterSize; i++) {
			auto word = dictionary.words[i];
			size_t wordIdx = i;
			while (word.size() > 1) {
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
			for (size_t j=0; j<dictionary.marlinAlphabet.size(); j++)
				if (jump(&ret->front(),i,j) == CompressorTableIdx(-1)) // words that are not parent of anyone else.
					jump(&ret->front(),i,j) = positions[i%NumChapters][Word(1,j)] + FLAG_NEXT_WORD;
										
	return ret;
}

template<typename TSource, typename MarlinIdx>
ssize_t TMarlinCompress<TSource,MarlinIdx>::compress(View<const TSource> src, View<uint8_t> dst) const {

	
	//memcpy(dst.start,src.start,src.nBytes()); return src.nBytes();

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

	size_t residualSize = src.nElements()*shift/8;


	std::vector<size_t> unrepresentedSymbols;		
	// This part, we encode the number of unrepresented symbols in a byte.
	// We are optimistic and we hope that no unrepresented symbols are required.
	*dst.start = 0;
	
	// Valid portion available to encode the marlin message.
	View<uint8_t> marlinDst = marlin::make_view(dst.start+1,dst.end-residualSize);
	ssize_t marlinSize = -1;
	if (false) {
		//marlinSize = compressMarlinReference(src, marlinDst, unrepresentedSymbols);
	} else if (K==8) {
		marlinSize = compressMarlin8(*this, src, marlinDst, unrepresentedSymbols);
	} else {
		marlinSize = compressMarlinFast(*this, src, marlinDst, unrepresentedSymbols);
	}
	
	size_t unrepresentedSize = unrepresentedSymbols.size() * ( sizeof(TSource) + (
		src.nElements() < 0x100 ? sizeof(uint8_t) :
		src.nElements() < 0x10000 ? sizeof(uint16_t) :
		src.nElements() < 0x100000000ULL ? sizeof(uint32_t) :sizeof(uint64_t)
		));
	
	
	//if (unrepresentedSize) printf("%d \n", unrepresentedSize);
	// If not worth encoding, we store raw.
	if (marlinSize < 0 	// If the encoded size is negative means that Marlin could not provide any meaningful compression, and the whole stream will be copied.
		or unrepresentedSymbols.size() > 255 
		or 1 + marlinSize + unrepresentedSize + residualSize > src.nBytes()) {

		memcpy(dst.start,src.start,src.nBytes());
		return padding + src.nBytes();
	}
	
	
	*dst.start++ = unrepresentedSymbols.size();
	dst.start += marlinSize;
	
	
	// Encode unrepresented symbols
	for (auto &s : unrepresentedSymbols) {
		if (src.nElements() < 0x100) {
			*reinterpret_cast<uint8_t *&>(dst.start)++ = s;	
		} else if (src.nElements() < 0x10000) {
			*reinterpret_cast<uint16_t *&>(dst.start)++ = s;	
		} else if (src.nElements() < 0x100000000ULL) {
			*reinterpret_cast<uint32_t *&>(dst.start)++ = s;	
		} else {
			*reinterpret_cast<uint64_t *&>(dst.start)++ = s;	
		}
		*reinterpret_cast<TSource *&>(dst.start)++ = src.start[s];	
	}
	
	// Encode residuals
	shift8(*this, src, dst);
	
	return padding + 1 + marlinSize + unrepresentedSize + residualSize; 
}

template<typename TSource, typename MarlinIdx>
std::array<MarlinIdx, 1U<<(sizeof(TSource)*8)> TMarlinCompress<TSource,MarlinIdx>::buildSource2marlin(
	const TMarlinDictionary<TSource,MarlinIdx> &dictionary) const {

	std::array<MarlinIdx, 1U<<(sizeof(TSource)*8)> source2marlin_;
	source2marlin_.fill(unrepresentedSymbolToken);
	for (size_t i=0; i<dictionary.marlinAlphabet.size(); i++)
		source2marlin_[dictionary.marlinAlphabet[i].sourceSymbol>>shift] = i;
	return source2marlin_;
}


////////////////////////////////////////////////////////////////////////
//
// Explicit Instantiations
#include "instantiations.h"
INSTANTIATE(TMarlinCompress)	

//INSTANTIATE_MEMBER(buildCompressorTableInit() const -> std::unique_ptr<std::vector<CompressorTableIdx>>)	
//INSTANTIATE_MEMBER(buildCompressorTable() const -> std::unique_ptr<std::vector<CompressorTableIdx>>)	
//INSTANTIATE_MEMBER(compress(View<const typename TMarlin::TSource_Type> src, View<uint8_t> dst) const -> ssize_t)	

