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


#ifndef MARLIN_HPP
#define MARLIN_HPP
#if defined (__cplusplus)

#include "marlin.h"

#include <iostream>
#include <vector>
#include <map>
#include <memory>

// TSource is a type that represents the data type of the source 
//  (i.e., uint8_t and uint16_t are supported now)

// MarlinIdx is a type that can hold all symbols that will be encoded by Marlin
//  (which usually will be less than the Source symbols). 
//  uint8_t should be enough for all practical cases.
namespace marlin {

template<typename T>
struct View {
	T *start, *end;
	View(T *start_, T *end_) : start(start_), end(end_) {}
	size_t nElements() const { return end - start; }
	size_t nBytes() const { return sizeof(T)*(end - start); }
	T &operator[](size_t v) { return start[v]; }
	const T &operator[](size_t v) const { return start[v]; }
};
template<typename T> static View<T> make_view(T *start, T *end)  { return View<T>(start,end); }
template<typename T> static View<T> make_view(std::vector<T> &v) { return View<T>(&v[0], &v[v.size()]); }
template<typename T> static View<const T> make_view(const std::vector<T> &v) { return View<const T>(&v[0], &v[v.size()]); }


template<typename TSource>
struct MarlinSymbol_ {
	TSource sourceSymbol;
	double p;
};

template<typename MarlinIdx>
struct Word_ : std::vector<MarlinIdx> {
	using std::vector<MarlinIdx>::vector;
	double p = 0;
	MarlinIdx state = 0;
};

typedef std::map<std::string, double> Configuration;
	
template<typename TSource, typename MarlinIdx>
struct TMarlinDictionary{
	
	typedef Word_<MarlinIdx>  Word;
	typedef MarlinSymbol_<TSource> MarlinSymbol;

	const Configuration conf;

	const size_t K                = conf.at("K");           // Non overlapping bits of codeword.
	const size_t O                = conf.at("O");           // Bits that overlap between codewprds.
	const size_t shift            = conf.at("shift");       // Bits that can be stored raw
	const size_t maxWordSize      = conf.at("maxWordSize"); // Maximum number of symbols per word.
	
	struct MarlinAlphabet : public std::vector<MarlinSymbol> {
		using std::vector<MarlinSymbol>::vector;
		double probabilityOfUnrepresentedSymbol;
	};
	
	/// ALPHABETS
	const std::vector<double> sourceAlphabet;
	const double sourceEntropy = calcSourceEntropy(sourceAlphabet);
	const MarlinAlphabet marlinAlphabet = buildMarlinAlphabet();
		
	/// DICTIONARY
	//Marlin only encodes a subset of the possible source symbols.
	//Marlin symbols are sorted by probability in descending order, 
	//so the Marlin Symbol 0 is always corresponds to the most probable alphabet symbol.
	const std::vector<Word> words = buildDictionary(); // All dictionary words.
	const double efficiency       = calcEfficiency();  // Expected efficiency of the dictionary.
	const double compressionRatio = efficiency/sourceEntropy; //Expected compression ratio
	const bool isSkip = calcSkip();        // If all words are small, we can do a faster decoding algorithm;

	/// CONSTRUCTOR
	TMarlinDictionary( 
		const std::vector<double> &sourceAlphabet_,
		Configuration conf_ = Configuration()) 
		: 
		conf(updateConf(sourceAlphabet_, conf_)), 
		sourceAlphabet(sourceAlphabet_)
		{}

private:
	// Sets default configurations
	static std::map<std::string, double> updateConf(const std::vector<double> &sourceAlphabet, Configuration conf);

	MarlinAlphabet buildMarlinAlphabet() const;
	
	std::vector<Word> buildDictionary() const;
	static double calcSourceEntropy(const std::vector<double> &sourceAlphabet);
	double calcEfficiency() const;
	bool calcSkip() const;
};


template<typename TSource, typename MarlinIdx>
struct TMarlinCompress {
	
	typedef Word_<MarlinIdx>  Word;
	typedef MarlinSymbol_<TSource> MarlinSymbol;
	const size_t K,O,shift,maxWordSize;
	double efficiency;
	
	// Structured as: FLAG_NEXT_WORD Where to jump next	
	typedef uint32_t CompressorTableIdx;      
	const MarlinIdx unrepresentedSymbolToken;
	const std::array<MarlinIdx, 1U<<(sizeof(TSource)*8)> source2marlin;
	const std::shared_ptr<std::vector<CompressorTableIdx>> compressorTableVector;	
	const CompressorTableIdx* const compressorTablePointer;	
	const std::shared_ptr<std::vector<CompressorTableIdx>> compressorTableInitVector;
	const CompressorTableIdx* const compressorTableInitPointer;	

	ssize_t compress(View<const TSource> src, View<uint8_t> dst) const;
	ssize_t compress(const std::vector<TSource> &src, std::vector<uint8_t> &dst) const {
		ssize_t r = compress(make_view(src), make_view(dst));
		if (r<0) return r;
		dst.resize(r);
		return dst.size();
	}
	
	TMarlinCompress(const TMarlinDictionary<TSource,MarlinIdx> &dictionary) :
		K(dictionary.K), O(dictionary.O), shift(dictionary.shift), maxWordSize(dictionary.maxWordSize), 
		efficiency(dictionary.efficiency),
		unrepresentedSymbolToken(dictionary.marlinAlphabet.size()),
		source2marlin(buildSource2marlin(dictionary)),
		compressorTableVector(buildCompressorTable(dictionary)),
		compressorTablePointer(compressorTableVector->data()),
		compressorTableInitVector(buildCompressorTableInit(dictionary)),
		compressorTableInitPointer(compressorTableInitVector->data())
	{}

	TMarlinCompress(
		size_t K_, size_t O_, size_t shift_, size_t maxWordSize_, double efficiency_,
		MarlinIdx unrepresentedSymbolToken_,
		const std::array<MarlinIdx, 1U<<(sizeof(TSource)*8)> &source2marlin_,
		const CompressorTableIdx* const compressorTablePointer_,
		const CompressorTableIdx* const compressorTableInitPointer_
		) :
		K(K_), O(O_), shift(shift_), maxWordSize(maxWordSize_), efficiency(efficiency_),
		unrepresentedSymbolToken(unrepresentedSymbolToken_),
		source2marlin(source2marlin_),
		compressorTableVector(),
		compressorTablePointer(compressorTablePointer_),
		compressorTableInitVector(),
		compressorTableInitPointer(compressorTableInitPointer_)
	{}

	constexpr static const size_t FLAG_NEXT_WORD = 1UL<<(8*sizeof(CompressorTableIdx)-1);

private:
	std::array<MarlinIdx, 1U<<(sizeof(TSource)*8)> buildSource2marlin(const TMarlinDictionary<TSource,MarlinIdx> &dictionary) const;
	std::unique_ptr<std::vector<CompressorTableIdx>> buildCompressorTable(const TMarlinDictionary<TSource,MarlinIdx> &dictionary) const;
	std::unique_ptr<std::vector<CompressorTableIdx>> buildCompressorTableInit(const TMarlinDictionary<TSource,MarlinIdx> &dictionary) const;
};

template<typename TSource, typename MarlinIdx>
struct TMarlinDecompress {

	typedef Word_<MarlinIdx>  Word;
	typedef MarlinSymbol_<TSource> MarlinSymbol;
	const size_t K,O,shift,maxWordSize;
	
	const std::unique_ptr<std::vector<TSource>> decompressorTableVector;	
	const TSource* const decompressorTablePointer;
	const TSource marlinMostCommonSymbol;
	const bool isSkip;
	
	ssize_t decompress(View<const uint8_t> src, View<TSource> dst) const;
	ssize_t decompress(const std::vector<uint8_t> &src, std::vector<TSource> &dst) const {
		return decompress(make_view(src), make_view(dst));
	}

	TMarlinDecompress(const TMarlinDictionary<TSource,MarlinIdx> &dictionary) :
		K(dictionary.K), O(dictionary.O), shift(dictionary.shift), maxWordSize(dictionary.maxWordSize),
		decompressorTableVector(buildDecompressorTable(dictionary)),
		decompressorTablePointer(decompressorTableVector->data()),
		marlinMostCommonSymbol(dictionary.marlinAlphabet.front().sourceSymbol),
		isSkip(dictionary.isSkip)
	{}
	
	TMarlinDecompress(
		size_t K_, size_t O_, size_t shift_, size_t maxWordSize_,
		const TSource* const decompressorTablePointer_,
		const TSource marlinMostCommonSymbol_,
		const bool isSkip_
		) :
		K(K_), O(O_), shift(shift_), maxWordSize(maxWordSize_),
		decompressorTableVector(),
		decompressorTablePointer(decompressorTablePointer_),
		marlinMostCommonSymbol(marlinMostCommonSymbol_),
		isSkip(isSkip_)
	{}	
private:
	std::unique_ptr<std::vector<TSource>> buildDecompressorTable(const TMarlinDictionary<TSource,MarlinIdx> &dictionary) const;
};

template<typename TSource, typename MarlinIdx>
struct TMarlin : 
	public TMarlinCompress<TSource,MarlinIdx>, 
	public TMarlinDecompress<TSource,MarlinIdx> {
		
	const std::string name;
	const size_t K,O,shift,maxWordSize;

	
	TMarlin( 
		std::string name_,
		const std::vector<double> &sourceAlphabet_,
		Configuration conf_ = Configuration() ) :
		TMarlin( name_, TMarlinDictionary<TSource,MarlinIdx>(sourceAlphabet_, conf_) ) {}
	
	
	TMarlin( 		
		std::string name_,
		TMarlinDictionary<TSource,MarlinIdx> dictionary ) :
		TMarlinCompress<TSource,MarlinIdx>(dictionary),
		TMarlinDecompress<TSource,MarlinIdx>(dictionary),
		name(name_),
		K(dictionary.K), O(dictionary.O), shift(dictionary.shift), maxWordSize(dictionary.maxWordSize)
		 {}

		
	TMarlin(
		std::string name_,
		size_t K_, size_t O_, size_t shift_, size_t maxWordSize_, double efficiency_,
		MarlinIdx unrepresentedSymbolToken_,
		const std::array<MarlinIdx, 1U<<(sizeof(TSource)*8)> source2marlin_,
		const typename TMarlinCompress<TSource,MarlinIdx>::CompressorTableIdx* compressorTablePointer_,
		const typename TMarlinCompress<TSource,MarlinIdx>::CompressorTableIdx* compressorTableInitPointer_,
		const TSource* decompressorTablePointer_,
		const TSource marlinMostCommonSymbol_,
		const bool isSkip_
		) :
		TMarlinCompress<TSource,MarlinIdx>(
			K_, O_, shift_, maxWordSize_, efficiency_,
			unrepresentedSymbolToken_, source2marlin_, compressorTablePointer_, compressorTableInitPointer_),
		TMarlinDecompress<TSource,MarlinIdx>(
			K_, O_, shift_, maxWordSize_, decompressorTablePointer_, marlinMostCommonSymbol_, isSkip_),
		name(name_),
		K(K_), O(O_), shift(shift_), maxWordSize(maxWordSize_) {}
};
}


typedef marlin::TMarlin<uint8_t,uint8_t> Marlin;

#endif
#endif

