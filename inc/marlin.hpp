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
};
template<typename T> static View<T> make_view(T *start, T *end)  { return View<T>(start,end); }
template<typename T> static View<T> make_view(std::vector<T> &v) { return View<T>(&v[0], &v[v.size()]); }
template<typename T> static View<const T> make_view(const std::vector<T> &v) { return View<const T>(&v[0], &v[v.size()]); }
	
	
template<typename TSource, typename MarlinIdx>
struct TMarlin {
	
	typedef TSource TSource_Type;

	/// Name and Configuration map
	typedef std::map<std::string, double> Configuration;
	const std::string name; // Name of this dictionary (optional)
	const uint32_t versionMajor = MARLIN_VERSION_MAJOR;
	const uint32_t versionMinor = MARLIN_VERSION_MINOR;
//	enum { SimpleDictionary = 0 } type = Simple;
	const Configuration conf;


	/// Main Typedefs
	//typedef uint8_t  TSource; // Type that can source the source symbols.
	struct MarlinSymbol {
		TSource sourceSymbol;
		double p;
	};

	//typedef uint8_t  MarlinIdx; // Type that can source the marlin symbols.
	struct Word : std::vector<MarlinIdx> {
		using std::vector<MarlinIdx>::vector;
		double p = 0;
		MarlinIdx state = 0;
	};
	
	/// Convenience configuration
	const size_t K                = conf.at("K");           // Non overlapping bits of codeword.
	const size_t O                = conf.at("O");           // Bits that overlap between codewprds.
	const size_t shift            = conf.at("shift");       // Bits that can be stored raw
	const size_t maxWordSize      = conf.at("maxWordSize"); // Maximum number of symbols per word.

	/// ALPHABETS
	const std::vector<double> sourceAlphabet;	
	const std::vector<MarlinSymbol> marlinAlphabet = buildMarlinAlphabet();
		
	/// DICTIONARY
	//Marlin only encodes a subset of the possible source symbols.
	//Marlin symbols are sorted by probability in descending order, 
	//so the Marlin Symbol 0 is always corresponds to the most probable alphabet symbol.
	const std::vector<Word> words = buildDictionary(); // All dictionary words.
	const double efficiency       = calcEfficiency();  // Expected efficiency of the dictionary.
	
	/// DECOMPRESSOR STUFF
	const std::unique_ptr<std::vector<TSource>> decompressorTableVector = buildDecompressorTable();	
	const TSource* const decompressorTablePointer = decompressorTableVector->data();
	const TSource marlinMostCommonSymbol = marlinAlphabet.front().sourceSymbol;
	const bool isSkip = calcSkip();        // If all words are small, we can do a faster decoding algorithm;
	ssize_t decompress(View<const uint8_t> src, View<TSource> dst) const;
	ssize_t decompress(const std::vector<uint8_t> &src, std::vector<TSource> &dst) const {
		return decompress(make_view(src), make_view(dst));
	}
	
	/// COMPRESSOR STUFF
	typedef uint32_t CompressorTableIdx;      // Structured as: FLAG_NEXT_WORD Where to jump next	
	const std::unique_ptr<std::vector<CompressorTableIdx>> compressorTableVector = buildCompressorTable();	
	const CompressorTableIdx* const compressorTablePointer = compressorTableVector->data();	
	const std::unique_ptr<std::vector<CompressorTableIdx>> compressorTableInitVector = buildCompressorTableInit();
	const CompressorTableIdx* const compressorTableInitPointer = compressorTableInitVector->data();	
	ssize_t compress(View<const TSource> src, View<uint8_t> dst) const;
	ssize_t compress(const std::vector<TSource> &src, std::vector<uint8_t> &dst) const {
		ssize_t r = compress(make_view(src), make_view(dst));
		if (r<0) return r;
		dst.resize(r);
		return dst.size();
	}

	/// CONSTRUCTOR
	TMarlin( 
		std::string name_,
		const std::vector<double> &sourceAlphabet_,
		Configuration conf_ = Configuration()) 
		: 
		name(name_), 
		conf(updateConf(sourceAlphabet_, conf_)), 
		sourceAlphabet(sourceAlphabet_) 
		{}

	TMarlin( 
		std::string name_,
		Configuration conf_,
		double efficiency_,
		const TSource* const decompressorTablePointer_,
		TSource marlinMostCommonSymbol_,
		bool isSkip_,
		const CompressorTableIdx* const compressorTablePointer_,
		const CompressorTableIdx* const compressorTableInitPointer_		
	) : 
		name(name_), 
		conf(conf_),
		efficiency(efficiency_),
		decompressorTableVector(),
		decompressorTablePointer(decompressorTablePointer_),
		marlinMostCommonSymbol(marlinMostCommonSymbol_),
		isSkip(isSkip_),
		compressorTableVector(),
		compressorTablePointer(compressorTablePointer_),
		compressorTableInitVector(),
		compressorTableInitPointer(compressorTableInitPointer_)
		{}
		
private:
	// Sets default configurations
	static std::map<std::string, double> updateConf(const std::vector<double> &sourceAlphabet, Configuration conf);

	std::vector<MarlinSymbol> buildMarlinAlphabet() const;
	
	std::vector<Word> buildDictionary() const;
	double calcEfficiency() const;
	bool calcSkip() const;

	std::unique_ptr<std::vector<TSource>> buildDecompressorTable() const;
	std::unique_ptr<std::vector<CompressorTableIdx>> buildCompressorTable() const;
	std::unique_ptr<std::vector<CompressorTableIdx>> buildCompressorTableInit() const;
};

}

typedef marlin::TMarlin<uint8_t,uint8_t> Marlin;

#endif
#endif

