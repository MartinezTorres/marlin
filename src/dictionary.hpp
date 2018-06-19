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

#ifndef MARLIN_DICTIONARY_H
#define MARLIN_DICTIONARY_H

#include <marlin.h>

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <cstring>
#include <queue>
#include <stack>

#include <memory>
#include <algorithm>

#include <cassert>
 
	
class MarlinDictionary {
public:

	/// Configuration map
	typedef std::map<std::string, double> Configuration;
	const Configuration conf;


	/// Main Typedefs
	typedef uint8_t  SourceSymbol; // Type taht store the raw symbols from the source.
	struct MarlinSymbol {
		SourceSymbol sourceSymbol;
		double p;
	};
	struct Word : std::vector<SourceSymbol> {
		using std::vector<SourceSymbol>::vector;
		double p = 0;
		SourceSymbol state = 0;
	};
	
	/// BASIC CONFIGURATION
	const size_t K                = conf.at("K");           // Non overlapping bits of codeword.
	const size_t O                = conf.at("O");           // Bits that overlap between codewprds.
	const size_t shift            = conf.at("shift");       // Bits that can be stored raw
	const size_t maxWordSize      = conf.at("maxWordSize"); // Maximum number of symbols per word.

	/// ALPHABETS
	const std::map<SourceSymbol, double> sourceAlphabet;	
	const std::vector<MarlinSymbol> marlinAlphabet = buildMarlinAlphabet();
		
	/// DICTIONARY
	//Marlin only encodes a subset of the possible source symbols.
	//Marlin symbols are sorted by probability in descending order, 
	//so the Marlin Symbol 0 is always corresponds to the most probable alphabet symbol.
	const std::vector<Word> words = buildDictionary(); // All dictionary words.
	const double efficiency       = calcEfficiency();  // Theoretical efficiency of the dictionary.
	
	/// DECOMPRESSOR STUFF
	const std::shared_ptr<std::vector<SourceSymbol>> decompressorTableVector = buildDecompressorTable();	
	const SourceSymbol* const decompressorTablePointer = decompressorTableVector->data();
//	const SourceSymbol mostCommonSourceSymbol = alphabet.marlinSymbols.front().sourceSymbol;
	ssize_t decompress(uint8_t* dst, size_t dstSize, const uint8_t* src, size_t srcSize) const;
	
	/// COMPRESSOR STUFF
	typedef uint32_t CompressorTableIdx;      // Structured as: FLAG_NEXT_WORD Where to jump next	
	const std::shared_ptr<std::vector<CompressorTableIdx>> compressorTableVector = buildCompressorTable();	
	const CompressorTableIdx* const compressorTablePointer = compressorTableVector->data();	

	ssize_t compress(uint8_t* dst, size_t dstCapacity, const uint8_t* src, size_t srcSize) const;


	/// CONSTRUCTORS
	MarlinDictionary( const std::map<SourceSymbol, double> &sourceAlphabet_,
		Configuration conf_ = Configuration()) 
		: conf(updateConf(sourceAlphabet_, conf_)), sourceAlphabet(sourceAlphabet_) {}
	
	// Turns the vector into a map and uses the previous constructor
	MarlinDictionary( const std::vector<double> &sourceAlphabet_,
		Configuration conf_ = Configuration()) 
		: MarlinDictionary(
		[&sourceAlphabet_](){
			std::map<SourceSymbol, double> ret;
			for (size_t i=0; i<sourceAlphabet_.size(); i++)
				ret.emplace(SourceSymbol(i), sourceAlphabet_[i]);
			return ret;
		}()
		, conf_) {}
		
private:
	// Sets default configurations
	static std::map<std::string, double> updateConf(const std::map<SourceSymbol, double> &sourceAlphabet, Configuration conf);

	std::vector<MarlinSymbol> buildMarlinAlphabet() const;
	
	std::vector<Word> buildDictionary() const;
	double calcEfficiency() const;

	std::shared_ptr<std::vector<SourceSymbol>> buildDecompressorTable() const;
	std::shared_ptr<std::vector<CompressorTableIdx>> buildCompressorTable() const;
};



#endif //MARLIN_DICTIONARY_H
