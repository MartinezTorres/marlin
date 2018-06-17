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
 
	
struct MarlinDictionary {

	// Main Typedefs
	typedef uint8_t  SourceSymbol; // Type taht store the raw symbols from the source.
	typedef uint8_t  MarlinSymbol; 
	
	// Configuration map
	typedef std::map<std::string, double> Configuration;
	const Configuration conf;

	static std::map<std::string, double> updateConf( // Sets default configurations
		const std::map<SourceSymbol, double> &symbols, 
		Configuration conf);


	// The Alphabet Class acts as a translation layer between SourceSymbols and MarlinSymbols.
	class Alphabet {

		static double calcEntropy(const std::map<SourceSymbol, double> &symbols) {
			
			double distEntropy=0;
			for (auto &&s : symbols)
				if (s.second>0.)
					distEntropy += -s.second*std::log2(s.second);
			return distEntropy;
		}
		
	public:
	
		const std::map<SourceSymbol, double> sourceSymbolProbability;
		const double sourceEntropy;
		
		double rareSymbolProbability;
		std::vector<double> marlinToSource;
		std::vector<SymbolAndProbability> marlinSymbols;

		
		Alphabet(std::map<SourceSymbol, double> sourceSymbols,const  Configuration &conf) : 
			sourceSymbols(sourceSymbols)
			sourceEntropy(calcEntropy(symbols)) {

			size_t shift = conf.at("shift"));

			
			// Group symbols by their high bits
			std::map<SourceSymbol, double> symbolsShifted;
			for (auto &&symbol : symbols)
				symbolsShifted[symbol.first>>shift] += symbol.second;
			
			for (auto &&symbol : symbolsShifted)
				marlinSymbols.push_back(SymbolAndProbability({SourceSymbol(symbol.first<<shift), symbol.second}));
				
			std::stable_sort(marlinSymbols.begin(),marlinSymbols.end());
			
			rareSymbolProbability = 0;
			while (marlinSymbols.size()>conf.at("minMarlinSymbols") and 
				  (marlinSymbols.size()>conf.at("maxMarlinSymbols") or
				  rareSymbolProbability<conf.at("purgeProbabilityThreshold"))) {
				
				rareSymbolProbability += marlinSymbols.back().p;
//					marlinSymbols.front().p +=  marlinSymbols.back().p;
				marlinSymbols.pop_back();
			}
		}
	};
				
	struct Word : std::vector<SourceSymbol> {
		using std::vector<SourceSymbol>::vector;
		double p = 0;
		SourceSymbol state = 0;

		friend std::ostream& operator<< (std::ostream& stream, const Word& word) {
			for (auto &&s : word) 
				if (s<=26) 
					stream << char('a'+s); 
				else 
					stream << " #" << uint(s);
			return stream;
		}
	};
		
	/// DICTIONARY
	//Marlin only encodes a subset of the possible source symbols.
	//Marlin symbols are sorted by probability in descending order, 
	//so the Marlin Symbol 0 is always corresponds to the most probable alphabet symbol.
	const Alphabet alphabet;
	const size_t K                = conf.at("K");           // Non overlapping bits of codeword.
	const size_t O                = conf.at("O");           // Bits that overlap between codewprds.
	const size_t shift            = conf.at("shift");       // Bits that can be stored raw
	const size_t maxWordSize      = conf.at("maxWordSize"); // Maximum number of symbols per word.
	const std::vector<Word> words = buildDictionary();      // All dictionary words.
	const double efficiency       = calcEfficiency(words);  // Theoretical efficiency of the dictionary.

	double calcEfficiency(std::vector<Word> dictionary) const;
	std::vector<Word> buildDictionary() const;

	// Decompressor
	const std::shared_ptr<std::vector<SourceSymbol>> decompressorTableVector = buildDecompressorTable();	
	const SourceSymbol* const decompressorTablePointer = decompressorTableVector->data();
	const SourceSymbol mostCommonSourceSymbol = alphabet.marlinSymbols.front().sourceSymbol;

	std::shared_ptr<std::vector<SourceSymbol>> buildDecompressorTable() const;
	ssize_t decompress(uint8_t* dst, size_t dstSize, const uint8_t* src, size_t srcSize) const;
	
	// Compressor
	typedef uint32_t WordIdx;      // Structured as: FLAG_NEXT_WORD Where to jump next	
	const std::shared_ptr<std::vector<WordIdx>> compressorTableVector = buildCompressorTable();	
	const WordIdx* const compressorTablePoointer = compressorTableVector->data();	

	std::shared_ptr<std::vector<WordIdx>> buildCompressorTable() const;
	ssize_t compress(uint8_t* dst, size_t dstCapacity, const uint8_t* src, size_t srcSize) const;


	MarlinDictionary( const std::map<SourceSymbol, double> &symbols,
		std::map<std::string, double> conf_ = std::map<std::string, double>()) 
		: conf(updateConf(symbols, conf_)), alphabet(symbols, conf) {}
	
	// Turns the vector into a map and uses the previous constructor
	MarlinDictionary( const std::vector<double> &symbols,
		std::map<std::string, double> conf_ = std::map<std::string, double>()) 
		: MarlinDictionary(
		[&symbols](){
			std::map<SourceSymbol, double> ret;
			for (size_t i=0; i<symbols.size(); i++)
				ret.emplace(SourceSymbol(i), symbols[i]);
			return ret;
		}()
		, conf_) {}
};



#endif //MARLIN_DICTIONARY_H
