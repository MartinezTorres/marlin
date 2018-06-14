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
	
	
	
};

namespace {
	
struct Dictionary : public MarlinDictionary {

	// Configuration map
	typedef std::map<std::string, double> Configuration;
	const Configuration conf;

	// Main Typedefs
	typedef uint8_t SourceSymbol;
	typedef uint8_t MarlinSymbol;
	typedef uint32_t WordIdx; // Structured as: FLAG_NEXT_WORD Where to jump next	

	// The Alphabet Class acts as a translation layer between SourceSymbols and MarlinSymbols.
	class Alphabet {

		struct SymbolAndProbability {
			SourceSymbol sourceSymbol;
			double p;
			constexpr bool operator<(const SymbolAndProbability &rhs) const {
				if (p!=rhs.p) return p>rhs.p; // Descending in probability
				return sourceSymbol<rhs.sourceSymbol; // Ascending in symbol index
			}
		};
		
		static double calcEntropy(const std::map<SourceSymbol, double> &symbols) {
			
			double distEntropy=0;
			for (auto &&s : symbols)
				if (s.second>0.)
					distEntropy += -s.second*std::log2(s.second);
			return distEntropy;
		}
		
	public:
	
		const std::map<SourceSymbol, double> symbols;
		const size_t shift;
		const double sourceEntropy;
		
		double rareSymbolProbability;
		std::vector<SymbolAndProbability> marlinSymbols;
		
		Alphabet(std::map<SourceSymbol, double> symbols_, Configuration conf) : 
			symbols(symbols_),
			shift(conf.at("S")),
			sourceEntropy(calcEntropy(symbols)) {
			
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
		MarlinSymbol state = 0;

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


	// Decoder
	std::vector<SourceSymbol> decoderTableVector;
	const SourceSymbol * const decoderTablePointer;
//	const SourceSymbol mostCommonSourceSymbol;
	
	// Encoder
//	const size_t nMarlinSymbols;
	std::vector<WordIdx> encoderTableVector;	//table(((1<<wordStride)+unalignment)*(1<<alphaStride),JumpIdx(-1))
	const WordIdx* const encoderTablePointer;	
//	std::array<MarlinSymbol, 1U<<(sizeof(SourceSymbol)*8)> Source2JumpTableShifted;
//constexpr static const size_t FLAG_NEXT_WORD = 1UL<<(8*sizeof(WordIdx)-1);

};

}


#endif //MARLIN_DICTIONARY_H
