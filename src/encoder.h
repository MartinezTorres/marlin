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

#ifndef MARLIN_ENCODER_H
#define MARLIN_ENCODER_H

#include "dictionary.h"

namespace {

	class Encoder { 
		
		using SourceSymbol = Dictionary::SourceSymbol;
		using MarlinSymbol = Dictionary::MarlinSymbol;
		using Word = Dictionary::Word;

		typedef uint32_t JumpIdx;
		// Structured as:
		// FLAG_NEXT_WORD
		// Where to jump next		
		
		constexpr static const size_t FLAG_NEXT_WORD = 1UL<<(8*sizeof(JumpIdx)-1);
		
		class JumpTable {

			constexpr static const size_t unalignment = 8; // Too much aligned reads break cache
			const size_t alphaStride;  // Bit stride of the jump table corresponding to the word dimension
			const size_t wordStride;  // Bit stride of the jump table corresponding to the word dimension
		public:

			std::vector<JumpIdx> table;		
		
			JumpTable(size_t keySize, size_t overlap, size_t nAlpha) :
				alphaStride(std::ceil(std::log2(nAlpha))),
				wordStride(keySize+overlap),
				table(((1<<wordStride)+unalignment)*(1<<alphaStride),JumpIdx(-1))
				{}
			
			template<typename T0, typename T1>
			JumpIdx &operator()(const T0 &word, const T1 &nextLetter) { 
				return table[(word&((1<<wordStride)-1))+(nextLetter*((1<<wordStride)+unalignment))];
			}

			template<typename T0, typename T1>
			constexpr JumpIdx operator()(const T0 &word, const T1 &nextLetter) const { 
				return table[(word&((1<<wordStride)-1))+(nextLetter*((1<<wordStride)+unalignment))];
			}
		};
		JumpTable jumpTable;

		const size_t shift;
		const size_t nMarlinSymbols;
		std::array<MarlinSymbol, 1U<<(sizeof(SourceSymbol)*8)> Source2JumpTableShifted;
		MarlinSymbol Source2JumpTable(SourceSymbol ss) const {
			return Source2JumpTableShifted[ss>>shift];
		}

	public:
		
		Encoder(const Dictionary &dict, const Configuration &) :
			jumpTable(dict.K, dict.O, dict.alphabet.marlinSymbols.size()),
			shift(dict.alphabet.shift),
			nMarlinSymbols(dict.alphabet.marlinSymbols.size()) { 

			Source2JumpTableShifted.fill(nMarlinSymbols);
			for (size_t i=0; i<dict.alphabet.marlinSymbols.size(); i++)
				Source2JumpTableShifted[dict.alphabet.marlinSymbols[i].sourceSymbol>>shift] = i;

			
			const size_t NumSections = 1<<dict.O;
			const size_t SectionSize = 1<<dict.K;
			std::vector<std::map<Word, size_t>> positions(NumSections);

			// Init the mapping (to know where each word goes)
			for (size_t k=0; k<NumSections; k++)
				for (size_t i=k*SectionSize; i<(k+1)*SectionSize; i++)
					positions[k][dict.words[i]] = i;

			// Link each possible word to its continuation
			for (size_t k=0; k<NumSections; k++) {
				for (size_t i=k*SectionSize; i<(k+1)*SectionSize; i++) {
					auto word = dict.words[i];
					size_t wordIdx = i;
					while (not word.empty()) {
						SourceSymbol lastSymbol = word.back();						
						word.pop_back();
						if (not positions[k].count(word)) throw(std::runtime_error("SHOULD NEVER HAPPEN"));
						size_t parentIdx = positions[k][word];
						jumpTable(parentIdx, Source2JumpTable(lastSymbol)) = wordIdx;
						wordIdx = parentIdx;
					}
				}
			}
						
			//Link between inner dictionaries
			for (size_t k=0; k<NumSections; k++)
				for (size_t i=k*SectionSize; i<(k+1)*SectionSize; i++)
					for (size_t j=0; j<dict.alphabet.marlinSymbols.size(); j++)
						if (jumpTable(i,j)==JumpIdx(-1)) // words that are not parent of anyone else.
							jumpTable(i,j) = positions[i%NumSections][Word(1,dict.alphabet.marlinSymbols[j].sourceSymbol)] + FLAG_NEXT_WORD;
							
		}
		
		__attribute__ ((target ("bmi2")))
		size_t operator()(const uint8_t * const i8start, const uint8_t * const i8end, uint8_t * const o8start, uint8_t * const o8end) const {
			
			assert(o8end-o8start >= i8end-i8start);
			assert( (i8end-i8start)%8 == 0 );
			if (i8start==i8end) return 0;

			// Fast check to see if all the block is made of a single symbol
			{
				const uint8_t *i8test = i8start+1;
				while (i8test!=i8end and *i8test==i8start[0]) i8test++;
				if (i8test==i8end) {
					*o8start = i8start[0];
					return 1;
				}
			}

			uint8_t *o8 = o8start;
			const uint8_t *i8 = i8start;
			
			// Encode Marlin, with rare symbols preceded by an empty word
			{
				
				ssize_t maxTargetSize = std::max(0UL, (i8end-i8start)-((i8end-i8start)*shift/8));
				
				JumpIdx j = 0;
				while (i8<i8end and maxTargetSize>8) {				
					
					SourceSymbol ss = *i8++;
					
					MarlinSymbol ms = Source2JumpTable(ss);
					bool rare = ms==nMarlinSymbols;
					if (rare) {
						if (j) *o8++ = j; // Finish current word, if any;
						*o8++ = j = 0;
						*o8++ = (ss>>shift)<<shift;
						maxTargetSize-=3;
						continue;
					}
					
					JumpIdx jOld = j;
					j = jumpTable(j, ms);
					
					if (j & FLAG_NEXT_WORD) {
						*o8++ = jOld & 0xFF;
						maxTargetSize--;
					}
				}
				if (j) *o8++ = j;
				if (maxTargetSize <= 8) { // Just encode the block uncompressed.
					memcpy(o8start, i8start, i8end-i8start);
					return i8end-i8start;
				}
			}
			
			// Encode residuals
			if (shift) {
				uint64_t mask=0;
				for (size_t i=0; i<8; i++)
					mask |= ((1ULL<<shift)-1)<<(8ULL*i);
				
				const uint64_t *i64    = (const uint64_t *)i8start;
				const uint64_t *i64end = (const uint64_t *)i8end;

				while (i64 != i64end) {
					*(uint64_t *)o8 = _pext_u64(*i64++, mask);
					o8 += shift;
				}
			}
			return o8-o8start;
		}
	};
	
}

#endif //MARLIN_DECODER_H
