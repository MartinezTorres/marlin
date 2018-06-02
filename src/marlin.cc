/***********************************************************************

Marlin: A Fast Entropy Codec

MIT License

Copyright (c) 2017 Manuel Martinez Torres

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

#include <x86intrin.h>
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




namespace {
	
	typedef std::map<std::string, double> Configuration;

	class Dictionary {

		friend class Encoder;
		friend class Decoder;
		
		//Marlin only encodes a subset of the possible source symbols.
		//Marlin symbols are sorted by probability in descending order, 
		//so the Marlin Symbol 0 is always corresponds to the most probable alphabet symbol.
		typedef uint8_t SourceSymbol;
		typedef uint8_t MarlinSymbol; 
		
		// Finds the best configuration for a given alphabet
		static Configuration updateConf( 
			const std::map<SourceSymbol, double> &symbols, 
			Configuration conf) {

			conf.emplace("K",8);
			conf.emplace("O",2);
			
			conf.emplace("debug",1);
			conf.emplace("purgeProbabilityThreshold",1e-5);
			conf.emplace("iterations",3);
			conf.emplace("minMarlinSymbols", std::max(1U<<size_t(conf.at("O")),8U));
			conf.emplace("maxMarlinSymbols",(1U<<size_t(conf.at("K")))-1);
				
			if (not conf.count("S")) {
				conf["S"] = 0;
				double best = Dictionary(symbols, conf).efficiency;
				for (int s=1; s<6; s++) {
					conf["S"] = s;
					double e = Dictionary(symbols, conf).efficiency;
					if (e<=best) {
						conf["S"] = s-1;
						break;
					}
					best = e;
				}
			}
			
			if (not conf.count("maxWordSize")) {
				conf["maxWordSize"] = 15;
				double e15 = Dictionary(symbols, conf).efficiency;
				conf["maxWordSize"] = 7;
				double e7 = Dictionary(symbols, conf).efficiency;
				conf["maxWordSize"] = 3;
				double e3 = Dictionary(symbols, conf).efficiency;
				if (e7>1.0001*e3) {
					conf["maxWordSize"] = 7;
				}
				if (e15>1.0001*e7) {
					conf["maxWordSize"] = 15;
				}
			}
			
			return conf;
		}	

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
		

		struct Node;
		typedef std::shared_ptr<Node> SNode;	
		struct Node : std::vector<SNode> {
			double p=0;
			size_t sz=0;
		};
		


		SNode buildTree(std::vector<double> Pstates) const {

			// Normalizing the probabilities makes the algorithm more stable
			double factor = 1e-10;
			for (auto &&p : Pstates) factor += p;
			for (auto &&p : Pstates) p/=factor;
			for (auto &&p : Pstates) if (std::abs(p-1.)<0.0001) p=1.;
			for (auto &&p : Pstates) if (std::abs(p-0.)<0.0001) p=0.;


			std::vector<double> PN;
			for (auto &&a : alphabet.marlinSymbols) PN.push_back(a.p);
			PN.back() += alphabet.rareSymbolProbability;
			for (size_t i=PN.size()-1; i; i--)
				PN[i-1] += PN[i];

			std::vector<double> Pchild(PN.size());
			for (size_t i=0; i<PN.size(); i++)
				Pchild[i] = alphabet.marlinSymbols[i].p/PN[i];
			
			auto cmp = [](const SNode &lhs, const SNode &rhs) { return lhs->p<rhs->p;};
			std::priority_queue<SNode, std::vector<SNode>, decltype(cmp)> pq(cmp);

			// DICTIONARY INITIALIZATION
			SNode root = std::make_shared<Node>();
			
			// Include empty word
			pq.push(root);
			root->p = 1;
			
			for (size_t c=0; c<alphabet.marlinSymbols.size(); c++) {			
					
				root->push_back(std::make_shared<Node>());
				double sum = 0;
				for (size_t t = 0; t<=c; t++) sum += Pstates[t]/PN[t];
				root->back()->p = sum * alphabet.marlinSymbols[c].p;
				root->p -= root->back()->p;
				root->back()->sz = 1;
				pq.push(root->back());
			}
				
			// DICTIONARY GROWING
			size_t retiredNodes=0;
			while (not pq.empty() and (pq.size() + retiredNodes < (1U<<K))) {
					
				SNode node = pq.top();
				pq.pop();
				
				// retire words larger than maxWordSize that are meant to be extended by a symbol different than zero.
				if (node->sz >= maxWordSize and not node->empty()) {
					retiredNodes++;
					continue;
				}

				if (node->sz == 255) {
					retiredNodes++;
					continue;
				}
				
				if (node->size() == alphabet.marlinSymbols.size()) {
					retiredNodes++;
					continue;					
				}
				
				double p = node->p * Pchild[node->size()];
				node->push_back(std::make_shared<Node>());
				node->back()->p = p;
				node->back()->sz = node->sz+1;
				pq.push(node->back());
				node->p -= p;
				pq.push(node);
			}

			// Renormalize probabilities.
			{
				std::queue<SNode> q(std::deque<SNode>{ root });
				double sum=0, num=0;
				while (not q.empty()) {
					sum += q.front()->p; num++;
					q.front()->p *= factor;
					for (auto &&child : *q.front()) 
						q.push(child);
					q.pop();
				}
				std::cerr << sum << " sum - num " << num << std::endl;
			}
			return root;
		}


		
		std::vector<Word> buildChapterWords( const SNode root) const {
		
			std::vector<Word> ret;
			
			std::stack<std::pair<SNode, Word>> q;
			Word rootWord;
			rootWord.p = root->p;
			q.emplace(root, rootWord);
			
			while (not q.empty()) {
				SNode n = q.top().first;
				Word w = q.top().second;
				q.pop();
				ret.push_back(w);
				for (size_t i = 0; i<n->size(); i++) {
					
					Word w2 = w;
					w2.push_back(alphabet.marlinSymbols[i].sourceSymbol);
					w2.p = n->at(i)->p;
					w2.state = n->at(i)->size();
					
					assert(n->at(i)->sz == w2.size());
					q.emplace(n->at(i), w2);
				}
			}
			
			std::cout << ret.size() << std::endl;
			return ret;
		}


		
		std::vector<Word> arrangeAndFuse( const std::vector<SNode> chapters) const {

			std::vector<Word> ret;
			for (auto &&chapter : chapters) {
				
				std::vector<Word> sortedDictionary = buildChapterWords(chapter);
				
				auto cmp = [](const Word &lhs, const Word &rhs) { 
					if (lhs.state != rhs.state) return lhs.state<rhs.state;
					if (std::abs(lhs.p-rhs.p)/(lhs.p+rhs.p) > 1e-10) return lhs.p > rhs.p;
					return lhs<rhs;
				};
				// Note the +1, we keep the empty word in the first position.
				std::stable_sort(sortedDictionary.begin()+1, sortedDictionary.end(), cmp);
				
				std::vector<Word> w(1U<<K,Word());
				for (size_t i=0,j=0,k=0; i<sortedDictionary.size(); j+=(1U<<O)) {
					
					if (j>=w.size()) 
						j=++k;

					w[j] = sortedDictionary[i++];
				}
				ret.insert(ret.end(), w.begin(), w.end());
			}
			return ret;
		}
		

		// Debug functions		
		void print(std::vector<Word> dictionary) const {

			if (conf.at("debug")<3) return;
			if (conf.at("debug")<4 and dictionary.size()/(1U<<O) > 40) return;

			for (size_t i=0; i<dictionary.size()/(1U<<O); i++) { 
				
				for (size_t k=0; k<(1U<<O); k++) {
					
					auto idx = i + (k* (dictionary.size()/(1U<<O)));
					auto &&w = dictionary[idx];
					printf(" %02lX %01ld %2d %01.2le ",idx,i%(1U<<O),w.state,w.p);
					for (size_t j=0; j<8; j++) {
						if (j<w.size()) {
							char a = 'a';
							for (size_t x=0; alphabet.marlinSymbols[x].sourceSymbol != w[j]; x++, a++);
							putchar(a);
						} else {
							putchar(' ');
						}
					}
				}
				putchar('\n');
			}		
			putchar('\n');
		}



		void print(std::vector<std::vector<double>> Pstates) const {
			
			if (conf.at("debug")<3) return;
			for (size_t i=0; i<Pstates[0].size() and i<4; i++) { 
				
				printf("S: %02ld",i);
				for (size_t k=0; k<Pstates.size() and k<8; k++) 
						 printf(" %01.3lf",Pstates[k][i]);
				putchar('\n');
			}		
			putchar('\n');
		}



		double calcEfficiency(std::vector<Word> dictionary) const {
		
			double meanLength = 0;
			for (auto &&w : dictionary)
					meanLength += w.p * w.size();
			
			double shannonLimit = alphabet.sourceEntropy;
			
			// The decoding algorithm has 4 steps:
			double meanBitsPerSymbol = 0;                           // a memset
			meanBitsPerSymbol += (K/meanLength)*(1-alphabet.rareSymbolProbability);                      // Marlin VF
			meanBitsPerSymbol += alphabet.shift;                    // Raw storing of lower bits
			meanBitsPerSymbol += 2*K*alphabet.rareSymbolProbability;// Recovering rare symbols

			return shannonLimit / meanBitsPerSymbol;
		}

		

		std::vector<Word> buildDictionary() const {
			
			std::vector<std::vector<double>> Pstates;
			for (size_t k=0; k<(1U<<O); k++) {
				std::vector<double> PstatesSingle(alphabet.marlinSymbols.size(), 0.);
				PstatesSingle[0] = 1./(1U<<O);
				Pstates.push_back(PstatesSingle);
			}
			
			std::vector<SNode> dictionaries;
			for (size_t k=0; k<(1U<<O); k++)
				dictionaries.push_back(buildTree(Pstates[k]));
				
			std::vector<Word> ret = arrangeAndFuse(dictionaries);
				
			print(ret);
			
			size_t iterations = conf.at("iterations");
				
			while (iterations--) {

				// UPDATING STATE PROBABILITIES
				{
					for (auto &&pk : Pstates)
						for (auto &&p : pk)
							p = 0.;

					for (size_t i=0; i<ret.size(); i++)
						Pstates[i%(1U<<O)][ret[i].state] += ret[i].p;
				}
				
				print(Pstates);

				dictionaries.clear();
				for (size_t k=0; k<(1U<<O); k++)
					dictionaries.push_back(buildTree(Pstates[k]));
				
				ret = arrangeAndFuse(dictionaries);
				
				print(ret);
				if (conf.at("debug")>2) printf("Efficiency: %3.4lf\n", calcEfficiency(ret));		
			}
			if (conf.at("debug")>1) for (auto &&c : conf) std::cout << c.first << ": " << c.second << std::endl;
			if (conf.at("debug")>0) printf("Efficiency: %3.4lf\n", calcEfficiency(ret));

			return ret;
		}

	public:
		
		const Configuration conf;
		const Alphabet alphabet;
		const size_t K                = conf.at("K");           // Non overlapping bits of codeword.
		const size_t O                = conf.at("O");           // Bits that overlap between codewprds.
		const size_t maxWordSize      = conf.at("maxWordSize"); // Maximum number of symbols per word.
		const std::vector<Word> words = buildDictionary();      // All dictionary words.
		const double efficiency       = calcEfficiency(words);  // Theoretical efficiency of the dictionary.

		Dictionary( const std::map<SourceSymbol, double> &symbols,
			Configuration conf_ = Configuration()) 
			: conf(updateConf(symbols, conf_)), alphabet(symbols, conf) {}
		
		
		// Turns the vector into a map and uses the previous constructor
		Dictionary( const std::vector<double> &symbols,
			Configuration conf_ = Configuration()) 
			: Dictionary(
			[&symbols](){
				std::map<SourceSymbol, double> ret;
				for (size_t i=0; i<symbols.size(); i++)
					ret.emplace(SourceSymbol(i), symbols[i]);
				return ret;
			}()
			, conf_) {}
	};

	
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
	

	class Decoder {

		using SourceSymbol = Dictionary::SourceSymbol;
		using MarlinSymbol = Dictionary::MarlinSymbol;
		
		const size_t shift;
		const size_t O;
		const size_t maxWordSize;
		
		std::vector<SourceSymbol> decoderTable;
		const SourceSymbol * const D;
		const SourceSymbol mostCommonSourceSymbol;

		
		template<typename T, size_t CO>
		size_t decode8(const uint8_t * const i8start, const uint8_t * const i8end, uint8_t * const o8start, uint8_t * const o8end) const {
			
			      uint8_t *o8 = o8start;
			const uint8_t *i8 = i8start;
			
			// Special case, same size! this means the block is uncompressed.
			if (i8end-i8start == o8end-o8start) {
				memcpy(o8start,i8start,i8end-i8start);
				return o8end-o8start;
			}

			// Special case, size 1! this means the block consists all of just one symbol.
			if (i8end-i8start == 1) {
				memset(o8start,*i8start,o8end-o8start);
				return o8end-o8start;
			}

			memset(o8start,mostCommonSourceSymbol,o8end-o8start);
//			return o8end-o8start;
			
			// Decode the Marlin Section
			{

				const uint8_t *endMarlin = i8end - (o8end-o8start)*shift/8;

//				const uint32_t overlappingMask = (1<<(8+O))-1;
				constexpr const uint32_t overlappingMask = (1<<(8+CO))-1;
//				constexpr const T clearSizeMask = T(-1)>>8;
				constexpr const T clearSizeMask = 0;
				uint64_t value = 0;

				while (i8<endMarlin-9) {
					
					uint32_t v32 = (*(const uint32_t *)i8);
/*					if (((v32 - 0x01010101UL) & ~v32 & 0x80808080UL)) { // Fast test for zero

						uint8_t in = *i8++;
						if (in==0) {
							*o8++ = *i8++;
							value = (value<<8) + 0;
						} else {
							value = (value<<8) + in;
							T v = ((const T *)D)[value & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}
						
						in = *i8++;
						if (in==0) {
							*o8++ = *i8++;
							value = (value<<8) + 0;
						} else {
							value = (value<<8) + in;
							T v = ((const T *)D)[value & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}
						
						in = *i8++;
						if (in==0) {
							*o8++ = *i8++;
							value = (value<<8) + 0;
						} else {
							value = (value<<8) + in;
							T v = ((const T *)D)[value & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}
						
						in = *i8++;
						if (in==0) {
							*o8++ = *i8++;
							value = (value<<8) + 0;
						} else {
							value = (value<<8) + in;
							T v = ((const T *)D)[value & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}
						
					} else { // Has no zeroes! hurray!*/
						i8+=4;
						//clearSizeMask = 0;
						value = (value<<32) +  v32; //__builtin_bswap32(v32);
						{
							T v = ((const T *)D)[(value>>24) & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
							
						}

						{
							T v = ((const T *)D)[(value>>16) & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}

						{
							T v = ((const T *)D)[(value>>8) & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}

						{
							T v = ((const T *)D)[value & overlappingMask];
							*((T *)o8) = v & clearSizeMask;
							o8 += v >> ((sizeof(T)-1)*8);
						}
					//}
				}
				
				while (i8<endMarlin) {
					uint8_t in = *i8++;
					if (in==0) {
						*o8++ = *i8++;
					} else {
						value = (value<<8) + in;
						const T *v = &((const T *)D)[value & overlappingMask];
						memcpy(o8, v, std::min(sizeof(T)-1,size_t(*v >> ((sizeof(T)-1)*8))));
						o8 += *v >> ((sizeof(T)-1)*8);
					}
				}				
				//if (endMarlin-i8 != 0) std::cerr << " {" << endMarlin-i8 << "} "; // SOLVED! PROBLEM IN THE CODE
				//if (o8end-o8 != 0) std::cerr << " [" << o8end-o8 << "] "; // SOLVED! PROBLEM IN THE CODE
			}

			// Decode residuals
			if (shift) {
				uint64_t mask=0;
				for (size_t i=0; i<8; i++)
					mask |= ((1ULL<<shift)-1)<<(8ULL*i);
				
				uint64_t *o64    = (uint64_t *)o8start;
				uint64_t *o64end = (uint64_t *)o8end;

				while (o64 != o64end) {
					*o64++ += _pdep_u64(*(const uint64_t *)i8, mask);
					i8 += shift;
				}
			}
			return o8end-o8start;
		}
		
	public:
		
		
		Decoder(const Dictionary &dict, const Configuration &) :
			shift(dict.alphabet.shift),
			O(dict.O),
			maxWordSize(dict.maxWordSize),
			decoderTable(dict.words.size()*(maxWordSize+1)),
			D(decoderTable.data()),
			mostCommonSourceSymbol(dict.alphabet.marlinSymbols.front().sourceSymbol) {
				
			
			SourceSymbol *d = &decoderTable.front();
			for (size_t i=0; i<dict.words.size(); i++) {
				for (size_t j=0; j<maxWordSize; j++)
					*d++ = (dict.words[i].size()>j ? dict.words[i][j] : SourceSymbol(0));
				*d++ = dict.words[i].size();
			}
		}
		
		size_t operator()(const uint8_t * const i8start, const uint8_t * const i8end, uint8_t * const o8start, uint8_t * const o8end) const {
			
			if (maxWordSize==3) {
				switch (O) {
					case   0: return decode8<uint32_t,0>(i8start, i8end, o8start, o8end);
					case   1: return decode8<uint32_t,1>(i8start, i8end, o8start, o8end);
					case   2: return decode8<uint32_t,2>(i8start, i8end, o8start, o8end);
					case   3: return decode8<uint32_t,3>(i8start, i8end, o8start, o8end);
					case   4: return decode8<uint32_t,4>(i8start, i8end, o8start, o8end);
				}
			}

			if (maxWordSize==7) {
				switch (O) {
					case   0: return decode8<uint64_t,0>(i8start, i8end, o8start, o8end);
					case   1: return decode8<uint64_t,1>(i8start, i8end, o8start, o8end);
					case   2: return decode8<uint64_t,2>(i8start, i8end, o8start, o8end);
					case   3: return decode8<uint64_t,3>(i8start, i8end, o8start, o8end);
					case   4: return decode8<uint64_t,4>(i8start, i8end, o8start, o8end);
				}
			}
			throw std::runtime_error ("unsupported maxWordSize");	
		}
	};
	

}

struct MarlinDictionary {
	
//	const std::string name;
//	const double hist[256];
//	const double bps[256]; // Expected bits per symbol
	
	const double efficiency;
	const Encoder encoder;
	const Decoder decoder;
	
	MarlinDictionary(const Dictionary &dict, const Configuration &configuration = Configuration()) : 
			efficiency(dict.efficiency),
			encoder(dict, configuration), 
			decoder(dict, configuration) {}
			
	MarlinDictionary(const std::vector<double> &pdf, const Configuration &configuration = Configuration()) : 
			MarlinDictionary(Dictionary(pdf, configuration), configuration) {}
	
};


////////////////////////////////////////////////////////////////////////
//
// Local Methods


////////////////////////////////////////////////////////////////////////
//
// Public Methods

size_t Marlin_compress(uint8_t* dst, size_t dstCapacity, const uint8_t* src, size_t srcSize, const MarlinDictionary *dict) {
	
	return dict->encoder(src, src+srcSize, dst, dst+dstCapacity);
}

size_t Marlin_decompress(uint8_t* dst, size_t dstSize, const uint8_t* src, size_t srcSize, const MarlinDictionary *dict) {
	
	return dict->decoder(src, src+srcSize, dst, dst+dstSize);
}

MarlinDictionary *Marlin_build_dictionary(const char *name, const double hist[256], size_t indexSizeBits, size_t indexOverlapBits, size_t maxWordSizeSymbols, size_t rawStorageBits) {
}

void Marlin_free_dictionary(MarlinDictionary *&dict) {
	
	free(dict);
	dict = nullptr;
}

double Marlin_estimate_size(const double hist[256], MarlinDictionary *dict) {
	
	double ret = 0;
	//for (int i=0; i<256; i++)
	//	ret += hist[i]*dict->bps[i];
	return ret;
}


