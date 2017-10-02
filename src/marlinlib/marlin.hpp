#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <stack>

#include <memory>
#include <algorithm>

#include <util/distribution.hpp>
#include <cassert>

class Marlin2018Simple {
	
	// Configuration		
	static const constexpr bool enableVictimDictionary = true;
	static const constexpr double purgeProbabilityThreshold = 1e-10;
	static const constexpr size_t iterationLimit = 3;
	static const constexpr bool debug = false;

	typedef uint8_t Symbol; // storage used to store an input symbol.
	typedef uint16_t WordIdx; // storage that suffices to store a word index.

	struct SymbolAndProbability {
		Symbol symbol;
		double p;
		bool operator<(const SymbolAndProbability &rhs) {
			if (p!=rhs.p) return p>rhs.p; // Descending in probability
			return symbol<rhs.symbol; // Ascending in symbol index
		}
	};
	
	struct Alphabet : public std::vector<SymbolAndProbability> {
		
		Alphabet(const std::map<Symbol, double> &symbols) {
			for (auto &&symbol : symbols)
				this->push_back(SymbolAndProbability({symbol.first, symbol.second}));
			std::sort(this->begin(),this->end());
		}
		Alphabet(const std::vector<double> &symbols) {
			for (size_t i=0; i<symbols.size(); i++)
				this->push_back(SymbolAndProbability({Symbol(i), symbols[i]}));
			std::sort(this->begin(),this->end());
		}
	};
	
	struct Word : std::vector<Symbol> {
		using std::vector<Symbol>::vector;
		double p = 0;
		Symbol state = 0;
		
		friend std::ostream& operator<< (std::ostream& stream, const Word& word) {
			stream << "{";
			for (size_t i=0; i<word.size(); i++) {
				if (i) stream << ",";
				stream << int(word[i]);
			}
			stream << "}";
			return stream;
        }
	};

	class Dictionary : public std::vector<Word> {
		
		struct Node;		
		struct Node : std::vector<std::shared_ptr<Node>> {
			double p=0;
			size_t sz=0;
			size_t erased=0;
		};

		std::shared_ptr<Node> buildTree(const std::vector<double> &Pstates, bool isVictim) const {

			std::vector<double> PN;
			for (auto &&a : alphabet) PN.push_back(a.p);
			for (size_t i=alphabet.size()-1; i; i--)
				PN[i-1] += PN[i];

			std::vector<double> Pchild(alphabet.size());
			for (size_t i=0; i<alphabet.size(); i++)
				Pchild[i] = alphabet[i].p/PN[i];

			
			auto cmp = [](const std::shared_ptr<Node> &lhs, const std::shared_ptr<Node> &rhs) { return lhs->p < rhs->p;};
		
//			auto cmp = [](const std::shared_ptr<Node> &lhs, const std::shared_ptr<Node> &rhs)
//				return lhs->p*(1+std::pow(lhs->sz,1)) < rhs->p*(1+std::pow(rhs->sz,1));	};

			std::priority_queue<std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>, decltype(cmp)> pq(cmp);
			size_t retiredNodes=0;

			auto pushAndPrune = [this,&pq,&retiredNodes,isVictim](std::shared_ptr<Node> node) {
				if (isVictim or (not Marlin2018Simple::enableVictimDictionary) or node->p>Marlin2018Simple::purgeProbabilityThreshold) {
					if (node->sz<=maxWordSize) {
						pq.push(node);
					} else {
						retiredNodes++;
					}
				} else {
					node->erased = true;
				}
			};

			// DICTIONARY INITIALIZATION
			std::shared_ptr<Node> root = std::make_shared<Node>();
			root->erased = true;
			
			// Include empty word
			pq.push(root); // Does not do anything, only uses a spot

			for (size_t c=0; c<alphabet.size(); c++) {			
					
				root->push_back(std::make_shared<Node>());
				double sum = 0;
				for (size_t t = 0; t<=c; t++) sum += Pstates[t]/PN[t];
				root->back()->p = sum * alphabet[c].p;
				root->back()->sz = 1;
				pushAndPrune(root->back());
			}
				
			// DICTIONARY GROWING
			while (not pq.empty() and (pq.size() + retiredNodes < (1U<<keySize))) {
					
				std::shared_ptr<Node> node = pq.top();
				pq.pop();
				
				//if (node->sz>2) 
					
				double p = node->p * Pchild[node->size()];
				node->push_back(std::make_shared<Node>());
				node->back()->p = p;
				node->back()->sz = node->sz+1;

				node->p -= p;
				pushAndPrune(node->back());
					
				if (node->size()<alphabet.size()-1) {

					pushAndPrune(node);
						
				} else {
					node->erased = true;
					node->p = 0;

					node->push_back(std::make_shared<Node>());
					node->back()->p = node->p;
					node->back()->sz = node->sz+1;
					pushAndPrune(node->back());
				}
			}
			return root;
		}
		
		std::vector<Word> buildWords( const std::shared_ptr<Node> &root) const {
		
			std::vector<Word> ret;
			
			std::stack<std::pair<std::shared_ptr<Node>, Word>> q;
			q.emplace(root, Word());
			while (not q.empty()) {
				std::shared_ptr<Node> n = q.top().first;
				Word w = q.top().second;
				q.pop();
				if (not n->erased) ret.push_back(w);
				for (size_t i = 0; i<n->size(); i++) {
					
					Word w2 = w;
					w2.push_back(alphabet[i].symbol);
					w2.p = n->at(i)->p;
					w2.state = n->at(i)->size();
					q.emplace(n->at(i), w2);
				}
			}
			return ret;
		}
		
		std::vector<Word> arrangeAndFuse( const std::vector<std::shared_ptr<Node>> &nodes, size_t victimIdx ) const {
			
			std::vector<Word> ret;
			for (size_t n = 0; n<nodes.size(); n++) {
				
				std::vector<Word> sortedDictionary = buildWords(nodes[n]);
				auto cmp = [](const Word &lhs, const Word &rhs) { 
					if (lhs.state != rhs.state) return lhs.state<rhs.state;
					return lhs.p > rhs.p;
				};
				std::sort(sortedDictionary.begin(), sortedDictionary.end(), cmp);
				
				std::vector<Word> w(1<<keySize);
				for (size_t i=0,j=0,k=0; i<sortedDictionary.size(); j+=(1<<overlap)) {
					
					if (j>=w.size()) 
						j=++k;
						
					if (victimIdx==j) {
						w[j] = Word();
					} else {
						w[j] = sortedDictionary[i++];
					}
				}
				ret.insert(ret.end(), w.begin(), w.end());
			}
			return ret;
		}
		
		// Debug functions
		
		void print(std::vector<Word> dictionary) {
			
			if (dictionary.size()>40) return;

			for (size_t i=0; i<dictionary.size()/(1U<<overlap); i++) { 
				
				for (size_t k=0; k<(1U<<overlap); k++) {
					
					auto idx = i + (k* (dictionary.size()/(1U<<overlap)));
					auto &&w = dictionary[idx];
					printf(" %02lX %01ld %2d %01.3lf ",idx,i%(1<<overlap),w.state,w.p);
					for (size_t j=0; j<16; j++) putchar("0123456789ABCDEF "[j<w.size()?w[j]:16]);
				}
				putchar('\n');
			}		
			putchar('\n');
		}

		static void print(std::vector<std::vector<double>> Pstates) {
			
			for (size_t i=0; i<Pstates[0].size() and i<4; i++) { 
				
				printf("S: %02ld",i);
				for (size_t k=0; k<Pstates.size() and k<8; k++) 
						 printf(" %01.3lf",Pstates[k][i]);
				putchar('\n');
			}		
			putchar('\n');
		}
			
		
	public:

		const Alphabet alphabet;
		const size_t keySize;     // Non overlapping bits of the word index in the big dictionary.
		const size_t overlap;     // Bits that overlap between keys.s
		const size_t maxWordSize; // Maximum number of symbols that a word in the dictionary can have.

		double calcEfficiency() const {
		
			double meanLength = 0;
			for (auto &&w : *this)
				meanLength += w.p * w.size();
			
			std::vector<double> P;
			for (auto &&a: alphabet) P.push_back(a.p);
			double shannonLimit = Distribution::entropy(P)/std::log2(P.size());

			return shannonLimit / (keySize / (meanLength*std::log2(P.size())));
		}
		
		Dictionary(const Alphabet &alphabet, size_t keySize, size_t overlap, size_t maxWordSize)
			: alphabet(alphabet), keySize(keySize), overlap(overlap), maxWordSize(maxWordSize) {
			
			std::vector<std::vector<double>> Pstates;
			for (auto k=0; k<(1<<overlap); k++) {
				std::vector<double> PstatesSingle(alphabet.size(), 0.);
				PstatesSingle[0] = 1./(1<<overlap);
				Pstates.push_back(PstatesSingle);
			}
			
			int victimDictionary = 0;
			
			std::vector<std::shared_ptr<Node>> dictionaries;
			for (auto k=0; k<(1<<overlap); k++)
				dictionaries.push_back(buildTree(Pstates[k], k==victimDictionary) );
				
			*(std::vector<Word> *)this = arrangeAndFuse(dictionaries,victimDictionary);
				
			if (Marlin2018Simple::debug) print(*this);
				
			for (size_t iteration=0; iteration<Marlin2018Simple::iterationLimit; iteration++) {

				// UPDATING STATE PROBABILITIES
				{
					for (auto k=0; k<(1<<overlap); k++)
						Pstates[k] = std::vector<double>(alphabet.size(), 0.);

					for (size_t i=0; i<size(); i++)
						Pstates[i%(1<<overlap)][(*this)[i].state] += (*this)[i].p;
				}
				
				// Find least probable subdictionary
				{
					double minP = 1.1;
					for (size_t i=0; i<Pstates.size(); i++) {
						double sumProb = 0.;
						for (auto &&ps : Pstates[i])
							sumProb += ps;
						if (sumProb > minP) continue;
						minP = sumProb;
						victimDictionary = i;
					}
				}
				if (Marlin2018Simple::debug) print(Pstates);

				dictionaries.clear();
				for (auto k=0; k<(1<<overlap); k++)
					dictionaries.push_back(buildTree(Pstates[k], k==victimDictionary) );
				
				*(std::vector<Word> *)this = arrangeAndFuse(dictionaries,victimDictionary);
				
				if (Marlin2018Simple::debug) print(*this);
				if (Marlin2018Simple::debug) printf("Efficiency: %3.4lf\n", calcEfficiency());		
			}
			if (Marlin2018Simple::debug) printf("Efficiency: %3.4lf\n", calcEfficiency());				
		}			
	};
	const Dictionary dictionary;
	
	struct Encoder { 

		typedef uint32_t JumpIdx;
		// Structured as:
		// FLAG_NEXT_WORD
		// FLAG_INSERT_EMPTY_WORD
		// Current Dictionary
		// Where to jump next
		
		constexpr static const size_t FLAG_NEXT_WORD = 1<<(8*sizeof(JumpIdx)-1);
		constexpr static const size_t FLAG_INSERT_EMPTY_WORD = 1<<(8*sizeof(JumpIdx)-2);
		
		class JumpTable {

			const size_t alphaStride;  // Bit stride of the jump table corresponding to the word dimension
			const size_t wordStride;  // Bit stride of the jump table corresponding to the word dimension
			size_t nextIntermediatePos = 1<<(wordStride-1);	
			std::vector<JumpIdx> table;
		
		public:
		
			JumpTable(size_t keySize, size_t overlap, size_t nAlpha) :
				alphaStride(std::ceil(std::log2(nAlpha))),
				wordStride(keySize+overlap+1), // Extra bit for intermediate nodes.
				table((1<<wordStride)*nAlpha,JumpIdx(-1))
				{}
			
			template<typename T0, typename T1>
			JumpIdx &operator()(const T0 &word, const T1 &nextLetter) { 
//				if ((word&((1<<wordStride)-1))+(nextLetter<<wordStride) < 0) std::cerr << "Underrun" << std::endl;
//				if ((word&((1<<wordStride)-1))+(nextLetter<<wordStride) >= table.size()) std::cerr << "Overrun" << std::endl;

//				return table[((word&((1<<wordStride)-1))<<alphaStride) + nextLetter];
				return table[(word&((1<<wordStride)-1))+(nextLetter<<wordStride)];
			}

			template<typename T0, typename T1>
			JumpIdx operator()(const T0 &word, const T1 &nextLetter) const { 
//				if ((word&((1<<wordStride)-1))+(nextLetter<<wordStride) < 0) std::cerr << "Underrun" << std::endl;
//				if ((word&((1<<wordStride)-1))+(nextLetter<<wordStride) >= table.size()) std::cerr << "Overrun" << std::endl;
//				return table[((word&((1<<wordStride)-1))<<alphaStride) + nextLetter];
				return table[(word&((1<<wordStride)-1))+(nextLetter<<wordStride)];
			}
			
			size_t getNewPos() { return nextIntermediatePos++; }
			
			bool isIntermediate(size_t pos) const { return pos & (1<<(wordStride-1)); }
		};
		JumpTable jumpTable;
		
		JumpIdx start;
		std::vector<JumpIdx> emptyWords; //emptyWords are pointers to the victim dictionary.
		const Dictionary dict;

		Encoder(const Dictionary &dict) :
			jumpTable(dict.keySize, dict.overlap, dict.alphabet.size()),
			dict(dict) { 
			
			size_t NumSections = 1<<dict.overlap;
			size_t SectionSize = 1<<dict.keySize;
			std::vector<std::map<Word, size_t>> positions(NumSections);

			// Init the mapping (to know where each word goes)
			for (size_t k=0; k<NumSections; k++)
				for (size_t i=k*SectionSize; i<(k+1)*SectionSize; i++)
					positions[k][dict[i]] = i;

			// Link each possible word to its continuation
			for (size_t k=0; k<NumSections; k++) {
				for (size_t i=k*SectionSize; i<(k+1)*SectionSize; i++) {
					Word word = dict[i];
					size_t wordIdx = i;
//					std::cerr << "start " << k << " " << i << " " << dict.size() << " " << word << std::endl;
					while (not word.empty()) {
						auto lastSymbol = word.back();						
						word.pop_back();
						size_t parentIdx;
						if (positions[k].count(word)) {
							parentIdx = positions[k][word];
						} else {
							std::cerr << "SHOULD NEVER HAPPEN" << std::endl;
							parentIdx = jumpTable.getNewPos();
						}
						jumpTable(parentIdx, lastSymbol) = wordIdx;
						wordIdx = parentIdx;
					}
				}
			}
						
			//Link between inner dictionaries
			for (size_t k=0; k<NumSections; k++) {
				for (size_t i=k*SectionSize; i<(k+1)*SectionSize; i++) {
					for (size_t j=0; j<dict.alphabet.size(); j++) {
						if (jumpTable(i,j)==JumpIdx(-1)) {
							if (positions[i%(1<<dict.overlap)].count(Word(1,Symbol(j)))) {
								jumpTable(i, j) = positions[i%(1<<dict.overlap)][Word(1,Symbol(j))] +
									FLAG_NEXT_WORD;
							} else { 
//								std::cerr << "k" << i << " " << j << std::endl;
								jumpTable(i, j) = positions[i%(1<<dict.overlap)][Word()] +
									FLAG_NEXT_WORD + 
									FLAG_INSERT_EMPTY_WORD;
							}
						}
					}
				}
			}
			
			// Fill list of empty words
			emptyWords.resize(NumSections,JumpIdx(-1));
			for (size_t k=0; k<NumSections; k++)
				emptyWords[k] = positions[k][Word()];

			// Get Starting Positions (not encoded)
			size_t victim = 0;
			while (not dict[victim].empty()) 
				victim++;
			victim = victim % (1<<dict.overlap);

			start = victim*SectionSize;
			while (not dict[start].empty()) 
				start++;
		}
		
		template<class TIN, typename TOUT, typename std::enable_if<sizeof(typename TIN::value_type)==1,int>::type = 0>		
		void encodeA(const TIN &in, TOUT &out) const {
			
			if (out.size() < in.size()) out.resize(in.size());
			
			uint32_t *o = (uint32_t *)&out.front();
			const uint8_t *i = (const uint8_t *)&in.front();
			const uint8_t *iend = i + in.size();
			
			uint64_t value=0; int32_t bits=0;
			if (i<iend) {
				
				JumpIdx j0 = jumpTable(start, *i++);
				while (i<iend) {

					JumpIdx j1 = jumpTable(j0, *i++);
					
					if (j1 & FLAG_NEXT_WORD) {
						value <<= dict.keySize;
						bits += dict.keySize;
						value += j0 & ((1<<dict.keySize)-1);
						if (bits>=32) {
							bits -= 32;
							*o++ = value>>bits;
						}

						if (j1 & FLAG_INSERT_EMPTY_WORD) i--;
					}
					j0=j1;
				}
				assert (not jumpTable.isIntermediate(j0)); //If we end in an intermediate node, we should roll back. Not implemented.
				value <<= dict.keySize;
				bits += dict.keySize;
				value += j0 & ((1<<dict.keySize)-1);

//std::cerr << (j0 & ((1<<dict.keySize)-1)) << " " << bits << std::endl;
				
				while (bits) {
					while (bits<32) {
						j0 = emptyWords[j0 % emptyWords.size()];
						value <<= dict.keySize;
						bits += dict.keySize;
						value += j0 & ((1<<dict.keySize)-1);

//std::cerr << (j0 & ((1<<dict.keySize)-1)) << " " << bits << std::endl;
					}
					bits -= 32;
					*o++ = value>>bits;
				}
			}
			out.resize((uint8_t *)o-(uint8_t *)&out.front());
		}

		template<class TIN, typename TOUT, typename std::enable_if<sizeof(typename TIN::value_type)==1,int>::type = 0>		
		void encode12(const TIN &in, TOUT &out) const {
			
			if (out.size() < in.size()) out.resize(in.size());
			
			uint32_t *o = (uint32_t *)&out.front();
			const uint8_t *i = (const uint8_t *)&in.front();
			const uint8_t *iend = i + in.size();
			
			uint64_t value=0; int32_t bits=0;
			if (i<iend) {
				
				JumpIdx j0 = jumpTable(start, *i++);
				while (i<iend) {
					
					JumpIdx j1 = jumpTable(j0, *i++);
					if (j1 & FLAG_NEXT_WORD) {
						if (j1 & FLAG_INSERT_EMPTY_WORD)  i--;
						value <<= 12;
						bits += 12;
						value += j0 & ((1<<12)-1);
						if (bits>=32) {
							bits -= 32;
							*o++ = value>>bits;
						}
					}
					j0=j1;
				}
				assert (not jumpTable.isIntermediate(j0)); //If we end in an intermediate node, we should roll back. Not implemented.
				value <<= dict.keySize;
				bits += dict.keySize;
				value += j0 & ((1<<dict.keySize)-1);

//std::cerr << (j0 & ((1<<dict.keySize)-1)) << " " << bits << std::endl;
				
				while (bits) {
					while (bits<32) {
						j0 = emptyWords[j0 % emptyWords.size()];
						value <<= dict.keySize;
						bits += dict.keySize;
						value += j0 & ((1<<dict.keySize)-1);

//std::cerr << (j0 & ((1<<dict.keySize)-1)) << " " << bits << std::endl;
					}
					bits -= 32;
					*o++ = value>>bits;
				}
			}
			out.resize((uint8_t *)o-(uint8_t *)&out.front());
		}

		template<class TIN, typename TOUT, typename std::enable_if<sizeof(typename TIN::value_type)==1,int>::type = 0>		
		void operator()(const TIN &in, TOUT &out) const {
			if (dict.keySize==12) 
				encode12(in,out);
			else
				encodeA(in,out);
		}
	};
	
	struct EncoderSlow {
		
		const Dictionary W;
		
		EncoderSlow(const Dictionary &dict) : W(dict) {}
		
		template<typename T, typename IT>
		static bool areEqual(const T &obj, const IT it) {
			for (size_t i=0; i<obj.size(); i++)
				if (obj[i] != it[i])
					return false;
			return true;
		}

		void operator()(const std::vector<uint8_t> in, std::vector<uint8_t> & out) const {

			out.clear();
			uint64_t value=0; int32_t bits=0;
			
			size_t lastWord = 0;
			while (not W[lastWord].empty()) 
				lastWord++;
			
			for (size_t i=0, longest=0; i<in.size(); i+=longest) {
				
				size_t remaining = in.size()-i;
				
				size_t best = 0;
				longest = 0;
				for (size_t j=0; j<W.size()>>W.overlap; j++) {
					size_t idx = (lastWord%(1<<W.overlap))*(W.size()>>W.overlap) + j;
					if (W[idx].size()>longest and W[idx].size()<= remaining and areEqual(W[idx],in.begin()+i) ) {
						best = idx;
						longest = W[idx].size();
					}
				}
				
				value <<= W.keySize;
				bits += W.keySize;
				value += best & ((1<<W.keySize)-1);
				lastWord = best;
				std::cerr << best << ":" << bits << std::endl;

				if (bits>=32) {
					bits-=32;
					out.resize(out.size()+4);
					uint32_t *v = (uint32_t *)&*out.end();
					*--v = value>>bits;
				}
			}
			while (bits) {
				while (bits<32) {
					for (size_t j=0; j<W.size()>>W.overlap; j++) {
						size_t idx = (lastWord%(1<<W.overlap))*(W.size()>>W.overlap) + j;
						if (W[idx].empty()) {
							value <<= W.keySize;
							bits += W.keySize;
							value += idx & ((1<<W.keySize)-1);
							lastWord = idx;
							std::cerr << idx << ":" << bits << std::endl;
							break;
						}
					}
				}
				bits-=32;
				out.resize(out.size()+4);
				uint32_t *v = (uint32_t *)&*out.end();
				*--v = value>>bits;
			}
		}
	};
	const Encoder encoder = Encoder(dictionary);
	
	struct Decoder {

		const size_t keySize;     // Non overlapping bits of the word index in the big dictionary
		const size_t overlap;     // Bits that overlap between keys
		const size_t maxWordSize;
		
		size_t start;
		
		std::vector<Symbol> decoderTable;
		template<typename T, size_t N, typename TIN, typename TOUT>
		void decodeA(const TIN &in, TOUT &out) const  {
			
			uint8_t *o = (uint8_t *)&out.front();
			const uint32_t *i = (const uint32_t *)in.data();
	
			uint64_t mask = (1<<(keySize+overlap))-1;
			const std::array<T,N>  *DD = (const std::array<T,N> *)decoderTable.data();
			uint64_t v32 = start; int32_t c=-keySize;
			
			while (c>=0 or i<(const uint32_t *)&*in.end()) {
				
//				std::cerr << ((v32>>c) & mask) << ":" << c << std::endl;
				//endianmess
				if (c<0) {
					v32 = (v32<<32) + *i++;
					c   += 32;
//				std::cerr << ((v32>>c) & mask) << ":" << c << std::endl;
				}
				{				
					
//				std::cerr << ((v32>>c) & mask) << ":" << c << std::endl;
					const uint8_t *&&v = (const uint8_t *)&DD[(v32>>c) & mask];
					c -= keySize;
					for (size_t n=0; n<N; n++)
						*(((T *)o)+n) = *(((const T *)v)+n);
					o += v[N*sizeof(T)-1];
				}
			}
			out.resize(o-(uint8_t *)&out.front());	
		}
		
		template<typename T, typename TIN, typename TOUT>
		void decode12(const TIN &in, TOUT &out) const {
			
			uint8_t *o = (uint8_t *)&out.front();
			const uint8_t *i = (const uint8_t *)in.data();
			
			uint64_t mask = (1<<(keySize+overlap))-1;
			const T *D = (const T *)decoderTable.data();
			uint64_t v64 = start;
			while (i<(const uint8_t *)&*in.end()) {

				v64 <<= 6;
				v64 += (*(const uint64_t *)i) & 0x0000FFFFFFFFFFFFULL;
				i+=6;
				
				{				
					T v = D[(v64>>0 ) & mask];
					*((T *)o) = v;
					o += v >> ((sizeof(T)-1)*8);
				}
				{				
					T v = D[(v64>>12) & mask];
					*((T *)o) = v;
					o += v >> ((sizeof(T)-1)*8);
				}
				{				
					T v = D[(v64>>24) & mask];
					*((T *)o) = v;
					o += v >> ((sizeof(T)-1)*8);
				}
				{				
					T v = D[(v64>>36) & mask];
					*((T *)o) = v;
					o += v >> ((sizeof(T)-1)*8);
				}
			}
		}
	
		Decoder(const Dictionary &dict) :
			keySize(std::log2(dict.size())-overlap),
			overlap(dict.overlap),
			maxWordSize(dict.maxWordSize) {
				
			start = 0;
			while (not dict[start].empty()) 
				start++;
				
			decoderTable.resize(dict.size()*(maxWordSize+1));
			for (size_t i=0; i<dict.size(); i++) {

				Symbol *d = &decoderTable[i*(maxWordSize+1)];
				d[maxWordSize] = dict[i].size();
				assert(dict[i].size()<=maxWordSize-3);
				for (auto c : dict[i])
					*d++ = c;
			}
		}
		
		template<typename TIN, typename TOUT>
		void operator()(const TIN &in, TOUT &out) const {
			
			if (out.size() < in.size()) out.resize(in.size());

			if (keySize==12) {
				switch (maxWordSize+1) {
					case   4: return decode12<uint32_t>(in, out);
					case   8: return decode12<uint64_t>(in, out);
				}
			} 
			switch (maxWordSize+1) {
				case   4: return decodeA<uint32_t, 1>(in, out);
				case   8: return decodeA<uint64_t, 1>(in, out);
				case  16: return decodeA<uint64_t, 2>(in, out);
				case  32: return decodeA<uint64_t, 4>(in, out);
				case  64: return decodeA<uint64_t, 8>(in, out);
				case 128: return decodeA<uint64_t,16>(in, out);
				default: throw std::runtime_error ("unsupported maxWordSize");
			}
		}
	};
	
	struct DecoderSlow {
		
		const Dictionary &W;
		DecoderSlow(const Dictionary &dict) : W(dict) {}
		
/*		void operator()(const std::vector<uint8_t> in, std::vector<uint8_t> & out) const {	
			
			out.clear();
			
			uint64_t value = 0;
			while (not W[value].empty()) 
				value++;
			uint32_t bits = 0;
			
			const uint32_t *it = (const uint32_t *)&*in.begin();
			while (bits>0 or it != (const uint32_t *)&*in.end()) {
				if (bits < W.keySize) {
					bits += 32;
					value = (value<<32) + *it++;
				}
				auto idx = (value>>(bits-W.keySize)) & ((1<<(W.keySize+W.overlap))-1);
//				std::cerr << idx << ":" << bits << std::endl;
				bits -= W.keySize;
				out.insert(out.end(),W[idx].begin(),W[idx].end());
			}
		}*/
		
		template<typename TIN, typename TOUT>
		void operator()(const TIN &in, TOUT &out) const  {

			uint8_t *o = (uint8_t *)&out.front();
			const uint32_t *i = (const uint32_t *)&in.front();
			const uint32_t *iend = i + in.size();
			
			uint64_t value = 0;
			while (not W[value].empty()) 
				value++;
			uint32_t bits = 0;
			
			while (bits>0 or i<iend) {
				if (bits < W.keySize) {
					bits += 32;
					value = (value<<32) + *i++;
				}
				auto idx = (value>>(bits-W.keySize)) & ((1<<(W.keySize+W.overlap))-1);
//				std::cerr << idx << ":" << bits << std::endl;
				bits -= W.keySize;
				for (auto c : W[idx]) *o++ = c;
			}
			out.resize((uint8_t *)o-(uint8_t *)&out.front());
		}
		
	};	
	const DecoderSlow decoder = DecoderSlow(dictionary);
	
	
struct TestTimer {
	timespec c_start, c_end;
	void start() { clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &c_start); };
	void stop () { clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &c_end); };
	double operator()() { return (c_end.tv_sec-c_start.tv_sec) + 1.E-9*(c_end.tv_nsec-c_start.tv_nsec); }
};

public:

	const double efficiency;

	Marlin2018Simple (const std::vector<double> &pdf, size_t keySize, size_t overlap, size_t maxWordSize)
		: dictionary(pdf, keySize, overlap, maxWordSize),
		  efficiency(dictionary.calcEfficiency())  {
			  /*
		auto testData = Distribution::getResiduals(pdf, 1<<24);
		
		auto compressedData = testData;
		
		compressedData.resize(8*testData.size());
		encode(testData,compressedData);
		compressedData.resize(8*testData.size());
		TestTimer tt;
		tt.start();
		encode(testData,compressedData);
		tt.stop();
		
		double shannonLimit = Distribution::entropy(pdf)/std::log2(pdf.size());
		double efficiency = shannonLimit / (compressedData.size()/double(testData.size()));
		std::cerr << testData.size() << " " << compressedData.size() << " " << efficiency <<  " " << (testData.size()/tt()/(1<<20)) << std::endl;

		auto uncompressedData = testData;
		uncompressedData.resize(8*testData.size());
		decode(compressedData,uncompressedData);
		uncompressedData.resize(8*testData.size());
		tt.start();
		decode(compressedData,uncompressedData);
		tt.stop();
		
		std::cerr << testData.size() << " " << uncompressedData.size() << std::endl;
		for (size_t i=0; i<10; i++) std::cerr << int(testData[i]) << " | "; std::cerr << std::endl;
		for (size_t i=0; i<10; i++) std::cerr << int(uncompressedData[i]) << " | "; std::cerr << std::endl;
		
		for (size_t i=0; i<1000 and i<testData.size(); i++) std::cerr << int(testData[i]==uncompressedData[i]); std::cerr << std::endl;
		
		
		
		std::cerr << (testData==uncompressedData) <<  " " << (testData.size()/tt()/(1<<20)) << std::endl;
		*/
	}
		  
	template<typename TIN, typename TOUT>
	void encode(const TIN &in, TOUT &out) const { encoder(in, out); }

	template<typename TIN, typename TOUT>
	void decode(const TIN &in, TOUT &out) const { decoder(in, out); }
};
