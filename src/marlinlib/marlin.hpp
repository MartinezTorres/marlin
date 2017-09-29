#pragma once
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
				if (isVictim or node->p>Marlin2018Simple::purgeProbabilityThreshold) {
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
			
			// Include empty?
			if (isVictim) {
				retiredNodes++;
			} else {
				pq.push(root);
			}
			root->erased = true;

			for (size_t c=0; c<alphabet.size(); c++) {			
					
				root->push_back(std::make_shared<Node>());
				double sum = 0;
				for (size_t t = 0; t<=c; t++) sum += Pstates[t]/PN[t];
				root->back()->p = sum * alphabet[c].p;
				root->back()->sz = 1;
				pushAndPrune(root->back());
			}
				
			// DICTIONARY GROWING
			while (pq.size() + retiredNodes < (1U<<keySize)) {
					
				auto node = pq.top();
				pq.pop();
					
				double p = node->p * Pchild[node->size()];
				node->push_back(std::make_shared<Node>());
				node->back()->p = p;
				root->back()->sz = node->sz+1;
				node->p -= p;
				pushAndPrune(node->back());
					
				if (node->size()<alphabet.size()-1) {

					pushAndPrune(node);
						
				} else {
					node->erased = true;
					node->p = 0;

					node->push_back(std::make_shared<Node>());
					node->back()->p = node->p;
					root->back()->sz = node->sz+1;
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
			for (auto &&node : nodes) {

				std::vector<Word> sortedDictionary = buildWords(node);
				auto cmp = [](const Word &lhs, const Word &rhs) { 
					if (lhs.state != rhs.state) return lhs.state<rhs.state;
					return lhs.p > rhs.p;
				};
				std::sort(sortedDictionary.begin(), sortedDictionary.end(), cmp);
				
				std::vector<Word> w(sortedDictionary.size());
				for (size_t i=0,j=0,k=0; i<sortedDictionary.size(); j+=(1<<overlap)) {
					
					if (j>=ret.size()) 
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
				
			double efficiency = shannonLimit / keySize / (meanLength*std::log2(P.size()));

			//printf("Compress Ratio: %3.4lf\n", (std::log2(dictSize)-overlapping)/(meanLength*std::log2(P.size())));
			
			
			return efficiency;
		}
		
//		Dictionary(const std::map<Symbol, double> &symbols, size_t dictSize, size_t overlap, size_t maxWordSize)
		Dictionary(const Alphabet &alphabet, size_t keySize, size_t overlap, size_t maxWordSize)
			: alphabet(alphabet), keySize(keySize), overlap(overlap), maxWordSize(maxWordSize) {
			
			std::vector<std::vector<double>> Pstates;
			for (auto k=0; k<(1<<overlap); k++) {
				std::vector<double> PstatesSingle(alphabet.size(), 0.);
				PstatesSingle[0] = 1./(1<<overlap);
				Pstates.push_back(PstatesSingle);
			}
			
			int leastProbable = 0;
			
			std::vector<std::shared_ptr<Node>> dictionaries;
			for (auto k=0; k<(1<<overlap); k++)
				dictionaries.push_back(buildTree(Pstates[k], k==leastProbable) );
				
			std::vector<Word> words = arrangeAndFuse(dictionaries,leastProbable);
				
			//if (Marlin2018Simple::debug) print(dictionary);
				
			for (size_t iteration=0; iteration<Marlin2018Simple::iterationLimit; iteration++) {

				// UPDATING STATE PROBABILITIES
				{
					for (auto k=0; k<(1<<overlap); k++)
						Pstates[k] = std::vector<double>(alphabet.size(), 0.);

					for (size_t i=0; i<words.size(); i++)
						Pstates[i%(1<<overlap)][words[i].state] += words[i].p;
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
						leastProbable = i;
					}
				}
				//if (Marlin2018Simple::debug) print(Pstates);

				dictionaries.clear();
				for (auto k=0; k<(1<<overlap); k++)
					dictionaries.push_back(buildTree(Pstates[k], k==leastProbable) );
				
				words = arrangeAndFuse(dictionaries,leastProbable);
				
				//if (Marlin2018Simple::debug) print(dictionary);
				//if (Marlin2018Simple::debug) printf("Efficiency: %3.4lf\n", efficiency);
			}
			this->insert(this->end(), words.begin(), words.end());
			//if (Marlin2018Simple::debug) printf("Efficiency: %3.4lf\n", efficiency);				
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

			const size_t wordStride;  // Bit stride of the jump table corresponding to the word dimension
			size_t nextIntermediatePos = 1<<(wordStride-1);	
			std::vector<JumpIdx> table;
		
		public:
		
			JumpTable(size_t keySize, size_t overlap, size_t nAlpha) :
				wordStride(keySize+overlap+1), // Extra bit for intermediate nodes.
				table(wordStride*nAlpha,JumpIdx(-1))
				{}
			
			JumpIdx &operator()(const size_t &word, const size_t &nextLetter) { 
				return table[(word&(1<<wordStride))+(nextLetter<<wordStride)];
			}
			
			JumpIdx  operator()(const size_t &word, const size_t &nextLetter) const { 
				return table[(word&(1<<wordStride))+(nextLetter<<wordStride)];
			}
			
			size_t getNewPos() { return nextIntermediatePos++; }
			
			bool isIntermediate(size_t pos) const { return pos & (1<<(wordStride-1)); }
		};
		JumpTable jumpTable;
		
		std::vector<JumpIdx> start;
		std::vector<JumpIdx> emptyWords;
		const size_t keySize;

		Encoder(const Dictionary &dict) :
			jumpTable(dict.keySize, dict.overlap, dict.alphabet.size()),
			keySize(dict.keySize) { 
			
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
					while (not word.empty()) {
						auto lastSymbol = word.back();						
						word.pop_back();
						size_t parentIdx;
						if (positions[k].count(word)) {
							parentIdx = positions[k][word];
						} else {
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
							if (positions[i>>dict.keySize].count(Word(1,Symbol(j)))) {
								jumpTable(i, j) = positions[i>>dict.keySize][Word(1,Symbol(j))] +
									FLAG_NEXT_WORD;
							} else { 
								jumpTable(i, j) = positions[i>>dict.keySize][Word(1,Symbol(j))] +
									FLAG_NEXT_WORD + 
									FLAG_INSERT_EMPTY_WORD;
							}
						}
					}
				}
			}
			
			// Fill list of empty words
			emptyWords.resize(dict.alphabet.size(),JumpIdx(-1));
			for (size_t j=0; j<dict.alphabet.size(); j++)
				if (positions[j].count(Word()))
				emptyWords[j] = positions[j][Word()];
			
			// Get Starting Positions
			start.resize(dict.alphabet.size(),JumpIdx(-1));
			for (size_t j=0; j<dict.alphabet.size(); j++) {
				for (size_t k=0; k<NumSections; k++) {
					if (positions[k].count(Word(1,Symbol(j)))) {
						start[j] = positions[k][Word(1,Symbol(j))];
						break;
					}
				}
			}
		}
		
		template<class TIN, typename TOUT, typename std::enable_if<sizeof(typename TIN::value_type)==1,int>::type = 0>		
		void operator()(const TIN &in, TOUT &out) const {
			
			if (out.size() < in.size()) out.resize(in.size());
			
			uint8_t *o = (uint8_t *)&out.front();
			const uint8_t *i = (const uint8_t *)&in.front();
			const uint8_t *iend = i + in.size();
			
			uint64_t value=0; int32_t bits=0;
			if (i<iend) {
				
				JumpIdx j0 = start[*i++];					
				while (i<iend) {
					
					JumpIdx j1 = jumpTable(j0, *i++);

					if (j1 & FLAG_NEXT_WORD) {
						value <<= keySize;
						bits += keySize;
						value += j0 & ((1<<keySize)-1);
						if (bits>32) {
							bits -= 32;
							*((uint32_t *)o) = uint32_t(value>>bits);
							o+=4;
						}

						if (j1 & FLAG_INSERT_EMPTY_WORD)  i--;
					}
					j0=j1;
				}
				assert (not jumpTable.isIntermediate(j0)); //If we end in an intermediate node, we should roll back. Not implemented.
				value <<= keySize;
				bits += keySize;
				value += j0 & ((1<<keySize)-1);
				if (bits) {
					*o++ = uint8_t((value<<8)>>bits);
					bits -=8;
				}
			}
		}
	};
	const Encoder encoder;
	
	struct Decoder { //Only valid for dictionaries with sizes multiple of 2

		const size_t keySize;     // Non overlapping bits of the word index in the big dictionary
		const size_t overlap;     // Bits that overlap between keys
		const size_t maxWordSize;
		
		std::vector<uint8_t> decoderTable;
		template<typename T, size_t N, typename TIN, typename TOUT>
		void decodeA(const TIN &in, TOUT &out) const  {
			
			uint8_t *o = (uint8_t *)out.data();
			const uint32_t *i = (const uint32_t *)in.data();
	
			uint64_t mask = (1<<(keySize+overlap))-1;
			const std::array<T,N>  *DD = (const std::array<T,N> *)decoderTable.data();
			uint64_t v32=0; uint32_t c=0;
			while (i<(const uint32_t *)in.end()) {
				
				if (c<32) {
					v32 += uint64_t(*i++) << c;
					c   += 32;
				}
				{				
					const uint8_t *&&v = (const uint8_t *)&DD[v32 & mask];
					v32 >>= keySize;
					c -= keySize;
					for (size_t n=0; n<N; n++)
						*(((T *)o)+n) = *(((const T *)v)+n);
					o += v[N*sizeof(T)-1];
				}
				{				
					const uint8_t *&&v = (const uint8_t *)&DD[v32 & mask];
					v32 >>= keySize;
					c -= keySize;
					for (size_t n=0; n<N; n++)
						*(((T *)o)+n) = *(((const T *)v)+n);
					o += v[N*sizeof(T)-1];
				}
			}				
		}
		
		template<typename T, typename TIN, typename TOUT>
		void decode12(const TIN &in, TOUT &out) const {
			
			uint8_t *o = (uint8_t *)out.data();
			const uint8_t *i = (const uint8_t *)in.data();
			
			uint64_t mask = (1<<(keySize+overlap))-1;
			const T *D = (const T *)decoderTable.data();
			uint64_t v64=0;
			while (i<(const uint8_t *)in.end()) {

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
				
			decoderTable.resize(dict.size()*maxWordSize);
			for (size_t i=0; i<dict.size(); i++) {

				uint8_t *d = &decoderTable[i*maxWordSize];
				d[maxWordSize] = dict[i].size();
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
	const Decoder decoder;	
	
public:

	const double efficiency;

	Marlin2018Simple (const Alphabet &alphabet, size_t keySize, size_t overlap, size_t maxWordSize)
		: dictionary(alphabet, keySize, overlap, maxWordSize),
		  encoder(dictionary),
		  decoder(dictionary),
		  efficiency(dictionary.calcEfficiency())  {}
		  
	template<typename TIN, typename TOUT>
	void encode(const TIN &in, TOUT &out) const { encoder(in, out); }

	template<typename TIN, typename TOUT>
	void decode(const TIN &in, TOUT &out) const { decoder(in, out); }
};
