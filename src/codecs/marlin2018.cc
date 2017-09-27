#include <codecs/marlin2018.hpp>
#include <util/distribution.hpp>

#include <functional>
#include <string>
#include <iostream>
#include <cassert>
#include <cstring>
#include <stack>
#include <queue>
#include <map>
#include <bitset>
#include <unordered_map>
#include <algorithm>
#include <memory>

struct Marlin2018Pimpl : public CODEC8Z {
	
	class Marlin2018Simple {
	public:

		// Configuration		
		static const constexpr bool enableVictimDictionary = true;
		static const constexpr double purgeProbabilityThreshold = 1e-10;
		static const constexpr size_t iterationLimit = 3;
		static const constexpr bool debug = false;

		typedef uint8_t Symbol; // storage used to store an input symbol.
		typedef uint16_t WordIdx; // storage that suffices to store a word index.

		// Member Types and Variables

		struct Word : std::vector<Symbol> {
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
			
//				auto cmp = [](const std::shared_ptr<Node> &lhs, const std::shared_ptr<Node> &rhs)
//					return lhs->p*(1+std::pow(lhs->sz,1)) < rhs->p*(1+std::pow(rhs->sz,1));	};

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
				while (pq.size() + retiredNodes < dictSize) {
						
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
			
			std::vector<Word> buildWords(const std::shared_ptr<Node> &root) const {
			
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
			
			std::vector<Word> arrangeAndFuse( const std::vector<std::shared_ptr<Node>> &nodes, size_t victimIndex ) const {
				
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
							
						if (victimIdx==int(j)) {
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

			const size_t dictSize; // Size of the big dictionary in words
			const size_t overlap; // Overlap, in bytes
			const size_t maxWordSize; // Maximum number of symbols that a word in the dictionary can have.

			struct SymbolWithP {
				Symbol symbol;
				double p;
				bool operator<(const SymbolWithP &rhs) {
					if (p!=rhs.p) return p>rhs.p; // Descending in probability
					return symbol<rhs.symbol; // Ascending in symbol index
				}
			};
					
			struct Alphabet : std::vector<SymbolWithP> {
				Alphabet(const std::map<Symbol, double> &symbols) {
					for (auto &&symbol : symbols)
						this->push_back(SymbolWithP({symbol.first, symbol.second}));
					std::sort(this->begin(),this->end());
				}
			};
			const Alphabet alphabet;

			double calcEfficiency() const {
			
				double meanLength = 0;
				for (auto &&w : *this)
					meanLength += w.p * w.size();
				
				std::vector<double> P;
				for (auto &&a: alphabet) P.push_back(a.p);
				double shannonLimit = Distribution::entropy(P)/std::log2(P.size());
					
				double efficiency = shannonLimit / ((std::log2(dictSize)-overlap)/(meanLength*std::log2(P.size())));

				//printf("Compress Ratio: %3.4lf\n", (std::log2(dictSize)-overlapping)/(meanLength*std::log2(P.size())));
				
				
				return efficiency;
			}
			
			Dictionary(const std::map<Symbol, double> &symbols, size_t dictSize, size_t overlap, size_t maxWordSize)
				: dictSize(dictSize), overlap(overlap), maxWordSize(maxWordSize),
				  alphabet(symbols) {
				
				std::vector<std::vector<double>> Pstates;
				for (auto k=0; k<(1<<overlap); k++) {
					std::vector<double> PstatesSingle(alphabet.size(), 0.);
					PstatesSingle[0] = 1./(1<<overlap);
					Pstates.push_back(PstatesSingle);
				}
				
				int leastProbable = 0;
				
				std::vector<std::vector<Word>> dictionaries;
				for (auto k=0; k<(1<<overlap); k++)
					dictionaries.push_back( arrangeWords( buildWords( buildTree(Pstates[k], k!=leastProbable) ), (k==leastProbable?-1:k) ) );
					
				std::vector<Word> words = concatenate(dictionaries);
					
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

					for (auto k=0; k<(1<<overlap); k++)
						dictionaries.push_back( arrangeWords( buildWords( buildTree(Pstates[k], k!=leastProbable) ), (k==leastProbable?-1:k) ) );
					
					words = concatenate(dictionaries);
					
					//if (Marlin2018Simple::debug) print(dictionary);
					//if (Marlin2018Simple::debug) printf("Efficiency: %3.4lf\n", efficiency);
				}
				
				//if (Marlin2018Simple::debug) printf("Efficiency: %3.4lf\n", efficiency);				
			}			
		};
		const Dictionary dictionary;
		
		struct Encoder {

			typedef uint32_t JumpIdx;
			size_t keySize,overlap;
			size_t wordStride,alphaStride;
			std::vector<JumpIdx> jumpTable;
			std::vector<JumpIdx> start;

			Encoder(const Dictionary &dictionary) {
				
				overlap(dictionary.overlap);
				keySize = 0;
				while ((1<<keySize)<dictionary.size()) keySize++;
				wordStride = keySize+1; // Extra bit for intermediate nodes.
				
				alphaStride = 0;
				while ((1<<alphaStride)<dictionary.alphabet.size()) alphaStride++;
				
				jumpTable.resize((1<<wordStride)*(1<<alphaStride),JumpIdx(-1));
				
				std::vector<std::map<Word, size_t>> positions(nDict);
				size_t nDict = 1<<dictionary.overlap;				
				// Init the mapping (to know where each word goes)
				for (size_t k=0; k<nDict; k++)				
					for (size_t i=k*(dictionary.size()/nDict); i<(k+1)*(dictionary.size()/nDict); i++)
						positions[k][dictionary[i]] = i;
				
				// Link each possible word to its continuation
				size_t nextIntermediatePos = 1<<(wordStride-1);
				for (size_t k=0; k<nDict; k++) {
					for (size_t i=k*(dictionary.size()/nDict); i<(k+1)*(dictionary.size()/nDict); i++) {
						Word parent = dictionary[i];
						size_t pos = i;
						while (not dictionary[i].empty()) {
							auto lastSymbol = parent.back();						
							parent.pop_back();
							size_t newPos;
							if (positions[k].count(parent)) {
								newPos = positions[k][parent];
							} else {
								newPos = nextIntermediatePos++;
							}
							jumpTable[newPos+(lastSymbol<<wordStride)] = pos;
							pos = newPos;
						}
					}
				}
				
				//Link between inner dictionaries
				for (size_t k=0; k<nDict; k++) {
					for (size_t i=k*(dictionary.size()/nDict); i<(k+1)*(dictionary.size()/nDict); i++) {
						Word &&w = dictionary[i];
						for (size_t j=0; j<(1<<alphaStride); j++) {
							if (jumpTable[i+(j<<wordStride)]==JumpIdx(-1)) {
								if (positions[w.state].count(Word(1,Symbol(j)))) {
									jumpTable[i+(j<<wordStride)] = 
										positions[w.state][Word(1,Symbol(j))] +
										(1<<(2*wordStride)) // Marker to advance step;
								} else if (positions[w.state].count(Word())) {
									jumpTable[i+(j<<wordStride)] =
										 positions[w.state][Word(1,Symbol(j))] +
										(positions[w.state][Word()] << wordStride) +
										(2<<(2*wordStride)); // Marker to advance step twice.
								}
							}
						}
					}
				}
				
				// Get Starting Positions
				start.resize(1<<alphaStride,JumpIdx(-1));
				if (positions[0].count(Word())) { // 0 is not victim
					for (size_t j=0; j<(1<<alphaStride); j++)
						start[j] = jumpTable[positions[0][positions[0][Word()]]+(j<<wordStride)];
				} else { // 0 is victim and must have representation of all words of a single letter
					for (size_t j=0; j<(1<<alphaStride); j++)
						if (positions[0].count(Word(1,Symbol(j))))
							start[j] = positions[0][Word(1,Symbol(j))];
				}
			}
			
			template<typename TIN, typename TOUT>
			void encode(const TIN &in, TOUT &out) const {
				
				if (out.size() < in.size()) out.resize(in.size());
				
				uint8_t *o = (uint8_t *)out.data();
				const uint8_t *i = (const uint8_t *)in.data();
				
				uint64_t mask = (1<<(keySize-overlap))-1;
				uint64_t v=0; uint32_t c=0;
				if (i<(const uint8_t *)in.end()) {
					
					JumpIdx j = start[*i++];					
					while (i<(const uint8_t *)in.end()) {
						
						
					} 
				// TODO: Take care of rolling back. Do it using a recursive function.
				
				
				
			}
		};
		const Encoder encoder;
		
		struct Decoder {

	// Decoder functions
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
	
				Decoder() {
				decoderTable.resize(words.size()*maxWordSize);
			for (size_t i=0; i<words.size(); i++) {

				uint8_t *d = &decoderTable[i*maxWordSize];
				d[maxWordSize] = words[i].symbols.size();
				for (auto c : words[i].symbols)
					*d++ = distSorted[c].second;
			}
			}
			template<typename TIN, typename TOUT>
			void decode(const TIN &in, TOUT &out) const {
				
				
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

		Marlin2018Simple (size_t dictSize, size_t overlap, size_t maxWordSize, const std::map<Symbol, double> &symbols)
			: dictSize(dictSize), overlap(overlap), maxWordSize(maxWordSize),
			  alphabet(symbols),
			  dictionary(),
			  efficiency(calcEfficiency()),
			  encoder(dictionary,overlap),
			  decoder(dictionary,maxWordSize,overlap)  {}
			  
		template<typename TIN, typename TOUT>
		void encode(const TIN &in, TOUT &out) const { encoder.encode(in, out); }
	
		template<typename TIN, typename TOUT>
		void decode(const TIN &in, TOUT &out) const { decoder.decode(in, out); }
	};


	std::vector<std::shared_ptr<Marlin2018Simple>> dictionaries;
	
	std::string coderName;
	std::string name() const { return coderName; }
	
	Marlin2018Pimpl(Distribution::Type distType, size_t keySize, size_t overlap, size_t numDict) {

		{
			std::ostringstream oss;
			oss << "Marlin2018 " << (distType==Distribution::Laplace?"Lap:":"Exp:") <<  ":" << keySize << ":" << overlap << ":" << numDict;
			coderName = oss.str();
		}

		std::vector<std::shared_ptr<Dictionary>> builtDictionaries(numDict);

		#pragma omp parallel for
		for (size_t p=0; p<numDict; p++) {
			
			std::array<double,256> dist; dist.fill(0.);
			for (double i=0.05; i<1; i+=0.1) {
				
				std::array<double,256> pdf = Distribution::pdf(distType, (p+i)/numDict);
				for (size_t j=0; j<256; j++)
					dist[j] += pdf[j]/9.;
			}
			
			std::vector<std::pair<double, uint8_t>> distSorted; distSorted.resize(256); 
			for (size_t i=0; i<256; i++) distSorted[i] = std::make_pair(dist[i],i);
			std::sort   (distSorted.begin(), distSorted.end());
			std::reverse(distSorted.begin(), distSorted.end());
			
			std::vector<double> P;
			for (auto &ds : distSorted)
				P.emplace_back(ds.first);

			builtDictionaries[p] = std::make_shared<Dictionary>(P, keySize, overlap, 4, distSorted);
			for (int maxWordLength=8; maxWordLength <= 128; maxWordLength*=2) {

				std::shared_ptr<Dictionary> dict = std::make_shared<Dictionary>(P, keySize, overlap, maxWordLength, distSorted);
				if (dict->efficiency > builtDictionaries[p]->efficiency+0.005)
					builtDictionaries[p] = dict;
			}
		}
		
		dictionaries.resize(256);
		
		#pragma omp parallel for
		for (size_t h=0; h<256; h+=4) {
			
			auto testData = Distribution::getResiduals(Distribution::pdf(distType, (h+2)/256.), 1<<18);
			
			double lowestSize = testData.size()*0.99; // If efficiency is not enough to compress 1%, skip compression
			for (auto &&dict : builtDictionaries) {
				std::string out;
				dict->encode(testData, out);
				if (out.size() < lowestSize) {
					lowestSize = out.size();
					for (size_t hh = 0; hh<4; hh++)
						dictionaries[h+hh] = dict;
				}
			}	
		}
	}

	
	void   compress(
		const std::vector<std::reference_wrapper<const AlignedArray8>> &in,
		      std::vector<std::reference_wrapper<      AlignedArray8>> &out,
		      std::vector<std::reference_wrapper<      uint8_t      >> &entropy) const { 
		
		for (size_t i=0; i<in.size(); i++)
			if (dictionaries[entropy[i]])
				dictionaries[entropy[i]]->encode(in[i].get(), out[i].get());
			else
				out[i].get().resize(in[i].get().size());
	}

	void uncompress(
		const std::vector<std::reference_wrapper<const AlignedArray8>> &in,
		      std::vector<std::reference_wrapper<      AlignedArray8>> &out,
		      std::vector<std::reference_wrapper<const uint8_t      >> &entropy) const {
		
		for (size_t i=0; i<in.size(); i++)
			if (dictionaries[entropy[i]])
				dictionaries[entropy[i]]->decode(in[i].get(), out[i].get());
			else
				out[i].get().resize(in[i].get().size());
	}

};


Marlin2018::Marlin2018(Distribution::Type distType, size_t keySize, size_t overlap, size_t numDict) 
	: CODEC8withPimpl( new Marlin2018Pimpl(distType, keySize, overlap, numDict) ) {}

