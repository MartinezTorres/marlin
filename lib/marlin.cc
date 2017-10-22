#include <marlin.h>
#include <bitstream.h>
#include <cx.h>
#include <distribution.h>

#include <algorithm>
#include <cassert>
#include <memory>
#include <cstring>
#include <vector>
#include <set>


namespace {
	
	// Drop it, only 8 bit, 4096 dictionary, 8 byte entry length.
	
	constexpr const size_t ALPHABET_SIZE = 256;
	constexpr const size_t NUM_WORDS = 4096; // Word 0 is reserved to be size 0
	constexpr const size_t WORD_SIZE = 7;
	typedef uint8_t Symbol;
	
	template<class Key, class Compare = std::less<Key>, class Allocator = std::allocator<Key>>
	struct depq : public std::multiset<Key,Compare,Allocator> {
		
		const Key& min() const { return *this->begin(); }
		const Key& max() const { return *this->rbegin(); }
		void push(const Key& e) { this->insert(e); }
		void removeMin() { this->erase(this->begin()); }
		void removeMax() { this->erase(std::prev(this->end())); }
	};
	
	struct Word : public cx::vector<Symbol,WORD_SIZE> {

		int state = 0; // Target State
		double p = 0.0; // Probability
		
		constexpr bool operator<(const Word& rhs) const {
						
			if (this->size() == WORD_SIZE and rhs.size() != WORD_SIZE ) return true;
			if (this->size() != WORD_SIZE and rhs.size() == WORD_SIZE ) return false;
			return p < rhs.p;
		}
	};
	
	struct SingleDictionaryEncoder {
		
		typedef uint16_t WordIdx;
		typedef uint16_t NodeIdx;
		
		struct Node {
			std::array<NodeIdx,ALPHABET_SIZE> child = {};
			WordIdx code = {};
		};            
		cx::vector<Node, 2*NUM_WORDS> nodes = {};

		// Words[0] is empty
		constexpr SingleDictionaryEncoder(cx::vector<Word,NUM_WORDS> words) {
			
			Node blank = {};
			blank.code = 0;

			for (size_t i=0; i<ALPHABET_SIZE; ++i)
				blank.child[i] = i;

			for (size_t i=0; i<ALPHABET_SIZE; ++i)
				nodes.push_back(blank);
			
			for (size_t i=1; i<words.size(); ++i) {
				
				auto &&w = words[i];
				
				NodeIdx nodeId = w[0];
				
				for (size_t j=1; j<w.size(); ++j) {
					
					auto c = w[j];
					if (nodes[nodeId].child[c] == c) {
						nodes[nodeId].child[c] = nodes.size();
						nodes.push_back(blank);
					}
					nodeId = nodes[nodeId].child[c];
				}
				nodes[nodeId].code = i;
			}
		}
		
		size_t encode(const uint8_t *src, size_t srcSize, uint8_t *dst) const {

			uint8_t *dst0 = dst;
			
			while (srcSize and src[srcSize-1]==0) srcSize--;
			const uint8_t *srcEnd = src + srcSize;
			
			if (srcSize==0) return dst0-dst;

			NodeIdx nodeId = *src++;
			if (srcSize > 16*WORD_SIZE+2) while (src < srcEnd-16*WORD_SIZE+2) {

				NodeIdx oldNodeId1 = 0;
				do { 
					oldNodeId1 = nodeId;
					nodeId = nodes[nodeId].child[*src++]; 
				} while (nodeId > ALPHABET_SIZE - 1);
				
				//printf("NN %ld\n", nodeId);
				
				NodeIdx oldNodeId2 = 0;
				do { 
					oldNodeId2 = nodeId;
					nodeId = nodes[nodeId].child[*src++]; 
				} while (nodeId > ALPHABET_SIZE - 1);
				
				*((uint32_t *)dst) = nodes[oldNodeId1].code + (nodes[oldNodeId2].code << 12);
				dst += 3;
			}
			while (true) {
				// nodeId must be encoded.
				NodeIdx oldNodeId1 = 0;
				do { 
					if (src==srcEnd) {
						*((uint32_t *)dst) = nodes[nodeId].code;
						dst += 3;
						return dst-dst0;
					}					
					oldNodeId1 = nodeId;
					nodeId = nodes[nodeId].child[*src++]; 
				} while (nodeId > ALPHABET_SIZE - 1);
				
				NodeIdx oldNodeId2 = 0;
				do { 

					if (src==srcEnd) {
						*((uint32_t *)dst) = nodes[oldNodeId1].code + (nodes[nodeId].code << 12);
						dst += 3;
						return dst-dst0;
					}

					oldNodeId2 = nodeId;
					nodeId = nodes[nodeId].child[*src++]; 
				} while (nodeId > ALPHABET_SIZE - 1);
				
				*((uint32_t *)dst) = nodes[oldNodeId1].code + (nodes[oldNodeId2].code << 12);
				dst += 3;
			}
			assert(false);
		}
	};
		
	struct SingleDictionaryDecoder {
		
		typedef cx::array<uint8_t, WORD_SIZE+1>  Entry;
		cx::vector<Entry, NUM_WORDS> decoderTable = {};

		constexpr SingleDictionaryDecoder(const cx::vector<Word ,NUM_WORDS> &words) {
			
			decoderTable.resize(words.size());
			for (size_t i=0; i<words.size(); ++i) {
				for (size_t j=0; j<words[i].size(); ++j)
					decoderTable[i][j] = words[i][j];
				
				decoderTable[i].back() = words[i].size();
			}
		}
		
		size_t decode(const uint8_t *src, size_t srcSize, uint8_t *dst, size_t dstSize) const {
			
			memset(dst, 0, dstSize);
			assert(srcSize%3 == 0);
			
			uint8_t *dst0 = dst;		
			
			while (srcSize) {
				
				uint32_t v32=0;
				v32 = *(const uint32_t *)src;
				src+=3;
				srcSize-=3;
				
				{				
					Entry v = decoderTable[(v32>>0 ) & 0xFFF];
					*((Entry *)dst) = v;
					dst += v.back();
				}
				{				
					Entry v = decoderTable[(v32>>12) & 0xFFF];
					*((Entry *)dst) = v;
					dst += v.back();
				}
			}
			memset(dst, 0, dstSize - (dst-dst0));
			return dst-dst0;
		}
	};

	struct Codec {
		virtual ~Codec() {}
		virtual size_t encode(const uint8_t *src, size_t srcSize, uint8_t *dst) const = 0;
		virtual size_t decode(const uint8_t *src, size_t srcSize, uint8_t *dst, size_t dstSize) const = 0;
		virtual double predict(const std::array<uint16_t, ALPHABET_SIZE>&, size_t) const = 0;
	};

	struct MarlinV1 : public Codec {
		
		
		const cx::array<double, ALPHABET_SIZE> P;
		const cx::vector<Word,NUM_WORDS> words;
		const SingleDictionaryEncoder encoder;
		const SingleDictionaryDecoder decoder;
		virtual size_t encode(const uint8_t *src, size_t srcSize, uint8_t *dst                ) const { return encoder.encode(src, srcSize, dst);          }
		virtual size_t decode(const uint8_t *src, size_t srcSize, uint8_t *dst, size_t dstSize) const { return decoder.decode(src, srcSize, dst, dstSize); }
		virtual double predict(const std::array<uint16_t, ALPHABET_SIZE>& hist, size_t) const {
			
			double ret = 0.0;
			for (size_t i=0; i<ALPHABET_SIZE; i++)
				ret += -std::log2(P[i])*hist[i];
		
			return ret/8;
		}
		
		static cx::vector<Word,NUM_WORDS> getWords(const cx::array<double, ALPHABET_SIZE> &P) {
			
			printf("DD %lf\n", Distribution::entropy(P));
				
			cx::vector<Word ,NUM_WORDS> words;
			
			// The formulas on the paper expect the alphabet to consist of symbols with decreasing order of probability.
			struct A {
				double p;
				Symbol symbol;
				constexpr bool operator>(const A& rhs) const { return p > rhs.p; }
			};
			
			// States sorted by state number (most probable state first)
			cx::array<A,ALPHABET_SIZE> a = {}; 
			for (size_t i=0; i<ALPHABET_SIZE; ++i)
				a[i] = { P[i]+1e-100 , Symbol(i) };

			std::sort(a.begin(), a.end(), std::greater<A>());			
			
			
			// Initial State Probability (Ps_i)
			cx::array<double,ALPHABET_SIZE> Ps = {}; 
			Ps.front() = 1.;
			
			// 1st Preprocessing
			// Reverse cumulative symbol probability (divisor of the first term in Eq1)
			cx::array<double,ALPHABET_SIZE> P1 = {};
			P1[ALPHABET_SIZE-1] = a[ALPHABET_SIZE-1].p;
			for (size_t i=ALPHABET_SIZE-1; i; i--)
				P1[i-1] = P1[i] + a[i-1].p;

			for (size_t StateUpdateIterations = 3; StateUpdateIterations; --StateUpdateIterations) {

				// 2nd Preprocessing
				// Parts of eq2 that depend on the state
				cx::array<double,ALPHABET_SIZE> P2 = {};
				{
					double sum = 0.0;
					for (size_t i=0; i<ALPHABET_SIZE; ++i)
						P2[i] = sum = sum + Ps[i]/P1[i];
				}				

				// Build the Dictionary
				depq<Word> depqNodes;
				{
					// DICTIONARY INITIALIZATION
					for (size_t n=0; n<ALPHABET_SIZE; ++n) {
						
						Word node = {};
						node.push_back(n);
						node.state = 0;
						node.p = P2[n] * a[n].p;
						
						depqNodes.push(node);
					}
					
					// GROW THE DICTIONARY
					while (depqNodes.size() < NUM_WORDS-1 and depqNodes.max().size() < WORD_SIZE) {
						
						Word node = depqNodes.max();
						depqNodes.removeMax();
						
						Word newNode = {};
						newNode = node;
						newNode.push_back(node.state);

						newNode.state = 0;
						
						newNode.p = node.p * a[node.state].p / P1[node.state];
						depqNodes.push(newNode);
						
						node.p -= newNode.p;
						node.state++;
						
						// We drop this optimization to ease the encoding process.. almost never happened anyways.
						if (node.state == ALPHABET_SIZE-1) {
							printf(". %ld ", depqNodes.size());
							node.push_back(node.state);
							node.state = 0;
						}
						
						depqNodes.push(node);
						
						/// Test if removing a few improbable words would affect performance... minimally.
						//while (depqNodes.min().p<.02/NUM_WORDS) depqNodes.removeMin();
					}
				}
				
				// Update probability state
				{
					cx::array<cx::array<double,ALPHABET_SIZE>,ALPHABET_SIZE> T = {};
					for (auto &&node : depqNodes) 
						for (size_t state=0, w0 = node.front(); state<=w0; ++state)
							 T[state][node.state] += node.p /(P2[w0] * P1[state]);

					for (size_t MaxwellIt = 3; MaxwellIt; --MaxwellIt) {
						
						cx::array<cx::array<double,ALPHABET_SIZE>,ALPHABET_SIZE> T2 = {};
								
						for (size_t i=0; i<ALPHABET_SIZE; ++i)
							for (size_t j = 0; j<ALPHABET_SIZE; ++j)
								for (size_t k=0; k<ALPHABET_SIZE; ++k)
									T2[i][j] += T[i][k] * T[k][j];
						
						T = T2;
					}
					Ps = T[0];					
				}
				
				cx::array<double,ALPHABET_SIZE> P3 = {};
				{
					double sum = 0.0;
					for (size_t i=0; i<ALPHABET_SIZE; ++i)
						P3[i] = sum = sum + Ps[i]/P1[i];
				}				

				

				double meanLength = 0, sump=0;
				words.clear();
				words.push_back(Word());
				for (auto&& node : depqNodes) {
					Word word = {};
					for (auto&& s : node)
						word.push_back(a[s].symbol);
					words.push_back(word);
					
					double p = node.p/P2[node.front()]*P3[node.front()];
					meanLength += p*node.size();
					sump += p;
				}
				printf("Eff: %lf %lf\n", 12./meanLength, sump);
				
/*				for (size_t i=0; i<ALPHABET_SIZE; ++i) {
					craz += P[i]*symbolLength[i]/symbolP[i];
					craz2 += -std::log2(P[i])*P[i];
				}

				for (size_t i=0; i<10; ++i) {
					std::cout << "i " << P[i] << " " << symbolLength[i] << " " << symbolP[i] <<  " " << -std::log2(P[i]) << std::endl;
				}
					
				std::cout << "Mean length: " << RequiredBits<NUM_WORDS-1>::value / meanLength << " " <<  RequiredBits<NUM_WORDS-1>::value / meanLengthApprox << " " << RequiredBits<NUM_WORDS-1>::value * craz << " " << craz2 << std::endl;*/
			}
			
			
/*			struct Compare {
				constexpr bool operator()(const Word &lhs, const Word &rhs) {
					
					if (lhs.size()==rhs.size()) 
						return lhs.size()<rhs.size();
						
					for (size_t i=0; i<lhs.size(); ++i)
						if (lhs[i]!=rhs[i])
							return lhs[i]<rhs[i];
						
					return false;
				}
			};
			sort(words.begin(), words.end());*/
			
			
			printf("WW %ld\n", words.size());
			return words;
		}	

		MarlinV1(const cx::array<double, ALPHABET_SIZE> &P) :
			P(P), words(getWords(P)),
			encoder(words), decoder(words) {}
	};
	
	struct DoubleSharedDictionaryEncoder {
		
		typedef uint16_t WordIdx;
		typedef uint16_t NodeIdx;
		
		struct Node {
			std::array<NodeIdx,ALPHABET_SIZE> child = {};
			WordIdx code = {};
		};            
		cx::vector<Node, 2*NUM_WORDS> nodes = {};

		// Words[0] is empty
		constexpr DoubleSharedDictionaryEncoder(cx::array<cx::vector<Word,NUM_WORDS>,2>) {
		}
		
		size_t encode(const uint8_t *src, size_t srcSize, uint8_t *dst) const {
		}
	};
		
	struct DoubleSharedDictionaryDecoder {
		
		typedef cx::array<uint8_t, WORD_SIZE+1>  Entry;
		cx::vector<Entry, NUM_WORDS> decoderTable = {};

		constexpr DoubleSharedDictionaryDecoder(cx::array<cx::vector<Word,NUM_WORDS>,2> ) {
		}
		
		size_t decode(const uint8_t *src, size_t srcSize, uint8_t *dst, size_t dstSize) const {
		}
	};

	struct MarlinV2 : public Codec {
		
		
		const cx::array<double, ALPHABET_SIZE> P;
		const cx::array<cx::vector<Word,NUM_WORDS>,2> words;
		const DoubleSharedDictionaryEncoder encoder;
		const DoubleSharedDictionaryDecoder decoder;
		virtual size_t encode(const uint8_t *src, size_t srcSize, uint8_t *dst                ) const { return encoder.encode(src, srcSize, dst);          }
		virtual size_t decode(const uint8_t *src, size_t srcSize, uint8_t *dst, size_t dstSize) const { return decoder.decode(src, srcSize, dst, dstSize); }
		virtual double predict(const std::array<uint16_t, ALPHABET_SIZE>& hist, size_t) const {
			
			double ret = 0.0;
			for (size_t i=0; i<ALPHABET_SIZE; i++)
				ret += -std::log2(P[i])*hist[i];
		
			return ret/8;
		}
		
		static cx::array<cx::vector<Word,NUM_WORDS>,2> getWords(const cx::array<double, ALPHABET_SIZE> &P) {
			
			printf("DD %lf\n", Distribution::entropy(P));
				
			cx::array<cx::vector<Word,NUM_WORDS>,2> words;
			
			// The formulas on the paper expect the alphabet to consist of symbols with decreasing order of probability.
			struct A {
				double p;
				Symbol symbol;
				constexpr bool operator>(const A& rhs) const { return p > rhs.p; }
			};
			
			// States sorted by state number (most probable state first)
			cx::array<A,ALPHABET_SIZE> a = {}; 
			for (size_t i=0; i<ALPHABET_SIZE; ++i)
				a[i] = { P[i]+1e-100 , Symbol(i) };

			std::sort(a.begin(), a.end(), std::greater<A>());			
			
			
			// Initial State Probability (Ps_i)
			cx::array<double,ALPHABET_SIZE> Ps = {}; 
			Ps.front() = 1.;
			
			// 1st Preprocessing
			// Reverse cumulative symbol probability (divisor of the first term in Eq1)
			cx::array<double,ALPHABET_SIZE> P1 = {};
			P1[ALPHABET_SIZE-1] = a[ALPHABET_SIZE-1].p;
			for (size_t i=ALPHABET_SIZE-1; i; i--)
				P1[i-1] = P1[i] + a[i-1].p;

			for (size_t StateUpdateIterations = 3; StateUpdateIterations; --StateUpdateIterations) {
				
				// 2nd Preprocessing
				// Parts of eq2 that depend on the state
				cx::array<double,ALPHABET_SIZE> P2 = {};
				{
					double sum = 0.0;
					for (size_t i=0; i<ALPHABET_SIZE; ++i)
						P2[i] = sum = sum + Ps[i]/P1[i];
				}
				
				depq<Word> depqNodes; // Where 
				// DICTIONARY INITIALIZATION
				for (size_t n=0; n<ALPHABET_SIZE; ++n) {
					
					Word node = {};
					node.push_back(n);
					node.state = 0;
					node.p = P2[n] * a[n].p;
					
					depqNodes.push(node);
				}

				/// EXPAND DICTIONARY FOR EACH STATE
				for (size_t state = 0; state < words.size(); ++state) {
					
					// GROW THE DICTIONARY
					while (depqNodes.size() < NUM_WORDS-1 and depqNodes.max().size() < WORD_SIZE) {
						
						Word node = depqNodes.max();
						depqNodes.removeMax();
						
						Word newNode = {};
						newNode = node;
						newNode.push_back(node.state);

						newNode.state = 0;
						
						newNode.p = node.p * a[node.state].p / P1[node.state];
						depqNodes.push(newNode);
						
						node.p -= newNode.p;
						node.state++;
						
						depqNodes.push(node);
					}
					
					// Init words
					words[state].clear();
					
					// Add words that start with zero
					for (auto&& w : depqNodes)
						if (w[0]==state or state==(words.size()-1))
							words[state].push_back(w);
					
					// Add an empty word
					words[state].push_back(Word());
					
					//Init with all nodes that do not start with zero:
					depq<Word> depqNodes2;
					for (auto&& w : depqNodes)
						if (w[0]>state)
							depqNodes2.push(w);
					depqNodes = depqNodes2;
				}

				// Update probability state
				{
					cx::array<cx::array<double,ALPHABET_SIZE>,ALPHABET_SIZE> T = {};
					
					for (size_t state = 0; state < ALPHABET_SIZE; ++state)
						for (auto &&word : words[std::min(words.size()-1,state)])
							if (state <= word.front())
								T[state][word.state] += word.p /(P2[word.front()] * P1[state]);

					for (size_t MaxwellIt = 3; MaxwellIt; --MaxwellIt) {
						
						cx::array<cx::array<double,ALPHABET_SIZE>,ALPHABET_SIZE> T2 = {};
								
						for (size_t i=0; i<ALPHABET_SIZE; ++i)
							for (size_t j = 0; j<ALPHABET_SIZE; ++j)
								for (size_t k=0; k<ALPHABET_SIZE; ++k)
									T2[i][j] += T[i][k] * T[k][j];
						
						T = T2;
					}
					Ps = T[0];					
				}
				
				// Estimate efficiency
				{
					
					cx::array<double,ALPHABET_SIZE> P3 = {};
					{
						double sum = 0.0;
						for (size_t i=0; i<ALPHABET_SIZE; ++i)
							P3[i] = sum = sum + Ps[i]/P1[i];
					}
				
					double meanLength = 0, sump = 0;
					for (auto && wordStates : words) {
						for (auto && word : wordStates) {
					
							//double p = word.p/P2[word.front()]*P3[word.front()];
							double p = word.p;
							meanLength += p*word.size();
							sump += p;
						}
					}
					
					printf("Eff: %lf %lf\n", 12./meanLength, sump);
				}
			}
			
			// Swap symbols order for actual symbols
			for (auto && wordStates : words)
				for (auto && word : wordStates)
					for (auto&& c : word)
						c = a[c].symbol;
			return words;
		}	

		MarlinV2(const cx::array<double, ALPHABET_SIZE> &P) :
			P(P), words(getWords(P)),
			encoder(words), decoder(words) {}
	};
	



	struct CopyCodec : public Codec {
		
		virtual size_t encode(const uint8_t *src, size_t srcSize, uint8_t *dst) const {
			memcpy(dst, src, srcSize);
			return srcSize;
		}
		virtual size_t decode(const uint8_t *src, size_t srcSize, uint8_t *dst, size_t dstSize) const {
			assert(srcSize == dstSize);
			memcpy(dst, src, srcSize);
			return srcSize;
		}
		virtual double predict(const std::array<uint16_t, ALPHABET_SIZE>&, size_t sz) const { return sz; }
	};

	struct SkipCodec : public Codec {
		
		virtual size_t encode(const uint8_t *, size_t, uint8_t *) const { return 0; }
		virtual size_t decode(const uint8_t *, size_t, uint8_t *dst, size_t dstSize) const {
			memset(dst,0,dstSize);
			return dstSize;
		}
		virtual double predict(const std::array<uint16_t, ALPHABET_SIZE>& hist, size_t sz) const {
			return hist[0]==sz?0.0:1E100;
		}
	};
	
	template<typename T>
	static inline std::vector<std::shared_ptr<Codec>> getDictionaries() {
		
		std::vector<std::shared_ptr<Codec>> ret;
		ret.emplace_back( std::make_shared<CopyCodec>() );
		ret.emplace_back( std::make_shared<SkipCodec>() );
		
		//return ret;
		double e = 0.9;
		for (size_t i=0; i<6; i++) {
			printf("G %lf", e); ret.emplace_back( std::make_shared<T>(Distribution::getWithEntropy(Distribution::Gaussian<256>,e) ) );
			printf("G2 %lf", e); MarlinV2(Distribution::getWithEntropy(Distribution::Gaussian<256>,e) );

			printf("L %lf", e); ret.emplace_back( std::make_shared<T>(Distribution::getWithEntropy(Distribution::Laplace<256>,e) ) );
			printf("E %lf", e); ret.emplace_back( std::make_shared<T>(Distribution::getWithEntropy(Distribution::Exponential<256>,e) ) );
			printf("P %lf", e); ret.emplace_back( std::make_shared<T>(Distribution::getWithEntropy(Distribution::Poisson<256>,e) ) );
			e*=0.8;
		}
		return ret;
	}
		
	const auto dictionariesMarlin = getDictionaries<MarlinV1>();
	
}

constexpr const static size_t BS = 4096;
constexpr const static size_t NB = 16;

// Streaming is important...
// In compression: 
// If dstCapacity > srcSize -> rearrange and send
// If dstCapacity < srcSize -> do not rearrange

// Always ENCODE blocks in order anyways.
// Header: isBlockSize(1) isPredicted(1) size(14) {if isBlockSize, uncompSize(16)} dict(8) data(size)


size_t MarlinEncode(const uint8_t*&& src, size_t srcSize, uint8_t* dst, size_t dstCapacity) {
	
	
	struct Block {
		
		cx::array<uint8_t,BS  > bufferFiltered;
		cx::array<uint8_t,BS*4> buffer;
		const uint8_t *src;
		size_t sizeUncompressed;

		size_t sizeCompressed;
		
		size_t dictionaryIdx;
		bool isFiltered;
		bool isCompressed;
	};
	
	cx::array<Block, NB> blocks = {};
	size_t start=0, end=0;
		
	const uint8_t* dst0 = dst;
	
	while (srcSize or start!=end) {
		
		// Prepare
		while (srcSize and end-start<NB) {
			
			Block &block = blocks[end % NB]; 
			end++;
			
			block.src = src;
			block.sizeUncompressed = std::min(srcSize,BS);

			block.dictionaryIdx = 0;
			block.isFiltered = 0;
			block.isCompressed = 0;
			
			if (block.sizeUncompressed > 64) {

				double expectedSize = block.sizeUncompressed;
				
				std::array<uint16_t, 256> hist = {};
				for (size_t i=0; i<block.sizeUncompressed; i++) 
					hist[src[i]]++;
				
				// We first try to compress it without prediction
				for (size_t i=0; i<dictionariesMarlin.size(); i++) {
					
					double expected = dictionariesMarlin[i]->predict(hist, block.sizeUncompressed);
					//printf("%ld %lf\n", i, expected);
					
					if (expected < expectedSize) {
						expectedSize = expected;
						block.dictionaryIdx = i;
					}
				}
				
			
				// Now with filtering
				hist.fill(0);
				for (size_t i=16; i<block.sizeUncompressed; i++) 
					hist[uint8_t(int8_t(src[i])-int8_t(src[i-16]))]++;

				for (size_t i=0; i<dictionariesMarlin.size(); i++) {
					
					double expected = 16 + dictionariesMarlin[i]->predict(hist, block.sizeUncompressed);
					//printf("%ld %lf\n", i, expected);
						
					if (expected < expectedSize) {
						expectedSize = expected;
						block.dictionaryIdx = i;
						block.isFiltered = true;
					}				
				}
				//printf("KK %ld %lf\n", block.dictionaryIdx, expectedSize);
				
			}
			
			srcSize -= block.sizeUncompressed;
			src += block.sizeUncompressed;
		}
		
		// Encode		
		for (size_t i=start; i<end; i++) {
			
			if (blocks[i%NB].dictionaryIdx != blocks[start%NB].dictionaryIdx) 
				continue;
				
			Block &block = blocks[i % NB]; 
			
			if (block.isFiltered) {
				
				memcpy(&block.buffer[0], block.src, 16);
				for (size_t i=0; i<block.sizeUncompressed-16; i++)
					block.bufferFiltered[i] = uint8_t(int8_t(src[i+16])-int8_t(src[i]));
					
				block.sizeCompressed = 16 + dictionariesMarlin[block.dictionaryIdx]->encode(block.src+16, block.sizeUncompressed-16, &block.buffer[16]);
			} else {
				
				block.sizeCompressed =      dictionariesMarlin[block.dictionaryIdx]->encode(block.src, block.sizeUncompressed, &block.buffer[0]);
			}
			
			printf("KK %ld %ld\n", block.dictionaryIdx, block.sizeCompressed);
			
			if (block.sizeCompressed > 0.98*block.sizeUncompressed) {

				block.dictionaryIdx = 0;
				block.isFiltered = false;

				block.sizeCompressed =      dictionariesMarlin[block.dictionaryIdx]->encode(block.src, block.sizeUncompressed, &block.buffer[0]);
			}
			
			block.isCompressed = true;
		}
		
		//Write
		while (start<end and blocks[start%NB].isCompressed) {
			
			Block &block = blocks[start%NB]; start++;
			if (dstCapacity < block.sizeCompressed+3) {
				src = block.src;
				return dst-dst0;
			}
			
			if (block.sizeUncompressed==BS) { // Normal header
				*((uint16_t *)dst) =  0x8000*0 + 0x4000*block.isFiltered + block.sizeCompressed;
				dst += sizeof(uint16_t);
				*((uint8_t *)dst) =  block.dictionaryIdx;
				dst += sizeof(uint8_t);
				dstCapacity -= sizeof(uint16_t) + sizeof(uint8_t);
			} else {
				*((uint16_t *)dst) =  0x8000*1 + 0x4000*block.isFiltered + block.sizeCompressed;
				dst += sizeof(uint16_t);
				*((uint8_t *)dst) =  block.dictionaryIdx;
				dst += sizeof(uint8_t);
				*((uint16_t *)dst) =  block.sizeUncompressed;
				dst += sizeof(uint16_t);
				dstCapacity -= sizeof(uint16_t) + sizeof(uint16_t)+ sizeof(uint8_t);
			}
			
			std::memcpy(dst, block.buffer.data(), block.sizeCompressed);

			dstCapacity -= block.sizeCompressed;
			dst += block.sizeCompressed;
		}
	}	
	return dst-dst0;
}


size_t MarlinDecode(const uint8_t*&& src, size_t srcSize, uint8_t* dst, size_t dstCapacity) {
	
	
	struct Block {
		
		const uint8_t *src;
		size_t sizeCompressed;

		uint8_t *dst;
		size_t sizeUncompressed;
		
		
		size_t dictionaryIdx;
		bool isFiltered;
		bool isProcessed;
	};
	
	cx::array<Block, NB> blocks = {};
	size_t start=0, end=0;
	
	const uint8_t *dst0 = dst, *dstEnd = dst+dstCapacity;
	
	while (srcSize > 5 or start!=end) {
		
		// Prepare
		while (srcSize and end-start<NB) {
			
			Block &block = blocks[end % NB]; 
			
			block.src = src;
			
			bool isShort = !!(*((const uint16_t *)block.src) & 0x8000);
			block.isFiltered = !!(*((const uint16_t *)block.src) & 0x4000);
			block.sizeCompressed = *((const uint16_t *)block.src) & 0x3FFF;
			block.src += sizeof(uint16_t);
			block.dictionaryIdx = *((const uint8_t *)block.src);
			block.src += sizeof(uint8_t);

			if (srcSize < block.sizeCompressed + 3 + (isShort?2:0)) break;

			if (isShort) {
				block.sizeUncompressed = *((const uint16_t *)block.src);
				block.src += sizeof(uint16_t);
			} else {
				block.sizeUncompressed = BS;
			}
			
			block.dst = dst;
			dst += block.sizeUncompressed;
			
			block.isProcessed = false;
			
			srcSize -= block.sizeCompressed + 3 + (isShort?2:0);
			src     += block.sizeCompressed + 3 + (isShort?2:0);
			end++;
		}
		
		// Decode and Write		
		for (size_t i=start; i<end; i++) {
			
			if (blocks[i%NB].dictionaryIdx != blocks[start%NB].dictionaryIdx) 
				continue;
				
			Block &block = blocks[i % NB]; 
			
			if (block.isProcessed) continue;
			
			if (block.dst + block.sizeUncompressed > dstEnd) {
				if (i==start) {
					return block.dst - dst0;
				} else {
					break;
				}
			}

			std::memset(block.dst, 0, block.sizeUncompressed);
			if (block.isFiltered) {
				
				printf("F");
				memcpy(block.dst, block.src, 16);

				dictionariesMarlin[block.dictionaryIdx]->decode(block.src+16, block.sizeCompressed-16, block.dst+16, block.sizeUncompressed-16);
				
				for (size_t i=16; i<block.sizeUncompressed; i++)
					block.dst[i] = uint8_t(int8_t(block.dst[i])+int8_t(block.dst[i-16]));
				
			} else {
				printf("N");

				dictionariesMarlin[block.dictionaryIdx]->decode(block.src, block.sizeCompressed, block.dst, block.sizeUncompressed);
				
			}
			
			block.isProcessed = true;
			
			if (i==start) start++;
		}
	}
	return dst - dst0;
}


