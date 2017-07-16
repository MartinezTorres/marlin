#include <marlin.h>
#include <util.h>
#include <distribution.h>

#include <set>


namespace {
	
	template<uint64_t N> struct RequiredBits    { enum { value = 1 + RequiredBits<(N>>1)>::value }; };
	template<>           struct RequiredBits<0> { enum { value = 0 }; };
	template<>           struct RequiredBits<1> { enum { value = 1 }; };

	template<uint64_t Max, uint8_t requiredBits = RequiredBits<Max>::value>
	struct SmallestUint {typedef typename SmallestUint<Max, requiredBits+1>::type type; };	
	template<uint64_t Max> struct SmallestUint<Max, 8> { typedef uint8_t  type; };
	template<uint64_t Max> struct SmallestUint<Max,16> { typedef uint16_t type; };
	template<uint64_t Max> struct SmallestUint<Max,32> { typedef uint32_t type; };
	template<uint64_t Max> struct SmallestUint<Max,64> { typedef uint64_t type; };

	struct obitstream {
		
		uint8_t *start;
		uint8_t *p;
		uint8_t *end;
		uint64_t val = 0;
		uint64_t roll = 0;
		
		constexpr obitstream(uint8_t *p, size_t sz) : start(p), p(p), end(p+sz) {}
		
		constexpr void write(size_t bits, uint64_t data) {
			
			if (p + 16 > end) { // End is near
				
				while (roll >= 8 and p!= end) { 
					*(uint8_t *)p = val + (data << roll);
					p += sizeof(uint8_t);
					roll -= 8;
				}
				
				val += (data << roll);
				roll += bits;
				
				while (roll >= 8 and p!= end) { 
					*(uint8_t *)p = val + (data << roll);
					p += sizeof(uint8_t);
					roll -= 8;
				}
				
				
			} else if (roll >= 64 - bits ) {
				
				*(uint64_t *)p = val + (data << roll);
				p += sizeof(uint64_t);
				val = (data >> (64-roll));
				roll = bits + roll - 64;
				
			} else {
				
				val += (data << roll);
				roll += bits;
			}
		}
		
		constexpr void sync() {}
		
		constexpr size_t size() { return p-start; }
		
		constexpr void writeBytes(size_t bytes, uint64_t data) {
			
			*(uint64_t *)p = data;
			p += bytes;
		}
		
		constexpr size_t bytesLeft() const {
			return end - p;
		}
	};

	struct ibitstream {
		
		const uint8_t *start;
		const uint8_t *p;
		const uint8_t *end;
		uint64_t val = 0;
		uint64_t roll = 0;
		
		constexpr ibitstream(const uint8_t *p, size_t sz) : start(p), p(p), end(p+sz) {}
		
		constexpr uint64_t read(size_t bits) {

			if (roll < bits) {
				if (p + 4 < end) {
					val += ( (*(uint32_t *)p) << roll);
					p += sizeof(uint32_t);
					roll += sizeof(uint32_t)*8;
				} else while (p < end) {
					val += ( (*(uint8_t *)p) << roll);
					p += sizeof(uint8_t);
					roll += sizeof(uint8_t)*8;
				}
			}

			uint64_t ret = val;
			roll -= bits;
			val = (val >> bits);
			return ret & ((1<<bits)-1);
		}
		
		constexpr uint64_t readBytes(size_t bytes) {
			
			uint64_t ret = *(uint64_t *)p;
			p += bytes;
			return ret;
		}
		
		constexpr operator bool() const {
			return not roll or p < end;
		}
		
		constexpr void sync() {}
		
		constexpr size_t bytesLeft() const {
			return end - p;
		}
	};

	// Uncompressed symbols may be anything from 1 bit to N bits, and it shold even be variable, and it might be padded, and and and... s**t.
	
	// Compressed symbols may even be variable.
	// Let's keep Marlin being Marlin: fixed sized uncompressed
	
	// Template parameters:
	// ALPHABET_SIZE  -> number of symbols in the alphabet
	// WORD_SIZE -> maximum number of symbols in a word
	// NUM_WORDS -> maximum number of words in the dictionary
	template<size_t ALPHABET_SIZE, size_t WORD_SIZE, size_t NUM_WORDS>
	class Dictionary {
		
		// Generic
		typedef typename SmallestUint<  ALPHABET_SIZE-1>::type Symbol;
		typedef cx::vector<Symbol, WORD_SIZE> Word;
		
		cx::vector<Word,NUM_WORDS> words;
		
		
		// Encoder
		typedef typename SmallestUint<2*NUM_WORDS-1>::type NodeIdx;
		typedef typename SmallestUint<  NUM_WORDS-1>::type WordIdx;
		
		struct Node {
			cx::array<NodeIdx,ALPHABET_SIZE> child = {};
			WordIdx code = {};
		};            
		cx::vector<Node, 2*NUM_WORDS> nodes = {};
		
		constexpr void initEncoder() {
			
			Node blank = {};
			blank.code = 0;

			for (size_t i=0; i<ALPHABET_SIZE; ++i)
				blank.child[i] = i;

			for (size_t i=0; i<ALPHABET_SIZE; ++i)
				nodes.push_back(blank);
			
			for (size_t i=0; i<words.size(); ++i) {
				
				auto &&w = words[i];
				// No word should be empty
				NodeIdx nodeId = Symbol(w[0]);
				
				for (size_t j=1; j<w.size(); ++j) {
					
					Symbol c = w[j];
					if (nodes[nodeId].child[c] == c) {
						nodes[nodeId].child[c] = nodes.size();
						nodes.push_back(blank);
					}
					nodeId = nodes[nodeId].child[c];
				}
				nodes[nodeId].code = i;
			}
		}

		
		// Decoder
		typedef typename SmallestUint<0, 
			   WORD_SIZE*RequiredBits<ALPHABET_SIZE-1>::value
			 + RequiredBits<WORD_SIZE>::value >::type Entry;
		
		constexpr static const size_t sizeShift = sizeof(Entry)*8 - RequiredBits<WORD_SIZE>::value;
		
		cx::array<Entry, NUM_WORDS> table = {};

		constexpr void initDecoder() {
			
			std::cout << (WORD_SIZE*RequiredBits<ALPHABET_SIZE-1>::value
			 + RequiredBits<WORD_SIZE) << " " << sizeof(Entry) << std::endl;
			
			for (size_t i=0; i<words.size(); ++i) {
				
				for (size_t j=0; j<words[i].size(); ++j)
					table[i] += uint64_t(words[i][j]) << (j * RequiredBits<ALPHABET_SIZE-1>::value);
				
				table[i] +=  words[i].size() << sizeShift;
			}
		}
		
	public:
	
		constexpr Dictionary(const cx::vector<Word,NUM_WORDS> words) : words(words) { 
		
			initEncoder(); 
			initDecoder(); 
		}
		 
		// It takes forever to compile as constexpr. Hardly recommended.
		constexpr Dictionary(const cx::array<double,ALPHABET_SIZE> &P) {
				
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

			cx::sort(a.begin(), a.end(), std::greater<A>());
			//cx::priority_queue<A,ALPHABET_SIZE,std::greater<A>> pq1;
			//for (size_t i=0; i<ALPHABET_SIZE; ++i) pq1.push(a[i]);
			//for (size_t i=0; i<ALPHABET_SIZE; ++i) { a[i] = pq1.top(); pq1.pop(); }
			
			
			
			//for (size_t i=0; i<ALPHABET_SIZE; ++i)
			//	std::cout << i << ": " << a[i].p << " " << uint32_t(a[i].symbol) << std::endl;

			// Initial State Probability (Ps_i)
			cx::array<double,ALPHABET_SIZE> Ps = {}; 
			Ps.front() = 1.;
			
			// 1st Preprocessing
			// Reverse cumulative symbol probability (divisor of the first term in Eq1)
			cx::array<double,ALPHABET_SIZE> P1 = {};
			P1[ALPHABET_SIZE-1] = a[ALPHABET_SIZE-1].p;
			for (size_t i=ALPHABET_SIZE; i; i--)
				P1[i-1] = P1[i] + a[i-1].p;


			for (size_t StateUpdateIterations = 10; StateUpdateIterations; --StateUpdateIterations) {

				// 2nd Preprocessing
				// Parts of eq2 that depend on the state
				cx::array<double,ALPHABET_SIZE> P2 = {};
				{
					double sum = 0.0;
					for (size_t i=0; i<ALPHABET_SIZE; ++i)
						P2[i] = sum = sum + Ps[i]/P1[i];
				}
				
				//for (size_t i=0; i<10; ++i)
				//	std::cout << i << ": " << Ps[i] << " " <<  std::endl;


				struct Node {
					
					cx::vector<Symbol, WORD_SIZE> symbols = {};
					int state = 0;
					double p = 0.0;
					constexpr bool operator>(const Node& rhs) const {
						
						if (symbols.size() == rhs.symbols.size()) return p > rhs.p;
						if (symbols.size() == WORD_SIZE) return false;
						if (rhs.symbols.size() == WORD_SIZE) return true;
						return p > rhs.p;
					}
				};

				// Build the Dictionary
				cx::priority_queue<Node,NUM_WORDS,std::greater<Node>> pq;
				{
					// DICTIONARY INITIALIZATION
					for (size_t n=0; n<ALPHABET_SIZE; ++n) {
						
						Node node = {};
						node.symbols.push_back(n);
						node.state = 0;
						node.p = P2[n] * a[n].p;
						
						pq.push(node);
					}
					
					// GROW THE DICTIONARY
					while (pq.size() < NUM_WORDS and pq.top().symbols.size() < WORD_SIZE) {
						
						Node node = pq.top();
						pq.pop();
						
						Node newNode = {};
						newNode.symbols = node.symbols;
						newNode.symbols.push_back(node.state);
	
						newNode.state = 0;
						
						newNode.p = node.p * a[node.state].p / P1[node.state];
						pq.push(newNode);
						
						node.p -= newNode.p;
						node.state++;
						
						if (node.state == ALPHABET_SIZE-1) {
							node.symbols.push_back(node.state);
							node.state = 0;
						}
						pq.push(node);
					}
					std::cout << "pq.size" << pq.size() <<  std::endl;
				}
				
				double meanLength = 0, meanLengthApprox = 0, sump=0;

				words.clear();
				for (auto&& node : pq) {
					Word word = {};
					for (auto&& s : node.symbols)
						word.push_back(a[s].symbol);
					words.push_back(word);
					meanLength += node.p*node.symbols.size();
					meanLengthApprox += node.symbols.size()/double(NUM_WORDS);
					sump += node.p;
					//std::cout << "kk " << node.p << std::endl;
				}
				std::cout << "Mean length: " << RequiredBits<NUM_WORDS-1>::value / meanLength << " " <<  RequiredBits<NUM_WORDS-1>::value / meanLengthApprox << " " << sump << std::endl;
				
				
				// Update probability state
				{
					cx::array<cx::array<double,ALPHABET_SIZE>,ALPHABET_SIZE> T = {};
					for (auto &&node : pq) 
						for (size_t state=0, w0 = node.symbols.front(); state<=w0; ++state)
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
			}
			
			initEncoder(); 
			initDecoder();
		}	


		constexpr void decode(ibitstream &src, obitstream &dst) const {
			
			constexpr size_t IN  = RequiredBits<NUM_WORDS    -1>::value;
			constexpr size_t OUT = RequiredBits<ALPHABET_SIZE-1>::value;
			while (src) {
				Entry e = table[src.read(IN)];
				dst.write((e >> sizeShift)*OUT,e);
			};
			
			
			/*
			while (src + 6 <= end) {
				
				uint64_t v64=0;
				v64 = *(const uint64_t *)src;
				src+=6;
				
				{
					uint64_t v = data[(v64>>0 ) & 0xFFF];
					*((uint64_t *)dst) = v;
					dst += v >> ((sizeof(uint64_t)-1)*8);
				}
				{				
					uint64_t v = data[(v64>>0 ) & 0xFFF];
					*((uint64_t *)dst) = v;
					dst += v >> ((sizeof(uint64_t)-1)*8);
				}
				{				
					uint64_t v = data[(v64>>0 ) & 0xFFF];
					*((uint64_t *)dst) = v;
					dst += v >> ((sizeof(uint64_t)-1)*8);
				}
				{				
					uint64_t v = data[(v64>>0 ) & 0xFFF];
					*((uint64_t *)dst) = v;
					dst += v >> ((sizeof(uint64_t)-1)*8);
				}
			}

			uint64_t v64=0, c=0;
			while (src < end) {
				while (c<12) {
					v64 += *src++ << c;
					c+= 8;
				}
				
				{
					const uint8_t *v = (const uint8_t *)&data[(v64>>0 ) & 0xFFF];
					size_t sz = v[sizeof(uint64_t)-1];
					for (size_t i=0; i<sz and dst<dstEnd; i++)
						*dst++ = *v++;
				}
				
				v64 = v64 >> 12;
				c-= 12;
			}*/
		}				

		constexpr void encode(ibitstream &src, obitstream &dst) const {

			constexpr size_t IN  = RequiredBits<ALPHABET_SIZE-1>::value;
			constexpr size_t OUT = RequiredBits<NUM_WORDS-1>::value;
			
			std::cout << "IN  " << IN << std::endl;
			std::cout << "OUT " << OUT << std::endl;
			NodeIdx nodeId = src.read(IN);
			do {
				NodeIdx oldNodeId = 0;
				do { 
					oldNodeId = nodeId;
					nodeId = nodes[nodeId].child[src.read(IN)]; 
				} while (nodeId >= ALPHABET_SIZE);
				dst.write(OUT, nodes[oldNodeId].code);
			} while (src);
		}
	};
	
	
	
}



/*
size_t Marlin_compress  (uint8_t* dst, size_t dstCapacity, const uint8_t* src, size_t srcSize, size_t bs = 4096) {
	
	std::vector<Block> blocks(1 + srcSize/bs);

	if (not distCapacity < blocks.size()*(bs + 4)) 
		throw std::runtime_error("Insufficient Output Capacity");
		
	for (size_t i=0; i<blocks.size(); i++) block.idx = i;
	
	for (auto &&block : blocks) {
		
		block.uncompressedSize = std::min(srcSize - src, Block.maxCapacity());
		block.src = src;
		src += block.uncompressedSize;
				
		// Skip compression of very small blocks
		if (sz < 256) continue;
		
		// Find the best dictionary
		block.expectedLength = block.uncompressedSize*.99;		
		std::array<size_t, 256> hist; 
		
		hist.fill(0);
		for (size_t j = 0; j<sz; j++) hist[src[j]]++;
		
		for (size_t j = 0; j<dictionaries.size(); j++) {
			double el = dictionaries[j].expectedLength(hist);
			if (el < block.expectedLength) {
				block.expectedLength = el;
				block.dictIndex = j;
			}
		}

		hist.fill(0);
		for (size_t j = 1; j<sz and j<16; j++) hist[uint8_t(int8_t(src[j])-int8_t(src[j-1]))]++;
		for (size_t j = 16; j<sz; j++) hist[uint8_t(int8_t(src[j])-int8_t(src[j-16]))]++;
		
		for (size_t j = 0; j<dictionaries.size(); j++) {
			double el = dictionaries[j].expectedLength(hist);
			if (el < block.expectedLength) {
				block.expectedLength = el;
				block.delta = true;
				block.dictIndex = j;
			}
		}
	}
	
	std::sort(
		blocks.begin(), blocks.end(), 
		[](const Block &a, const Block &b) { 
			return a.dictIndex == b.dictIndex ? a.delta < b.delta : a.dictIndex < b.dictIndex; 
		}
	);
	
	size_t outSize = 0;
	for (auto  it1 = blocks.begin(); i != blocks.end();) {
		auto it2 = it1;
		while (it2 != blocks.end() and it1->dictIndex == it2->dictIndex and it1->delta == it2->delta) 
			it2++;
		
		outSize += dictionaries[it1->dictIndex].encode(it1, it2, &dst);
	}

	return outSize;
}*/
/*
size_t Marlin_decompress(int8_t* dst, size_t dstCapacity, const int8_t* Src, size_t SrcSize) {
	
	out.resize(in.size()-1);
		assert(in.back().size()==out.size());
		uint8_t *head = in.back().begin();
		
		std::vector<std::pair<std::pair<int64_t, int64_t>, size_t>> packets;
		for (size_t i=0; i<out.size(); i++) {

			if        (head[i] == 255) {
			
				out[i] = in[i];
			} else if (head[i] == 0  ) {
				
				memset(out[i].begin(), 0, out[i].size());
				out[i][0] = in[i][0];
			} else  {
			
				packets.emplace_back(std::make_pair(head[i], -out[i].size()), i);
			}
		}
		std::sort(packets.begin(), packets.end());
		
		std::vector<std::reference_wrapper<const AlignedArray8>> rIn;
		std::vector<std::reference_wrapper<      AlignedArray8>> rOut;
		std::vector<std::reference_wrapper<const uint8_t      >> zeroCounts;

		for (size_t i=0; i<packets.size(); i++) {
			rIn       .emplace_back(std::cref(in  [packets[i].second]));
			rOut      .emplace_back(std:: ref(out [packets[i].second]));
			zeroCounts.emplace_back(std:: ref(head[packets[i].second]));
		}
		
		uncompress(rIn, rOut, zeroCounts);

		return out.nBytes();		
}
*/
