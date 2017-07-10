
namespace {
	
    template<uint64_t N> struct bitcount    { enum { value = 1 + bitcount<(N>>1)>::value }; };
	template<>           struct bitcount<0> { enum { value = 0 }; };
	template<>           struct bitcount<1> { enum { value = 1 }; };
    template<uint64_t N> struct bytecount   { enum { value = (bitcount<N>::value + 7) / 8 }; };

	
	template<uint64_t Max>
	struct RequiredBits { enum { value = 
			Max < 0x100         ?  8 :
			Max < 0x10000       ? 16 :
			Max < 0x100000000LL ? 32 :
								  64 };
	};

	template<int bits> struct SelectInteger_;
	template<> struct SelectInteger_ <8> { typedef uint8_t type; };
	template<> struct SelectInteger_<16> { typedef uint16_t type; };
	template<> struct SelectInteger_<32> { typedef uint32_t type; };
	template<> struct SelectInteger_<64> { typedef uint64_t type; };

	template<uint64_t Max>
	struct SelectInteger : SelectInteger_<RequiredBits<Max>::value> {};


	// Template parameters:
	// A  -> number of letters in the alphabet, maximum is 2^16
	// WS -> maximum size of the word, maximum is 128
	// NW -> number of words in the dictionary, maximum is 2^16
	template<size_t A, size_t WS, size_t NW>
	class Dictionary {
		
		struct Word : cx::vector<uint16_t, WS> {
			
			double probability;
		};
		
		struct AllWords : cx::vector<Word, NW> {
		};
		
		
		struct EncoderNode {
			cx::array<uint16_t,N> child;
			uint16_t code;
		};            
		cv:array<EncoderNode, 2*Words::capacity()> nodes = {};
			
		constexpr void prepareEncoder(const Words &words) {
			
			Node blank;
			blank.code = 0;
			blank.increment = 0;
			for (size_t i=0; i<N; i++)
				blank.child[i] = i;

			for (size_t i=0; i<N; i++)
				nodes.push_back(blank);
			
			for (size_t i=0; i<words.size(); ++i) {
				
				auto &&w = words[i];
				// No word should be empty
				uint32_t nodeId = uint8_t(w[0]);
				
				for (size_t j=1; j<w.size(); ++j) {
					
					uint8_t c = w[j];
					if (nodeId and nodes[nodeId].child[c] == c) {
						nodes[nodeId].child[c] = nodes.size();
						nodes.push_back(blank);
					}
					nodeId = nodes[nodeId].child[c];
				}
				nodes[nodeId].code = i;
			}
		}
			
		size_t encode(uint8_t* dst, const uint8_t* src, size_t srcSize) const { 
				
			uint8_t* dst0 = dst;
			
			const uint8_t *end = src+srcSize;
			if (src == end) return;

			uint32_t nodeId = *src++, oldNodeId = 0;
			uint64_t v64 = 0;
			while (src + 4*Word::capacity() < end) {
				
				do { 
					oldNodeId = nodeId;
					nodeId = nodes[nodeId].child[*src++]; 
				} while (nodeId < N);
				v64 = nodes[oldNodeId].code;
				do { 
					oldNodeId = nodeId;
					nodeId = nodes[nodeId].child[*src++]; 
				} while (nodeId < N);
				v64 += nodes[oldNodeId].code << 12;
				do { 
					oldNodeId = nodeId;
					nodeId = nodes[nodeId].child[*src++]; 
				} while (nodeId < N);
				v64 += nodes[oldNodeId].code << 24;
				do { 
					oldNodeId = nodeId;
					nodeId = nodes[nodeId].child[*src++]; 
				} while (nodeId < N);
				v64 += nodes[oldNodeId].code << 36;

				*(uint64_t *)dst = v64;
				dst += 6;
			}
			
			v64 = 0;
			uint32_t c = 0;
			while (src < end) {
				oldNodeId = nodeId;
				nodeId = nodes[nodeId].child[*src++]; 
				if (nodeId < N) {
					v64 += nodes[oldNodeId].code << c;
					c += 12;
					if (c==48) {
						*(uint64_t *)dst = v64;
						dst += 6;
						c    =  0;                            
					}
				}
			}

			do {
				oldNodeId = nodeId;
				nodeId = nodes[nodeId].child[0]; 
			} while (not (nodeId < N));
			v64 += nodes[oldNodeId].code << c;
			c += 12;
			*(uint64_t *)dst = v64;
			dst += (c+7)/8;
			
			return dst - dst0;
		}
			
		cx::array<cx::array<uint8_t,8>> decodeData;

		void decode(uint8_t* dst, size_t dstSize, const uint8_t* src, size_t srcSize) const { 
							
			uint8_t* dst0 = dst;
			
			const uint8_t *end = src+srcSize;
			if (src == end) return;
			
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
			}
		}

		constexpr prepareDecoder(const Dictionary &dict, const std::array<std::pair<double, uint8_t>,256> &distSorted) 
			: maxWordSize(dict.maxWordSize), dictSize2(dict.dictSize2) {
			
			data.resize(dict.W.size()*maxWordSize);
			for (size_t i=0; i<dict.W.size(); i++) {

				uint8_t *d = &data[i*maxWordSize];
				d[maxWordSize-1] = dict.W[i].size();
				for (auto c : dict.W[i])
					*d++ = distSorted[c].second;
			}
		}
		
		
		cx::array<N> meanLengthPerSymbol = {};
		
	public:
			
		constexpr Dictionary(cx::array<double,A> Punsorted) {
			
			struct SymbolWithProbability {
				double p;
				Tsymbol symbol;
				constexpr bool operator>(const &rhs) { return p > rhs.p; }
			};
			
			cx::array<SymbolWithProbability,A> symbols;
			
			for (size_t i=0; i<symbols.size(); i++)
				symbols[i] = { Punsorted[i] , i };
			
			cx::sort<std::greater>(symbols);
			
						
			cx::array<double,A> Pstate; 
			Pstate.fill(0.);
			Pstate.front() = 1.;
			
			
			cx::array<double,A> PN;
			PN.back() = symbols.back().p;
			for (size_t i=PN.size()-1; i; i--)
				PN[i-1] = PN[i] + symbols[i-1].p;


			cx::array<double,A> Pchild;
			for (size_t i=0; i<P.size(); i++)
				Pchild[i] = symbols[i].p / PN[i];
				
			for (size_t StateUpdateIterations = 3; StateUpdateIterations; --StateUpdateIterations) {

				struct Node {
					Word word; //Not symbols but symbol orders
					int nChildren;
					double p;
					constexpr bool operator>(const &rhs) { 
						return (not (word.size() == WS)) and (p > rhs.p); 
					}
				};
				
				cx::priority_queue<Node,NW,std::greater> pq;

				// DICTIONARY INITIALIZATION
				for (size_t c=0; c<A; ++c) {
					
					double sum = 0;
					for (size_t t = 0; t<=c; ++t) sum += Pstate[t]/PN[t];
					
					Node node;
					node.word.push_back(c);
					node.nChildren = 0;
					node.p = sum * symbols[c].p;
					
					pq.push(node);
				}
				
				// GROW THE DICTIONARY
				while (pq.size() < NW and pw.top().word.size() < WS) {
					
					Node node = pq.top();
					pq.pop();
					
					Node newNode;
					newNode.word = node.word;
					newNode.word.push_back(node.word.nChildren);
					
					newNode.nChildren = 0;
					
					newNode.p = node.p * Pchild[node.word.nChildren];
					pq.push(newNode);
					
					node.p -= newWord.p;
					node.nChildren++;
					
					if (node.nChildren == A-1) {

						node.word.push_back(node.word.nChildren);
						node.p = node.p * Pchild[node.word.nChildren];
						node.nChildren = 0;
					}					
					pq.push(node);
				}
								
				// UPDATE STATE PROBABILITIES
				{
					cx::array<cx::array<double,A>,A> T;
					for (auto &&t : T)
						t.fill(0.);

					cx::array<double,A> sums;
					{
						double sum = 0.0;
						for (size_t state = 0; state < A; ++state)
							sums[state] = sum = sum + Pstate[state]/PN[state];
					}
			
					for (auto &&node : pq)
						for (size_t state=0; state<=node.word.front(); ++state)
							T[state][node.nChildren] += node.p /(sum[node.word.front()] * PN[state]);

					// Solving Maxwell
					for (size_t MaxwellIterations = 10; MaxwellIterations; --MaxwellIterations) {
						
						auto T2 = T;
						for (auto &&t : T)
							t.fill(0.);
						for (size_t i=0; i<A; ++i)
							for (size_t j=0; j<A; ++j)
								for (size_t k=0; k<A; ++k)
									T[i][j] += T2[i][k] * T2[k][j];
					}					
					Pstate = T[0];
				}
				
				W.clear();
				for (auto &&node : pq) {
					W.push_back(pq.word);
					for (auto &&c: W.back().word)
						c = symbols[c].symbol;
				}					
			}		
		}
				
		constexpr double averageBitsPerSymbol() const {
			
			double meanLength = 0;
			for (auto &w : W) 
				meanLength += prob(w)*w.size();

			return dictSize2/meanLength;		
		}
	
		constexpr double expectedLength(const std::array<size_t, N> &hist) const {
			
			double ret = 0;
			for (size_t i=0; i<N; i++)
				ret += meanLengthPerSymbol[i]*hist[i];
			return ret;
		}
		
	};
	
	Dictionary<256> dictionaries[] = {};

	void encode(const Block *begin, const Block *end, uint8_t *&dst) {
		
		
	}
		
	
	class Block {
		uint8_t  dictIndex = 0;
		uint16_t uncompressedSize = 0;
		double expectedLength = 0;
		bool delta = false;
	};

};


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
}
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
