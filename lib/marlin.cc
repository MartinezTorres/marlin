
#include <util.h>
#include <distribution.h>
#include <marlin.h>


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


	// Template parameters:
	// A  -> number of symbols in the alphabet
	// WS -> maximum number of symbols in a word
	// NW -> maximum number of words in the dictionary
	template<size_t A, size_t WS, size_t NW>
	class Dictionary {
		
		typedef typename SmallestUint<A-1>::type TSymbol;

/*
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
		*/
	public:
			
		constexpr Dictionary(cx::array<double,A> Punsorted) {
			
			struct SymbolWithProbability {
				double p;
				TSymbol symbol;
				constexpr bool operator>(const SymbolWithProbability& rhs) const { return p > rhs.p; }
			};
			
			cx::array<SymbolWithProbability,A> symbols = {};
			
			for (size_t i=0; i<symbols.size(); i++)
				symbols[i] = { Punsorted[i] , TSymbol(i) };
			
			cx::sort(symbols.begin(), symbols.end(), std::greater<SymbolWithProbability>());
			
						
			cx::array<double,A> Pstate = {}; 
			Pstate.fill(0.);
			Pstate.front() = 1.;
			
			
			cx::array<double,A> PN = {};
			PN.back() = symbols.back().p;
			for (size_t i=A-1; i; i--)
				PN[i-1] = PN[i] + symbols[i-1].p;


			cx::array<double,A> Pchild = {};
			for (size_t i=0; i<A; i++)
				Pchild[i] = symbols[i].p / PN[i];
				
				
			for (size_t StateUpdateIterations = 3; StateUpdateIterations; --StateUpdateIterations) {

				struct Node {
					cx::vector<TSymbol,WS> symbolIndexes; //Not symbols but symbol idn
					int nChildren;
					double p;
					constexpr bool operator>(const Node& rhs) const {
						return ((symbolIndexes.size() == WS) xor (rhs.symbolIndexes.size() == WS)) ?
							p > rhs.p  : symbolIndexes.size() != WS;
					}
				};
				
				cx::priority_queue<Node,NW,std::greater<Node>> pq;

				// DICTIONARY INITIALIZATION
				for (size_t c=0; c<A; ++c) {
					
					double sum = 0;
					for (size_t t = 0; t<=c; ++t) sum += Pstate[t]/PN[t];
					
					Node node = {};
					node.symbolIndexes.push_back(c);
					node.nChildren = 0;
					node.p = sum * symbols[c].p;
					
					pq.push(node);
				}
				
				// GROW THE DICTIONARY
				while (pq.size() < NW and pq.top().symbolIndexes.size() < WS) {
					
					Node node = pq.top();
					pq.pop();
					
					Node newNode = {};
					newNode.symbolIndexes = node.symbolIndexes;
					newNode.symbolIndexes.push_back(node.nChildren);
					
					newNode.nChildren = 0;
					
					newNode.p = node.p * Pchild[node.nChildren];
					pq.push(newNode);
					
					node.p -= newNode.p;
					node.nChildren++;
					
					if (node.nChildren == A-1) {

						node.symbolIndexes.push_back(node.nChildren);
						node.p = node.p * Pchild[node.nChildren];
						node.nChildren = 0;
					}					
					pq.push(node);
				}
								
				// UPDATE STATE PROBABILITIES
				{
					
					cx::array<double,A> sums = {};
					{
						double sum = 0.0;
						for (size_t state = 0; state < A; ++state)
							sums[state] = sum = sum + Pstate[state]/PN[state];
					}
					
					cx::matrix<double,A,A> T = {};
			
					for (auto &&node : pq)
						for (size_t state=0; state<=node.symbolIndexes.front(); ++state)
							T[state][node.nChildren] += node.p /(sums[node.symbolIndexes.front()] * PN[state]);

					// Solving Maxwell
					for (size_t MaxwellIterations = 1; MaxwellIterations; --MaxwellIterations)
						T = T* T;

/*					Pstate = T[0];*/
				}
				
				/*{
					W.clear();
					for (auto &&node : pq) {
						W.push_back(pq.word);
						for (auto &&c: W.back().word)
							c = symbols[c].symbol;
					}					
				}*/
			}		
		}
				
		/*constexpr double expectedLength(const std::array<size_t, N> &hist) const {
			
			double ret = 0;
			for (size_t i=0; i<N; i++)
				ret += meanLengthPerSymbol[i]*hist[i];
			return ret;
		}*/
	};
	/*
	Dictionary<256> dictionaries[] = {};

	void encode(const Block *begin, const Block *end, uint8_t *&dst) {
		
		
	}
		
	
	class Block {
		uint8_t  dictIndex = 0;
		uint16_t uncompressedSize = 0;
		double expectedLength = 0;
		bool delta = false;
	};*/

}

int main() {
	
	constexpr auto dictionary  = Dictionary<256,7,4096>( Distribution::getWithEntropy(Distribution::Gaussian<256>,.5) );
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
