#include <marlin.h>
#include <distribution.h>


namespace {
	
	
	class Encoder {
		
		size_t ALPHABET_SIZE;
		size_t WORD_SIZE;
		size_t NUM_WORDS;
	
		struct EncoderNode {
			cx::array<uint16_t,N> child;
			uint16_t code;
		};            
		cv:array<EncoderNode, 2*Words::capacity()> nodes = {};
			
		
			
		size_t encode8to12(uint8_t* dst, const uint8_t* src, size_t srcSize) const { 
				
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
		
	public:
		constexpr Encoder(const Words &words) {
			
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
		
		size_t operator()(const void* src, size_t srcSize, void* dst) const {
			
		}
	};
	
/*	constexpr Encoder encoders[] = {
#include "dictionaries.inc"
	};*/
	

}

int main() {
	
	//constexpr auto dictionary  = Dictionary<256,7,4096>( Distribution::getWithEntropy(Distribution::Gaussian<256>,.5) );
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
