
namespace {
	
	class Block {
		uint8_t  dictIndex = 0;
		uint16_t uncompressedSize = 0;
		double expectedLength = 0;
		bool delta = false;
	};

	template<size_t N>
	class Dictionary {
        
        typedef cx::vector<uint8_t, 7> Word;
        typedef cx::vector<Word, 1<<12> Words;
        
        
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
		
		std::vector<double> PnextState, PcurrState;
		
/*		struct Word : public std::vector<uint8_t> {
			uint8_t nextState = 0;
			using vector::vector;
		};*/
        
        constexpr double phi(const Word &w) const { return prob(w); }
		
		constexpr std::vector<Word> split(const Word &w) const {
		
			std::vector<Word> ret(2,w);
			ret[0].push_back(w.nextState); ret[0].nextState = 0;
										   ret[1].nextState++;
										   
			if (ret[1].nextState == P.size()-1) {
				ret[1].push_back(P.size()-1);
				ret[1].nextState = 0;
			}
			return ret;
		}	
		
		// This calculates the probability, for a word to be used to encode a sequence.
		constexpr double prob(const Word &w) const {

			double ret = 0;
			for (size_t t = 0; t<=w[0]; t++) {
				
				double p = PcurrState[t]; // Start with the probability of being in such state
				
				p *= P[w[0]]/PnextState[t];
				
				ret += p;
			}
				
			for (size_t i=1; i<w.size(); i++)
				ret *= P[w[i]];

			ret *= PnextState[w.nextState];

			return ret;
		}
		
		// For each symbol, it updates the probabilities of the state chain.
		constexpr void updatePcurrState() {
			
			std::vector<std::vector<double>> T(P.size(),std::vector<double>(P.size(),0.));

			for (auto &w : W) {
				
				double pw = 1.;
				for (auto &c : w) pw *= P[c];
				pw *=  PnextState[w.nextState];
				
				for (size_t s=0; s<=w[0]; s++)
					T[s][w.nextState] += pw/PnextState[s];
					
			}
			
			int t = 10;
			double diff = 0;
			do {
				
				auto T2 = T;
				for (size_t i=0; i<P.size(); i++) {
					for (size_t j=0; j<P.size(); j++) {
						T2[i][j]=0;
						for (size_t k=0; k<P.size(); k++)
							T2[i][j] += T[i][k] * T[k][j];
					}
				}
				diff = 0;
				for (size_t i=0; i<P.size(); i++)
					for (size_t j=0; j<P.size(); j++)
						diff += std::abs(T[i][j]-T2[i][j]);

				//std::cerr << diff << std::endl;
				T = T2;

			} while (t-- and diff>.00001);
			
			PcurrState = T[0];
		}
		
	public:
	
		const std::vector<double>  P;
		const size_t dictSize2, maxWordSize;

		std::vector<Word> W;
		
		Dictionary(const std::vector<double> &P, size_t dictSize2, size_t maxWordSize) :
			P(P), dictSize2(dictSize2), maxWordSize(maxWordSize) {


			PnextState = P;  
			for (size_t i=P.size()-1; i; i--)
				PnextState[i-1] += PnextState[i];
				
			PcurrState = std::vector<double>(P.size(), 0.);
			PcurrState[0] = 1.;
		}

		Dictionary &tune() {

			std::cerr << "Tune: " << dictSize2 << " and " << maxWordSize << std::endl;

			// Compares using phi operation. It also ensures that big words will be sorted last. 
			auto phiComp = [this](const Word &lhs, const Word &rhs) { 
				
				if (lhs.size()==maxWordSize-1) return true;
				if (rhs.size()==maxWordSize-1) return false;
				return this->phi(lhs) < this->phi(rhs); 
			};

			// we update the P0 probability iteratively
			for (int tries = 3; tries; tries--) {
				
				W.clear();
			
				for (size_t i=0; i<P.size(); i++)
					W.emplace_back(1,i);

				std::make_heap(W.begin(), W.end(), phiComp);
								
				while (W.size() < (1U<<dictSize2) and W.front().size() < maxWordSize-1) {
					
					auto child = split(W.front());
					
					if (child.size() + W.size() > (1U<<dictSize2)) break;
					
					std::pop_heap(W.begin(), W.end(), phiComp);					
					W.pop_back();
					for (auto &w : child) {
						W.push_back(w);
						std::push_heap(W.begin(), W.end(), phiComp);
					}
				}
				
				updatePcurrState();
				
			}
			
			std::sort(W.begin(), W.end());
			
			return *this;
		}
		
		
		Dictionary &tune(std::vector<Word> &_W) {
			
			W = _W;
			updatePcurrState();
			return *this;
		}

				
		double averageBitsPerSymbol() const {
			
			double meanLength = 0;
			for (auto &w : W) 
				meanLength += prob(w)*w.size();

			return dictSize2/meanLength;		
		}
		
		double averageBitsPerSymbolEmpirically() const {
			
			static const size_t TEST_SIZE = 1<<16;
			
			std::vector<uint8_t> testData(TEST_SIZE);
			// create data
			{
				size_t i=0;
				double ap=0;
				for (size_t j=0; j<P.size(); j++) {
					while (i<(ap+P[j])*testData.size())
						testData[i++]=j;
					ap += P[j];
				}
			}
			std::random_shuffle(testData.begin(), testData.end());
			
			// pad with 0s
			for (size_t i=0; i<maxWordSize; i++)
				testData.push_back(0);

			double meanLengthE = 0;
			size_t nWords = 0;
			for (size_t i=0, longest=0; i<TEST_SIZE; i+=longest) {
				
				Word const* lw = nullptr;
				for (auto &w : W)
					if (w[0] == testData[i])
						if (w[w.size()-1] == testData[i+w.size()-1])
							if (w==std::vector<uint8_t>(&testData[i], &testData[i+w.size()]))
								if (lw == nullptr or w.size()>lw->size())
									lw = &w;
									
				longest = lw->size();
				meanLengthE += lw->size();
				nWords++;
			}
			
			double meanLength = 0;
			for (auto &w : W) 
				meanLength += prob(w)*w.size();

			printf("Empirical Meanlength: %3.2lf\n", (meanLengthE/nWords)/meanLength);

			return double(nWords*dictSize2)/TEST_SIZE;
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
