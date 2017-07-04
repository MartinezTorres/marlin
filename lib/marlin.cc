
namespace {

	class Dictionary {
	protected:

		std::vector<double> PnextState, PcurrState;
		
		struct Word : public std::vector<uint8_t> {
			uint8_t nextState = 0;
			using vector::vector;
		};
		
		virtual double phi(const Word &) const = 0;
		virtual std::vector<Word> split(const Word &) const = 0;
		
		// This calculates the probability, for a word to be used to encode a sequence.
		double prob(const Word &w) const {

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
		void updatePcurrState() {
			
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
	};

	class TunstallDictionary : public Dictionary {
		
		using Dictionary::Dictionary;
		
		virtual double phi(const Word &w) const { return prob(w); }
		
		virtual std::vector<Word> split(const Word &w) const {
		
			std::vector<Word> ret(P.size(), w);
			for (size_t i=0; i<P.size(); i++)
				ret[i].push_back(i);

			return ret;
		}	
	};

	class MarlinDictionary : public Dictionary {

		using Dictionary::Dictionary;
		
		virtual double phi(const Word &w) const { return prob(w); }
		
		virtual std::vector<Word> split(const Word &w) const {
		
			std::vector<Word> ret(2,w);
			ret[0].push_back(w.nextState); ret[0].nextState = 0;
										   ret[1].nextState++;
										   
			if (ret[1].nextState == P.size()-1) {
				ret[1].push_back(P.size()-1);
				ret[1].nextState = 0;
			}
			return ret;
		}	
	};

	class Encoder {

		struct Node {
			std::array<uint32_t,256> child;
			uint32_t code;
			uint32_t increment;
		};
		std::vector<Node> nodes;
		
		void add(const std::vector<uint8_t> &s, uint32_t code) {
			
			Node blank;
			blank.code = 0;
			blank.increment = 0;
			for (size_t i=0; i<256; i++)
				blank.child[i] = i+1;

			if (nodes.empty()) {
				nodes.push_back(blank);
				for (size_t i=0; i<256; i++) {
					nodes.push_back(blank);
					nodes.back().increment = 1;
				}
			}
				
			uint32_t nodeId = 0;
			for (uint8_t c : s) {
				if (nodeId and nodes[nodeId].child[c] == uint32_t(c)+1) {
					nodes[nodeId].child[c] = nodes.size();
					nodes.push_back(blank);
				}
				nodeId = nodes[nodeId].child[c];
			}
			nodes[nodeId].code = code;
		}
		
		void bitCrush(AlignedArray8 &out) const {
		
			if (dictSize2==16) return;
			if (dictSize2 >16) std::runtime_error ("dictSize not supported");
			if (dictSize2 <8 ) std::runtime_error ("dictSize not supported");
			
			const uint16_t *i = (const uint16_t *)out.begin();
			uint8_t *o = out.begin();
			uint64_t v = 0; uint32_t c = 0;
			while (i!=(const uint16_t *)out.end()) {
				v += *i++ << c;
				c += dictSize2;
				while (c>=8) {
					*o++ = v & 0xFF;
					v >>= 8;
					c -= 8;
				}
			}
			if (c) *o++ = v & 0xFF;

			out.resize(o-out.begin());
		}

	public: 

		uint maxWordSize;
		uint dictSize2;
	
		Encoder() {}
	
		Encoder(const Dictionary &dict, const std::array<std::pair<double, uint8_t>,256> &distSorted) 
			: maxWordSize(dict.maxWordSize), dictSize2(dict.dictSize2) {
			
			for (size_t i=0; i<dict.W.size(); i++) {
				std::vector<uint8_t> w = dict.W[i];
				for (auto &c : w) c = distSorted[c].second;
				add (w, i);
			}
		}
		
		void operator()(const AlignedArray8 &in, AlignedArray8 &out) const { // Slower and I don't know why.
			
			const uint8_t *i = in.begin();
			const uint8_t *end = in.end();
			uint16_t *o = (uint16_t *)out.begin();

			uint32_t nodeId = 1 + *i++;
			while (i<end) {
				*o = nodes[nodeId].code;
				nodeId = nodes[nodeId].child[*i++];
				o += nodes[nodeId].increment;
			};
			do {
				*o = nodes[nodeId].code;
				nodeId = nodes[nodeId].child[0];
			} while (not nodes[nodeId].increment);
			o++;
			
			out.resize( ((uint8_t *)o)-out.begin() );
			
			bitCrush(out);
		}
	};

	class Decoder {

		AlignedArray<uint8_t,8<<20> data;
		uint maxWordSize;
		uint dictSize2;

		template<typename T, size_t N>
		void uncompressA(const AlignedArray8 &in, AlignedArray8 &out) const  {
			
			uint8_t *o = out.begin();
			const uint32_t *i = (const uint32_t *)in.begin();
			
			std::array<T,N>  *DD = (std::array<T,N> *)&data[0];
			uint64_t v32=0; uint32_t c=0;
			while (o<out.end()) {
				
				if (c<32) {
					v32 += uint64_t(*i++) << c;
					c   += 32;
				}
				{				
					uint8_t *&&v = (uint8_t *)&DD[v32 & ((1<<dictSize2)-1)];
					v32 >>= dictSize2;
					c -= dictSize2;
					for (size_t n=0; n<N; n++)
						*(((T *)o)+n) = *(((T *)v)+n);
					o += v[N*sizeof(T)-1];
				}
				{				
					uint8_t *&&v = (uint8_t *)&DD[v32 & ((1<<dictSize2)-1)];
					v32 >>= dictSize2;
					c -= dictSize2;
					for (size_t n=0; n<N; n++)
						*(((T *)o)+n) = *(((T *)v)+n);
					o += v[N*sizeof(T)-1];
				}
			}		
		}
		
		template<typename T>
		void uncompress12(const AlignedArray8 &in, AlignedArray8 &out) const  {
			
			uint8_t *o = out.begin();
			const uint8_t *i = (const uint8_t *)in.begin();
			
			T *D = (T *)data.begin();
			while (o<out.end()) {
				
				uint64_t v64=0;
				v64 = *(const uint64_t *)i;
				i+=6;
				
				{				
					T v = D[(v64>>0 ) & 0xFFF];
					*((T *)o) = v;
					o += v >> ((sizeof(T)-1)*8);
				}
				{				
					T v = D[(v64>>12) & 0xFFF];
					*((T *)o) = v;
					o += v >> ((sizeof(T)-1)*8);
				}
				{				
					T v = D[(v64>>24) & 0xFFF];
					*((T *)o) = v;
					o += v >> ((sizeof(T)-1)*8);
				}
				{				
					T v = D[(v64>>36) & 0xFFF];
					*((T *)o) = v;
					o += v >> ((sizeof(T)-1)*8);
				}
			}
		}

	public:

		Decoder() {}
		
		Decoder(const Dictionary &dict, const std::array<std::pair<double, uint8_t>,256> &distSorted) 
			: maxWordSize(dict.maxWordSize), dictSize2(dict.dictSize2) {
			
			data.resize(dict.W.size()*maxWordSize);
			for (size_t i=0; i<dict.W.size(); i++) {

				uint8_t *d = &data[i*maxWordSize];
				d[maxWordSize-1] = dict.W[i].size();
				for (auto c : dict.W[i])
					*d++ = distSorted[c].second;
			}
		}

		void operator()(const AlignedArray8 &in, AlignedArray8 &out) const {

			if (dictSize2==12) {
				switch (maxWordSize) {
					case   4: return uncompress12<uint32_t>(in, out);
					case   8: return uncompress12<uint64_t>(in, out);
				}
			} 
			switch (maxWordSize) {
				case   4: return uncompressA<uint32_t, 1>(in, out);
				case   8: return uncompressA<uint64_t, 1>(in, out);
				case  16: return uncompressA<uint64_t, 2>(in, out);
				case  32: return uncompressA<uint64_t, 4>(in, out);
				case  64: return uncompressA<uint64_t, 8>(in, out);
				case 128: return uncompressA<uint64_t,16>(in, out);
				default: throw std::runtime_error ("unsupported maxWordSize");
			}
		}

	};

	std::vector<std::shared_ptr<Encoder>> encoders;
	std::vector<std::shared_ptr<Decoder>> decoders;
	
	std::string coderName;
	std::string name() const { return coderName; }
	
	MarlinPimpl(Distribution::Type distType, Marlin::Type dictType, size_t dictSize2, size_t numDict) {

		{
			std::ostringstream oss;
			oss << "Marlin " << (distType==Distribution::Laplace?"Lap:":"Exp:") <<  (dictType == Marlin::TUNSTALL?"T":"M") << ":" << dictSize2 << ":" << numDict;
			coderName = oss.str();
		}

		std::vector<std::shared_ptr<Dictionary>> builtDictionaries(numDict);
		std::vector<std::shared_ptr<Encoder   >> builtEncoders    (numDict);
		std::vector<std::shared_ptr<Decoder   >> builtDecoders    (numDict);

//		#pragma omp parallel for
		for (size_t p=0; p<numDict; p++) {
			
			std::array<double,256> dist; dist.fill(0.);
			for (double i=0.05; i<1; i+=0.1) {
				
				std::array<double,256> pdf = Distribution::pdf(distType, (p+0.5)/numDict);
				for (size_t j=0; j<256; j++)
					dist[j] += pdf[j]/10.;
			}
			
			std::array<std::pair<double, uint8_t>,256> distSorted; 
			for (size_t i=0; i<256; i++) distSorted[i] = std::make_pair(dist[i],i);
			std::sort   (distSorted.begin(), distSorted.end());
			std::reverse(distSorted.begin(), distSorted.end());
			
			std::vector<double> P;
			for (auto &ds : distSorted)
				P.emplace_back(ds.first);

			double bestBitsPerSymbol = 1e10;
			for (int maxWordLength=4; maxWordLength <= 128; maxWordLength*=2) {

				std::shared_ptr<Dictionary> dict;
				if (dictType == Marlin::TUNSTALL) dict.reset(new TunstallDictionary(P, dictSize2, maxWordLength));
				if (dictType == Marlin::MARLIN) dict.reset(new MarlinDictionary(P, dictSize2, maxWordLength));

				dict->tune();
				double averageBitsPerSymbol = dict->averageBitsPerSymbol();
				
//				std::cerr << "P: " << p << " " << maxWordLength << " " << averageBitsPerSymbol << std::endl;

				if (averageBitsPerSymbol*1.005 > bestBitsPerSymbol) break;
				
				builtDictionaries[p] = dict;
				bestBitsPerSymbol = averageBitsPerSymbol;
			}
			
			builtEncoders[p] = std::make_shared<Encoder>(*builtDictionaries[p], distSorted);
			builtDecoders[p] = std::make_shared<Decoder>(*builtDictionaries[p], distSorted);
			
//			#pragma omp critical
//			std::cerr << "P: " << p << " " << " " << builtDictionaries[p]->averageBitsPerSymbol() << " " << builtDictionaries[p]->averageBitsPerSymbolEmpirically() << std::endl;
		}
		
		encoders.resize(256);
		decoders.resize(256);
		
		size_t bestDictionary = 0;
		for (size_t h=0; h<256; h+=4) {
			
			std::array<double,256> pdf = Distribution::pdf(distType, (h+2)/256.);
			std::sort(pdf.rbegin(), pdf.rend());
			std::vector<double> P;
			for (auto &p : pdf) P.push_back(p);
			
			double bps0 = MarlinDictionary(P, 
				builtDictionaries[bestDictionary  ]->dictSize2, 
				builtDictionaries[bestDictionary  ]->maxWordSize)
				.tune(builtDictionaries[bestDictionary  ]->W).averageBitsPerSymbol();
				
			double bps1 = 1e10;
			if (bestDictionary < builtDictionaries.size()-1)
				bps1= MarlinDictionary(P, 
				builtDictionaries[bestDictionary+1]->dictSize2, 
				builtDictionaries[bestDictionary+1]->maxWordSize)
				.tune(builtDictionaries[bestDictionary+1]->W).averageBitsPerSymbol();
				
			if (bps1<bps0) bestDictionary = bestDictionary+1;
			
			if (std::min(bps0, bps1)>8*.99) break;
						
//			encoders[h] = builtEncoders[h*builtEncoders.size()/256];
//			decoders[h] = builtDecoders[h*builtDecoders.size()/256];
			
			for (size_t hh = 0; hh<4; hh++) {
				encoders[h+hh] = builtEncoders[bestDictionary];
				decoders[h+hh] = builtDecoders[bestDictionary];
			}
		}
	}


	
	
	class Block {
		uint8_t  dictIndex = 0;
		uint16_t uncompressedSize = 0;
		double expectedLength = 0;
		bool delta = false;
	};

	template<size_t N, size_t W, size_t WL>
	class Dictionary {
        
        class Encoder {
            
            struct Node {
                cx::array<uint16_t,256> child;
                uint16_t code;
                uint16_t increment;
            };
            
            cv:array<Node, W*2> nodes = {};
            
            
					
			void add(const cx::array<uint8_t, N> &translation, uint16_t code) {
				
				Node blank;
				blank.code = 0;
				blank.increment = 0;
				for (size_t i=0; i<N; i++)
					blank.child[i] = i+1;

				if (nodes.empty()) {
					nodes.push_back(blank);
					for (size_t i=0; i<256; i++) {
						nodes.push_back(blank);
						nodes.back().increment = 1;
					}
				}
					
				uint32_t nodeId = 0;
				for (uint8_t c : s) {
					if (nodeId and nodes[nodeId].child[c] == uint32_t(c)+1) {
						nodes[nodeId].child[c] = nodes.size();
						nodes.push_back(blank);
					}
					nodeId = nodes[nodeId].child[c];
				}
				nodes[nodeId].code = code;
			}
			
		public: 

			uint maxWordSize;
			uint dictSize2;
		
			Encoder() {}
		
			Encoder(const cx::array<cx::array<N,WL>,W> &words) {
				
                Node &root = nodes[0];
				root.code = 0;
				root.increment = 0;
				for (size_t i=0; i<N; i++)
					root.child[i] = i+1;

                for (size_t i=0; i<N; i++) {
                    nodes[i+1] = root;
                    nodes[i+1].increment = 1;
				}
                
				for (auto &&w : words) {
                    
                    uint32_t nodeId = 0;
                    for (uint8_t c : s) {
                        if (nodeId and nodes[nodeId].child[c] == uint32_t(c)+1) {
                            nodes[nodeId].child[c] = nodes.size();
                            nodes.push_back(blank);
                        }
                        nodeId = nodes[nodeId].child[c];
                    }
                    nodes[nodeId].code = code;
				}
			}
			
			void operator()(const AlignedArray8 &in, AlignedArray8 &out) const { // Slower and I don't know why.
				
				const uint8_t *i = in.begin();
				const uint8_t *end = in.end();
				uint16_t *o = (uint16_t *)out.begin();

				uint32_t nodeId = 1 + *i++;
				while (i<end) {
					*o = nodes[nodeId].code;
					nodeId = nodes[nodeId].child[*i++];
					o += nodes[nodeId].increment;
				};
				do {
					*o = nodes[nodeId].code;
					nodeId = nodes[nodeId].child[0];
				} while (not nodes[nodeId].increment);
				o++;
				
				out.resize( ((uint8_t *)o)-out.begin() );
				
				bitCrush(out);
			}
		};
		
		cx::array<N> meanLengthPerSymbol = {};
		
	public:
	
		constexpr double expectedLength(const std::array<size_t, N> &hist) const {
			
			double ret = 0;
			for (size_t i=0; i<N; i++)
				ret += meanLengthPerSymbol[i]*hist[i];
			return ret;
		}
		
		void encode(const Block *begin, const Block *end, uint8_t *&dst) {
			
			
		}
		
	};
	
	Dictionary<256> dictionaries[] = {};
	
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
