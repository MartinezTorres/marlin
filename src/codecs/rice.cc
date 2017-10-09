#include <codecs/rice.hpp>
#include <util/distribution.hpp>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <iostream>

#define UNROLL1(c, sz_, a)  { typeof(sz_) sz = sz_; while(sz>c      + 0U) { sz-= 1; static const int fast = 1; (void)fast; a; }                               while(sz--) { static const int fast = 0; (void)fast; a; } }
#define UNROLL4(c, sz_, a)  { typeof(sz_) sz = sz_; while(sz>c * 4U + 3U) { sz-= 4; static const int fast = 1; (void)fast; a;a;a;a; }                         while(sz--) { static const int fast = 0; (void)fast; a; } }
#define UNROLL8(c, sz_, a)  { typeof(sz_) sz = sz_; while(sz>c * 8U + 7U) { sz-= 8; static const int fast = 1; (void)fast; a;a;a;a;a;a;a;a; }                 while(sz--) { static const int fast = 0; (void)fast; a; } }
#define UNROLL16(c, sz_,a)  { typeof(sz_) sz = sz_; while(sz>c *16U +15U) { sz-=16; static const int fast = 1; (void)fast; a;a;a;a;a;a;a;a;a;a;a;a;a;a;a;a; } while(sz--) { static const int fast = 0; (void)fast; a; } }

//using namespace std;


struct RicePimpl : public CODEC8Z {

	static const size_t SPLITS = 64;

	Distribution::Type distType; // LAPLACE or EXPONENTIAL
	
	uint8_t  T[SPLITS][256]; // Transforms from Input space to Rice space (highest probability symbols first)
	uint8_t iT[SPLITS][256]; // Transforms from Rice space (highest probability symbols first) to Input space
	
	uint32_t bestm[SPLITS];

	int32_t cQ[SPLITS][256];
	uint64_t cR[SPLITS][256];

	static const int RICE_UTABLE_SIZE = 2*1024;
	struct RiceUTable {
		union {
			uint64_t d64;
			uint8_t d8[8];
		}; 
		uint32_t nData;
		uint32_t bDisp;
	};
	RiceUTable uT[SPLITS][RICE_UTABLE_SIZE];

	static double efficiency(const std::array<double,256> &pdf, const uint8_t *T, int m) {

		double riceCost = 1+m;
		for (int i=0;i<256;i++)
			riceCost += pdf[i] * (T[i]/(1<<m));
		
		return Distribution::entropy(pdf)/riceCost;
	}

	std::string name() const { return "Rice"; }

	RicePimpl(Distribution::Type distType_) : distType(distType_) {

		memset(&uT, 0, sizeof(uT));
		for (size_t split=0; split<SPLITS; split++) {
			
			double p = (split+0.5)/SPLITS;
			
			// Get the PDF
			std::array<double,256> pdf = Distribution::pdf(distType, p);
		
			// Calculate PDF, sorted from most probable symbol to least probable symbol
			{
				std::array<std::pair<double,uint8_t>,256> pdfSorted;
				for (uint i=0;i<256;i++) pdfSorted[i] = {pdf[i],i};
				std::sort(pdfSorted.begin(), pdfSorted.end());
				std::reverse(pdfSorted.begin(), pdfSorted.end());
				
				for (uint i=0;i<256;i++)  T[split][pdfSorted[i].second] = i; 
				for (uint i=0;i<256;i++) iT[split][i] = pdfSorted[i].second; 
			}
		
			// Calculate best M parameter for compression
			double bestEfficiency = 0.;
			for (int m=0; m<8; m++) {
				double eff = efficiency( pdf, T[split], m);
				if (eff>bestEfficiency) {
					bestEfficiency=eff;
					bestm[split]=m;
				}
			}
			
			// Precalculate compression tables
			for (uint in=0; in<256; in++) {
				int m = bestm[split];
				int M = 1<<m;
				int rice = T [split][in];
				cQ[split][in] = -(1+m+(rice>>m));
				cR[split][in] = M|(rice&(M-1));
			}

			// Precalculate uncompression tables
			for (int x=0; x<RICE_UTABLE_SIZE; x++) {
		
				uint m = bestm[split];
				uint M = 1<<m;
					
				RiceUTable &ut = uT[split][x];
				uint bit = RICE_UTABLE_SIZE/2;
				for (int d=0; d<8; d++) {
					int q = 0;
					while (bit and not (x&bit)) { q++; bit>>=1; }
					if (bit<M) break;
					bit >>= 1+m;
					if ((q<<m) + (x/(bit*2+!bit))%M > 255) break;
					ut.d8[d] = iT[split][(q<<m) + (x/(bit*2+!bit))%M];
					ut.nData ++;
					ut.bDisp += 1+m+q;
				}
			}
		}
	}

	void   compress(
		const std::vector<std::reference_wrapper<const AlignedArray8>> &in,
		      std::vector<std::reference_wrapper<      AlignedArray8>> &out,
		      std::vector<std::reference_wrapper<      uint8_t      >> &zeroCounts) const { 
			
		for (size_t j=0; j<in.size(); j++) {

			int split = zeroCounts[j]*SPLITS/256;
						
			const uint8_t *i8 = in[j].get().begin();
			uint32_t *o32 =(uint32_t *)out[j].get().begin();

			const int32_t  *Q = cQ[split];
			const uint64_t *R = cR[split];

			uint64_t st = 0;
			int32_t sts = 64;

			UNROLL16(0, in[j].get().size(), {
				
				auto i = *i8++;
				sts += Q[i];
				
				while (sts<0) {
					*o32++ = st>>32U; 
					st <<= 32U;
					sts += 32U;
				}
				st |= R[i]<<sts;
			})
			while (sts<64) {
				
				*o32++ = st>>32U;
				st <<= 32U;
				sts += 32U;
			}
			
			out[j].get().resize((uint8_t *)o32-(uint8_t *)out[j].get().begin());
		}
	}

	void uncompress(
		const std::vector<std::reference_wrapper<const AlignedArray8>> &in,
		      std::vector<std::reference_wrapper<      AlignedArray8>> &out,
		      std::vector<std::reference_wrapper<const uint8_t      >> &zeroCounts) const {
		
		for (size_t j=0; j<out.size(); j++) {
						
			uint8_t *o8 = out[j].get().begin();
			const uint32_t *i32 = (const uint32_t *)in[j].get().begin();
			uint split = zeroCounts[j]*SPLITS/256;
						
			const RiceUTable *ut = uT[split];
			uint uts = 64-std::round(std::log2(RICE_UTABLE_SIZE));
				
			uint m = bestm[split];
			uint m64 = 64-m;
			
			uint64_t st = 0;
			uint32_t sts = 64;
			
			UNROLL4(8, out[j].get().size(), {
				if (sts>=32U) {
					
					uint64_t v = *i32++;
					sts -= 32U;
					st |= v << sts;
				}

				const RiceUTable &t = ut[st>>uts];
				if (fast and t.nData) {

					*(uint64_t *)o8 = t.d64;
					o8  += t.nData;
					st <<= t.bDisp;
					sts += t.bDisp;
					sz  -= t.nData-1;

				} else {
					
					uint q=0;
					while (st<0x100000000ULL) {

						uint64_t v = *i32++;
						st <<= 32U;
						q += 32U;
						st |= v << sts;
					}

					uint leadZ = __builtin_clzll(st);
					q += leadZ;
					st <<= leadZ+1;
					sts += leadZ+1;

					if (m) {

						if (sts>m64) {
							uint64_t v = *i32++;
							sts -= 32;
							st |= v<<sts;
						}

						*o8++ = iT[split][(q<<m) + (st>>m64)];

						st <<= m;
						sts += m;
					} else {
						*o8++ = iT[split][q];
					}
				}
			})

		}
	}
};

Rice::Rice(Distribution::Type distType) : CODEC8withPimpl( new RicePimpl(distType) ) {}
