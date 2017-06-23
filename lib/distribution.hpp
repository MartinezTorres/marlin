//#pragmonce
#include <vector>
#include <array>
#include <string>
#include <limits>
#include <cmath>
#include <ieee754.h>
#include <numeric>

#include <iostream>

#include <util.hpp>

namespace { // internal linkage
namespace Distribution {
	
    
    template<size_t N>
    constexpr cx::array<N> Gaussian(double b) {
        
        cx::array<N> arr;
        for (int64_t i=-10*N+1; i<10*int(N); i++)
            arr[(10*N+i) % N] += cx::exp(-i*i/b);
            
        return arr.norm1();
    }

    template<size_t N>
    constexpr cx::array<N> Laplace(double b) {
        
        cx::array<N> arr;
        for (int64_t i=1; i<10*N; i++) {
            arr[      i  % N] += cx::exp(-i/b );
            arr[(10*N-i) % N] += cx::exp(-i/b );
        }    
        return arr.norm1();
    }

    template<size_t N>
    constexpr cx::array<N> Exponential(double b) {
        
        cx::array<N> arr;
        arr[0] += 1.;
        for (int64_t i=1; i<10*N; i++)
            arr[      i  % N] += cx::exp(-i/b );

        return arr.norm1();
    }

    template<size_t N>
    constexpr cx::array<N> Poisson(double l) {
        
        cx::array<N> arr;
        double lp=0,kf=0;
        for (int64_t k=0; k<10*N and kf>k*std::numeric_limits<double>::min(); k++) {
            kf = (k?kf*k:1);
            lp = (k?lp*l:1);
            arr[      k  % N] += lp*cx::exp(-l)/kf;
        }
        return arr.norm1();
    }
	
    template<typename F>
    constexpr long double entropy(const F &pmf) {

		long double distEntropy=0;
        for (auto &p : pmf)
            if (p)
                distEntropy += -p*cx::log2(p);
		return distEntropy;        
    }
    
    
    template<typename F>
	constexpr auto getWithEntropySlow(const F &pmf, double h) {

		double b=1<<16;
		// Estimate parameter b from p using dicotomic search
		double stepSize = 1<<15;
		while (stepSize>1e-10) {
			if (h > entropy(pmf(b))/cx::log2( pmf(b).size()) ) b+=stepSize;
			else b-=stepSize;
			stepSize/=2.;
		}

		return pmf(b);
	}

    template<typename F>
	constexpr auto getWithEntropy(const F &pmf, double h) {

		double bMax = 1<<16, bMin = 1./(1<<30);
        double hMax = entropy(pmf(bMax)), hMin = entropy(pmf(bMin));
        double hGoal = h * cx::log2(pmf(bMax).size());
        
        while (hMax-hMin>(1e-7)) {
            cout << "A " << bMax << " " << bMin << " " << hMax << " " << hMin << " " << ((hGoal-hMin)/(hMax-hMin)) << endl;          

            double bHinge = bMin+(bMax-bMin)*((hGoal-hMin)/(hMax-hMin));
            double hHinge = entropy(pmf(bHinge));

            cout << "B " << bMax << " " << bMin << " " << bHinge << " " << ((hHinge-hMin)/(hHinge-hMin)) << endl;          

            if ( (hGoal-hHinge)/(hMax-hHinge) >=0. ) {
                hMin = hHinge;
                bMin = bHinge;
            } else {
                hMax = hHinge;
                bMax = bHinge;
            } 
        }

		return pmf((bMax+bMin)/2);
	}

	template<typename T, typename F>
	std::vector<uint8_t> getResiduals(const F &pmf, size_t S) {

		int8_t cdf[0x100000];
		uint j=0;
		double lim=0;
		for (uint i=0; i<pmf.size(); i++) {
			lim += pmf[i]*0x100000;
			uint ilim = round(lim);
			while (j<ilim)
				cdf[j++]=i;
		}

		std::vector<uint8_t> ret;
		uint32_t rnd =  135154;
		for (size_t i=0; i<S; i++) {
			rnd = 36969 * (rnd & 65535) + (rnd >> 16);
			ret.push_back(cdf[rnd&0xFFFF]);
		}

		return ret;
	}
}
}

