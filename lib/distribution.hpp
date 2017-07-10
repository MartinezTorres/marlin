#pragma once
#include <util.hpp>

namespace { // internal linkage
namespace Distribution {
	
    template<size_t N>
    constexpr cx::array<double,N> Gaussian(double b) {
        
        cx::array<double,N> arr = {}; arr.fill(std::numeric_limits<double>::min());
        arr[0] += 1.;
        for (int64_t i=1; i<10*int(N) and cx::exp(-i*i/b )>std::numeric_limits<double>::min(); i++) {
            arr[      i  % N] += cx::exp(-i*i/b);
            arr[(10*N-i) % N] += cx::exp(-i*i/b);
        }
            
        return cx::norm1(arr);
    }

    template<size_t N>
    constexpr cx::array<double,N> Laplace(double b) {
        
        cx::array<double,N> arr = {}; arr.fill(std::numeric_limits<double>::min());
        arr[0] += 1.;
        for (int64_t i=1; i<10*int(N) and cx::exp(-i/b )>std::numeric_limits<double>::min(); i++) {
            arr[      i  % N] += cx::exp(-i/b );
            arr[(10*N-i) % N] += cx::exp(-i/b );
        }    
        return cx::norm1(arr);
    }

    template<size_t N>
    constexpr cx::array<double,N> Exponential(double b) {
        
        cx::array<double,N> arr = {}; arr.fill(std::numeric_limits<double>::min());
        arr[0] += 1.;
        for (int64_t i=1; i<10*int(N) and cx::exp(-i/b )>std::numeric_limits<double>::min(); i++)
            arr[      i  % N] += cx::exp(-i/b );

        return cx::norm1(arr);
    }

    template<size_t N>
    constexpr cx::array<double,N> Poisson(double l) {
        
        cx::array<double,N> arr = {}; arr.fill(std::numeric_limits<double>::min());
        double lp=0,kf=0;
        for (int64_t k=0; k<10*N and kf>k*std::numeric_limits<double>::min(); k++) {
            kf = (k?kf*k:1);
            lp = (k?lp*l:1);
            arr[      k  % N] += lp*cx::exp(-l)/kf;
        }
        return cx::norm1(arr);
    }
	
    template<typename F>
    constexpr double entropy(const F &pmf) {

		double distEntropy=0;
        for (auto &p : pmf) {
            if (p) {
                distEntropy += -p*cx::log2(p);
            }
        }
		return distEntropy;        
    }

    template<typename F>
	constexpr auto getWithEntropy(const F &pmf, double h) {

        h = std::min(std::max(h,1e-5),1.-1e-5);
		double bMax = 1<<1, bMin = 1./(1<<30);
        double hMax = entropy(pmf(bMax)), hMin = entropy(pmf(bMin));
        double hGoal = h * cx::log2(pmf(bMax).size());
        
        while (hMax<hGoal) {
            bMax *= 2;
            hMax = entropy(pmf(bMax));
        }
        
        while (hMax-hMin>(1e-6)) {

            double bHinge = bMin+(bMax-bMin)*( .95*((hGoal-hMin)/(hMax-hMin)) + .025);
            double hHinge = entropy(pmf(bHinge));

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
/*
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
	}*/
}
}

