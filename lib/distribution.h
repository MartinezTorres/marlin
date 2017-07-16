#pragma once
#include <array>
#include <cmath>
#include <limits>
#include <cstddef>


namespace {
namespace Distribution {
	
	template<typename T, size_t N> 
	std::array<T,N> norm1(const std::array<T,N> &a) { 
			
		std::array<T,N> ret = {};
		double sum = 0;
		for (size_t n=0; n<N; ++n) sum += a[n];
		for (size_t n=0; n<N; ++n) ret[n] = a[n]/sum;
		return ret;
	}
	
    template<size_t N>
    std::array<double,N> Gaussian(double b) {
        
        std::array<double,N> arr = {}; 
        for (auto&& v : arr) v = std::numeric_limits<double>::min();

        arr[0] += 1.;
        for (int64_t i=1; i<10*int(N) and std::exp(-i*i/b )>std::numeric_limits<double>::min(); i++) {
            arr[      i  % N] += std::exp(-i*i/b);
            arr[(10*N-i) % N] += std::exp(-i*i/b);
        }
            
        return norm1(arr);
    }

    template<size_t N>
    std::array<double,N> Laplace(double b) {
        
        std::array<double,N> arr = {}; 
        for (auto&& v : arr) v = std::numeric_limits<double>::min();

        arr[0] += 1.;
        for (int64_t i=1; i<10*int(N) and std::exp(-i/b )>std::numeric_limits<double>::min(); i++) {
            arr[      i  % N] += std::exp(-i/b );
            arr[(10*N-i) % N] += std::exp(-i/b );
        }    
        return norm1(arr);
    }

    template<size_t N>
    std::array<double,N> Exponential(double b) {
        
        std::array<double,N> arr = {}; 
        for (auto&& v : arr) v = std::numeric_limits<double>::min();

        arr[0] += 1.;
        for (int64_t i=1; i<10*int(N) and std::exp(-i/b )>std::numeric_limits<double>::min(); i++)
            arr[      i  % N] += std::exp(-i/b );

        return norm1(arr);
    }

    template<size_t N>
    std::array<double,N> Poisson(double l) {
        
        std::array<double,N> arr = {}; 
        for (auto&& v : arr) v = std::numeric_limits<double>::min();

        double lp=0,kf=0;
        for (int64_t k=0; k<10*N and kf>k*std::numeric_limits<double>::min(); k++) {
            kf = (k?kf*k:1);
            lp = (k?lp*l:1);
            arr[      k  % N] += lp*std::exp(-l)/kf;
        }
        return norm1(arr);
    }
	
    template<typename F>
    double entropy(const F &pmf) {

		double distEntropy=0;
        for (auto &p : pmf) {
            if (p) {
                distEntropy += -p*std::log2(p);
            }
        }
		return distEntropy;        
    }

    template<typename F>
	auto getWithEntropy(const F &pmf, double h) {

        h = std::min(std::max(h,1e-5),1.-1e-5);
		double bMax = 1<<1, bMin = 1./(1<<30);
        double hMax = entropy(pmf(bMax)), hMin = entropy(pmf(bMin));
        double hGoal = h * std::log2(pmf(bMax).size());
        
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
}
}

