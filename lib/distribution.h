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
		
		l = std::min(l, 200.);
		l = std::max(l,1./200.);
		        
        std::array<double,N> arr = {}; 
        for (auto&& v : arr) v = std::numeric_limits<double>::min();

        double lp=1;
        for (int64_t k=0; k<10*int(N); k++) {
            lp = (k?lp*l/k:1);
            if (lp > 1e100 or lp < 1e-100) break;
            arr[      k  % N] += lp*std::exp(-l);
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

        h = std::min(std::max(h,1e-5),1.-1e-5) * std::log2(pmf(0.0).size());

		double b=1<<16;
		// Estimate parameter b from p using dicotomic search
		double stepSize = 1<<15;
		while (stepSize>1E-12) {
			if (h > entropy(pmf(b))) b+=stepSize;
			else b-=stepSize;
			stepSize/=2.;
		}
		return pmf(b);
	}
}
}

