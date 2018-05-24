#pragma once
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <numeric>

namespace Distribution {

	inline std::vector<double> norm1(const std::vector<double> &pmf) {

		double sum  = std::accumulate(pmf.rbegin(), pmf.rend(), 0.);
		for (auto &&v : pmf) v/=sum;
		return pmf;
	}

	inline std::vector<double> PMFNormal(size_t N, double b) {

		std::vector<double> pmf(N, std::numeric_limits<double>::min());
		for (int i=-10*int(N)+1; i<int(10*N); i++)
			pmf[(10*N+i) % N] += std::exp(-double(  i*i)/b );
		return norm1(pmf);
	}

	inline std::vector<double> PMFLaplace(size_t N, double b) {

		std::vector<double> pmf(N, std::numeric_limits<double>::min());
		pmf[0] += 1.;
		for (size_t i=1; i<10*N; i++) {
			pmf[      i  % N] += std::exp(-double(  i)/b );
			pmf[(10*N-i) % N] += std::exp(-double(  i)/b );
		}
		return norm1(pmf);
	}

	inline std::vector<double> PMFExponential(size_t N, double b) {

		std::vector<double> pmf(N, std::numeric_limits<double>::min());
		pmf[0] += 1.;
		for (size_t i=1; i<10*N; i++)
			pmf[      i  % N] += std::exp(-double(  i)/b );

		return norm1(pmf);
	}

	inline double entropy(const std::vector<double> &pmf) {

		double distEntropy=0;
		for (size_t i=0;i<pmf.size();i++)
			if (pmf[i]>0.)
				distEntropy += -pmf[i]*std::log2(pmf[i]);

		return distEntropy;
	}
/*
	inline std::vector<double> pmf(size_t N, Type type, double h) {

		auto *dist = &PDFGaussian;
		if      (type == Gaussian   ) dist = &PDFGaussian;
		else if (type == Laplace    ) dist = &PDFLaplace;
		else if (type == Exponential) dist = &PDFExponential;
		else throw std::runtime_error("Unsupported distribution");

		double b=1<<16;
		// Estimate parameter b from p using dicotomic search
		double stepSize = 1<<15;
		while (stepSize>1E-12) {
			if (h > entropy(dist(N,b))/std::log2(N) ) b+=stepSize;
			else b-=stepSize;
			stepSize/=2.;
		}

		//std::cerr << "b: " << b << std::endl;

		return dist(N,b);
	}*/
}

int main() {
	
	
	
}
