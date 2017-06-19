#pragma once
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <numeric>

class Distribution {


	static inline std::vector<double> norm1(std::vector<double> pdf) {

		double sum  = std::accumulate(pdf.rbegin(), pdf.rend(), 0.);
		for (auto &v : pdf) v/=sum;
		return pdf;
	}

	static inline std::vector<double> PDFGaussian(size_t N, double b) {

		std::vector<double> pdf(N, 1e-100);
		for (int i=-10*int(N)+1; i<int(10*N); i++)
			pdf[(10*N+i) % N] += std::exp(-double(  i*i)/b );
		return norm1(pdf);
	}

	static inline std::vector<double> PDFLaplace(size_t N, double b) {

		std::vector<double> pdf(N, 1e-100);
		pdf[0] += 1.;
		for (size_t i=1; i<10*N; i++) {
			pdf[      i  % N] += std::exp(-double(  i)/b );
			pdf[(10*N-i) % N] += std::exp(-double(  i)/b );
		}
		return norm1(pdf);
	}

	static inline std::vector<double> PDFExponential(size_t N, double b) {

		std::vector<double> pdf(N, 1e-100);
		pdf[0] += 1.;
		for (size_t i=1; i<10*N; i++)
			pdf[      i  % N] += std::exp(-double(  i)/b );

		return norm1(pdf);
	}

public:

	enum Type { Gaussian, Laplace, Exponential };


	static inline double entropy(const std::vector<double> &pdf) {

		double distEntropy=0;
		for (size_t i=0;i<pdf.size();i++)
			if (pdf[i]>0.)
				distEntropy += -pdf[i]*std::log2(pdf[i]); //Should'n I use log2?

		return distEntropy;
	}

	template<size_t N>
	static inline double entropy(const std::array<double,N> &pdf) {
		return entropy(std::vector<double>(pdf.begin(), pdf.end()));
	}

	static inline std::vector<double> pdf(size_t N, Type type, double h) {

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
	}

	static inline std::array<double,256> pdf(Type type, double h) {

		auto P = pdf(256, type, h);
		std::array<double,256> A;
		for (size_t i=0; i<256; i++) A[i]=P[i];
		return A;
	}

	static inline std::vector<uint8_t> getResiduals(const std::vector<double> &pdf, size_t S) {

		int8_t cdf[0x10000];
		uint j=0;
		double lim=0;
		for (uint i=0; i<pdf.size(); i++) {
			lim += pdf[i]*0x10000;
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

	static inline std::vector<uint8_t> getResiduals(const std::array<double,256> &pdf, size_t S) {
		return getResiduals(std::vector<double>(pdf.begin(), pdf.end()), S);
	}

};
