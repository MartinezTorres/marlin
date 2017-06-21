//#pragmonce
#include <vector>
#include <array>
#include <string>
#include <limits>
#include <numeric>

#include <constexpr/src/include/cx_math.h>
#include <constexpr/src/include/cx_array.h>


class Distribution {
public:
	
	constexpr static long double constexpr_exp(long double d) {

		if ( d >= 600) d = 600;
		if (-d >= 600) return std::numeric_limits<double>::min();
	
		if (d>1 or d<-1) return constexpr_exp(d/2)*constexpr_exp(d/2);

		long double t=0, f=1;
		for (size_t i=1; i<20; i++) {
			f*=d/i;
			t+=f;
		}
		return t+1.;
	}


	constexpr static long double constexpr_log(long double d) {
		return constexpr_log2(d) * 0.693147180559945309417232121458176568075500134360255254120680009;
	}

	constexpr static long double constexpr_log2(long double d) {
        
        if ( d < std::numeric_limits<double>::min() ) return std::numeric_limits<double>::lowest();
        int exp = ((ieee854_long_double *)&d)->ieee.exponent - IEEE854_LONG_DOUBLE_BIAS+1;
        ((ieee854_long_double *)&d)->ieee.exponent = IEEE854_LONG_DOUBLE_BIAS-1;
        
        long double x = 1 - d, r = 0,  xp = 1;
		for (size_t i=1; i<1024 and xp>std::numeric_limits<double>::min(); i++)
            r -= (xp *= x)/i;

        return exp + r*1.442695040888963407359924681001892137426645954152985934135449407;
    }

	template<size_t N> 
	class constexpr_array {
	protected:
		double arr[N];
	public:
		constexpr constexpr_array() : arr() {
			for (size_t n = 0; n<N; n++)
				arr[n] = std::numeric_limits<double>::min();
		};

		constexpr constexpr_array(const constexpr_array &src) : arr() { 
			for (size_t n = 0; n<N; n++)
				arr[n] = src[n];
		}
		
		constexpr double operator[](size_t i) const { return arr[i]; }

		constexpr size_t size() const { return N; }
        
        constexpr const double *begin() const { return &arr[0]; }
        constexpr const double *end() const { return &arr[N]; }
	};
	
	template<size_t N> 
	struct Gaussian : constexpr_array<N> {
		using constexpr_array<N>::arr;
		
		constexpr Gaussian(double b) : constexpr_array<N>() {
			for (int64_t i=-10*N+1; i<10*N; i++)
				arr[(10*N+i) % N] += constexpr_exp(-i*i/b);
			
			double sum = 0;
			for (size_t n=N; n; n--) sum += arr[n-1];
			for (size_t n=N; n; n--) arr[n-1] /= sum;
		};
	};

	template<size_t N> 
	struct Laplace : constexpr_array<N> {
		using constexpr_array<N>::arr;
		
		constexpr Laplace(double b) : constexpr_array<N>() {

            arr[0] += 1.;
            for (int64_t i=1; i<10*N; i++) {
                arr[      i  % N] += constexpr_exp(-i/b );
                arr[(10*N-i) % N] += constexpr_exp(-i/b );
            }
			
			double sum = 0;
			for (size_t n=N; n; n--) sum += arr[n-1];
			for (size_t n=N; n; n--) arr[n-1] /= sum;
		};
	};

	template<size_t N> 
	struct Exponential : constexpr_array<N> {
		using constexpr_array<N>::arr;
		
		constexpr Exponential(double b) : constexpr_array<N>() {

            arr[0] += 1.;
            for (int64_t i=1; i<10*N; i++)
                arr[      i  % N] += constexpr_exp(-i/b );
			
			double sum = 0;
			for (size_t n=N; n; n--) sum += arr[n-1];
			for (size_t n=N; n; n--) arr[n-1] /= sum;
		};
	};
    
    template<size_t N> 
	struct Poisson : constexpr_array<N> {
		using constexpr_array<N>::arr;
		
		constexpr Poisson(double l) : constexpr_array<N>() {

            double lp=0,kf=0;
            for (int64_t k=0; k<10*N and kf<1e300; k++) {
                kf = (k?kf*k:1);
                lp = (k?lp*l:1);
                arr[      k  % N] += lp*constexpr_exp(-l)/kf;
            }
			
			double sum = 0;
			for (size_t n=N; n; n--) sum += arr[n-1];
			for (size_t n=N; n; n--) arr[n-1] /= sum;
		};
	};
	
    template<size_t N> 
    static constexpr long double entropy(constexpr_array<N> const &pmf) {

		long double distEntropy=0;
        for (auto &p : pmf)
            if (p)
                distEntropy += -p*constexpr_log2(p);
		return distEntropy;        
    }
/*	
	

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
*/
};

