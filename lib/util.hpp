#include <string>
#include <limits>
#include <cmath>
#include <ieee754.h>
#include <numeric>

namespace {
namespace cx {
	
	// MATH
	
	constexpr long double exp(long double d) {

		if ( d >= 600) d = 600;
		if (-d >= 600) return std::numeric_limits<double>::min();
	
		if (d>1 or d<-1) return exp(d/2)*exp(d/2);

		long double t=0, f=1;
		for (size_t i=1; i<20; i++) {
			f*=d/i;
			t+=f;
		}
		return t+1.;
	}

	constexpr long double log_2 = 0.693147180559945309417232121458176568075500134360255254120680009L;

	constexpr long double log2(long double d) {
		
		if ( d < std::numeric_limits<double>::min() ) return std::numeric_limits<double>::lowest();
		if ( d > 1 ) return -log2(1./d);
		
		int exp = 0;
		while (d < 1./(1<<30)) { exp -= 30; d*= (1<<30); }
		while (d < 1./(1<< 6)) { exp -=  6; d*= (1<< 6); }
		while (d < 1./(1<< 1)) { exp -=  1; d*= (1<< 1); }
		
		long double x = 1 - d, r = 0,  xp = 1;
		for (size_t i=1; i<1024 and xp>std::numeric_limits<double>::min(); i++)
			r -= (xp *= x)/i;

		return exp + r*(1/log_2);
	}
	
	constexpr long double log(long double d) { return log2(d) * log_2; }

	// ARRAY

	template<size_t N> 
	class array {
	protected:
		double arr[N];
	public:		
		constexpr array() : arr() {
			for (auto &a : arr) a = std::numeric_limits<double>::min();
		}

		constexpr array(const double(&a)[N]) : arr() {
			for (size_t n = 0; n<N; n++)
				arr[n] = a[n];
		}
		
		constexpr double &operator[](size_t i)       { return arr[i]; }
		constexpr double  operator[](size_t i) const { return arr[i]; }

		static constexpr size_t size() { return N; }
		
		constexpr       double *begin()       { return &arr[0]; }
		constexpr const double *begin() const { return &arr[0]; }
		constexpr       double *end()         { return &arr[N]; }
		constexpr const double *end()  const { return &arr[N]; }
	
		constexpr array norm1() const { 
			
			array a;
		    double sum = 0;
			for (size_t n=N; n; n--) sum += arr[n-1];
			for (size_t n=N; n; n--) a[n-1] = arr[n-1]/sum;
			return a;
		}
	};
}
}
