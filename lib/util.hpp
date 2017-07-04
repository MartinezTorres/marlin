#pragma once
#include <cstddef>
#include <limits>

namespace {
namespace cx {
	
	//using std::array;
	
	
	// ARRAY
	template<typename T, size_t N> 
	class array {
	protected:
		T arr[N];
	public:
//		constexpr array() {}
		
//		constexpr array(const T &val) { for (auto &a : arr) a = val; }

//		constexpr array(const T(&a)[N]) { for (size_t n = 0; n<N; n++) arr[n] = a[n]; }
	
		constexpr void fill( const T& value ) { for (auto &a : arr) a = value; }
		
		constexpr T &operator[](size_t i)       { return arr[i]; }
		constexpr T  operator[](size_t i) const { return arr[i]; }

		static constexpr size_t size() { return N; }
		
		constexpr       T *begin()       { return &arr[0]; }
		constexpr const T *begin() const { return &arr[0]; }
		constexpr       T *end()         { return &arr[N]; }
		constexpr const T *end()   const { return &arr[N]; }
        
        constexpr array<T,N+1> operator+(const T &rho ) const {
            array<T,N+1> ret;
			for (size_t n = 0; n<N; n++)
				ret[n] = arr[n];
            ret[N] = rho;
            return ret;
        }

        template<size_t M> 
        constexpr array<T,N+M> append(const array<T,M> &rho) const {
            array<T,N+M> ret;
			for (size_t n = 0; n<N; n++)
				ret[n] = arr[n];
			for (size_t n = 0; n<M; n++)
				ret[N+n] = rho[n];
            return ret;            
        }
	};
    
	// ARRAY
	template<typename T> 
	class constvector {
	protected:
		T *arr = nullptr;
		size_t sz = 0;
		size_t cap = 0;
	public:
		
		constexpr T &&operator[](size_t i)       && { return arr[i]; }
		constexpr T   operator[](size_t i) const && { return arr[i]; }

		constexpr size_t size()     const && { return sz; }
		constexpr size_t capacity() const && { return capacity; }
		
		constexpr       T *begin()       { return &arr[0]; }
		constexpr const T *begin() const { return &arr[0]; }
		constexpr       T *end()         { return &arr[sz]; }
		constexpr const T *end()   const { return &arr[sz]; }
		
		constexpr constvector && operator+(const T &rho ) const && {
			
            array<T,N+1> ret;
			for (size_t n = 0; n<N; n++)
				ret[n] = arr[n];
            ret[N] = rho;
            return ret;
        }
	};

	// MATH
	
	constexpr double exp(double d) {

		if ( d >= 600) d = 600;
		if (-d >= 600) return std::numeric_limits<double>::min();
	
		if (d>1 or d<-1) return exp(d/2)*exp(d/2);

		double t=0, f=1;
		for (size_t i=1; i<1024; i++) {
			f*=d/i;
			if (t+f == t) break;
			t+=f;
		}
		return t+1.;
	}

	constexpr double log_2 = 0.693147180559945309417232121458176568075500134360255254120680009L;

	constexpr double log2(double d) {
		
		if ( d > 1 ) return -log2(1./d);
		
		int exp = 0;
		while (d < 1./(1<<30)) { exp -= 30; d*= (1<<30); }
		while (d < 1./(1<< 6)) { exp -=  6; d*= (1<< 6); }
		while (d < 1./(1<< 1)) { exp -=  1; d*= (1<< 1); }
		
		double x = 1 - d, r = 0,  xp = 1;
		for (size_t i=1; i<1024; i++) {
			double update = (xp *= x)/i;
			if (r-update == r) break;
			r -= update;
		}

		return exp + r*(1/log_2);
	}
	
	constexpr double log(double d) { return log2(d) * log_2; }

	template<typename T, size_t N> 
	constexpr array<T,N> norm1(const array<T,N> &a) { 
			
		array<T,N> ret = {};
		double sum = 0;
		for (size_t n=0; n<N; ++n) sum += a[n];
		for (size_t n=0; n<N; ++n) ret[n] = a[n]/sum;
		return ret;
	}

}
}
