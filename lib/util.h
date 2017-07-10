#pragma once
#include <cstddef>
#include <functional>
#include <limits>

namespace {
namespace cx {
	
	// ARRAY
	template<typename T, size_t N> 
	class array {
	public:
		T arr[N];
	public:

		constexpr void fill( const T& value ) { for (auto &a : arr) a = value; }
		
		constexpr T &operator[](size_t i)       { return arr[i]; }
		constexpr T  operator[](size_t i) const { return arr[i]; }

		static constexpr size_t size() { return N; }
		
		constexpr       T* begin()       { return &arr[0]; }
		constexpr const T* begin() const { return &arr[0]; }
		constexpr       T& front()       { return  arr[0]; }
		constexpr const T& front() const { return  arr[0]; }
		constexpr       T* end()         { return &arr[N]; }
		constexpr const T* end()   const { return &arr[N]; }
		constexpr       T& back()        { return  arr[N-1]; }
		constexpr const T& back()  const { return  arr[N-1]; }
		
		constexpr static cx::array<T,N> mulLine(
			size_t i,
			const cx::array<cx::array<T,N>,N> M) {
		
			cx::array<T,N> L = {};
			
			for (size_t j=0; j<N; ++j)
				for (size_t k=0; k<N; ++k)
					L.arr[j] += M.arr[i].arr[k] * M.arr[k].arr[j];

			return L;
		}		
    };

    template<typename T, size_t N>
    struct Liner {
		
		cx::array<T,N> L = {};
		
		constexpr Liner( size_t i, const cx::array<cx::array<T,N>,N> MM) {

			L = array<T,N>::mulLine(i, MM);
		}
		
	};

    template<typename T, size_t N>
    struct Squarer {
		
		cx::array<cx::array<T,N>,N> M = {};
		cx::array<T,N> L = {};
		
		constexpr Squarer( const cx::array<cx::array<T,N>,N> MM) {

			for (size_t i=0; i<N; i++) {
				L = Liner<T,N>(i, MM).L;
				M[i] = L;
			}
		}
		
	};
	
	

    
	// VECTOR
	template<typename T, size_t C> 
	class vector {
	protected:
		T arr[C];
        size_t sz;
	public:
    		
		constexpr T &operator[](size_t i)       { return arr[i]; }
		constexpr T  operator[](size_t i) const { return arr[i]; }

		constexpr size_t size()     const { return sz; }
		constexpr bool   empty()    const { return sz==0; }
		constexpr size_t capacity() const { return C; }

		// Todo: unallocate unused objects
        constexpr void   resize(size_t newSz) { sz = newSz; }
        constexpr void   clear() { sz = 0; }
		
		constexpr       T* begin()        { return &arr[0]; }
		constexpr const T* begin()  const { return &arr[0]; }
		constexpr       T& front()        { return  arr[0]; }
		constexpr const T& front()  const { return  arr[0]; }
		constexpr       T* end()          { return &arr[sz]; }
		constexpr const T* end()    const { return &arr[sz]; }
		constexpr       T& back()         { return  arr[sz-1]; }
		constexpr const T& back()   const { return  arr[sz-1]; }
        
        constexpr void push_back(const T &item) { arr[sz++] = item; }
        constexpr void pop_back() { if (sz!=C) arr[sz] = T(); sz--; }
	};
    
    // MATRIX
	template<typename T, size_t N1, size_t N2> 
	class matrix {
	protected:
		T arr[N1][N2];

		template<size_t N3>
		constexpr static double vecmul(size_t i, size_t j, const matrix<T,N2,N3> lhs, const matrix<T,N2,N3> rhs) {

			double ret = 0.0;
			for (size_t k=0; k<N2; ++k)
				ret += lhs.arr[i][k] * rhs.arr[k][j];
			return ret;
		}	

	public:

//		constexpr void fill( const T& value ) { for (auto &a : arr) a = value; }
		
		constexpr       T* operator[](size_t i)       { return arr[i]; }
		constexpr const T* operator[](size_t i) const { return arr[i]; }
		
		
		template<size_t N3>
		constexpr matrix<T,N1,N3> operator*(matrix<T,N2,N3> rhs) const {
			
			matrix<T,N1,N3> ret = {};
			for (size_t i=0; i<N1; ++i)
				for (size_t j=0; j<N3; ++j)
					ret.arr[i][j] = vecmul(i,j,*this,rhs);
			return ret;
		}
    };
        
    // PRIORITY_QUEUE
    template<typename T, size_t C, typename Compare = std::less<T>>
    class priority_queue {
    protected:
        vector<T,C> container = {};

    public:
    
        constexpr void push(const T&  item) {

            size_t pos = container.size();
            container.resize(pos + 1);
            
            size_t parent = (pos-1)>>1;
            while (pos and Compare()(container[parent], item)) {
                container[pos] = std::move(container[parent]);
                pos = parent;
                parent = (pos-1)>>1;
            }
            container[pos] = std::move(item);
        }
        
        constexpr void pop() {
            
            if (container.empty()) return;
            
            auto && item = container[size()-1];
            size_t pos = 0;
            while (true) {
                size_t l = (pos<<1)+1;
                if (l < size()-1 and not Compare()(container[l], item)) {
                    container[pos] = std::move(container[l]);
                    pos = l;
                    continue;
                }
                
                size_t r = (pos<<1)+2;
                if (r < size()-1 and not Compare()(container[r], item)) {
                    container[pos] = std::move(container[r]);
                    pos = r;
                    continue;
                }
                break;
            }
            container[pos] = std::move(item);
            container.pop_back();
        }

        constexpr const T& top()   const { return container.front(); }
        constexpr size_t   size()  const { return container.size();  }
        constexpr const T* begin() const { return container.begin(); }
		constexpr const T* end()   const { return container.end();   }       
    };
    
    // SORT ALGORITHM
    template<
        typename RandomAccessIterator, 
        typename Compare = std::less<typename RandomAccessIterator::value_type> 
    >
    constexpr void sort(
        RandomAccessIterator first, 
        RandomAccessIterator last,
        Compare compare = Compare()) {
        
        // Create heap
        size_t sz = 0;
        while (first + sz != last) {

            size_t pos = sz;
            size_t parent = (pos-1)>>1;
            while (pos and compare(first[parent], first[sz])) {
                first[pos] = std::move(first[parent]);
                pos = parent;
                parent = (pos-1)>>1;
            }
            first[pos] = first[sz];
            sz++;
        }
            
        // Empty heap
        while (sz) {

            auto larger = std::move(first[0]);
            auto && item = first[sz-1];
            size_t pos = 0;
            while (true) {
                size_t l = (pos<<1)+1;
                if (l < sz-1 and not compare(first[l], item)) {
                    first[pos] = std::move(first[l]);
                    pos = l;
                    continue;
                }
                
                size_t r = (pos<<1)+2;
                if (r < sz-1 and not compare(first[r], item)) {
                    first[pos] = std::move(first[r]);
                    pos = r;
                    continue;
                }
                break;
            }
            first[pos] = std::move(item);
            first[sz-1] = std::move(larger);
            sz--;
        } 
    }
    
      
    
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
