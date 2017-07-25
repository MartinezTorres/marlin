#pragma once
#include <cstddef>
#include <functional>
#include <limits>

namespace {
namespace cx {
	
	using std::array; // Can be used when array is finally constexpr compliant
	
	template<typename T>
	constexpr inline void swap(T &a, T &b) { T c = std::move(a); a = std::move(b); b = std::move(c); }
	
	// ARRAY
	template<typename T, size_t N> 
	class arrayCx {
	public:
		T arr[N] = {};
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
    };
    
	// VECTOR
	template<typename T, size_t C> 
	class vector {
	protected:
		T arr[C] = {};
        size_t sz = 0;
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
        
        void push_back2(const T &item) { std::cout << C << " " << sz << std::endl; arr[sz++] = item; }
        constexpr void push_back(const T &item) { arr[sz++] = item; }
        constexpr void pop_back() { if (sz!=C) arr[sz] = T(); sz--; }
	};
    
    // PRIORITY_QUEUE
    template<typename T, size_t C, typename Compare = std::less<T>>
    class minmaxHeap {
    protected:
        vector<T,C> container = {};
        
        constexpr size_t parent(size_t pos) { return (pos-1)>>1; };
        constexpr size_t left  (size_t pos) { return (pos<<1)+1; };
        constexpr size_t right (size_t pos) { return (pos<<1)+2; };

        constexpr size_t level (size_t pos) { return  pos?1+level(parent(pos)):0; };
        
        
        constexpr bool test(size_t a, size_t b, bool reverse) {
			if (reverse) swap(a,b);
			return a < size() and b<size() and Compare()(container[a], container[b]);
		}
        
		constexpr void sink(size_t pos) { 
			
			bool found = true;
            while (found) {

				found = false;
				size_t newPos = pos, testPos = 0;
               
				testPos = left (pos);        if (test(testPos, newPos, level(testPos) & 1)) newPos=testPos;
                testPos = right(pos);        if (test(testPos, newPos, level(testPos) & 1)) newPos=testPos;
                testPos = left (left (pos)); if (test(testPos, newPos, level(testPos) & 1)) newPos=testPos;
                testPos = left (right(pos)); if (test(testPos, newPos, level(testPos) & 1)) newPos=testPos;
                testPos = right(left (pos)); if (test(testPos, newPos, level(testPos) & 1)) newPos=testPos;
                testPos = right(right(pos)); if (test(testPos, newPos, level(testPos) & 1)) newPos=testPos;

				if (newPos != pos) {
					swap(container[pos], container[newPos]);
					pos = newPos;				
					found = true;
				}
			}
        }
        
        constexpr void bubble(size_t pos) {
			
            for ( size_t p = pos; p;  p = parent(p) ) {
				
				if ( test(p,pos, level(p) & 1) ) {
					
					swap(container[pos], container[p]);
					pos = p;
				}
            }
		}


    public:
    
        constexpr void push(const T&  item) {

			container.push_back(item);
			bubble(container.size()-1);
        }
        
        constexpr void popMin() {
            
            if (container.empty()) return;
            
            swap(container.back(), container.front());
            container.pop_back();
            sink(0);
		}


        constexpr void popMax() {
            
            if (container.empty()) return;
            
            if (size()<2) {
				container.pop_back();
			} else if (Compare()(container[1], container[2])) {
				
				swap(container.back(), container[1]);
				container.pop_back();
				sink(1);
			} else {
				swap(container.back(), container[2]);
				container.pop_back();
				sink(2);
			}
        }

        constexpr const T min()   const { return container.front(); }
        constexpr const T max()   const { 
			
			if (size()<2) return container.back();
			if (Compare()(container[1], container[2])) {
				return container[1];
			} else {
				return container[2];
			}
		}
				

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
        
        for (size_t sz = 0; first + sz != last; sz++) {

			auto item = std::move(first[sz]);
            size_t pos = sz;
            size_t parent = (pos-1)>>1;
            while (pos and compare(first[parent], item)) {
				
                first[pos] = std::move(first[parent]);
                pos = parent;
                parent = (pos-1)>>1;
            }
            first[pos] = item;
        }
                    
        // Empty heap
        for (size_t sz = last - first; sz>1; sz--) {

			size_t pos = sz-1, newPos = 0;            
            while (pos!=newPos) {

				swap(first[pos], first[newPos]);
				pos = newPos;				
                size_t l = (pos<<1)+1, r = (pos<<1)+2;
                if (l < sz-2 and Compare()(first[newPos], first[l])) newPos = l;
                if (r < sz-2 and Compare()(first[newPos], first[r])) newPos = r;
			}
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
/*
	template<typename T, size_t N> 
	constexpr array<T,N> norm1(const array<T,N> &a) { 
			
		array<T,N> ret = {};
		double sum = 0;
		for (size_t n=0; n<N; ++n) sum += a[n];
		for (size_t n=0; n<N; ++n) ret[n] = a[n]/sum;
		return ret;
	}*/
}
}
