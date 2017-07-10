
#include <array>
#include <iostream>
#include <cmath>
#include <bitset>

using namespace std;
#include "distribution.h"

template<typename T, size_t N, typename Compare = std::less<T>>
static constexpr cx::array<T,N> getSorted2(cx::array<T,N> arr, Compare = std::less<T>()) {
    
    cx::sort(arr.begin(), arr.end(), Compare());
    return arr;
}
	template<uint64_t N> struct RequiredBits    { enum { value = 1 + RequiredBits<(N>>1)>::value }; };
	template<>           struct RequiredBits<0> { enum { value = 0 }; };
	template<>           struct RequiredBits<1> { enum { value = 1 }; };

	template<uint64_t Max, uint8_t requiredBits = RequiredBits<Max>::value>
	struct SmallestUint {typedef typename SmallestUint<Max, requiredBits+1>::type type; };	
	template<uint64_t Max> struct SmallestUint<Max, 8> { typedef uint8_t  type; };
	template<uint64_t Max> struct SmallestUint<Max,16> { typedef uint16_t type; };
	template<uint64_t Max> struct SmallestUint<Max,32> { typedef uint32_t type; };
	template<uint64_t Max> struct SmallestUint<Max,64> { typedef uint64_t type; };
	
	
int main() {

    
    constexpr double v = 0.20001;
    std::cout << (cx::exp(v)-std::exp(v))/std::exp(v) << std::endl;

    
    constexpr auto p = Distribution::Gaussian<256>(100000);
    //for (auto &&i : p) std::cout << i << std::endl;
    constexpr auto a  = Distribution::entropy(p);
    std::cout << a << std::endl;
    
    constexpr auto e = Distribution::getWithEntropy(Distribution::Gaussian<256>,.5);
    std::cout << Distribution::entropy(e) << std::endl;

    std::cout << std::log2(1e+307) << endl;
    std::cout << std::log2(1e+308) << endl;
    
    constexpr auto d = Distribution::getWithEntropy(Distribution::Gaussian<5>,.5);
    for (auto &&dd : d) {
        std::cout << dd << " "; 
    }
    std::cout << std::endl;

    constexpr auto k = getSorted2(d);

    for (auto &&dd : k) {
        std::cout << dd << " "; 
    }
    std::cout << std::endl;
    
    std::cout << sizeof(SmallestUint<255>::type) << std::endl;

    
/*
    constexpr Distribution::Gaussian<256> t1(3);
    constexpr Distribution::Laplace<256> t2(3);
    constexpr Distribution::Exponential<256> t3(3);
    constexpr Distribution::Poisson<256> t4(10);

    std::cout << Distribution::entropy(t4) << std::endl;
    
    
    
    
    constexpr double k = Distribution::constexpr_log(v);

    std::cout << (Distribution::constexpr_log(v)-std::log(v))/std::log(v) << std::endl;

    {    
        double v = 237847823;
        std::cout << Distribution::constexpr_log2(v) << " " << std::log2(v) << std::endl;

        std::cout << (Distribution::constexpr_log2(v)-std::log2(v))/std::log2(v) << std::endl;
    }*/

//    for (auto &t : t4)
//        std::cout << t << std::endl;
    
}
