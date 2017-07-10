
#include <array>
#include <iostream>
#include <cmath>
#include <bitset>

using namespace std;
#include "distribution.hpp"

template<typename T, size_t N, typename Compare = std::less<T>>
static constexpr cx::array<T,N> getSorted2(cx::array<T,N> arr, Compare = std::less<T>()) {
    
    cx::sort(arr.begin(), arr.end(), Compare());
    return arr;
}

    template<uint64_t N> struct bitcount    { enum { value = 1 + bitcount<(N>>1)>::value }; };
	template<>           struct bitcount<0> { enum { value = 0 }; };
	template<>           struct bitcount<1> { enum { value = 1 }; };
    template<uint64_t N> struct bytecount   { enum { value = (bitcount<N>::value + 7) / 8 }; };

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
    
    std::cout << bytecount<0>::value << std::endl;
    std::cout << bytecount<1>::value << std::endl;
    std::cout << bytecount<2>::value << std::endl;
    std::cout << bytecount<3>::value << std::endl;
    std::cout << bytecount<4>::value << std::endl;
    std::cout << bytecount<5>::value << std::endl;
    std::cout << bytecount<6>::value << std::endl;
    std::cout << bytecount<7>::value << std::endl;
    std::cout << bytecount<8>::value << std::endl;
    std::cout << bytecount<16>::value << std::endl;
    std::cout << bytecount<32>::value << std::endl;
    std::cout << bytecount<0x100000000LL>::value << std::endl;

    
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
