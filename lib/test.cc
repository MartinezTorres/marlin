
#include <array>
#include <iostream>

using namespace std;
#include "distribution.hpp"

int main() {

    constexpr double v = 0.20001;
    std::cout << (cx::exp(v)-std::exp(v))/std::exp(v) << std::endl;

    
    constexpr auto p = Distribution::Gaussian<256>(100000);
    //for (auto &&i : p) std::cout << i << std::endl;
    constexpr auto a  = Distribution::entropy(Distribution::Gaussian<256>(10000));
    std::cout << a << std::endl;

    auto e = Distribution::getWithEntropy(Distribution::Gaussian<256>,.5);
    std::cout << Distribution::entropy(e) << std::endl;

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
