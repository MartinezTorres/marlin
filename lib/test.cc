#include "distribution.hpp"

#include <iostream>

int main() {
	
    constexpr Distribution::Gaussian<256> t1(3);
    constexpr Distribution::Laplace<256> t2(3);
    constexpr Distribution::Exponential<256> t3(3);
    constexpr Distribution::Poisson<256> t4(10);

    std::cout << Distribution::entropy(t4) << std::endl;
    
    
    
    constexpr double v = 0.000000000237847823;
    std::cout << Distribution::constexpr_log(v) << " " << std::log(v) << std::endl;
    
    constexpr double k = Distribution::constexpr_log(v);

    std::cout << (Distribution::constexpr_log(v)-std::log(v))/std::log(v) << std::endl;

    {    
        double v = 237847823;
        std::cout << Distribution::constexpr_log2(v) << " " << std::log2(v) << std::endl;

        std::cout << (Distribution::constexpr_log2(v)-std::log2(v))/std::log2(v) << std::endl;
    }

//    for (auto &t : t4)
//        std::cout << t << std::endl;
    
}
