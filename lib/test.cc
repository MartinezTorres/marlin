#include "distribution.hpp"

#include <iostream>

int main() {
	
		constexpr Distribution::Gaussian<256> test(3),test2(3);
		constexpr Distribution::Gaussian<65536> test3(3);
        std::cout << test[2] << std::endl;
}
