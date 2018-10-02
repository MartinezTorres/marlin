#include "marlin.h"
#include "distribution.hpp"
#include <iostream>

static bool testMini() {
	
	size_t sz = 1<<5;
	
	auto distribution = Distribution::pdf(8, Distribution::Laplace, 0.5);
	for (int i=0; i<1000; i++)
		std::swap(distribution[rand()%8], distribution[rand()%8]);
	
	std::vector<uint8_t> original(Distribution::getResiduals(distribution,sz));
	std::vector<uint8_t> compressed(sz);
	std::vector<uint8_t> uncompressed(sz);
	
	MarlinDictionary::Configuration conf;
	conf["K"] = 4;
	conf["O"] = 0;
//	conf["debug"] = 99;
	conf["purgeProbabilityThreshold"] = 1e-99;
	
	
	
	MarlinDictionary dict("test", distribution, conf);
	
	dict.compress(original, compressed);
	dict.decompress(compressed, uncompressed);
	
	std::cout << "Compressed Size: " << compressed.size() << std::endl;

	for (size_t i=0; i<original.size(); i++)
		printf("%c", char('a'+original[i]));
		printf("\n");

	for (size_t i=0; i<compressed.size(); i++)
		printf("%02X", compressed[i]);
		printf("\n");

	
	if (original != uncompressed) {
		
		std::cout << "P: " << 0.5 << " " << "FAIL!     sizes(" << original.size() << "," << uncompressed.size() << ")" << std::endl;
		for (size_t i=0; i<10; i++)
			printf("%02X:%02X ", original[i], uncompressed[i]);
		std::cout << std::endl;
		
		{
			int c = 0;
			for (size_t i=0; i<original.size(); i++) {
				if (original[i] != uncompressed[i]) {
					printf("Pos %04X = %02X:%02X\n", uint(i), original[i], uncompressed[i]);
					if (c++==4) break;
				}
			}
		}
		
		return false;
	} else {
		
		
		std::cout << "Original == uncompressed!" << " " << compressed.size() << std::endl;
	}
	return compressed.size() < 0.55*uncompressed.size();
}

static bool testLaplace() {
	
	size_t sz = 1<<20;
	
	for (double p=0.1; p<.995; p+=0.1) {
		
		std::vector<uint8_t> original(Distribution::getResiduals(Distribution::pdf(Distribution::Laplace, p),sz));
		std::vector<uint8_t> compressed(sz);
		std::vector<uint8_t> uncompressed(sz);
		
		std::cout << "IO!" << std::endl;
		
		MarlinDictionary dict("test", Distribution::pdf(256, Distribution::Laplace, p));
		
		dict.compress(original, compressed);
		dict.decompress(compressed, uncompressed);
		
		
		
		if (original != uncompressed) {
			
			std::cout << "P: " << p << " " << "FAIL!     sizes(" << original.size() << "," << uncompressed.size() << ")" << std::endl;
			for (size_t i=0; i<10; i++)
				printf("%02X:%02X ", original[i], uncompressed[i]);
			std::cout << std::endl;
			
			{
				int c = 0;
				for (size_t i=0; i<original.size(); i++) {
					if (original[i] != uncompressed[i]) {
						printf("Pos %04X = %02X:%02X\n", uint(i), original[i], uncompressed[i]);
						if (c++==4) break;
					}
				}
			}
			
			return false;
		} else {
			
			
			std::cout << "Original == uncompressed!" << " " << compressed.size() << std::endl;
		}
	}
	return true;
}


int main() {

	return 
		testMini() and
		testLaplace() and
		true?0:-1;
}
