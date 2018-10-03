#include "marlin.h"
#include "distribution.hpp"
#include <iostream>

static void printAlpha(std::vector<uint8_t> msg) {
	for (size_t i=0; i<msg.size(); i++)
		printf("%c", char('a'+msg[i]));
		printf("\n");
}

static void printHex(std::vector<uint8_t> msg) {
	for (size_t i=0; i<msg.size(); i++)
		printf("%02X", msg[i]);
		printf("\n");
}

static bool testMini() {
	
	std::cout << "Test Mini" << std::endl;
	
	size_t sz = 1<<9;
	
	auto distribution = Distribution::pdf(8, Distribution::Laplace, 0.5);
	for (int i=0; i<1000; i++)
		std::swap(distribution[rand()%8], distribution[rand()%8]);
	
	std::vector<uint8_t> original(Distribution::getResiduals(distribution,sz));
	std::vector<uint8_t> compressed(sz);
	std::vector<uint8_t> uncompressed(sz);
	
	MarlinDictionary::Configuration conf;
	conf["K"] = 8;
	conf["O"] = 0;
	//conf["debug"] = 99;
	conf["purgeProbabilityThreshold"] = 1e-99;	
	
	MarlinDictionary dict("test", distribution, conf);
	
	dict.compress(original, compressed);
	dict.decompress(compressed, uncompressed);
	
	std::cout << "Compressed Size: " << compressed.size() << std::endl;

	printAlpha(original);
	printHex(compressed);
	printAlpha(uncompressed);
	
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
		
		std::cout << "Get Dictionary!" << std::endl;		
		MarlinDictionary dict("test", Distribution::pdf(256, Distribution::Laplace, p));
		std::cout << "Compress:" << std::endl;
		dict.compress(original, compressed);
		std::cout << "Compressed to:" << double(compressed.size()) / original.size() << std::endl;
		std::cout << "Decompress:" << std::endl;
		//dict.decompress(compressed, uncompressed);
		std::cout << "Done!" << std::endl;
		
		
		
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
			
			//return false;
		} else {
			
			
			std::cout  << "P: " << p << " " << "Original == uncompressed!" << " " << compressed.size() << std::endl;
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
