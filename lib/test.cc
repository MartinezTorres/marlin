#include <iostream>

#include <fstream>

#include "marlin.cc"
#include "distribution.h"

#include <iostream>
#include <vector>
template<typename F>
inline std::vector<uint8_t> getResiduals(const F &pmf, size_t S) {

	int8_t cdf[0x10000];
	uint j=0;
	double lim=0;
	for (uint i=0; i<pmf.size(); i++) {
		lim += pmf[i]*0x10000;
		uint ilim = round(lim);
		while (j<ilim)
			cdf[j++]=i;
	}

	std::vector<uint8_t> ret;
	uint32_t rnd =  135154;
	for (size_t i=0; i<S; i++) {
		rnd = 36969 * (rnd & 65535) + (rnd >> 16);
		ret.push_back(cdf[rnd&0xFFFF]);
	}

	return ret;
}

    auto dist = Distribution::getWithEntropy(Distribution::Gaussian<256>,2./8);
	
	auto dictionary  = Dictionary<256,15,4096>( dist );

int main() {
    
	
	auto in = getResiduals( dist, 1<<20);
	ibitstream ibs1(in.data(), in.size());
	
/*	for (size_t i = 0; i<in.size(); ++i) {
		if (ibs1.read(8) != in[i])
			std::cerr << "ii " << i << std::endl;		
	}*/
		
	
	
	std::vector<uint8_t> compressed, decompressed;
	compressed.resize(in.size()*2);
	obitstream obs1(compressed.data(), compressed.size());
	
	dictionary.encode(ibs1, obs1); obs1.sync();

	std::cout << in.size() << " " << obs1.size() << " " << obs1.size()*8./in.size() << std::endl;
    
    obs1.sync();
    compressed.resize(obs1.size());
    ibitstream ibs2(compressed.data(), compressed.size());

	decompressed.resize(in.size()*2);
	obitstream obs2(decompressed.data(), decompressed.size());
    
    dictionary.decode(ibs2, obs2); obs2.sync();
    
	std::cout << in.size() << " " << obs2.size() << std::endl;
    
    

    for (size_t i = 0, count = 10; i<in.size() and count; ++i) {
		if (in[i] != decompressed[i]) {
			std::cerr << "i6 " << i << " " << uint64_t(in[i]) << " " << uint64_t(decompressed[i]) << std::endl;
            count--;
        }
	}
}
