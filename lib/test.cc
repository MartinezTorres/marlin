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

    auto dist = Distribution::getWithEntropy(Distribution::Gaussian<256>,1./8);
/*	
	auto marlinV1  = MarlinV1( dist );
    */

using namespace std;

int main() {

//    cout << "SIZE: " << sizeof(marlinV1) << " " << sizeof(marlinV1) << endl;
    
	
	auto in = getResiduals( dist, 1<<20);
/*	ibitstream ibs1(in.data(), in.size());
	
	
	std::vector<uint8_t> compressed, decompressed;
	compressed.resize(in.size()*2);
	obitstream obs1(compressed.data(), compressed.size());
	
	marlinV1.encode(ibs1, obs1); obs1.sync();

	std::cout << in.size() << " " << obs1.size() << " " << obs1.size()*8./in.size() << std::endl;
    
    obs1.sync();
    compressed.resize(obs1.size());
    ibitstream ibs2(compressed.data(), compressed.size());

	decompressed.resize(in.size()*2);
	obitstream obs2(decompressed.data(), decompressed.size());
    
    marlinV1.decode(ibs2, obs2); obs2.sync();
    
	std::cout << in.size() << " " << obs2.size() << std::endl;
    
    

    for (size_t i = 0, count = 10; i<in.size() and count; ++i) {
		if (in[i] != decompressed[i]) {
			std::cerr << "i6 " << i << " " << uint64_t(in[i]) << " " << uint64_t(decompressed[i]) << std::endl;
            count--;
        }
	}
	
	
	depq<int> pq;
	for (auto &&e : {1,2,234,32,1,2,3}) pq.insert(e);
	for (auto && e: pq) std::cout << e << " ";
	std::cout << std::endl;
	*/
    std::cout << "Test2" << std::endl;
    
    std::vector<uint8_t> compressed, decompressed;
    compressed.resize(4*in.size());
    size_t sz1 = MarlinEncode(in.data(), in.size(), compressed.data(), compressed.size());
    compressed.resize(sz1);

    decompressed.resize(4*in.size());
    size_t sz2 = MarlinDecode(compressed.data(), compressed.size(), decompressed.data(), decompressed.size());
    decompressed.resize(sz2);
    
    std::cout << in.size() << " " << compressed.size() << " " << decompressed.size() << std::endl;
    
    for (size_t i = 0, count = 10; i<in.size() and count; ++i) {
		if (in[i] != decompressed[i]) {
			std::cerr << "AA " << i << " " << i%4096 << " " << uint64_t(in[i]) << " " << uint64_t(decompressed[i]) << std::endl;
            count--;
        }
	}
}
