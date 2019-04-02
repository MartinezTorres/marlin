/***********************************************************************

Marlin: A Fast Entropy Codec

MIT License

Copyright (c) 2018 Manuel Martinez Torres

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

***********************************************************************/

#include <marlin.h>

using namespace marlin;

template<typename TSource, typename MarlinIdx>
std::map<std::string, double> TMarlinDictionary<TSource,MarlinIdx>::updateConf( 
	const std::vector<double> &sourceAlphabet, 
	Configuration conf) {
	
	conf.emplace("K",8);
	conf.emplace("O",4);
	
	conf.emplace("debug",1);
//	conf.emplace("purgeProbabilityThreshold",1e-99);
//	conf.emplace("purgeProbabilityThreshold",1e-6);
//	conf.emplace("purgeProbabilityThreshold",0.5/4096/32);
	conf.emplace("purgeProbabilityThreshold",0.5/4096/32);
	conf.emplace("iterations",3);
//	conf.emplace("minMarlinSymbols", std::max(1U<<size_t(conf.at("O")),2U));
	conf.emplace("minMarlinSymbols", 2U);
	conf.emplace("maxMarlinSymbols",(1U<<size_t(conf.at("K"))));
	conf["maxMarlinSymbols"] = std::min(conf["maxMarlinSymbols"], double((1U<<size_t(conf.at("K")))));
	conf["maxMarlinSymbols"] = std::min(conf["maxMarlinSymbols"], double((1U<<(8*sizeof(MarlinIdx)))-1));

	double maxWordSize = conf["maxWordSize"];
	conf["maxWordSize"] = 255;

	double sourceEntropy = TMarlinDictionary<TSource,MarlinIdx>::calcSourceEntropy(sourceAlphabet);

	if (not conf.count("shift")) {
		
		double best = 0;
		size_t shift = 0;
		for (int i=0; i<6; i++) {
			
			auto testConf = conf;
			testConf["shift"] = shift;
			double test = TMarlinDictionary<TSource,MarlinIdx>(sourceAlphabet, testConf).compressionRatio;
			if (test > 1.0001*best) {
				best = test;
				conf = testConf;
			};
			shift ++;
		}
	}

	conf["maxWordSize"] = maxWordSize;
	
	if (conf["maxWordSize"]==0) {
		
		conf.emplace("autoMaxWordSize",64);
		
		double best = 0.;
		size_t sz = 4;
		while (sz <= conf["autoMaxWordSize"]) {
			
			auto testConf = conf;
			testConf["maxWordSize"] = sz-1;
			double test = TMarlinDictionary<TSource,MarlinIdx>(sourceAlphabet, testConf).compressionRatio;
			if (test > 1.0001*best) {
				best = test;
				conf = testConf;
			} else break;
			sz*=2;
		}
	}
	
	//printf("%lf %lf\n", conf["maxWordSize"], conf["shift"]);
	
	return conf;
}


////////////////////////////////////////////////////////////////////////
//
// Explicit Instantiations
#include "instantiations.h"
INSTANTIATE(TMarlinDictionary)
//typedef std::map<std::string, double> phonyMap; // Commas do not fit well within macros
//INSTANTIATE_MEMBER(TMarlinDictionary, updateConf(const std::vector<double> &sourceAlphabet, Configuration conf) -> phonyMap)
