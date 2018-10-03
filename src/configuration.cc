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

template<typename TSource, typename MarlinIdx>
std::map<std::string, double> TMarlinDictionary<TSource,MarlinIdx>::updateConf( 
	const std::vector<double> &sourceAlphabet, 
	TMarlinDictionary<TSource,MarlinIdx>::Configuration conf) {
	
	conf.emplace("K",8);
	conf.emplace("O",2);
	
	conf.emplace("debug",1);
	conf.emplace("purgeProbabilityThreshold",1e-99);
	conf.emplace("iterations",5);
	conf.emplace("minMarlinSymbols", std::max(1U<<size_t(conf.at("O")),8U));
	conf.emplace("maxMarlinSymbols",(1U<<size_t(conf.at("K")))-1);
		
	if (not conf.count("shift")) {
		conf["shift"] = 0;
		double best = MarlinDictionary("", sourceAlphabet, conf).efficiency;
		for (int s=1; s<6; s++) {
			conf["shift"] = s;
			double e = MarlinDictionary("", sourceAlphabet, conf).efficiency;
			if (e<=best) {
				conf["shift"] = s-1;
				break;
			}
			best = e;
		}
	}
	
	if (not conf.count("maxWordSize")) {
		conf["maxWordSize"] = 15;
		double e15 = MarlinDictionary("", sourceAlphabet, conf).efficiency;
		conf["maxWordSize"] = 7;
		double e7 = MarlinDictionary("", sourceAlphabet, conf).efficiency;
		conf["maxWordSize"] = 3;
		double e3 = MarlinDictionary("", sourceAlphabet, conf).efficiency;
		if (e7>1.0001*e3) {
			conf["maxWordSize"] = 7;
		}
		if (e15>1.0001*e7) {
			conf["maxWordSize"] = 15;
		}
	}
	
	//printf("%lf %lf\n", conf["maxWordSize"], conf["shift"]);
	
	return conf;
}
