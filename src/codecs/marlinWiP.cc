#include <marlinlib/marlin.hpp>

#include <codecs/marlinWiP.hpp>
#include <util/distribution.hpp>

#include <functional>
#include <string>
#include <iostream>
#include <cassert>
#include <cstring>
#include <stack>
#include <queue>
#include <map>
#include <bitset>
#include <unordered_map>
#include <algorithm>
#include <memory>

struct MarlinWiPPimpl : public CODEC8Z {
	
	std::vector<std::shared_ptr<MarlinWiP::SingleDictionaryCodec>> dictionaries;
	
	std::string coderName;
	std::string name() const { return coderName; }
	
	MarlinWiPPimpl(Distribution::Type distType, size_t overlap, size_t numDict) {

		{
			std::ostringstream oss;
			oss << "MarlinWiP " << (distType==Distribution::Laplace?"Lap:":"Exp:") <<  ":" << overlap << ":" << numDict;
			coderName = oss.str();
		}

		std::vector<std::shared_ptr<MarlinWiP::SingleDictionaryCodec>> builtDictionaries(numDict);

//		#pragma omp parallel for
		for (size_t p=0; p<numDict; p++) {
			
			std::vector<double> pdf(256,0.);
			for (double i=0.05; i<0.99; i+=0.1) {
				
				auto pdf0 = Distribution::pdf(distType, (p+i)/numDict);
				for (size_t j=0; j<pdf.size(); j++)
					pdf[j] += pdf0[j]/10.;
			}
			std::map<std::string, double> conf; conf["O"]=overlap;
			builtDictionaries[p] = std::make_shared<MarlinWiP::SingleDictionaryCodec>(pdf, conf );
		}
		
		dictionaries.resize(256);
		
//		#pragma omp parallel for
		for (size_t h=0; h<256; h+=4) {
			
			auto testData = Distribution::getResiduals(Distribution::pdf(distType, (h+2)/256.), 1<<16);
			
			double lowestSize = testData.size()*0.99; // If efficiency is not enough to compress 1%, skip compression
			for (auto &&dict : builtDictionaries) {
				std::string out;
				dict->encode(testData, out);
				if (out.size() < lowestSize) {
					lowestSize = out.size();
					for (size_t hh = 0; hh<4; hh++)
						dictionaries[h+hh] = dict;
				}
			}	
		}
	}

	
	void   compress(
		const std::vector<std::reference_wrapper<const AlignedArray8>> &in,
		      std::vector<std::reference_wrapper<      AlignedArray8>> &out,
		      std::vector<std::reference_wrapper<      uint8_t      >> &entropy) const { 
		
		for (size_t i=0; i<in.size(); i++)
			if (dictionaries[entropy[i]])
				dictionaries[entropy[i]]->encode(in[i].get(), out[i].get());
			else
				out[i].get().resize(in[i].get().size());
	}

	void uncompress(
		const std::vector<std::reference_wrapper<const AlignedArray8>> &in,
		      std::vector<std::reference_wrapper<      AlignedArray8>> &out,
		      std::vector<std::reference_wrapper<const uint8_t      >> &entropy) const {
		
		for (size_t i=0; i<in.size(); i++)
			if (dictionaries[entropy[i]])
				dictionaries[entropy[i]]->decode(in[i].get(), out[i].get());
			else
				out[i].get().resize(in[i].get().size());
	}

};


MarlinCodecWiP::MarlinCodecWiP(Distribution::Type distType, size_t overlap, size_t numDict) 
	: CODEC8withPimpl( new MarlinWiPPimpl(distType, overlap, numDict) ) {}

