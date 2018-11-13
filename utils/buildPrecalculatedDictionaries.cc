#include <marlin.h>
#include <distribution.hpp>
#include <sstream>

static  void buildDictionaries(
	std::map<std::string, std::shared_ptr<Marlin>> &dictionaries,
	size_t numDict, 
	Distribution::Type type, 
	std::map<std::string, double> conf = std::map<std::string, double>()
	) {


	for (size_t p=0; p<numDict; p++) {
		
		std::ostringstream oss;
		oss << "Marlin";
		if      (type == Distribution::Gaussian   ) oss << "_Gaussian";
		else if (type == Distribution::Laplace    ) oss << "_Laplace";
		else if (type == Distribution::Exponential) oss << "_Exponential";
		oss << "_" << (p<10?"0":"") << p << "_" << numDict; 
		for (auto &&c : conf)
			oss << "_" << c.first << "_" << c.second;
			
		std::vector<double> pdf(256,0.);
		
		size_t nSamples = 10;
		for (double i=0.5/nSamples; i<0.9999999; i+=1./nSamples) {
			auto pdf0 = Distribution::pdf(type, (p+i)/double(numDict));
			for (size_t j=0; j<pdf.size(); j++)
				pdf[j] += pdf0[j]/nSamples;
		}

		dictionaries[oss.str()] = std::make_shared<Marlin>(oss.str(),pdf,conf);
	}
}


int main() {

	auto &&out = std::cout;
	
	std::map<std::string, std::shared_ptr<Marlin>> builtDictionaries;
	
	buildDictionaries(builtDictionaries,16,Distribution::Laplace);
	buildDictionaries(builtDictionaries,16,Distribution::Gaussian);
	buildDictionaries(builtDictionaries,16,Distribution::Exponential);
		
	out << "#include <marlin.h>" << std::endl;
	for (auto &&dict : builtDictionaries) {
		
		
		out << "static const std::array<uint8_t, 256> prebuilt_dictionary_" << dict.first << "_source2marlin = {";
		for (auto &&p: dict.second->source2marlin) out << uint64_t(p) << ","; out << "};" << std::endl; 

		out << "static const uint32_t prebuilt_dictionary_" << dict.first << "_compressorTableVector[] = {";
		for (auto &&p: *dict.second->compressorTableVector) out << uint64_t(p) << ","; out << "};" << std::endl; 

		out << "static const uint32_t prebuilt_dictionary_" << dict.first << "_compressorTableInitVector[] = {";
		for (auto &&p: *dict.second->compressorTableInitVector) out << uint64_t(p) << ","; out << "};" << std::endl; 

		out << "static const uint8_t prebuilt_dictionary_" << dict.first << "_decompressorTableVector[] = {";
		for (auto &&p: *dict.second->decompressorTableVector) out << uint64_t(p) << ","; out << "};" << std::endl; 
		
		out << "static const Marlin prebuilt_dictionary_" << dict.first << "(" << std::endl;
		out << "    \"" << dict.second->name << "\", // name" << std::endl; 
		out << "    " << dict.second->K << ", // K" << std::endl; 
		out << "    " << dict.second->O << ", // O" << std::endl; 
		out << "    " << dict.second->shift << ", // shift" << std::endl; 
		out << "    " << dict.second->maxWordSize << ", // maxWordSize" << std::endl; 
		out << "    " << dict.second->efficiency << ", // Efficiency" << std::endl; 
		out << "    " << uint64_t(dict.second->unrepresentedSymbolToken) << ", // unrepresentedSymbolToken" << std::endl; 

		out << "    prebuilt_dictionary_" << dict.first << "_source2marlin, " << std::endl; 
		out << "   &prebuilt_dictionary_" << dict.first << "_compressorTableVector[0], " << std::endl; 
		out << "   &prebuilt_dictionary_" << dict.first << "_compressorTableInitVector[0], " << std::endl; 
		out << "   &prebuilt_dictionary_" << dict.first << "_decompressorTableVector[0], " << std::endl; 
		
		out << "    " << uint32_t(dict.second->marlinMostCommonSymbol) << ", // marlinMostCommonSymbol" << std::endl; 
		out << "    " << dict.second->isSkip << " // isSkip" << std::endl; 
		out << ");" << std::endl;
	};
	
	out << "static const Marlin *Marlin_all_prebuilt_dictionaries[] = {" << std::endl;
	for (auto &&dict : builtDictionaries)
		out << "    &prebuilt_dictionary_" << dict.first << ", " << std::endl;
	out << "nullptr};" << std::endl;

	out << "const Marlin **Marlin_get_prebuilt_dictionaries() {" << std::endl;
	out << "    return Marlin_all_prebuilt_dictionaries;" << std::endl;
	out << "}" << std::endl;
	
}
