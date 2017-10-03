#include <iostream>

#include <marlinlib/marlin.hpp>

using std::cin;
using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char **argv) {

	// Parse command line optins
	std::map<std::string,double> options;
	options["--overlap"]=2;
	options["--maxWordSize"]=1024;

	for (int i=1; i<argc; i++) {
		
		std::string name;
		while (isalpha(*argv[i]) or *argv[i]=='-') name += *argv[i]++;
		options[name] = 1; 
		if (*argv[i] and *argv[i]=='=') 
			options[name] = atof(argv[i]+1);
	}
	
	for (auto &op : options)
		cerr << op.first << " " << op.second << endl;
	
	auto pdf = Distribution::pdf(256,Distribution::Laplace,0.25);

	for (size_t keySize = 9; keySize<17; keySize++) {
		std::cout << keySize << " " << 
			Marlin2018Simple(pdf,keySize,0,options["--maxWordSize"]).efficiency << " " <<
			Marlin2018Simple(pdf,keySize,1,options["--maxWordSize"]).efficiency << " " <<
			Marlin2018Simple(pdf,keySize,2,options["--maxWordSize"]).efficiency << " " << 
			Marlin2018Simple(pdf,keySize,3,options["--maxWordSize"]).efficiency << " " <<
			Marlin2018Simple(pdf,keySize,4,options["--maxWordSize"]).efficiency << std::endl;
	}
	return 0;
}
