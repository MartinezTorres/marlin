#include <iostream>

#include <marlinlib/marlin.hpp>

using std::cin;
using std::cout;
using std::cerr;
using std::endl;

void usage() {
	
	cout << "Usage: testDictionary [options]" << endl;
	cout << "Options:" << endl;
	cout << "   -h or --help: show this help" << endl;
	cout << "   --custom  <N-1 probabilities>: use a custom distribution" << endl;
	cout << "   --exp  <size> <entropy>: use an exponential distribution with <size> alphabet and <entropy> entropy" << endl;
	cout << "   --lap  <size> <entropy>: use an Laplace distribution with <size> alphabet and <entropy> entropy" << endl;
	cout << "   --norm <size> <entropy>: use an normal distribution with <size> alphabet and <entropy> entropy" << endl;
	exit(-1);
}


int main(int argc, char **argv) {

	if (argc==1) usage();
	
	// Parse command line optins
	std::map<std::string,double> options;
	options["--keySize"]=12;
	options["--overlap"]=2;
	options["--maxWordSize"]=127;
	options["--testSize"]=1<<20;

	std::vector<double> P;
	for (int i=1; i<argc; i++) {
		if (argv[i][0]=='-') {
			std::string name;
			while (isalpha(*argv[i]) or *argv[i]=='-') name += *argv[i]++;
			options[name] = 1; 
			if (*argv[i] and *argv[i]=='=') 
				options[name] = atof(argv[i]+1);
		} else {
			P.push_back(atof(argv[i]));
		}
	}
	
	for (auto &op : options)
		cerr << op.first << " " << op.second << endl;
	
	if (options["-h"] or options["--help"] or options.empty())
		usage();
	
	//cerr << "P: "; for (auto p : P) cerr << p << " "; cerr  << endl;


	if (options["--custom"]) {
		// Adding probability for the last symbol, ensuring that all probabilities sum 1
		{
			double lp=1;
			for (auto p : P)
				lp -= p;
			P.push_back(lp);	
		}

		// Check that all probabilities are positive
		for (auto p : P)
			if (p<0)
				usage();
				
	} else if (options["--exp"]) {
		
		if (P.size()!=2 or P[0] < 1 or P[1]<0 or P[1]>1) usage(); 
		P = Distribution::pdf(P[0],Distribution::Exponential,P[1]);
	
	} else if (options["--lap"]) {
		
		if (P.size()!=2 or P[0] < 1 or P[1]<0 or P[1]>1) usage(); 
		P = Distribution::pdf(P[0],Distribution::Laplace,P[1]);
	} else {
		usage();
	}
		

	if (true) {
		std::cerr << "MarlinWiP" << std::endl;
		MarlinWiP::SingleDictionaryCodec(P).benchmark(P,options["--testSize"]);
	}


	if (true) {
		std::cerr << "MarlinDCC2018" << std::endl;
		Marlin2018Simple::setConfiguration("debug",true);
		Marlin2018Simple(P,options["--keySize"],options["--overlap"],options["--maxWordSize"]).benchmark(P,options["--testSize"]);
	}
		
	return 0;
}
