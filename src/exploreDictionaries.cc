#include <queue>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>
#include <stack>
#include <list>
#include <map>

#include <util/distribution.hpp>

using std::cin;
using std::cout;
using std::cerr;
using std::endl;

struct Marlin2Dictionary {
	
	size_t overlapping = 0;

	struct Word {
		
		std::string symbols = "";
		double p = 0;
		uint8_t state = 0;
	};
	
	std::vector<Word> dictionary;
	
	struct Node;
	struct Node : std::vector<std::shared_ptr<Node>> {
		double p=0;
		size_t sz=0;
		size_t erased=0;
	};

	std::shared_ptr<Node> buildTree(const std::vector<double> &P, const std::vector<double> &Pstates, size_t dictSize, bool prune) {

//		prune = false;
		
		std::vector<double> PN = P;
		for (size_t i=P.size()-1; i; i--)
			PN[i-1] += PN[i];

		std::vector<double> Pchild = P;
		for (size_t i=0; i<P.size(); i++)
			Pchild[i] = P[i]/PN[i];


//		auto cmp = [](const std::shared_ptr<Node> &lhs, const std::shared_ptr<Node> &rhs) { return lhs->p < rhs->p;};

		auto cmp = [](const std::shared_ptr<Node> &lhs, const std::shared_ptr<Node> &rhs) { 
			
			return lhs->p*(1+std::pow(lhs->sz,1)) < rhs->p*(1+std::pow(rhs->sz,1));
			//if (std::abs(lhs->p - rhs->p) > 0.0000001)
			//	return lhs->p < rhs->p;
			//return lhs->sz < rhs->sz;
		};

		std::priority_queue<std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>, decltype(cmp)> pq(cmp);
		auto pushAndPrune = [&pq,prune](std::shared_ptr<Node> node) {
			if (not prune or node->p>1e-10)
				pq.push(node);
			else
				node->erased = true;
		};

		// DICTIONARY INITIALIZATION
		std::shared_ptr<Node> root = std::make_shared<Node>();
		
		// Include empty?
		if (prune) {
			pq.push(root);
		} else {
			root->erased = true;
		}

		for (size_t c=0; c<P.size(); c++) {			
				
			root->push_back(std::make_shared<Node>());
			double sum = 0;
			for (size_t t = 0; t<=c; t++) sum += Pstates[t]/PN[t];
			root->back()->p = sum * P[c];
			root->back()->sz = 1;
			pushAndPrune(root->back());
		}
			
		// DICTIONARY GROWING
		while (pq.size()<dictSize) {
				
			auto node = pq.top();
			pq.pop();
				
			double p = node->p * Pchild[node->size()];
			node->push_back(std::make_shared<Node>());
			node->back()->p = p;
			root->back()->sz = node->sz+1;
			node->p -= p;
			pushAndPrune(node->back());
				
			if (node->size()<P.size()-1) {

				pushAndPrune(node);
					
			} else {
				node->erased = true;
				node->p = 0;

				node->push_back(std::make_shared<Node>());
				node->back()->p = node->p;
				root->back()->sz = node->sz+1;
				pushAndPrune(node->back());
			}
		}
		
		return root;
	}
	
	std::vector<Word> buildWords(const std::shared_ptr<Node> &root) {
	
		std::vector<Word> ret;
		
		std::stack<std::pair<std::shared_ptr<Node>, Word>> q;
		q.emplace(root, Word());
		while (not q.empty()) {
			std::shared_ptr<Node> n = q.top().first;
			Word w = q.top().second;
			q.pop();
			if (not n->erased) ret.push_back(w);
			for (size_t i = 0; i<n->size(); i++) {
				
				Word w2 = w;
				w2.symbols += char(i);
				w2.p = n->at(i)->p;
				w2.state = n->at(i)->size();
				q.emplace(n->at(i), w2);
			}
		}
		//std::cerr << "NW: " << ret.size() << std::endl;
		return ret;
	}
	
	std::vector<Word> arrangeWords(const std::vector<Word> &words) {
		
		std::vector<Word> sortedWords = words;
		auto cmp = [](const Word &lhs, const Word &rhs) { 
			if (lhs.state != rhs.state) return lhs.state<rhs.state;
			return lhs.p > rhs.p;
		};
		std::sort(sortedWords.begin(), sortedWords.end(), cmp);
		
		std::vector<Word> ret = sortedWords;
		for (size_t i=0,j=0,k=0; i<sortedWords.size(); i++,j+=(1<<overlapping)) {
			if (j>=ret.size()) j=++k;
			ret[j] = sortedWords[i];
		}
		return ret;
	}
	
	std::vector<Word> arrangeWordsAlt(const std::vector<Word> &words) {
		
		std::vector<Word> sortedWords = words;
		auto cmp = [](const Word &lhs, const Word &rhs) { 
			return lhs.p < rhs.p;
		};
		std::sort(sortedWords.begin(), sortedWords.end(), cmp);
		
		std::vector<Word> ret = sortedWords;
		for (size_t i=0,j=0,k=0; i<sortedWords.size(); i++,j+=(1<<overlapping)) {
			if (j>=ret.size()) j=++k;
			ret[j] = sortedWords[i];
		}
		return ret;
	}
	
	template<typename T>
	std::vector<T> concatenate( const std::vector<std::vector<T>> &W) {

		std::vector<T> A;
		for (auto &&w : W)
			A.insert(A.end(), w.begin(), w.end());
		
		countUnique(A);		
		return A;
	}

	void countUnique( std::vector<Word> W ) {
		
		std::vector<std::string> C;
		for (auto &&w :W) C.push_back(w.symbols);
		std::sort(C.begin(), C.end());
		std::cerr << "Total: " << C.size() << " Unique: " << unique(C.begin(),C.end())-C.begin() << std::endl;
	}

	void print(std::vector<Word> dictionary) {
		if (dictionary.size()>40) return;

		for (size_t i=0; i<dictionary.size()/(1<<overlapping); i++) { 
			
			for (size_t k=0; k<(1U<<overlapping); k++) {
				
				auto idx = i + (k* (dictionary.size()/(1<<overlapping)));
				auto &&w = dictionary[idx];
				printf(" %02lX %01ld %2d %01.3lf ",idx,i%(1<<overlapping),w.state,w.p);
				for (size_t j=0; j<16; j++) putchar("0123456789ABCDEF "[j<w.symbols.size()?w.symbols[j]:16]);
			}
			putchar('\n');
		}		
		putchar('\n');
	}

	void print(std::vector<std::vector<double>> Pstates) {
		
		for (size_t i=0; i<Pstates[0].size() and i<4; i++) { 
			
			printf("S: %02ld",i);
			for (size_t k=0; k<Pstates.size(); k++) 
					 printf(" %01.3lf",Pstates[k][i]);
			putchar('\n');
		}		
		putchar('\n');
	}


	Marlin2Dictionary(const std::vector<double> &P, size_t dictSize, size_t tries=7, size_t overlapping=0) 
		: overlapping(overlapping) {
		
		std::vector<std::vector<double>> Pstates;
		for (auto k=0; k<(1<<overlapping); k++) {
			std::vector<double> PstatesSingle(P.size(), 0.);
			PstatesSingle[0] = 1./(1<<overlapping);
			Pstates.push_back(PstatesSingle);
		}
		
		int leastProbable = 0;
		
		std::vector<std::vector<Word>> dictionaries;
		for (auto k=0; k<(1<<overlapping); k++)
			dictionaries.push_back( arrangeWords( buildWords( buildTree(P, Pstates[k], dictSize/(1<<overlapping), k!=leastProbable) ) ) );
			
		dictionary = concatenate(dictionaries);
			
		print(dictionary);
			
		while (tries--) {

			// UPDATING STATE PROBABILITIES
			{
				for (auto k=0; k<(1<<overlapping); k++)
					Pstates[k] = std::vector<double>(P.size(), 0.);

				for (size_t i=0; i<dictionary.size(); i++)
					Pstates[i%(1<<overlapping)][dictionary[i].state] += dictionary[i].p;
			}
			
			// Find least probable subdictionary
			{
				double minP = 1.1;
				for (size_t i=0; i<Pstates.size(); i++) {
					double sumProb = 0.;
					for (auto &&ps : Pstates[i])
						sumProb += ps;
					if (sumProb > minP) continue;
					minP = sumProb;
					leastProbable = i;
				}
			}
			print(Pstates);

			for (auto k=0; k<(1<<overlapping); k++)
				dictionaries[k] = arrangeWords( buildWords( buildTree(P, Pstates[k], dictSize/(1<<overlapping), k!=leastProbable) ) );
			
			dictionary = concatenate(dictionaries);
			
			print(dictionary);

			averageBitsPerSymbol(P,dictSize);
		}
		
		//test(P,dictSize);
		//averageBitsPerSymbolEmpirically(P,dictSize);
		
		std::ofstream off(std::string("probabilities") + "01234"[overlapping] + ".txt");
		for (auto &&w : dictionary)
			off << w.p << std::endl;
	}
		
	double averageBitsPerSymbol(const std::vector<double> &P, size_t dictSize) const {
	
		double meanLength = 0;
		for (auto &&w : dictionary)
			meanLength += w.p * w.symbols.size();

		//printf("Meanlength: %3.2lf\n", meanLength);
		printf("Compress Ratio: %3.4lf\n", (std::log2(dictSize)-overlapping)/(meanLength*std::log2(P.size())));
		printf("Efficiency: %3.4lf\n", 
			Distribution::entropy(P)/8./
			((std::log2(dictSize)-overlapping)/(meanLength*std::log2(P.size()))));

		return 0.;
	}
	
	double averageBitsPerSymbolEmpirically(const std::vector<double> &P, size_t dictSize) const {
		
		static const size_t TEST_SIZE = 1<<18;
		
		std::string testData(TEST_SIZE,0);
		// create data
		{
			size_t i=0;
			double ap=0;
			for (size_t j=0; j<P.size(); j++) {
				while (i<(ap+P[j])*testData.size())
					testData[i++]=j;
				ap += P[j];
			}
		}
		std::random_shuffle(testData.begin(), testData.end());
		
		// pad with 0s
		for (size_t i=0; i<10; i++)
			testData.push_back(0);

		auto &&W = dictionary;
		double meanLengthE = 0;
		size_t nWords = 0;
		size_t lastWord = 0;
		for (size_t i=0, longest=0; i<TEST_SIZE; i+=longest) {
			
			size_t best = 0;
			longest = 0;
			for (size_t j=0; j<W.size()/(1<<overlapping); j++) {
				auto idx = (lastWord%(1<<overlapping))*(W.size()/(1<<overlapping)) + j;
				auto &&w = W[idx].symbols;
				if (w.size()>longest and testData[i] == w[0] and testData[i+w.size()-1] == w[w.size()-1] and testData.compare(i,w.size(),w)==0 ) {
					best = idx;
					longest = w.size();
				}
			}
			lastWord = best;
			meanLengthE += longest;
			nWords++;
		}

		printf("Empirical Meanlength: %3.2lf\n", meanLengthE/nWords);
		printf("Empirical Compress Ratio: %3.4lf\n", 
			(nWords*(std::log2(dictSize)-overlapping))/(double(TEST_SIZE)*std::log2(P.size())));


		return 0.;
	}

	double test(const std::vector<double> &P, size_t dictSize) const {
		
		static const size_t TEST_SIZE = 1<<20;
		
		std::string testData(TEST_SIZE,0);
		// create data
		{
			size_t i=0;
			double ap=0;
			for (size_t j=0; j<P.size(); j++) {
				while (i<(ap+P[j])*testData.size())
					testData[i++]=j;
				ap += P[j];
			}
		}
		std::random_shuffle(testData.begin(), testData.end());
		
		auto encoded = encode(testData);
		auto decoded = decode(encoded);
		std::cout << testData.size() << " " << decoded.size() << std::endl;
		printf("Test ok: %s\n", (decoded==testData?"true":"false"));
		

//		printf("Compression Ratio: %2.3lf\n", );

//		return double(nWords*std::log2(dictSize))/TEST_SIZE;
		return 0.0;
	}

	std::vector<uint16_t> encode(std::string testData) const {

		std::string reconstructed;

		std::vector<uint16_t> ret;
		auto W = dictionary;
		size_t lastWord = 0;
		for (size_t i=0, longest=0; i<testData.size(); i+=longest) {
			
			size_t remaining = testData.size()-i;
			
			
			size_t best = 0;
			longest = 0;
			for (size_t j=0; j<W.size()/2; j++) {
				auto &&w = W[(lastWord%2)*(W.size()/2) + j].symbols;
				if (w.size()>longest and w.size()<= remaining and testData.substr(i,w.size()).compare(w)==0 ) {
					best = (lastWord%2)*(W.size()/2) + j;
					longest = w.size();
				}
			}
			ret.push_back(best);
			reconstructed += W[best].symbols;
			
			//std::cout << best << " " << i << " " << reconstructed.size() << std::endl;
			lastWord = best;
		}
		return ret;
	}
	
	std::string decode(std::vector<uint16_t> data) const {

		std::string ret;
		auto W = dictionary;
		for (auto &&d : data)
			ret += W[d].symbols;

		return ret;
	}


};

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
	options["--tries"]=3;
	options["--size"]=256;
//	options["--maxWordSize"]=256;

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
	
		cerr << "P: "; for (auto p : P) cerr << p << " "; cerr  << endl;


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
		
	cerr << "P: "; for (auto p : P) cerr << p << " "; cerr << Distribution::entropy(P) << endl;
		
	// Ensure that symbols are sorted in order of decreasing probabilities
	std::sort(P.begin(), P.end(), std::greater<double>());
	
	
	std::cerr << "Marlin" << std::endl;
//	MarlinDictionary(P,options["--size"],options["--tries"]);
	std::cerr << "Tunstall" << std::endl;
	//TunstallDictionary(P,options["--size"]);
	std::cerr << "Marlin2" << std::endl;
	Marlin2Dictionary(P,options["--size"],options["--tries"]);
	Marlin2Dictionary(P,options["--size"]*2,options["--tries"],1);
	Marlin2Dictionary(P,options["--size"]*4,options["--tries"],2);
	Marlin2Dictionary(P,options["--size"]*8,options["--tries"],3);
	Marlin2Dictionary(P,options["--size"]*16,options["--tries"],4);
	Marlin2Dictionary(P,options["--size"]*32,options["--tries"],5);
		
	return 0;
}
