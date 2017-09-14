#include <queue>
#include <iostream>
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


class Dictionary {
protected:

	struct Node;
	struct Node : std::vector<std::shared_ptr<Node>> {
		double p=0;
	};

	std::shared_ptr<Node> root;

public:

	const std::vector<double>  P;
	const size_t dictSize;

	Dictionary(const std::vector<double> &P, size_t dictSize) :
		P(P), dictSize(dictSize) {}
		
	void print( const std::shared_ptr<Node> &node, const std::string &s = "" ) const {
		
		if (node==root) {
			
			if (dictSize>9) return;
			cout << endl;
		}
		
		if (node->size() < P.size())
			cout << s << " " << node->p << endl;

		for (size_t i=0; i<node->size(); i++)
			print( (*node)[i], s+char('A'+i) );
	}

	double meanLength(const std::shared_ptr<Node> &node, size_t length = 0) const {

		double ret = 0;
		
		for (auto &c : *node)
			ret += meanLength(c, length+1);
		
		if (node->size() < P.size())
			ret += node->p * double(length);
			
		return ret;		
	}
	
	double averageBitsPerSymbol() const {
			
		return std::log2(dictSize)/meanLength(root);		
	}
	
	std::list<std::vector<uint8_t>> getWords(const std::shared_ptr<Node> &node, std::vector<uint8_t> prefix = std::vector<uint8_t>()) const {
		
		std::list<std::vector<uint8_t>> ret;
		if (node->size() < P.size())
			ret.push_back(prefix);
		
		for (size_t i=0; i<node->size(); i++) {
			prefix.push_back(i);
			ret.splice(ret.end(), getWords((*node)[i], prefix));
			prefix.pop_back();
		}
		return ret;
	}

	double averageBitsPerSymbolEmpirically() const {
		
		static const size_t TEST_SIZE = 1<<20;
		
		std::vector<uint8_t> testData(TEST_SIZE);
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

		auto W = getWords(root);
		double meanLengthE = 0;
		size_t nWords = 0;
		for (size_t i=0, longest=0; i<TEST_SIZE; i+=longest) {
			
			std::vector<uint8_t> const* lw = nullptr;
			for (auto &w : W)
				if (w[0] == testData[i])
					if (w[w.size()-1] == testData[i+w.size()-1])
						if (w==std::vector<uint8_t>(&testData[i], &testData[i+w.size()]))
							if (lw == nullptr or w.size()>lw->size())
								lw = &w;
								
			longest = lw->size();
			meanLengthE += lw->size();
			nWords++;
		}

		printf("Empirical Meanlength: %3.2lf\n", meanLengthE/nWords);
		printf("Empirical Compress Ratio: %3.4lf\n", (nWords*std::log2(dictSize))/(double(TEST_SIZE)*std::log2(P.size())));

		return double(nWords*std::log2(dictSize))/TEST_SIZE;
	}

};

struct MarlinDictionary : public Dictionary {

	MarlinDictionary(const std::vector<double> &P, size_t dictSize, size_t tries=3) : Dictionary(P, dictSize) {

		std::vector<double> Pstate(P.size(), 0.);
		Pstate[0] = 1.;

//		std::vector<double> Pstate(P.size(), 1./(P.size()-1));
//		Pstate.back() = 0.;

		std::vector<double> PN = P;
		for (size_t i=P.size()-1; i; i--)
			PN[i-1] += PN[i];

		std::vector<double> Pchild = P;
		for (size_t i=0; i<P.size(); i++)
			Pchild[i] = P[i]/PN[i];
	
		while (tries--) {
			
			auto cmp = [](const std::shared_ptr<Node> &lhs, const std::shared_ptr<Node> &rhs) { return lhs->p < rhs->p;};
			std::priority_queue<std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>, decltype(cmp)> pq(cmp);

			// DICTIONARY INITIALIZATION
			root = std::make_shared<Node>();

			for (size_t c=0; c<P.size(); c++) {
				
				root->push_back(std::make_shared<Node>());

				double sum = 0;
				for (size_t t = 0; t<=c; t++) sum += Pstate[t]/PN[t];

				root->back()->p = sum * P[c];
				
				pq.push(root->back());
			}
			
			print(root);
			
			// DICTIONARY GROWING
			while (pq.size()<dictSize) {
				
				auto node = pq.top();
				pq.pop();
				
				double p = node->p * Pchild[node->size()];
				node->push_back(std::make_shared<Node>());
				node->back()->p = p;
				node->p -= p;
				pq.push(node->back());
				
				if (node->size()<P.size()-1) {

					pq.push(node);
				} else {

					node->push_back(std::make_shared<Node>());
					node->back()->p = node->p;
					node->p = 0;
					pq.push(node->back());
				}
				
				print(root);
			}

			
			
			auto oldPstate = Pstate;
			// UPDATING STATE PROBABILITIES
			{
				std::vector<std::vector<double>> T(P.size(),std::vector<double>(P.size(),0.));
		
				for (size_t w0 = 0; w0<P.size(); w0++) {
					
					std::stack<std::shared_ptr<Node>> st;
					st.push((*root)[w0]);

					double sum = 0;
					for (size_t t = 0; t<=w0; t++) sum += Pstate[t]/PN[t];

					
					while (not st.empty()) {
					
						auto node = st.top();
						st.pop();
						
						for (auto &c : *node) st.push(c);
						
						if (node->size()<P.size())
							for (size_t s=0; s<=w0; s++)
								T[s][node->size()] += node->p /(sum * PN[s]);
					}
				}
								
				int t = 10;
				double diff = 0;
				do {
					
					auto T2 = T;
					for (size_t i=0; i<P.size(); i++) {
						for (size_t j=0; j<P.size(); j++) {
							T2[i][j]=0;
							for (size_t k=0; k<P.size(); k++)
								T2[i][j] += T[i][k] * T[k][j];
						}
					}
					diff = 0;
					for (size_t i=0; i<P.size(); i++)
						for (size_t j=0; j<P.size(); j++)
							diff += std::abs(T[i][j]-T2[i][j]);

					//std::cerr << diff << std::endl;
					T = T2;

				} while (t-- and diff>.00001);
				
				Pstate = T[0];
			}

			// UPDATE NODE PROBABILITIES
			{
				cerr << "ST: "; for (auto ps : Pstate) cerr << ps << " ";
				cerr  << averageBitsPerSymbol() << " " << meanLength(root) << endl;

				for (size_t w0 = 0; w0<P.size(); w0++) {
					
					std::stack<std::shared_ptr<Node>> st;
					st.push((*root)[w0]);

					double oldsum = 0, newsum = 0;
					for (size_t t = 0; t<=w0; t++) {
						oldsum += oldPstate[t]/PN[t];
						newsum +=    Pstate[t]/PN[t];
					}

					while (not st.empty()) {
					
						auto node = st.top();
						st.pop();						
						for (auto &c : *node) st.push(c);
						
						if (node->size()<P.size())
							node->p *= newsum / oldsum;
					}
				}
				cerr << "ST: "; for (auto ps : Pstate) cerr << ps << " ";
				cerr  << averageBitsPerSymbol() << " " << meanLength(root) << " " << averageBitsPerSymbolEmpirically() << endl;

			}
		}
	}
};

struct TunstallDictionary : public Dictionary {

	TunstallDictionary(const std::vector<double> &P, size_t dictSize) : Dictionary(P, dictSize) {
			
		auto cmp = [](const std::shared_ptr<Node> &lhs, const std::shared_ptr<Node> &rhs) { return lhs->p < rhs->p;};
		std::priority_queue<std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>, decltype(cmp)> pq(cmp);

		// DICTIONARY INITIALIZATION
		root = std::make_shared<Node>();

		for (size_t c=0; c<P.size(); c++) {
				
			root->push_back(std::make_shared<Node>());
			root->back()->p = P[c];
			pq.push(root->back());
		}
			
		print(root);
			
		// DICTIONARY GROWING
		while (pq.size()+P.size()-1 <= dictSize) {
			
			auto node = pq.top();
			pq.pop();
			
			for (size_t c=0; c<P.size(); c++) {
			
				node->push_back(std::make_shared<Node>());
				node->back()->p = node->p * P[c];
				pq.push(node->back());
			}
		}

		print(root);
		cerr  << averageBitsPerSymbol() << " " << meanLength(root) << endl;
	}
};

// NAH! We simply have more states!!! MATH!
struct Marlin2Dictionary {

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
	};

	std::shared_ptr<Node> buildTree(const std::vector<double> &P, const std::vector<double> &Pstates, size_t dictSize) {

		std::vector<double> PN = P;
		for (size_t i=P.size()-1; i; i--)
			PN[i-1] += PN[i];

		std::vector<double> Pchild = P;
		for (size_t i=0; i<P.size(); i++)
			Pchild[i] = P[i]/PN[i];


		auto cmp = [](const std::shared_ptr<Node> &lhs, const std::shared_ptr<Node> &rhs) { return lhs->p < rhs->p;};
		std::priority_queue<std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>, decltype(cmp)> pq(cmp);

		// DICTIONARY INITIALIZATION
		std::shared_ptr<Node> root = std::make_shared<Node>();

		for (size_t c=0; c<P.size(); c++) {
				
			root->push_back(std::make_shared<Node>());
			double sum = 0;
			for (size_t t = 0; t<=c; t++) sum += Pstates[t]/PN[t];
			root->back()->p = sum * P[c];
			root->back()->sz = 1;
			pq.push(root->back());
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
//			if (node->back()->sz==1 or node->back()->p > 1e-5)
				pq.push(node->back());
				
			if (node->size()<P.size()-1) {

//				if (node->sz==1 or node->p > 1e+10)
					pq.push(node);
					
			} else {

				node->push_back(std::make_shared<Node>());
				node->back()->p = node->p;
				root->back()->sz = node->sz+1;
				node->p = 0;
//				if (node->back()->sz==1 or node->back()->p > 1e-5)
					pq.push(node->back());
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
			if (n->p > 0.) ret.push_back(w);
			for (size_t i = 0; i<n->size(); i++) {
				
				Word w2 = w;
				w2.symbols += char(i);
				w2.p = n->at(i)->p;
				w2.state = n->at(i)->size();
				q.emplace(n->at(i), w2);
			}
		}
		std::cerr << "NW: " << ret.size() << std::endl;
		return ret;
	}
	
	std::vector<Word> arrangeWords(const std::vector<Word> &words) {
		
		std::vector<Word> sortedWords = words;
		auto cmp = [](const Word &lhs, const Word &rhs) { 
			if (lhs.state != rhs.state) return lhs.state<rhs.state;
			return lhs.symbols.size() < rhs.symbols.size();
		};
		std::sort(sortedWords.begin(), sortedWords.end(), cmp);
		
		std::vector<Word> ret = sortedWords;
		for (size_t i=0,j=0; i<sortedWords.size(); i++,j+=2) {
			if (j>=ret.size()) j=1;
			ret[j] = sortedWords[i];
		}
		return ret;
	}
	
	template<typename T>
	std::vector<T> concatenate( std::vector<T> A, std::vector<T> B) {

		A.insert(A.end(), B.begin(), B.end());
		return A;
	}

	void print(std::vector<Word> dictionary) {
		if (dictionary.size()>40) return;

		for (size_t i=0; i<dictionary.size()/2; i++) { 
			{
				auto &&w = dictionary[i]; 
				printf(" %02lX %01ld %2d %01.3lf ",i,i%2,w.state,w.p);
				for (size_t j=0; j<16; j++) putchar("0123456789ABCDEF "[j<w.symbols.size()?w.symbols[j]:16]);
			}
			{
				auto &&w = dictionary[i+dictionary.size()/2]; 
				printf(" %02lX %01ld %2d %01.3lf ",i+dictionary.size()/2,i%2,w.state,w.p);
				for (size_t j=0; j<16; j++) putchar("0123456789ABCDEF "[j<w.symbols.size()?w.symbols[j]:16]);
			}
			putchar('\n');
		}		
		putchar('\n');
	}

	Marlin2Dictionary(const std::vector<double> &P, size_t dictSize, size_t tries=4) {
		
		std::vector<double> PstatesA(P.size(), 0.); PstatesA[0] = .5;
		std::vector<double> PstatesB(P.size(), 0.); PstatesB[0] = .5;
		
		dictionary = concatenate( 
			arrangeWords( buildWords( buildTree(P, PstatesA, dictSize/2) ) ),
			arrangeWords( buildWords( buildTree(P, PstatesB, dictSize/2) ) ) );
			
		print(dictionary);
			
		while (tries--) {

			// UPDATING STATE PROBABILITIES
			{
				PstatesA = std::vector<double>(P.size(), 0.);
				PstatesB = std::vector<double>(P.size(), 0.);

				for (size_t i=0; i<dictionary.size(); i++) {
					
					//Last bit, determines next dictionary
					if (i%2)
						PstatesB[dictionary[i].state] += dictionary[i].p;
					else
						PstatesA[dictionary[i].state] += dictionary[i].p;
				}
			}
			
			if (P.size()<16)
				for (size_t i=0; i<P.size(); i++)
					printf("S: %02ld %01.3lf %01.3lf\n",i,PstatesA[i],PstatesB[i]);
			putchar('\n');

			
			dictionary = concatenate( 
				arrangeWords( buildWords( buildTree(P, PstatesA, dictSize/2) ) ),
				arrangeWords( buildWords( buildTree(P, PstatesB, dictSize/2) ) ) );

			print(dictionary);
		}
		
		//test(P,dictSize);
		averageBitsPerSymbolEmpirically(P,dictSize);
	}
	
	double averageBitsPerSymbolEmpirically(const std::vector<double> &P, size_t dictSize) const {
		
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
		
		// pad with 0s
		for (size_t i=0; i<10; i++)
			testData.push_back(0);

		auto W = dictionary;
		double meanLengthE = 0;
		size_t nWords = 0;
		size_t lastWord = 0;
		for (size_t i=0, longest=0; i<TEST_SIZE; i+=longest) {
			
			size_t best = 0;
			longest = 0;
			for (size_t j=0; j<W.size()/2; j++) {
				auto &&w = W[(lastWord%2)*(W.size()/2) + j].symbols;
				if (w.size()>longest and testData.compare(i,w.size(),w)==0 ) {
					best = (lastWord%2)*(W.size()/2) + j;
					longest = w.size();
				}
			}
			lastWord = best;
			meanLengthE += longest;
			nWords++;
		}

		printf("Empirical Meanlength: %3.2lf\n", meanLengthE/nWords);
		printf("Empirical Compress Ratio: %3.4lf\n", 
			(nWords*(std::log2(dictSize)-1))/(double(TEST_SIZE)*std::log2(P.size())));


		return double(nWords*(std::log2(dictSize)-1))/TEST_SIZE;
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
	MarlinDictionary(P,options["--size"],options["--tries"]);
	std::cerr << "Tunstall" << std::endl;
	//TunstallDictionary(P,options["--size"]);
	std::cerr << "Marlin2" << std::endl;
	Marlin2Dictionary(P,options["--size"]);
		
	return 0;
}
