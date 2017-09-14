#include <queue>
#include <iostream>
#include <algorithm>
#include <memory>
#include <stack>
#include <list>
#include <map>

#include <distribution.hpp>

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

	Dictionary buildSingleDictionary(const std::vector<double> &P, const std::vector<double> &Pstates, size_t dictSize, bool allSymbols) {

		Dictionary ret(P,dictSize);
		
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
			ret.root = std::make_shared<Node>();

			for (size_t c=0; c<P.size(); c++) {
				
				ret.root->push_back(std::make_shared<Node>());

				double sum = 0;
				for (size_t t = 0; t<=c; t++) sum += Pstates[t]/PN[t];

				ret.root->back()->p = sum * P[c];
				
				pq.push(ret.root->back());
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

	Marlin2Dictionary(const std::vector<double> &P, size_t dictSize, size_t tries=3) {
		
		std::vector<double> PstatesInit(P.size(), 0.); PstatesInit[0] = 1.;
		
		Dictionary dictionaryA = buildSingleDictionary(P, PstatesInit, dictSize/2, true);
		Dictionary dictionaryB = dictionaryA;
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
	
	MarlinDictionary(P,options["--size"],options["--tries"]);
	TunstallDictionary(P,options["--size"]);
		
	return 0;
}
