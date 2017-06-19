#include <codecs/huf.hpp>
#include <FiniteStateEntropy/lib/huf.h>

class Huff0Pimpl : public CODEC8AA {
	
	std::string name() const { return "Huff0"; };

	void   compress(const AlignedArray8 &in, AlignedArray8 &out) const {
		
		size_t ret =  HUF_compress(out.begin(), out.capacity(), in.begin(), in.size());
		if (ret==0) {
			out=in;
			return;
		} else if (ret==1) { 
			out[0] = in[0];
			out.resize(1);
		}
		out.resize(ret);
	}

	void uncompress(const AlignedArray8 &in, AlignedArray8 &out) const {

		if (in.size()==out.size()) {
			out = in;
		} else if (in.size()==1) {
			memset(out.begin(), in[0], in.size());
		} else {
			HUF_decompress(out.begin(), out.size(), in.begin(), in.size());
		}
	}
};

Huff0::Huff0() : CODEC8withPimpl(new Huff0Pimpl()) {}

