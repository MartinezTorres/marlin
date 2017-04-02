#include <codecs/lzo.hpp>
#include <lzo/lzo1x.h>

class LzoPimpl : public CODEC8AA {
	
	int level;
	
	std::string name() const { return std::string("Lzo1-15 ")+char('0'+level); };

	
	void   compress(const AlignedArray8 &in, AlignedArray8 &out) const {

		uint8_t workMem[LZO1X_1_15_MEM_COMPRESS];

		lzo_uint sz = out.capacity();
		lzo1x_1_15_compress( in.begin(), in.size(), out.begin(), &sz, workMem );
		out.resize(sz);
	}

	void uncompress(const AlignedArray8 &in, AlignedArray8 &out) const {
		
		lzo_uint sz = out.capacity();
		lzo1x_decompress( in.begin(), in.size(), out.begin(), &sz, nullptr);
		out.resize(sz);
	}

public:
	
	LzoPimpl(int level) : level(level) {}

};

Lzo::Lzo(int level) : CODEC8withPimpl(new LzoPimpl(level)) {}
