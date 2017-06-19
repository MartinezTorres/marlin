#include <codecs/snappy.hpp>
#include <snappy.h>


class SnappyPimpl : public CODEC8AA {
	
	std::string name() const { return "Snappy"; }

	void   compress(const AlignedArray8 &in, AlignedArray8 &out) const {

		size_t compressed_length;
		snappy::RawCompress((const char *)in.begin(), in.size(), (char *)out.begin(), &compressed_length);
		out.resize(compressed_length);
	}

	void uncompress(const AlignedArray8 &in, AlignedArray8 &out) const {

		size_t uncompressed_length;
		snappy::GetUncompressedLength((const char *)in.begin(), in.size(), &uncompressed_length);
		out.resize(uncompressed_length);
		snappy::RawUncompress((const char *)in.begin(), in.size(), (char *)out.begin());
	}
};

Snappy::Snappy() : CODEC8withPimpl(new SnappyPimpl()) {}
