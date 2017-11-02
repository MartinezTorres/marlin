#include <codecs/gzip.hpp>
#include <zlib.h>


class GzipPimpl : public CODEC8AA {
	
	int level;
	
	std::string name() const { return std::string("Gzip")+char('0'+level); };

	void   compress(const AlignedArray8 &in, AlignedArray8 &out) const {

		z_stream strm;
		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;
		deflateInit(&strm, level);
		
		strm.avail_in = in.size();
		strm.next_in = const_cast<Bytef*>(&in[0]);
		strm.avail_out = out.capacity();
		strm.next_out = (Bytef*)out.begin();

		deflate( &strm, Z_FINISH);
		assert( strm.avail_out != 0 );

		out.resize( out.capacity()-strm.avail_out );
		deflateEnd(&strm);		
	}

	void uncompress(const AlignedArray8 &in, AlignedArray8 &out) const {
		
		z_stream strm;
		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;
		inflateInit(&strm);
		
		strm.avail_in = in.size();
		strm.next_in = const_cast<Bytef*>(&in[0]);
		strm.avail_out = out.capacity();
		strm.next_out = (Bytef*)out.begin();

		inflate( &strm, Z_NO_FLUSH);
		assert( strm.avail_out != 0 );
		out.resize( out.capacity()-strm.avail_out );
		
		inflateEnd(&strm);		
	}

public:
	
	GzipPimpl(int level_) : level(level_) {}

};

Gzip::Gzip(int level) : CODEC8withPimpl(new GzipPimpl(level)) {}
