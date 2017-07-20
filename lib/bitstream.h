namespace {
	
struct obitstream {
	
	uint8_t *start;
	uint8_t *p;
	uint8_t *end;
	uint64_t val = 0;
	int16_t roll = 0;
	
	constexpr 
	void emitByte() {
		
	//	std::cout << "B " << roll << std::endl;

		*p++ = (val & 0xFF);
		val = (val >> 8);
		roll -= 8;
	}
	
	constexpr obitstream(uint8_t *p, size_t sz) : start(p), p(p), end(p+sz) {}
	
	template<typename T>		
	constexpr 
	void write(size_t bits, const T& data) {
		
		for (size_t b = 0, c = 0; b < bits; b+=8, ++c) {
			//std::cout << "B " << b << " " << bits << " " << data << " " << uint32_t(((const uint8_t *)&data)[c]) << " " << int(bits-b) << std::endl;
			write8(std::min(8,int(bits-b)), ((const uint8_t *)&data)[c]);
		}
	}
	
	constexpr 
	void write8(size_t bits, uint8_t data) {

		
		data = (data & ((1<<bits)-1));
//			std::cout << "K1 " << data << " " << val << std::endl;
//			std::cout << "e " << bits <<  " " << roll << std::endl;
		val += (uint64_t(data) << roll);
//			std::cout << "K2 " << val << std::endl;
		roll += bits;

		while (roll >= 8 and p!= end) 
			emitByte();
	}
	
	constexpr void sync() {
		
		while (roll+7 >= 8 and p!= end) 
			emitByte();
	}
	
	constexpr size_t size() { return p-start; }
};

struct ibitstream {
	
	const uint8_t *start;
	const uint8_t *p;
	const uint8_t *end;
	uint64_t val = 0;
	int64_t roll = 0;
	
	constexpr ibitstream(const uint8_t *p, size_t sz) : start(p), p(p), end(p+sz) {}
	
	constexpr uint64_t read(size_t bits) {

		if (roll < int(bits)) {
			if (p + 4 < end) {
				
				val += ( uint64_t(*(uint32_t *)p) << roll);
				p += sizeof(uint32_t);
				roll += sizeof(uint32_t)*8;
			} else while (p < end) {
				
				val += ( uint64_t(*(uint8_t *)p) << roll);
				p += sizeof(uint8_t);
				roll += sizeof(uint8_t)*8;
			}
		}

		uint64_t ret = val & uint64_t((1<<bits)-1);
		roll -= bits;
		val = (val >> bits);
		return ret;
	}
	
	constexpr uint64_t readBytes(size_t bytes) {
		
		uint64_t ret = *(uint64_t *)p;
		p += bytes;
		return ret;
	}
	
	constexpr operator bool() const {
		return (roll>0) or (p!=end);
	}
	
	constexpr void sync() {}
	
	constexpr size_t bytesLeft() const {
		return end - p;
	}
};
}
