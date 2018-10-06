#define INSTANTIATE() \
	template class marlin::TMarlin<uint8_t,uint8_t>; \
	template class marlin::TMarlin<uint16_t,uint8_t>; \
	template class marlin::TMarlin<uint8_t,uint16_t>; 

#define INSTANTIATE_MEMBER(A) \
	template auto marlin::TMarlin<uint8_t,uint8_t>::A; \
	template auto marlin::TMarlin<uint16_t,uint8_t>::A; \
	template auto marlin::TMarlin<uint8_t,uint16_t>::A;
