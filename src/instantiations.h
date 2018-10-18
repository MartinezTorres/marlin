#define INSTANTIATE(A) \
	template class marlin::A<uint16_t,uint8_t>; 
//	template class marlin::A<uint16_t,uint8_t>; 

#define INSTANTIATE_MEMBER(A,B) \
	template auto marlin::A<uint8_t,uint8_t>::B;
//	template auto marlin::A<uint16_t,uint8_t>::B;


