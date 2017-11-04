#pragma once
#include <util/codec.hpp>
#include <util/distribution.hpp>

struct MarlinCodecWiP : public CODEC8withPimpl { 

	MarlinCodecWiP(
		Distribution::Type distType = Distribution::Laplace, 
		size_t overlap = 2,
		size_t numDict = 11);	
};
