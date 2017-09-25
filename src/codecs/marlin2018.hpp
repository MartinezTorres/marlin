#pragma once
#include <util/codec.hpp>
#include <util/distribution.hpp>

struct Marlin2018 : public CODEC8withPimpl { 

	Marlin2018(
		Distribution::Type distType = Distribution::Laplace, 
		size_t keySize = 12, 
		size_t overlap = 2,
		size_t numDict = 11);	
};
