#pragma once
#include <util/codec.hpp>
#include <util/distribution.hpp>

struct Marlin2 : public CODEC8withPimpl { 

	enum Type { TUNSTALL, MARLIN };	
	Marlin2(Distribution::Type distType = Distribution::Laplace, Type dictType = MARLIN, size_t dictSize = 12, size_t numDict = 11);	
};
