#pragma once
#include <util/codec.hpp>
#include <util/distribution.hpp>

struct Rice : public CODEC8withPimpl { Rice(Distribution::Type distType = Distribution::Laplace); };

