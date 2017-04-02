#pragma once
#include <util/codec.hpp>

struct Gzip : public CODEC8withPimpl { Gzip(int level=1); };
