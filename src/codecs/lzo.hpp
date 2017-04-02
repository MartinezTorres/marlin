#pragma once
#include <util/codec.hpp>

struct Lzo : public CODEC8withPimpl { Lzo(int level=1); };
