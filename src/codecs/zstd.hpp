#pragma once
#include <util/codec.hpp>

struct Zstd : public CODEC8withPimpl { Zstd(int level = 1); };
