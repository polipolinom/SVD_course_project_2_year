#pragma once

#include <string>

#include "complex.h"

namespace svd_computation {
namespace details {
class Parser {
   public:
    Complex parse(::std::string);
};
}  // namespace details
}  // namespace svd_computation