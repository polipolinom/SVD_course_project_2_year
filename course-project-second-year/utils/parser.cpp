#include "parser.h"

#include <cstdlib>
#include <iomanip>

namespace svd_computation {
namespace details {
Complex Parser::parse(std::string buf) {
    std::string::size_type sz;

    long double num1, num2;

    try {
        num1 = std::stold(buf, &sz);
    } catch (std::invalid_argument& e) {
        throw std::invalid_argument("Incorrect format of complex number " + buf);
    } catch (std::out_of_range& e) {
        throw std::out_of_range("Out of range in complex number " + buf);
    }
    if (sz == buf.size()) {
        Complex result = num1;
        return result;
    }

    if (buf[sz] == '*' || buf[sz] == 'i') {
        Complex result = Complex(0.0l, num1);
        return result;
    }

    try {
        num2 = std::stold(buf.substr(sz), &sz);
    } catch (std::invalid_argument& e) {
        throw std::invalid_argument("Incorrect format of complex number " + buf);
    } catch (std::out_of_range& e) {
        throw std::out_of_range("Out of range in complex number " + buf);
    }

    Complex result = Complex(num1, num2);

    return result;
}
}  // namespace details
}  // namespace svd_computation