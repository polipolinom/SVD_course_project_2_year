#include <utility>

namespace svd_computation {
std::pair<long double, long double> get_givens_rotation(long double a, long double b,
                                                        long double eps = constants::DEFAULT_EPSILON) {
    if (abs(b) < eps) {
        return {1.0, 0.0};
    }
    if (abs(b) > abs(a)) {
        long double t = -a / b;
        return {t / sqrtl(1 + t * t), 1 / sqrtl(1 + t * t)};
    }
    long double t = -b / a;
    return {1 / sqrtl(1 + t * t), t / sqrtl(1 + t * t)};
}
}  // namespace svd_computation