#include <cmath>
#include <limits>
#include <stdexcept>

#ifndef ELLINT_H
#define ELLINT_H

#include <cmath>

namespace ellint
{
    // Complete elliptic integral of the first kind K(k) for real |k| <= 1.
    float comp_ellint_1(float k);
    double comp_ellint_1(double k);
    long double comp_ellint_1(long double k);

    // Complete elliptic integral of the second kind E(k) for real |k| <= 1.
    float comp_ellint_2(float k);
    double comp_ellint_2(double k);
    long double comp_ellint_2(long double k);

} // namespace ellint

namespace ellint_F
{
    // Complete elliptic integral of the first kind K(k) for real |k| <= 1.
    float comp_ellint_1(float k);
    double comp_ellint_1(double k);
    long double comp_ellint_1(long double k);

    // Complete elliptic integral of the second kind E(k) for real |k| <= 1.
    float comp_ellint_2(float k);
    double comp_ellint_2(double k);
    long double comp_ellint_2(long double k);

} // namespace ellint

#endif // ELLINT_H