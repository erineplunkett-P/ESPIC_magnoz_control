#include "comp_ellint.h"
#include <limits>
#include <stdexcept>
#include <algorithm>

namespace ellint
{

    // ------------------------------
    // Internal: Carlson RF and RD
    // ------------------------------

    template <class T>
    static T carlson_RF(T x, T y, T z)
    {
        // SLATEC / Carlson duplication with 5th-order correction (DRF)
        // Coefficients from SLATEC RF: C1=1/24, C2=3/44, C3=1/14
        using std::abs;
        using std::sqrt;

        if (x < T(0) || y < T(0) || z < T(0))
            throw std::domain_error("RF: arguments must be >= 0");

        const T ERRTOL = std::pow(T(4) * std::numeric_limits<T>::epsilon(), T(1) / T(6));

        T XN = x, YN = y, ZN = z;
        for (;;)
        {
            const T MU = (XN + YN + ZN) / T(3);
            const T XNDEV = T(2) - (MU + XN) / MU;
            const T YNDEV = T(2) - (MU + YN) / MU;
            const T ZNDEV = T(2) - (MU + ZN) / MU;

            const T EPS = std::max(abs(XNDEV), std::max(abs(YNDEV), abs(ZNDEV)));
            if (EPS < ERRTOL)
            {
                const T E2 = XNDEV * YNDEV - ZNDEV * ZNDEV;
                const T E3 = XNDEV * YNDEV * ZNDEV;
                const T C1 = T(1) / T(24);
                const T C2 = T(3) / T(44);
                const T C3 = T(1) / T(14);
                const T S = T(1) + (C1 * E2 - T(0.1) - C2 * E3) * E2 + C3 * E3;
                return S / sqrt(MU);
            }

            const T XNROOT = sqrt(XN);
            const T YNROOT = sqrt(YN);
            const T ZNROOT = sqrt(ZN);
            const T LAM = XNROOT * (YNROOT + ZNROOT) + YNROOT * ZNROOT;

            XN = (XN + LAM) / T(4);
            YN = (YN + LAM) / T(4);
            ZN = (ZN + LAM) / T(4);
        }
    }

    template <class T>
    static T carlson_RD(T x, T y, T z)
    {
        // SLATEC / Carlson duplication for RD with 5th-order correction (DRD)
        // Coefficients: C1=3/14, C2=1/6, C3=9/22, C4=3/26
        using std::abs;
        using std::sqrt;

        if (x < T(0) || y < T(0) || z <= T(0))
            throw std::domain_error("RD: need x,y >= 0 and z > 0");

        const T ERRTOL = std::pow(std::numeric_limits<T>::epsilon() / T(3), T(1) / T(6));

        T XN = x, YN = y, ZN = z;
        T SIGMA = T(0);
        T POWER4 = T(1);

        for (;;)
        {
            const T MU = (XN + YN + T(3) * ZN) / T(5);
            const T XNDEV = (MU - XN) / MU;
            const T YNDEV = (MU - YN) / MU;
            const T ZNDEV = (MU - ZN) / MU;

            const T EPS = std::max(abs(XNDEV), std::max(abs(YNDEV), abs(ZNDEV)));
            if (EPS < ERRTOL)
            {
                const T EA = XNDEV * YNDEV;
                const T EB = ZNDEV * ZNDEV;
                const T EC = EA - EB;
                const T ED = EA - T(6) * EB;
                const T EF = ED + EC + EC;

                const T C1 = T(3) / T(14);
                const T C2 = T(1) / T(6);
                const T C3 = T(9) / T(22);
                const T C4 = T(3) / T(26);

                const T S1 = ED * (-C1 + T(0.25) * C3 * ED - T(1.5) * C4 * ZNDEV * EF);
                const T S2 = ZNDEV * (C2 * EF + ZNDEV * (-C3 * EC + ZNDEV * C4 * EA));

                return T(3) * SIGMA + POWER4 * (T(1) + S1 + S2) / (MU * sqrt(MU));
            }

            const T XNROOT = sqrt(XN);
            const T YNROOT = sqrt(YN);
            const T ZNROOT = sqrt(ZN);
            const T LAM = XNROOT * (YNROOT + ZNROOT) + YNROOT * ZNROOT;

            SIGMA += POWER4 / (ZNROOT * (ZN + LAM));
            POWER4 *= T(0.25);

            XN = (XN + LAM) / T(4);
            YN = (YN + LAM) / T(4);
            ZN = (ZN + LAM) / T(4);
        }
    }

    // ------------------------------
    // Public: comp_ellint_1 and _2
    // ------------------------------

    template <class T>
    static T comp_ellint_1_impl(T k)
    {
        using std::abs;
        using std::numeric_limits;

        if (std::isnan(k))
            return numeric_limits<T>::quiet_NaN();
        const T kk = abs(k);
        if (kk > T(1))
            throw std::domain_error("comp_ellint_1: |k| > 1 gives complex result");
        if (kk == T(1))
            return numeric_limits<T>::infinity();

        const T m = kk * kk;
        return carlson_RF(T(0), T(1) - m, T(1));
    }

    template <class T>
    static T comp_ellint_2_impl(T k)
    {
        using std::abs;
        using std::numeric_limits;

        if (std::isnan(k))
            return numeric_limits<T>::quiet_NaN();
        const T kk = abs(k);
        if (kk > T(1))
            throw std::domain_error("comp_ellint_2: |k| > 1 gives complex result");
        if (kk == T(1))
            return T(1); // E(1) = 1

        const T m = kk * kk;
        const T rf = carlson_RF(T(0), T(1) - m, T(1));
        const T rd = carlson_RD(T(0), T(1) - m, T(1));
        return rf - (m / T(3)) * rd;
    }

    // float/double/long double overloads

    float comp_ellint_1(float k) { return comp_ellint_1_impl(k); }
    double comp_ellint_1(double k) { return comp_ellint_1_impl(k); }
    long double comp_ellint_1(long double k) { return comp_ellint_1_impl(k); }

    float comp_ellint_2(float k) { return comp_ellint_2_impl(k); }
    double comp_ellint_2(double k) { return comp_ellint_2_impl(k); }
    long double comp_ellint_2(long double k) { return comp_ellint_2_impl(k); }

} // namespace ellint

// namespace ellint
// {

//     template <class T>
//     struct KE
//     {
//         T K, E;
//     };

//     template <class T>
//     static KE<T> agm_complete_KE(T k)
//     {
//         using std::abs;
//         using std::sqrt;

//         const T NaN = std::numeric_limits<T>::quiet_NaN();
//         if (std::isnan(k))
//             return {NaN, NaN};

//         const T kk = abs(k);
//         if (kk > T(1))
//         {
//             throw std::domain_error("comp_ellint_*: |k| > 1 yields complex result");
//         }

//         // Special cases
//         if (kk == T(1))
//         {
//             return {std::numeric_limits<T>::infinity(), T(1)}; // K(1)=∞, E(1)=1
//         }

//         // AGM iteration
//         T a = T(1);
//         T b = sqrt(T(1) - kk * kk);

//         T sum = T(0); // Σ 2^n c_n^2
//         T pow2 = T(1);

//         const T tol = std::numeric_limits<T>::epsilon() * T(8);

//         for (int iter = 0; iter < 64; ++iter)
//         {
//             // c_n from current (a_n, b_n)
//             const T c = (a - b) / T(2);
//             sum += pow2 * c * c;
//             pow2 *= T(2);

//             // next (a_{n+1}, b_{n+1})
//             const T an = (a + b) / T(2);
//             const T bn = sqrt(a * b);
//             a = an;
//             b = bn;

//             if (abs(a - b) <= tol * a)
//                 break;
//         }

//         const T pi = T(3.141592653589793238462643383279502884L);
//         const T K = pi / (T(2) * a);
//         const T E = K * (a * a - sum);
//         return {K, E};
//     }

//     template <class T>
//     static T comp_ellint_1_impl(T k) { return agm_complete_KE(k).K; }

//     template <class T>
//     static T comp_ellint_2_impl(T k) { return agm_complete_KE(k).E; }

//     float comp_ellint_1(float k) { return comp_ellint_1_impl(k); }
//     double comp_ellint_1(double k) { return comp_ellint_1_impl(k); }
//     long double comp_ellint_1(long double k) { return comp_ellint_1_impl(k); }

//     float comp_ellint_2(float k) { return comp_ellint_2_impl(k); }
//     double comp_ellint_2(double k) { return comp_ellint_2_impl(k); }
//     long double comp_ellint_2(long double k) { return comp_ellint_2_impl(k); }

// } // namespace ellint