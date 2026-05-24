// Smoke test for fp80 on aarch64 (or any non-x86 platform).
//
// Build:  aarch64-linux-gnu-g++ -std=c++20 -O2 -static -I.. \
//             -o x87aarch64_smoke x87aarch64_smoke.cpp
// Run:    qemu-aarch64 ./x87aarch64_smoke
//
// Computes a handful of fp80 operations against precomputed golden values
// (captured on a real Intel x87 via the main test suite). Each line prints
// PASS / FAIL plus the result bits. A non-zero exit code if any FAILs.

#include <cstdio>
#include <cstdint>
#include <cstring>

#include "../x87fp80.h"
#include "../x87fp64.h"

// Library implementation files included as a unit (matches the main test
// program's organization).
#include "../x87fp64trans.cpp"
#include "../x87fp80.cpp"
#include "../x87fp80trans.cpp"

using namespace x87;

static int failures = 0;

static void check(const char *label, fp80_t actual, uint16_t expect_se, uint64_t expect_mant)
{
    bool ok = (actual.sign_exp() == expect_se) && (actual.mantissa() == expect_mant);
    std::printf("%-30s %s  %04X:%016lX %s %04X:%016lX\n",
                label, ok ? "PASS" : "FAIL",
                actual.sign_exp(), actual.mantissa(),
                ok ? "==" : "!=",
                expect_se, expect_mant);
    if (!ok) failures++;
}

int main()
{
    // Constants
    check("const_one",  fp80_t::const_one(),  0x3FFF, 0x8000000000000000ull);
    check("const_pi",   fp80_t::const_pi(),   0x4000, 0xc90fdaa22168c235ull);
    check("const_l2t",  fp80_t::const_l2t(),  0x4000, 0xd49a784bcd1b8afeull);
    check("const_ln2",  fp80_t::const_ln2(),  0x3FFE, 0xb17217f7d1cf79acull);

    // Predicates
    auto qnan_v = fp80_t::const_qnan();
    auto nan_v  = fp80_t::const_snan();
    auto inf_v  = fp80_t::const_pinf();
    auto zero_v = fp80_t::const_zero();
    std::printf("%-30s %s\n", "isnan(qnan)",
                qnan_v.isnan() ? "PASS" : "FAIL"); failures += !qnan_v.isnan();
    std::printf("%-30s %s\n", "issnan(snan)",
                nan_v.issnan() ? "PASS" : "FAIL"); failures += !nan_v.issnan();
    std::printf("%-30s %s\n", "isinf(pinf)",
                inf_v.isinf() ? "PASS" : "FAIL"); failures += !inf_v.isinf();
    std::printf("%-30s %s\n", "iszero(zero)",
                zero_v.iszero() ? "PASS" : "FAIL"); failures += !zero_v.iszero();

    // Arithmetic: 1.0 + 1.0 = 2.0
    {
        fp80_t a = fp80_t::const_one();
        fp80_t b = fp80_t::const_one();
        fp80_t r;
        fp80_t::x87_fadd(a, b, r);
        check("1.0 + 1.0", r, 0x4000, 0x8000000000000000ull);
    }
    // 2.0 * 3.0 = 6.0  (2.0 = 0x4000:0x8000.., 3.0 = 0x4000:0xC000..)
    {
        fp80_t a(0x8000000000000000ull, 0x4000); // 2.0
        fp80_t b(0xC000000000000000ull, 0x4000); // 3.0
        fp80_t r;
        fp80_t::x87_fmul(a, b, r);
        check("2.0 * 3.0 = 6.0", r, 0x4001, 0xC000000000000000ull); // 6.0
    }
    // 1.0 / 2.0 = 0.5
    {
        fp80_t a = fp80_t::const_one();                                  // 1.0
        fp80_t b(0x8000000000000000ull, 0x4000);                          // 2.0
        // Asm convention: x87_fdiv(a,b) computes b/a; for 1/2 we want a/b
        fp80_t r;
        fp80_t::x87_fdivr(a, b, r);  // = 1.0 / 2.0
        check("1.0 / 2.0 = 0.5", r, 0x3FFE, 0x8000000000000000ull);     // 0.5
    }
    // sqrt(4.0) = 2.0
    {
        fp80_t four(0x8000000000000000ull, 0x4001);                       // 4.0
        fp80_t r;
        fp80_t::x87_fsqrt(four, r);
        check("sqrt(4.0) = 2.0", r, 0x4000, 0x8000000000000000ull);
    }
    // sqrt(2.0) ~= 1.4142135623730950488
    // On aarch64 the lib uses HW fp64 fsqrt (X87_MATCH_XEFU-class
    // precision) so the lower 11 mantissa bits are zero. The C
    // shift-and-subtract path produces the bit-exact fp80 value.
    {
        fp80_t two(0x8000000000000000ull, 0x4000);                        // 2.0
        fp80_t r;
        fp80_t::x87_fsqrt(two, r);
#if defined(__aarch64__)
        check("sqrt(2.0)", r, 0x3FFF, 0xB504F333F9DE6800ull);              // fp64-precision
#else
        check("sqrt(2.0)", r, 0x3FFF, 0xB504F333F9DE6484ull);              // fp80-precision
#endif
    }

    // Comparisons
    {
        fp80_t one  = fp80_t::const_one();
        fp80_t two(0x8000000000000000ull, 0x4000);
        bool ok = (one < two) && !(two < one) && !(one == two) && (one != two);
        std::printf("%-30s %s\n", "comparisons", ok ? "PASS" : "FAIL");
        if (!ok) failures++;
    }

    std::printf("\n%s: %d failure(s)\n",
                failures == 0 ? "OVERALL PASS" : "OVERALL FAIL",
                failures);
    return failures != 0;
}
