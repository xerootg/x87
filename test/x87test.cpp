// Build (Windows):  cl /EHsc /Zi /std:c++20 x87test.cpp x87testasm.obj
// Build (Linux):    see Makefile (g++ -std=c++20 x87test.cpp x87testasm.o)

#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <format>
#include <assert.h>
#include <chrono>

//
// references:
//    https://www.researchgate.net/publication/3612479_The_K5_transcendental_functions
//    https://www.researchgate.net/publication/3798381_New_algorithms_for_improved_transcendental_functions_on_IA-64
//

#include "../x87fp64trans.cpp"
#include "../x87fp80.cpp"
#include "../x87fp80trans.cpp"
#include "../x87fpext.h"
#undef print_val

using namespace x87;

//
// external assembly stubs with actual implementations to compare against
//
extern "C"
{
    uint16_t x87getsw();
    void x87consts80(fp80_t *dst);
    void x87consts64(fp64_t *dst);
    void x87setcw(uint16_t *val);
    void x87test1(fp64_t *src, fp80_t *dst);
    void x87test2(fp80_t *src, fp64_t *dst);

    uint16_t fld8080(fp80_t const *src, fp80_t *dst);
    uint16_t fld6480(fp64_t const *src, fp80_t *dst);
    uint16_t fld3280(float const *src, fp80_t *dst);
    uint16_t fild6480(int64_t const *src, fp80_t *dst);
    uint16_t fild3280(int32_t const *src, fp80_t *dst);
    uint16_t fild1680(int16_t const *src, fp80_t *dst);
    uint16_t fst8080(fp80_t const *src, fp80_t *dst);
    uint16_t fst8064(fp80_t const *src, fp64_t *dst);
    uint16_t fst8032(fp80_t const *src, float *dst);
    uint16_t fist8064(fp80_t const *src, int64_t *dst);
    uint16_t fist8032(fp80_t const *src, int32_t *dst);
    uint16_t fist8016(fp80_t const *src, int16_t *dst);
    uint16_t fadd80(fp80_t *src1, fp80_t *src2, fp80_t *dst);
    uint16_t fsub80(fp80_t *src1, fp80_t *src2, fp80_t *dst);
    uint16_t fmul80(fp80_t *src1, fp80_t *src2, fp80_t *dst);
    uint16_t fdiv80(fp80_t *src1, fp80_t *src2, fp80_t *dst);
    uint16_t fsqrt80(fp80_t *src, fp80_t *dst);
    uint16_t f2xm180(fp80_t *src, fp80_t *dst);
    uint16_t fyl2x80(fp80_t *src1, fp80_t *src2, fp80_t *dst);
    uint16_t fptan80(fp80_t *src, fp80_t *dst1, fp80_t *dst2);
    uint16_t fsincos80(fp80_t *src, fp80_t *dst1, fp80_t *dst2);
    uint16_t fpatan80(fp80_t *src1, fp80_t *src2, fp80_t *dst);
    uint16_t fxtract80(fp80_t *src, fp80_t *dst1, fp80_t *dst2);
    uint16_t fprem180(fp80_t *src1, fp80_t *src2, fp80_t *dst);
    uint16_t fprem80(fp80_t *src1, fp80_t *src2, fp80_t *dst);
    uint16_t fyl2xp180(fp80_t *src1, fp80_t *src2, fp80_t *dst);
    uint16_t frndint80(fp80_t *src, fp80_t *dst);
    uint16_t fscale80(fp80_t *src1, fp80_t *src2, fp80_t *dst);
    uint16_t fsin80(fp80_t *src, fp80_t *dst);
    uint16_t fcos80(fp80_t *src, fp80_t *dst);
    uint16_t fcmp80(fp80_t const *src1, fp80_t const *src2);
    uint16_t fxam80(fp80_t const *src);
    uint16_t fabs80(fp80_t const *src, fp80_t *dst);
    uint16_t fchs80(fp80_t const *src, fp80_t *dst);
    uint16_t fsubr80(fp80_t const *src1, fp80_t const *src2, fp80_t *dst);
    uint16_t fdivr80(fp80_t const *src1, fp80_t const *src2, fp80_t *dst);
    uint16_t ftst80(fp80_t const *src);
    uint16_t fcom80(fp80_t const *src1, fp80_t const *src2);
    uint16_t fcomi80(fp80_t const *src1, fp80_t const *src2);
    uint16_t fucomi80(fp80_t const *src1, fp80_t const *src2);

    uint16_t fadd64(fp64_t *src1, fp64_t *src2, fp64_t *dst);
    uint16_t fsub64(fp64_t *src1, fp64_t *src2, fp64_t *dst);
    uint16_t fmul64(fp64_t *src1, fp64_t *src2, fp64_t *dst);
    uint16_t fdiv64(fp64_t *src1, fp64_t *src2, fp64_t *dst);
    uint16_t fsqrt64(fp64_t *src, fp64_t *dst);
    uint16_t f2xm164(fp64_t *src, fp64_t *dst);
    uint16_t fyl2x64(fp64_t *src1, fp64_t *src2, fp64_t *dst);
    uint16_t fptan64(fp64_t *src, fp64_t *dst1, fp64_t *dst2);
    uint16_t fsincos64(fp64_t *src, fp64_t *dst1, fp64_t *dst2);
    uint16_t fpatan64(fp64_t *src1, fp64_t *src2, fp64_t *dst);
    uint16_t fxtract64(fp64_t *src, fp64_t *dst1, fp64_t *dst2);
    uint16_t fprem164(fp64_t *src1, fp64_t *src2, fp64_t *dst);
    uint16_t fprem64(fp64_t *src1, fp64_t *src2, fp64_t *dst);
    uint16_t fyl2xp164(fp64_t *src1, fp64_t *src2, fp64_t *dst);
    uint16_t frndint64(fp64_t *src, fp64_t *dst);
    uint16_t fscale64(fp64_t *src1, fp64_t *src2, fp64_t *dst);
    uint16_t fsin64(fp64_t *src, fp64_t *dst);
    uint16_t fcos64(fp64_t *src, fp64_t *dst);
}

//
// global variables
//
std::chrono::nanoseconds min_timing_duration;
std::vector<fp80_t> values80;
std::vector<fp64_t> values64;
std::vector<float> values32;
std::vector<int64_t> valuesi64;
std::vector<int32_t> valuesi32;
std::vector<int16_t> valuesi16;

#define MAX_PRINT_ERRORS 1000

//
// print helper since we're C++20
//
template<typename... Args>
void print(std::format_string<Args...> fmt, Args &&...args)
{
    puts(std::format(fmt, std::forward<Args>(args)...).c_str());
}

//
// stderr print helper since we're C++20
//
template<typename... Args>
void eprint(std::format_string<Args...> fmt, Args &&...args)
{
    fputs(std::format(fmt, std::forward<Args>(args)...).c_str(), stderr);
}

//
// formatter for fp64_t
//
template <> class std::formatter<fp64_t>
{
public:
    constexpr auto parse (format_parse_context& ctx) { return ctx.begin(); }
    template <typename Context>
    constexpr auto format (fp64_t const &src, Context &ctx) const
    {
        return format_to(ctx.out(), "{:016X} [{:+.12e}]", src.as_fpbits64(), src.as_double());
    }
};

//
// formatter for fp80_t
//
template <> class std::formatter<fp80_t>
{
public:
    constexpr auto parse (format_parse_context& ctx) { return ctx.begin(); }
    template <typename Context>
    constexpr auto format (fp80_t const &src, Context &ctx) const
    {
        return format_to(ctx.out(), "{:04X}:{:016X} [{:+.12e}]", src.sign_exp(), src.mantissa(), src.as_double());
    }
};

//
// helper class for reporting and tracking errors
//
struct errors_t
{
    //
    // core error counts
    //
    uint32_t count = 0;
    uint32_t printed = 0;
    uint32_t matches = 0;
    uint32_t signerrors = 0;
    uint32_t experrors = 0;
    uint32_t infinities = 0;
    uint32_t swerrors = 0;
    std::array<uint32_t, 64> errors = {{ 0 }};
    char const *name = nullptr;
    int print_thresh;

    //
    // constructor
    //
    errors_t(char const *_name, int _thresh) :
        name(_name), print_thresh(_thresh)
    {
    }

    //
    // check an integer value against a result from the real FPU
    //
    template<typename PrintFuncType, typename ReExecFuncType, typename IntType>
    void check_value(IntType const &ourdst, IntType const &x87dst, uint16_t oursw, uint16_t x87sw, PrintFuncType printname, ReExecFuncType reex)
    {
        bool print = false;
        count++;

        oursw &= ~X87SW_TOP_MASK;
        x87sw &= ~X87SW_TOP_MASK;

        if (x87dst != ourdst)
            experrors++, print = true;
        else if (x87sw != oursw)
            swerrors++, print = true;
        else
            matches++;

        if (print && printed++ < MAX_PRINT_ERRORS)
        {
            printname();
            ::print(" = {:08X} {{{:04X}}} (should be {:08X} {{{:04X}}})\n", ourdst, oursw, x87dst, x87sw);
            reex();
        }
    }

    //
    // check a 32-bit floating-point value against a result from the real FPU
    //
    template<typename PrintFuncType, typename ReExecFuncType>
    void check_value(float const &ourdstx, float const &x87dstx, uint16_t oursw, uint16_t x87sw, PrintFuncType printname, ReExecFuncType reex)
    {
        bool print = false;
        count++;

        oursw &= ~X87SW_TOP_MASK;
        x87sw &= ~X87SW_TOP_MASK;

        float_int32_t ourdst{ourdstx};
        float_int32_t x87dst{x87dstx};
        uint32_t mandiff = std::abs(int32_t(x87dst.i - ourdst.i));
        if ((x87dst.i & FP32_EXPONENT_MASK) != (ourdst.i & FP32_EXPONENT_MASK) && mandiff > 16)
            experrors++, print = true;
        else if ((x87dst.i & FP32_MANTISSA_MASK) != (ourdst.i & FP32_MANTISSA_MASK))
        {
            int index = 0;
            for ( ; mandiff != 0; index++, mandiff >>= 1) ;
            errors[index]++;
            print = (index >= print_thresh);
        }
        else if ((x87dst.i & FP32_SIGN_MASK) != (ourdst.i & FP32_SIGN_MASK))
            signerrors++, print = true;
        else if (x87sw != oursw)
            swerrors++, print = true;
        else
            matches++;

        if (print && printed++ < MAX_PRINT_ERRORS)
        {
            printname();
            ::print(" = {:08X} [{:+.12e}] {{{:04X}}} (should be {:08X} [{:+.12e}] {{{:04X}}})\n",
                ourdst.i, double(ourdst.d), oursw,
                x87dst.i, double(x87dst.d), x87sw);
            reex();
        }
    }

    //
    // check a 64-bit floating-point value against a result from the real FPU
    //
    template<typename PrintFuncType, typename ReExecFuncType>
    void check_value(fp64_t const &ourdst, fp64_t const &x87dst, uint16_t oursw, uint16_t x87sw, PrintFuncType printname, ReExecFuncType reex, bool secondary = false)
    {
        bool print = false;
        if (!secondary)
            count++;

        oursw &= ~X87SW_TOP_MASK;
        x87sw &= ~X87SW_TOP_MASK;

        uint64_t mandiff = std::abs(int64_t(x87dst.as_fpbits64() - ourdst.as_fpbits64()));
        if (x87dst.exponent() != ourdst.exponent() && mandiff > 16)
            experrors++, print = true;
        else if (x87dst.mantissa() != ourdst.mantissa())
        {
            int index = 0;
            for ( ; mandiff != 0; index++, mandiff >>= 1) ;
            errors[index]++;
            print = (index >= print_thresh);
        }
        else if (x87dst.sign() != ourdst.sign())
            signerrors++, print = true;
        else if (!secondary && x87sw != oursw)
            swerrors++, print = true;
        else
            matches++;

        if (print && printed++ < MAX_PRINT_ERRORS)
        {
            printname();
            ::print(" = {:016X} [{:+.12e}] {{{:04X}}} (should be {:016X} [{:+.12e}] {{{:04X}}})\n",
                ourdst.as_fpbits64(), ourdst.as_double(), oursw,
                x87dst.as_fpbits64(), x87dst.as_double(), x87sw);
            reex();
        }
    }

    //
    // check an 80-bit floating-point value against a result from the real FPU
    //
    template<typename PrintFuncType, typename ReExecFuncType>
    void check_value(fp80_t const &ourdst, fp80_t const &x87dst, uint16_t oursw, uint16_t x87sw, PrintFuncType printname, ReExecFuncType reex, bool secondary = false)
    {
        bool print = false;
        if (!secondary)
            count++;

        oursw &= ~X87SW_TOP_MASK;
        x87sw &= ~X87SW_TOP_MASK;

        if (x87dst.sign_exp() != ourdst.sign_exp())
            experrors++, print = true;
        else if (x87dst.mantissa() != ourdst.mantissa())
        {
            uint64_t mandiff = std::abs(int64_t(x87dst.mantissa() - ourdst.mantissa()));
            int index = 0;
            for ( ; mandiff != 0; index++, mandiff >>= 1) ;
            errors[index]++;
            print = (index >= print_thresh);
        }
        else if (x87dst.sign() != ourdst.sign())
            signerrors++, print = true;
        else if (!secondary && x87sw != oursw)
            swerrors++, print = true;
        else
            matches++;

        if (print && printed++ < MAX_PRINT_ERRORS)
        {
            printname();
            ::print(" = {:04X}:{:016X} [{:+.12e}] {{{:04X}}} (should be {:04X}:{:016X} [{:+.12e}] {{{:04X}}})\n",
                ourdst.sign_exp(), ourdst.mantissa(), ourdst.as_double(), uint16_t(oursw),
                x87dst.sign_exp(), x87dst.mantissa(), x87dst.as_double(), uint16_t(x87sw));
            reex();
        }
    }

    //
    // print a summary report
    //
    void print_report(char const *name)
    {
        eprint("{} results:\n", name);
        eprint("   {:9} matches [{:.2f}%]\n", matches, double(matches) * 100.0 / double(count));
        for (int index = 1; index < 64; index++)
            if (errors[index] != 0)
                eprint("   {:9} off by {} bits [{:.2f}%]\n", errors[index], index, double(errors[index]) * 100.0 / double(count));
        if (experrors != 0)
            eprint("   {:9} differ by exponent [{:.2f}%]\n", experrors, double(experrors) * 100.0 / double(count));
        if (signerrors != 0)
            eprint("   {:9} differ by sign only [{:.2f}%]\n", signerrors, double(signerrors) * 100.0 / double(count));
        if (infinities != 0)
            eprint("   {:9} pseudo infinities [{:.2f}%]\n", infinities, double(infinities) * 100.0 / double(count));
        if (swerrors != 0)
            eprint("   {:9} status word errors [{:.2f}%]\n", swerrors, double(swerrors) * 100.0 / double(count));
        eprint("\n");
    }
};

//
// create a set of interesting 80-bit float values
//
void make_values80()
{
    // +/-0
    values80.emplace_back(0x0000000000000000ull, 0x0000);
    values80.emplace_back(0x0000000000000000ull, 0x8000);

    // denorms
    for (uint32_t bits = 2; bits < 63; bits++)
    {
        if (bits != 63)
        {
            values80.emplace_back((1ull << bits), 0x0000);
            values80.emplace_back((1ull << bits), 0x8000);
        }
        values80.emplace_back((1ull << bits) - 1, 0x0000);
        values80.emplace_back((1ull << bits) - 1, 0x8000);
        values80.emplace_back(0x123456789abcdefull & ((1ull << bits) - 1), 0x0000);
        values80.emplace_back(0x123456789abcdefull & ((1ull << bits) - 1), 0x8000);
        values80.emplace_back(0xfedcba987654321ull & ((1ull << bits) - 1), 0x0000);
        values80.emplace_back(0xfedcba987654321ull & ((1ull << bits) - 1), 0x8000);
    }

    // normal values80
    static int const s_exponents[] =
    {
                0x0001, 0x0002, 0x0003, 0x0004, 0x0005, 0x0006, 0x0007,
        0x0ffc, 0x0ffd, 0x0ffe, 0x0fff, 0x1000, 0x1001, 0x1002, 0x1003,
        0x1ffc, 0x1ffd, 0x1ffe, 0x1fff, 0x2000, 0x2001, 0x2002, 0x2003,
        0x3fec, 0x3fed, 0x3fee, 0x3fef, 0x3ff0, 0x3ff1, 0x3ff2, 0x3ff3,
        0x3ff4, 0x3ff5, 0x3ff6, 0x3ff7, 0x3ff8, 0x3ff9, 0x3ffa, 0x3ffb,
        0x3ffc, 0x3ffd, 0x3ffe, 0x3fff, 0x4000, 0x4001, 0x4002, 0x4003,
        0x4004, 0x4005, 0x4006, 0x4007, 0x4008, 0x4009, 0x400a, 0x400b,
        0x400c, 0x400d, 0x400e, 0x400f, 0x4010, 0x4011, 0x4012, 0x4013,
        0x4ff8, 0x4ff9, 0x4ffa, 0x4ffb, 0x4ffc, 0x4ffd, 0x4ffe, 0x4fff,
        0x5ffc, 0x5ffd, 0x5ffe, 0x5fff, 0x6000, 0x6001, 0x6002, 0x6003,
        0x7ff8, 0x7ff9, 0x7ffa, 0x7ffb, 0x7ffc, 0x7ffd, 0x7ffe, 0x7fff
    };
    for (int iexp = 0; iexp < (sizeof(s_exponents) / sizeof(s_exponents[0])); iexp++)
    {
        auto exp = s_exponents[iexp];
        values80.emplace_back(0x8000000000000000ull, uint16_t(exp | 0x0000));
        values80.emplace_back(0x8000000000000000ull, uint16_t(exp | 0x8000));
        for (uint32_t bits = 2; bits < 63; bits++)
        {
            if (bits != 63)
            {
                values80.emplace_back(0x8000000000000000ull | (1ull << bits), uint16_t(exp | 0x0000));
                values80.emplace_back(0x8000000000000000ull | (1ull << bits), uint16_t(exp | 0x8000));
            }
            values80.emplace_back(0x8000000000000000ull | ((1ull << bits) - 1), uint16_t(exp | 0x0000));
            values80.emplace_back(0x8000000000000000ull | ((1ull << bits) - 1), uint16_t(exp | 0x8000));
            values80.emplace_back(0x8000000000000000ull | (0x123456789abcdefull & ((1ull << bits) - 1)), uint16_t(exp | 0x0000));
            values80.emplace_back(0x8000000000000000ull | (0x123456789abcdefull & ((1ull << bits) - 1)), uint16_t(exp | 0x8000));
            values80.emplace_back(0x8000000000000000ull | (0xfedcba987654321ull & ((1ull << bits) - 1)), uint16_t(exp | 0x0000));
            values80.emplace_back(0x8000000000000000ull | (0xfedcba987654321ull & ((1ull << bits) - 1)), uint16_t(exp | 0x8000));
        }
    }

    // +/- infinity
    values80.emplace_back(0x8000000000000000ull, 0x7fff);
    values80.emplace_back(0x8000000000000000ull, 0xffff);

    // NaNs
    for (uint32_t bits = 2; bits < 63; bits++)
    {
        if (bits != 63)
        {
            values80.emplace_back(0x8000000000000000ull | (1ull << bits), 0x7fff);
            values80.emplace_back(0x8000000000000000ull | (1ull << bits), 0xffff);
        }
        values80.emplace_back(0x8000000000000000ull | ((1ull << bits) - 1), 0x7fff);
        values80.emplace_back(0x8000000000000000ull | ((1ull << bits) - 1), 0xffff);
        values80.emplace_back(0x8000000000000000ull | (0x123456789abcdefull & ((1ull << bits) - 1)), 0x7fff);
        values80.emplace_back(0x8000000000000000ull | (0x123456789abcdefull & ((1ull << bits) - 1)), 0xffff);
        values80.emplace_back(0x8000000000000000ull | (0xfedcba987654321ull & ((1ull << bits) - 1)), 0x7fff);
        values80.emplace_back(0x8000000000000000ull | (0xfedcba987654321ull & ((1ull << bits) - 1)), 0xffff);
    }
}

//
// create a set of interesting 64-bit float values
//
void make_values64()
{
    constexpr uint64_t p = 0x0000000000000000ull;
    constexpr uint64_t m = 0x8000000000000000ull;

    // +/-0
    values64.push_back(fp64_t::from_fpbits64(p | 0x0000000000000000ull));
    values64.push_back(fp64_t::from_fpbits64(m | 0x0000000000000000ull));

    // denorms
    for (uint32_t bits = 2; bits < 52; bits++)
    {
        if (bits != 52)
        {
            values64.push_back(fp64_t::from_fpbits64(p | (1ull << bits)));
            values64.push_back(fp64_t::from_fpbits64(m | (1ull << bits)));
        }
        values64.push_back(fp64_t::from_fpbits64(p | ((1ull << bits) - 1)));
        values64.push_back(fp64_t::from_fpbits64(m | ((1ull << bits) - 1)));
        values64.push_back(fp64_t::from_fpbits64(p | (0x456789abcdefull & ((1ull << bits) - 1))));
        values64.push_back(fp64_t::from_fpbits64(m | (0x456789abcdefull & ((1ull << bits) - 1))));
        values64.push_back(fp64_t::from_fpbits64(p | (0xfedcba987654ull & ((1ull << bits) - 1))));
        values64.push_back(fp64_t::from_fpbits64(m | (0xfedcba987654ull & ((1ull << bits) - 1))));
    }

    // normal values
    static int const s_exponents[] =
    {
               0x001, 0x002, 0x003, 0x004, 0x005, 0x006, 0x007,
        0x0fc, 0x0fd, 0x0fe, 0x0ff, 0x100, 0x101, 0x102, 0x103,
        0x1fc, 0x1fd, 0x1fe, 0x1ff, 0x200, 0x201, 0x202, 0x203,
        0x3cc, 0x3cd, 0x3ce, 0x3cf, 0x3d0, 0x3d1, 0x3d2, 0x3d3,
        0x3d4, 0x3d5, 0x3d6, 0x3d7, 0x3d8, 0x3d9, 0x3da, 0x3db,
        0x3dc, 0x3dd, 0x3de, 0x3df, 0x3e0, 0x3e1, 0x3e2, 0x3e3,
        0x3e4, 0x3e5, 0x3e6, 0x3e7, 0x3e8, 0x3e9, 0x3ea, 0x3eb,
        0x3ec, 0x3ed, 0x3ee, 0x3ef, 0x3f0, 0x3f1, 0x3f2, 0x3f3,
        0x3f4, 0x3f5, 0x3f6, 0x3f7, 0x3f8, 0x3f9, 0x3fa, 0x3fb,
        0x3fc, 0x3fd, 0x3fe, 0x3ff, 0x400, 0x401, 0x402, 0x403,
        0x404, 0x405, 0x406, 0x407, 0x408, 0x409, 0x40a, 0x40b,
        0x40c, 0x40d, 0x40e, 0x40f, 0x410, 0x411, 0x412, 0x413,
        0x414, 0x415, 0x416, 0x417, 0x418, 0x419, 0x41a, 0x41b,
        0x41c, 0x41d, 0x41e, 0x41f, 0x420, 0x421, 0x422, 0x423,
        0x424, 0x425, 0x426, 0x427, 0x428, 0x429, 0x42a, 0x42b,
        0x42c, 0x42d, 0x42e, 0x42f, 0x430, 0x431, 0x432, 0x433,
        0x43d, 0x43e, 0x43f, 0x440, 0x441, 0x442, 0x443, 0x444,  // for sin/cos/tan
        0x4f8, 0x4f9, 0x4fa, 0x4fb, 0x4fc, 0x4fd, 0x4fe, 0x4ff,
        0x5fc, 0x5fd, 0x5fe, 0x5ff, 0x600, 0x601, 0x602, 0x603,
        0x7f8, 0x7f9, 0x7fa, 0x7fb, 0x7fc, 0x7fd, 0x7fe, 0x7ff
    };
    for (int iexp = 0; iexp < (sizeof(s_exponents) / sizeof(s_exponents[0])); iexp++)
    {
        uint64_t e = uint64_t(s_exponents[iexp]) << 52;
        values64.push_back(fp64_t::from_fpbits64(p | e));
        values64.push_back(fp64_t::from_fpbits64(m | e));
        for (uint32_t bits = 2; bits < 52; bits++)
        {
            if (bits != 52)
            {
                values64.push_back(fp64_t::from_fpbits64(p | e | (1ull << bits)));
                values64.push_back(fp64_t::from_fpbits64(m | e | (1ull << bits)));
            }
            values64.push_back(fp64_t::from_fpbits64(p | e | ((1ull << bits) - 1)));
            values64.push_back(fp64_t::from_fpbits64(m | e | ((1ull << bits) - 1)));
            values64.push_back(fp64_t::from_fpbits64(p | e | (0x3456789abcdefull & ((1ull << bits) - 1))));
            values64.push_back(fp64_t::from_fpbits64(m | e | (0x3456789abcdefull & ((1ull << bits) - 1))));
            values64.push_back(fp64_t::from_fpbits64(p | e | (0xfedcba9876543ull & ((1ull << bits) - 1))));
            values64.push_back(fp64_t::from_fpbits64(m | e | (0xfedcba9876543ull & ((1ull << bits) - 1))));
        }
    }

    // +/- infinity
    values64.push_back(fp64_t::from_fpbits64(p | 0x7ff0000000000000ull));
    values64.push_back(fp64_t::from_fpbits64(m | 0x7ff0000000000000ull));

    // NaNs
    {
        uint64_t e = 0x7ffull << 52;
        values64.push_back(fp64_t::from_fpbits64(p | e | 0x0000000000000001ull));
        values64.push_back(fp64_t::from_fpbits64(m | e | 0x0000000000000001ull));
        values64.push_back(fp64_t::from_fpbits64(p | e | 0x0003ffffffffffffull));
        values64.push_back(fp64_t::from_fpbits64(m | e | 0x0003ffffffffffffull));
        values64.push_back(fp64_t::from_fpbits64(p | e | 0x0004000000000000ull));
        values64.push_back(fp64_t::from_fpbits64(m | e | 0x0004000000000000ull));
        values64.push_back(fp64_t::from_fpbits64(p | e | 0x0007ffffffffffffull));
        values64.push_back(fp64_t::from_fpbits64(m | e | 0x0007ffffffffffffull));
        values64.push_back(fp64_t::from_fpbits64(p | e | 0x0008000000000000ull));
        values64.push_back(fp64_t::from_fpbits64(m | e | 0x0008000000000000ull));
        values64.push_back(fp64_t::from_fpbits64(p | e | 0x000bffffffffffffull));
        values64.push_back(fp64_t::from_fpbits64(m | e | 0x000bffffffffffffull));
        values64.push_back(fp64_t::from_fpbits64(p | e | 0x000c000000000000ull));
        values64.push_back(fp64_t::from_fpbits64(m | e | 0x000c000000000000ull));
        values64.push_back(fp64_t::from_fpbits64(p | e | 0x000fffffffffffffull));
        values64.push_back(fp64_t::from_fpbits64(m | e | 0x000fffffffffffffull));
        values64.push_back(fp64_t::from_fpbits64(p | e | 0x0003456789abcdefull));
        values64.push_back(fp64_t::from_fpbits64(m | e | 0x0003456789abcdefull));
        values64.push_back(fp64_t::from_fpbits64(p | e | 0x000fedcba9876543ull));
        values64.push_back(fp64_t::from_fpbits64(m | e | 0x000fedcba9876543ull));
    }
}

//
// create a set of interesting 32-bit float values
//
void make_values32()
{
    auto from_fpbits32 = [](uint32_t i) { int32_float_t val{i}; return val.d; };

    constexpr uint32_t p = 0x00000000;
    constexpr uint32_t m = 0x80000000;

    // +/-0
    values32.push_back(from_fpbits32(p | 0x00000000));
    values32.push_back(from_fpbits32(m | 0x00000000));

    // denorms
    for (uint32_t bits = 2; bits < 23; bits++)
    {
        if (bits != 23)
        {
            values32.push_back(from_fpbits32(p | (1 << bits)));
            values32.push_back(from_fpbits32(m | (1 << bits)));
        }
        values32.push_back(from_fpbits32(p | ((1 << bits) - 1)));
        values32.push_back(from_fpbits32(m | ((1 << bits) - 1)));
        values32.push_back(from_fpbits32(p | (0x456789ab & ((1 << bits) - 1))));
        values32.push_back(from_fpbits32(m | (0x456789ab & ((1 << bits) - 1))));
        values32.push_back(from_fpbits32(p | (0xfedcba98 & ((1 << bits) - 1))));
        values32.push_back(from_fpbits32(m | (0xfedcba98 & ((1 << bits) - 1))));
    }

    // normal values
    static int const s_exponents[] =
    {
               0x001, 0x002, 0x003, 0x004, 0x005, 0x006, 0x007,
        0x0fc, 0x0fd, 0x0fe, 0x0ff
    };
    for (int iexp = 0; iexp < (sizeof(s_exponents) / sizeof(s_exponents[0])); iexp++)
    {
        uint64_t e = s_exponents[iexp] << 23;
        values32.push_back(from_fpbits32(p | e));
        values32.push_back(from_fpbits32(m | e));
        for (uint32_t bits = 2; bits < 23; bits++)
        {
            if (bits != 52)
            {
                values32.push_back(from_fpbits32(p | e | (1 << bits)));
                values32.push_back(from_fpbits32(m | e | (1 << bits)));
            }
            values32.push_back(from_fpbits32(p | e | ((1 << bits) - 1)));
            values32.push_back(from_fpbits32(m | e | ((1 << bits) - 1)));
            values32.push_back(from_fpbits32(p | e | (0x3456789a & ((1 << bits) - 1))));
            values32.push_back(from_fpbits32(m | e | (0x3456789a & ((1 << bits) - 1))));
            values32.push_back(from_fpbits32(p | e | (0xfedcba98 & ((1 << bits) - 1))));
            values32.push_back(from_fpbits32(m | e | (0xfedcba98 & ((1 << bits) - 1))));
        }
    }

    // +/- infinity
    values32.push_back(from_fpbits32(p | 0x7f800000));
    values32.push_back(from_fpbits32(m | 0x7f800000));

    // NaNs
    {
        uint32_t e = 0xff << 23;
        values32.push_back(from_fpbits32(p | e | 0x00000000));
        values32.push_back(from_fpbits32(m | e | 0x00000000));
        values32.push_back(from_fpbits32(p | e | 0x003fffff));
        values32.push_back(from_fpbits32(m | e | 0x003fffff));
        values32.push_back(from_fpbits32(p | e | 0x00400000));
        values32.push_back(from_fpbits32(m | e | 0x00400000));
        values32.push_back(from_fpbits32(p | e | 0x007fffff));
        values32.push_back(from_fpbits32(m | e | 0x007fffff));
        values32.push_back(from_fpbits32(p | e | 0x00000000));
        values32.push_back(from_fpbits32(m | e | 0x00000000));
        values32.push_back(from_fpbits32(p | e | 0x003fffff));
        values32.push_back(from_fpbits32(m | e | 0x003fffff));
        values32.push_back(from_fpbits32(p | e | 0x00400000));
        values32.push_back(from_fpbits32(m | e | 0x00400000));
        values32.push_back(from_fpbits32(p | e | 0x007fffff));
        values32.push_back(from_fpbits32(m | e | 0x007fffff));
        values32.push_back(from_fpbits32(p | e | 0x00345678));
        values32.push_back(from_fpbits32(m | e | 0x00345678));
        values32.push_back(from_fpbits32(p | e | 0x007edcba));
        values32.push_back(from_fpbits32(m | e | 0x007edcba));
    }
}

//
// create a set of interesting n-bit integer values
//
template<typename Type>
void make_valuesi(std::vector<Type> &values)
{
    values.push_back(0);
    for (uint32_t index = 0; index < 32; index++)
    {
        values.push_back(index + 1);
        values.push_back(-Type(index + 1));
    }
    for (uint32_t bit = 6; bit < 8 * sizeof(Type); bit++)
    {
        Type base = 1 << bit;
        for (int delta = -3; delta <= 3; delta++)
        {
            values.push_back(base + delta);
            values.push_back(-(base + delta));
        }
    }
}

//
// helpers to print out values of various types
//
void print_val(char const *name, fp64_t const &val)
{
    print("{} = {:c}{:013X}e{:+05d} ({:+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.exponent(), val.as_double());
}

void print_val(char const *name, fp80_t const &val)
{
    print("{} = {:c}{:016X}e{:+05d} ({:+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.exponent(), val.as_double());
}

void print_val(char const *name, fpext52_t const &val)
{
    print("{} = {:c}{:013X}e{:+05d} ({:+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.exponent(), val.as_double());
}

void print_val(char const *name, fpext64_t const &val)
{
    print("{} = {:c}{:016X}e{:+05d} ({:+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.exponent(), val.as_double());
}

void print_val(char const *name, fpext96_t const &val)
{
    print("{} = {:c}{:016X}`{:08X}e{:+05d} ({:+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.extend(), val.exponent(), val.as_double());
}

//
// validate conversions between different types
//
void validate_conversions()
{
    for (int index = 0; index < values64.size(); index++)
    {
        auto src = values64[index];
        if (src.isnormal() || src.isdenorm())
        {
            fpext64_t temp80(src);
            auto src80 = temp80.as_fp64();
            if (src.as_fpbits64() != src80.as_fpbits64())
                print("64-bit: {:016X} -> {}.{}:{:016X} -> {:016X}\n",
                    src.as_fpbits64(),
                    temp80.sign(), temp80.exponent(), temp80.mantissa(),
                    src80.as_fpbits64());

            fpext96_t temp96(src);
            auto src96 = temp96.as_fp64();
            if (src.as_fpbits64() != src96.as_fpbits64())
                print("64-bit: {:016X} -> {}.{}:{:016X}{:04X} -> {:016X}\n",
                    src.as_fpbits64(),
                    temp96.sign(), temp96.exponent(), temp96.mantissa(), temp96.extend(),
                    src96.as_fpbits64());

            fpext96_t temp112(src);
            auto src112 = temp112.as_fp64();
            if (src.as_fpbits64() != src112.as_fpbits64())
                print("64-bit: {:016X} -> {}.{}:{:016X}{:04X} -> {:016X}\n",
                    src.as_fpbits64(),
                    temp112.sign(), temp112.exponent(), temp112.mantissa(), temp112.extend(),
                    src112.as_fpbits64());
        }
    }

    for (int index = 0; index < values80.size(); index++)
    {
        auto src = values80[index];
        if (src.isnormal() || src.isdenorm())
        {
            fpext64_t temp80(src);
            auto src80 = temp80.as_fp80();
            if (src.sign_exp() != src80.sign_exp() || src.mantissa() != src80.mantissa())
                print("80-bit: {:04X}:{:016X} -> {}.{}:{:016X} -> {:04X}:{:016X}\n",
                    src.sign_exp(), src.mantissa(),
                    temp80.sign(), temp80.exponent(), temp80.mantissa(),
                    src80.sign_exp(), src80.mantissa());

            fpext96_t temp96(src);
            auto src96 = temp96.as_fp80();
            if (src.sign_exp() != src96.sign_exp() || src.mantissa() != src96.mantissa())
                print("80-bit: {:04X}:{:016X} -> {}.{}:{:016X}{:04X} -> {:04X}:{:016X}\n",
                    src.sign_exp(), src.mantissa(),
                    temp96.sign(), temp96.exponent(), temp96.mantissa(), temp96.extend(),
                    src96.sign_exp(), src96.mantissa());

            fpext96_t temp112(src);
            auto src112 = temp112.as_fp80();
            if (src.sign_exp() != src112.sign_exp() || src.mantissa() != src112.mantissa())
                print("80-bit: {:04X}:{:016X} -> {}.{}:{:016X}{:04X} -> {:04X}:{:016X}\n",
                    src.sign_exp(), src.mantissa(),
                    temp112.sign(), temp112.exponent(), temp112.mantissa(), temp112.extend(),
                    src112.sign_exp(), src112.mantissa());
        }
    }

    // conversion
    int errors = 0;
    for (int sign = 0; sign < 2; sign++)
        for (int exp = 0; exp <= 0x7ff; exp++)
            for (uint64_t man = 0; man < 0x10000000000000ull; man = (man == 0) ? 1 : ((man << 2) | 1))
            {
                auto src = fp64_t::from_fpbits64((uint64_t(sign) << 63) | (uint64_t(exp) << 52) | man);

                fp80_t ourext(src);

                fp80_t x87ext;
                x87test1(&src, &x87ext);

                if (ourext != x87ext)
                {
                    if (++errors < MAX_PRINT_ERRORS)
                        print("Src = {:016X} ({:+.9e}) -> {:04X}:{:016X} (should be {:04X}:{:016X})\n",
                            src.as_fpbits64(), src.as_double(), ourext.sign_exp(), ourext.mantissa(), x87ext.sign_exp(), x87ext.mantissa());
                }

                fp64_t ourret(ourext);

                fp64_t x87ret;
                x87test2(&x87ext, &x87ret);

                if (ourret.as_fpbits64() != x87ret.as_fpbits64())
                {
                    if (++errors < MAX_PRINT_ERRORS)
                        print("Src = {:04X}:{:016X} -> {:016X} ({:f}) (should be {:016X} ({:f}))\n",
                            ourext.sign_exp(), ourext.mantissa(), ourret.as_fpbits64(), ourret.as_double(), x87ret.as_fpbits64(), x87ret.as_double());
                }
            }
}

//
// test a unary 64-bit operation
//
template<typename FpFuncType, typename X87FuncType>
void test_unary64(FpFuncType fpfunc, X87FuncType x87func, char const *name, int print_thresh)
{
    errors_t errs(name, print_thresh);
    for (int isrc1 = 0; isrc1 < values64.size(); isrc1++)
    {
        fp64_t src1(values64[isrc1]);

        fp64_t x87dst;
        auto x87sw = x87func(&src1, &x87dst) & ~X87SW_TOP_MASK;

        fp64_t ourdst;
        auto oursw = fpfunc(src1, ourdst) & ~X87SW_TOP_MASK;

        errs.check_value(ourdst, x87dst, oursw, x87sw,
            [&]() { print("{}({:016X} [{:+.12e}])", name, src1.as_fpbits64(), src1.as_double()); },
            [&]()
            {
                fp64_t res;
                fpfunc(src1, res);
            });
    }
    auto start = std::chrono::steady_clock::now();
    size_t reps = 0;
    std::chrono::nanoseconds elapsed;
    do
    {
        for (int isrc1 = 0; isrc1 < values64.size(); isrc1++)
        {
            fp64_t ourdst;
            fpfunc(values64[isrc1], ourdst);
        }
        reps += values64.size();
        elapsed = std::chrono::steady_clock::now() - start;
    } while (elapsed < min_timing_duration);
    eprint("{}: ns/op = {:.2f}\n", name, double(elapsed.count()) / double(reps));
    errs.print_report(name);
}

//
// test a unary 64-bit operation with two results
//
template<typename FpFuncType, typename X87FuncType>
void test_unary64_2(FpFuncType fpfunc, X87FuncType x87func, char const *name, int print_thresh)
{
    errors_t errs(name, print_thresh);
    for (int isrc1 = 0; isrc1 < values64.size(); isrc1++)
    {
        fp64_t src1(values64[isrc1]);

        fp64_t x87dst1, x87dst2;
        auto x87sw = x87func(&src1, &x87dst1, &x87dst2) & ~X87SW_TOP_MASK;

        fp64_t ourdst1, ourdst2;
        auto oursw = fpfunc(src1, ourdst1, ourdst2) & ~X87SW_TOP_MASK;

        errs.check_value(ourdst1, x87dst1, oursw, x87sw,
            [&]() { print("{}({:016X} [{:+.12e}])[1]", name, src1.as_fpbits64(), src1.as_double()); },
            [&]()
            {
                fp64_t res1, res2;
                fpfunc(src1, res1, res2);
            });
        errs.check_value(ourdst2, x87dst2, oursw, oursw,
            [&]() { print("{}({:016X} [{:+.12e}])[2]", name, src1.as_fpbits64(), src1.as_double()); },
            [&]()
            {
                fp64_t res1, res2;
                fpfunc(src1, res1, res2);
            }, true);
    }
    auto start = std::chrono::steady_clock::now();
    size_t reps = 0;
    std::chrono::nanoseconds elapsed;
    do
    {
        for (int isrc1 = 0; isrc1 < values64.size(); isrc1++)
        {
            fp64_t ourdst1, ourdst2;
            fpfunc(values64[isrc1], ourdst1, ourdst2);
        }
        reps += values64.size();
        elapsed = std::chrono::steady_clock::now() - start;
    } while (elapsed < min_timing_duration);
    eprint("{}: ns/op = {:.2f}\n", name, double(elapsed.count()) / double(reps));
    errs.print_report(name);
}

//
// test a binary 64-bit operation
//
template<typename FpFuncType, typename X87FuncType>
void test_binary64(FpFuncType fpfunc, X87FuncType x87func, char const *name, int print_thresh)
{
    errors_t errs(name, print_thresh);
    for (int isrc1 = 0; isrc1 < values64.size(); isrc1 += 5)
    {
        fp64_t src1(values64[isrc1]);
        for (int isrc2 = 0; isrc2 < values64.size(); isrc2 += 5)
        {
            fp64_t src2(values64[isrc2]);

            fp64_t x87dst;
            auto x87sw = x87func(&src2, &src1, &x87dst) & ~X87SW_TOP_MASK;

            fp64_t ourdst;
            auto oursw = fpfunc(src2, src1, ourdst) & ~X87SW_TOP_MASK;

            errs.check_value(ourdst, x87dst, oursw, x87sw,
                [&]() { print("{}({:016X} [{:+.12e}], {:016X} [{:+.12e}])", name, src2.as_fpbits64(), src2.as_double(), src1.as_fpbits64(), src1.as_double()); },
                [&]()
                {
                    fp64_t res;
                    fpfunc(src2, src1, res);
                });
        }
    }
    auto start = std::chrono::steady_clock::now();
    size_t reps = 0;
    int initial = 0;
    std::chrono::nanoseconds elapsed;
    do
    {
        for (int isrc1 = ++initial; isrc1 < values64.size(); isrc1 += 23)
            for (int isrc2 = initial + 1; isrc2 < values64.size(); isrc2 += 17)
            {
                fp64_t ourdst;
                fpfunc(values64[isrc2], values64[isrc1], ourdst);
            }
        reps += ((values64.size() - initial) / 23) * ((values64.size() - initial - 1) / 17);
        elapsed = std::chrono::steady_clock::now() - start;
    } while (elapsed < min_timing_duration);
    eprint("{}: ns/op = {:.2f}\n", name, double(elapsed.count()) / double(reps));
    errs.print_report(name);
}

//
// test a unary 80-bit operation
//
template<typename FpFuncType, typename X87FuncType>
void test_unary80(FpFuncType fpfunc, X87FuncType x87func, char const *name, int print_thresh)
{
    errors_t errs(name, print_thresh);
    for (int isrc1 = 0; isrc1 < values80.size(); isrc1++)
    {
        fp80_t src1(values80[isrc1]);

        fp80_t x87dst;
        auto x87sw = x87func(&src1, &x87dst) & ~X87SW_TOP_MASK;

        fp80_t ourdst;
        auto oursw = fpfunc(src1, ourdst) & ~X87SW_TOP_MASK;

        errs.check_value(ourdst, x87dst, oursw, x87sw,
            [&]() { print("{}({:04X}:{:016X} [{:+.12e}])", name, src1.sign_exp(), src1.mantissa(), src1.as_double()); },
            [&]() { fp80_t res; fpfunc(src1, res); });
    }
    auto start = std::chrono::steady_clock::now();
    size_t reps = 0;
    std::chrono::nanoseconds elapsed;
    do
    {
        for (int isrc1 = 0; isrc1 < values80.size(); isrc1++)
        {
            fp80_t ourdst;
            fpfunc(values80[isrc1], ourdst);
        }
        reps += values80.size();
        elapsed = std::chrono::steady_clock::now() - start;
    } while (elapsed < min_timing_duration);
    eprint("{}: ns/op = {:.2f}\n", name, double(elapsed.count()) / double(reps));
    errs.print_report(name);
}

//
// test an 80-bit unary operation with two results
//
template<typename FpFuncType, typename X87FuncType>
void test_unary80_2(FpFuncType fpfunc, X87FuncType x87func, char const *name, int print_thresh)
{
    errors_t errs(name, print_thresh);
    for (int isrc1 = 0; isrc1 < values80.size(); isrc1++)
    {
        fp80_t src1(values80[isrc1]);

        fp80_t x87dst1, x87dst2;
        auto x87sw = x87func(&src1, &x87dst1, &x87dst2) & ~X87SW_TOP_MASK;

        fp80_t ourdst1, ourdst2;
        auto oursw = fpfunc(src1, ourdst1, ourdst2) & ~X87SW_TOP_MASK;

        errs.check_value(ourdst1, x87dst1, oursw, x87sw,
            [&]() { print("{}({:04X}:{:016X} [{:+.12e}])[1]", name, src1.sign_exp(), src1.mantissa(), src1.as_double()); },
            [&]() { fp80_t r1, r2; fpfunc(src1, r1, r2); });
        errs.check_value(ourdst2, x87dst2, oursw, oursw,
            [&]() { print("{}({:04X}:{:016X} [{:+.12e}])[2]", name, src1.sign_exp(), src1.mantissa(), src1.as_double()); },
            [&]() { fp80_t r1, r2; fpfunc(src1, r1, r2); }, true);
    }
    auto start = std::chrono::steady_clock::now();
    size_t reps = 0;
    std::chrono::nanoseconds elapsed;
    do
    {
        for (int isrc1 = 0; isrc1 < values80.size(); isrc1++)
        {
            fp80_t ourdst1, ourdst2;
            fpfunc(values80[isrc1], ourdst1, ourdst2);
        }
        reps += values80.size();
        elapsed = std::chrono::steady_clock::now() - start;
    } while (elapsed < min_timing_duration);
    eprint("{}: ns/op = {:.2f}\n", name, double(elapsed.count()) / double(reps));
    errs.print_report(name);
}

//
// test an 80-bit binary operation
//
template<typename FpFuncType, typename X87FuncType>
void test_binary80(FpFuncType fpfunc, X87FuncType x87func, char const *name, int print_thresh)
{
    errors_t errs(name, print_thresh);
    for (int isrc1 = 0; isrc1 < values80.size(); isrc1 += 5)
    {
        fp80_t src1(values80[isrc1]);
        for (int isrc2 = 0; isrc2 < values80.size(); isrc2 += 5)
        {
            fp80_t src2(values80[isrc2]);

            fp80_t x87dst;
            auto x87sw = x87func(&src2, &src1, &x87dst) & ~X87SW_TOP_MASK;

            fp80_t ourdst;
            auto oursw = fpfunc(src2, src1, ourdst) & ~X87SW_TOP_MASK;

            errs.check_value(ourdst, x87dst, oursw, x87sw,
                [&]() { print("{}({:04X}:{:016X} [{:+.12e}], {:04X}:{:016X} [{:+.12e}])",
                    name,
                    src2.sign_exp(), src2.mantissa(), src2.as_double(),
                    src1.sign_exp(), src1.mantissa(), src1.as_double()); },
                [&]() { fp80_t r; fpfunc(src2, src1, r); });
        }
    }
    auto start = std::chrono::steady_clock::now();
    size_t reps = 0;
    int initial = 0;
    std::chrono::nanoseconds elapsed;
    do
    {
        for (int isrc1 = ++initial; isrc1 < values80.size(); isrc1 += 23)
            for (int isrc2 = initial + 1; isrc2 < values80.size(); isrc2 += 17)
            {
                fp80_t ourdst;
                fpfunc(values80[isrc2], values80[isrc1], ourdst);
            }
        reps += ((values80.size() - initial) / 23) * ((values80.size() - initial - 1) / 17);
        elapsed = std::chrono::steady_clock::now() - start;
    } while (elapsed < min_timing_duration);
    eprint("{}: ns/op = {:.2f}\n", name, double(elapsed.count()) / double(reps));
    errs.print_report(name);
}

//
// 80-bit wrappers around the free-function arithmetic operators so they can
// be passed to test_binary80 (which expects a uniform (a, b, &dst) -> sw API).
//
inline uint16_t fp80_add_wrap(fp80_t const &a, fp80_t const &b, fp80_t &dst) { dst = a + b; return 0; }
inline uint16_t fp80_sub_wrap(fp80_t const &a, fp80_t const &b, fp80_t &dst) { dst = a - b; return 0; }
inline uint16_t fp80_mul_wrap(fp80_t const &a, fp80_t const &b, fp80_t &dst) { dst = a * b; return 0; }
inline uint16_t fp80_div_wrap(fp80_t const &a, fp80_t const &b, fp80_t &dst) { dst = a / b; return 0; }
inline uint16_t fp80_sqrt_wrap(fp80_t const &a, fp80_t &dst) { dst = fp80_t::sqrt(a); return 0; }

//
// Test that fp80_t's built-in static constants match the real FPU. The first
// seven slots come from fld1/fldl2t/fldl2e/fldpi/fldlg2/fldln2/fldz; the next
// six are derived (nzero via fldz+fchs, pinf/ninf/indef via 1/0,-1/0,0/0) or
// memory-roundtripped (snan, qnan) — see x87consts80 in x87testasm.asm.
//
void test_consts80()
{
    errors_t errs("consts(80)", 0);
    fp80_t fpu[13];
    x87consts80(fpu);
    fp80_t ours[13] = {
        fp80_t::const_one(),
        fp80_t::const_l2t(),
        fp80_t::const_l2e(),
        fp80_t::const_pi(),
        fp80_t::const_lg2(),
        fp80_t::const_ln2(),
        fp80_t::const_zero(),
        fp80_t::const_nzero(),
        fp80_t::const_pinf(),
        fp80_t::const_ninf(),
        fp80_t::const_snan(),
        fp80_t::const_qnan(),
        fp80_t::const_indef(),
    };
    static char const * const names[13] = {
        "one", "l2t", "l2e", "pi", "lg2", "ln2", "zero",
        "nzero", "pinf", "ninf", "snan", "qnan", "indef"
    };
    for (int i = 0; i < 13; i++)
    {
        errs.check_value(ours[i], fpu[i], 0, 0,
            [&]() { print("const_{}", names[i]); },
            [&]() {});
    }
    errs.print_report("consts(80)");
}

//
// Test fp80_t's classification predicates against FXAM.
//
// Per Intel x87 spec, Table 8-1, FXAM returns:
//   C3:C2:C0 = 000  Unsupported
//             001  NaN
//             010  Normal finite
//             011  Infinity
//             100  Zero
//             101  Empty (stack slot)
//             110  Denormal
//             111  Empty
// and C1 = sign of the operand.
//
// We test each predicate is consistent with the FPU's classification. isqnan
// vs issnan can't be distinguished by FXAM (both are NaN to FXAM), so we
// derive the expected from the bit pattern (mantissa MSB after the explicit
// one) — but only over values FXAM agrees are NaN, so a mis-classifying NaN
// would still be flagged.
//
void test_predicates80()
{
    errors_t errs("preds(80)", 0);
    int const C3 = X87SW_C3, C2 = X87SW_C2, C0 = X87SW_C0, C1 = X87SW_C1;
    for (int i = 0; i < int(values80.size()); i++)
    {
        fp80_t const &v = values80[i];
        uint16_t sw = fxam80(&v);
        bool fxam_nan       = (sw & (C3|C2|C0)) == C0;
        bool fxam_normfin   = (sw & (C3|C2|C0)) == C2;
        bool fxam_inf       = (sw & (C3|C2|C0)) == (C2|C0);
        bool fxam_zero      = (sw & (C3|C2|C0)) == C3;
        bool fxam_denorm    = (sw & (C3|C2|C0)) == (C3|C2);
        bool fxam_unsup     = (sw & (C3|C2|C0)) == 0;
        bool fxam_sign      = (sw & C1) != 0;

        // Build a list of (name, ours, expected) tuples and check each.
        struct entry { char const *name; bool ours; bool expected; };
        // fp80_t::isnormal() returns true for normal *or* zero (i.e., bits are
        // already in canonical form). The FPU has no single class for this;
        // it's the union of "normal finite" and "zero".
        bool exp_isnormal = fxam_normfin || fxam_zero;
        bool exp_ispinf   = fxam_inf && !fxam_sign;
        bool exp_isninf   = fxam_inf && fxam_sign;
        // qnan vs snan: only meaningful when FXAM agrees this is a NaN.
        bool exp_isqnan = fxam_nan && ((v.mantissa() & FP80_MANTISSA_MASK) >= 0x4000000000000000ull);
        bool exp_issnan = fxam_nan && !exp_isqnan;

        entry checks[] = {
            { "isnan",    v.isnan(),    fxam_nan    },
            { "isinf",    v.isinf(),    fxam_inf    },
            { "iszero",   v.iszero(),   fxam_zero   },
            { "isdenorm", v.isdenorm(), fxam_denorm },
            { "isnormal", v.isnormal(), exp_isnormal },
            { "ispinf",   v.ispinf(),   exp_ispinf  },
            { "isninf",   v.isninf(),   exp_isninf  },
            { "isqnan",   v.isqnan(),   exp_isqnan  },
            { "issnan",   v.issnan(),   exp_issnan  },
            // "sign()" is exposed as a uint8_t getter; verify against FXAM C1.
            // Unsupported encodings don't have a meaningful FXAM sign, so we
            // only check when FXAM gave us a classification.
        };

        for (auto const &c : checks)
        {
            errs.count++;
            if (c.ours != c.expected)
            {
                errs.experrors++;
                if (errs.printed++ < MAX_PRINT_ERRORS)
                    print("{}({:04X}:{:016X}) = {} (FXAM sw={:04X}, expected {})\n",
                        c.name, v.sign_exp(), v.mantissa(),
                        c.ours ? "true" : "false", sw,
                        c.expected ? "true" : "false");
            }
            else
                errs.matches++;
        }

        // sign() consistency, but only for non-unsupported encodings
        if (!fxam_unsup)
        {
            bool exp_sign = fxam_sign;
            bool our_sign = (v.sign() != 0);
            errs.count++;
            if (our_sign != exp_sign)
            {
                errs.experrors++;
                if (errs.printed++ < MAX_PRINT_ERRORS)
                    print("sign({:04X}:{:016X}) = {} (FXAM C1 says {})\n",
                        v.sign_exp(), v.mantissa(),
                        our_sign ? "true" : "false",
                        exp_sign ? "true" : "false");
            }
            else
                errs.matches++;
        }
    }
    errs.print_report("preds(80)");
}

//
// Sign manipulation: abs() clears the sign, chs() flips it, copysign(a,b)
// transplants b's sign onto a, samesign() reports whether two values agree.
//
// fabs / fchs map to FPU instructions of the same name. copysign / samesign
// have no x87 instruction; the oracle is bit-pattern manipulation of the
// fp80_t encoding, which is well-defined regardless of the FPU.
//
inline uint16_t fp80_abs_wrap(fp80_t const &a, fp80_t &dst) { dst = fp80_t::abs(a); return 0; }
inline uint16_t fp80_chs_wrap(fp80_t const &a, fp80_t &dst) { dst = fp80_t::chs(a); return 0; }

//
// floor() and ceil() are defined as round-to-int with rounding mode Down
// and Up respectively. The asm side uses frndint with an appropriately
// preconfigured control word.
//
inline uint16_t fp80_floor_wrap(fp80_t const &a, fp80_t &dst) { dst = fp80_t::floor(a); return 0; }
inline uint16_t fp80_ceil_wrap (fp80_t const &a, fp80_t &dst) { dst = fp80_t::ceil(a);  return 0; }

void test_copysign80()
{
    errors_t errs("copysign(80)", 0);
    // Sample sparsely; copysign is purely bit ops so we don't need full N^2.
    for (int i = 0; i < int(values80.size()); i += 11)
    {
        fp80_t const &a = values80[i];
        for (int j = 0; j < int(values80.size()); j += 11)
        {
            fp80_t const &b = values80[j];

            // Oracle: |a| with b's sign bit. Reconstruct from raw fields so
            // we don't rely on the very function under test.
            fp80_t expected(a.mantissa(),
                            (a.sign_exp() & ~FP80_SIGN_MASK) | (b.sign_exp() & FP80_SIGN_MASK));

            fp80_t ours = a;
            ours.copysign(b);

            errs.check_value(ours, expected, 0, 0,
                [&]() { print("copysign({:04X}:{:016X}, {:04X}:{:016X})",
                    a.sign_exp(), a.mantissa(),
                    b.sign_exp(), b.mantissa()); },
                [&]() {});
        }
    }
    errs.print_report("copysign(80)");
}

void test_samesign80()
{
    errors_t errs("samesign(80)", 0);
    for (int i = 0; i < int(values80.size()); i += 11)
    {
        fp80_t const &a = values80[i];
        for (int j = 0; j < int(values80.size()); j += 11)
        {
            fp80_t const &b = values80[j];
            bool expected = ((a.sign_exp() ^ b.sign_exp()) & FP80_SIGN_MASK) == 0;
            bool ours = fp80_t::samesign(a, b);
            errs.count++;
            if (ours != expected)
            {
                errs.experrors++;
                if (errs.printed++ < MAX_PRINT_ERRORS)
                    print("samesign({:04X}:{:016X}, {:04X}:{:016X}) = {} (expected {})\n",
                        a.sign_exp(), a.mantissa(),
                        b.sign_exp(), b.mantissa(),
                        ours ? "true" : "false",
                        expected ? "true" : "false");
            }
            else
                errs.matches++;
        }
    }
    errs.print_report("samesign(80)");
}

void test_make_qnan80()
{
    errors_t errs("make_qnan(80)", 0);
    for (int i = 0; i < int(values80.size()); i++)
    {
        fp80_t const &v = values80[i];
        if (!v.isnan()) continue;  // make_qnan precondition: input must be a NaN
        fp80_t q = fp80_t::make_qnan(v);
        errs.count++;
        // Spec for fp80_t::make_qnan: keep all other bits, force the "quiet"
        // bit (mantissa MSB below the explicit 1) on.
        if (!q.isqnan()) {
            errs.experrors++;
            if (errs.printed++ < MAX_PRINT_ERRORS)
                print("make_qnan({:04X}:{:016X}) = {:04X}:{:016X} (not qnan)\n",
                    v.sign_exp(), v.mantissa(), q.sign_exp(), q.mantissa());
            continue;
        }
        // Bits other than the quiet bit should be preserved.
        uint64_t expected_man = v.mantissa() | 0xc000000000000000ull;
        if (q.mantissa() != expected_man || q.sign_exp() != v.sign_exp())
        {
            errs.experrors++;
            if (errs.printed++ < MAX_PRINT_ERRORS)
                print("make_qnan({:04X}:{:016X}) = {:04X}:{:016X} (expected {:04X}:{:016X})\n",
                    v.sign_exp(), v.mantissa(),
                    q.sign_exp(), q.mantissa(),
                    v.sign_exp(), expected_man);
        }
        else
            errs.matches++;
    }
    errs.print_report("make_qnan(80)");
}

//
// Test an SW-only fp80 op (FXAM-style — no result value, only flags). For
// each value in values80 the C++ implementation's returned SW is compared
// directly against what the real CPU's FPU produces.
//
template<typename FpFuncType, typename X87FuncType>
void test_classify80(FpFuncType fpfunc, X87FuncType x87func, char const *name)
{
    errors_t errs(name, 0);
    for (int i = 0; i < int(values80.size()); i++)
    {
        fp80_t const &v = values80[i];
        uint16_t x87sw = x87func(&v) & ~X87SW_TOP_MASK;
        uint16_t oursw = fpfunc(v) & ~X87SW_TOP_MASK;
        errs.count++;
        if (oursw != x87sw)
        {
            errs.swerrors++;
            if (errs.printed++ < MAX_PRINT_ERRORS)
                print("{}({:04X}:{:016X}) sw = {:04X} (should be {:04X})",
                    name, v.sign_exp(), v.mantissa(), oursw, x87sw);
        }
        else
            errs.matches++;
    }
    errs.print_report(name);
}

//
// Test an SW-returning binary compare op (fcom/fucom/fcomi/fucomi). Like
// test_classify80 but takes two operands. Strides through the value table
// to keep the run time bounded.
//
template<typename FpFuncType, typename X87FuncType>
void test_compare_sw80(FpFuncType fpfunc, X87FuncType x87func, char const *name, int stride = 7)
{
    errors_t errs(name, 0);
    for (int i = 0; i < int(values80.size()); i += stride)
    {
        fp80_t const &a = values80[i];
        for (int j = 0; j < int(values80.size()); j += stride)
        {
            fp80_t const &b = values80[j];
            uint16_t x87sw = x87func(&a, &b) & ~X87SW_TOP_MASK;
            uint16_t oursw = fpfunc(a, b) & ~X87SW_TOP_MASK;
            errs.count++;
            if (oursw != x87sw)
            {
                errs.swerrors++;
                if (errs.printed++ < MAX_PRINT_ERRORS)
                    print("{}({:04X}:{:016X}, {:04X}:{:016X}) sw = {:04X} (should be {:04X})",
                        name,
                        a.sign_exp(), a.mantissa(),
                        b.sign_exp(), b.mantissa(),
                        oursw, x87sw);
            }
            else
                errs.matches++;
        }
    }
    errs.print_report(name);
}

//
// Compare two fp80 status words for the cmp(80) test. Currently fp80_t has
// no x87_fcom that returns a SW, so this validates the *test harness's* model
// (test_compare80 derives booleans from C0/C2/C3 — verify those are the only
// bits fucompp produces from the FPU on these inputs).
//
// Specifically: fucompp's only meaningful output bits beyond TOP and ES are
// C0, C2, C3 (and C1=0 unless stack fault, which can't happen since we
// always reinit). Anything else in the low SW bits would indicate an issue
// in the asm bridge.
//
void test_fcmp_sw_shape()
{
    errors_t errs("cmp_sw_shape(80)", 0);
    // Bits we DON'T expect fucompp to set: precision/underflow/overflow/divzero;
    // and we do NOT expect C1 since we always reinit before each call.
    // (Note: fucompp legitimately sets DE on denormal operands and IE on SNaN
    // inputs, so those are excluded from FORBIDDEN.)
    uint16_t const FORBIDDEN =
        X87SW_PRECISION_EX | X87SW_UNDERFLOW_EX | X87SW_OVERFLOW_EX |
        X87SW_DIVZERO_EX   | X87SW_C1;
    for (int i = 0; i < int(values80.size()); i += 7)
    {
        for (int j = 0; j < int(values80.size()); j += 7)
        {
            uint16_t sw = fcmp80(&values80[i], &values80[j]) & ~X87SW_TOP_MASK;
            errs.count++;
            // FUCOMPP sets #IA only when an operand is an SNaN (QNaNs are
            // "unordered" without raising IE). Skip IE checks: they're
            // expected for SNaN inputs.
            uint16_t leak = sw & FORBIDDEN;
            if (leak != 0)
            {
                errs.experrors++;
                if (errs.printed++ < MAX_PRINT_ERRORS)
                    print("cmp_sw_shape({:04X}:{:016X}, {:04X}:{:016X}) sw={:04X} (unexpected bits: {:04X})\n",
                        values80[i].sign_exp(), values80[i].mantissa(),
                        values80[j].sign_exp(), values80[j].mantissa(),
                        sw, leak);
            }
            else
                errs.matches++;
        }
    }
    errs.print_report("cmp_sw_shape(80)");
}

//
// Test the six fp80 comparison operators against fucompp from the real FPU.
//
// The FPU encodes the comparison result in C3:C2:C0:
//    a > b      -> 0:0:0
//    a < b      -> 0:0:1
//    a == b     -> 1:0:0
//    unordered  -> 1:1:1   (either operand is NaN)
//
void test_compare80()
{
    errors_t errs("cmp(80)", 0);
    // O(N^2) over values80 is too expensive at the full set; stride to keep
    // the count tractable while still hitting every kind of value.
    int const stride = 7;
    for (int i = 0; i < int(values80.size()); i += stride)
    {
        fp80_t const &a = values80[i];
        for (int j = 0; j < int(values80.size()); j += stride)
        {
            fp80_t const &b = values80[j];

            uint16_t sw = fcmp80(&a, &b);
            bool x87_gt = (sw & (X87SW_C3 | X87SW_C2 | X87SW_C0)) == 0;
            bool x87_lt = (sw & (X87SW_C3 | X87SW_C2 | X87SW_C0)) == X87SW_C0;
            bool x87_eq = (sw & (X87SW_C3 | X87SW_C2 | X87SW_C0)) == X87SW_C3;
            bool x87_unord = (sw & X87SW_C2) != 0;

            // Derived expected boolean for each C++ operator. Standard
            // IEEE comparison semantics: every ordered relation is false
            // when either operand is NaN; != is the lone exception.
            struct { char const *name; bool ours; bool expected; } cases[6] = {
                { "==", a == b, x87_eq                          },
                { "!=", a != b, !x87_eq                         },
                { "<",  a <  b, x87_lt && !x87_unord            },
                { "<=", a <= b, (x87_lt || x87_eq) && !x87_unord},
                { ">",  a >  b, x87_gt && !x87_unord            },
                { ">=", a >= b, (x87_gt || x87_eq) && !x87_unord},
            };

            for (auto const &c : cases)
            {
                if (c.ours != c.expected)
                {
                    errs.experrors++;
                    if (errs.printed++ < MAX_PRINT_ERRORS)
                        print("cmp({:04X}:{:016X}, {:04X}:{:016X}) {} = {} (should be {})\n",
                            a.sign_exp(), a.mantissa(),
                            b.sign_exp(), b.mantissa(),
                            c.name,
                            c.ours ? "true" : "false",
                            c.expected ? "true" : "false");
                }
                else
                    errs.matches++;
                errs.count++;
            }
        }
    }
    errs.print_report("cmp(80)");
}

//
// test a load operation
//
template<typename DstType, typename FpFuncType, typename X87FuncType, typename SrcType>
void test_load(FpFuncType fpfunc, X87FuncType x87func, std::vector<SrcType> const &vals, char const *name, int print_thresh = 0)
{
    errors_t errs(name, print_thresh);
    for (int isrc1 = 0; isrc1 < vals.size(); isrc1++)
    {
        DstType x87dst;
        auto x87sw = x87func(&vals[isrc1], &x87dst);

        DstType ourdst;
        auto oursw = fpfunc(&vals[isrc1], &ourdst);

        errs.check_value(ourdst, x87dst, oursw, x87sw,
            [&]() { print("{}(", name); for (uint32_t index = 0; index < sizeof(vals[isrc1]); index++) print("{:02X}", ((uint8_t const *)&vals[isrc1])[sizeof(vals[isrc1]) - 1 - index]); print(")"); },
            [&]() { DstType res; fpfunc(&vals[isrc1], &res); });
    }
    errs.print_report(name);
}

//
// test a store operation
//
template<typename DstType, typename FpFuncType, typename X87FuncType, typename SrcType>
void test_store(FpFuncType fpfunc, X87FuncType x87func, std::vector<SrcType> const &vals, char const *name, int print_thresh = 0)
{
    errors_t errs(name, print_thresh);
    for (int isrc1 = 0; isrc1 < vals.size(); isrc1++)
    {
        DstType x87dst;
        auto x87sw = x87func(&vals[isrc1], &x87dst);

        DstType ourdst;
        auto oursw = fpfunc(&vals[isrc1], &ourdst);

        errs.check_value(ourdst, x87dst, oursw, x87sw,
            [&]() { print("{}({})", name, vals[isrc1]); },
            [&]() { DstType res; fpfunc(&vals[isrc1], &res); });
    }
    errs.print_report(name);
}

//
// Control word:
//   bits 11-10 = rounding control
//      00 = nearest
//      01 = -infinity
//      10 = +infinity
//      11 = toward zero
//   bits 9-8 = precision control
//      00 = single
//      01 = reserved
//      10 = double
//      11 = extended
//

int main(int argc, char *argv[])
{
    int errors = 0;

    min_timing_duration = std::chrono::milliseconds(500);

    uint16_t cw = 0xf3f;
    x87setcw(&cw);

    make_values32();
    make_values64();
    make_values80();
    make_valuesi(valuesi64);
    make_valuesi(valuesi32);
    make_valuesi(valuesi16);

    validate_conversions();

    static std::array<x87cw_t, 3> const s_precision = { X87CW_PRECISION_EXTENDED, X87CW_PRECISION_DOUBLE, X87CW_PRECISION_SINGLE };
    static std::array<x87cw_t, 4> const s_round = { X87CW_ROUNDING_NEAREST, X87CW_ROUNDING_DOWN, X87CW_ROUNDING_UP, X87CW_ROUNDING_ZERO };

    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            x87cw_t cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);

            print("Testing precision {} round {}\n", prec, round);

            test_load<fp80_t>(
                [&](auto const *src, auto *dst) { x87sw_t sw = 0; fp80_t::x87_fld80(cw, sw, *dst, src); return sw; },
                fld8080, values80, "fld80");
        }

    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            x87cw_t cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);

            print("Testing precision {} round {}\n", prec, round);

            test_load<fp80_t>(
                [&](auto const *src, auto *dst) { x87sw_t sw = 0; fp80_t::x87_fld64(cw, sw, *dst, src); return sw; },
                fld6480, values64, "fld64");
        }

    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            x87cw_t cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);

            print("Testing precision {} round {}\n", prec, round);

            test_load<fp80_t>(
                [&](auto const *src, auto *dst) { x87sw_t sw = 0; fp80_t::x87_fld32(cw, sw, *dst, src); return sw; },
                fld3280, values32, "fld32");
        }

    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            x87cw_t cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);

            print("Testing precision {} round {}\n", prec, round);

            test_load<fp80_t>(
                [&](auto const *src, auto *dst) { x87sw_t sw = 0; fp80_t::x87_fild64(cw, sw, *dst, src); return sw; },
                fild6480, valuesi64, "fild64");
        }

    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            x87cw_t cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);

            print("Testing precision {} round {}\n", prec, round);

            test_load<fp80_t>(
                [&](auto const *src, auto *dst) { x87sw_t sw = 0; fp80_t::x87_fild32(cw, sw, *dst, src); return sw; },
                fild3280, valuesi32, "fild32");
        }

    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            x87cw_t cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);

            print("Testing precision {} round {}\n", prec, round);

            test_load<fp80_t>(
                [&](auto const *src, auto *dst) { x87sw_t sw = 0; fp80_t::x87_fild16(cw, sw, *dst, src); return sw; },
                fild1680, valuesi16, "fild16");
        }


    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            x87cw_t cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);

            print("Testing precision {} round {}\n", prec, round);

            test_store<fp80_t>(
                [&](auto const *src, auto *dst) { x87sw_t sw = 0; fp80_t::x87_fst80(cw, sw, dst, *src); return sw; },
                fst8080, values80, "fst80");
        }

    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            x87cw_t cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);

            print("Testing precision {} round {}\n", prec, round);

            test_store<fp64_t>(
                [&](auto const *src, auto *dst) { x87sw_t sw = 0; fp80_t::x87_fst64(cw, sw, dst, *src); return sw; },
                fst8064, values80, "fst64");
        }

    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            x87cw_t cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);

            print("Testing precision {} round {}\n", prec, round);

            test_store<float>(
                [&](auto const *src, auto *dst) { x87sw_t sw = 0; fp80_t::x87_fst32(cw, sw, dst, *src); return sw; },
                fst8032, values80, "fst32");
        }

    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            x87cw_t cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);

            print("Testing precision {} round {}\n", prec, round);

            test_store<int64_t>(
                [&](auto const *src, auto *dst) { x87sw_t sw = 0; fp80_t::x87_fist64(cw, sw, dst, *src); return sw; },
                fist8064, values80, "fist64");
        }

    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            x87cw_t cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);

            print("Testing precision {} round {}\n", prec, round);

            test_store<int32_t>(
                [&](auto const *src, auto *dst) { x87sw_t sw = 0; fp80_t::x87_fist32(cw, sw, dst, *src); return sw; },
                fist8032, values80, "fist32");
        }

    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            x87cw_t cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);

            print("Testing precision {} round {}\n", prec, round);

            test_store<int16_t>(
                [&](auto const *src, auto *dst) { x87sw_t sw = 0; fp80_t::x87_fist16(cw, sw, dst, *src); return sw; },
                fist8016, values80, "fist16");
        }

    // set round: to zero, precision: 53 bits
    cw = X87CW_MASK_ALL_EX | X87CW_ROUNDING_ZERO | X87CW_PRECISION_DOUBLE;
    x87setcw(&cw);

    test_unary64_2(&fp64_t::x87_fxtract, &fxtract64, "fxtract(64)", 1);
    test_unary64(&fp64_t::x87_f2xm1, &f2xm164, "f2xm1(64)", 2);
    test_unary64(&fp64_t::x87_fsin, &fsin64, "fsin(64)", 3);
    test_unary64(&fp64_t::x87_fcos, &fcos64, "fcos(64)", 3);
    test_unary64_2(&fp64_t::x87_fsincos, &fsincos64, "fsincos(64)", 3);
    test_unary64_2(&fp64_t::x87_fptan, &fptan64, "fptan(64)", 3);

    test_binary64(&fp64_t::x87_fscale, &fscale64, "fscale(64)", 1);
    test_binary64(&fp64_t::x87_fprem, &fprem64, "fprem(64)", 1);
    test_binary64(&fp64_t::x87_fprem1, &fprem164, "fprem1(64)", 1);
    test_binary64(&fp64_t::x87_fyl2xp1, &fyl2xp164, "fyl2xp1(64)", 3);
    test_binary64(&fp64_t::x87_fyl2x, &fyl2x64, "fyl2x(64)", 2);
    test_binary64(&fp64_t::x87_fpatan, &fpatan64, "fpatan(64)", 3);

    // set round: to nearest, precision: 64 bits
    cw = X87CW_MASK_ALL_EX | X87CW_ROUNDING_NEAREST | X87CW_PRECISION_EXTENDED;
    x87setcw(&cw);

    // Built-in constants vs. the real FPU (all 13 constants)
    test_consts80();

    // FXAM-driven classification: isnan/isinf/iszero/isdenorm/isnormal/
    // ispinf/isninf/isqnan/issnan/sign
    test_predicates80();

    // Sign manipulation: abs, chs, copysign, samesign, make_qnan
    test_unary80(&fp80_abs_wrap, &fabs80, "fabs(80)", 1);
    test_unary80(&fp80_chs_wrap, &fchs80, "fchs(80)", 1);
    test_copysign80();
    test_samesign80();
    test_make_qnan80();

    // Comparison operators (==, !=, <, <=, >, >=) against fucompp
    test_compare80();
    // Verify nothing leaks into the fucompp status word besides C0/C2/C3
    test_fcmp_sw_shape();

    // SW-returning classification and compare entry points: validate that
    // C++ x87_fxam / x87_ftst / x87_fcom / x87_fucom / x87_fcomi / x87_fucomi
    // each produce exactly the SW the live FPU does. (Currently all stubs.)
    test_classify80(&fp80_t::x87_fxam, &fxam80, "x87_fxam(80)");
    test_classify80(&fp80_t::x87_ftst, &ftst80, "x87_ftst(80)");
    test_compare_sw80(&fp80_t::x87_fcom,   &fcom80,   "x87_fcom(80)");
    test_compare_sw80(&fp80_t::x87_fucom,  &fcmp80,   "x87_fucom(80)");
    test_compare_sw80(&fp80_t::x87_fcomi,  &fcomi80,  "x87_fcomi(80)");
    test_compare_sw80(&fp80_t::x87_fucomi, &fucomi80, "x87_fucomi(80)");

    // 80-bit arithmetic: sweep across all 4 rounding modes x 3 precision
    // modes, because per the spec (Table 8-2) PC affects FADD/FSUB/FMUL/FDIV
    // (+ reverse variants) and FSQRT. All other transcendentals use full
    // 80-bit internally and are unaffected by PC, so we run them at default
    // CW only.
    //
    // These tests use the SW-returning entry points (x87_fadd etc.) so
    // C1=roundup and the precision/underflow/overflow/invalid exception
    // bits are validated against the real FPU rather than ignored.
    for (uint32_t prec = 0; prec < s_precision.size(); prec++)
        for (uint32_t round = 0; round < s_round.size(); round++)
        {
            cw = X87CW_MASK_ALL_EX | s_precision[prec] | s_round[round];
            x87setcw(&cw);
            print("80-bit arithmetic: precision {} round {}\n", prec, round);
            test_binary80(&fp80_t::x87_fadd,  &fadd80,  "fadd(80)",  1);
            test_binary80(&fp80_t::x87_fsub,  &fsub80,  "fsub(80)",  1);
            test_binary80(&fp80_t::x87_fsubr, &fsubr80, "fsubr(80)", 1);
            test_binary80(&fp80_t::x87_fmul,  &fmul80,  "fmul(80)",  1);
            test_binary80(&fp80_t::x87_fdiv,  &fdiv80,  "fdiv(80)",  1);
            test_binary80(&fp80_t::x87_fdivr, &fdivr80, "fdivr(80)", 1);
            test_unary80 (&fp80_t::x87_fsqrt, &fsqrt80, "fsqrt(80)", 1);
        }

    // Back to default CW (round-to-nearest, extended precision) for the
    // remaining tests.
    cw = X87CW_MASK_ALL_EX | X87CW_ROUNDING_NEAREST | X87CW_PRECISION_EXTENDED;
    x87setcw(&cw);

    // floor/ceil are defined as frndint with rounding mode Down/Up respectively.
    cw = X87CW_MASK_ALL_EX | X87CW_ROUNDING_DOWN | X87CW_PRECISION_EXTENDED;
    x87setcw(&cw);
    test_unary80(&fp80_floor_wrap, &frndint80, "floor(80)", 1);
    cw = X87CW_MASK_ALL_EX | X87CW_ROUNDING_UP | X87CW_PRECISION_EXTENDED;
    x87setcw(&cw);
    test_unary80(&fp80_ceil_wrap, &frndint80, "ceil(80)", 1);
    cw = X87CW_MASK_ALL_EX | X87CW_ROUNDING_NEAREST | X87CW_PRECISION_EXTENDED;
    x87setcw(&cw);

    // 80-bit transcendentals (PC has no effect — single mode is sufficient)
    test_unary80_2(&fp80_t::x87_fxtract, &fxtract80, "fxtract(80)", 1);
    test_unary80(&fp80_t::x87_f2xm1, &f2xm180, "f2xm1(80)", 2);
    test_unary80(&fp80_t::x87_fsin, &fsin80, "fsin(80)", 3);
    test_unary80(&fp80_t::x87_fcos, &fcos80, "fcos(80)", 3);
    test_unary80_2(&fp80_t::x87_fsincos, &fsincos80, "fsincos(80)", 3);
    test_unary80_2(&fp80_t::x87_fptan, &fptan80, "fptan(80)", 3);
    test_unary80(&fp80_t::x87_frndint, &frndint80, "frndint(80)", 1);

    test_binary80(&fp80_t::x87_fscale, &fscale80, "fscale(80)", 1);
    test_binary80(&fp80_t::x87_fprem, &fprem80, "fprem(80)", 1);
    test_binary80(&fp80_t::x87_fprem1, &fprem180, "fprem1(80)", 1);
    test_binary80(&fp80_t::x87_fyl2xp1, &fyl2xp180, "fyl2xp1(80)", 3);
    test_binary80(&fp80_t::x87_fyl2x, &fyl2x80, "fyl2x(80)", 2);
    test_binary80(&fp80_t::x87_fpatan, &fpatan80, "fpatan(80)", 3);

    return 0;
}

#if 0
{
    int sign = 0;
    for (int srcval = 0; srcval < 128; srcval++)
    {
        int srcexp = 0x7fe;
        fp64_t src1 = fp64_t::from_fpbits64((uint64_t(sign) << 63) | (uint64_t(srcexp) << 52) | srcval);
        int over = -1, under = 1;
        for (int exp = 1; exp < 1000000 && (over < 0 || under > 0); exp++)
        {
            if (over < 0)
            {
                fp64_t src2(exp), dst;
                auto sw = fscale64(&src1, &src2, &dst);
                if (sw & X87SW_OVERFLOW_EX)
                    over = exp;
            }
            if (under > 0)
            {
                fp64_t src2(-exp), dst;
                auto sw = fscale64(&src1, &src2, &dst);
                if (sw & X87SW_UNDERFLOW_EX)
                    under = -exp;
            }
        }
        int thresh = -16394 - ((srcval == 0) ? 52 : ctz64(srcval));
        print("{:016X}: over={} under={} ({})\n", src1.as_fpbits64(), over + (srcexp - 0x3ff), under + (srcexp - 0x3ff), thresh);
    }
    return 0;

    for (int srcexp = 0; srcexp < 0x1000; srcexp++)
    {
        fp64_t src1 = fp64_t::from_fpbits64((uint64_t(srcexp) << 52) | 0x1000);
        int over = -1, under = 1;
        for (int exp = 1; exp < 1000000 && (over < 0 || under > 0); exp++)
        {
            if (over < 0)
            {
                fp64_t src2(exp), dst;
                auto sw = fscale64(&src1, &src2, &dst);
                if (sw & X87SW_OVERFLOW_EX)
                    over = exp;
            }
            if (under > 0)
            {
                fp64_t src2(-exp), dst;
                auto sw = fscale64(&src1, &src2, &dst);
                if (sw & X87SW_UNDERFLOW_EX)
                    under = -exp;
            }
        }
        print("{:03X}: over={} under={}\n", srcexp, over, under);
    }
    return 0;

    for (int exp1 = 0x100; exp1 < 0x500; exp1++)
    {
        fp64_t src2 = fp64_t::from_fpbits64((uint64_t(exp1) << 52) | 0x123456789abull);
        fp64_t src1 = fp64_t::from_fpbits64((uint64_t(0x400) << 52) | 0xba987654321ull);
        fp64_t dst;
        uint16_t sw = fscale64(&src1, &src2, &dst);
        fp64_t ourdst;
        uint16_t oursw = fp64_t::x87_fscale(src1, src2, ourdst);
        uint64_t idst = dst.as_fpbits64();
        print("({}) {:016X} / {:016X} = {:016X} [{:+.16e}] {{{:04X}}} = {:016X} [{:+.16e}] {{{:04X}}}\n", (sw >> X87SW_C2_BIT) & 1, src1.as_fpbits64(), src2.as_fpbits64(), dst.as_fpbits64(), dst.as_double(), sw, ourdst.as_fpbits64(), ourdst.as_double(), oursw);
    }
    print("=============================\n");
    for (int exp1 = 0; exp1 < 0x500; exp1++)
    {
        fp64_t src2 = fp64_t::from_fpbits64((uint64_t(exp1) << 52) | 0x123456789abull);
        fp64_t src1 = fp64_t::from_fpbits64((uint64_t(0) << 52) | 0x8ull);
        fp64_t dst;
        uint16_t sw = fscale64(&src1, &src2, &dst);
        fp64_t ourdst;
        uint16_t oursw = fp64_t::x87_fscale(src1, src2, ourdst);
        uint64_t idst = dst.as_fpbits64();
        print("({}) {:016X} / {:016X} = {:016X} [{:+.16e}] {{{:04X}}} = {:016X} [{:+.16e}] {{{:04X}}}\n", (sw >> X87SW_C2_BIT) & 1, src1.as_fpbits64(), src2.as_fpbits64(), dst.as_fpbits64(), dst.as_double(), sw, ourdst.as_fpbits64(), ourdst.as_double(), oursw);
    }


    return 0;
}
#endif

#if 0
{
    auto val = fp64_t::from_fpbits64(0x3FEA876CCC000000ull);
    fp64_t res1, res2;
    fp64_t::x87_fsin(val, res1);
    return 0;
}
#endif

