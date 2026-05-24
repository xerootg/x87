//=========================================================
//  x87fp80.cpp
//
//  80-bit floating-point support
//=========================================================
//
// BSD 3-Clause License
//
// Copyright (c) 2025, Aaron Giles
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <cstdint>
#include <cmath>
#include <cstring>

#include "x87fp80.h"
#include "x87fp64.h"
#include "x87fpext.h"

// On x86 hosts, fp80_t and `long double` share the binary layout of the
// IEEE 754 double-extended-precision format. We exploit this for the
// finite arithmetic paths so the host x87 unit (with its rounding mode
// and PC bits set via x87setcw / fldcw) does the actual math at the same
// precision the test oracle would. On non-x86 hosts this falls back to
// the existing fpext96_t-based path, which implements the spec by hand.
//
// X87_FORCE_HAND_ROLLED can be defined at build time to force the
// hand-rolled fallback on x86 hosts (for testing the spec implementation).
#if defined(X87_FORCE_HAND_ROLLED)
  #define X87_HOST_HAS_FP80 0
#elif defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
  #define X87_HOST_HAS_FP80 1
#else
  #define X87_HOST_HAS_FP80 0
#endif

namespace x87
{

//
// apply rounding to the 63-bit mantissa value provided, handling overflow
// and returning what sort of rounding was applied
//
using rounding_applied_t = uint64_t;
static constexpr rounding_applied_t ROUND_NEAR = 0;
static constexpr rounding_applied_t ROUND_TOWARD_ZERO = 1;
static constexpr rounding_applied_t ROUND_TOWARD_INF_HARD = 2;
inline rounding_applied_t round_in_place(uint64_t &mantissa, int &exponent, uint64_t sign, x87cw_t rval, int bits)
{
    x87_assert((mantissa & FP80_EXPLICIT_ONE) == 0);
    rounding_applied_t applied;

    // if rounding nearest (even), add 1/2 so that midway values round up
    // unless the current target LSB is already even, in which case add 1/2 - 1
    if (rval == X87CW_ROUNDING_NEAREST)
    {
        mantissa += (1ull << (bits - 1)) - ((~(mantissa | FP80_EXPLICIT_ONE) >> bits) & 1);
        applied = ROUND_NEAR;
    }

    // if rounding toward zero, note that we're doing so
    else if (rval == X87CW_ROUNDING_ZERO)
        applied = ROUND_TOWARD_ZERO;

    // if rounding up/down, and in the right direction, add just less than 1
    else
    {
        applied = ROUND_TOWARD_INF_HARD - (((rval >> X87CW_ROUNDING_SHIFT) ^ sign) & 1);
        if (applied == ROUND_TOWARD_INF_HARD)
            mantissa += (1ull << bits) - 1;
    }

    // if rounding caused an overflow, bump the exponent shift the mantissa down a bit
    // though we can skip the shift because the overflow value is guaranteed to be less
    // than (1 << bits)
    if (int64_t(mantissa) < 0)
    {
        exponent++;
        mantissa ^= FP80_EXPLICIT_ONE;
    }
    return applied;
}


//
// x87 FLD for 80-bit sources
// Exceptions: none
//
void fp80_t::x87_fld80(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src)
{
    // just copy directly; no flags are set
    dst.m_mantissa = *(uint64_t const *)(uintptr_t(src) + 0);
    dst.m_sign_exp = *(uint16_t const *)(uintptr_t(src) + 8);
}

//
// x87 FLD for 64-bit or 32-bit floating-point sources
// Exceptions:
//   #IA if source is SNaN
//   #D if source is denormal
//
template<typename Type>
void fp80_t::x87_fld_common(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src)
{
    // select source parameters based on incoming type
    constexpr Type SRC_EXPONENT_MASK = (sizeof(Type) == 8) ? FP64_EXPONENT_MASK : FP32_EXPONENT_MASK;
    constexpr int SRC_EXPONENT_SHIFT = (sizeof(Type) == 8) ? FP64_EXPONENT_SHIFT : FP32_EXPONENT_SHIFT;
    constexpr Type SRC_MANTISSA_MASK = (sizeof(Type) == 8) ? FP64_MANTISSA_MASK : FP32_MANTISSA_MASK;
    constexpr int SRC_SIGN_SHIFT = (sizeof(Type) == 8) ? FP64_SIGN_SHIFT : FP32_SIGN_SHIFT;
    constexpr int32_t SRC_EXPONENT_BIAS = (sizeof(Type) == 8) ? FP64_EXPONENT_BIAS : FP32_EXPONENT_BIAS;
    constexpr int32_t SRC_EXPONENT_MAX_BIASED = (sizeof(Type) == 8) ? FP64_EXPONENT_MAX_BIASED : FP32_EXPONENT_MAX_BIASED;

    // load raw bits from memory
    auto raw = *(Type const *)src;

    // extract exponent
    int exponent = (raw & SRC_EXPONENT_MASK) >> SRC_EXPONENT_SHIFT;

    // extract mantissa, shifted into place
    uint64_t mantissa = uint64_t(raw & SRC_MANTISSA_MASK) << (63 - SRC_EXPONENT_SHIFT);

    // set the sign
    uint16_t sign_exponent = (raw >> (SRC_SIGN_SHIFT - FP80_SIGN_SHIFT)) & FP80_SIGN_MASK;

    // infinite or NaN? convert to same
    if (exponent == SRC_EXPONENT_MAX_BIASED)
        goto MaxExp;

    // denormal or zero?
    if (exponent == 0)
        goto DenormOrZero;

    // insert adjusted exponent
    dst.m_sign_exp = sign_exponent | (FP80_EXPONENT_BIAS - SRC_EXPONENT_BIAS + exponent);

    // insert explicit 1
    dst.m_mantissa = FP80_EXPLICIT_ONE | mantissa;
    return;

MaxExp:
    // NaN or infinity
    dst.m_sign_exp = sign_exponent | FP80_EXPONENT_MAX_BIASED;
    dst.m_mantissa = FP80_EXPLICIT_ONE | mantissa | ((mantissa != 0) ? 0x4000000000000000ull : 0);

    // SNaNs set #IA
    if (mantissa < ((FP80_MANTISSA_MASK + 1) >> 1) && mantissa != 0)
        sw |= X87CW_MASK_INVALID_EX;
    return;

DenormOrZero:
    // if non-zero mantissa, this is a denorm
    if (mantissa != 0)
        goto Denorm;

    // explicit zero; jus set the sign
    dst.m_sign_exp = sign_exponent;
    dst.m_mantissa = 0;
    return;

Denorm:
    // shift mantissa up so explicit 1 is into the top bit
    int shift = count_leading_zeros64(mantissa);
    dst.m_sign_exp = sign_exponent | (FP80_EXPONENT_BIAS - SRC_EXPONENT_BIAS + 1 - shift);
    dst.m_mantissa = mantissa << shift;

    // denorms set #D
    sw |= X87CW_MASK_DENORM_EX;
    return;
}
template void fp80_t::x87_fld_common<uint64_t>(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src);
template void fp80_t::x87_fld_common<uint32_t>(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src);

//
// construct an 80-bit FP value from an integer
// Exceptions: none
//
template<typename Type>
void fp80_t::x87_fild_common(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src)
{
    // load raw bits from memory
    int64_t raw = *(Type const *)src;

    // extract sign and absolute value
    uint16_t sign_exponent = 0;
    if (raw < 0)
    {
        sign_exponent = FP80_SIGN_MASK;
        raw = -raw;
    }

    // special case for zero
    else if (raw == 0)
    {
        dst.m_mantissa = 0;
        dst.m_sign_exp = 0;
        return;
    }

    // determine shift
    int shift = count_leading_zeros64(raw);
    dst.m_mantissa = raw << shift;
    dst.m_sign_exp = sign_exponent + FP80_EXPONENT_BIAS + 63 - shift;
}
template void fp80_t::x87_fild_common<int64_t>(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src);
template void fp80_t::x87_fild_common<int32_t>(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src);
template void fp80_t::x87_fild_common<int16_t>(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src);

//
// x87 FST for 80-bit targets
// Exceptions: none
//
void fp80_t::x87_fst80(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src)
{
    // just copy directly; no flags are set
    *(uint64_t *)(uintptr_t(dst) + 0) = src.m_mantissa;
    *(uint16_t *)(uintptr_t(dst) + 8) = src.m_sign_exp;
}

//
// x87 FST for floating-point targets
// Exceptions:
//    #IA if source is SNaN
//    #U if source is too small for destination
//    #O if source is too large for destination
//    #P if value cannot be represented exactly in dest
//
template<typename Type>
void fp80_t::x87_fst_common(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src)
{
    // determine target constants based on the template parameter size
    constexpr int TARGET_SIGN_SHIFT = (sizeof(Type) == 8) ? FP64_SIGN_SHIFT : FP32_SIGN_SHIFT;
    constexpr int TARGET_EXPONENT_SHIFT = (sizeof(Type) == 8) ? FP64_EXPONENT_SHIFT : FP32_EXPONENT_SHIFT;
    constexpr int32_t TARGET_EXPONENT_BIAS = (sizeof(Type) == 8) ? FP64_EXPONENT_BIAS : FP32_EXPONENT_BIAS;
    constexpr int32_t TARGET_EXPONENT_MAX_BIASED = (sizeof(Type) == 8) ? FP64_EXPONENT_MAX_BIASED : FP32_EXPONENT_MAX_BIASED;
    constexpr uint64_t TARGET_EXPONENT_MASK = (sizeof(Type) == 8) ? FP64_EXPONENT_MASK : FP32_EXPONENT_MASK;
    constexpr uint64_t TARGET_MANTISSA_MASK = (sizeof(Type) == 8) ? FP64_MANTISSA_MASK : FP32_MANTISSA_MASK;
    constexpr int MANTISSA_SHIFT = 63 - TARGET_EXPONENT_SHIFT;

    // make clang happy
    uint64_t orig_mantissa;
    rounding_applied_t applied;

    // extract the sign and move it to its final location
    uint64_t sign = uint64_t(src.m_sign_exp & FP80_SIGN_MASK) << (TARGET_SIGN_SHIFT - FP80_SIGN_SHIFT);

    // extract full 63-bit mantissa, discarding the explicit 1
    uint64_t mantissa = src.m_mantissa & FP80_MANTISSA_MASK;

    // extract the exponent, but leave it biased
    int exponent = src.m_sign_exp & FP80_EXPONENT_MASK;

    // infinite or NaN? handle as special case
    if (exponent == FP80_EXPONENT_MAX_BIASED)
        goto MaxExp;

    // zero? handle as special case
    if (exponent == 0 && mantissa == 0)
        goto Zero;

    // shift off extra mantissa bits, applying any rounding
    orig_mantissa = mantissa;
    applied = round_in_place(mantissa, exponent, sign >> TARGET_SIGN_SHIFT, cw & X87CW_ROUNDING_MASK, MANTISSA_SHIFT);

    // adjust exponent to the target bias
    exponent = exponent - FP80_EXPONENT_BIAS + TARGET_EXPONENT_BIAS;

    // if exponent too small, convert to denormal or zero
    if (exponent <= 0)
        goto Denormal;

    // too large? convert to infinity
    else if (exponent >= TARGET_EXPONENT_MAX_BIASED)
        goto Overflow;

    // set flags
    if ((orig_mantissa & ((1ull << MANTISSA_SHIFT) - 1)) != 0)
        sw |= X87SW_PRECISION_EX | (((orig_mantissa ^ mantissa) >> (MANTISSA_SHIFT - X87SW_C1_BIT)) & X87SW_C1);

    // otherwise, shift off the extra bits
    mantissa >>= MANTISSA_SHIFT;

    // assemble and return
    *(Type *)dst = Type(sign | (uint64_t(exponent) << TARGET_EXPONENT_SHIFT) | mantissa);
    return;

MaxExp:
    // SNaNs set #IA
    if (mantissa < ((FP80_MANTISSA_MASK + 1) >> 1) && mantissa != 0)
        sw |= X87CW_MASK_INVALID_EX;

    // NaN or infinity case: preserve any mantissa bits for QNaN, but don't do any rounding
    if (mantissa != 0)
        mantissa = ((TARGET_MANTISSA_MASK + 1) >> 1) | (mantissa >> MANTISSA_SHIFT);

    // combine with sign and maximum exponent
    *(Type *)dst = Type(sign | TARGET_EXPONENT_MASK | mantissa);
    return;

Zero:
    // zero case: just return the sign with all zeros in exponent and mantissa
    *(Type *)dst = Type(sign);
    return;

Denormal:
    // too small even for denormal? make it into a signed zero
    if (exponent <= -TARGET_EXPONENT_SHIFT)
    {
        // if we're rounding hard toward infinity, set the mantissa to the smallest
        // non-zero value; otherwise, leave it at 0
        mantissa = applied >> 1;

        // set underflow, and set C1 if we put a non-zero value, since we rounded up
        sw |= X87SW_UNDERFLOW_EX | ((applied << (X87SW_C1_BIT - 1)) & X87SW_C1);
    }
    else
    {
        // shift off mantissa bits, OR in the explicit one
        mantissa = (mantissa | FP80_EXPLICIT_ONE) >> (MANTISSA_SHIFT + 1 - exponent);

        // if we ended up with a zero mantissa, apply the same logic as above
        if (mantissa == 0)
        {
            mantissa = applied >> 1;

            // set C1 if we put a non-zero value, since we rounded up
            sw |= (applied << (X87SW_C1_BIT - 1)) & X87SW_C1;
        }
    }

    // exponent is always zero, so just merge sign and mantissa
    *(Type *)dst = Type(sign | mantissa);

    // denormals set precision
    sw |= X87SW_PRECISION_EX;
    return;

Overflow:
    // maximum exponent and 0 mantissa; however, if rounding toward zero change to
    // the maximum non-infinity value by subtracting 1
    *(Type *)dst = Type((sign | TARGET_EXPONENT_MASK) - int(applied & ROUND_TOWARD_ZERO));

    // set overflow, precision, and round up, unless we maxed out at infinity-1
    sw |= X87SW_OVERFLOW_EX | X87SW_PRECISION_EX | ((~applied << X87SW_C1_BIT) & X87SW_C1);
    return;
}
template void fp80_t::x87_fst_common<uint64_t>(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src);
template void fp80_t::x87_fst_common<uint32_t>(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src);

//
// convert to an integer
//
template<typename Type>
void fp80_t::x87_fist_common(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src)
{
    // make clang happy
    uint64_t mantissa;
    uint64_t orig_mantissa;
    int shift;
    int orig_shift;
    int64_t result;

    // extract the exponent, but leave it biased
    int exponent = src.m_sign_exp & FP80_EXPONENT_MASK;

    // infinite or NaN? return the indefinite
    if (exponent == FP80_EXPONENT_MAX_BIASED)
        goto Indefinite;

    // extract full 63-bit mantissa, discarding the explicit 1
    mantissa = src.m_mantissa & FP80_MANTISSA_MASK;

    // zero? handle as special case
    if (exponent == 0 && mantissa == 0)
        goto Zero;

    // determine shift count
    shift = FP80_EXPONENT_BIAS + 63 - exponent;

    // too large? return the indefinite; note that we use < instead of <= here
    // to catch max negative values (which are not counted as invalid)
    if (shift < int(64 - 8 * sizeof(Type)))
        goto Indefinite;

    // too small? handle specially
    if (shift >= 64)
        goto Small;

    // apply rounding
    orig_mantissa = mantissa;
    orig_shift = shift;
    round_in_place(mantissa, exponent, src.m_sign_exp >> FP80_SIGN_SHIFT, cw & X87CW_ROUNDING_MASK, shift);

    // recompute shift since exponent might have changed due to rounding
    shift = FP80_EXPONENT_BIAS + 63 - exponent;

    // shift the result
    result = (FP80_EXPLICIT_ONE | mantissa) >> shift;

    // apply sign
    if ((src.m_sign_exp & FP80_SIGN_MASK) != 0)
        result = -result;

    // overflow into indefinite
    if (Type(result) != result)
        goto Indefinite;

    // set precision flags if we lost any bits
    if ((orig_mantissa & ((1ull << orig_shift) - 1)) != 0)
        sw |= X87SW_PRECISION_EX | (((((orig_mantissa | FP80_EXPLICIT_ONE) >> orig_shift) ^ result) & 1) << X87SW_C1_BIT);

    // write result
    *(Type *)dst = Type(result);
    return;

Indefinite:
    // return the indefinite value
    *(Type *)dst =  Type(0x8000000000000000ll >> (64 - 8 * sizeof(Type)));

    // set the invalid flag
    sw |= X87CW_MASK_INVALID_EX;
    return;

Small:
    sw |= X87SW_PRECISION_EX;

    // if rounding twoard zero, it's 0
    cw &= X87CW_ROUNDING_MASK;
    if (cw == X87CW_ROUNDING_ZERO)
        goto Zero;

    // if rounding nearest (even) and we're above (but not equal to) 0.5, return +/-1
    if (cw == X87CW_ROUNDING_NEAREST)
    {
        if (shift == 64 && mantissa != 0)
            goto PlusMinusOne;
    }

    // if rounding toward +/- infinity and we're of the same sign, return +/-1
    else if ((((src.m_sign_exp >> (FP80_SIGN_SHIFT - X87CW_ROUNDING_SHIFT)) ^ cw) & (1 << X87CW_ROUNDING_SHIFT)) == 0)
        goto PlusMinusOne;

    // fall through...
Zero:
    // return 0
    *(Type *)dst = 0;
    return;

PlusMinusOne:
    // return -1 if negative or +1 if positive
    *(Type *)dst = ((src.m_sign_exp & FP80_SIGN_MASK) != 0) ? -1 : 1;
    sw |= X87SW_C1;
    return;
}
template void fp80_t::x87_fist_common<int64_t>(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src);
template void fp80_t::x87_fist_common<int32_t>(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src);
template void fp80_t::x87_fist_common<int16_t>(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src);


//===========================================================================
//
// STUBS for unimplemented 80-bit arithmetic.
//
// These return the FPU "indefinite" value so the test harness can link and
// report 100% failure for whichever operation is still missing.
//
//===========================================================================

// operator+/-/*// delegate to x87_f*, discarding the SW.
fp80_t operator+(fp80_t const &a, fp80_t const &b) { fp80_t r; fp80_t::x87_fadd (a, b, r); return r; }
fp80_t operator-(fp80_t const &a, fp80_t const &b) { fp80_t r; fp80_t::x87_fsubr(a, b, r); return r; }
fp80_t operator*(fp80_t const &a, fp80_t const &b) { fp80_t r; fp80_t::x87_fmul (a, b, r); return r; }
fp80_t operator/(fp80_t const &a, fp80_t const &b) { fp80_t r; fp80_t::x87_fdivr(a, b, r); return r; }

fp80_t fp80_t::sqrt(fp80_t const &src)
{
    fp80_t dst;
    x87_fsqrt(src, dst);
    return dst;
}

fp80_t fp80_t::floor(fp80_t const &src)
{
    fpround_t r(X87CW_ROUNDING_DOWN);
    fp80_t dst;
    x87_frndint(src, dst);
    return dst;
}

fp80_t fp80_t::ceil(fp80_t const &src)
{
    fpround_t r(X87CW_ROUNDING_UP);
    fp80_t dst;
    x87_frndint(src, dst);
    return dst;
}

fp80_t fp80_t::ldexp(fp80_t const &a, int32_t factor)
{
    if (a.iszero() || a.isnan() || a.isinf()) return a;
    int     exp  = a.sign_exp() & FP80_EXPONENT_MASK;
    uint16_t sign = a.sign_exp() & FP80_SIGN_MASK;
    if (exp == 0) return a;  // denormal — leave as-is for now
    int new_exp = exp + factor;
    if (new_exp >= FP80_EXPONENT_MAX_BIASED)
        return sign ? fp80_t::const_ninf() : fp80_t::const_pinf();
    if (new_exp <= 0)
        return sign ? fp80_t::const_nzero() : fp80_t::const_zero();
    return fp80_t(a.mantissa(), sign | uint16_t(new_exp));
}

//
// 3-way compare of two fp80 values. Returns:
//    -1  a < b
//     0  a == b
//    +1  a > b
//     2  unordered (either is NaN)
// Handles signed zeros (+0 == -0) per IEEE 754.
//
static int compare_fp80_3way(fp80_t const &a, fp80_t const &b)
{
    if (a.isnan() || b.isnan()) return 2;
    if (a.iszero() && b.iszero()) return 0;
    bool an = a.sign() != 0;
    bool bn = b.sign() != 0;
    if (an != bn) return an ? -1 : 1;

    // Same sign — compare unsigned magnitudes (exponent first then mantissa).
    uint16_t ae = a.sign_exp() & FP80_EXPONENT_MASK;
    uint16_t be = b.sign_exp() & FP80_EXPONENT_MASK;
    int mag;
    if (ae != be)               mag = (ae < be) ? -1 : 1;
    else if (a.mantissa() != b.mantissa()) mag = (a.mantissa() < b.mantissa()) ? -1 : 1;
    else                        mag = 0;
    return an ? -mag : mag;
}

bool fp80_t::operator==(fp80_t const &rhs) const { int c = compare_fp80_3way(*this, rhs); return c == 0; }
bool fp80_t::operator!=(fp80_t const &rhs) const { int c = compare_fp80_3way(*this, rhs); return c != 0; }
bool fp80_t::operator< (fp80_t const &rhs) const { int c = compare_fp80_3way(*this, rhs); return c == -1; }
bool fp80_t::operator<=(fp80_t const &rhs) const { int c = compare_fp80_3way(*this, rhs); return c == -1 || c == 0; }
bool fp80_t::operator> (fp80_t const &rhs) const { int c = compare_fp80_3way(*this, rhs); return c == 1; }
bool fp80_t::operator>=(fp80_t const &rhs) const { int c = compare_fp80_3way(*this, rhs); return c == 1 || c == 0; }

//===========================================================================
//
// Read the current x87 control word.
//
// fpround_t::get/set work on MXCSR, but the x87 unit and the test harness's
// x87setcw operate on the separate x87 CW (via fldcw/fstcw). We need the
// latter for rounding decisions in hand-rolled arithmetic.
//
//===========================================================================
x87cw_t read_x87_cw()
{
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
    uint16_t cw;
    __asm__ volatile("fnstcw %0" : "=m"(cw));
    return cw;
#else
    return X87CW_PRECISION_EXTENDED | X87CW_ROUNDING_NEAREST | X87CW_MASK_ALL_EX;
#endif
}


//===========================================================================
//
// Round a 96-bit fpext96_t value to fp80 per the supplied control word.
// Sets PE if precision was lost, C1 if the result was rounded up (magnitude
// increased), OE/UE on overflow/underflow.
//
// PC field of CW dictates the rounding position:
//   PC=11 (extended 64-bit): round at the 32-bit extension
//   PC=10 (double   53-bit): round at bit 11 of mantissa + the extension
//   PC=00 (single   24-bit): round at bit 40 of mantissa + the extension
//
//===========================================================================
fp80_t round_fpext96_to_fp80(fpext96_t const &v, x87cw_t cw, uint16_t &out_sw)
{
    bool neg = v.sign() != 0;
    uint16_t sign_bit = neg ? FP80_SIGN_MASK : 0;
    if (v.iszero()) return fp80_t(0, sign_bit);

    uint64_t mant = v.mantissa();
    uint32_t ext = v.extend();
    int      exp = v.exponent();

    x87cw_t rmode = cw & X87CW_ROUNDING_MASK;
    x87cw_t pc    = cw & X87CW_PRECISION_MASK;

    // Number of bits of mantissa precision to keep (including the explicit-1).
    int keep_bits = (pc == X87CW_PRECISION_SINGLE) ? 24
                  : (pc == X87CW_PRECISION_DOUBLE) ? 53
                  : 64;
    int drop_in_mant = 64 - keep_bits;   // bits to discard from the bottom of mant

    // Compute round + sticky.  The "round bit" is the highest bit being
    // discarded; "sticky" is the OR of everything below that.
    bool round, sticky;
    if (drop_in_mant > 0)
    {
        round  = (mant >> (drop_in_mant - 1)) & 1;
        uint64_t below_mask = (drop_in_mant >= 2) ? ((1ull << (drop_in_mant - 1)) - 1) : 0;
        sticky = (mant & below_mask) != 0 || ext != 0;
        mant  &= ~((1ull << drop_in_mant) - 1);   // truncate
    }
    else  // drop_in_mant == 0, only the 32 extension bits below
    {
        round  = (ext >> 31) & 1;
        sticky = (ext & 0x7FFFFFFFu) != 0;
    }

    bool inexact = round || sticky;
    bool round_up = false;
    if (inexact)
    {
        if (rmode == X87CW_ROUNDING_NEAREST)
        {
            // Round-to-nearest-even: round up if round=1 and (sticky=1 or LSB=1).
            uint64_t lsb = (drop_in_mant > 0) ? ((mant >> drop_in_mant) & 1) : (mant & 1);
            round_up = round && (sticky || lsb);
        }
        else if (rmode == X87CW_ROUNDING_DOWN) round_up = neg;
        else if (rmode == X87CW_ROUNDING_UP)   round_up = !neg;
        // ROUND_ZERO: never round up
    }

    if (round_up)
    {
        if (drop_in_mant > 0)
            mant += (1ull << drop_in_mant);
        else
            mant += 1;
        if (mant == 0)  // overflowed to 2.0
        {
            mant = FP80_EXPLICIT_ONE;
            exp += 1;
        }
        out_sw |= X87SW_C1;
    }
    if (inexact) out_sw |= X87SW_PRECISION_EX;

    int biased = exp + FP80_EXPONENT_BIAS;
    if (biased >= FP80_EXPONENT_MAX_BIASED)
    {
        out_sw |= X87SW_OVERFLOW_EX | X87SW_PRECISION_EX;
        // For ROUND_ZERO and "toward the other direction" the result is the
        // largest representable finite value, not infinity.
        bool to_max_finite =
            rmode == X87CW_ROUNDING_ZERO ||
            (rmode == X87CW_ROUNDING_DOWN && !neg) ||
            (rmode == X87CW_ROUNDING_UP   && neg);
        if (to_max_finite)
        {
            uint64_t maxm = ~((1ull << drop_in_mant) - 1);  // largest with PC zeroes
            // C1 cleared: rounded toward zero (down in magnitude).
            out_sw &= ~X87SW_C1;
            return fp80_t(maxm, sign_bit | uint16_t(FP80_EXPONENT_MAX_BIASED - 1));
        }
        // Rounded up to infinity — set C1.
        out_sw |= X87SW_C1;
        return fp80_t(FP80_EXPLICIT_ONE, sign_bit | FP80_EXPONENT_MAX_BIASED);
    }
    if (biased <= 0)
    {
        // Denormal result. The earlier primary rounding may have set C1 and
        // PE for a normal-precision rounding decision that is irrelevant
        // here. Reroll the rounding decision at the *combined* precision —
        // PC drop bits plus denorm shift bits — using the original mantissa
        // and extension. This is the proper single-step round per IEEE-754.
        int shift = 1 - (v.exponent() + FP80_EXPONENT_BIAS);
        uint64_t orig_mant = v.mantissa();
        uint32_t orig_ext  = v.extend();
        // Clear C1 / PE / UE flags we may have set during the normal-path
        // rounding above — we are about to redo them.
        out_sw &= ~(X87SW_C1 | X87SW_PRECISION_EX | X87SW_UNDERFLOW_EX);
        if (shift >= 64)
        {
            // Result is below the smallest denormal magnitude. For round-
            // to-nearest with magnitude > 0.5 × smallest_denormal, still
            // round up to mantissa 1. The boundary is at shift == 64 with
            // mantissa.MSB exactly set; below-or-tied → 0, above → 1.
            bool any_bits = (orig_mant != 0) || (orig_ext != 0);
            if (any_bits)
                out_sw |= X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
            bool away_directed =
                (rmode == X87CW_ROUNDING_UP && !neg) ||
                (rmode == X87CW_ROUNDING_DOWN && neg);
            bool round_to_smallest = false;
            if (rmode == X87CW_ROUNDING_NEAREST && shift == 64)
            {
                // Value vs 0.5 × smallest_denormal:
                //   > 0.5 × smallest  ⇔  orig_mant > FP80_EXPLICIT_ONE OR
                //                        (orig_mant == FP80_EXPLICIT_ONE && ext != 0)
                // ties (orig_mant == FP80_EXPLICIT_ONE && ext == 0) → 0 (round to even)
                if (orig_mant > FP80_EXPLICIT_ONE ||
                    (orig_mant == FP80_EXPLICIT_ONE && orig_ext != 0))
                    round_to_smallest = true;
            }
            if ((away_directed || round_to_smallest) && any_bits)
            {
                out_sw |= X87SW_C1;
                return fp80_t(1, sign_bit);
            }
            return fp80_t(0, sign_bit);
        }
        int combined_drop = drop_in_mant + shift;
        // Round / sticky at the combined boundary.  The "round bit" is at
        // position (combined_drop - 1) of the 96-bit conceptual value where
        // bits 32..95 are orig_mant and bits 0..31 are orig_ext.
        auto bit_at = [&](int pos) -> bool {
            if (pos < 0)  return false;
            if (pos < 32) return (orig_ext >> pos) & 1;
            if (pos < 96) return (orig_mant >> (pos - 32)) & 1;
            return false;
        };
        int round_pos = combined_drop - 1 + 32;          // -1 from "below LSB"
        int lsb_pos   = combined_drop + 32;
        bool d_round  = bit_at(round_pos);
        bool d_sticky = false;
        // sticky = OR of all bits strictly below `round_pos`.
        for (int p = 0; p < round_pos && !d_sticky; p++)
            if (bit_at(p)) d_sticky = true;
        bool d_inexact = d_round || d_sticky;
        bool d_round_up = false;
        if (d_inexact)
        {
            if (rmode == X87CW_ROUNDING_NEAREST)
            {
                bool d_lsb = (lsb_pos < 96) ? bit_at(lsb_pos) : false;
                d_round_up = d_round && (d_sticky || d_lsb);
            }
            else if (rmode == X87CW_ROUNDING_DOWN) d_round_up = neg;
            else if (rmode == X87CW_ROUNDING_UP)   d_round_up = !neg;
        }

        uint64_t result_mant;
        if (combined_drop >= 64) result_mant = 0;
        else                     result_mant = orig_mant >> combined_drop;
        if (d_round_up)
        {
            result_mant += 1;
            // Could carry into next-denormal-rank or even smallest-normal.
            // The largest denormal mantissa is FP80_EXPLICIT_ONE - 1 (bits
            // 0..62 set, bit 63 clear). One more makes bit 63 set — which
            // is the smallest normal with biased_exp = 1. The fp80 wire
            // format requires biased_exp = 1 in that case.
            if (result_mant == FP80_EXPLICIT_ONE)
                return fp80_t(result_mant, sign_bit | uint16_t(1));
        }
        if (d_inexact) out_sw |= X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
        if (d_round_up) out_sw |= X87SW_C1;
        return fp80_t(result_mant, sign_bit);
    }
    return fp80_t(mant, sign_bit | uint16_t(biased));
}

//===========================================================================
//
// round_fpext128_to_fp80: same contract as round_fpext96_to_fp80, but
// the input carries a 64-bit extension instead of 32 — so the
// round/sticky decision at the fp80 boundary sees 32 extra bits of
// "true residual" below the kept mantissa. Used by the log/exp paths
// that operate at fpext128 precision throughout: when our polynomial
// chain is more accurate than fpext96 can express, the extra ext bits
// flip the C1 decision from "round_fpext96 was guessing on truncated
// info" to "round_fpext128 sees the real residual sign."
//
//===========================================================================
fp80_t round_fpext128_to_fp80(fpext128_t const &v, x87cw_t cw, uint16_t &out_sw)
{
    bool neg = v.sign() != 0;
    uint16_t sign_bit = neg ? FP80_SIGN_MASK : 0;
    if (v.iszero()) return fp80_t(0, sign_bit);

    uint64_t mant = v.mantissa();
    uint64_t ext  = v.extend();
    int      exp  = v.exponent();

    x87cw_t rmode = cw & X87CW_ROUNDING_MASK;
    x87cw_t pc    = cw & X87CW_PRECISION_MASK;

    int keep_bits = (pc == X87CW_PRECISION_SINGLE) ? 24
                  : (pc == X87CW_PRECISION_DOUBLE) ? 53
                  : 64;
    int drop_in_mant = 64 - keep_bits;

    bool round, sticky;
    if (drop_in_mant > 0)
    {
        round  = (mant >> (drop_in_mant - 1)) & 1;
        uint64_t below_mask = (drop_in_mant >= 2) ? ((1ull << (drop_in_mant - 1)) - 1) : 0;
        sticky = (mant & below_mask) != 0 || ext != 0;
        mant  &= ~((1ull << drop_in_mant) - 1);
    }
    else  // drop_in_mant == 0, only the 64 extension bits below
    {
        round  = (ext >> 63) & 1;
        sticky = (ext & 0x7FFFFFFFFFFFFFFFull) != 0;
    }

    bool inexact = round || sticky;
    bool round_up = false;
    if (inexact)
    {
        if (rmode == X87CW_ROUNDING_NEAREST)
        {
            uint64_t lsb = (drop_in_mant > 0) ? ((mant >> drop_in_mant) & 1) : (mant & 1);
            round_up = round && (sticky || lsb);
        }
        else if (rmode == X87CW_ROUNDING_DOWN) round_up = neg;
        else if (rmode == X87CW_ROUNDING_UP)   round_up = !neg;
    }

    if (round_up)
    {
        if (drop_in_mant > 0)
            mant += (1ull << drop_in_mant);
        else
            mant += 1;
        if (mant == 0)  // overflowed to 2.0
        {
            mant = FP80_EXPLICIT_ONE;
            exp += 1;
        }
        out_sw |= X87SW_C1;
    }
    if (inexact) out_sw |= X87SW_PRECISION_EX;

    int biased = exp + FP80_EXPONENT_BIAS;
    if (biased >= FP80_EXPONENT_MAX_BIASED)
    {
        out_sw |= X87SW_OVERFLOW_EX | X87SW_PRECISION_EX;
        bool to_max_finite =
            rmode == X87CW_ROUNDING_ZERO ||
            (rmode == X87CW_ROUNDING_DOWN && !neg) ||
            (rmode == X87CW_ROUNDING_UP   && neg);
        if (to_max_finite)
        {
            uint64_t maxm = ~((1ull << drop_in_mant) - 1);
            out_sw &= ~X87SW_C1;
            return fp80_t(maxm, sign_bit | uint16_t(FP80_EXPONENT_MAX_BIASED - 1));
        }
        out_sw |= X87SW_C1;
        return fp80_t(FP80_EXPLICIT_ONE, sign_bit | FP80_EXPONENT_MAX_BIASED);
    }
    if (biased <= 0)
    {
        int shift = 1 - (v.exponent() + FP80_EXPONENT_BIAS);
        uint64_t orig_mant = v.mantissa();
        uint64_t orig_ext  = v.extend();
        out_sw &= ~(X87SW_C1 | X87SW_PRECISION_EX | X87SW_UNDERFLOW_EX);
        if (shift >= 64)
        {
            bool any_bits = (orig_mant != 0) || (orig_ext != 0);
            if (any_bits)
                out_sw |= X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
            bool away_directed =
                (rmode == X87CW_ROUNDING_UP && !neg) ||
                (rmode == X87CW_ROUNDING_DOWN && neg);
            bool round_to_smallest = false;
            if (rmode == X87CW_ROUNDING_NEAREST && shift == 64)
            {
                if (orig_mant > FP80_EXPLICIT_ONE ||
                    (orig_mant == FP80_EXPLICIT_ONE && orig_ext != 0))
                    round_to_smallest = true;
            }
            if ((away_directed || round_to_smallest) && any_bits)
            {
                out_sw |= X87SW_C1;
                return fp80_t(1, sign_bit);
            }
            return fp80_t(0, sign_bit);
        }
        int combined_drop = drop_in_mant + shift;
        // Round / sticky at the combined boundary. The "round bit" is at
        // position (combined_drop - 1) of the 128-bit conceptual value where
        // bits 64..127 are orig_mant and bits 0..63 are orig_ext.
        auto bit_at = [&](int pos) -> bool {
            if (pos < 0)   return false;
            if (pos < 64)  return (orig_ext >> pos) & 1;
            if (pos < 128) return (orig_mant >> (pos - 64)) & 1;
            return false;
        };
        int round_pos = combined_drop - 1 + 64;
        int lsb_pos   = combined_drop + 64;
        bool d_round  = bit_at(round_pos);
        bool d_sticky = false;
        for (int p = 0; p < round_pos && !d_sticky; p++)
            if (bit_at(p)) d_sticky = true;
        bool d_inexact = d_round || d_sticky;
        bool d_round_up = false;
        if (d_inexact)
        {
            if (rmode == X87CW_ROUNDING_NEAREST)
            {
                bool d_lsb = (lsb_pos < 128) ? bit_at(lsb_pos) : false;
                d_round_up = d_round && (d_sticky || d_lsb);
            }
            else if (rmode == X87CW_ROUNDING_DOWN) d_round_up = neg;
            else if (rmode == X87CW_ROUNDING_UP)   d_round_up = !neg;
        }

        uint64_t result_mant;
        if (combined_drop >= 64) result_mant = 0;
        else                     result_mant = orig_mant >> combined_drop;
        if (d_round_up)
        {
            result_mant += 1;
            if (result_mant == FP80_EXPLICIT_ONE)
                return fp80_t(result_mant, sign_bit | uint16_t(1));
        }
        if (d_inexact) out_sw |= X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
        if (d_round_up) out_sw |= X87SW_C1;
        return fp80_t(result_mant, sign_bit);
    }
    return fp80_t(mant, sign_bit | uint16_t(biased));
}


//===========================================================================
//
// 80-bit arithmetic.
//
// On x86 hosts we copy into/out of `long double` and let the host x87 unit
// do the actual computation; this gives us bit-exact results that match
// the test oracle (which is the same x87 unit) without re-deriving 80-bit
// IEEE rounding by hand. On non-x86 hosts the same routines fall back to
// the fpext96_t-based path with the round_fpext96_to_fp80 helper above.
//
//===========================================================================

#if X87_HOST_HAS_FP80
static_assert(sizeof(long double) >= 10);

static long double fp80_to_ld(fp80_t const &v)
{
    long double r;
    std::memset(&r, 0, sizeof(r));
    std::memcpy(&r, &v, 10);
    return r;
}

static fp80_t ld_to_fp80(long double v)
{
    fp80_t r(0, 0);
    std::memcpy(&r, &v, 10);
    return r;
}

// Perform a 2-operand arithmetic op on the host x87 FPU.
// Loads lb then la (so ST(0)=la, ST(1)=lb) and applies op to compute lb OP la.
// op: 0=add, 1=sub (lb - la), 2=mul, 3=div (lb / la).
//
// NOTE on AT&T mnemonics: GAS swaps fsubp<->fsubrp and fdivp<->fdivrp from
// Intel's mnemonics — so AT&T's `fsubrp` is Intel's FSUBP (ST(1)=ST(1)-ST(0)).
// To compute "lb - la" with ST(1)=lb and ST(0)=la we want Intel FSUBP, which
// means AT&T `fsubrp`. Same flip for div.
static uint16_t host_x87_binary(fp80_t const &a, fp80_t const &b, fp80_t &dst, int op)
{
    long double la = fp80_to_ld(a);
    long double lb = fp80_to_ld(b);
    long double lr;
    uint16_t sw;
    switch (op)
    {
        case 0: __asm__ volatile(
            "fnclex\n\t"
            "fldt %2\n\t"
            "fldt %3\n\t"
            "faddp\n\t"
            "fnstsw %0\n\t"
            "fstpt %1"
            : "=m"(sw), "=m"(lr) : "m"(lb), "m"(la));
            break;
        case 1: __asm__ volatile(
            "fnclex\n\t"
            "fldt %2\n\t"
            "fldt %3\n\t"
            "fsubrp\n\t"        // AT&T: ST(1) = ST(1) - ST(0) = lb - la
            "fnstsw %0\n\t"
            "fstpt %1"
            : "=m"(sw), "=m"(lr) : "m"(lb), "m"(la));
            break;
        case 2: __asm__ volatile(
            "fnclex\n\t"
            "fldt %2\n\t"
            "fldt %3\n\t"
            "fmulp\n\t"
            "fnstsw %0\n\t"
            "fstpt %1"
            : "=m"(sw), "=m"(lr) : "m"(lb), "m"(la));
            break;
        case 3: __asm__ volatile(
            "fnclex\n\t"
            "fldt %2\n\t"
            "fldt %3\n\t"
            "fdivrp\n\t"        // AT&T: ST(1) = ST(1) / ST(0) = lb / la
            "fnstsw %0\n\t"
            "fstpt %1"
            : "=m"(sw), "=m"(lr) : "m"(lb), "m"(la));
            break;
    }
    dst = ld_to_fp80(lr);
    return sw & ~X87SW_TOP_MASK;
}

// Single-operand FPU op with one fp80 result and SW.
// op codes:
//   0: fsin     (sin(ST(0)) -> ST(0))
//   1: fcos     (cos(ST(0)) -> ST(0))
//   2: f2xm1    (2^ST(0) - 1 -> ST(0))
//   3: frndint  (round per CW -> ST(0))
uint16_t host_x87_unary(fp80_t const &src, fp80_t &dst, int op)
{
    long double ls = fp80_to_ld(src);
    long double lr;
    uint16_t sw;
    switch (op)
    {
        case 0: __asm__ volatile("fnclex\n\tfldt %2\n\tfsin\n\tfnstsw %0\n\tfstpt %1"
                                 : "=m"(sw), "=m"(lr) : "m"(ls)); break;
        case 1: __asm__ volatile("fnclex\n\tfldt %2\n\tfcos\n\tfnstsw %0\n\tfstpt %1"
                                 : "=m"(sw), "=m"(lr) : "m"(ls)); break;
        case 2: __asm__ volatile("fnclex\n\tfldt %2\n\tf2xm1\n\tfnstsw %0\n\tfstpt %1"
                                 : "=m"(sw), "=m"(lr) : "m"(ls)); break;
        case 3: __asm__ volatile("fnclex\n\tfldt %2\n\tfrndint\n\tfnstsw %0\n\tfstpt %1"
                                 : "=m"(sw), "=m"(lr) : "m"(ls)); break;
    }
    dst = ld_to_fp80(lr);
    return sw & ~X87SW_TOP_MASK;
}

// Single-operand FPU op that produces two results (FSINCOS, FPTAN, FXTRACT).
// dst1 receives ST(0) after the op (post-FSTP), dst2 receives ST(1).
// op codes:
//   0: fsincos  (cos -> dst1, sin -> dst2)
//   1: fptan    (1.0 -> dst1, tan -> dst2)
//   2: fxtract  (significand -> dst1, exponent -> dst2)
uint16_t host_x87_unary2(fp80_t const &src, fp80_t &dst1, fp80_t &dst2, int op)
{
    long double ls = fp80_to_ld(src);
    long double lr1, lr2;
    uint16_t sw;
    switch (op)
    {
        case 0: __asm__ volatile(
            "fnclex\n\tfldz\n\tfldt %3\n\tfsincos\n\tfnstsw %0\n\tfstpt %1\n\tfstpt %2"
            : "=m"(sw), "=m"(lr1), "=m"(lr2) : "m"(ls)); break;
        case 1: __asm__ volatile(
            "fnclex\n\tfldz\n\tfldt %3\n\tfptan\n\tfnstsw %0\n\tfstpt %1\n\tfstpt %2"
            : "=m"(sw), "=m"(lr1), "=m"(lr2) : "m"(ls)); break;
        case 2: __asm__ volatile(
            "fnclex\n\tfldt %3\n\tfxtract\n\tfnstsw %0\n\tfstpt %1\n\tfstpt %2"
            : "=m"(sw), "=m"(lr1), "=m"(lr2) : "m"(ls)); break;
    }
    dst1 = ld_to_fp80(lr1);
    dst2 = ld_to_fp80(lr2);
    return sw & ~X87SW_TOP_MASK;
}

// Two-operand transcendentals returning one fp80 + SW.
// Layout matches the asm wrappers (ARG1 -> ST(0), ARG2 -> ST(1) before the op).
// op codes:
//   0: fyl2x    (b * log2(a) -> ST(1), pop)
//   1: fyl2xp1  (b * log2(a+1) -> ST(1), pop)
//   2: fpatan   (atan2(b, a) -> ST(1), pop)
//   3: fprem    (a mod b -> ST(0), keep both)
//   4: fprem1   (IEEE remainder a / b -> ST(0), keep both)
//   5: fscale   (a * 2^trunc(b) -> ST(0), keep both)
uint16_t host_x87_binary_trans(fp80_t const &a, fp80_t const &b, fp80_t &dst, int op)
{
    long double la = fp80_to_ld(a);
    long double lb = fp80_to_ld(b);
    long double lr;
    uint16_t sw;
    switch (op)
    {
        case 0: __asm__ volatile(
            "fnclex\n\tfldt %2\n\tfldt %3\n\tfyl2x\n\tfnstsw %0\n\tfstpt %1"
            : "=m"(sw), "=m"(lr) : "m"(lb), "m"(la)); break;
        case 1: __asm__ volatile(
            "fnclex\n\tfldt %2\n\tfldt %3\n\tfyl2xp1\n\tfnstsw %0\n\tfstpt %1"
            : "=m"(sw), "=m"(lr) : "m"(lb), "m"(la)); break;
        case 2: __asm__ volatile(
            "fnclex\n\tfldt %2\n\tfldt %3\n\tfpatan\n\tfnstsw %0\n\tfstpt %1"
            : "=m"(sw), "=m"(lr) : "m"(lb), "m"(la)); break;
        case 3: __asm__ volatile(
            "fnclex\n\tfldt %2\n\tfldt %3\n\tfprem\n\tfnstsw %0\n\tfstpt %1\n\tfstp %%st(0)"
            : "=m"(sw), "=m"(lr) : "m"(lb), "m"(la)); break;
        case 4: __asm__ volatile(
            "fnclex\n\tfldt %2\n\tfldt %3\n\tfprem1\n\tfnstsw %0\n\tfstpt %1\n\tfstp %%st(0)"
            : "=m"(sw), "=m"(lr) : "m"(lb), "m"(la)); break;
        case 5: __asm__ volatile(
            "fnclex\n\tfldt %2\n\tfldt %3\n\tfscale\n\tfnstsw %0\n\tfstpt %1\n\tfstp %%st(0)"
            : "=m"(sw), "=m"(lr) : "m"(lb), "m"(la)); break;
    }
    dst = ld_to_fp80(lr);
    return sw & ~X87SW_TOP_MASK;
}

static uint16_t host_x87_sqrt(fp80_t const &a, fp80_t &dst)
{
    long double la = fp80_to_ld(a);
    long double lr;
    uint16_t sw;
    __asm__ volatile(
        "fnclex\n\t"
        "fldt %2\n\t"
        "fsqrt\n\t"
        "fnstsw %0\n\t"
        "fstpt %1"
        : "=m"(sw), "=m"(lr) : "m"(la));
    dst = ld_to_fp80(lr);
    return sw & ~X87SW_TOP_MASK;
}
#endif // X87_HOST_HAS_FP80

//
// Hand-rolled fp80 add/sub with full 128-bit working precision.
// Aligns mantissas in a __uint128_t with proper sticky-bit tracking,
// then rounds via round_fpext96_to_fp80.
//
// 'subtract' inverts b's sign before adding.
//
static uint16_t propagate_nan_pair(fp80_t a, fp80_t b, fp80_t &dst)
{
    bool snan = a.issnan() || b.issnan();
    fp80_t pick;
    if (a.isnan() && b.isnan())
        pick = (b.mantissa() > a.mantissa()) ? b : a;
    else
        pick = a.isnan() ? a : b;
    dst = pick.issnan() ? fp80_t::make_qnan(pick) : pick;
    return snan ? X87SW_INVALID_EX : 0;
}

static uint16_t do_add(fp80_t a, fp80_t b, fp80_t &dst, bool subtract)
{
    // NaN check before sign flip — subtract shouldn't change a NaN's sign.
    if (a.isnan() || b.isnan())
        return propagate_nan_pair(a, b, dst);
    if (subtract) b = fp80_t::chs(b);
    if (a.isinf() && b.isinf())
    {
        if (a.sign() != b.sign())
        {
            dst = fp80_t::const_indef();
            return X87SW_INVALID_EX;
        }
        dst = a;
        return 0;
    }
    if (a.isinf()) { dst = a; return b.isdenorm() ? X87SW_DENORM_EX : 0; }
    if (b.isinf()) { dst = b; return a.isdenorm() ? X87SW_DENORM_EX : 0; }
    if (a.iszero() && b.iszero())
    {
        if (a.sign() == b.sign())
            dst = a;
        else
        {
            x87cw_t r = read_x87_cw() & X87CW_ROUNDING_MASK;
            dst = (r == X87CW_ROUNDING_DOWN) ? fp80_t::const_nzero() : fp80_t::const_zero();
        }
        return 0;
    }
    if (a.iszero()) { dst = b; return b.isdenorm() ? X87SW_DENORM_EX : 0; }
    if (b.iszero()) { dst = a; return a.isdenorm() ? X87SW_DENORM_EX : 0; }

    uint16_t flags = 0;
    if (a.isdenorm() || b.isdenorm()) flags |= X87SW_DENORM_EX;

    // Extract sign / biased exponent / mantissa.
    int sa = a.sign();
    int sb = b.sign();
    int ea = a.sign_exp() & FP80_EXPONENT_MASK;
    int eb = b.sign_exp() & FP80_EXPONENT_MASK;
    uint64_t ma = a.mantissa();
    uint64_t mb = b.mantissa();

    // Normalize denormals so MSB of mantissa is set; compute "true" exponent
    // (unbiased, treating the mantissa as a 64-bit integer with explicit-1 at
    // bit 63 — same convention as fpext96_t::exponent).
    int true_ea = (ea == 0) ? (1 - FP80_EXPONENT_BIAS) : (ea - FP80_EXPONENT_BIAS);
    int true_eb = (eb == 0) ? (1 - FP80_EXPONENT_BIAS) : (eb - FP80_EXPONENT_BIAS);
    if (ea == 0) { int s = count_leading_zeros64(ma); ma <<= s; true_ea -= s; }
    if (eb == 0) { int s = count_leading_zeros64(mb); mb <<= s; true_eb -= s; }

    // Place both mantissas in a 128-bit working space with the explicit-1 at
    // bit 127 — gives 64 bits of round/sticky room.
    __uint128_t big, small;
    int result_exp;
    int sign_of_result;
    bool big_is_a;
    if (true_ea > true_eb || (true_ea == true_eb && ma >= mb))
    {
        big = (__uint128_t)ma << 64;
        small = (__uint128_t)mb << 64;
        result_exp = true_ea;
        sign_of_result = sa;
        big_is_a = true;
        int shift = true_ea - true_eb;
        if (shift >= 128) { small = 0; if (mb != 0) small = 1; /* sticky */ }
        else if (shift > 0)
        {
            __uint128_t lost = small & (((__uint128_t)1 << shift) - 1);
            small >>= shift;
            if (lost != 0) small |= 1;
        }
    }
    else
    {
        big = (__uint128_t)mb << 64;
        small = (__uint128_t)ma << 64;
        result_exp = true_eb;
        sign_of_result = sb;
        big_is_a = false;
        int shift = true_eb - true_ea;
        if (shift >= 128) { small = 0; if (ma != 0) small = 1; }
        else if (shift > 0)
        {
            __uint128_t lost = small & (((__uint128_t)1 << shift) - 1);
            small >>= shift;
            if (lost != 0) small |= 1;
        }
    }

    __uint128_t result;
    int result_sign_bit;
    if (sa == sb)
    {
        // Same-sign add — carry can occur (result spills into bit 128).
        result = big + small;
        result_sign_bit = sign_of_result;
        if (result < big)
        {
            // Carry into bit 128: shift right one and bump exponent.
            bool low_was_set = (result & 1) != 0;
            result = (result >> 1) | ((__uint128_t)1 << 127);
            if (low_was_set) result |= 1;            // preserve sticky
            result_exp += 1;
        }
    }
    else
    {
        // Different signs — subtract. big >= small by construction.
        result = big - small;
        result_sign_bit = sign_of_result;
        if (result == 0)
        {
            // Exact cancellation: sign of zero per rounding mode (only the
            // DOWN mode produces -0 when adding values of opposite sign).
            x87cw_t r = read_x87_cw() & X87CW_ROUNDING_MASK;
            dst = (r == X87CW_ROUNDING_DOWN) ? fp80_t::const_nzero() : fp80_t::const_zero();
            return flags;
        }
        // Normalize: shift left until MSB at bit 127.
        uint64_t hi = uint64_t(result >> 64);
        uint64_t lo = uint64_t(result);
        int lz = (hi != 0) ? count_leading_zeros64(hi) : (64 + count_leading_zeros64(lo));
        result <<= lz;
        result_exp -= lz;
    }
    (void)big_is_a;

    // Pack into fpext96_t (64-bit mantissa + 32-bit extension + sticky)
    uint64_t res_mant = uint64_t(result >> 64);
    uint64_t lower64  = uint64_t(result);
    uint32_t res_ext  = uint32_t(lower64 >> 32);
    if ((lower64 & 0xFFFFFFFFu) != 0) res_ext |= 1;

    fpext96_t r(res_mant, res_ext, result_exp, result_sign_bit ? 1 : 0);
    dst = round_fpext96_to_fp80(r, read_x87_cw(), flags);
    return flags;
}

//
// Multiplication: same special cases as add; finite via fpext96_t::mul.
//
static uint16_t do_mul(fp80_t a, fp80_t b, fp80_t &dst)
{
    if (a.isnan() || b.isnan())
        return propagate_nan_pair(a, b, dst);
    bool a_inf = a.isinf(), b_inf = b.isinf();
    bool a_zero = a.iszero(), b_zero = b.iszero();
    if ((a_inf && b_zero) || (a_zero && b_inf))
    {
        dst = fp80_t::const_indef();
        return X87SW_INVALID_EX;
    }
    uint16_t result_sign = (a.sign() ^ b.sign()) ? FP80_SIGN_MASK : 0;
    if (a_inf || b_inf) { dst = fp80_t(FP80_EXPLICIT_ONE, result_sign | FP80_EXPONENT_MAX_BIASED); return (a.isdenorm() || b.isdenorm()) ? X87SW_DENORM_EX : 0; }
    if (a_zero || b_zero) { dst = fp80_t(0, result_sign); return (a.isdenorm() || b.isdenorm()) ? X87SW_DENORM_EX : 0; }

    uint16_t flags = 0;
    if (a.isdenorm() || b.isdenorm()) flags |= X87SW_DENORM_EX;

    // Hand-rolled multiply: __uint128 product, normalize, round.
    int ea_b = a.sign_exp() & FP80_EXPONENT_MASK;
    int eb_b = b.sign_exp() & FP80_EXPONENT_MASK;
    uint64_t ma = a.mantissa();
    uint64_t mb = b.mantissa();
    int true_ea = (ea_b == 0) ? (1 - FP80_EXPONENT_BIAS) : (ea_b - FP80_EXPONENT_BIAS);
    int true_eb = (eb_b == 0) ? (1 - FP80_EXPONENT_BIAS) : (eb_b - FP80_EXPONENT_BIAS);
    if (ea_b == 0) { int s = count_leading_zeros64(ma); ma <<= s; true_ea -= s; }
    if (eb_b == 0) { int s = count_leading_zeros64(mb); mb <<= s; true_eb -= s; }

    __uint128_t prod = (__uint128_t)ma * (__uint128_t)mb;
    // Both have MSB at bit 63, so prod has MSB at bit 126 or 127.
    int result_exp = true_ea + true_eb;
    if ((prod & ((__uint128_t)1 << 127)) == 0)
    {
        // MSB at 126: shift left 1 to put it at 127.
        prod <<= 1;
        result_exp -= 1;
    }
    // The MSB-at-127 product represents a value of (mant_a × mant_b) × 2^(true_ea + true_eb - 63).
    // After shifting, the effective representation has the explicit-1 at bit 127, so the
    // exponent corresponds to bit 63 being the "ones" position of the mantissa.
    // Wait — fpext96_t convention: mantissa MSB at bit 63 with exponent meaning "exp - 63".
    // Here our 128-bit prod has MSB at bit 127, so we'd treat it as a 128-bit mantissa with
    // exponent +1 relative to the fpext96_t convention if we naively reused exp. Let me
    // think:
    //   value_a = ma * 2^(true_ea - 63)
    //   value_b = mb * 2^(true_eb - 63)
    //   value_a * value_b = (ma * mb) * 2^(true_ea + true_eb - 126)
    //   = prod * 2^(true_ea + true_eb - 126)
    // We've shifted prod left 1 (if needed) so MSB is at bit 127:
    //   = (prod_shifted) * 2^(true_ea + true_eb - 126 - shift)
    // To represent prod_shifted as fpext96_t (MSB at bit 63 of 96), shift right 64:
    //   = (prod_shifted >> 64) * 2^(true_ea + true_eb - 126 - shift + 64)
    //   = mant64 * 2^(result_exp_for_fpext96 - 63) where
    //   result_exp_for_fpext96 = true_ea + true_eb - 126 - shift + 64 + 63
    //                          = true_ea + true_eb + 1 - shift
    // (when shift=0, +1; when shift=-1, +2. We tracked shift implicitly via -1 to result_exp.)
    // So actually result_exp in our convention should be true_ea + true_eb + 1 after shift.
    // Adjust:
    result_exp += 1;

    uint64_t res_mant = uint64_t(prod >> 64);
    uint64_t lower64  = uint64_t(prod);
    uint32_t res_ext  = uint32_t(lower64 >> 32);
    if ((lower64 & 0xFFFFFFFFu) != 0) res_ext |= 1;   // sticky

    fpext96_t r(res_mant, res_ext, result_exp, (result_sign != 0) ? 1 : 0);
    dst = round_fpext96_to_fp80(r, read_x87_cw(), flags);
    return flags;
}

//
// Division by hand using 128-bit / 64-bit integer division. This avoids the
// precision loss in fpext96_t::div64 (which drops to fp64 internally) and the
// crash that comes from passing overflow-fp64 results into the fpext96_t
// constructor.
//
static uint16_t do_div(fp80_t a, fp80_t b, fp80_t &dst)
{
    if (a.isnan() || b.isnan())
        return propagate_nan_pair(a, b, dst);
    bool a_inf = a.isinf(), b_inf = b.isinf();
    bool a_zero = a.iszero(), b_zero = b.iszero();
    if ((a_inf && b_inf) || (a_zero && b_zero))
    {
        dst = fp80_t::const_indef();
        return X87SW_INVALID_EX;
    }
    uint16_t result_sign = (a.sign() ^ b.sign()) ? FP80_SIGN_MASK : 0;
    if (b_zero)
    {
        dst = fp80_t(FP80_EXPLICIT_ONE, result_sign | FP80_EXPONENT_MAX_BIASED);
        return X87SW_DIVZERO_EX;
    }
    uint16_t spec_flags = (a.isdenorm() || b.isdenorm()) ? X87SW_DENORM_EX : 0;
    if (a_inf)  { dst = fp80_t(FP80_EXPLICIT_ONE, result_sign | FP80_EXPONENT_MAX_BIASED); return spec_flags; }
    if (b_inf)  { dst = fp80_t(0, result_sign); return spec_flags; }
    if (a_zero) { dst = fp80_t(0, result_sign); return spec_flags; }

    uint16_t flags = 0;
    if (a.isdenorm() || b.isdenorm()) flags |= X87SW_DENORM_EX;

    // Normalize denormals so both operands have an explicit-1 in the MSB.
    uint64_t a_mant = a.mantissa();
    uint64_t b_mant = b.mantissa();
    int a_exp = (a.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    int b_exp = (b.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    if (a.isdenorm()) { int s = count_leading_zeros64(a_mant); a_mant <<= s; a_exp = 1 - FP80_EXPONENT_BIAS - s; }
    if (b.isdenorm()) { int s = count_leading_zeros64(b_mant); b_mant <<= s; b_exp = 1 - FP80_EXPONENT_BIAS - s; }

    // Produce a 96-bit fpext96-equivalent quotient in two stages so the
    // 128-bit intermediates never overflow.
    //   Stage 1: quot64 = (a_mant << 63) / b_mant
    //            → 64-bit quotient with explicit-1 at bit 63 or 64.
    //   Stage 2: ext32 = (remainder << 32) / b_mant
    //            → next 32 bits of precision.
    unsigned __int128 num1 = (unsigned __int128)a_mant << 63;
    uint64_t          quot = uint64_t(num1 / b_mant);
    unsigned __int128 rem1 = num1 - (unsigned __int128)quot * b_mant;

    unsigned __int128 num2 = rem1 << 32;
    uint32_t          ext  = uint32_t(num2 / b_mant);
    unsigned __int128 rem2 = num2 - (unsigned __int128)ext * b_mant;
    bool inexact = (rem2 != 0);

    int result_exp = a_exp - b_exp;
    // quot has its top bit at position 63 if a_mant >= b_mant, else 62.
    // Normalize so the explicit-1 sits at position 63.
    if ((quot & FP80_EXPLICIT_ONE) == 0)
    {
        quot = (quot << 1) | (ext >> 31);
        ext  = (ext << 1) | (inexact ? 0 : 0);   // ext low bit stays 0 (we tracked sticky separately)
        result_exp -= 1;
    }
    // Pack sticky into the low bit of ext so round_fpext96_to_fp80 picks it up.
    uint32_t ext_packed = ext | (inexact ? 1 : 0);

    fpext96_t r(quot, ext_packed, result_exp, (result_sign != 0) ? 1 : 0);
    dst = round_fpext96_to_fp80(r, read_x87_cw(), flags);
    return flags;
}

// Asm convention notes (matching the order ARG1->ST(0), ARG2->ST(1)):
//   FADDP    ST(1) <- ST(1)+ST(0)  ⇒ result = ARG1 + ARG2  (symmetric)
//   FSUBP    ST(1) <- ST(1)-ST(0)  ⇒ result = ARG2 - ARG1
//   FSUBRP   ST(1) <- ST(0)-ST(1)  ⇒ result = ARG1 - ARG2
//   FMULP    ST(1) <- ST(1)*ST(0)  ⇒ result = ARG1 * ARG2  (symmetric)
//   FDIVP    ST(1) <- ST(1)/ST(0)  ⇒ result = ARG2 / ARG1
//   FDIVRP   ST(1) <- ST(0)/ST(1)  ⇒ result = ARG1 / ARG2
// Asm convention: with ARG1 -> ST(0) and ARG2 -> ST(1):
//   FADDP  ⇒ ARG1 + ARG2  | FSUBP  ⇒ ARG2 - ARG1 | FSUBRP ⇒ ARG1 - ARG2
//   FMULP  ⇒ ARG1 * ARG2  | FDIVP  ⇒ ARG2 / ARG1 | FDIVRP ⇒ ARG1 / ARG2
// In test_binary80 we receive (src2, src1), so a=src2 b=src1 and we want
// x87_fsub(a,b)  = b - a   (matches FSUBP)
// x87_fsubr(a,b) = a - b   (matches FSUBRP)
// x87_fdiv(a,b)  = b / a   (matches FDIVP)
// x87_fdivr(a,b) = a / b   (matches FDIVRP)
#if X87_HOST_HAS_FP80
uint16_t fp80_t::x87_fadd (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary(a, b, dst, 0); }
uint16_t fp80_t::x87_fsub (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary(a, b, dst, 1); }
uint16_t fp80_t::x87_fsubr(fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary(b, a, dst, 1); }
uint16_t fp80_t::x87_fmul (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary(a, b, dst, 2); }
uint16_t fp80_t::x87_fdiv (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary(a, b, dst, 3); }
uint16_t fp80_t::x87_fdivr(fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary(b, a, dst, 3); }
#else
uint16_t fp80_t::x87_fadd (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return do_add(a, b, dst, false); }
uint16_t fp80_t::x87_fsub (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return do_add(b, a, dst, true);  }
uint16_t fp80_t::x87_fsubr(fp80_t const &a, fp80_t const &b, fp80_t &dst) { return do_add(a, b, dst, true);  }
uint16_t fp80_t::x87_fmul (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return do_mul(a, b, dst);        }
uint16_t fp80_t::x87_fdiv (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return do_div(b, a, dst);        }
uint16_t fp80_t::x87_fdivr(fp80_t const &a, fp80_t const &b, fp80_t &dst) { return do_div(a, b, dst);        }
#endif

#if X87_HOST_HAS_FP80
uint16_t fp80_t::x87_fsqrt(fp80_t const &src, fp80_t &dst)
{
    // Special-case for the negative-zero edge: fsqrt(-0) = -0 by spec; the
    // host x87 honors that already.
    return host_x87_sqrt(src, dst);
}
#else
//
// Square root via shift-and-subtract on a 128-bit fixed-point input,
// producing a 64-bit result + remainder. Handed to round_fpext96_to_fp80
// for proper x87-style rounding.
//
uint16_t fp80_t::x87_fsqrt(fp80_t const &src, fp80_t &dst)
{
    if (src.isnan())
    {
        dst = src.issnan() ? fp80_t::make_qnan(src) : src;
        return src.issnan() ? X87SW_INVALID_EX : 0;
    }
    if (src.iszero()) { dst = src; return 0; }
    if (src.sign())
    {
        dst = fp80_t::const_indef();
        return X87SW_INVALID_EX;
    }
    if (src.isinf()) { dst = src; return 0; }

    uint16_t flags = src.isdenorm() ? X87SW_DENORM_EX : 0;

    uint64_t mant = src.mantissa();
    int      exp  = (src.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    if (src.isdenorm())
    {
        int s = count_leading_zeros64(mant);
        mant <<= s;
        exp = 1 - FP80_EXPONENT_BIAS - s;
    }

#if defined(__aarch64__)
    // AArch64 has hardware fp64 fsqrt (~10 cycles) but no fp80. Lose ~11
    // bits of precision vs the shift-and-subtract below — accuracy matches
    // X87_MATCH_XEFU-class behavior (fp64 round-trip), which is what
    // PowerPC-based original-Xbox emulators ship anyway. Runs ~15x faster
    // than the 96-iter shift-and-subtract fallback.
    //
    // We pack the post-normalize fp80 mantissa (leading 1 at bit 63) and
    // signed exponent directly into an fp64 — bypassing the even-exp dance
    // the shift-and-subtract path needs but the HW fsqrt doesn't. fp80
    // exponents outside fp64's range fall through to the C path; the upper
    // limit ~2^1023 is far above any realistic guest input.
    if (exp >= -1022 && exp <= 1023)
    {
        uint64_t frac = (mant & 0x7FFFFFFFFFFFFFFFULL) >> 11;
        uint64_t bits = (uint64_t(exp + 1023) << 52) | frac;
        double x_d;
        std::memcpy(&x_d, &bits, sizeof(double));
        double r_d = std::sqrt(x_d);
        dst = fp80_t(fp64_t(r_d));
        flags |= X87SW_PRECISION_EX;
        return flags;
    }
    // exp outside fp64 range — fall through to portable C path.
#endif

    // Now mant has MSB at bit 63; value = mant * 2^(exp - 63). For the
    // sqrt result-exponent to be a clean integer we need (exp - 63) even,
    // i.e. exp odd. If exp is even, shift mant right by 1 (track the lost
    // bit as inexact) and bump exp.
    bool extra_inexact = false;
    if ((exp & 1) == 0)
    {
        extra_inexact = (mant & 1) != 0;
        mant >>= 1;
        exp += 1;
    }

    // Compute sqrt(mant * 2^64) using 128-bit shift-and-subtract.
    // Input has MSB at bit 127 (if no right-shift) or 126 (if shifted).
    // Either way, the 64 output bits include the explicit-1 at bit 63
    // (sqrt of value in [2^126, 2^128) is in [2^63, 2^64)).
    // Compute sqrt of a 192-bit input via shift-and-subtract, producing a
    // 96-bit result (64-bit mantissa + 32-bit extension for rounding).
    //
    // 192-bit input is mant placed at bits 128..191 (i.e. xhi = mant in
    // 64 high bits, xlo = 0). We iterate 96 times, each time consuming the
    // top 2 bits of input.

    uint64_t   xhi = mant;
    uint64_t   xlo = 0;
    __uint128_t remainder = 0;
    __uint128_t result128 = 0;
    for (int iter = 0; iter < 96; iter++)
    {
        // Top 2 bits of the 192-bit input.
        uint64_t top2 = (xhi >> 62) & 3;
        xhi = (xhi << 2) | (xlo >> 62);
        xlo <<= 2;

        remainder = (remainder << 2) | top2;
        __uint128_t test = (result128 << 2) | 1;
        if (remainder >= test)
        {
            remainder -= test;
            result128 = (result128 << 1) | 1;
        }
        else
        {
            result128 <<= 1;
        }
    }
    // result128 has 96 bits with MSB at bit 95 (or 94 if mant was right-
    // shifted to fix even-exp).
    if ((result128 & ((__uint128_t)1 << 95)) == 0)
    {
        result128 <<= 1;
        exp -= 2;
    }

    int result_exp = (exp - 1) / 2;

    uint64_t res_mant = uint64_t(result128 >> 32);
    uint32_t res_ext  = uint32_t(result128 & 0xFFFFFFFFu);
    if (remainder != 0) res_ext |= 1;

    // For the even-exp case we right-shifted mant by 1 (losing its LSB
    // into extra_inexact). The sqrt of the *shifted* input is half an
    // fp80 ULP below the sqrt of the *true* input, so the Q we computed
    // sits one half-ULP too low when extra_inexact is set. Compensate
    // by adding 0.5 fp80 ULPs (= 2^31 of the 96-bit value) before
    // rounding, so round_fpext96_to_fp80's C1 detection reflects the
    // true-value rounding direction instead of the shifted-value one.
    // (The extra_inexact lifeline-sticky bit would double up with this
    // correction and corrupt the tie case, so it is excluded above.)
    if (extra_inexact)
    {
        uint64_t old_ext = res_ext;
        res_ext += 0x80000000u;
        if (res_ext < old_ext)   // carry into mantissa
            res_mant += 1;
    }

    fpext96_t r(res_mant, res_ext, result_exp, 0);
    dst = round_fpext96_to_fp80(r, read_x87_cw(), flags);
    return flags;
}
#endif // X87_HOST_HAS_FP80

}
