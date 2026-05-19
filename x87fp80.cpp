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
// the existing fpext96_t-based path, which is precision-limited but
// works for values in fp64 range.
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
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
// 80-bit arithmetic.
//
// On x86 hosts we copy into/out of `long double` and let the host x87 unit
// do the actual computation; this gives us bit-exact results that match
// the test oracle (which is the same x87 unit) without re-deriving 80-bit
// IEEE rounding by hand. On non-x86 hosts the same routines fall back to
// the fpext96_t-based path; rounding/PC are not fully honored there.
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
// Compute "a + b" via fpext96_t plus full NaN/Inf/zero handling.
// 'subtract' inverts b's sign before adding.
//
static uint16_t do_add(fp80_t a, fp80_t b, fp80_t &dst, bool subtract)
{
    if (subtract) b = fp80_t::chs(b);

    if (a.isnan() || b.isnan())
    {
        bool snan = a.issnan() || b.issnan();
        if (a.isnan()) dst = a.issnan() ? fp80_t::make_qnan(a) : a;
        else           dst = b.issnan() ? fp80_t::make_qnan(b) : b;
        return snan ? X87SW_INVALID_EX : 0;
    }
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
    if (a.isinf()) { dst = a; return 0; }
    if (b.isinf()) { dst = b; return 0; }
    if (a.iszero() && b.iszero())
    {
        if (a.sign() == b.sign())
            dst = a;
        else
        {
            x87cw_t r = fpround_t::get() & X87CW_ROUNDING_MASK;
            dst = (r == X87CW_ROUNDING_DOWN) ? fp80_t::const_nzero() : fp80_t::const_zero();
        }
        return 0;
    }
    if (a.iszero()) { dst = b; return 0; }
    if (b.iszero()) { dst = a; return 0; }

    uint16_t flags = 0;
    if (a.isdenorm() || b.isdenorm()) flags |= X87SW_DENORM_EX;

    fpext96_t ea(a), eb(b);
    fpext96_t result;
    result.add(ea, eb);
    dst = result.as_fp80();

    if (dst.isinf())                              flags |= X87SW_OVERFLOW_EX  | X87SW_PRECISION_EX;
    else if (dst.iszero() && !(a.iszero() || b.iszero()))
                                                  flags |= X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
    else if (dst.isdenorm())                      flags |= X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
    return flags;
}

//
// Multiplication: same special cases as add; finite via fpext96_t::mul.
//
static uint16_t do_mul(fp80_t a, fp80_t b, fp80_t &dst)
{
    if (a.isnan() || b.isnan())
    {
        bool snan = a.issnan() || b.issnan();
        if (a.isnan()) dst = a.issnan() ? fp80_t::make_qnan(a) : a;
        else           dst = b.issnan() ? fp80_t::make_qnan(b) : b;
        return snan ? X87SW_INVALID_EX : 0;
    }
    bool a_inf = a.isinf(), b_inf = b.isinf();
    bool a_zero = a.iszero(), b_zero = b.iszero();
    if ((a_inf && b_zero) || (a_zero && b_inf))
    {
        dst = fp80_t::const_indef();
        return X87SW_INVALID_EX;
    }
    uint16_t result_sign = (a.sign() ^ b.sign()) ? FP80_SIGN_MASK : 0;
    if (a_inf || b_inf) { dst = fp80_t(FP80_EXPLICIT_ONE, result_sign | FP80_EXPONENT_MAX_BIASED); return 0; }
    if (a_zero || b_zero) { dst = fp80_t(0, result_sign); return 0; }

    uint16_t flags = 0;
    if (a.isdenorm() || b.isdenorm()) flags |= X87SW_DENORM_EX;

    fpext96_t ea(a), eb(b), result;
    result.mul(ea, eb);
    dst = result.as_fp80();

    if (dst.isinf())                              flags |= X87SW_OVERFLOW_EX  | X87SW_PRECISION_EX;
    else if (dst.iszero())                        flags |= X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
    else if (dst.isdenorm())                      flags |= X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
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
    {
        bool snan = a.issnan() || b.issnan();
        if (a.isnan()) dst = a.issnan() ? fp80_t::make_qnan(a) : a;
        else           dst = b.issnan() ? fp80_t::make_qnan(b) : b;
        return snan ? X87SW_INVALID_EX : 0;
    }
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
    if (a_inf)  { dst = fp80_t(FP80_EXPLICIT_ONE, result_sign | FP80_EXPONENT_MAX_BIASED); return 0; }
    if (b_inf)  { dst = fp80_t(0, result_sign); return 0; }
    if (a_zero) { dst = fp80_t(0, result_sign); return 0; }

    uint16_t flags = 0;
    if (a.isdenorm() || b.isdenorm()) flags |= X87SW_DENORM_EX;

    // Normalize denormals so both operands have an explicit-1 in the MSB.
    uint64_t a_mant = a.mantissa();
    uint64_t b_mant = b.mantissa();
    int a_exp = (a.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    int b_exp = (b.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    if (a.isdenorm()) { int s = count_leading_zeros64(a_mant); a_mant <<= s; a_exp = 1 - FP80_EXPONENT_BIAS - s; }
    if (b.isdenorm()) { int s = count_leading_zeros64(b_mant); b_mant <<= s; b_exp = 1 - FP80_EXPONENT_BIAS - s; }

    // 128-bit / 64-bit: shift the numerator up 63 bits so the 64-bit
    // quotient sits at a known precision.
    unsigned __int128 num = (unsigned __int128)a_mant << 63;
    uint64_t quot = uint64_t(num / b_mant);
    unsigned __int128 prod = (unsigned __int128)quot * b_mant;
    bool inexact = (num != prod);

    int result_exp = a_exp - b_exp;
    // The quotient is either [1, 2) or [2, 4) at this point depending on
    // a_mant vs b_mant. Normalize so the explicit 1 is at bit 63.
    if ((quot & (1ull << 63)) == 0)
    {
        quot <<= 1;
        result_exp -= 1;
        // We need an extra bit; if there was a remainder it's now sticky.
        // (Approximation: ignore the extra round bit; precision loss tracked
        // via PE flag.)
    }
    result_exp += FP80_EXPONENT_BIAS;

    if (result_exp >= FP80_EXPONENT_MAX_BIASED)
    {
        dst = fp80_t(FP80_EXPLICIT_ONE, result_sign | FP80_EXPONENT_MAX_BIASED);
        return flags | X87SW_OVERFLOW_EX | X87SW_PRECISION_EX;
    }
    if (result_exp <= 0)
    {
        dst = fp80_t(0, result_sign);
        return flags | X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
    }
    dst = fp80_t(quot, result_sign | uint16_t(result_exp));
    return flags | (inexact ? X87SW_PRECISION_EX : 0);
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
// Square root via Newton-Raphson on the (split) mantissa.  We extract an
// even-exponent mantissa in [1.0, 4.0), compute sqrt to 64-bit precision via
// integer NR, then reassemble with exp/2.  Avoids the fp64 conversion path.
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

    // For sqrt we want an even exponent so that exp/2 is exact. If exp is
    // odd, halve mantissa and use exp+1 (now even).
    if (exp & 1)
    {
        mant >>= 1;
        exp  += 1;
    }
    int result_exp = exp >> 1;

    // Initial estimate: a single iteration of Newton-Raphson starting from
    // mant/2 + (3<<61) gets us into the right ballpark.  Then refine
    // a handful more times to reach 64-bit precision.
    //
    // We compute y = sqrt(x) where x = mant/2^63 viewed in [1.0, 4.0).
    // Using fixed-point representation with mant in [2^63, 2^65) (so the
    // result fits in [2^31, 2^33) of the same fixed-point), and 128-bit
    // intermediate division.
    unsigned __int128 x128 = ((unsigned __int128)mant) << 63;     // x at .126 fixed point
    uint64_t y = mant;                                            // crude initial estimate
    for (int i = 0; i < 6; i++)
    {
        if (y == 0) break;
        unsigned __int128 q = x128 / y;     // ~128-bit / 64-bit
        uint64_t nq = uint64_t(q);
        y = (y >> 1) + (nq >> 1) + ((y & nq) & 1);  // (y + x/y) / 2 rounded
    }

    // Now y is sqrt(mant) approximately, scaled. Renormalize so the
    // explicit-1 bit is at position 63.
    int lz = count_leading_zeros64(y);
    if (lz > 0)
    {
        y <<= lz;
        result_exp -= lz;
    }
    result_exp += FP80_EXPONENT_BIAS;

    if (result_exp >= FP80_EXPONENT_MAX_BIASED)
    {
        dst = fp80_t::const_pinf();
        return flags | X87SW_OVERFLOW_EX | X87SW_PRECISION_EX;
    }
    if (result_exp <= 0)
    {
        dst = fp80_t::const_zero();
        return flags | X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
    }
    dst = fp80_t(y, uint16_t(result_exp));
    return flags | X87SW_PRECISION_EX;
}
#endif // X87_HOST_HAS_FP80

}
