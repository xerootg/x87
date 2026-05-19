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

#include "x87fp80.h"
#include "x87fp64.h"

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

fp80_t operator+(fp80_t const &, fp80_t const &) { return fp80_t::const_indef(); }
fp80_t operator-(fp80_t const &, fp80_t const &) { return fp80_t::const_indef(); }
fp80_t operator*(fp80_t const &, fp80_t const &) { return fp80_t::const_indef(); }
fp80_t operator/(fp80_t const &, fp80_t const &) { return fp80_t::const_indef(); }

fp80_t fp80_t::sqrt(fp80_t const &)  { return fp80_t::const_indef(); }
fp80_t fp80_t::floor(fp80_t const &) { return fp80_t::const_indef(); }
fp80_t fp80_t::ceil(fp80_t const &)  { return fp80_t::const_indef(); }

fp80_t fp80_t::ldexp(fp80_t const &, int32_t) { return fp80_t::const_indef(); }

bool fp80_t::operator< (fp80_t const &) const { return false; }
bool fp80_t::operator<=(fp80_t const &) const { return false; }
bool fp80_t::operator> (fp80_t const &) const { return false; }
bool fp80_t::operator>=(fp80_t const &) const { return false; }

uint16_t fp80_t::x87_fadd (fp80_t const &, fp80_t const &, fp80_t &dst) { dst = fp80_t::const_indef(); return X87SW_INVALID_EX; }
uint16_t fp80_t::x87_fsub (fp80_t const &, fp80_t const &, fp80_t &dst) { dst = fp80_t::const_indef(); return X87SW_INVALID_EX; }
uint16_t fp80_t::x87_fsubr(fp80_t const &, fp80_t const &, fp80_t &dst) { dst = fp80_t::const_indef(); return X87SW_INVALID_EX; }
uint16_t fp80_t::x87_fmul (fp80_t const &, fp80_t const &, fp80_t &dst) { dst = fp80_t::const_indef(); return X87SW_INVALID_EX; }
uint16_t fp80_t::x87_fdiv (fp80_t const &, fp80_t const &, fp80_t &dst) { dst = fp80_t::const_indef(); return X87SW_INVALID_EX; }
uint16_t fp80_t::x87_fdivr(fp80_t const &, fp80_t const &, fp80_t &dst) { dst = fp80_t::const_indef(); return X87SW_INVALID_EX; }
uint16_t fp80_t::x87_fsqrt(fp80_t const &, fp80_t &dst)                 { dst = fp80_t::const_indef(); return X87SW_INVALID_EX; }

}
