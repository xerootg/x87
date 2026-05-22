//=========================================================
//  x87fpext.h
//
//  Helpers for extended-precision temporary calculations.
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

#pragma once

#ifndef X87FPEXT_H
#define X87FPEXT_H

#include "x87fp64.h"
#include "x87fp80.h"


namespace x87
{

template<typename ExtendedType> class fpextxx_t;

//===========================================================================
//
// fpext52_t
//
// This class implements a non-exploded 64-bit floating point value in a way
// that is compatible with the fpextxx_t class below.
//
//===========================================================================

class fpext52_t : public fp64_t
{
public:
    //
    // public constants
    //
    static constexpr int32_t EXPONENT_MIN = 0 - FP64_EXPONENT_BIAS;
    static constexpr uint64_t EXPLICIT_ONE = 0x8000000000000000ull;

    //
    // constructor for the common extended form
    //
    explicit fpext52_t() { }
    explicit fpext52_t(uint64_t high, uint32_t low, int32_t exponent, uint16_t sign);

    //
    // converting constructors
    //
    template<typename SrcExtendedType> explicit fpext52_t(fpextxx_t<SrcExtendedType> const &src);
    explicit fpext52_t(fp64_t const &src) : fp64_t(src) { }
    explicit fpext52_t(fp80_t const &src) : fp64_t(src) { }
    explicit fpext52_t(double src) : fp64_t(src) { }

    //
    // raw parts
    //
    constexpr bool extended() const { return false; }
    uint16_t sign() const { return this->fp64_t::sign(); }
    int32_t exponent() const { return this->fp64_t::exponent(); }
    uint64_t mantissa() const { return (this->fp64_t::mantissa() << 11) | EXPLICIT_ONE; }
    uint8_t extend() const { return 0; }

    //
    // raw setters
    //
    void set_sign(uint8_t sign) { m_value.i = (m_value.i & ~(1ull << FP64_SIGN_SHIFT)) | (uint64_t(sign) << FP64_SIGN_SHIFT); }
    void set_exponent(int32_t exp) { m_value.i = (m_value.i & ~FP64_EXPONENT_MASK) | (uint64_t(exp + FP64_EXPONENT_BIAS) << FP64_EXPONENT_SHIFT); }

    //
    // conversions
    //
    fp64_t as_fp64() const { return *this; }
    double as_double() const { return m_value.d; }
    fp80_t as_fp80() const { return this->fp64_t::as_fp80(); }

    //
    // queries
    //
    bool iszero() const { return this->fp64_t::iszero(); }

    //
    // unary self operations
    //
    fpext52_t &abs() { m_value.d = std::abs(m_value.d); return *this; }
    fpext52_t &chs() { m_value.d = -m_value.d; return *this; }

    //
    // operators
    //
    fpext52_t &operator+=(fpext52_t const &rhs) { static_cast<fp64_t &>(*this) += static_cast<fp64_t const &>(rhs); return *this; }
    fpext52_t &operator-=(fpext52_t const &rhs) { static_cast<fp64_t &>(*this) -= static_cast<fp64_t const &>(rhs); return *this; }
    fpext52_t &operator*=(fpext52_t const &rhs) { static_cast<fp64_t &>(*this) *= static_cast<fp64_t const &>(rhs); return *this; }

    //
    // comparison operators
    //
    bool operator==(fpext52_t const &rhs) const { return static_cast<fp64_t const &>(*this) == static_cast<fp64_t const &>(rhs); }
    bool operator!=(fpext52_t const &rhs) const { return static_cast<fp64_t const &>(*this) != static_cast<fp64_t const &>(rhs); }
    bool operator>(fpext52_t const &rhs) const { return static_cast<fp64_t const &>(*this) > static_cast<fp64_t const &>(rhs); }
    bool operator>=(fpext52_t const &rhs) const { return static_cast<fp64_t const &>(*this) >= static_cast<fp64_t const &>(rhs); }
    bool operator<(fpext52_t const &rhs) const { return static_cast<fp64_t const &>(*this) < static_cast<fp64_t const &>(rhs); }
    bool operator<=(fpext52_t const &rhs) const { return static_cast<fp64_t const &>(*this) <= static_cast<fp64_t const &>(rhs); }

    //
    // friends
    //
    friend fpext52_t operator+(fpext52_t const &a, fpext52_t const &b) { return fpext52_t(static_cast<fp64_t const &>(a) + static_cast<fp64_t const &>(b)); }
    friend fpext52_t operator-(fpext52_t const &a, fpext52_t const &b) { return fpext52_t(static_cast<fp64_t const &>(a) - static_cast<fp64_t const &>(b)); }
    friend fpext52_t operator*(fpext52_t const &a, fpext52_t const &b) { return fpext52_t(static_cast<fp64_t const &>(a) * static_cast<fp64_t const &>(b)); }

    //
    // core operations
    //
    fpext52_t div64(fpext52_t const &b) const { return fpext52_t(static_cast<fp64_t const &>(*this) / static_cast<fp64_t const &>(b)); }

    //
    // static helpers
    //
    static fpext52_t ldexp(fpext52_t const &a, int32_t dexp) { return fpext52_t(fp64_t::ldexp(a, dexp)); }
    static fpext52_t floor(fpext52_t const &a) { return fpext52_t(fp64_t::floor(a)); }
    static fpext52_t floor_abs_loint(fpext52_t const &a, uint64_t &intbits);

    //
    // constant values
    //
    static fpext52_t const zero;
    static fpext52_t const nzero;
    static fpext52_t const one;
    static fpext52_t const none;
    static fpext52_t const l2t;
    static fpext52_t const l2e;
    static fpext52_t const pi;
    static fpext52_t const pio2;
    static fpext52_t const pio4;
    static fpext52_t const lg2;
    static fpext52_t const ln2;
};



//===========================================================================
//
// fpextxx_t
//
// "Exploded" form of a floating-point value, with a 64-bit or 80-bit
// mantissa, 32-bit exponent, and separated sign. These classes are
// used for intermediate values during floating-point calculations and have
// the following limitations:
//
//   * Denormals are not supported; instead we have an enormous exponent
//      range. Incoming denormals are converted to this range, and when
//      collapsing to fp64_t or fp80_t, we create denormals as needed.
//
//   * NaNs/infinities are not supported: these should be filtered out
//      ahead of time. Infinities can be produced when collapsing huge
//      values to fp64_t or fp80_t.
//
//   * Full-precision divides are not supported.
//
// Multiple precisions are supported, based on the mantissa size:
//
//   fpext64_t: 64-bit mantissa; fastest by far, so use if possible
//
//   fpext96_t: 64-bit mantissa plus 32 bits of extension; add/sub are
//      a bit more expensive, multiply is a lot more expensive
//
//===========================================================================

template<typename ExtendedType>
class fpextxx_t
{
    //
    // ExtendedType can be uint8_t (no extension), uint16_t, uint32_t, or
    // uint64_t (the last giving 128-bit total precision: fpext128_t).
    //
    static_assert(sizeof(ExtendedType) <= 8);

    //
    // bit of a kludge, but treat ExtendedType of less than 16 bits as none
    //
    static constexpr bool EXTENDED = (sizeof(ExtendedType) >= 2);

    //
    // types of our components
    //
    using sign_t = uint16_t;
    using mantissa_t = uint64_t;
    using extend_t = ExtendedType;
    using exponent_t = int32_t;

public:
    //
    // public constants
    //
    static constexpr int EXTEND_BITS = 8 * sizeof(extend_t);
    static constexpr int MANTISSA_BITS = 8 * sizeof(mantissa_t) + (EXTENDED ? EXTEND_BITS : 0);
    static constexpr exponent_t EXPONENT_MIN = -10000000;
    static constexpr mantissa_t EXPLICIT_ONE = 0x8000000000000000ull;

    //
    // default constructor
    //
    explicit fpextxx_t()
    {
    }

    //
    // constexpr constructor from raw parts for constants
    //
    constexpr explicit fpextxx_t(mantissa_t high, uint32_t low, exponent_t exp, sign_t sign) :
        m_mantissa(high + (EXTENDED ? 0 : (low >> 31))),
        m_extend(EXTEND_BITS == 64 ? (extend_t(low) << 32) :
                 EXTEND_BITS == 32 ? extend_t(low) :
                 EXTEND_BITS == 16 ? extend_t((low >> 16) + ((low >> 15) & 1)) : 0),
        m_exponent(exp),
        m_sign(sign)
    {
    }

    // Tag for 64-bit-low constructor.
    struct low64_tag {};

    // For fpext128_t constants where the bottom 64 bits matter, accept a
    // 64-bit low via the tagged form. (Existing callers using uint32_t low
    // are unaffected.)
    constexpr explicit fpextxx_t(mantissa_t high, uint64_t low64, exponent_t exp, sign_t sign, low64_tag) :
        m_mantissa(high + (EXTENDED ? 0 : (low64 >> 63))),
        m_extend(EXTEND_BITS == 64 ? extend_t(low64) :
                 EXTEND_BITS == 32 ? extend_t(low64 >> 32) :
                 EXTEND_BITS == 16 ? extend_t(low64 >> 48) : 0),
        m_exponent(exp),
        m_sign(sign)
    {
    }

    //
    // converting constructors
    //
    template<typename SrcExtendedType> explicit fpextxx_t(fpextxx_t<SrcExtendedType> const &src, bool round = false);
    explicit fpextxx_t(fp64_t const &src);
    explicit fpextxx_t(fp80_t const &src);
    explicit fpextxx_t(double src) : fpextxx_t(fp64_t(src)) { }

    //
    // raw parts
    //
    constexpr bool extended() const { return EXTENDED; }
    sign_t sign() const { return m_sign; }
    exponent_t exponent() const { return m_exponent; }
    mantissa_t mantissa() const { return m_mantissa; }
    extend_t extend() const { return EXTENDED ? m_extend : 0; }

    //
    // raw setters
    //
    void set_sign(sign_t sign) { m_sign = sign; }
    void set_exponent(exponent_t exp) { m_exponent = exp; }

    //
    // conversions
    //
    fp64_t as_fp64() const;
    double as_double() const { return this->as_fp64().as_double(); }
    fp80_t as_fp80() const;

    //
    // queries
    //
    bool iszero() const { return (m_mantissa == 0 && (!EXTENDED || m_extend == 0)); }

    //
    // unary self operations
    //
    fpextxx_t &abs() { m_sign = 0; return *this; }
    fpextxx_t &chs() { m_sign ^= 1; return *this; }

    //
    // operators
    //
    fpextxx_t &operator+=(fpextxx_t const &rhs) { this->add(*this, rhs); return *this; }
    fpextxx_t &operator-=(fpextxx_t const &rhs) { this->sub(*this, rhs); return *this; }
    fpextxx_t &operator*=(fpextxx_t const &rhs) { this->mul(*this, rhs); return *this; }

    //
    // comparison operators
    //
    bool operator==(fpextxx_t const &rhs) const;
    bool operator!=(fpextxx_t const &rhs) const;
    bool operator>(fpextxx_t const &rhs) const;
    bool operator>=(fpextxx_t const &rhs) const;
    bool operator<(fpextxx_t const &rhs) const;
    bool operator<=(fpextxx_t const &rhs) const;

    //
    // friends
    //
    friend fpextxx_t operator+(fpextxx_t const &a, fpextxx_t const &b) { fpextxx_t res; res.add(a, b); return res; }
    friend fpextxx_t operator-(fpextxx_t const &a, fpextxx_t const &b) { fpextxx_t res; res.sub(a, b); return res; }
    friend fpextxx_t operator*(fpextxx_t const &a, fpextxx_t const &b) { fpextxx_t res; res.mul(a, b); return res; }

    //
    // core operations
    //
    void add(fpextxx_t const &a, fpextxx_t const &b);
    void sub(fpextxx_t const &a, fpextxx_t const &b);
    void mul(fpextxx_t const &a, fpextxx_t const &b);
    fpextxx_t div64(fpextxx_t const &b) const { return fpextxx_t(this->as_fp64() / b.as_fp64()); }

    //
    // static helpers
    //
    static fpextxx_t ldexp(fpextxx_t const &a, int32_t dexp) { fpextxx_t res = a; res.m_exponent += dexp; return res; }
    static fpextxx_t floor(fpextxx_t const &a);
    static fpextxx_t floor_abs_loint(fpextxx_t const &a, uint64_t &intbits);

    //
    // constant values
    //
    static fpextxx_t const zero;
    static fpextxx_t const nzero;
    static fpextxx_t const one;
    static fpextxx_t const none;
    static fpextxx_t const l2t;
    static fpextxx_t const l2e;
    static fpextxx_t const pi;
    static fpextxx_t const pio2;
    static fpextxx_t const pio4;
    static fpextxx_t const lg2;
    static fpextxx_t const ln2;

private:
    //
    // internal primitives
    //
    bool mantissa_eq(fpextxx_t const &a) const;
    bool mantissa_gt(fpextxx_t const &a) const;
    bool mantissa_lt(fpextxx_t const &a) const;
    void round_mantissa_up();
    void round_extend_up();
    void shift_mantissa_right(int count);
    void normalize();
    void add_values(fpextxx_t const &src1, fpextxx_t const &src2, int src2shift);
    void sub_values(fpextxx_t const &src1, fpextxx_t const &src2, int src2shift);

    //
    // internal state
    //
    mantissa_t m_mantissa;
    extend_t m_extend;
    sign_t m_sign;
    exponent_t m_exponent;
};

//
// derived types
//
using fpext64_t = fpextxx_t<uint8_t>;
using fpext96_t = fpextxx_t<uint32_t>;
using fpext128_t = fpextxx_t<uint64_t>;



//
// constants — `inline` (C++17) gives these vague linkage so the linker
// collapses duplicate definitions across translation units that include
// this header.
//
inline fpext52_t const fpext52_t::zero (0x0000000000000000ull, 0x00000000, fpext52_t::EXPONENT_MIN, 0);
inline fpext52_t const fpext52_t::nzero(0x0000000000000000ull, 0x00000000, fpext52_t::EXPONENT_MIN, 1);
inline fpext52_t const fpext52_t::one  (0x8000000000000000ull, 0x00000000,  0, 0);
inline fpext52_t const fpext52_t::none (0x8000000000000000ull, 0x00000000,  0, 1);
inline fpext52_t const fpext52_t::l2t  (0xd49a784bcd1b8afeull, 0x492bf6ff,  1, 0);
inline fpext52_t const fpext52_t::l2e  (0xb8aa3b295c17f0bbull, 0xbe87fed0,  0, 0);
inline fpext52_t const fpext52_t::pi   (0xc90fdaa22168c234ull, 0xc4c6628c,  1, 0);
inline fpext52_t const fpext52_t::pio2 (0xc90fdaa22168c234ull, 0xc4c6628c,  0, 0);
inline fpext52_t const fpext52_t::pio4 (0xc90fdaa22168c234ull, 0xc4c6628c, -1, 0);
inline fpext52_t const fpext52_t::lg2  (0x9a209a84fbcff798ull, 0x8f8959ac, -2, 0);
inline fpext52_t const fpext52_t::ln2  (0xb17217f7d1cf79abull, 0xc9e3b398, -1, 0);

template<> inline fpext64_t const fpext64_t::zero (0x0000000000000000ull, 0x00000000, fpext64_t::EXPONENT_MIN, 0);
template<> inline fpext64_t const fpext64_t::nzero(0x0000000000000000ull, 0x00000000, fpext64_t::EXPONENT_MIN, 1);
template<> inline fpext64_t const fpext64_t::one  (0x8000000000000000ull, 0x00000000,  0, 0);
template<> inline fpext64_t const fpext64_t::none (0x8000000000000000ull, 0x00000000,  0, 1);
template<> inline fpext64_t const fpext64_t::l2t  (0xd49a784bcd1b8afeull, 0x492bf6ff,  1, 0);
template<> inline fpext64_t const fpext64_t::l2e  (0xb8aa3b295c17f0bbull, 0xbe87fed0,  0, 0);
template<> inline fpext64_t const fpext64_t::pi   (0xc90fdaa22168c234ull, 0xc4c6628c,  1, 0);
template<> inline fpext64_t const fpext64_t::pio2 (0xc90fdaa22168c234ull, 0xc4c6628c,  0, 0);
template<> inline fpext64_t const fpext64_t::pio4 (0xc90fdaa22168c234ull, 0xc4c6628c, -1, 0);
template<> inline fpext64_t const fpext64_t::lg2  (0x9a209a84fbcff798ull, 0x8f8959ac, -2, 0);
template<> inline fpext64_t const fpext64_t::ln2  (0xb17217f7d1cf79abull, 0xc9e3b398, -1, 0);

template<> inline fpext96_t const fpext96_t::zero (0x0000000000000000ull, 0x00000000, fpext96_t::EXPONENT_MIN, 0);
template<> inline fpext96_t const fpext96_t::nzero(0x0000000000000000ull, 0x00000000, fpext96_t::EXPONENT_MIN, 1);
template<> inline fpext96_t const fpext96_t::one  (0x8000000000000000ull, 0x00000000,  0, 0);
template<> inline fpext96_t const fpext96_t::none (0x8000000000000000ull, 0x00000000,  0, 1);
template<> inline fpext96_t const fpext96_t::l2t  (0xd49a784bcd1b8afeull, 0x492bf6ff,  1, 0);
template<> inline fpext96_t const fpext96_t::l2e  (0xb8aa3b295c17f0bbull, 0xbe87fed0,  0, 0);
template<> inline fpext96_t const fpext96_t::pi   (0xc90fdaa22168c234ull, 0xc4c6628c,  1, 0);
template<> inline fpext96_t const fpext96_t::pio2 (0xc90fdaa22168c234ull, 0xc4c6628c,  0, 0);
template<> inline fpext96_t const fpext96_t::pio4 (0xc90fdaa22168c234ull, 0xc4c6628c, -1, 0);
template<> inline fpext96_t const fpext96_t::lg2  (0x9a209a84fbcff798ull, 0x8f8959ac, -2, 0);
template<> inline fpext96_t const fpext96_t::ln2  (0xb17217f7d1cf79abull, 0xc9e3b398, -1, 0);

template<> inline fpext128_t const fpext128_t::zero (0x0000000000000000ull, 0x0000000000000000ull, fpext128_t::EXPONENT_MIN, 0, fpext128_t::low64_tag{});
template<> inline fpext128_t const fpext128_t::nzero(0x0000000000000000ull, 0x0000000000000000ull, fpext128_t::EXPONENT_MIN, 1, fpext128_t::low64_tag{});
template<> inline fpext128_t const fpext128_t::one  (0x8000000000000000ull, 0x0000000000000000ull,  0, 0, fpext128_t::low64_tag{});
template<> inline fpext128_t const fpext128_t::none (0x8000000000000000ull, 0x0000000000000000ull,  0, 1, fpext128_t::low64_tag{});
template<> inline fpext128_t const fpext128_t::l2t  (0xd49a784bcd1b8afeull, 0x492bf6ff0e7e1f80ull,  1, 0, fpext128_t::low64_tag{});
template<> inline fpext128_t const fpext128_t::l2e  (0xb8aa3b295c17f0bbull, 0xbe87fed0691d3e88ull,  0, 0, fpext128_t::low64_tag{});
template<> inline fpext128_t const fpext128_t::pi   (0xc90fdaa22168c234ull, 0xc4c6628b80dc1cd1ull,  1, 0, fpext128_t::low64_tag{});
template<> inline fpext128_t const fpext128_t::pio2 (0xc90fdaa22168c234ull, 0xc4c6628b80dc1cd1ull,  0, 0, fpext128_t::low64_tag{});
template<> inline fpext128_t const fpext128_t::pio4 (0xc90fdaa22168c234ull, 0xc4c6628b80dc1cd1ull, -1, 0, fpext128_t::low64_tag{});
template<> inline fpext128_t const fpext128_t::lg2  (0x9a209a84fbcff798ull, 0x8f8959ac0b7c9178ull, -2, 0, fpext128_t::low64_tag{});
template<> inline fpext128_t const fpext128_t::ln2  (0xb17217f7d1cf79abull, 0xc9e3b39803f2f6afull, -1, 0, fpext128_t::low64_tag{});



//
// construct an fpex52_t from high-precision components
//
inline fpext52_t::fpext52_t(uint64_t high, uint32_t low, int32_t exponent, uint16_t sign)
{
    int32_t exp = exponent + FP64_EXPONENT_BIAS;

    // compute a signed zero as the default case
    m_value.i = uint64_t(sign) << FP64_SIGN_SHIFT;

    // too big to fit? return signed infinity
    if (exp >= FP64_EXPONENT_MAX_BIASED)
        m_value.i |= FP64_EXPONENT_MASK;

    // normal case
    else if (exp > 0)
    {
        m_value.i |= (uint64_t(exp) << FP64_EXPONENT_SHIFT) | ((high >> (63 - FP64_EXPONENT_SHIFT)) & FP64_MANTISSA_MASK);
        m_value.i += (high >> (62 - FP64_EXPONENT_SHIFT)) & 1;
    }

    // denormal case
    else if (exp > -52)
    {
        m_value.i |= high >> (64 - FP64_EXPONENT_SHIFT - exp);
        m_value.i += (63 - FP64_EXPONENT_SHIFT - exp);
    }
}



//
// construct an fpext52_t from a higher precision type
//
template<typename SrcExtendedType>
fpext52_t::fpext52_t(fpextxx_t<SrcExtendedType> const &src) :
    fpext52_t(src.mantissa(), src.extend(), src.exponent(), src.sign())
{
}



//
// copy constructor
//
template<typename ExtendedType>
template<typename SrcExtendedType>
inline fpextxx_t<ExtendedType>::fpextxx_t(fpextxx_t<SrcExtendedType> const &src, bool round) :
    m_mantissa(src.mantissa()),
    m_exponent(src.exponent()),
    m_sign(src.sign())
{
    constexpr bool SRC_EXTENDED = (sizeof(SrcExtendedType) >= 2);
    constexpr int SRC_EXTEND_BITS = 8 * sizeof(SrcExtendedType);

    // set extend to 0 if either src/dest is not extended
    if (!EXTENDED || !SRC_EXTENDED)
    {
        m_extend = 0;

        // if the source was extended, optionally round based on top extension bit
        if (SRC_EXTENDED && round)
            if ((src.extend() & (1 << (SRC_EXTEND_BITS - 1))) != 0)
                this->round_mantissa_up();
    }

    // if we have fewer extension bits, shift the source bits down
    else if (EXTEND_BITS < SRC_EXTEND_BITS)
    {
        constexpr int SRCSHIFT = (SRC_EXTEND_BITS - EXTEND_BITS) % SRC_EXTEND_BITS;
        constexpr int SRCSHIFTM1 = (SRCSHIFT + SRC_EXTEND_BITS - 1) % SRC_EXTEND_BITS;
        m_extend = src.extend() >> SRCSHIFT;

        // if rounding, check the bit we shifted out
        if (round && (src.extend() & (1 << SRCSHIFTM1)) != 0)
            this->round_extend_up();
    }

    // otherwise, shift the source up to our number of bits
    else
    {
        constexpr int SRCSHIFT = (EXTEND_BITS - SRC_EXTEND_BITS) % EXTEND_BITS;
        m_extend = src.extend() << SRCSHIFT;
    }
}



//
// construct from an fp64 type
//
template<typename ExtendedType>
inline fpextxx_t<ExtendedType>::fpextxx_t(fp64_t const &src) :
    m_mantissa(src.mantissa() << (63 - FP64_EXPONENT_SHIFT)),
    m_extend(0),
    m_exponent(src.exponent()),
    m_sign(src.sign())
{
    x87_assert(!src.ismaxexp());

    // insert the explicit one for normal numbers
    if (m_exponent != 0x000 - FP64_EXPONENT_BIAS)
        m_mantissa |= EXPLICIT_ONE;

    // normalize if we have a denorm or zero
    else
    {
        m_exponent += 1;
        this->normalize();
    }
}



//
// construct from an fp80 type
//
template<typename ExtendedType>
inline fpextxx_t<ExtendedType>::fpextxx_t(fp80_t const &src) :
    m_mantissa(src.mantissa()),
    m_extend(0),
    m_exponent(src.exponent()),
    m_sign(src.sign())
{
    x87_assert(!src.ismaxexp());

    // normalize if we have a denorm or zero
    if (int64_t(m_mantissa) >= 0)
    {
        m_exponent += 1;
        this->normalize();
    }
    x87_assert(m_sign == 0 || m_sign == 1);
}



//
// convert to an fp64
//
template<typename ExtendedType>
inline fp64_t fpextxx_t<ExtendedType>::as_fp64() const
{
    uint64_t result = uint64_t(m_sign) << FP64_SIGN_SHIFT;
    exponent_t exp = m_exponent + FP64_EXPONENT_BIAS;

    // too big to fit? return signed infinity
    if (exp >= FP64_EXPONENT_MAX_BIASED)
        return fp64_t::from_fpbits64(result | FP64_EXPONENT_MASK);

    // normal case
    else if (exp > 0)
        return fp64_t::from_fpbits64(result | (uint64_t(exp) << FP64_EXPONENT_SHIFT) | ((m_mantissa >> (63 - FP64_EXPONENT_SHIFT)) & FP64_MANTISSA_MASK));

    // denormal case
    else if (exp > -52)
        return fp64_t::from_fpbits64(result | (m_mantissa >> (64 - FP64_EXPONENT_SHIFT - exp)));

    // zero case
    return fp64_t::from_fpbits64(result);
}



//
// convert to an fp80
//
template<typename ExtendedType>
inline fp80_t fpextxx_t<ExtendedType>::as_fp80() const
{
    uint16_t sign_exp = m_sign << FP80_SIGN_SHIFT;
    exponent_t exp = m_exponent + FP80_EXPONENT_BIAS;

    // too big to fit? return signed infinity
    if (exp >= FP80_EXPONENT_MAX_BIASED)
        return fp80_t(0, sign_exp | FP80_EXPONENT_MASK);

    // normal case
    else if (exp > 0)
        return fp80_t(m_mantissa, sign_exp | exp);

    // denormal case
    else if (exp > -63)
        return fp80_t(m_mantissa >> (1 - exp), sign_exp);

    // zero case
    return fp80_t(0, sign_exp);
}



//
// compare two mantissas in various ways
//
template<typename ExtendedType>
inline bool fpextxx_t<ExtendedType>::mantissa_eq(fpextxx_t const &a) const
{
    return (m_mantissa == a.m_mantissa && (!EXTENDED || m_extend == a.m_extend));
}

template<typename ExtendedType>
inline bool fpextxx_t<ExtendedType>::mantissa_gt(fpextxx_t const &a) const
{
    return (m_mantissa > a.m_mantissa || (EXTENDED && m_mantissa == a.m_mantissa && m_extend > a.m_extend));
}

template<typename ExtendedType>
inline bool fpextxx_t<ExtendedType>::mantissa_lt(fpextxx_t const &a) const
{
    return (m_mantissa < a.m_mantissa || (EXTENDED && m_mantissa == a.m_mantissa && m_extend < a.m_extend));
}



//
// round the mantissa up, overflowing into the exponent
//
template<typename ExtendedType>
inline void fpextxx_t<ExtendedType>::round_mantissa_up()
{
    if (++m_mantissa == 0)
    {
        m_mantissa = EXPLICIT_ONE;
        m_exponent += 1;
    }
}



//
// round the extension bits up, overflowing into the mantissa
//
template<typename ExtendedType>
inline void fpextxx_t<ExtendedType>::round_extend_up()
{
    if (++m_extend == 0 && ++m_mantissa == 0)
    {
        m_mantissa = EXPLICIT_ONE;
        m_exponent += 1;
    }
}



//
// shift the mantissa 'count' bits to the right
//
template<typename ExtendedType>
inline void fpextxx_t<ExtendedType>::shift_mantissa_right(int count)
{
    if (count < EXTEND_BITS)
    {
        m_extend = extend_t((m_extend >> count) | (m_mantissa << (EXTEND_BITS - count)));
        m_mantissa >>= count;
    }
    else
    {
        m_extend = extend_t(m_mantissa >> (count - EXTEND_BITS));
        if (count < 64)
            m_mantissa >>= count;
        else
            m_mantissa = EXPLICIT_ONE;
    }
}

template<>
inline void fpextxx_t<uint8_t>::shift_mantissa_right(int count)
{
    m_mantissa >>= count;
}



//
// normalize a denormalized or zero value
//
template<typename ExtendedType>
inline void fpextxx_t<ExtendedType>::normalize()
{
    // if mantissa is all zeros, set exponent to minimum
    if (this->iszero())
        m_exponent = EXPONENT_MIN;

    // if we have bits in the mantissa, shift them high
    else if (!EXTENDED || m_mantissa != 0)
    {
        // determine shift amount
        int shift = count_leading_zeros64(m_mantissa);
        if (shift == 0)
            return;

        // shift it up
        m_mantissa <<= shift;
        m_exponent -= shift;

        // shift in extended bits if we have some
        if (EXTENDED)
        {
            if (shift < EXTEND_BITS)
            {
                m_mantissa |= m_extend >> (EXTEND_BITS - shift);
                m_extend <<= shift;
            }
            else
            {
                m_mantissa |= mantissa_t(m_extend) << (shift - EXTEND_BITS);
                m_extend = 0;
            }
        }
    }

    // if we only have bits in the extension, move them high
    else if (EXTENDED)
    {
        // determine shift amount
        int shift = count_leading_zeros64(m_extend);

        // shift it up
        m_mantissa = mantissa_t(m_extend) << shift;
        m_extend = 0;
        m_exponent -= 64 + shift - (64 - EXTEND_BITS);
    }

    // verify
    x87_assert(this->iszero() || (m_mantissa & EXPLICIT_ONE) != 0);
}



//
// add two values of the same sign, assuming src1 is the larger value
//
template<typename ExtendedType>
inline void fpextxx_t<ExtendedType>::add_values(fpextxx_t const &src1, fpextxx_t const &src2, int src2shift)
{
    // if src2 is way too small, treat as zero
    if (src2shift >= MANTISSA_BITS)
    {
        m_mantissa = src1.m_mantissa;
        m_exponent = src1.m_exponent;
        m_extend = src1.m_extend;
        return;
    }

    // shift the second value
    extend_t src2e;
    mantissa_t src2m;
    if (src2shift == 0)
    {
        src2e = src2.m_extend;
        src2m = src2.m_mantissa;
    }
    else if (src2shift < EXTEND_BITS)
    {
        src2e = extend_t((src2.m_extend >> src2shift) | (src2.m_mantissa << (EXTEND_BITS - src2shift)));
        src2m = src2.m_mantissa >> src2shift;
        if (((src2.m_extend & (1 << (src2shift - 1))) != 0) && ++src2e == 0) src2m++;
    }
    else
    {
        src2e = extend_t(src2.m_mantissa >> (src2shift - EXTEND_BITS));
        src2m = (src2shift < 64) ? (src2.m_mantissa >> src2shift) : 0;
        if (src2shift != EXTEND_BITS && ((src2.m_mantissa & (1ull << (src2shift - EXTEND_BITS - 1))) != 0) && ++src2e == 0) src2m++;
    }

    // add the main mantissa and note carry out
    m_mantissa = src1.m_mantissa + src2m;
    m_exponent = src1.m_exponent;
    bool carry = (m_mantissa < src2m);

    // add the extension and handle carry out
    m_extend = src1.m_extend + src2e;
    if (m_extend < src2e)
        this->round_mantissa_up();

    // handle carry out
    if (carry)
    {
        this->shift_mantissa_right(1);
        m_mantissa |= EXPLICIT_ONE;
        m_exponent += 1;
    }
}

template<>
inline void fpextxx_t<uint8_t>::add_values(fpextxx_t const &src1, fpextxx_t const &src2, int src2shift)
{
    // if src2 is way too small, treat as zero
    if (src2shift >= MANTISSA_BITS)
    {
        m_mantissa = src1.m_mantissa;
        m_exponent = src1.m_exponent;
        return;
    }

    // shift the second value
    mantissa_t src2m = src2.m_mantissa >> src2shift;
    if (src2shift != 0 && (src2.m_mantissa & (1ull << (src2shift - 1))) != 0) src2m++;

    // add and handle overflow
    m_mantissa = src1.m_mantissa + src2m;
    m_exponent = src1.m_exponent;
    if (m_mantissa < src2m)
    {
        this->shift_mantissa_right(1);
        m_mantissa |= EXPLICIT_ONE;
        m_exponent += 1;
    }
}



//
// subtract the source mantissa, normalizing
//
template<typename ExtendedType>
inline void fpextxx_t<ExtendedType>::sub_values(fpextxx_t const &src1, fpextxx_t const &src2, int src2shift)
{
    // if src2 is way too small, treat as zero
    if (src2shift >= MANTISSA_BITS)
    {
        m_mantissa = src1.m_mantissa;
        m_exponent = src1.m_exponent;
        m_extend = src1.m_extend;
        return;
    }

    // shift the second value
    extend_t src2e;
    mantissa_t src2m;
    if (src2shift == 0)
    {
        src2e = src2.m_extend;
        src2m = src2.m_mantissa;
    }
    else if (src2shift < EXTEND_BITS)
    {
        src2e = extend_t((src2.m_extend >> src2shift) | (src2.m_mantissa << (EXTEND_BITS - src2shift)));
        src2m = src2.m_mantissa >> src2shift;
        if (((src2.m_extend & (1 << (src2shift - 1))) != 0) && ++src2e == 0) src2m++;
    }
    else
    {
        src2e = extend_t(src2.m_mantissa >> (src2shift - EXTEND_BITS));
        src2m = (src2shift < 64) ? (src2.m_mantissa >> src2shift) : 0;
        if (src2shift != EXTEND_BITS && ((src2.m_mantissa & (1ull << (src2shift - EXTEND_BITS - 1))) != 0) && ++src2e == 0) src2m++;
    }

    // subtract the main mantissa
    m_mantissa = src1.m_mantissa - src2m;
    m_exponent = src1.m_exponent;

    // subtract the extension and borrow
    auto orig = src1.m_extend;
    m_extend = orig - src2e;
    if (m_extend > orig)
        m_mantissa--;

    // normalize
    this->normalize();
}

template<>
inline void fpextxx_t<uint8_t>::sub_values(fpextxx_t const &src1, fpextxx_t const &src2, int src2shift)
{
    // if src2 is way too small, treat as zero
    if (src2shift >= MANTISSA_BITS)
    {
        m_mantissa = src1.m_mantissa;
        m_exponent = src1.m_exponent;
        return;
    }

    // shift the second value
    mantissa_t src2m = src2.m_mantissa >> src2shift;
    if (src2shift != 0 && (src2.m_mantissa & (1ull << (src2shift - 1))) != 0) src2m++;

    // subtract
    m_mantissa = src1.m_mantissa - src2m;
    m_exponent = src1.m_exponent;

    // normalize
    this->normalize();
}



//
// equality comparisons
//
template<typename ExtendedType>
inline bool fpextxx_t<ExtendedType>::operator==(fpextxx_t const &rhs) const
{
    return (m_mantissa == rhs.m_mantissa &&
            (!EXTENDED || m_extend == rhs.m_extend) &&
            m_exponent == rhs.m_exponent &&
            m_sign == rhs.m_sign);
}

template<typename ExtendedType>
inline bool fpextxx_t<ExtendedType>::operator!=(fpextxx_t const &rhs) const
{
    return (m_mantissa != rhs.m_mantissa ||
            (EXTENDED && m_extend != rhs.m_extend) ||
            m_exponent != rhs.m_exponent ||
            m_sign != rhs.m_sign);
}



//
// less than comparisons
//
template<typename ExtendedType>
inline bool fpextxx_t<ExtendedType>::operator<(fpextxx_t const &rhs) const
{
    if ((m_sign ^ rhs.m_sign) != 0)
        return m_sign;
    if (m_exponent < rhs.m_exponent)
        return true;
    if (m_exponent > rhs.m_exponent)
        return false;
    if (!EXTENDED)
        return (m_mantissa < rhs.m_mantissa);
    if (m_mantissa < rhs.m_mantissa)
        return true;
    if (m_mantissa > rhs.m_mantissa)
        return false;
    return (m_extend < rhs.m_extend);
}

template<typename ExtendedType>
inline bool fpextxx_t<ExtendedType>::operator<=(fpextxx_t const &rhs) const
{
    if ((m_sign ^ rhs.m_sign) != 0)
        return m_sign;
    if (m_exponent < rhs.m_exponent)
        return true;
    if (m_exponent > rhs.m_exponent)
        return false;
    if (!EXTENDED)
        return (m_mantissa <= rhs.m_mantissa);
    if (m_mantissa < rhs.m_mantissa)
        return true;
    if (m_mantissa > rhs.m_mantissa)
        return false;
    return (m_extend <= rhs.m_extend);
}



//
// greater than comparisons
//
template<typename ExtendedType>
inline bool fpextxx_t<ExtendedType>::operator>(fpextxx_t const &rhs) const
{
    if ((m_sign ^ rhs.m_sign) != 0)
        return rhs.m_sign;
    if (m_exponent > rhs.m_exponent)
        return true;
    if (m_exponent < rhs.m_exponent)
        return false;
    if (!EXTENDED)
        return (m_mantissa > rhs.m_mantissa);
    if (m_mantissa > rhs.m_mantissa)
        return true;
    if (m_mantissa < rhs.m_mantissa)
        return false;
    return (m_extend > rhs.m_extend);
}

template<typename ExtendedType>
inline bool fpextxx_t<ExtendedType>::operator>=(fpextxx_t const &rhs) const
{
    if ((m_sign ^ rhs.m_sign) != 0)
        return rhs.m_sign;
    if (m_exponent > rhs.m_exponent)
        return true;
    if (m_exponent < rhs.m_exponent)
        return false;
    if (!EXTENDED)
        return (m_mantissa >= rhs.m_mantissa);
    if (m_mantissa > rhs.m_mantissa)
        return true;
    if (m_mantissa < rhs.m_mantissa)
        return false;
    return (m_extend >= rhs.m_extend);
}



//
// perform addition between two source values
//
template<typename ExtendedType>
inline void fpextxx_t<ExtendedType>::add(fpextxx_t const &a, fpextxx_t const &b)
{
    // get difference in signs
    sign_t signdiff = a.m_sign ^ b.m_sign;

    // compute difference in exponents
    exponent_t dexp = a.m_exponent - b.m_exponent;

    // handle same-sign
    if (signdiff == 0)
    {
        m_sign = a.m_sign;

        // A >= B
        if (dexp >= 0)
            this->add_values(a, b, dexp);

        // A < B
        else
            this->add_values(b, a, -dexp);
    }

    // handle opposite sign
    else
    {
        // A >= B
        if (dexp > 0 || (dexp == 0 && !a.mantissa_lt(b)))
        {
            m_sign = a.m_sign;
            this->sub_values(a, b, dexp);
        }

        // A < B
        else
        {
            m_sign = b.m_sign;
            this->sub_values(b, a, -dexp);
        }
    }

    // verify
    x87_assert(this->iszero() || (m_mantissa & EXPLICIT_ONE) != 0);
}



//
// perform subtraction between two source values
//
template<typename ExtendedType>
inline void fpextxx_t<ExtendedType>::sub(fpextxx_t const &a, fpextxx_t const &b)
{
    // get difference in signs
    sign_t signdiff = a.m_sign ^ b.m_sign;

    // compute difference in exponents
    exponent_t dexp = a.m_exponent - b.m_exponent;

    // handle same-sign
    if (signdiff != 0)
    {
        m_sign = a.m_sign;

        // A >= B
        if (dexp >= 0)
            this->add_values(a, b, dexp);

        // A < B
        else
            this->add_values(b, a, -dexp);
    }

    // handle opposite sign
    else
    {
        // A >= B
        if (dexp > 0 || (dexp == 0 && !a.mantissa_lt(b)))
        {
            m_sign = a.m_sign;
            this->sub_values(a, b, dexp);
        }

        // A < B
        else
        {
            m_sign = a.m_sign ^ 1;
            this->sub_values(b, a, -dexp);
        }
    }

    // verify
    x87_assert(this->iszero() || (m_mantissa & EXPLICIT_ONE) != 0);
}



//
// perform multiplication between two source values
//
template<typename ExtendedType>
inline void fpextxx_t<ExtendedType>::mul(fpextxx_t const &a, fpextxx_t const &b)
{
    // compute final sign
    m_sign = a.m_sign ^ b.m_sign;

    // check for 0
    if (a.iszero() || b.iszero())
    {
        m_exponent = EXPONENT_MIN;
        m_mantissa = 0;
        m_extend = 0;
        return;
    }

    // compute 64x64 mantissa multiplication
    auto [lo, hi] = multiply_64x64(a.m_mantissa, b.m_mantissa);

    // compute A.hi * B.lo and B.hi * A.lo
    auto [lo1, hi1] = multiply_64x64(a.m_mantissa, uint32_t(b.m_extend << (32 - EXTEND_BITS)));
    auto [lo2, hi2] = multiply_64x64(b.m_mantissa, uint32_t(a.m_extend << (32 - EXTEND_BITS)));

    // add upper parts; guaranteed not to overflow (max 64+32 bits each)
    uint64_t hiadd = hi1 + hi2;

    // add lower parts; carry into upper part
    uint64_t loadd = lo1 + lo2;
    if (loadd < lo2)
        hiadd++;

    // also compute A.lo * B.lo
    uint64_t lo3 = (uint64_t(a.m_extend) * uint64_t(b.m_extend)) >> (EXTEND_BITS*2 - 32);

    // add in this contribution, carrying into hi
    loadd += lo3;
    if (loadd < lo3)
        hiadd++;

    // shift remaining bits down and borrow 32 from hiadd
    loadd = (loadd >> 32) | (hiadd << 32);
    hiadd >>= 32;

    // add to the lower part and carry if necessary
    lo += loadd;
    if (lo < loadd)
        hi++;
    hi += hiadd;

    // compute final exponent
    m_exponent = a.m_exponent + b.m_exponent;

    // adjust for overflow
    if ((hi & EXPLICIT_ONE) == 0)
    {
        m_mantissa = (hi << 1) | (lo >> 63);
        m_extend = extend_t(lo >> (63 - EXTEND_BITS));
        if ((lo & (1ull << (63 - EXTEND_BITS - 1))) != 0) this->round_extend_up();
    }
    else
    {
        m_mantissa = hi;
        m_extend = extend_t(lo >> (64 - EXTEND_BITS));
        m_exponent += 1;
        if ((lo & (1ull << (64 - EXTEND_BITS - 1))) != 0) this->round_extend_up();
    }

    // double check to be sure we ended up as expected
    x87_assert((m_mantissa & EXPLICIT_ONE) != 0);
}

template<>
inline void fpextxx_t<uint8_t>::mul(fpextxx_t const &a, fpextxx_t const &b)
{
    // compute final sign
    m_sign = a.m_sign ^ b.m_sign;

    // check for 0
    if (a.iszero() || b.iszero())
    {
        m_exponent = EXPONENT_MIN;
        m_mantissa = 0;
        return;
    }

    // compute 64x64 mantissa multiplication
    auto [lo, hi] = multiply_64x64(a.m_mantissa, b.m_mantissa);

    // compute final exponent
    m_exponent = a.m_exponent + b.m_exponent;

    // adjust for overflow
    if ((hi & EXPLICIT_ONE) == 0)
        m_mantissa = ((hi << 1) | (lo >> 63)) + ((lo >> 62) & 1);
    else
        m_mantissa = hi + ((lo >> 63) & 1), m_exponent += 1;

    // double check to be sure we ended up as expected
    x87_assert((m_mantissa & EXPLICIT_ONE) != 0);
}

//
// fpext128_t multiply: full 128x128 -> 256-bit intermediate, retaining
// the top 128 bits as result and using the next bits for proper round-
// to-nearest-even.
//
template<>
inline void fpextxx_t<uint64_t>::mul(fpextxx_t const &a, fpextxx_t const &b)
{
    m_sign = a.m_sign ^ b.m_sign;
    if (a.iszero() || b.iszero())
    {
        m_exponent = EXPONENT_MIN;
        m_mantissa = 0;
        m_extend = 0;
        return;
    }

    // Schoolbook 128x128 -> 256 with four 64x64 partial products.
    // A = a.m_mantissa : a.m_extend  (high : low)
    // B = b.m_mantissa : b.m_extend
    auto [pp_lo, pp_hi] = multiply_64x64(a.m_mantissa, b.m_mantissa); // top
    auto [mE_lo, mE_hi] = multiply_64x64(a.m_mantissa, b.m_extend);
    auto [eM_lo, eM_hi] = multiply_64x64(a.m_extend,   b.m_mantissa);
    auto [ee_lo, ee_hi] = multiply_64x64(a.m_extend,   b.m_extend);   // bottom

    // Accumulator layout (each tN is 64 bits):
    //   t3 (bits 192..255) = pp_hi + carries
    //   t2 (bits 128..191) = pp_lo + mE_hi + eM_hi
    //   t1 (bits  64..127) = mE_lo + eM_lo + ee_hi
    //   t0 (bits   0..63 ) = ee_lo
    uint64_t t0 = ee_lo;
    uint64_t t1 = mE_lo;
    uint64_t carry_t1 = 0;
    t1 += eM_lo; if (t1 < eM_lo) carry_t1++;
    t1 += ee_hi; if (t1 < ee_hi) carry_t1++;

    uint64_t t2 = pp_lo;
    uint64_t carry_t2 = 0;
    t2 += mE_hi;     if (t2 < mE_hi)     carry_t2++;
    t2 += eM_hi;     if (t2 < eM_hi)     carry_t2++;
    t2 += carry_t1;  if (t2 < carry_t1)  carry_t2++;

    uint64_t t3 = pp_hi + carry_t2;

    m_exponent = a.m_exponent + b.m_exponent;

    // Normalize: result starts in t3:t2 (with t1:t0 as residue).
    // Mantissas in [2^63, 2^64) → 128-bit values in [2^127, 2^128). Their
    // product is in [2^254, 2^256), so t3's bit 63 is either set (no
    // normalization) or one below (need to shift left by 1).
    bool round_bit, sticky;
    if ((t3 & EXPLICIT_ONE) == 0)
    {
        // Shift everything left by 1.
        m_mantissa = (t3 << 1) | (t2 >> 63);
        m_extend   = (t2 << 1) | (t1 >> 63);
        round_bit  = (t1 >> 62) & 1;
        sticky     = ((t1 & ((1ull << 62) - 1)) != 0) || (t0 != 0);
    }
    else
    {
        m_mantissa = t3;
        m_extend   = t2;
        m_exponent += 1;
        round_bit  = (t1 >> 63) & 1;
        sticky     = ((t1 & ((1ull << 63) - 1)) != 0) || (t0 != 0);
    }
    // Round-to-nearest-even at the boundary below m_extend's LSB.
    if (round_bit && (sticky || (m_extend & 1)))
        this->round_extend_up();

    x87_assert((m_mantissa & EXPLICIT_ONE) != 0);
}



//
// compute the floor of a value
//
template<typename ExtendedType>
inline fpextxx_t<ExtendedType> fpextxx_t<ExtendedType>::floor(fpextxx_t const &a)
{
    mantissa_t mantissa = a.mantissa();
    extend_t extend = a.extend();
    exponent_t exp = a.exponent();

    // positive case
    if (a.sign() == 0)
    {
        // negative exponents are < 1.0, and thus floor to 0
        if (exp < 0)
            return fpextxx_t::zero;

        // small exponents mask below
        if (exp <= MANTISSA_BITS - 1)
        {
            int shift = MANTISSA_BITS - 1 - exp;
            extend_t extend_mask = (shift < EXTEND_BITS) ? ~((1 << shift) - 1) : 0;
            mantissa_t mantissa_mask = (shift > EXTEND_BITS) ? ~((1ull << (shift - EXTEND_BITS)) - 1) : ~0ull;
            return fpextxx_t(mantissa & mantissa_mask, extend & extend_mask, exp, 0);
        }

        // large exponents have nothing to floor
        return a;
    }

    // negative case
    else
    {
        // negative exponents are > -1.0, and thus floor to -1
        if (exp < 0)
            return fpextxx_t::none;

        // small exponents mask below after adding the maximal amount
        if (exp <= MANTISSA_BITS - 1)
        {
            int shift = MANTISSA_BITS - 1 - exp;
            extend_t extend_mask = (shift < EXTEND_BITS) ? ~((1 << shift) - 1) : 0;
            mantissa_t mantissa_mask = (shift > EXTEND_BITS) ? ~((1ull << (shift - EXTEND_BITS)) - 1) : ~0ull;
            extend_t extend_sum = extend + ~extend_mask;
            mantissa_t mantissa_sum = mantissa + ~mantissa_mask + (extend_sum < extend);
            if (mantissa_sum < mantissa)
                return fpextxx_t<uint8_t>(((mantissa_sum & mantissa_mask) >> 1) | EXPLICIT_ONE, extend_t((extend_sum >> 1) | (mantissa_sum << (EXTEND_BITS - 1))), exp + 1, 1);
            else
                return fpextxx_t(mantissa_sum & mantissa_mask, extend_sum & extend_mask, exp, 1);
        }

        // large exponents have nothing to floor
        return a;
    }
}

template<>
inline fpextxx_t<uint8_t> fpextxx_t<uint8_t>::floor(fpextxx_t const &a)
{
    mantissa_t mantissa = a.mantissa();
    exponent_t exp = a.exponent();

    // positive case
    if (a.sign() == 0)
    {
        // negative exponents are < 1.0, and thus floor to 0
        if (exp < 0)
            return fpextxx_t::zero;

        // small exponents mask below
        if (exp <= MANTISSA_BITS - 1)
        {
            int shift = MANTISSA_BITS - 1 - exp;
            mantissa_t mantissa_mask = ~((1ull << shift) - 1);
            return fpextxx_t<uint8_t>(mantissa & mantissa_mask, 0, exp, 0);
        }

        // large exponents have nothing to floor
        return a;
    }

    // negative case
    else
    {
        // negative exponents are > -1.0, and thus floor to -1
        if (exp < 0)
            return fpextxx_t<uint8_t>::none;

        // small exponents mask below after adding the maximal amount
        if (exp <= MANTISSA_BITS - 1)
        {
            int shift = MANTISSA_BITS - 1 - exp;
            mantissa_t mantissa_mask = ~((1ull << shift) - 1);
            mantissa_t mantissa_sum = mantissa + ~mantissa_mask;
            if (mantissa_sum < mantissa)
                return fpextxx_t<uint8_t>(((mantissa_sum & mantissa_mask) >> 1) | EXPLICIT_ONE, 0, exp + 1, 1);
            else
                return fpextxx_t<uint8_t>(mantissa_sum & mantissa_mask, 0, exp, 1);
        }

        // large exponents have nothing to floor
        return a;
    }
}



//
// compute the floor the absolute value of a number, plus return
// the low integral bits (used for trig functions)
//
template<typename ExtendedType>
inline fpextxx_t<ExtendedType> fpextxx_t<ExtendedType>::floor_abs_loint(fpextxx_t const &a, uint64_t &intbits)
{
    mantissa_t mantissa = a.mantissa();
    extend_t extend = a.extend();
    exponent_t exp = a.exponent();

    x87_assert(exp < 63);

    // negative exponents are < 1.0, and thus floor to 0
    if (exp < 0)
    {
        intbits = 0;
        return fpextxx_t::zero;
    }

    // mask mantissa bits, with extend == 0
    int shift = MANTISSA_BITS - 1 - exp;
    if (shift >= EXTEND_BITS)
    {
        shift -= EXTEND_BITS;
        mantissa_t mantissa_mask = ~((1ull << shift) - 1);
        intbits = mantissa >> shift;
        return fpextxx_t(mantissa & mantissa_mask, 0, exp, 0);
    }

    // mask only extend bits
    extend_t extend_mask = ~((1 << shift) - 1);
    intbits = int64_t((extend >> shift) | (mantissa << (EXTEND_BITS - shift)));
    return fpextxx_t(mantissa, extend & extend_mask, exp, 0);
}

template<>
inline fpextxx_t<uint8_t> fpextxx_t<uint8_t>::floor_abs_loint(fpextxx_t const &a, uint64_t &intbits)
{
    mantissa_t mantissa = a.mantissa();
    exponent_t exp = a.exponent();

    x87_assert(exp < 63);

    // negative exponents are < 1.0, and thus floor to 0
    if (exp < 0)
    {
        intbits = 0;
        return fpextxx_t::zero;
    }

    // small exponents mask below
    int shift = MANTISSA_BITS - 1 - exp;
    mantissa_t mantissa_mask = ~((1ull << shift) - 1);
    intbits = mantissa >> shift;
    return fpextxx_t<uint8_t>(mantissa & mantissa_mask, 0, exp, 0);
}

}

#endif
