//=========================================================
//  x87fp80.h
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

#pragma once

#ifndef X87FP80_H
#define X87FP80_H

#include "x87common.h"


//===========================================================================
//
// x87::fp80_t
//
// Class for handling x87 floating-point values using full 80-bit values.
// This code is mostly just a stub, with lots of stuff to be implemented.
//
// NOTE: much of this, apart from the loads/stores, is not yet implemented
//
//===========================================================================

namespace x87
{

struct fp64_t;

//
// packing needed to get this to come out as 10 bytes; most compilers will
// round to 16 bytes without the pack
//
#pragma pack(push)
#pragma pack(2)

//
// struct representing an 80-bit FPU value
//
struct fp80_t
{
public:
    //
    // construction/destruction
    //
    explicit fp80_t() { }
    explicit fp80_t(uint64_t man, uint16_t se) : m_mantissa(man), m_sign_exp(se) { }
    fp80_t(fp80_t const &v80) : m_mantissa(v80.m_mantissa), m_sign_exp(v80.m_sign_exp) { }
    explicit fp80_t(fp64_t const &v64) { x87sw_t sw; this->x87_fld64(fpround_t::get(), sw, *this, &v64); }
    explicit fp80_t(double _val) { x87sw_t sw; this->x87_fld64(fpround_t::get(), sw, *this, &_val); }
    explicit fp80_t(float _val) { x87sw_t sw; this->x87_fld32(fpround_t::get(), sw, *this, &_val); }
    explicit fp80_t(int64_t _val) { x87sw_t sw; this->x87_fild64(fpround_t::get(), sw, *this, &_val); }
    explicit fp80_t(int32_t _val) { x87sw_t sw; this->x87_fild32(fpround_t::get(), sw, *this, &_val); }
    explicit fp80_t(int16_t _val) { x87sw_t sw; this->x87_fild16(fpround_t::get(), sw, *this, &_val); }

    //
    // operators
    //
    fp80_t &operator=(fp80_t const &src) { m_mantissa = src.m_mantissa; m_sign_exp = src.m_sign_exp; return *this; }
    fp80_t &operator+=(fp80_t const &rhs);
    fp80_t &operator-=(fp80_t const &rhs);
    fp80_t &operator*=(fp80_t const &rhs);
    fp80_t &operator/=(fp80_t const &rhs);

    //
    // comparison operations (IEEE 754 semantics: NaN is unordered, ±0 are equal)
    //
    bool operator==(fp80_t const &rhs) const;
    bool operator!=(fp80_t const &rhs) const;
    bool operator< (fp80_t const &rhs) const;
    bool operator<=(fp80_t const &rhs) const;
    bool operator> (fp80_t const &rhs) const;
    bool operator>=(fp80_t const &rhs) const;

    //
    // pieces
    //
    uint64_t mantissa() const { return m_mantissa; }
    uint16_t sign_exp() const { return m_sign_exp; }
    int32_t exponent() const { return (m_sign_exp & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS; }
    uint8_t sign() const { return m_sign_exp >> FP80_SIGN_SHIFT; }
    uint64_t as_fpbits64() const;

    //
    // conversion
    //
    int16_t as_int16(x87cw_t round = fpround_t::get()) const { x87sw_t sw; int16_t res; this->x87_fist16(fpround_t::get(), sw, &res, *this); return res; }
    int32_t as_int32(x87cw_t round = fpround_t::get()) const { x87sw_t sw; int32_t res; this->x87_fist32(fpround_t::get(), sw, &res, *this); return res; }
    int64_t as_int64(x87cw_t round = fpround_t::get()) const { x87sw_t sw; int64_t res; this->x87_fist64(fpround_t::get(), sw, &res, *this); return res; }
    float as_float(x87cw_t round = fpround_t::get()) const { x87sw_t sw; float res; this->x87_fst32(fpround_t::get(), sw, &res, *this); return res; }
    double as_double(x87cw_t round = fpround_t::get()) const { x87sw_t sw; double res; this->x87_fst64(fpround_t::get(), sw, &res, *this); return res; }
    fp80_t as_fp80() const { return *this; }

    //
    // misc
    //
    bool isnormal() const { return (((m_sign_exp + 1) & 0x7ffe) != 0 || this->iszero()); }
    bool isminexp() const { return ((m_sign_exp & FP80_EXPONENT_MASK) == 0); }
    bool ismaxexp() const { return ((m_sign_exp & FP80_EXPONENT_MASK) == FP80_EXPONENT_MASK); }
    bool isnan() const { return (this->ismaxexp() && (m_mantissa & FP80_MANTISSA_MASK) != 0); }
    bool isqnan() const { return (this->ismaxexp() && (m_mantissa & FP80_MANTISSA_MASK) >= 0x4000000000000000ull); }
    bool issnan() const { return this->isnan() && !isqnan(); }
    bool isinf() const { return (this->ismaxexp() && (m_mantissa & FP80_MANTISSA_MASK) == 0); }
    bool ispinf() const { return (m_sign_exp == 0x7fff && (m_mantissa & FP80_MANTISSA_MASK) == 0); }
    bool isninf() const { return (m_sign_exp == 0xffff && (m_mantissa & FP80_MANTISSA_MASK) == 0); }
    bool iszero() const { return (this->isminexp() && m_mantissa == 0); }
    bool isdenorm() const { return (this->isminexp() && m_mantissa != 0); }
    fp80_t &copysign(fp80_t const &src) { m_sign_exp = (m_sign_exp & FP80_EXPONENT_MASK) | (src.m_sign_exp & FP80_SIGN_MASK); return *this; }

    //
    // static constants
    //
    static fp80_t const_zero()  { static fp80_t const c(0x0000000000000000, 0x0000); return c; }
    static fp80_t const_nzero() { static fp80_t const c(0x0000000000000000, 0x8000); return c; }
    static fp80_t const_one()   { static fp80_t const c(0x8000000000000000, 0x3fff); return c; }
    static fp80_t const_l2t()   { static fp80_t const c(0xd49a784bcd1b8afe, 0x4000); return c; }
    static fp80_t const_l2e()   { static fp80_t const c(0xb8aa3b295c17f0bc, 0x3fff); return c; }
    static fp80_t const_pi()    { static fp80_t const c(0xc90fdaa22168c235, 0x4000); return c; }
    static fp80_t const_lg2()   { static fp80_t const c(0x9a209a84fbcff799, 0x3ffd); return c; }
    static fp80_t const_ln2()   { static fp80_t const c(0xb17217f7d1cf79ac, 0x3ffe); return c; }
    static fp80_t const_snan()  { static fp80_t const c(0x8000000000000001, 0x7fff); return c; }
    static fp80_t const_qnan()  { static fp80_t const c(0xc000000000000001, 0x7fff); return c; }
    static fp80_t const_pinf()  { static fp80_t const c(0x8000000000000000, 0x7fff); return c; }
    static fp80_t const_ninf()  { static fp80_t const c(0x8000000000000000, 0xffff); return c; }
    static fp80_t const_indef() { static fp80_t const c(0xc000000000000000, 0xffff); return c; }

    //
    // static unary ops
    //
    static fp80_t abs(fp80_t const &src) { fp80_t res = src; res.m_sign_exp &= ~FP80_SIGN_MASK; return res; }
    static fp80_t chs(fp80_t const &src) { fp80_t res = src; res.m_sign_exp ^= FP80_SIGN_MASK; return res; }
    static fp80_t sqrt(fp80_t const &src);   // STUB: not yet implemented
    static fp80_t floor(fp80_t const &src);  // STUB: not yet implemented
    static fp80_t ceil(fp80_t const &src);   // STUB: not yet implemented

    //
    // static transcendental ops
    //
    static fp80_t ldexp(fp80_t const &a, int32_t factor);  // STUB: not yet implemented

    //
    // x87 ops
    //
    static uint16_t x87_fxtract(fp80_t const &src, fp80_t &dst1, fp80_t &dst2);                    // STUB
    static uint16_t x87_fscale(fp80_t const &src1, fp80_t const &src2, fp80_t &dst);               // STUB
    static uint16_t x87_fprem(fp80_t const &src1, fp80_t const &src2, fp80_t &dst);                // STUB
    static uint16_t x87_fprem1(fp80_t const &src1, fp80_t const &src2, fp80_t &dst);               // STUB
    static uint16_t x87_f2xm1(fp80_t const &src, fp80_t &dst);                                     // implemented
    static uint16_t x87_fyl2x(fp80_t const &src1, fp80_t const &src2, fp80_t &dst);                // STUB
    static uint16_t x87_fyl2xp1(fp80_t const &src1, fp80_t const &src2, fp80_t &dst);              // STUB
    static uint16_t x87_fsin(fp80_t const &src, fp80_t &dst);                                      // STUB
    static uint16_t x87_fcos(fp80_t const &src, fp80_t &dst);                                      // STUB
    static uint16_t x87_fsincos(fp80_t const &src, fp80_t &dst1, fp80_t &dst2);                    // STUB
    static uint16_t x87_fptan(fp80_t const &src, fp80_t &dst1, fp80_t &dst2);                      // STUB
    static uint16_t x87_fpatan(fp80_t const &src1, fp80_t const &src2, fp80_t &dst);               // STUB
    static uint16_t x87_frndint(fp80_t const &src, fp80_t &dst);                                   // STUB

    // SW-returning arithmetic. operator+/-/*// don't carry a status word;
    // these entry points expose what the FPU actually puts in C1 (roundup),
    // and the precision/under/over/invalid exception bits. STUBS for now.
    static uint16_t x87_fadd (fp80_t const &src1, fp80_t const &src2, fp80_t &dst);
    static uint16_t x87_fsub (fp80_t const &src1, fp80_t const &src2, fp80_t &dst);
    static uint16_t x87_fsubr(fp80_t const &src1, fp80_t const &src2, fp80_t &dst);
    static uint16_t x87_fmul (fp80_t const &src1, fp80_t const &src2, fp80_t &dst);
    static uint16_t x87_fdiv (fp80_t const &src1, fp80_t const &src2, fp80_t &dst);
    static uint16_t x87_fdivr(fp80_t const &src1, fp80_t const &src2, fp80_t &dst);
    static uint16_t x87_fsqrt(fp80_t const &src, fp80_t &dst);

    // Classification / compare entry points. Each returns the resulting SW;
    // x87_fxam sets only C0/C1/C2/C3 (no other side effect), x87_ftst is the
    // compare-against-0 form, x87_fcom signals on QNaN, x87_fucom doesn't.
    // x87_fcomi is the P6+ variant that writes EFLAGS-style result bits;
    // we return them packed into the same C0/C2/C3 layout for test parity.
    static uint16_t x87_fxam  (fp80_t const &src);
    static uint16_t x87_ftst  (fp80_t const &src);
    static uint16_t x87_fcom  (fp80_t const &src1, fp80_t const &src2);
    static uint16_t x87_fucom (fp80_t const &src1, fp80_t const &src2);
    static uint16_t x87_fcomi (fp80_t const &src1, fp80_t const &src2);
    static uint16_t x87_fucomi(fp80_t const &src1, fp80_t const &src2);

    //
    // floating point load helpers
    //
    template<typename Type> static void x87_fld_common(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src);
    static void x87_fld80(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src);
    static void x87_fld64(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src) { x87_fld_common<uint64_t>(cw, sw, dst, src); }
    static void x87_fld32(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src) { x87_fld_common<uint32_t>(cw, sw, dst, src); }

    //
    // integral load helpers
    //
    template<typename Type> static void x87_fild_common(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src);
    static void x87_fild64(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src) { x87_fild_common<int64_t>(cw, sw, dst, src); }
    static void x87_fild32(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src) { x87_fild_common<int32_t>(cw, sw, dst, src); }
    static void x87_fild16(x87cw_t cw, x87sw_t &sw, fp80_t &dst, void const *src) { x87_fild_common<int16_t>(cw, sw, dst, src); }

    //
    // floating point store helpers
    //
    template<typename Type> static void x87_fst_common(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src);
    static void x87_fst80(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src);
    static void x87_fst64(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src) { x87_fst_common<uint64_t>(cw, sw, dst, src); }
    static void x87_fst32(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src) { x87_fst_common<uint32_t>(cw, sw, dst, src); }

    //
    // integral store helpers
    //
    template<typename Type> static void x87_fist_common(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src);
    static void x87_fist64(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src) { x87_fist_common<int64_t>(cw, sw, dst, src); }
    static void x87_fist32(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src) { x87_fist_common<int32_t>(cw, sw, dst, src); }
    static void x87_fist16(x87cw_t cw, x87sw_t &sw, void *dst, fp80_t const &src) { x87_fist_common<int16_t>(cw, sw, dst, src); }

    //
    // static misc ops
    //
    static fp80_t make_qnan(fp80_t const &src) { x87_assert(src.isnan()); fp80_t result(src); result.m_mantissa |= 0xc000000000000000ull; return result; }
    // NYI static fp80_t from_fpbits32(uint32_t bits);
    // NYI static fp80_t from_fpbits64(uint64_t bits);
    static bool samesign(fp80_t const &src1, fp80_t const &src2) { return (((src1.m_sign_exp ^ src2.m_sign_exp) & FP80_SIGN_MASK) == 0); }

protected:
    //
    // internal state
    //
    uint64_t m_mantissa;
    uint16_t m_sign_exp;
};
static_assert(sizeof(fp80_t) == 10);

#pragma pack(pop)

//
// core math operations
//
fp80_t operator+(fp80_t const &a, fp80_t const &b);
fp80_t operator-(fp80_t const &a, fp80_t const &b);
fp80_t operator*(fp80_t const &a, fp80_t const &b);
fp80_t operator/(fp80_t const &a, fp80_t const &b);

inline fp80_t &fp80_t::operator+=(fp80_t const &rhs) { *this = *this + rhs; return *this; }
inline fp80_t &fp80_t::operator-=(fp80_t const &rhs) { *this = *this - rhs; return *this; }
inline fp80_t &fp80_t::operator*=(fp80_t const &rhs) { *this = *this * rhs; return *this; }
inline fp80_t &fp80_t::operator/=(fp80_t const &rhs) { *this = *this / rhs; return *this; }

}

#endif
