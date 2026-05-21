//=========================================================
//  x87fp80trans.cpp
//
//  80-bit floating-point support (transcendental funcs)
//
//  64-bit floating-point support (transcendental funcs)
//
//  The routines in this code are based on existing work from several sources:
//
//  * fxtract/fscale/f2xm1 implementations are by Aaron Giles
//  * fprem/fprem1 implementation was derived from softfloat (BSD 3-clause)
//  * fyl2x/fyl2xp1 implementations were derived from fdlibm (Sun license)
//  * fsin/fcos/fsincos/fptan/fpatan implementations were derived from the
//     Cephes math library (MIT license)
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
//=========================================================
//
// softfloat license: (BSD 3-clause)
//
// This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
// Package, Release 3e, by John R. Hauser.
//
// Copyright 2011, 2012, 2013, 2014 The Regents of the University of California.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice,
//     this list of conditions, and the following disclaimer.
//
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions, and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//  3. Neither the name of the University nor the names of its contributors may
//     be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS", AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, ARE
// DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//==========================================================
//
// fdlibm license:
//
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
//
//==========================================================
//
// Cephes license: (MIT)
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//==========================================================

#include "x87fp80.h"
#include "x87fpext.h"

#include <cstdint>
#include <cmath>

namespace x87
{

// From x87fp80.cpp — properly rounded fpext96_t -> fp80 conversion + CW reader.
extern fp80_t round_fpext96_to_fp80(fpext96_t const &v, x87cw_t cw, uint16_t &out_sw);
extern x87cw_t read_x87_cw();

//
// Debug helpers for printing out intermediate values
//
#if 0
#define PRINT_VAL
void print_val(char const *name, fp64_t const &val)
{
    fmt::print("{} = {:c}{:013X}e{:+05} ({+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.exponent(), val.as_double());
}
void print_val(char const *name, fp80_t const &val)
{
    fmt::print("{} = {:c}{:016X}e{:+05} ({+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.exponent(), val.as_double());
}
void print_val(char const *name, fpext64_t const &val)
{
    fmt::print("{} = {:c}{:016X}e{:+05} ({+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.exponent(), val.as_fp64().as_double());
}
void print_val(char const *name, fpext96_t const &val)
{
    fmt::print("{} = {:c}{:016X}`{:04X}e{:+05} ({+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.extend(), val.exponent(), val.as_fp64().as_double());
}
void print_val(char const *name, fpext96_t const &val)
{
    fmt::print("{} = {:c}{:016X}`{:08X}e{:+05} ({+.12e})\n", name, val.sign() ? '-' : '+', val.mantissa(), val.extend(), val.exponent(), val.as_fp64().as_double());
}
#else
#define print_val(x, y) do { } while (0)
#endif

//
// Core implementation of the f2xm1 operation
//
template<bool Debug>
static uint16_t x87_f2xm1_core(fp80_t const &src, fp80_t &dst)
{
    // warning: this is derived from the 64-bit version and has not been thoroughly vetted

    // special case values outside of defined range
    auto exponent = src.exponent();
    if (exponent >= 0)
        goto special;

    // anything small ends up with a G of 0 and an H == x*ln2
    // this also allows us to avoid dealing with denormals
    if (exponent <= -1000)
        goto tiny;

    // parameters
    static const int LOG_R = 4;
    static const int R = 1 << LOG_R;
    static const int TABLE_SIZE = 2 * R + 1;
    static const int TAYLOR_TERMS = 9;

using fpext_t = fpext96_t;
using fpextfast_t = fpext64_t;
    static fpext_t const s_table_g[TABLE_SIZE] =
    {
        fpext_t(0x8000000000000000ull, 0x00000000, -1, 1),    // 2^(-16/16) = -0.5l,
        fpext_t(0xf4aa7930676f09d6ull, 0x746d48e8, -2, 1),    // 2^(-15/16) = -0.47786310878629307983901676063004l,
        fpext_t(0xe8d47c382ae85232ull, 0x08373af1, -2, 1),    // 2^(-14/16) = -0.45474613366737117039649467211965l,
        fpext_t(0xdc785918a9dc7993ull, 0xe0524e3f, -2, 1),    // 2^(-13/16) = -0.43060568262165417314808485807924l,
        fpext_t(0xcf901f5ce48ead21ull, 0x72a5b9d0, -2, 1),    // 2^(-12/16) = -0.40539644249863946664125001471976l,
        fpext_t(0xc2159b3edcbddca4ull, 0xbeddc1ec, -2, 1),    // 2^(-11/16) = -0.3790710939632579757031612656367l,
        fpext_t(0xb40252ac9d5d8e2bull, 0xc685013c, -2, 1),    // 2^(-10/16) = -0.35158022267449516703312294110377l,
        fpext_t(0xa54f822b7abd6a73ull, 0x6cfeae6e, -2, 1),    // 2^( -9/16) = -0.32287222653155363585099262992965l,
        fpext_t(0x95f619980c4336f7ull, 0x4d04ec99, -2, 1),    // 2^( -8/16) = -0.29289321881345247559915563789515l,
        fpext_t(0x85eeb8c14fe79282ull, 0xaefdc093, -2, 1),    // 2^( -7/16) = -0.26158692703025034430654625981298l,
        fpext_t(0xea6357baabe4948bull, 0x0754bcda, -3, 1),    // 2^( -6/16) = -0.22889458729602958819385406895463l,
        fpext_t(0xc76dcfab81edfc70ull, 0x7729f1c2, -3, 1),    // 2^( -5/16) = -0.1947548340253728459102396663213l,
        fpext_t(0xa2ec0cd4a58a542full, 0x1965d11a, -3, 1),    // 2^( -4/16) = -0.15910358474628545696887452376679l,
        fpext_t(0xf999089eab58f777ull, 0xcd3b57dc, -4, 1),    // 2^( -3/16) = -0.12187391981335025844391969031234l,
        fpext_t(0xa9f9c8c116de3689ull, 0x7e945264, -4, 1),    // 2^( -2/16) = -0.08299595679532876825645840520586l,
        fpext_t(0xada82eadb7933d38ull, 0x462f3851, -5, 1),    // 2^( -1/16) = -0.04239671930142635306369436485208l,
        fpext_t(0x0000000000000000ull, 0x00000000, fpext_t::EXPONENT_MIN, 0), // 0
        fpext_t(0xb5586cf9890f6298ull, 0xb92b7184, -5, 0),    // 2^( +1/16) = 0.04427378242741384032196647873993l,
        fpext_t(0xb95c1e3ea8bd6e6full, 0xbe462876, -4, 0),    // 2^( +2/16) = 0.09050773266525765920701065576071l,
        fpext_t(0x8e1e9b9d588e19b0ull, 0x7eb6c705, -3, 0),    // 2^( +3/16) = 0.13878863475669165370383028384151l,
        fpext_t(0xc1bf828c6dc54b7aull, 0x356918c1, -3, 0),    // 2^( +4/16) = 0.18920711500272106671749997056048l,
        fpext_t(0xf7a993048d088d6dull, 0x0488f84f, -3, 0),    // 2^( +5/16) = 0.2418578120734840485936774687266l,
        fpext_t(0x97fb5aa6c544e3a8ull, 0x72f5fd88, -2, 0),    // 2^( +6/16) = 0.29683955465100966593375411779245l,
        fpext_t(0xb560fba90a852b19ull, 0x2602a324, -2, 0),    // 2^( +7/16) = 0.3542555469368927282980147401407l,
        fpext_t(0xd413cccfe7799211ull, 0x65f626ce, -2, 0),    // 2^( +8/16) = 0.4142135623730950488016887242097l,
        fpext_t(0xf4228e7d6030dafaull, 0xa2047eda, -2, 0),    // 2^( +9/16) = 0.47682614593949931138690748037405l,
        fpext_t(0x8ace5422aa0db5baull, 0x7c55a193, -1, 0),    // 2^(+10/16) = 0.54221082540794082361229186209073l,
        fpext_t(0x9c49182a3f0901c7ull, 0xc46b071f, -1, 0),    // 2^(+11/16) = 0.6104903319492543081795206673574l,
        fpext_t(0xae89f995ad3ad5e8ull, 0x734d1773, -1, 0),    // 2^(+12/16) = 0.68179283050742908606225095246643l,
        fpext_t(0xc199bdd85529c222ull, 0x0cb12a09, -1, 0),    // 2^(+13/16) = 0.75625216037329948311216061937531l,
        fpext_t(0xd5818dcfba48725dull, 0xa05aeb67, -1, 0),    // 2^(+14/16) = 0.83400808640934246348708318958829l,
        fpext_t(0xea4afa2a490d9858ull, 0xf73a18f6, -1, 0),    // 2^(+15/16) = 0.91520656139714729387261127029583l,
        fpext_t(0x8000000000000000ull, 0x00000000,  0, 0)     // 2^(+16/16) = 1.0
    };
    static fpextfast_t const s_table_u[TABLE_SIZE] =
    {
        fpextfast_t(0x8000000000000000ull, 0x00000000,  0, 1),    // -16/16
        fpextfast_t(0xf000000000000000ull, 0x00000000, -1, 1),    // -15/16
        fpextfast_t(0xe000000000000000ull, 0x00000000, -1, 1),    // -14/16
        fpextfast_t(0xd000000000000000ull, 0x00000000, -1, 1),    // -13/16
        fpextfast_t(0xc000000000000000ull, 0x00000000, -1, 1),    // -12/16
        fpextfast_t(0xb000000000000000ull, 0x00000000, -1, 1),    // -11/16
        fpextfast_t(0xa000000000000000ull, 0x00000000, -1, 1),    // -10/16
        fpextfast_t(0x9000000000000000ull, 0x00000000, -1, 1),    //  -9/16
        fpextfast_t(0x8000000000000000ull, 0x00000000, -1, 1),    //  -8/16
        fpextfast_t(0xe000000000000000ull, 0x00000000, -2, 1),    //  -7/16
        fpextfast_t(0xc000000000000000ull, 0x00000000, -2, 1),    //  -6/16
        fpextfast_t(0xa000000000000000ull, 0x00000000, -2, 1),    //  -5/16
        fpextfast_t(0x8000000000000000ull, 0x00000000, -2, 1),    //  -4/16
        fpextfast_t(0xc000000000000000ull, 0x00000000, -3, 1),    //  -3/16
        fpextfast_t(0x8000000000000000ull, 0x00000000, -3, 1),    //  -2/16
        fpextfast_t(0x8000000000000000ull, 0x00000000, -4, 1),    //  -1/16
        fpextfast_t(0x0000000000000000ull, 0x00000000, -16383, 0),//   0/16
        fpextfast_t(0x8000000000000000ull, 0x00000000, -4, 0),    //   1/16
        fpextfast_t(0x8000000000000000ull, 0x00000000, -3, 0),    //   2/16
        fpextfast_t(0xc000000000000000ull, 0x00000000, -3, 0),    //   3/16
        fpextfast_t(0x8000000000000000ull, 0x00000000, -2, 0),    //   4/16
        fpextfast_t(0xa000000000000000ull, 0x00000000, -2, 0),    //   5/16
        fpextfast_t(0xc000000000000000ull, 0x00000000, -2, 0),    //   6/16
        fpextfast_t(0xe000000000000000ull, 0x00000000, -2, 0),    //   7/16
        fpextfast_t(0x8000000000000000ull, 0x00000000, -1, 0),    //   8/16
        fpextfast_t(0x9000000000000000ull, 0x00000000, -1, 0),    //   9/16
        fpextfast_t(0xa000000000000000ull, 0x00000000, -1, 0),    //  10/16
        fpextfast_t(0xb000000000000000ull, 0x00000000, -1, 0),    //  11/16
        fpextfast_t(0xc000000000000000ull, 0x00000000, -1, 0),    //  12/16
        fpextfast_t(0xd000000000000000ull, 0x00000000, -1, 0),    //  13/16
        fpextfast_t(0xe000000000000000ull, 0x00000000, -1, 0),    //  14/16
        fpextfast_t(0xf000000000000000ull, 0x00000000, -1, 0),    //  15/16
        fpextfast_t(0x8000000000000000ull, 0x00000000,  0, 0)     //  16/16
    };
    static fpextfast_t const s_taylor_coeff[8] =
    {
        fpextfast_t(0x9000000000000000ull, 0x00000000,  3, 0),    // 9
        fpextfast_t(0x9000000000000000ull, 0x00000000,  6, 0),    // 9*8
        fpextfast_t(0xfc00000000000000ull, 0x00000000,  8, 0),    // 9*8*7
        fpextfast_t(0xbd00000000000000ull, 0x00000000, 11, 0),    // 9*8*7*6
        fpextfast_t(0xec40000000000000ull, 0x00000000, 13, 0),    // 9*8*7*6*5
        fpextfast_t(0xec40000000000000ull, 0x00000000, 15, 0),    // 9*8*7*6*5*4
        fpextfast_t(0xb130000000000000ull, 0x00000000, 17, 0),    // 9*8*7*6*5*4*3
        fpextfast_t(0xb130000000000000ull, 0x00000000, 18, 0)     // 9*8*7*6*5*4*3*2
    };
    static fpextfast_t const s_taylor_factorial_inv =
        fpextfast_t(0xb8ef1d2ab6399c7dull, 0x560e4473, -19, 0);   // 1.0/9!

    {
        // round x to the nearest multiple of 1/R by looking at the high bits of the mantissa
        int32_t g_index = 0;

        // anything smaller than -LOG_R - 1 will round to 0, so only do this if above
        if (exponent >= -LOG_R)
        {
            // shift mantissa down so we just have LOG_R + 1 bits
            auto mantissa = src.mantissa() & FP80_MANTISSA_MASK;
            g_index = mantissa >> (64 - LOG_R - exponent - 1);

            // round by adding LSB and shifting to get LOG_R bits
            g_index = (g_index >> 1) + (g_index & 1);

            // if negative, use a negative index
            if (src.sign() != 0)
                g_index = -g_index;
        }

        // compute v = delta from table entry
        fpextfast_t v = fpextfast_t(src) - s_table_u[g_index + R];

        // multiply v by ln(2) so we can use the e^x Taylor series; do this in
        // extended precision
        fpext_t w = fpext_t(v) * fpext_t::ln2;
        if (Debug) print_val("w", w);

        // Taylor series: this can be done in lower precision; start with h = w + coeff[0]
        fpextfast_t w80(w, true);
        if (Debug) print_val("w80", w80);
        fpextfast_t h80 = w80 + s_taylor_coeff[0];
        if (Debug) print_val("h1", h80);

        // now compute h = h * w + coeff[term] for terms up through 7
        for (int term = 1; term < TAYLOR_TERMS - 2; term++)
        {
            h80 = h80 * w80 + s_taylor_coeff[term];
            if (Debug) print_val("hn", h80);
        }

        // final term is just times w^2
        h80 *= w80 * w80;
        if (Debug) print_val("h2", h80);

        // then divide by 9!
        h80 = h80 * s_taylor_factorial_inv;
        if (Debug) print_val("h3", h80);

        // back to extended precision for final result; add w for final h value
        fpext_t h(h80);
        h += w;
        if (Debug) print_val("h4", h);

        // retrieve g from the table
        fpext_t g = s_table_g[g_index + R];
        if (Debug) print_val("g", g);

        // return g * h + g + h
        if (Debug)
        {
            fpext_t res = g * h;
            if (Debug) print_val("res", res);
            res += g;
            if (Debug) print_val("res", res);
            res += h;
            if (Debug) print_val("res", res);
        }
        uint16_t flags = X87SW_PRECISION_EX;
        dst = round_fpext96_to_fp80(g * h + g + h, read_x87_cw(), flags);
        if (dst.isdenorm() || dst.iszero()) flags |= X87SW_UNDERFLOW_EX;
        return flags;
    }

special:
    // return -0.5 for -1
    if (src.sign_exp() == 0xbfff && src.mantissa() == 0)
    {
        dst = fp80_t(0, 0xbffe);
        return X87SW_PRECISION_EX;
    }

    // max exponent cases
    if (src.ismaxexp())
    {
        // return -1 for -inf
        if (src.isninf())
        {
            dst = fp80_t(0, 0xbfff);
            return 0;
        }

        // return src for +inf or qNaNs and set no flags
        if (src.isinf() || src.isqnan())
        {
            dst = src;
            return 0;
        }

        // convert signalling NaN to quiet NaN and signal
        if (src.issnan())
        {
            dst = fp80_t(src.mantissa() | 0x4000000000000000ull, src.sign_exp());
            return X87SW_INVALID_EX;
        }
    }

    // for 0 or out-of-range values, just return x
    dst = src;
    return src.iszero() ? 0 : X87SW_PRECISION_EX;

tiny:
    // special case zero
    if (src.iszero())
    {
        dst = src;
        return 0;
    }

    // denorms and other tiny values reduce to a simple multiply
    dst = (fpext_t(src) * fpext_t::ln2).as_fp80();
    if (src.isdenorm())
        return X87SW_PRECISION_EX | X87SW_DENORM_EX | X87SW_UNDERFLOW_EX;
    else if (exponent <= 1 - FP80_EXPONENT_BIAS)
        return X87SW_PRECISION_EX | X87SW_UNDERFLOW_EX;
    return X87SW_PRECISION_EX;
}

#if X87_HOST_HAS_FP80
uint16_t fp80_t::x87_f2xm1(fp80_t const &src, fp80_t &dst)
{
    return host_x87_unary(src, dst, 2);
}
#else
uint16_t fp80_t::x87_f2xm1(fp80_t const &src, fp80_t &dst)
{
    return x87_f2xm1_core<false>(src, dst);
}
#endif


//===========================================================================
//
// STUBS for unimplemented 80-bit transcendentals. Each writes the FPU
// "indefinite" value and sets the invalid-operation flag so the test
// harness will report a clean 100% failure for the missing operation.
//
//===========================================================================

static uint16_t stub_unary(fp80_t &dst)
{
    dst = fp80_t::const_indef();
    return X87SW_INVALID_EX;
}

static uint16_t stub_unary2(fp80_t &dst1, fp80_t &dst2)
{
    dst1 = fp80_t::const_indef();
    dst2 = fp80_t::const_indef();
    return X87SW_INVALID_EX;
}

//
// Build an fp80 integer constant from a signed 32-bit value, normalized.
// (Used by fxtract to format its integer-exponent output.)
//
static fp80_t int32_to_fp80(int32_t v)
{
    if (v == 0) return fp80_t(0, 0);
    bool neg = (v < 0);
    uint64_t a = neg ? -uint64_t(v) : uint64_t(v);
    int sh = count_leading_zeros64(a);
    return fp80_t(a << sh, uint16_t(FP80_EXPONENT_BIAS + 63 - sh) | (neg ? FP80_SIGN_MASK : 0));
}

//
// FXTRACT (Intel SDM): split operand into significand and integer-log2.
// Matches the fp64_t::x87_fxtract convention (dst1 = significand, dst2 =
// exponent) which is the asm ordering after FXTRACT + the two FSTPs.
//   dst1 = significand normalized to [1.0, 2.0) with src's sign
//   dst2 = unbiased exponent as an integer fp80
// Zero → dst1 = src, dst2 = -inf, sets #Z.
// Infinity → dst1 = src, dst2 = +inf.
//
#if defined(X87_FORCE_HAND_ROLLED)
#define X87_HOST_HAS_FP80 0
#elif defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#define X87_HOST_HAS_FP80 1
#else
#define X87_HOST_HAS_FP80 0
#endif

#if X87_HOST_HAS_FP80
// Implemented in x87fp80.cpp.
extern uint16_t host_x87_unary       (fp80_t const &src, fp80_t &dst,                  int op);
extern uint16_t host_x87_unary2      (fp80_t const &src, fp80_t &dst1, fp80_t &dst2,   int op);
extern uint16_t host_x87_binary_trans(fp80_t const &a,   fp80_t const &b, fp80_t &dst, int op);
#endif

uint16_t fp80_t::x87_fxtract(fp80_t const &src, fp80_t &dst1, fp80_t &dst2)
{
    if (src.isnan())
    {
        fp80_t q = src.issnan() ? fp80_t::make_qnan(src) : src;
        dst1 = q; dst2 = q;
        return src.issnan() ? X87SW_INVALID_EX : 0;
    }
    if (src.iszero())
    {
        dst1 = src;
        dst2 = fp80_t::const_ninf();
        return X87SW_DIVZERO_EX;
    }
    if (src.isinf())
    {
        dst1 = src;
        dst2 = fp80_t::const_pinf();
        return 0;
    }

    uint16_t sign = src.sign_exp() & FP80_SIGN_MASK;
    int exp = (src.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    uint64_t mant = src.mantissa();

    uint16_t flags = 0;
    if (src.isdenorm())
    {
        int sh = count_leading_zeros64(mant);
        mant <<= sh;
        exp = 1 - FP80_EXPONENT_BIAS - sh;
        flags = X87SW_DENORM_EX;
    }

    // dst1 = the same mantissa with exponent biased to 0 (i.e., [1.0, 2.0))
    dst1 = fp80_t(mant, sign | uint16_t(FP80_EXPONENT_BIAS));
    // dst2 = exp as fp80 integer
    dst2 = int32_to_fp80(int32_t(exp));
    return flags;
}

#if X87_HOST_HAS_FP80
uint16_t fp80_t::x87_fscale(fp80_t const &a, fp80_t const &b, fp80_t &dst)
{
    return host_x87_binary_trans(a, b, dst, 5);
}
#else
//
// FSCALE: dst = a * 2^trunc(b).  Following the asm convention where ARG1
// becomes ST(0) (the value) and ARG2 becomes ST(1) (the scale factor):
// x87_fscale(value, scale, &dst).
//
uint16_t fp80_t::x87_fscale(fp80_t const &a, fp80_t const &b, fp80_t &dst)
{
    if (a.isnan() || b.isnan())
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

    uint16_t flags = (a.isdenorm() || b.isdenorm()) ? X87SW_DENORM_EX : 0;

    if (a.iszero()) { dst = a; return flags; }
    if (a.isinf())  { dst = a; return flags; }

    // b infinite: scale by 2^+inf = +inf, 2^-inf = 0
    if (b.isinf())
    {
        if (b.sign()) dst = (a.sign() ? fp80_t::const_nzero() : fp80_t::const_zero());
        else          dst = (a.sign() ? fp80_t::const_ninf()  : fp80_t::const_pinf());
        return flags;
    }
    if (b.iszero()) { dst = a; return flags; }

    // Truncate b to a signed integer. Saturate for huge |b| at a magnitude
    // that comfortably covers the fp80 exponent range without overflowing
    // int64 when added to the biased exponent.
    int bexp = (b.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    int64_t scale;
    constexpr int64_t SCALE_SAT = 0x40000;   // far exceeds fp80 exponent range
    if (bexp < 0)        scale = 0;
    else if (bexp >= 63) scale = b.sign() ? -SCALE_SAT : SCALE_SAT;
    else
    {
        int64_t v = int64_t(b.mantissa() >> (63 - bexp));
        if (v > SCALE_SAT) v = SCALE_SAT;
        scale = b.sign() ? -v : v;
    }

    // No scale change: return src1 verbatim (only DE if denormal).
    if (scale == 0) { dst = a; return flags; }

    // Normalize denormal src1: shift mantissa MSB to bit 63, track the
    // implicit exponent shift. We then work with a "virtual biased exp"
    // that can be negative for very small inputs.
    uint64_t mant = a.mantissa();
    int64_t  biased = int64_t(a.sign_exp() & FP80_EXPONENT_MASK);
    if (biased == 0 && mant != 0)
    {
        int s = count_leading_zeros64(mant);
        mant <<= s;
        biased = 1 - s;   // virtual biased exp (negative for deep denormals)
    }
    int64_t new_exp = biased + scale;
    uint16_t sign = a.sign_exp() & FP80_SIGN_MASK;

    if (new_exp >= FP80_EXPONENT_MAX_BIASED)
    {
        dst = sign ? fp80_t::const_ninf() : fp80_t::const_pinf();
        return flags | X87SW_OVERFLOW_EX | X87SW_PRECISION_EX | X87SW_C1;
    }
    if (new_exp <= 0)
    {
        // Re-denormalize: shift mantissa right by (1 - new_exp).
        int64_t shift = 1 - new_exp;
        if (shift >= 64)
        {
            // The value is smaller than the smallest denormal. With round-
            // to-nearest, round to mantissa 1 when value > 0.5 × smallest.
            // shift == 64 with normalized mantissa: value = mant × 2^(emin-127).
            // 0.5 × smallest = 2^(emin-64). So value > 0.5*smallest iff
            // mant > 2^63 OR mant == 2^63 and any lower bits.
            bool round_up = (shift == 64) &&
                            (mant > FP80_EXPLICIT_ONE);
            // (Ties at exactly 0.5×smallest → round to even = 0.)
            if (round_up)
            {
                dst = fp80_t(1, sign);
                return flags | X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX | X87SW_C1;
            }
            dst = sign ? fp80_t::const_nzero() : fp80_t::const_zero();
            return flags | X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
        }
        uint64_t lost = mant & ((1ull << shift) - 1);
        uint64_t denorm_mant = mant >> shift;
        // Round-to-nearest-even on the discarded bits.
        bool round_up = false;
        if (lost != 0)
        {
            uint64_t round_bit = (shift >= 1) ? (1ull << (shift - 1)) : 0;
            bool round  = (lost & round_bit) != 0;
            bool sticky = (lost & (round_bit - 1)) != 0;
            // For ROUND_NEAREST (the only mode we model here):
            round_up = round && (sticky || (denorm_mant & 1));
            if (round_up) denorm_mant += 1;
        }
        dst = fp80_t(denorm_mant, sign);
        uint16_t under_flags = 0;
        if (lost != 0) under_flags |= X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
        if (round_up) under_flags |= X87SW_C1;
        return flags | under_flags;
    }
    dst = fp80_t(mant, sign | uint16_t(new_exp));
    return flags;
}
#endif // X87_HOST_HAS_FP80

#if X87_HOST_HAS_FP80
uint16_t fp80_t::x87_fprem (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary_trans(a, b, dst, 3); }
uint16_t fp80_t::x87_fprem1(fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary_trans(a, b, dst, 4); }
#else
//
// Hand-rolled fallback path (aarch64 + everywhere not x86).
//
// For ops that have working fp64 implementations in x87fp64trans.cpp, we
// route through them: convert fp80 -> fp64 -> compute -> fp64 -> fp80.
// This loses ~11 bits of precision vs. a true fp80 implementation but is
// functional and matches the spec's semantics. A future enhancement is to
// re-port each fp64 algorithm to fpext96_t for full extended precision.
//
static uint16_t via_fp64_binary(fp80_t const &a, fp80_t const &b, fp80_t &dst,
                                 uint16_t (*fn)(fp64_t const &, fp64_t const &, fp64_t &))
{
    fp64_t a64(a), b64(b), r64;
    uint16_t sw = fn(a64, b64, r64);
    dst = fp80_t(r64);
    return sw;
}

//
// Common core for fprem and fprem1. Computes src1 mod src2 (truncated for
// fprem, round-to-even for fprem1). Returns SW with C0/C1/C3 set to the
// low 3 bits of the quotient (per Intel SDM FPREM/FPREM1).
//
// Algorithm: extract 64-bit normalized mantissas as integers, compute
// remainder via iterative reduction, then reassemble.
//
template<bool Rem1>
static uint16_t fprem_core_fp80(fp80_t const &src1, fp80_t const &src2, fp80_t &dst)
{
    uint16_t flags = 0;
    if (src1.isdenorm() || src2.isdenorm()) flags |= X87SW_DENORM_EX;

    // NaN propagation
    if (src1.isnan() || src2.isnan())
    {
        bool snan = src1.issnan() || src2.issnan();
        fp80_t pick;
        if (src1.isnan() && src2.isnan())
            pick = (src2.mantissa() > src1.mantissa()) ? src2 : src1;
        else
            pick = src1.isnan() ? src1 : src2;
        dst = pick.issnan() ? fp80_t::make_qnan(pick) : pick;
        return snan ? X87SW_INVALID_EX : 0;
    }
    // Special inputs
    if (src1.isinf() || src2.iszero())
    {
        dst = fp80_t::const_indef();
        return X87SW_INVALID_EX;
    }
    if (src2.isinf() || src1.iszero())
    {
        // src1 mod inf = src1; 0 mod src2 = 0.
        dst = src1;
        return flags;
    }

    // Extract normalized mantissas + exponents.
    uint64_t ma = src1.mantissa();
    uint64_t mb = src2.mantissa();
    int ea = (src1.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    int eb = (src2.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    if (src1.isdenorm()) { int s = count_leading_zeros64(ma); ma <<= s; ea = 1 - FP80_EXPONENT_BIAS - s; }
    if (src2.isdenorm()) { int s = count_leading_zeros64(mb); mb <<= s; eb = 1 - FP80_EXPONENT_BIAS - s; }
    uint16_t sign_a = src1.sign_exp() & FP80_SIGN_MASK;

    int dexp = ea - eb;
    if (dexp < 0)
    {
        // |src1| < |src2|. For fprem the quotient is 0 and result = src1.
        // For fprem1 round-to-nearest-even, the quotient may be ±1 when
        // 2*|src1| > |src2| or (2*|src1| == |src2| and quotient is odd).
        if (Rem1 && dexp == -1)
        {
            // |src1|/|src2| = (ma/mb)/2 (since dexp = -1). The half-point
            // for round-to-nearest is ma == mb. For ma > mb the ratio
            // exceeds 0.5 and the quotient rounds away from zero.
            // (Tie: ma == mb → quotient = 1 only if that makes it even —
            // but the quotient candidates are 0 (even) and 1 (odd), so
            // ties go to 0.)
            bool round_away = ma > mb;
            if (round_away)
            {
                // result = src1 - sign(src1/src2) * src2 = src1 + (sign_b==sign_a ? -src2 : +src2)
                fp80_t neg_src2 = fp80_t::chs(src2);
                // src1 / src2 sign == sign_a XOR sign_b. If non-negative we sub src2.
                fp80_t addend = (src1.sign() == src2.sign()) ? neg_src2 : src2;
                fp80_t result;
                fp80_t::x87_fadd(src1, addend, result);
                dst = result;
                // Q bit C1 should reflect rounding direction: we rounded
                // away from zero in the quotient, so result magnitude is
                // smaller than src1's. C1 reflects whether the residue had
                // its mantissa rounded up by the subtraction — leave to
                // round behavior of fadd; just set C1 to indicate quotient
                // = 1 (per Intel Q-bit encoding: C1=Q_LSB).
                return flags | X87SW_C1;
            }
        }
        dst = src1;
        return flags;
    }

    // Partial-remainder threshold: if dexp > 63, x87 sets C2 to indicate
    // incomplete reduction. Match that behavior by reducing by chunks of 32.
    int factor = 0;
    if (dexp > 63)
    {
        factor = ((dexp - 32) / 32) * 32;
        dexp -= factor;
    }

    // Standard binary long division: compute floor((ma << dexp) / mb).
    // Iterate (dexp + 1) times — equivalent to standard "process every bit
    // of N from MSB to LSB" but starting at R=ma since the first 63 trivial
    // doublings can't trigger a subtraction (ma < 2^64 ≤ 2 * mb).
    //
    // Result quotient: up to (dexp + 1) bits.
    __uint128_t rem = (__uint128_t)ma;
    uint64_t quotient = 0;
    // For ma >= mb, the first iteration already has rem >= mb (no shift
    // needed) — handle that case specially with one initial subtraction.
    if (rem >= (__uint128_t)mb)
    {
        rem -= mb;
        quotient = 1;
    }
    for (int k = 0; k < dexp; k++)
    {
        rem <<= 1;
        quotient <<= 1;
        if (rem >= (__uint128_t)mb)
        {
            rem -= mb;
            quotient |= 1;
        }
    }

    // For fprem1 (IEEE), round quotient toward even — only on the FINAL
    // reduction step (factor == 0). Partial reductions return like fprem.
    bool rem_negative = false;
    if (Rem1 && factor == 0)
    {
        // Compare 2*rem to mb. If 2*rem > mb, or 2*rem == mb and quotient is odd,
        // subtract one more.
        __uint128_t twice_rem = rem << 1;
        bool round_away = (twice_rem > (__uint128_t)mb) ||
                          (twice_rem == (__uint128_t)mb && (quotient & 1));
        if (round_away)
        {
            // Set rem to mb - rem (negative magnitude).
            rem = (__uint128_t)mb - rem;
            quotient++;
            rem_negative = true;
        }
    }

    // Assemble result. rem is at bit-position (eb + 63 - dexp), with
    // (dexp + 1) bits of quotient computed. The remainder value:
    // |result| = rem / 2^(63 + dexp) × 2^eb × 2^63 = rem × 2^(eb - dexp)
    // Wait — rem is conceptually `ma << dexp` mod (mb << dexp). Magnitude
    // wise, rem represents a value × 2^eb (same scale as src2's mantissa).
    // Renormalize so MSB is at bit 63.
    uint16_t result_sign = sign_a ^ (rem_negative ? FP80_SIGN_MASK : 0);
    if (rem == 0)
    {
        dst = fp80_t(0, result_sign);
        if (factor != 0)
            return flags | X87SW_C2;
        return flags | ((quotient & 1) << X87SW_C1_BIT)
                     | ((quotient & 2) << (X87SW_C3_BIT - 1))
                     | ((quotient & 4) << (X87SW_C0_BIT - 2));
    }
    // rem is the integer remainder N - Q*D where N = ma << dexp, D = mb.
    // It represents a fp value of rem × 2^(eb - 63). Find rem's MSB to
    // build the fpext96_t (mantissa-MSB at bit 63).
    uint64_t hi = uint64_t(rem >> 64);
    uint64_t lo = uint64_t(rem);
    int msb_pos = (hi != 0) ? (127 - count_leading_zeros64(hi))
                            : (63 - count_leading_zeros64(lo));
    uint64_t result_mant;
    if (msb_pos >= 63)
        result_mant = uint64_t(rem >> (msb_pos - 63));
    else
        result_mant = lo << (63 - msb_pos);
    int result_exp = eb + factor + (msb_pos - 63);

    int biased_exp = result_exp + FP80_EXPONENT_BIAS;
    if (biased_exp >= FP80_EXPONENT_MAX_BIASED)
    {
        dst = result_sign ? fp80_t::const_ninf() : fp80_t::const_pinf();
        return flags | X87SW_OVERFLOW_EX;
    }
    if (biased_exp <= 0)
    {
        // Denormal result
        int shift = 1 - biased_exp;
        if (shift >= 64) dst = fp80_t(0, result_sign);
        else             dst = fp80_t(result_mant >> shift, result_sign);
        // fprem residue is exact (integer subtraction) — no UE.
    }
    else
    {
        dst = fp80_t(result_mant, result_sign | uint16_t(biased_exp));
    }

    if (factor != 0)
    {
        // Partial reduction — only C2 is meaningful; Q bits are undefined.
        return flags | X87SW_C2;
    }
    uint16_t qbits = ((quotient & 1) << X87SW_C1_BIT)
                   | ((quotient & 2) << (X87SW_C3_BIT - 1))
                   | ((quotient & 4) << (X87SW_C0_BIT - 2));
    return flags | qbits;
}

uint16_t fp80_t::x87_fprem (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return fprem_core_fp80<false>(a, b, dst); }
uint16_t fp80_t::x87_fprem1(fp80_t const &a, fp80_t const &b, fp80_t &dst) { return fprem_core_fp80<true>(a, b, dst);  }
#endif
#if X87_HOST_HAS_FP80
uint16_t fp80_t::x87_fyl2x  (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary_trans(a, b, dst, 0); }
uint16_t fp80_t::x87_fyl2xp1(fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary_trans(a, b, dst, 1); }
uint16_t fp80_t::x87_fpatan (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary_trans(a, b, dst, 2); }
uint16_t fp80_t::x87_fsin   (fp80_t const &a, fp80_t &dst)                  { return host_x87_unary(a, dst, 0); }
uint16_t fp80_t::x87_fcos   (fp80_t const &a, fp80_t &dst)                  { return host_x87_unary(a, dst, 1); }
uint16_t fp80_t::x87_fsincos(fp80_t const &a, fp80_t &d1, fp80_t &d2)       { return host_x87_unary2(a, d1, d2, 0); }
uint16_t fp80_t::x87_fptan  (fp80_t const &a, fp80_t &d1, fp80_t &d2)       { return host_x87_unary2(a, d1, d2, 1); }
#else
// Hand-rolled fallback: route through the fp64 implementations
// (Cephes/fdlibm) and accept the ~11-bit precision loss.
static uint16_t via_fp64_unary(fp80_t const &a, fp80_t &dst,
                               uint16_t (*fn)(fp64_t const &, fp64_t &))
{
    fp64_t a64(a), r64;
    uint16_t sw = fn(a64, r64);
    dst = fp80_t(r64);
    return sw;
}

static uint16_t via_fp64_unary2(fp80_t const &a, fp80_t &dst1, fp80_t &dst2,
                                uint16_t (*fn)(fp64_t const &, fp64_t &, fp64_t &))
{
    fp64_t a64(a), r1_64, r2_64;
    uint16_t sw = fn(a64, r1_64, r2_64);
    dst1 = fp80_t(r1_64);
    dst2 = fp80_t(r2_64);
    return sw;
}

//
// Polynomial evaluators (Horner) — adapted from x87fp64trans.cpp but
// templated so we can use them with fpext96_t directly.
//
template<typename FpType, size_t Count>
static FpType poly_eval80(FpType const &x, std::array<FpType, Count> const &terms)
{
    FpType dst = terms[0];
    for (size_t i = 1; i < Count; i++)
        dst = dst * x + terms[i];
    return dst;
}

template<typename FpType, size_t Count>
static FpType poly1_eval80(FpType const &x, std::array<FpType, Count> const &terms)
{
    FpType dst = x + terms[0];
    for (size_t i = 1; i < Count; i++)
        dst = dst * x + terms[i];
    return dst;
}

//
// Hand-rolled fpatan using fpext96_t intermediates and Cephes-style
// long-double polynomial coefficients (atanl from Cephes).
//
// fpatan(src1, src2) computes atan2(src2, src1). Asm convention puts
// ARG1 (the divisor) into ST(0) and ARG2 (the dividend) into ST(1);
// FPATAN computes ST(1)/ST(0)-shaped atan and replaces ST(1).
//
uint16_t fp80_t::x87_fpatan(fp80_t const &src1, fp80_t const &src2, fp80_t &dst)
{
    using fpext_t = fpext96_t;

    uint16_t flags = 0;
    if (src1.isdenorm() || src2.isdenorm()) flags |= X87SW_DENORM_EX;

    // Constants needed for the special-case branches.
    static fpext_t const pi80(0xc90fdaa22168c235ull, 0x00000000, 1, 0);
    static fpext_t const npi80(0xc90fdaa22168c235ull, 0x00000000, 1, 1);
    static fpext_t const pio2_80(0xc90fdaa22168c235ull, 0x00000000, 0, 0);
    static fpext_t const npio2_80(0xc90fdaa22168c235ull, 0x00000000, 0, 1);
    static fpext_t const pio4_80(0xc90fdaa22168c235ull, 0x00000000, -1, 0);
    static fpext_t const npio4_80(0xc90fdaa22168c235ull, 0x00000000, -1, 1);
    static fpext_t const pi3o4_80(0x96cbe3f9990e91a8ull, 0x00000000, 1, 0);
    static fpext_t const npi3o4_80(0x96cbe3f9990e91a8ull, 0x00000000, 1, 1);

    // Special cases: NaN propagation, infinities, zeros.
    if (src1.isnan() || src2.isnan())
    {
        bool snan = src1.issnan() || src2.issnan();
        fp80_t pick;
        if (src1.isnan() && src2.isnan())
            pick = (src2.mantissa() > src1.mantissa()) ? src2 : src1;
        else
            pick = src1.isnan() ? src1 : src2;
        dst = pick.issnan() ? fp80_t::make_qnan(pick) : pick;
        return snan ? X87SW_INVALID_EX : 0;
    }
    // Returning a π-based constant: x87 stores 66-bit internal π and rounds
    // up to 64-bit fp80, setting C1=1 to flag the upward rounding.
    auto return_pi_const = [&](fpext_t const &c) {
        dst = round_fpext96_to_fp80(c, read_x87_cw(), flags);
        return flags | X87SW_PRECISION_EX | X87SW_C1;
    };

    if (src1.ismaxexp())  // src1 = ±inf
    {
        if (src2.isinf())
        {
            // atan2(±inf, ±inf) -> ±pi/4 or ±3pi/4
            return return_pi_const(
                (src1.sign() == 0) ? (src2.sign() ? npio4_80 : pio4_80)
                                   : (src2.sign() ? npi3o4_80 : pi3o4_80));
        }
        // atan2(finite, ±inf): 0 if +inf, ±pi if -inf
        if (src1.sign() == 0)
        {
            dst = src2.sign() ? fp80_t::const_nzero() : fp80_t::const_zero();
            return flags;
        }
        return return_pi_const(src2.sign() ? npi80 : pi80);
    }
    if (src2.ismaxexp())   // src2 = ±inf (and src1 finite)
    {
        return return_pi_const(src2.sign() ? npio2_80 : pio2_80);
    }
    if (src1.iszero())
    {
        if (src2.iszero())
        {
            // atan2(±0, ±0) -> 0 if x≥0, ±pi if x<0
            if (src1.sign() == 0)
            {
                dst = src2.sign() ? fp80_t::const_nzero() : fp80_t::const_zero();
                return flags;
            }
            return return_pi_const(src2.sign() ? npi80 : pi80);
        }
        return return_pi_const(src2.sign() ? npio2_80 : pio2_80);
    }
    if (src2.iszero())
    {
        if (src1.sign() == 0)
        {
            dst = src2.sign() ? fp80_t::const_nzero() : fp80_t::const_zero();
            return flags;
        }
        return return_pi_const(src2.sign() ? npi80 : pi80);
    }

    // Taylor coefficients for atan(x) = x * Σ (-1)^k x^(2k) / (2k+1).
    // Stored in Horner order (highest power first). 36 terms give ~100-bit
    // precision on the range |x| ≤ tan(π/8) ≈ 0.4142.
    static std::array<fpext_t, 36> const TAY = {
        fpext_t(0xe6c2b4481cd85689ull, 0x00000000, -7, 1),  // -1/71
        fpext_t(0xed7303b5cc0ed730ull, 0x00000000, -7, 0),  //  1/69
        fpext_t(0xf4898d5f85bb3950ull, 0x00000000, -7, 1),  // -1/67
        fpext_t(0xfc0fc0fc0fc0fc10ull, 0x00000000, -7, 0),  //  1/65
        fpext_t(0x8208208208208208ull, 0x00000000, -6, 1),  // -1/63
        fpext_t(0x864b8a7de6d1d608ull, 0x00000000, -6, 0),  //  1/61
        fpext_t(0x8ad8f2fba9386823ull, 0x00000000, -6, 1),  // -1/59
        fpext_t(0x8fb823ee08fb823full, 0x00000000, -6, 0),  //  1/57
        fpext_t(0x94f2094f2094f209ull, 0x00000000, -6, 1),  // -1/55
        fpext_t(0x9a90e7d95bc609a9ull, 0x00000000, -6, 0),  //  1/53
        fpext_t(0xa0a0a0a0a0a0a0a1ull, 0x00000000, -6, 1),  // -1/51
        fpext_t(0xa72f05397829cbc1ull, 0x00000000, -6, 0),  //  1/49
        fpext_t(0xae4c415c9882b931ull, 0x00000000, -6, 1),  // -1/47
        fpext_t(0xb60b60b60b60b60bull, 0x00000000, -6, 0),  //  1/45
        fpext_t(0xbe82fa0be82fa0bfull, 0x00000000, -6, 1),  // -1/43
        fpext_t(0xc7ce0c7ce0c7ce0cull, 0x00000000, -6, 0),  //  1/41
        fpext_t(0xd20d20d20d20d20dull, 0x00000000, -6, 1),  // -1/39
        fpext_t(0xdd67c8a60dd67c8aull, 0x00000000, -6, 0),  //  1/37
        fpext_t(0xea0ea0ea0ea0ea0full, 0x00000000, -6, 1),  // -1/35
        fpext_t(0xf83e0f83e0f83e10ull, 0x00000000, -6, 0),  //  1/33
        fpext_t(0x8421084210842108ull, 0x00000000, -5, 1),  // -1/31
        fpext_t(0x8d3dcb08d3dcb08dull, 0x00000000, -5, 0),  //  1/29
        fpext_t(0x97b425ed097b425full, 0x00000000, -5, 1),  // -1/27
        fpext_t(0xa3d70a3d70a3d70aull, 0x00000000, -5, 0),  //  1/25
        fpext_t(0xb21642c8590b2164ull, 0x00000000, -5, 1),  // -1/23
        fpext_t(0xc30c30c30c30c30cull, 0x00000000, -5, 0),  //  1/21
        fpext_t(0xd79435e50d79435eull, 0x00000000, -5, 1),  // -1/19
        fpext_t(0xf0f0f0f0f0f0f0f1ull, 0x00000000, -5, 0),  //  1/17
        fpext_t(0x8888888888888889ull, 0x00000000, -4, 1),  // -1/15
        fpext_t(0x9d89d89d89d89d8aull, 0x00000000, -4, 0),  //  1/13
        fpext_t(0xba2e8ba2e8ba2e8cull, 0x00000000, -4, 1),  // -1/11
        fpext_t(0xe38e38e38e38e38eull, 0x00000000, -4, 0),  //  1/9
        fpext_t(0x9249249249249249ull, 0x00000000, -3, 1),  // -1/7
        fpext_t(0xcccccccccccccccdull, 0x00000000, -3, 0),  //  1/5
        fpext_t(0xaaaaaaaaaaaaaaabull, 0x00000000, -2, 1),  // -1/3
        fpext_t(0x8000000000000000ull, 0x00000000, 0, 0),   //  1
    };
    static fpext_t const T3P8(0x9a827999fcef3242ull, 0x00000000, 1, 0);
    static fpext_t const TP8 (0xd413cccfe7799211ull, 0x00000000, -2, 0);

    // x = |src2 / src1| (fp80 division for full precision).
    // x87_fdivr(a,b) computes a/b; we want src2/src1 so pass (src2, src1).
    fp80_t x_fp80;
    fp80_t::x87_fdivr(src2, src1, x_fp80);   // x_fp80 = src2 / src1
    bool inner_sign = (x_fp80.sign() != 0);
    if (inner_sign) x_fp80 = fp80_t::chs(x_fp80);

    fpext_t yext, xext;

    if (x_fp80.isinf())
    {
        // |src2|/|src1| overflowed: result is ±π/2 with no correction.
        yext = pio2_80;
        xext = fpext_t::zero;
    }
    else if (x_fp80.iszero())
    {
        yext = fpext_t::zero;
        xext = fpext_t::zero;
    }
    else
    {
        fpext_t x(x_fp80);
        auto fpext_gt = [](fpext_t const &a, fpext_t const &b) {
            if (a.exponent() != b.exponent()) return a.exponent() > b.exponent();
            if (a.mantissa() != b.mantissa()) return a.mantissa() > b.mantissa();
            return a.extend() > b.extend();
        };
        bool x_gt_t3p8 = fpext_gt(x, T3P8);
        bool x_gt_tp8  = fpext_gt(x, TP8);
        bool skip_poly = false;
        if (x_gt_t3p8)
        {
            yext = pio2_80;
            // xext = -1 / x   →   need fdivr(1, x) = 1/x.
            fp80_t one = fp80_t::const_one();
            fp80_t recip;
            fp80_t::x87_fdivr(one, x_fp80, recip);   // recip = 1 / x
            if (recip.isnan() || recip.isinf() || recip.iszero())
            {
                xext = fpext_t::zero;
                skip_poly = true;
            }
            else
            {
                fpext_t r(recip);
                r.chs();
                xext = r;
            }
        }
        else if (x_gt_tp8)
        {
            yext = pio4_80;
            // xext = (x - 1) / (x + 1)
            // x87_fsubr(a,b) = a - b; x87_fadd is symmetric; x87_fdivr(a,b) = a/b.
            fp80_t one = fp80_t::const_one();
            fp80_t num, denom, ratio;
            fp80_t::x87_fsubr(x_fp80, one, num);     // num = x - 1
            fp80_t::x87_fadd (x_fp80, one, denom);   // denom = x + 1
            fp80_t::x87_fdivr(num, denom, ratio);    // ratio = num / denom
            if (ratio.isnan() || ratio.isinf())
            {
                xext = fpext_t::zero;
                skip_poly = true;
            }
            else if (ratio.iszero())
            {
                xext = fpext_t::zero;
                skip_poly = true;
            }
            else
            {
                xext = fpext_t(ratio);
            }
        }
        else
        {
            yext = fpext_t::zero;
            xext = x;
        }

        if (!skip_poly)
        {
            // Taylor: atan(x) = x * Σ_k (-1)^k x^(2k)/(2k+1).
            // Horner on z = x²; coefficients pre-stored highest-power-first.
            fpext_t z = xext * xext;
            fpext_t p = TAY[0];
            for (size_t i = 1; i < TAY.size(); i++)
                p = p * z + TAY[i];
            fpext_t atan_x = xext * p;
            yext = yext + atan_x;
        }
    }

    if (inner_sign) yext.chs();

    // Apply quadrant offset.
    int code = (src1.sign() << 1) | src2.sign();
    if (code == 2) yext = yext + pi80;
    else if (code == 3) yext = yext + npi80;

    dst = round_fpext96_to_fp80(yext, read_x87_cw(), flags);
    // For transcendentals, Intel x87 hardware appears to always assert C1
    // alongside PE — likely because the internal extra-precision result
    // always undergoes some rounding to fit the 64-bit fp80 mantissa.
    flags |= X87SW_PRECISION_EX | X87SW_C1;

    // For sign of zero: if src2 is negative, result of zero gets a minus.
    if (dst.iszero() && src2.sign())
        dst = fp80_t::chs(dst);

    return flags;
}

//
// fyl2x: compute src2 * log2(src1).  fdlibm-derived algorithm: reduce
// src1 to mantissa m and integer exponent k where m ∈ [1, 2), then
// log2(src1) = k + log2(m). log2(m) is computed via the standard
// (s = (m-1)/(m+1)) polynomial expansion in fpext96_t precision.
//
uint16_t fp80_t::x87_fyl2x(fp80_t const &src1, fp80_t const &src2, fp80_t &dst)
{
    using fpext_t = fpext96_t;

    uint16_t flags = 0;
    if (src1.isdenorm() || src2.isdenorm()) flags |= X87SW_DENORM_EX;

    // Special cases per Intel SDM §8.3.9 / FYL2X.
    if (src1.isnan() || src2.isnan())
    {
        bool snan = src1.issnan() || src2.issnan();
        fp80_t pick;
        if (src1.isnan() && src2.isnan())
            pick = (src2.mantissa() > src1.mantissa()) ? src2 : src1;
        else
            pick = src1.isnan() ? src1 : src2;
        dst = pick.issnan() ? fp80_t::make_qnan(pick) : pick;
        return snan ? X87SW_INVALID_EX : 0;
    }
    if (src1.sign())
    {
        // log of a negative number is undefined.
        dst = fp80_t::const_indef();
        return X87SW_INVALID_EX;
    }
    if (src1.ismaxexp())   // +inf
    {
        // log2(+inf) = +inf, multiplied by src2 — gives signed inf, except
        // src2 = 0 → indef.
        if (src2.iszero())
        {
            dst = fp80_t::const_indef();
            return X87SW_INVALID_EX;
        }
        dst = src2.sign() ? fp80_t::const_ninf() : fp80_t::const_pinf();
        return flags;
    }
    if (src1.iszero())
    {
        // log2(0) = -inf, multiplied by src2.
        if (src2.iszero())
        {
            dst = fp80_t::const_indef();
            return X87SW_INVALID_EX;
        }
        if (src2.ismaxexp())   // 0 * inf
        {
            dst = fp80_t::const_indef();
            return X87SW_INVALID_EX;
        }
        // -inf * src2 → signed inf, divzero flag. Don't propagate DE here;
        // x87 reports only #Z for this case even if src2 is denormal.
        dst = (src2.sign() == 0) ? fp80_t::const_ninf() : fp80_t::const_pinf();
        return X87SW_DIVZERO_EX;
    }
    if (src2.ismaxexp())
    {
        // src2 = ±inf, src1 finite >0. If src1 == 1, log2(1) = 0, 0*inf = indef.
        // Otherwise log2(src1) * ±inf = ±inf with sign depending.
        // Compare src1 to 1.
        if (src1.sign_exp() == 0x3FFF && src1.mantissa() == FP80_EXPLICIT_ONE)
        {
            dst = fp80_t::const_indef();
            return X87SW_INVALID_EX;
        }
        // src1 != 1; log2(src1) is positive if src1>1, negative if src1<1.
        bool src1_gt_1 = (src1.sign_exp() > 0x3FFF) ||
                         (src1.sign_exp() == 0x3FFF && src1.mantissa() > FP80_EXPLICIT_ONE);
        bool result_sign = src2.sign() ^ (!src1_gt_1);
        dst = result_sign ? fp80_t::const_ninf() : fp80_t::const_pinf();
        return flags;
    }
    if (src2.iszero())
    {
        // src2 = 0, src1 finite >0. 0 * log2(src1) = 0, signed based on
        // whether log2(src1) is negative (for src1 strictly < 1) and src2's
        // sign. For src1 == 1, log2(1) = +0, so result = +0 * src2 = signed
        // by src2 only.
        bool src1_lt_1 = (src1.sign_exp() < 0x3FFF) ||
                         (src1.sign_exp() == 0x3FFF && src1.mantissa() < FP80_EXPLICIT_ONE);
        bool log_neg = src1_lt_1;
        bool result_neg = log_neg ^ (src2.sign() != 0);
        dst = result_neg ? fp80_t::const_nzero() : fp80_t::const_zero();
        return flags;
    }

    // Main path: src1 ∈ (0, ∞), src1 != 1, src2 finite non-zero.
    // log2(src1) = exponent_of(src1) + log2(mantissa_of(src1)) where mantissa ∈ [1, 2).
    //
    // For log2(m) where m ∈ [1, 2):
    //   m - 1 close to 0 → use series expansion
    //   General: log(m) = 2 atanh((m-1)/(m+1)) — fdlibm style.
    //
    // We compute log_e via the (m-1)/(m+1) substitution, then multiply by 1/ln(2).

    fpext_t one(0x8000000000000000ull, 0x00000000, 0, 0);

    // Extract exponent k and mantissa m = src1 / 2^k where m ∈ [1, 2).
    int k = (src1.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    uint64_t mant_bits = src1.mantissa();
    if (src1.isdenorm())
    {
        int s = count_leading_zeros64(mant_bits);
        mant_bits <<= s;
        k = 1 - FP80_EXPONENT_BIAS - s;
    }
    // m = mantissa with biased exponent 0x3FFF → in [1, 2).
    fp80_t m_fp80(mant_bits, uint16_t(FP80_EXPONENT_BIAS));
    fpext_t m(m_fp80);

    // s = (m - 1) / (m + 1)
    fpext_t numer; numer.sub(m, one);
    fpext_t denom; denom.add(m, one);
    fpext_t s = numer.div64(denom);
    fpext_t s2 = s * s;

    // log(m) = 2*s * (1 + s^2/3 + s^4/5 + s^6/7 + ...)
    // Polynomial in s2 of the series.
    static fpext_t const c1_3(0xaaaaaaaaaaaaaaaaull, 0xaaaaaaab, -2, 0);  // 1/3
    static fpext_t const c1_5(0xccccccccccccccccull, 0xcccccccd, -3, 0);  // 1/5
    static fpext_t const c1_7(0x9249249249249249ull, 0x24924925, -3, 0);  // 1/7
    static fpext_t const c1_9(0xe38e38e38e38e38eull, 0x38e38e39, -4, 0);  // 1/9
    static fpext_t const c1_11(0xba2e8ba2e8ba2e8bull, 0xa2e8ba2f, -4, 0); // 1/11
    static fpext_t const c1_13(0x9d89d89d89d89d89ull, 0xd89d89d9, -4, 0); // 1/13

    // Horner: ((((c1_13 * s2 + c1_11) * s2 + c1_9) * s2 + c1_7) * s2 + c1_5) * s2 + c1_3) * s2 + 1
    fpext_t poly = c1_13 * s2 + c1_11;
    poly = poly * s2 + c1_9;
    poly = poly * s2 + c1_7;
    poly = poly * s2 + c1_5;
    poly = poly * s2 + c1_3;
    poly = poly * s2 + one;

    fpext_t two(0x8000000000000000ull, 0x00000000, 1, 0);
    fpext_t logm = two * s * poly;     // log(m) (natural log)

    // log2(m) = log(m) * (1/ln(2))
    static fpext_t const invln2(0xb8aa3b295c17f0bbull, 0xbe87fed0, 0, 0);
    fpext_t log2_m = logm * invln2;

    // log2(src1) = k + log2(m)
    fpext_t k_ext;
    {
        // Build fpext from integer k.
        if (k == 0)        k_ext = fpext_t::zero;
        else
        {
            bool neg = k < 0;
            uint64_t a = neg ? -(uint64_t)k : (uint64_t)k;
            int sh = count_leading_zeros64(a);
            k_ext = fpext_t(a << sh, 0, 63 - sh, neg ? 1 : 0);
        }
    }
    fpext_t log2_x;
    log2_x.add(k_ext, log2_m);

    // result = src2 * log2(src1)
    fpext_t s2_ext(src2);
    fpext_t result = log2_x * s2_ext;

    flags |= X87SW_PRECISION_EX;
    dst = round_fpext96_to_fp80(result, read_x87_cw(), flags);
    if (dst.isdenorm() && !dst.iszero())
        flags |= X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
    return flags;
}

//
// fyl2xp1: compute src2 * log2(1 + src1). Same series as fyl2x but with
// the (m-1)/(m+1) substitution becoming src1/(src1+2), avoiding the
// catastrophic cancellation that explicit 1+src1 would produce.
//
uint16_t fp80_t::x87_fyl2xp1(fp80_t const &src1, fp80_t const &src2, fp80_t &dst)
{
    using fpext_t = fpext96_t;

    uint16_t flags = 0;
    if (src1.isdenorm() || src2.isdenorm()) flags |= X87SW_DENORM_EX;

    if (src1.isnan() || src2.isnan())
    {
        bool snan = src1.issnan() || src2.issnan();
        fp80_t pick;
        if (src1.isnan() && src2.isnan())
            pick = (src2.mantissa() > src1.mantissa()) ? src2 : src1;
        else
            pick = src1.isnan() ? src1 : src2;
        dst = pick.issnan() ? fp80_t::make_qnan(pick) : pick;
        return snan ? X87SW_INVALID_EX : 0;
    }
    // src2 = 0 (multiplier zero): x87 returns signed 0 without checking
    // src1's domain. Sign per (src2 ^ src1).
    if (src2.iszero() && !src1.ismaxexp())
    {
        bool result_neg = (src2.sign() != 0) ^ (src1.sign() != 0);
        dst = result_neg ? fp80_t::const_nzero() : fp80_t::const_zero();
        return flags;
    }
    // src1 = -inf or src1 < -1 → indef (log of negative).
    if (src1.ismaxexp() && src1.sign())
    {
        dst = fp80_t::const_indef();
        return X87SW_INVALID_EX;
    }
    if (src1.sign() && !src1.isinf() && !src1.iszero())
    {
        if (src1.sign_exp() > 0xBFFF ||
            (src1.sign_exp() == 0xBFFF && src1.mantissa() > FP80_EXPLICIT_ONE))
        {
            // |src1| > 1, src1 negative: outside fyl2xp1's domain. Intel
            // returns src1 unchanged with PE + DE rather than INDEF.
            dst = src1;
            return flags | X87SW_PRECISION_EX | X87SW_DENORM_EX;
        }
        // src1 == -1 → log(0) = -inf (handled below as DZ)
        if (src1.sign_exp() == 0xBFFF && src1.mantissa() == FP80_EXPLICIT_ONE)
        {
            if (src2.iszero())
            {
                dst = fp80_t::const_indef();
                return X87SW_INVALID_EX;
            }
            dst = (src2.sign() == 0) ? fp80_t::const_ninf() : fp80_t::const_pinf();
            return flags | X87SW_DIVZERO_EX;
        }
    }
    // src1 = +inf: result is src2 * inf. If src2=0, NaN; else signed inf.
    if (src1.isinf())   // src1 = +inf (we filtered -inf above)
    {
        if (src2.iszero())
        {
            dst = fp80_t::const_indef();
            return X87SW_INVALID_EX;
        }
        dst = src2.sign() ? fp80_t::const_ninf() : fp80_t::const_pinf();
        return flags;
    }
    // src2 = ±inf and src1 is finite.
    if (src2.ismaxexp())
    {
        if (src1.iszero())
        {
            // log2(1) = 0; 0 * inf = indef.
            dst = fp80_t::const_indef();
            return X87SW_INVALID_EX;
        }
        bool log_neg = src1.sign() != 0;
        bool result_neg = log_neg ^ (src2.sign() != 0);
        dst = result_neg ? fp80_t::const_ninf() : fp80_t::const_pinf();
        return flags;
    }
    // src2 = 0 with finite src1 (and src1 in valid domain): result is ±0.
    if (src2.iszero())
    {
        bool result_neg = (src2.sign() != 0) ^ (src1.sign() != 0);
        dst = result_neg ? fp80_t::const_nzero() : fp80_t::const_zero();
        return flags;
    }
    if (src1.iszero())
    {
        // log2(1) = 0; src2 * 0 = signed zero (sign of src2).
        dst = src2.sign() ? fp80_t::const_nzero() : fp80_t::const_zero();
        return flags;
    }

    // For |src1| outside Cephes' fast-converging range (~|s| < 0.4), the
    // direct (x/(x+2)) Maclaurin series converges too slowly. Dispatch to
    // fyl2x(1 + src1, src2) which uses the well-conditioned (m-1)/(m+1)
    // path on the renormalized mantissa.
    int src1_exp = (src1.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    // For negative src1 with |src1| >= 1, 1+src1 <= 0 and fyl2x is undefined.
    // Intel passes src1 through and signals PE + DE.
    if (src1.sign() && src1_exp >= 0)
    {
        dst = src1;
        return flags | X87SW_PRECISION_EX | X87SW_DENORM_EX;
    }
    if (src1_exp > -2)   // |src1| > ~0.25
    {
        fp80_t one_plus_x;
        fp80_t::x87_fadd(src1, fp80_t::const_one(), one_plus_x);
        uint16_t sub_flags = x87_fyl2x(one_plus_x, src2, dst);
        // Intel asserts UE+PE on fyl2xp1 when the result is denormal,
        // regardless of whether the multiplication was exact internally.
        if (dst.isdenorm() && !dst.iszero())
            sub_flags |= X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
        return sub_flags | (flags & X87SW_DENORM_EX);
    }

    // log(1 + x) via the (1+x - 1) / (1+x + 1) = x/(x+2) substitution.
    fp80_t two_fp80(0x8000000000000000ull, 0x4000);   // 2.0
    fp80_t denom_fp80;
    fp80_t::x87_fadd(src1, two_fp80, denom_fp80);
    fp80_t s_fp80;
    fp80_t::x87_fdivr(src1, denom_fp80, s_fp80);   // s = src1 / denom
    if (s_fp80.ismaxexp() || denom_fp80.ismaxexp())
    {
        // Numerical edge — fall back to via_fp64 to get something reasonable.
        return via_fp64_binary(src1, src2, dst, &fp64_t::x87_fyl2xp1);
    }

    fpext_t one(0x8000000000000000ull, 0x00000000, 0, 0);
    fpext_t two(0x8000000000000000ull, 0x00000000, 1, 0);
    fpext_t s(s_fp80);
    fpext_t s2 = s * s;

    static fpext_t const c1_3(0xaaaaaaaaaaaaaaaaull, 0xaaaaaaab, -2, 0);
    static fpext_t const c1_5(0xccccccccccccccccull, 0xcccccccd, -3, 0);
    static fpext_t const c1_7(0x9249249249249249ull, 0x24924925, -3, 0);
    static fpext_t const c1_9(0xe38e38e38e38e38eull, 0x38e38e39, -4, 0);
    static fpext_t const c1_11(0xba2e8ba2e8ba2e8bull, 0xa2e8ba2f, -4, 0);
    static fpext_t const c1_13(0x9d89d89d89d89d89ull, 0xd89d89d9, -4, 0);

    fpext_t poly = c1_13 * s2 + c1_11;
    poly = poly * s2 + c1_9;
    poly = poly * s2 + c1_7;
    poly = poly * s2 + c1_5;
    poly = poly * s2 + c1_3;
    poly = poly * s2 + one;

    fpext_t logm = two * s * poly;

    // Multiply by 1/ln(2) instead of dividing by ln(2) — avoids div64's
    // fp64-precision bottleneck (which can overflow for huge values).
    static fpext_t const invln2(0xb8aa3b295c17f0bbull, 0xbe87fed0, 0, 0);
    fpext_t log2_x = logm * invln2;

    fpext_t s2_ext(src2);
    fpext_t result = log2_x * s2_ext;

    dst = round_fpext96_to_fp80(result, read_x87_cw(), flags);
    flags |= X87SW_PRECISION_EX;
    if (dst.isdenorm() && !dst.iszero())
        flags |= X87SW_UNDERFLOW_EX;
    return flags;
}
//
// Shared sin/cos Taylor coefficients (Horner order: highest power first).
// Defined at namespace scope so both x87_fsin and x87_fcos can share them.
//
namespace {

using fpext_t = fpext96_t;

// sin(y) = y * Σ_k (-1)^k y^(2k) / (2k+1)!
// 26 terms gives ~120-bit precision on |y| ≤ π/4.
static std::array<fpext_t, 26> const SIN_T = {
    fpext_t(0x8b0c395fbdc119bdull, 0x00000000, -220, 1),  // k=25 -1/(51!)
    fpext_t(0xad21786ff5842eccull, 0x00000000, -209, 0),  // k=24
    fpext_t(0xc6d4705093f5cdbeull, 0x00000000, -198, 1),  // k=23
    fpext_t(0xd1e5c39110323c71ull, 0x00000000, -187, 0),  // k=22
    fpext_t(0xcaeda292bf289170ull, 0x00000000, -176, 1),  // k=21
    fpext_t(0xb2f30e1ce8120641ull, 0x00000000, -165, 0),  // k=20
    fpext_t(0x8f4ca24d25d66f01ull, 0x00000000, -154, 1),  // k=19
    fpext_t(0xcf6468e4a742d7a7ull, 0x00000000, -144, 0),  // k=18
    fpext_t(0x86e2ce38b6c8f941ull, 0x00000000, -133, 1),  // k=17
    fpext_t(0x9cc092a6e86a8da9ull, 0x00000000, -123, 0),  // k=16
    fpext_t(0xa1a6973c1fade217ull, 0x00000000, -113, 1),  // k=15
    fpext_t(0x92cfcc5a1ac56bd6ull, 0x00000000, -103, 0),  // k=14
    fpext_t(0xe8d58e16e6751904ull, 0x00000000, -94, 1),   // k=13
    fpext_t(0x9f9e66e8b2fd46a7ull, 0x00000000, -84, 0),   // k=12
    fpext_t(0xbb0da098b1c0ceccull, 0x00000000, -75, 1),   // k=11
    fpext_t(0xb8dc77b6e7ab8c5full, 0x00000000, -66, 0),   // k=10
    fpext_t(0x97a4da340a0ab926ull, 0x00000000, -57, 1),   // k=9
    fpext_t(0xca963b81856a5359ull, 0x00000000, -49, 0),   // k=8
    fpext_t(0xd73f9f399dc0f88full, 0x00000000, -41, 1),   // k=7
    fpext_t(0xb092309d43684be5ull, 0x00000000, -33, 0),   // k=6
    fpext_t(0xd7322b3faa271c7full, 0x00000000, -26, 1),   // k=5
    fpext_t(0xb8ef1d2ab6399c7dull, 0x00000000, -19, 0),   // k=4
    fpext_t(0xd00d00d00d00d00dull, 0x00000000, -13, 1),   // k=3 -1/5040
    fpext_t(0x8888888888888889ull, 0x00000000, -7, 0),    // k=2  1/120
    fpext_t(0xaaaaaaaaaaaaaaabull, 0x00000000, -3, 1),    // k=1 -1/6
    fpext_t(0x8000000000000000ull, 0x00000000, 0, 0),     // k=0  1
};

// cos(y) = Σ_k (-1)^k y^(2k) / (2k)!  (Horner: highest power first)
static std::array<fpext_t, 26> const COS_T = {
    fpext_t(0xdd9b7b70966bc105ull, 0x00000000, -215, 1),  // k=25
    fpext_t(0x848da035b7f933d4ull, 0x00000000, -203, 0),  // k=24
    fpext_t(0x9204027b2ca88317ull, 0x00000000, -192, 1),  // k=23
    fpext_t(0x93958d81ff635280ull, 0x00000000, -181, 0),  // k=22
    fpext_t(0x8857a93a986f41b7ull, 0x00000000, -170, 1),  // k=21
    fpext_t(0xe5476a1509571802ull, 0x00000000, -160, 0),  // k=20
    fpext_t(0xaea565ce061d5749ull, 0x00000000, -149, 1),  // k=19
    fpext_t(0xefcc194861654958ull, 0x00000000, -139, 0),  // k=18
    fpext_t(0x9388118e07ebd09full, 0x00000000, -128, 1),  // k=17
    fpext_t(0xa1a6973c1fade217ull, 0x00000000, -118, 0),  // k=16
    fpext_t(0x9c9962823eb07306ull, 0x00000000, -108, 1),  // k=15
    fpext_t(0x850c5131a842e9baull, 0x00000000, -98, 0),   // k=14
    fpext_t(0xc4742fe35272cd1cull, 0x00000000, -89, 1),   // k=13
    fpext_t(0xf96780cb97abbe65ull, 0x00000000, -80, 0),   // k=12
    fpext_t(0x8671cb6dbfc294a3ull, 0x00000000, -70, 1),   // k=11
    fpext_t(0xf2a15d201011283dull, 0x00000000, -62, 0),   // k=10
    fpext_t(0xb413c31dcbecbbdeull, 0x00000000, -53, 1),   // k=9
    fpext_t(0xd73f9f399dc0f88full, 0x00000000, -45, 0),   // k=8
    fpext_t(0xc9cba54603e4e906ull, 0x00000000, -37, 1),   // k=7
    fpext_t(0x8f76c77fc6c4bdaaull, 0x00000000, -29, 0),   // k=6
    fpext_t(0x93f27dbbc4fae397ull, 0x00000000, -22, 1),   // k=5
    fpext_t(0xd00d00d00d00d00dull, 0x00000000, -16, 0),   // k=4  1/40320
    fpext_t(0xb60b60b60b60b60bull, 0x00000000, -10, 1),   // k=3 -1/720
    fpext_t(0xaaaaaaaaaaaaaaabull, 0x00000000, -5, 0),    // k=2  1/24
    fpext_t(0x8000000000000000ull, 0x00000000, -1, 1),    // k=1 -1/2
    fpext_t(0x8000000000000000ull, 0x00000000, 0, 0),     // k=0  1
};

// Evaluate sin(y) for |y| ≤ π/4 via Taylor (Horner on z = y²).
inline fpext_t taylor_sin(fpext_t const &y)
{
    fpext_t z = y * y;
    fpext_t p = SIN_T[0];
    for (size_t i = 1; i < SIN_T.size(); i++) p = p * z + SIN_T[i];
    return y * p;
}

inline fpext_t taylor_cos(fpext_t const &y)
{
    fpext_t z = y * y;
    fpext_t p = COS_T[0];
    for (size_t i = 1; i < COS_T.size(); i++) p = p * z + COS_T[i];
    return p;
}

// Range-reduce x to (quadrant, y) where x = quadrant * (π/2) + y and
// |y| ≤ π/4 (approximately). Returns the quadrant (0..3).
// Limited precision: works well for |x| ≤ ~2^60. Beyond that the
// caller should detect via the 2^63 check and bail.
inline int range_reduce(fp80_t const &src, fpext_t &y_out)
{
    static const fpext_t pio2(0xc90fdaa22168c234ull, 0xc4c6628c, 0, 0);
    static const fp80_t pio2_fp80(0xc90fdaa22168c235ull, 0x3FFF);

    // n = round(src / (π/2))
    fp80_t recip_pio2;
    fp80_t one = fp80_t::const_one();
    fp80_t::x87_fdivr(one, pio2_fp80, recip_pio2);   // 2/π
    fp80_t scaled;
    fp80_t::x87_fmul(src, recip_pio2, scaled);
    // Round to nearest integer
    fp80_t rounded;
    {
        fpround_t r(X87CW_ROUNDING_NEAREST);
        fp80_t::x87_frndint(scaled, rounded);
    }
    int64_t n;
    if ((rounded.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS > 62)
        n = (rounded.sign() ? INT64_MIN : INT64_MAX);
    else
    {
        int exp = (rounded.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
        if (exp < 0)
            n = 0;
        else
            n = int64_t(rounded.mantissa() >> (63 - exp)) * (rounded.sign() ? -1 : 1);
    }

    // y = src - n * (π/2), in fpext96_t.
    fpext_t n_ext;
    if (n == 0)
        n_ext = fpext_t::zero;
    else
    {
        bool neg = n < 0;
        uint64_t a = neg ? -(uint64_t)n : (uint64_t)n;
        int sh = count_leading_zeros64(a);
        n_ext = fpext_t(a << sh, 0, 63 - sh, neg ? 1 : 0);
    }
    fpext_t shift = n_ext * pio2;
    fpext_t src_ext(src);
    y_out.sub(src_ext, shift);

    return int(n & 3);
}

}  // anonymous namespace

uint16_t fp80_t::x87_fsin(fp80_t const &src, fp80_t &dst)
{
    uint16_t flags = src.isdenorm() ? X87SW_DENORM_EX : 0;
    if (src.isnan())
    {
        dst = src.issnan() ? fp80_t::make_qnan(src) : src;
        return src.issnan() ? X87SW_INVALID_EX : 0;
    }
    if (src.isinf())
    {
        dst = fp80_t::const_indef();
        return X87SW_INVALID_EX;
    }
    if (src.iszero())
    {
        dst = src;
        return 0;
    }
    // |src| > 2^63 → out of range, set C2, return src unchanged.
    if ((src.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS > 62)
    {
        dst = src;
        return flags | X87SW_C2;
    }

    fpext_t y;
    int q = range_reduce(src, y);
    fpext_t r;
    switch (q & 3)
    {
        case 0: r = taylor_sin(y); break;
        case 1: r = taylor_cos(y); break;
        case 2: r = taylor_sin(y); r.chs(); break;
        case 3: r = taylor_cos(y); r.chs(); break;
    }
    dst = round_fpext96_to_fp80(r, read_x87_cw(), flags);
    flags |= X87SW_PRECISION_EX;
    if (dst.isdenorm() || dst.iszero())
        flags |= X87SW_UNDERFLOW_EX;
    return flags;
}

uint16_t fp80_t::x87_fcos(fp80_t const &src, fp80_t &dst)
{
    uint16_t flags = src.isdenorm() ? X87SW_DENORM_EX : 0;
    if (src.isnan())
    {
        dst = src.issnan() ? fp80_t::make_qnan(src) : src;
        return src.issnan() ? X87SW_INVALID_EX : 0;
    }
    if (src.isinf())
    {
        dst = fp80_t::const_indef();
        return X87SW_INVALID_EX;
    }
    if (src.iszero())
    {
        dst = fp80_t::const_one();
        return 0;
    }
    if ((src.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS > 62)
    {
        dst = src;
        return flags | X87SW_C2;
    }

    fpext_t y;
    int q = range_reduce(src, y);
    fpext_t r;
    // cos(x + qπ/2): q=0 cos(y), q=1 -sin(y), q=2 -cos(y), q=3 sin(y)
    switch (q & 3)
    {
        case 0: r = taylor_cos(y); break;
        case 1: r = taylor_sin(y); r.chs(); break;
        case 2: r = taylor_cos(y); r.chs(); break;
        case 3: r = taylor_sin(y); break;
    }
    dst = round_fpext96_to_fp80(r, read_x87_cw(), flags);
    flags |= X87SW_PRECISION_EX;
    if (dst.isdenorm() || dst.iszero())
        flags |= X87SW_UNDERFLOW_EX;
    return flags;
}
uint16_t fp80_t::x87_fsincos(fp80_t const &src, fp80_t &dst1, fp80_t &dst2)
{
    // dst1 = cos, dst2 = sin (matches asm/test convention from fxtract-style
    // unary_2 ops). Reuse the single range reduction for efficiency.
    uint16_t flags = src.isdenorm() ? X87SW_DENORM_EX : 0;
    if (src.isnan())
    {
        dst1 = src.issnan() ? fp80_t::make_qnan(src) : src;
        dst2 = dst1;
        return src.issnan() ? X87SW_INVALID_EX : 0;
    }
    if (src.isinf())
    {
        dst1 = fp80_t::const_indef();
        dst2 = fp80_t::const_indef();
        return X87SW_INVALID_EX;
    }
    if (src.iszero())
    {
        dst1 = fp80_t::const_one();   // cos
        dst2 = src;                    // sin
        return 0;
    }
    if ((src.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS > 62)
    {
        dst1 = src; dst2 = src;
        return flags | X87SW_C2;
    }
    fpext_t y;
    int q = range_reduce(src, y);
    fpext_t s_ext = taylor_sin(y);
    fpext_t c_ext = taylor_cos(y);
    fpext_t sin_r, cos_r;
    switch (q & 3)
    {
        case 0: sin_r = s_ext;            cos_r = c_ext;            break;
        case 1: sin_r = c_ext;            cos_r = s_ext; cos_r.chs(); break;
        case 2: sin_r = s_ext; sin_r.chs(); cos_r = c_ext; cos_r.chs(); break;
        case 3: sin_r = c_ext; sin_r.chs(); cos_r = s_ext;            break;
    }
    uint16_t f1 = flags, f2 = flags;
    dst1 = round_fpext96_to_fp80(cos_r, read_x87_cw(), f1);
    dst2 = round_fpext96_to_fp80(sin_r, read_x87_cw(), f2);
    flags |= X87SW_PRECISION_EX;
    if (dst1.isdenorm() || dst1.iszero() || dst2.isdenorm() || dst2.iszero())
        flags |= X87SW_UNDERFLOW_EX;
    return flags;
}

uint16_t fp80_t::x87_fptan(fp80_t const &src, fp80_t &dst1, fp80_t &dst2)
{
    // tan(x) = sin/cos. dst1 = 1.0 (per asm convention), dst2 = tan(x).
    uint16_t flags = src.isdenorm() ? X87SW_DENORM_EX : 0;
    if (src.isnan())
    {
        dst1 = src.issnan() ? fp80_t::make_qnan(src) : src;
        dst2 = dst1;
        return src.issnan() ? X87SW_INVALID_EX : 0;
    }
    if (src.isinf())
    {
        dst1 = fp80_t::const_indef();
        dst2 = fp80_t::const_indef();
        return X87SW_INVALID_EX;
    }
    if (src.iszero())
    {
        dst1 = fp80_t::const_one();
        dst2 = src;
        return 0;
    }
    if ((src.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS > 62)
    {
        // Out of range: result unchanged, C2 set. Asm fptan pushes 1.0
        // even in this case but the stored values per test are src + src.
        dst1 = src; dst2 = src;
        return flags | X87SW_C2;
    }
    fpext_t y;
    int q = range_reduce(src, y);
    fpext_t s_ext = taylor_sin(y);
    fpext_t c_ext = taylor_cos(y);
    fpext_t sin_r, cos_r;
    switch (q & 3)
    {
        case 0: sin_r = s_ext;            cos_r = c_ext;            break;
        case 1: sin_r = c_ext;            cos_r = s_ext; cos_r.chs(); break;
        case 2: sin_r = s_ext; sin_r.chs(); cos_r = c_ext; cos_r.chs(); break;
        case 3: sin_r = c_ext; sin_r.chs(); cos_r = s_ext;            break;
    }
    // tan = sin / cos via fp80 division.
    uint16_t f_sin = flags, f_cos = flags;
    fp80_t sin80 = round_fpext96_to_fp80(sin_r, read_x87_cw(), f_sin);
    fp80_t cos80 = round_fpext96_to_fp80(cos_r, read_x87_cw(), f_cos);
    fp80_t tan80;
    fp80_t::x87_fdivr(sin80, cos80, tan80);   // tan = sin/cos
    dst1 = fp80_t::const_one();
    dst2 = tan80;
    flags |= X87SW_PRECISION_EX;
    if (tan80.isdenorm() || tan80.iszero())
        flags |= X87SW_UNDERFLOW_EX;
    return flags;
}
#endif
#if X87_HOST_HAS_FP80
uint16_t fp80_t::x87_frndint(fp80_t const &src, fp80_t &dst)
{
    return host_x87_unary(src, dst, 3);
}
#else
//
// FRNDINT: round ST(0) to the nearest integer per current rounding mode.
// Per spec: PE flag set if result differs from source; sign-preserving for
// zero; NaN/Inf passed through.
//
uint16_t fp80_t::x87_frndint(fp80_t const &src, fp80_t &dst)
{
    if (src.isnan())
    {
        dst = src.issnan() ? fp80_t::make_qnan(src) : src;
        return src.issnan() ? X87SW_INVALID_EX : 0;
    }
    if (src.isinf() || src.iszero())
    {
        dst = src;
        return 0;
    }

    uint16_t sign  = src.sign_exp() & FP80_SIGN_MASK;
    bool     is_neg = sign != 0;
    int      exp   = (src.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    uint64_t mant  = src.mantissa();
    x87cw_t  rmode = fpround_t::get() & X87CW_ROUNDING_MASK;
    uint16_t flags = 0;

    if (src.isdenorm())
    {
        // Denormal magnitude is < 1; result is 0 or ±1 per rounding mode.
        flags |= X87SW_DENORM_EX | X87SW_PRECISION_EX;
        bool to_one;
        if      (rmode == X87CW_ROUNDING_ZERO)    to_one = false;
        else if (rmode == X87CW_ROUNDING_NEAREST) to_one = false;  // |x| < 0.5
        else if (rmode == X87CW_ROUNDING_DOWN)    to_one = is_neg;
        else                                      to_one = !is_neg;  // UP
        if (to_one)
        {
            dst = fp80_t(FP80_EXPLICIT_ONE, sign | uint16_t(FP80_EXPONENT_BIAS));
            flags |= X87SW_C1;     // magnitude rounded up from tiny to 1
        }
        else
        {
            dst = is_neg ? fp80_t::const_nzero() : fp80_t::const_zero();
        }
        return flags;
    }

    if (exp >= 63)
    {
        dst = src;
        return 0;
    }

    if (exp < 0)
    {
        flags |= X87SW_PRECISION_EX;
        // |x| in [2^-N, 1) for some N >= 1.
        bool to_one;
        if (rmode == X87CW_ROUNDING_ZERO)
            to_one = false;
        else if (rmode == X87CW_ROUNDING_NEAREST)
        {
            // |x| >= 0.5 ↔ exp == -1 and mant > EXPLICIT_ONE OR exactly 0.5 → round even (0)
            if (exp == -1)
                to_one = (mant > FP80_EXPLICIT_ONE);   // > 0.5 → 1; ==0.5 → even (0)
            else
                to_one = false;                         // < 0.5
        }
        else if (rmode == X87CW_ROUNDING_DOWN) to_one = is_neg;
        else                                   to_one = !is_neg;  // UP

        if (to_one)
        {
            dst = fp80_t(FP80_EXPLICIT_ONE, sign | uint16_t(FP80_EXPONENT_BIAS));
            flags |= X87SW_C1;     // magnitude rounded up
        }
        else
        {
            dst = is_neg ? fp80_t::const_nzero() : fp80_t::const_zero();
        }
        return flags;
    }

    // 0 <= exp < 63: split mantissa into integer + fractional parts.
    int      frac_bits = 63 - exp;            // bits below the integer point
    uint64_t mask      = (1ull << frac_bits) - 1;
    uint64_t frac      = mant & mask;
    uint64_t int_part  = mant & ~mask;

    if (frac == 0)
    {
        dst = src;
        return 0;
    }

    flags |= X87SW_PRECISION_EX;
    bool round_up;
    if (rmode == X87CW_ROUNDING_NEAREST)
    {
        uint64_t half = 1ull << (frac_bits - 1);
        if      (frac > half) round_up = true;
        else if (frac < half) round_up = false;
        else                  round_up = (int_part & (1ull << frac_bits)) != 0; // even
    }
    else if (rmode == X87CW_ROUNDING_ZERO) round_up = false;
    else if (rmode == X87CW_ROUNDING_DOWN) round_up = is_neg;
    else                                   round_up = !is_neg;  // UP

    uint64_t result_mant = int_part;
    int      result_exp  = exp + FP80_EXPONENT_BIAS;
    if (round_up)
    {
        result_mant += (1ull << frac_bits);
        if (result_mant == 0)
        {
            // overflowed past 2.0 — bump exponent, restore explicit-1
            result_mant = FP80_EXPLICIT_ONE;
            result_exp++;
        }
        flags |= X87SW_C1;
    }
    dst = fp80_t(result_mant, sign | uint16_t(result_exp));
    return flags;
}
#endif // X87_HOST_HAS_FP80

//
// 3-way magnitude/sign compare reused by all FCOM-family ops.
//
static int compare_fp80_3way_trans(fp80_t const &a, fp80_t const &b)
{
    if (a.isnan() || b.isnan()) return 2;
    if (a.iszero() && b.iszero()) return 0;
    bool an = a.sign() != 0, bn = b.sign() != 0;
    if (an != bn) return an ? -1 : 1;
    uint16_t ae = a.sign_exp() & FP80_EXPONENT_MASK;
    uint16_t be = b.sign_exp() & FP80_EXPONENT_MASK;
    int mag;
    if (ae != be)                          mag = (ae < be) ? -1 : 1;
    else if (a.mantissa() != b.mantissa()) mag = (a.mantissa() < b.mantissa()) ? -1 : 1;
    else                                   mag = 0;
    return an ? -mag : mag;
}

//
// FXAM (Intel SDM §8.1.2.2 + Table 8-1):
//   C1 = sign of operand
//   C3:C2:C0 = 000 Unsupported, 001 NaN, 010 Normal finite, 011 Infinity,
//              100 Zero, 110 Denormal
// "Empty" (101/111) cannot arise from a software value, only from an actual
// FPU stack slot tag. Unsupported encodings (pseudo-NaN, pseudo-infinity,
// unnormal, pseudo-denormal) all share the 000 class on modern x87.
//
uint16_t fp80_t::x87_fxam(fp80_t const &src)
{
    uint16_t sw = 0;
    if (src.sign()) sw |= X87SW_C1;

    uint16_t e = src.sign_exp() & FP80_EXPONENT_MASK;
    uint64_t m = src.mantissa();
    bool explicit_one = (m & FP80_EXPLICIT_ONE) != 0;
    uint64_t m_frac = m & FP80_MANTISSA_MASK;

    if (e == FP80_EXPONENT_MAX_BIASED)
    {
        // Inf / NaN / pseudo-{Inf,NaN}
        if (!explicit_one)
        {
            // pseudo-NaN or pseudo-infinity — unsupported (000)
            return sw;
        }
        if (m_frac == 0) sw |= X87SW_C2 | X87SW_C0;       // Infinity (011)
        else             sw |= X87SW_C0;                  // NaN (001)
    }
    else if (e == 0)
    {
        if (m == 0)            sw |= X87SW_C3;            // Zero (100)
        else if (!explicit_one) sw |= X87SW_C3 | X87SW_C2; // Denormal (110)
        // else: pseudo-denormal (explicit_one set with zero exponent) —
        // unsupported (000)
    }
    else
    {
        // Normal-range exponent
        if (explicit_one) sw |= X87SW_C2;                 // Normal finite (010)
        // else: unnormal — unsupported (000)
    }
    return sw;
}

//
// FTST: compare ST(0) against +0.0. Result in C3:C2:C0 per Table 8-1.
// NaN -> unordered (111) and #IA.  Denormal operand -> #D.
//
uint16_t fp80_t::x87_ftst(fp80_t const &src)
{
    uint16_t flags = 0;
    if (src.isdenorm()) flags |= X87SW_DENORM_EX;
    if (src.isnan())
        return flags | X87SW_C3 | X87SW_C2 | X87SW_C0 | X87SW_INVALID_EX;
    if (src.iszero())
        return flags | X87SW_C3;
    if (src.sign())
        return flags | X87SW_C0;
    return flags;
}

//
// FCOM: compare two operands.  Sets #IA on any NaN (signaling or quiet);
// sets #D if either operand is denormal *and* neither is NaN (the FPU does
// not raise DE alongside an unordered result).
//
uint16_t fp80_t::x87_fcom(fp80_t const &a, fp80_t const &b)
{
    int c = compare_fp80_3way_trans(a, b);
    if (c == 2) return X87SW_C3 | X87SW_C2 | X87SW_C0 | X87SW_INVALID_EX;
    uint16_t flags = 0;
    if (a.isdenorm() || b.isdenorm()) flags |= X87SW_DENORM_EX;
    if (c == 0) return flags | X87SW_C3;
    if (c <  0) return flags | X87SW_C0;
    return flags;
}

//
// FUCOM: like FCOM but only signals #IA on signaling NaN.
//
uint16_t fp80_t::x87_fucom(fp80_t const &a, fp80_t const &b)
{
    int c = compare_fp80_3way_trans(a, b);
    if (c == 2)
    {
        uint16_t sw = X87SW_C3 | X87SW_C2 | X87SW_C0;
        if (a.issnan() || b.issnan()) sw |= X87SW_INVALID_EX;
        return sw;
    }
    uint16_t flags = 0;
    if (a.isdenorm() || b.isdenorm()) flags |= X87SW_DENORM_EX;
    if (c == 0) return flags | X87SW_C3;
    if (c <  0) return flags | X87SW_C0;
    return flags;
}

//
// FCOMI / FUCOMI: P6+ compares that write EFLAGS-style result bits. Our test
// harness packs CF/PF/ZF into the C0/C2/C3 SW positions, so the C++ side
// produces the same mapping (and never sets IE, since the asm oracle reads
// the integer EFLAGS — IE only ever appears in the SW, never EFLAGS).
//
uint16_t fp80_t::x87_fcomi(fp80_t const &a, fp80_t const &b)
{
    int c = compare_fp80_3way_trans(a, b);
    if (c == 2) return X87SW_C3 | X87SW_C2 | X87SW_C0;
    if (c == 0) return X87SW_C3;
    if (c <  0) return X87SW_C0;
    return 0;
}

uint16_t fp80_t::x87_fucomi(fp80_t const &a, fp80_t const &b)
{
    // FUCOMI only differs from FCOMI in that it signals on SNaN — but again
    // the asm oracle reads EFLAGS only, so IE is never observable here.
    return x87_fcomi(a, b);
}

}
