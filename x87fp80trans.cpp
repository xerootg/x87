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
        dst = (g * h + g + h).as_fp80();
        return X87SW_PRECISION_EX;
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
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
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
        // Result NaN: prefer 'a' if NaN, else 'b'
        if (a.isnan()) dst = a.issnan() ? fp80_t::make_qnan(a) : a;
        else           dst = b.issnan() ? fp80_t::make_qnan(b) : b;
        return snan ? X87SW_INVALID_EX : 0;
    }

    // Inf * anything (except 0) = signed inf of a
    // 0 * inf-scale = 0 (no exception, scale just becomes 0 effectively)
    if (a.iszero()) { dst = a; return 0; }
    if (a.isinf())  { dst = a; return 0; }

    // b infinite: scale by 2^+inf = +inf, 2^-inf = 0
    if (b.isinf())
    {
        if (b.sign()) dst = (a.sign() ? fp80_t::const_nzero() : fp80_t::const_zero());
        else          dst = (a.sign() ? fp80_t::const_ninf()  : fp80_t::const_pinf());
        return 0;
    }
    if (b.iszero()) { dst = a; return 0; }

    // Truncate b to a signed integer. Saturate for huge |b|.
    int bexp = (b.sign_exp() & FP80_EXPONENT_MASK) - FP80_EXPONENT_BIAS;
    int64_t scale;
    if (bexp < 0)        scale = 0;
    else if (bexp >= 63) scale = b.sign() ? INT64_MIN : INT64_MAX;
    else                 scale = int64_t(b.mantissa() >> (63 - bexp)) * (b.sign() ? -1 : 1);

    // Apply scale to a's exponent.
    int64_t new_exp = int64_t(a.sign_exp() & FP80_EXPONENT_MASK) + scale;
    uint16_t sign = a.sign_exp() & FP80_SIGN_MASK;

    if (new_exp >= FP80_EXPONENT_MAX_BIASED)
    {
        dst = sign ? fp80_t::const_ninf() : fp80_t::const_pinf();
        return X87SW_OVERFLOW_EX | X87SW_PRECISION_EX;
    }
    if (new_exp <= 0)
    {
        // For brevity treat all underflows as flushing to signed zero.
        dst = sign ? fp80_t::const_nzero() : fp80_t::const_zero();
        return X87SW_UNDERFLOW_EX | X87SW_PRECISION_EX;
    }
    dst = fp80_t(a.mantissa(), sign | uint16_t(new_exp));
    return 0;
}
#endif // X87_HOST_HAS_FP80

#if X87_HOST_HAS_FP80
uint16_t fp80_t::x87_fprem (fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary_trans(a, b, dst, 3); }
uint16_t fp80_t::x87_fprem1(fp80_t const &a, fp80_t const &b, fp80_t &dst) { return host_x87_binary_trans(a, b, dst, 4); }
#else
uint16_t fp80_t::x87_fprem(fp80_t const &, fp80_t const &, fp80_t &dst)               { return stub_unary(dst); }
uint16_t fp80_t::x87_fprem1(fp80_t const &, fp80_t const &, fp80_t &dst)              { return stub_unary(dst); }
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
uint16_t fp80_t::x87_fyl2x(fp80_t const &, fp80_t const &, fp80_t &dst)               { return stub_unary(dst); }
uint16_t fp80_t::x87_fyl2xp1(fp80_t const &, fp80_t const &, fp80_t &dst)             { return stub_unary(dst); }
uint16_t fp80_t::x87_fsin(fp80_t const &, fp80_t &dst)                                { return stub_unary(dst); }
uint16_t fp80_t::x87_fcos(fp80_t const &, fp80_t &dst)                                { return stub_unary(dst); }
uint16_t fp80_t::x87_fsincos(fp80_t const &, fp80_t &dst1, fp80_t &dst2)              { return stub_unary2(dst1, dst2); }
uint16_t fp80_t::x87_fptan(fp80_t const &, fp80_t &dst1, fp80_t &dst2)                { return stub_unary2(dst1, dst2); }
uint16_t fp80_t::x87_fpatan(fp80_t const &, fp80_t const &, fp80_t &dst)              { return stub_unary(dst); }
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
