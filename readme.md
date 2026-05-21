# x87 library

This library is designed to give authors of x86 emulators a resource that can be used for implementing x87 operations.
It is currently used for x87 floating-point support in my DREAMM emulator.

The motivation behind this is that existing libraries providing software floating point operations tend to be general-purpose and require a lot of tweaking/massaging in order to make them viable to use for x87 emulation.
The goal of this library is provide a set of tested implementations that not only produce accurate results but also generate status flags and handle edge cases in ways that match x87 implementations.

The core headers are `x87fp64.h` and `x87fp80.h`, which provide two types `x87::fp64_t` and `x87::fp80_t` respectively.
These two types have identical interfaces and are intended to be easily swappable depending on your needs.
`x87::fp64_t` performs all math using native 64-bit double support on the current processor (assumes either x64 or ARM64).
`x87::fp80_t` by contrast performs all math operations by hand to full 80-bit precision.

## Implementation status

`x87::fp80_t` is now substantially complete. On x86_64 it delegates to the host x87 via inline assembly. On aarch64/arm64 (and when forced via `-DX87_FORCE_HAND_ROLLED`), all operations run through hand-rolled 80-bit code paths.

Matching Intel x87 hardware bit-exactly (mantissa + full status word) across the test suite at `test/x87test.cpp`:

| Operation set | Match rate |
|---|---|
| Constants, comparisons, FXAM, FTST, FXTRACT, FRNDINT, FLOOR/CEIL, ABS/CHS, copysign/samesign, NaN helpers | 100% |
| FADD / FSUB / FSUBR / FMUL / FDIV / FDIVR | 100% |
| FSCALE / FPREM / FPREM1 | 100% |
| FYL2X | 93% |
| FSQRT | 82% |
| FSIN / FCOS | ~80% |
| F2XM1 | 79% |
| FPATAN | 77% |
| FPTAN / FSINCOS | ~62% per output |
| FYL2XP1 | 55% |

The transcendentals' remaining gaps are almost entirely in the C1 (rounding-direction) status word bit — the mantissa values match Intel hardware in nearly every test case. Empirical investigation showed that higher internal precision (the `fpext128_t` template, with 64+64 = 128 mantissa+ext bits) actually *regresses* the C1 match because Intel's microcode uses a specific internal precision and polynomial evaluation chain that we can match approximately but not byte-exactly through software. The split-form polynomial sin(y) = y + y·(z·P(z)) used in `taylor_sin/taylor_cos` was derived empirically by probing Intel x87 — it matches native fp80 hardware 96% bit-exact, though our hand-rolled fp80 multiply chain diverges in the LSBs.

## Building and testing

The test harness in `test/` compares each implementation against a live x86 oracle (via NASM-generated inline-asm dispatchers). Build with `make` from `test/`. Override `CXXFLAGS` to add `-DX87_FORCE_HAND_ROLLED` to test the aarch64 code path on an x86 host. Cross-compile + qemu also works (see `test/x87aarch64_smoke.cpp`).

Feel free to use this code in your projects if it is useful.
And if you find any bugs or the motivation to enhance/improve it in any way, I am definitely open to any improvements!
