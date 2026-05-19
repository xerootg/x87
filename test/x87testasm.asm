; Assemble (Windows): nasm -f win64 x87testasm.asm -o x87testasm.obj
; Assemble (Linux):   nasm -f elf64 x87testasm.asm -o x87testasm.o
;
; ABI argument register macros: Win64 uses rcx/rdx/r8 for the first 3
; integer args; SysV (Linux/macOS) uses rdi/rsi/rdx.
%ifidn __OUTPUT_FORMAT__,win64
    %define ARG1 rcx
    %define ARG2 rdx
    %define ARG3 r8
%else
    %define ARG1 rdi
    %define ARG2 rsi
    %define ARG3 rdx
%endif

    section .data

saved_cw dw 0

    ; Library-defined bit patterns for SNaN / QNaN. No x87 instruction
    ; directly generates these specific payloads; we load them via FLD TWORD
    ; (which is a verbatim 10-byte copy for the extended format) so the
    ; constants test round-trips them through the FPU.
    ;
    ; layout: 8 bytes mantissa (LE), then 2 bytes sign_exp (LE)
const_snan80: db 0x01, 0, 0, 0, 0, 0, 0, 0x80, 0xff, 0x7f  ; 0x8000_0000_0000_0001 / 0x7fff
const_qnan80: db 0x01, 0, 0, 0, 0, 0, 0, 0xc0, 0xff, 0x7f  ; 0xc000_0000_0000_0001 / 0x7fff

    section .text

    bits 64

    ; The destination buffer is an array of packed 10-byte fp80_t values
    ; (matches sizeof(fp80_t) under #pragma pack(2)), so stride is 10.
    ; Layout (13 slots):
    ;    0  one    (fld1)
    ;    1  l2t    (fldl2t)
    ;    2  l2e    (fldl2e)
    ;    3  pi     (fldpi)
    ;    4  lg2    (fldlg2)
    ;    5  ln2    (fldln2)
    ;    6  zero   (fldz)
    ;    7  nzero  (fldz; fchs)
    ;    8  pinf   (1.0 / 0.0)
    ;    9  ninf   (-1.0 / 0.0)
    ;   10  snan   (round-tripped from .data)
    ;   11  qnan   (round-tripped from .data)
    ;   12  indef  (0.0 / 0.0)
    global x87consts80
x87consts80:
    finit
    fldcw   [rel saved_cw]
    fld1
    fstp    tword [ARG1 + 0*10]
    fldl2t
    fstp    tword [ARG1 + 1*10]
    fldl2e
    fstp    tword [ARG1 + 2*10]
    fldpi
    fstp    tword [ARG1 + 3*10]
    fldlg2
    fstp    tword [ARG1 + 4*10]
    fldln2
    fstp    tword [ARG1 + 5*10]
    fldz
    fstp    tword [ARG1 + 6*10]

    fldz
    fchs
    fstp    tword [ARG1 + 7*10]

    fld1
    fldz
    fdivp                                ; +inf
    fstp    tword [ARG1 + 8*10]

    fld1
    fchs
    fldz
    fdivp                                ; -inf
    fstp    tword [ARG1 + 9*10]

    fld     tword [rel const_snan80]
    fstp    tword [ARG1 + 10*10]

    fld     tword [rel const_qnan80]
    fstp    tword [ARG1 + 11*10]

    fldz
    fldz
    fdivp                                ; indefinite
    fstp    tword [ARG1 + 12*10]
    ret

    global x87consts64
x87consts64:
    fld1
    fstp    qword [ARG1 + 0*8]
    fldl2t
    fstp    qword [ARG1 + 1*8]
    fldl2e
    fstp    qword [ARG1 + 2*8]
    fldpi
    fstp    qword [ARG1 + 3*8]
    fldlg2
    fstp    qword [ARG1 + 4*8]
    fldln2
    fstp    qword [ARG1 + 5*8]
    fldz
    fstp    qword [ARG1 + 6*8]
    ret

    global x87setcw
x87setcw:
    fldcw   [ARG1]
    fstcw   [rel saved_cw]
    ret

    global x87getsw
x87getsw:
    fstsw   ax
    movzx   eax, ax
    ret

    global x87test1
x87test1:
    fld     qword [ARG1]
    fstp    tword [ARG2]
    ret

    global x87test2
x87test2:
    fld     tword [ARG1]
    fstp    qword [ARG2]
    ret

    global fld8080
fld8080:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fstsw   ax
    fstp    tword [ARG2]
    ret

    global fld6480
fld6480:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG1]
    fstsw   ax
    fstp    tword [ARG2]
    ret

    global fld3280
fld3280:
    finit
    fldcw   [rel saved_cw]
    fld     dword [ARG1]
    fstsw   ax
    fstp    tword [ARG2]
    ret

    global fild6480
fild6480:
    finit
    fldcw   [rel saved_cw]
    fild    qword [ARG1]
    fstsw   ax
    fstp    tword [ARG2]
    ret

    global fild3280
fild3280:
    finit
    fldcw   [rel saved_cw]
    fild    dword [ARG1]
    fstsw   ax
    fstp    tword [ARG2]
    ret

    global fild1680
fild1680:
    finit
    fldcw   [rel saved_cw]
    fild    word [ARG1]
    fstsw   ax
    fstp    tword [ARG2]
    ret

    global fst8080
fst8080:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fstp    tword [ARG2]
    fstsw   ax
    ret

    global fst8064
fst8064:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fstp    qword [ARG2]
    fstsw   ax
    ret

    global fst8032
fst8032:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fstp    dword [ARG2]
    fstsw   ax
    ret

    global fist8064
fist8064:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fistp   qword [ARG2]
    fstsw   ax
    ret

    global fist8032
fist8032:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fistp   dword [ARG2]
    fstsw   ax
    ret

    global fist8016
fist8016:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fistp   word [ARG2]
    fstsw   ax
    ret

    global fadd80
fadd80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    faddp
    fstsw   ax
    fstp    tword [ARG3]
    ret

    global fsub80
fsub80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fsubp
    fstsw   ax
    fstp    tword [ARG3]
    ret

    global fmul80
fmul80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fmulp
    fstsw   ax
    fstp    tword [ARG3]
    ret

    global fdiv80
fdiv80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fdivp
    fstsw   ax
    fstp    tword [ARG3]
    ret

    global fsqrt80
fsqrt80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fsqrt
    fstsw   ax
    fstp    tword [ARG2]
    ret

    global f2xm180
f2xm180:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    f2xm1
    fstsw   ax
    fstp    tword [ARG2]
    ret

    global fyl2x80
fyl2x80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fyl2x
    fstsw   ax
    fstp    tword [ARG3]
    ret

    global fptan80
fptan80:
    finit
    fldcw   [rel saved_cw]
    fldz
    fld     tword [ARG1]
    fptan
    fstsw   ax
    fstp    tword [ARG2]
    fstp    tword [ARG3]
    ret

    global fsincos80
fsincos80:
    finit
    fldcw   [rel saved_cw]
    fldz
    fld     tword [ARG1]
    fsincos
    fstsw   ax
    fstp    tword [ARG2]
    fstp    tword [ARG3]
    ret

    global fpatan80
fpatan80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fpatan
    fstsw   ax
    fstp    tword [ARG3]
    ret

    global fxtract80
fxtract80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fxtract
    fstsw   ax
    fstp    tword [ARG2]
    fstp    tword [ARG3]
    ret

    global fprem180
fprem180:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fprem1
    fstsw   ax
    fstp    tword [ARG3]
    ret

    global fprem80
fprem80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fprem
    fstsw   ax
    fstp    tword [ARG3]
    ret

    global fyl2xp180
fyl2xp180:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fyl2xp1
    fstsw   ax
    fstp    tword [ARG3]
    ret

    global frndint80
frndint80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    frndint
    fstsw   ax
    fstp    tword [ARG2]
    ret

    global fscale80
fscale80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fscale
    fstsw   ax
    fstp    tword [ARG3]
    ret

    global fsin80
fsin80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fsin
    fstsw   ax
    fstp    tword [ARG2]
    ret

    global fcos80
fcos80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fcos
    fstsw   ax
    fstp    tword [ARG2]
    ret

    ; Unordered compare of two 80-bit values. Loads src2 then src1, executes
    ; fucompp (compare ST(0) to ST(1), pop both), returns the resulting status
    ; word in ax. Caller inspects C0/C2/C3 to decode the relation.
    global fcmp80
fcmp80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fucompp
    fstsw   ax
    ret

    ; Classify an 80-bit value via FXAM. Returns the resulting status word in
    ; ax; C3:C2:C0 encode the class and C1 is the sign bit.
    global fxam80
fxam80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fxam
    fstsw   ax
    fstp    st0
    ret

    global fabs80
fabs80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fabs
    fstsw   ax
    fstp    tword [ARG2]
    ret

    global fchs80
fchs80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    fchs
    fstsw   ax
    fstp    tword [ARG2]
    ret

    ; FSUBR / FDIVR: reverse-operand variants. Matching the convention of
    ; fsub80/fdiv80 above so test_binary80 passes them src2 then src1.
    global fsubr80
fsubr80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fsubrp
    fstsw   ax
    fstp    tword [ARG3]
    ret

    global fdivr80
fdivr80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fdivrp
    fstsw   ax
    fstp    tword [ARG3]
    ret

    ; FTST compares ST(0) against +0.0 and sets C0/C2/C3 per Table 8-1.
    global ftst80
ftst80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG1]
    ftst
    fstsw   ax
    fstp    st0
    ret

    ; FCOM (signals on QNaN as well as SNaN). Same pop pattern as fucompp.
    global fcom80
fcom80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fcompp
    fstsw   ax
    ret

    ; FCOMI / FUCOMI: P6+ compares that write CF/PF/ZF in EFLAGS rather than
    ; C0/C2/C3 in SW. We pack the EFLAGS result into a uint16_t at the same
    ; bit positions as the SW would use so the same test machinery applies:
    ;   CF -> C0, PF -> C2, ZF -> C3
    global fcomi80
fcomi80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fcomi   st1
    pushfq
    pop     rax                          ; full RFLAGS in rax
    mov     rcx, rax                     ; copy
    ; reassemble: bit 0 (CF) -> C0 (bit 8), bit 2 (PF) -> C2 (bit 10), bit 6 (ZF) -> C3 (bit 14)
    and     rax, 1                       ; CF
    shl     rax, 8                       ; -> C0
    mov     rdx, rcx
    and     rdx, 4                       ; PF
    shl     rdx, 8                       ; -> C2 (4 << 8 = 0x400 = bit 10)
    or      rax, rdx
    mov     rdx, rcx
    and     rdx, 0x40                    ; ZF
    shl     rdx, 8                       ; -> C3 (0x40 << 8 = 0x4000 = bit 14)
    or      rax, rdx
    push    rax                          ; stash result
    fstp    st0
    fstp    st0                          ; clear both operands from x87 stack
    pop     rax
    ret

    global fucomi80
fucomi80:
    finit
    fldcw   [rel saved_cw]
    fld     tword [ARG2]
    fld     tword [ARG1]
    fucomi  st1
    pushfq
    pop     rax
    mov     rcx, rax
    and     rax, 1
    shl     rax, 8
    mov     rdx, rcx
    and     rdx, 4
    shl     rdx, 8
    or      rax, rdx
    mov     rdx, rcx
    and     rdx, 0x40
    shl     rdx, 8
    or      rax, rdx
    push    rax
    fstp    st0
    fstp    st0
    pop     rax
    ret

    global fadd64
fadd64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG2]
    fld     qword [ARG1]
    faddp
    fstsw   ax
    fstp    qword [ARG3]
    ret

    global fsub64
fsub64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG2]
    fld     qword [ARG1]
    fsubp
    fstsw   ax
    fstp    qword [ARG3]
    ret

    global fmul64
fmul64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG2]
    fld     qword [ARG1]
    fmulp
    fstsw   ax
    fstp    qword [ARG3]
    ret

    global fdiv64
fdiv64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG2]
    fld     qword [ARG1]
    fdivp
    fstsw   ax
    fstp    qword [ARG3]
    ret

    global fsqrt64
fsqrt64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG1]
    fsqrt
    fstsw   ax
    fstp    qword [ARG2]
    ret

    global f2xm164
f2xm164:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG1]
    f2xm1
    fstsw   ax
    fstp    qword [ARG2]
    ret

    global fyl2x64
fyl2x64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG2]
    fld     qword [ARG1]
    fyl2x
    fstsw   ax
    fstp    qword [ARG3]
    ret

    global fptan64
fptan64:
    finit
    fldcw   [rel saved_cw]
    fldz
    fld     qword [ARG1]
    fptan
    fstsw   ax
    fstp    qword [ARG2]
    fstp    qword [ARG3]
    ret

    global fsincos64
fsincos64:
    finit
    fldcw   [rel saved_cw]
    fldz
    fld     qword [ARG1]
    fsincos
    fstsw   ax
    fstp    qword [ARG2]
    fstp    qword [ARG3]
    ret

    global fpatan64
fpatan64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG2]
    fld     qword [ARG1]
    fpatan
    fstsw   ax
    fstp    qword [ARG3]
    ret

    global fxtract64
fxtract64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG1]
    fxtract
    fstsw   ax
    fstp    qword [ARG2]
    fstp    qword [ARG3]
    ret

    global fprem164
fprem164:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG2]
    fld     qword [ARG1]
    fprem1
    fstsw   ax
    fstp    qword [ARG3]
    ret

    global fprem64
fprem64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG2]
    fld     qword [ARG1]
    fprem
    fstsw   ax
    fstp    qword [ARG3]
    ret

    global fyl2xp164
fyl2xp164:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG2]
    fld     qword [ARG1]
    fyl2xp1
    fstsw   ax
    fstp    qword [ARG3]
    ret

    global frndint64
frndint64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG1]
    frndint
    fstsw   ax
    fstp    qword [ARG2]
    ret

    global fscale64
fscale64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG2]
    fld     qword [ARG1]
    fscale
    fstsw   ax
    fstp    qword [ARG3]
    ret

    global fsin64
fsin64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG1]
    fsin
    fstsw   ax
    fstp    qword [ARG2]
    ret

    global fcos64
fcos64:
    finit
    fldcw   [rel saved_cw]
    fld     qword [ARG1]
    fcos
    fstsw   ax
    fstp    qword [ARG2]
    ret

%ifidn __OUTPUT_FORMAT__,elf64
    ; Mark stack as non-executable (ELF only).
    section .note.GNU-stack noalloc noexec nowrite progbits
%endif
