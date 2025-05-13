.text
.balign 4
.global ntt_inner
# fn ntt_inner(element: *mut u64, omega: *const u64, modulus: u64, precon_omega: *const u64, j1: u64, t: u64);
ntt_inner:
# here starts mod_mul_fast_const, vector extension version
# based on optimized version of mod_mul_fast_const
    li t0, 1            # parallel instructions = 1
    vsetvli t0, t0, e64, m8, tu, ma  # arg0 = items to be processed, arg1 = remaining items
    add t0, a4, a5
    slli t0, t0, 3
    add t0, t0, a0
    vle64.v v0, (t0)    # v0 = elements[0:t0]
    vle64.v v8, (a3)   # v8 = precon_omega[0:t0]
    vle64.v v16, (a1)   # v16 = omega[0:t0]


# Vectorial version of
# fn mod_mul_fast_const(a: u64, b: u64, modulus: u64, b_inv: u64) -> u64
    vmulhu.vv   v8, v8, v0
    vnot.v      v8, v8
    vmul.vv     v0, v16, v0
    vmul.vx     v16, v8, a2
    vadd.vv     v0, v0, v16
    li t1, 0x3f          # Load immediate 63 into t1
    vsra.vx     v16, v0, t1
    vand.vx     v16, v16, a2
    vadd.vv     v0, v0, v16
# mod_mul_fast_const finished


    vmv.x.s t1, v0         # v0 has new omega factor
#                  let mut lo_val = element[j1 as usize];
#                  let mut hi_val = lo_val + omega_factor;
#                  if hi_val >= modulus {
#                      hi_val -= modulus;
#                  }
#                  if lo_val < omega_factor {
#                      lo_val += modulus;
#                  }
#                  lo_val -= omega_factor;
#                  element[(j1 + 0) as usize] = hi_val;
#                  element[(j1 + t) as usize] = lo_val;
    slli    t2, a4, 3    # lo_val = element[j1];
    add     t2, t2, a0   # lo_val = element[j1];
    ld      t3, (t2)     # lo_val = element[j1];
    vle64.v v8, (t2)     # v8 lo_val[] = element[j1..];
    vadd.vv v0, v0, v8   # v0 hi_val = lo_val + omega
    vmv.x.s t3, v0
    sltu    t4, t3, a2   # t4 =  hi_val < modulus
    addi    t4, t4, -1   # t4 = 0 if hi_val < modulus, 0xffffffff if hi_val >= modulus
    and     t4, a2, t4   # t5 = modulus & t4
    sub     t4, t3, t4   # t5 = hi_val - modulus

    vmsltu.vx  v16, v0, a2   # t4 =  hi_val < modulus
    vmv.x.s t5, v16
    vadd.vi  v16, v16, -1   # t4 = 0 if hi_val < modulus, 0xffffffff if hi_val >= modulus
    vand.vx  v16, v16, a2  # v16: masked modulus
    vsub.vv  v16, v0, v16   # t5 = hi_val - modulus
    vmv.x.s t5, v16

    ld      t3, (t2)     # lo_val = element[j1];
    sd      t4, (t2)     # element[j1] = hi_val
    sltu    t4, t3, t1   # t4 =  lo_val < omega_factor
    sub     t4, zero, t4 # t4 = 0xffffffff if lo_val < omega_factor
    and     t4, a2, t4   # t4 = modulus & t4
    add     t4, t3, t4   # t4 = lo_val + masked modulus
    sub     t4, t4, t1   # t4 = lo_val - omega factor

    add t0, a4, a5
    slli t0, t0, 3
    add t0, t0, a0
    sd      t4, (t0)
    ret
