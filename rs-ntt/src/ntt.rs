use std::ffi::c_ulong;
use std::ops::{Index, IndexMut};
use std::ptr::{slice_from_raw_parts_mut};

#[repr(C)]
pub struct FfiVec {
    pub m_modulus: c_ulong,
    pub m_data: *mut c_ulong,
    pub m_size: c_ulong,
}

#[derive(Debug)]
pub struct FheVec<'t> {
    pub modulus: u64,
    pub data: &'t mut [u64]
}

impl From<&FfiVec> for FheVec<'_> {
    fn from(value: &FfiVec) -> Self {
        let modulus = value.m_modulus;
        let data = unsafe {&mut *slice_from_raw_parts_mut(value.m_data as *mut u64, value.m_size as usize)};
        Self {modulus, data}
    }
}

impl Index<usize> for FheVec<'_> {
    type Output = u64;
    fn index(&self, i: usize) -> & u64 {
        &self.data[i]
    }
}

impl IndexMut<usize> for FheVec<'_> {
    fn index_mut(&mut self, i: usize) -> &mut u64 {
        &mut self.data[i]
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn rust_ntt_transform(
    root_of_unity: *const FfiVec,
    precon_root_of_unity: *const FfiVec,
    element: *mut FfiVec,
) {
    let root_of_unity = FheVec::from(&*root_of_unity);
    let precon_root_of_unity = FheVec::from(&*precon_root_of_unity);
    let mut element = FheVec::from(&*element);
    ntt_transform_hybrid(&root_of_unity, &precon_root_of_unity, &mut element);
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn mod_mul_fast_wrap(a: u64, b: u64, modulus: u64, b_inv: u64) -> u64 {
    mod_mul_fast_const(a, b, modulus, b_inv)
}

#[cfg(target_arch = "riscv64")]
#[link(name="nttvec")]
unsafe extern "C" {
    fn ntt_inner(element: *mut u64, omega: *const u64, modulus: u64, precon_omega: *const u64, j1: u64, t: u64);
}

pub fn ntt_transform(root_of_unity: &FheVec, precon_root_of_unity: &FheVec, element: &mut FheVec) {
    let modulus = element.modulus;

    let n = element.data.len() as u64 >> 1;
    let mut m = 1;
    let mut t = n;
    let mut logt = get_msb(t);
    println!("Root of unity: {:?}", root_of_unity);
    println!("T: {:?}, log_t: {:?}", t, logt);
    while m < n {
        println!("Using m: {m}");
        println!("------------");
        for i in 0 .. m {
            let omega_idx = (i + m) as usize;
            let omega = root_of_unity[omega_idx];
            let precon_omega = precon_root_of_unity[omega_idx];
            let j1base = i << logt;
            let j2 = j1base + t;
            println!("i: {i}");
            println!("omega_idx: {omega_idx}, omega: {omega}, preccon_omega: {precon_omega}");
            println!("j1: {j1base}, j2: {j2}");
            for  j1 in j1base .. j2 {
                let omega_factor_old = element[(j1 + t) as usize];
                let omega_factor = mod_mul_fast_const(omega_factor_old, omega, modulus, precon_omega);
                let mut lo_val = element[j1 as usize];
                println!("Elements: {:?}", element.data);
                println!("Omega factor {omega_factor} - Lo val: {lo_val}");
                let mut hi_val = lo_val + omega_factor;
                if hi_val >= modulus {
                    hi_val -= modulus;
                }
                if lo_val < omega_factor {
                    lo_val += modulus;
                }
                lo_val -= omega_factor;
                element[(j1 + 0) as usize] = hi_val;
                element[(j1 + t) as usize] = lo_val;
                println!("element[{j1}] -> {hi_val}");
                println!("element[{}] -> {lo_val}", j1 + t);
            }
            println!("Element vector: {:?}", element);
            println!("-------------");
        }
        m = m << 1;
        t = t >> 1;
        logt -= 1;
    }
    // peeled off last ntt stage for performance
    let mut i = 0;
    while i < (n << 1) {
        let omega_factor = element[(i + t) as usize];
        let omega_idx = (i >> 1) + n;
        let omega = root_of_unity[omega_idx as usize];
        let precon_omega = precon_root_of_unity[omega_idx as usize];
        let omega_factor = mod_mul_fast_const(omega_factor, omega, modulus, precon_omega);
        let mut lo_val = (*element)[(i) as usize];
        let mut hi_val = lo_val + omega_factor;
        if hi_val >= modulus {
            hi_val -= modulus;
        }
        if lo_val < omega_factor {
            lo_val += modulus;
        }
        lo_val -= omega_factor;
        (*element)[(i + 0) as usize] = hi_val;
        (*element)[(i + 1) as usize] = lo_val;
        i += 2;
    }
    println!("Element vector: {:?}", element);
    println!("-------------");
}

const fn compute_mu(value: u64) -> u64 {
    if value == 0 {
        panic!("NativeIntegerT ComputeMu: Divide by zero");
    }
    let tmp = {1 << (2 * get_msb(value) + 3)};
    tmp / value
}

const fn get_msb(val: u64) -> u64 {
    64u64 - u64::leading_zeros(val) as u64
}

const fn mod_mul_fast_eq(a: u64, b: u64, modulus: u64, mu: u64) -> u64 {
    let n =  get_msb(modulus) - 2;
    let mut prod = a as u128 * b as u128;
    let orig = prod;
    let mv = modulus;
    prod = (prod >> n) * mu as u128;
    prod = (prod >> n + 7) * mv as u128;
    let rem = orig - prod;
    if rem as u64 >= mv {
        rem as u64 - mv
    } else {
        rem as u64
    }
}

const fn mod_mul_fast_const(a: u64, b: u64, modulus: u64, b_inv: u64) -> u64{
    let q = mul_hi(a, b_inv) + 1;
    let yprime = (a.overflowing_mul(b).0 as i64).overflowing_sub(q.overflowing_mul(modulus).0 as i64).0;
    if yprime >= 0 {
        yprime as u64
    } else {
        (yprime + modulus as i64) as u64
    }
}

const fn mul_hi(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) >> 64) as u64
}

pub fn ntt_transform_hybrid(root_of_unity: &FheVec, precon_root_of_unity: &FheVec, element: &mut FheVec) {
    let modulus = element.modulus;

    let n = element.data.len() as u64 >> 1;
    let mut m = 1;
    let mut t = n;
    let mut logt = get_msb(t);
    // println!("Input element: {:?}", element);
    // println!("Root of unity: {:?}", root_of_unity);
    // println!("Precon Root of unity: {:?}", precon_root_of_unity);
    // println!("T: {:?}, log_t: {:?}", t, logt);
    while m < n {
        //println!("Using m: {m}");
        //println!("------------");
        for i in 0 .. m {
            let omega_idx = (i + m) as usize;
            let omega = root_of_unity[omega_idx];
            let precon_omega = precon_root_of_unity[omega_idx];
            let j1base = i << logt;
            let j2 = j1base + t;
            // println!("i: {i}");
            // println!("omega_idx: {omega_idx}, omega: {omega}, precon_omega: {precon_omega}");
            // println!("j1: {j1base}, j2: {j2}");
            for  j1 in j1base .. j2 {
                #[cfg(target_arch = "riscv64")]
                {
                    let element_unsafe = element.data.as_mut_ptr();
                    let omega_unsafe = (&omega) as *const _;
                    let precon_omega_unsafe = (&precon_omega) as *const _;
                    unsafe { ntt_inner(element_unsafe, omega_unsafe, modulus, precon_omega_unsafe, j1, t) };
                }
                #[cfg(target_arch = "x86_64")]
                {
                    let omega_factor_old = element[(j1 + t) as usize];
                    let omega_factor = mod_mul_fast_const(omega_factor_old, omega, modulus, precon_omega);
                    let mut lo_val = element[j1 as usize];
                    let mut hi_val = lo_val + omega_factor;
                    if hi_val >= modulus {
                        hi_val -= modulus;
                    }
                    if lo_val < omega_factor {
                        lo_val += modulus;
                    }
                    lo_val -= omega_factor;
                    element[j1 as usize] = hi_val;
                    element[(j1 + t) as usize] = lo_val;
                }
            }
            // println!("Element vector: {:?}", element);
            // println!("-------------");
        }
        m = m << 1;
        t = t >> 1;
        logt -= 1;
    }
    // peeled off last ntt stage for performance
    let mut i = 0;
    while i < (n << 1) {
        let omega_factor = element[(i + t) as usize];
        let omega_idx = (i >> 1) + n;
        let omega = root_of_unity[omega_idx as usize];
        let precon_omega = precon_root_of_unity[omega_idx as usize];
        let omega_factor = mod_mul_fast_const(omega_factor, omega, modulus, precon_omega);
        let mut lo_val = (*element)[(i) as usize];
        let mut hi_val = lo_val + omega_factor;
        if hi_val >= modulus {
            hi_val -= modulus;
        }
        if lo_val < omega_factor {
            lo_val += modulus;
        }
        lo_val -= omega_factor;
        (*element)[(i + 0) as usize] = hi_val;
        (*element)[(i + 1) as usize] = lo_val;
        i += 2;
    }
    // println!("Element vector: {:?}", element);
    // println!("-------------");
}


#[test]
fn test_crt_small() {
    let mut rou_data = [1u64, 4, 2, 8];
    let mut prou_data = [1085102592571150095u64, 4340410370284600380, 2170205185142300190, 8680820740569200760];
    let mut element_data = [2, 1, 3, 5];
    let root_of_unity = FheVec { modulus: 17, data:  &mut rou_data};
    let precon_root_of_unity = FheVec { modulus: 17, data: &mut prou_data };
    let mut element = FheVec { modulus: 17, data: &mut element_data };
    ntt_transform(&root_of_unity, &precon_root_of_unity, &mut element);
    assert_eq!(element.data, [5,6,8,6]);

}
