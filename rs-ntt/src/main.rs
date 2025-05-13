use crate::ntt::{FheVec, ntt_transform_hybrid, ntt_transform};

pub mod ntt;

fn main() {
    println!("Start");
    let mut rou_data = [1u64, 4, 2, 8];
    let mut prou_data = [1085102592571150095u64, 4340410370284600380, 2170205185142300190, 8680820740569200760];
    let mut element_data = [2, 1, 3, 5];
    let root_of_unity = FheVec { modulus: 17, data:  &mut rou_data};
    let precon_root_of_unity = FheVec { modulus: 17, data: &mut prou_data };
    let mut element = FheVec { modulus: 17, data: &mut element_data };
    //ntt_transform(&root_of_unity, &precon_root_of_unity, &mut element);
    ntt_transform_hybrid(&root_of_unity, &precon_root_of_unity, &mut element);
    assert_eq!(element.data, [5,6,8,6]);
}