#![cfg(feature = "zeroize")]

use super::BigUint;

impl zeroize::Zeroize for BigUint {
    fn zeroize(&mut self) {
        self.data.zeroize();
    }
}
