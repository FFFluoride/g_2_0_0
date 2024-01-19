use std::ops::{Mul, MulAssign, Neg, Sub};

use std::cmp::PartialOrd;

use num_traits::{identities::One, Zero};

use radians::{Angle, Float, Radians};

#[derive(Clone, Debug, Copy, PartialEq, PartialOrd, Hash, Eq, Ord)]
pub struct MultiVector<S> {
    pub scalar: S,
    pub e1: S,
    pub e2: S,
    pub e12: S,
}

impl<S> MultiVector<S>
where
    S: MulAssign
        + Mul<Output = S>
        + Neg<Output = S>
        + Copy
        + Zero
        + Sub<Output = S>
        + One
        + From<i16>
        + PartialOrd,
{
    pub fn new(scalar: S, e1: S, e2: S, e12: S) -> Self {
        Self {
            scalar,
            e1,
            e2,
            e12,
        }
    }

    pub fn scale(self, scalar: S) -> Self {
        Self {
            scalar: self.scalar * scalar,
            e1: self.e1 * scalar,
            e2: self.e2 * scalar,
            e12: self.e12 * scalar,
        }
    }
    pub fn scale_mut(&mut self, scalar: S) {
        self.scalar *= scalar;
        self.e1 *= scalar;
        self.e2 *= scalar;
        self.e12 *= scalar;
    }
    pub fn flip_sign(self) -> Self {
        Self {
            scalar: -self.scalar,
            e1: -self.e1,
            e2: -self.e2,
            e12: -self.e12,
        }
    }
    pub fn flip_sign_mut(&mut self) {
        self.scalar = self.scalar.neg();
        self.e1 = self.e1.neg();
        self.e2 = self.e2.neg();
        self.e12 = self.e12.neg();
    }
    pub fn project(self, grade: usize) -> Self {
        match grade {
            0 => Self {
                scalar: self.scalar,
                e1: S::zero(),
                e2: S::zero(),
                e12: S::zero(),
            },
            1 => Self {
                scalar: S::zero(),
                e1: self.e1,
                e2: self.e2,
                e12: S::zero(),
            },
            2 => Self {
                scalar: S::zero(),
                e1: S::zero(),
                e2: S::zero(),
                e12: self.e12,
            },
            _ => Self {
                scalar: S::zero(),
                e1: S::zero(),
                e2: S::zero(),
                e12: S::zero(),
            },
        }
    }
    pub fn project_mut(&mut self, grade: usize) {
        match grade {
            0 => {
                self.e1 = S::zero();
                self.e2 = S::zero();
                self.e12 = S::zero();
            }
            1 => {
                self.scalar = S::zero();
                self.e12 = S::zero();
            }
            2 => {
                self.scalar = S::zero();
                self.e1 = S::zero();
                self.e2 = S::zero();
            }
            _ => {
                self.scalar = S::zero();
                self.e1 = S::zero();
                self.e2 = S::zero();
                self.e12 = S::zero();
            }
        }
    }
    pub fn reverse(self) -> Self {
        Self {
            scalar: self.scalar,
            e1: self.e1,
            e2: self.e2,
            e12: -self.e12,
        }
    }
    pub fn reverse_mut(&mut self) {
        self.e12 = self.e12.neg()
    }
    pub fn geometric_product(self, rhs: Self) -> Self {
        let a = self.scalar;
        let b = self.e1;
        let c = self.e2;
        let d = self.e12;
        let e = rhs.scalar;
        let f = rhs.e1;
        let g = rhs.e2;
        let h = rhs.e12;
        Self {
            scalar: (a * e) + (b * f) + (c * g) - (d * h),
            e1: (a * f) + (b * e) + (d * g) - (c * h),
            e2: (a * g) + (c * e) + (b * h) - (d * f),
            e12: (a * h) + (b * g) + (d * e) - (c * f),
        }
    }
    pub fn geometric_product_mut(&mut self, rhs: Self) {
        let a = self.scalar;
        let b = self.e1;
        let c = self.e2;
        let d = self.e12;
        let e = rhs.scalar;
        let f = rhs.e1;
        let g = rhs.e2;
        let h = rhs.e12;
        self.scalar = (a * e) + (b * f) + (c * g) - (d * h);
        self.e1 = (a * f) + (b * e) + (d * g) - (c * h);
        self.e2 = (a * g) + (c * e) + (b * h) - (d * f);
        self.e12 = (a * h) + (b * g) + (d * e) - (c * f);
    }
    pub fn geometric_division(mut self, rhs: Self) -> Self {
        self.scalar = self.scalar.neg();
        self.scale_mut(
            (self.scalar * self.scalar) + (self.e1 * self.e1) + (self.e2 * self.e2)
                - (self.e12 * self.e12),
        );
        self.geometric_product(rhs)
    }
    pub fn geometric_division_mut(&mut self, rhs: Self) {
        self.scalar = self.scalar.neg();
        self.scale_mut(
            (self.scalar * self.scalar) + (self.e1 * self.e1) + (self.e2 * self.e2)
                - (self.e12 * self.e12),
        );
        self.geometric_product_mut(rhs);
    }

    pub fn geometric_division_left(self, mut rhs: Self) -> Self {
        rhs.scalar = rhs.scalar.neg();
        rhs.scale_mut(
            (rhs.scalar * rhs.scalar) + (rhs.e1 * rhs.e1) + (rhs.e2 * rhs.e2)
                - (self.e12 * self.e12),
        );
        rhs.geometric_product(self)
    }

    pub fn geometric_division_left_mut(self, rhs: &mut Self) {
        rhs.scalar = rhs.scalar.neg();
        rhs.scale_mut(
            (rhs.scalar * rhs.scalar) + (rhs.e1 * rhs.e1) + (rhs.e2 * rhs.e2)
                - (self.e12 * self.e12),
        );
        rhs.geometric_product_mut(self)
    }

    pub fn dual(self) -> Self {
        let ps = Self {
            scalar: S::zero(),
            e1: S::zero(),
            e2: S::zero(),
            e12: S::one(),
        };
        ps.geometric_division(self)
    }
    pub fn dual_mut(&mut self) {
        let ps = Self {
            scalar: S::zero(),
            e1: S::zero(),
            e2: S::zero(),
            e12: S::one(),
        };
        ps.geometric_division_left_mut(self);
    }
    pub fn inverse_dual(self) -> Self {
        let ps = Self {
            scalar: S::zero(),
            e1: S::zero(),
            e2: S::zero(),
            e12: S::one(),
        };
        ps.geometric_division(self)
    }
    pub fn inverse_dual_mut(&mut self) {
        let ps = Self {
            scalar: S::zero(),
            e1: S::zero(),
            e2: S::zero(),
            e12: S::one(),
        };
        ps.geometric_division_left_mut(self)
    }
    pub fn grade(&self) -> usize {
        let mut grade = 0;
        if self.e1 > 0.into() || self.e2 > 0.into() {
            grade += 1;
        }
        if self.e2 > 0.into() {
            grade += 1;
        }
        grade
    }
    pub fn outer_product(self, rhs: Self) -> Self {
        let grade1 = self.grade();
        let grade2 = rhs.grade();
        self.geometric_product(rhs).project(grade1 + grade2)
    }
    pub fn outer_product_mut(&mut self, rhs: Self) {
        let grade1 = self.grade();
        let grade2 = rhs.grade();
        self.geometric_product_mut(rhs);
        self.project_mut(grade1 + grade2)
    }
    pub fn inner_product(self, rhs: Self) -> Self {
        let grade1 = self.grade();
        let grade2 = rhs.grade();
        self.geometric_product(rhs).project(grade1.abs_diff(grade2))
    }
    pub fn inner_product_mut(&mut self, rhs: Self) {
        let grade1 = self.grade();
        let grade2 = rhs.grade();
        self.geometric_product_mut(rhs);
        self.project_mut(grade1.abs_diff(grade2))
    }
    pub fn magnitude_squared(self) -> S {
        self.geometric_product(self.reverse()).scalar
    }
}

pub fn exp<F, S>(theta: Angle<F, Radians>) -> MultiVector<S>
where
    F: Float,
    S: MulAssign
        + Mul<Output = S>
        + Neg<Output = S>
        + Copy
        + Zero
        + Sub<Output = S>
        + One
        + From<F>
        + From<i16>
        + PartialOrd,
{
    let (sin, cos) = theta.sin_cos();
    MultiVector::<S>::new(cos.into(), S::zero(), S::zero(), sin.into())
}

#[cfg(test)]
mod tests {
    use super::MultiVector;
    use rand::prelude::*;

    fn rand_mvec() -> MultiVector<f32> {
        let mut rng = thread_rng();
        MultiVector {
            scalar: rng.gen(),
            e1: rng.gen(),
            e2: rng.gen(),
            e12: rng.gen(),
        }
    }

    #[test]
    fn reverse_twice() {
        let m = rand_mvec();
        assert_eq!(m, m.reverse().reverse());
    }

    #[test]
    fn outer_product_anticommutes() {
        let u = rand_mvec();
        let v = rand_mvec();
        assert_eq!(u.outer_product(v), v.outer_product(u).scale(-1f32));
    }

    #[test]
    fn inner_product_commutes() {
	let u = rand_mvec();
	let v = rand_mvec();
	assert_eq!(u.inner_product(v), v.inner_product(u));
    }

    #[test]
    fn inner_equals_inner_mut() {
	let u = rand_mvec();
	let v = rand_mvec();
	let mut u_clone = u.clone();
	let mut v_clone = v.clone();
	u_clone.inner_product_mut(v);
	v_clone.inner_product_mut(u);
	assert_eq!(u_clone, v_clone);
    }
}
