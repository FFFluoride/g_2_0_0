pub use std::ops::{Mul, MulAssign, Neg, Sub};

pub use std::cmp::PartialOrd;

pub use num_traits::{
    Zero,
    identities::One,
};

#[derive(Clone)]
pub struct MultiVector<S> {
    scalar: S,
    e1: S,
    e2: S,
    e12: S,
}

impl<S> MultiVector<S>
where
    S: MulAssign + Mul<Output = S> + Neg<Output = S> + Copy + Zero + Sub<Output = S> + One + PartialOrd<isize>,
{
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
        Self {
            scalar: (self.scalar * rhs.scalar) + (self.e1 * rhs.e1) + (self.e2 * rhs.e2)
                - (self.e12 * rhs.e12),
            e1: (self.scalar * rhs.e1) + (self.e1 * rhs.scalar) + (self.e12 * rhs.e2)
                - (self.e2 * rhs.e12),
            e2: (self.scalar * rhs.e2) + (self.e2 * rhs.scalar) + (self.e1 * rhs.e12)
                - (self.e12 * rhs.e1),
            e12: (self.e1 * rhs.e2) - (self.e2 * rhs.e1),
        }
    }
    pub fn geometric_product_mut(&mut self, rhs: Self) {
        self.scalar = (self.scalar * rhs.scalar) + (self.e1 * rhs.e1) + (self.e2 * rhs.e2)
            - (self.e12 * rhs.e12);
        self.e1 = (self.scalar * rhs.e1) + (self.e1 * rhs.scalar) + (self.e12 * rhs.e2)
            - (self.e2 * rhs.e12);
        self.e2 = (self.scalar * rhs.e2) + (self.e2 * rhs.scalar) + (self.e1 * rhs.e12)
            - (self.e12 * rhs.e1);
        self.e12 = (self.e1 * rhs.e2) - (self.e2 * rhs.e1);
    }
    pub fn geometric_division(mut self, rhs: Self) -> Self {
	self.scalar = self.scalar.neg();
	self.scale_mut((self.scalar * self.scalar) + (self.e1 * self.e1) + (self.e2 * self.e2) - (self.e12 * self.e12));
	self.geometric_product(rhs)
    }
    pub fn geometric_division_mut(&mut self, rhs: Self) {
	self.scalar = self.scalar.neg();
	self.scale_mut((self.scalar * self.scalar) + (self.e1 * self.e1) + (self.e2 * self.e2) - (self.e12 * self.e12));
	self.geometric_product_mut(rhs);
    }

    pub fn geometric_division_left(self, mut rhs: Self) -> Self {
	rhs.scalar = rhs.scalar.neg();
	rhs.scale_mut((rhs.scalar * rhs.scalar) + (rhs.e1 * rhs.e1) + (rhs.e2 * rhs.e2) - (self.e12 * self.e12));
	rhs.geometric_product(self)
    }

    pub fn geometric_division_left_mut(self, rhs: &mut Self) {
	rhs.scalar = rhs.scalar.neg();
	rhs.scale_mut((rhs.scalar * rhs.scalar) + (rhs.e1 * rhs.e1) + (rhs.e2 * rhs.e2) - (self.e12 * self.e12));
	rhs.geometric_product_mut(self)
    }
    
    pub fn dual(self) -> Self {
	let ps = Self { scalar: S::zero(), e1: S::zero(), e2: S::zero(), e12: S::one() };
	ps.geometric_division(self)
    }
    pub fn dual_mut(&mut self) {
	let ps = Self { scalar: S::zero(), e1: S::zero(), e2: S::zero(), e12: S::one() };
	ps.geometric_division_left_mut(self);
    }
    pub fn inverse_dual(self) -> Self {
	let ps = Self { scalar: S::zero(), e1: S::zero(), e2: S::zero(), e12: S::one() };
	ps.geometric_division(self)
    }
    pub fn inverse_dual_mut(&mut self) {
	let ps = Self { scalar: S::zero(), e1: S::zero(), e2: S::zero(), e12: S::one() };
	ps.geometric_division_left_mut(self)
    }
    pub fn grade(&self) -> usize {
	let mut grade = 0;
	if self.e1 > 0 || self.e2 > 0 {
	    grade += 1;
	}
	if self.e2 > 0 {
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
	self.clone().geometric_product(self.reverse()).scalar
    }
}
