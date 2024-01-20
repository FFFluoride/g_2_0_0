//! Library made for ppputilising geometric algebra with only two basis vectors each squaring to 1.

use std::ops::Add;
use std::ops::{Mul, MulAssign, Neg, Sub};

use std::cmp::PartialOrd;

use num_traits::{identities::One, Zero};

use radians::{Angle, Float, Radians};

/// In 2D GA, Multivectors have 4 components
#[derive(Clone, Debug, Copy, PartialEq, PartialOrd, Hash, Eq, Ord)]
pub struct MultiVector<S> {
    pub scalar: S,
    pub e1: S,
    pub e2: S,
    pub e12: S,
}

impl<S> Add for MultiVector<S>
where
    S: Add<Output = S>,
{
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            scalar: self.scalar + other.scalar,
            e1: self.e1 + other.e1,
            e2: self.e2 + other.e2,
            e12: self.e12 + other.e12,
        }
    }
}

impl<S> Neg for MultiVector<S>
where
    S: Neg<Output = S>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            scalar: -self.scalar,
            e1: -self.e1,
            e2: -self.e2,
            e12: -self.e12,
        }
    }
}

/// Most of these functions have an in place mutation counterpart, they have "_mut" at the end of it
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

    // todo: Mul implementation using the geometric project

    pub fn scale(self, scalar: S) -> Self {
        //! Fun fact: Multivectors are vectors because you can scale and add them
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

    /// Grade projection is when you take only the elements of a multivector that is a certain grade. Scalars (or 0-Vectors) have grade 0. Vectors (or 1-Vectors) have grade 1. Bivectors (or 2-Vectors) have grade 2. In general K-Vectors have grade K.
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

    /// This operation is distributive across addittion and involves reversing the order of products. e.g. reverse(abcd) = dcba. The reverse predictably changes the size of specifically graded elements. Grade => Swap | No Swap. 0 => No Swap. 1 => No Swap. 2 => Swap. 3 => Swap. 4 => No Swap. 5 => No Swap &c.
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

    /// The Geometric product is distributive across addittion and is defined as: uv = 1/2(u.v + u^v). If you apply this rule and several others, in two dimensions, you get this formula for calculating the geometric product for arbitrary Multivectors in two dimensions. The Geometric product for vectors is particularly nice as it commutes when the vectors are parallel and anticommutes if they are perpendicular (this is more useful theoretically)
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
            e1: (a * f) + (b * e) - (c * h) + (d * g),
            e2: (a * g) + (b * h) + (c * e) - (d * f),
            e12: (a * h) + (b * g) - (c * f) + (d * e),
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
	*self = Self {
            scalar: (a * e) + (b * f) + (c * g) - (d * h),
            e1: (a * f) + (b * e) - (c * h) + (d * g),
            e2: (a * g) + (b * h) + (c * e) - (d * f),
            e12: (a * h) + (b * g) - (c * f) + (d * e),
        }
    }
    /// This is the geometic product where you multiply on the right and divide by the the magnitude squared of rhs
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
    /// Multiply by rhs on the left instead of right then account for coefficient
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
    /// The dual multiplied by the the original gives the psuedoscalar. Therefore simply do geometric division by the psuedoscalar
    pub fn dual(self) -> Self {
        let ps = Self {
            scalar: S::zero(),
            e1: S::zero(),
            e2: S::zero(),
            e12: S::one(),
        };
        ps.geometric_division_left(self)
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
    /// Acquiring what multiplies by the input multivector to produce the psuedoscalar. Theoretically self.dual().inverse() = self
    pub fn inverse_dual(self) -> Self {
        let ps = Self {
            scalar: S::zero(),
            e1: S::zero(),
            e2: S::zero(),
            e12: S::one(),
        };
        ps.geometric_division_left(self)
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
    /// Gets the maximum grade element that =/= 0. There may be some imprecise floating point errors
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
    /// The outer product is zero if the multivectors are parrallel across some space. If they are perpendicular it joins the two into a higher grade subspace. In general, the outer product of a j-vector and a k-vector the grade projection of the geometric product of the two into the j+k grade vector.
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
    /// Contrary to the outer product, the product of the two vectors is projected into the absolute difference of the grades of the two multivectors
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
    /// Magnitude is sometimes negative in GA so the magintude sqaured is more useful. This function doesn't have a mut version because I think that would be weird.
    pub fn magnitude_squared(self) -> S {
        self.geometric_product(self.reverse()).scalar
    }
}

/// The exponential of the 2d psuedoscalar squares to -1 so it is also i. This means we can use eulers formula. Also if you multiply a vector on the left and right by the exponential of angle in 2d, rotates the vector anticlockwise by double the angle.
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
    /// Generates a random multivector for testing
    fn rand_mvec() -> MultiVector<f32> {
        let mut rng = thread_rng();
        MultiVector {
            scalar: rng.gen(),
            e1: rng.gen(),
            e2: rng.gen(),
            e12: rng.gen(),
        }
    }
    /// Testing for parity between the mut and non-mut versions of functions
    #[cfg(test)]
    mod parity {
	use super::rand_mvec;
	use rand::prelude::*;
	 
        #[test]
        fn scale_test() {
	    let m = rand_mvec();
	    let mut m_clone = m.clone();
	    
	    let mut rng = thread_rng();
	    let scalar: f32 = rng.gen();
	    
	    m_clone.scale_mut(scalar);
	    assert_eq!(m.scale(scalar), m_clone)
	}

	#[test]
	fn project_test() {
	    for scalar in 0..2 {
		let m = rand_mvec();
		let mut m_clone = m.clone();
		
		m_clone.project_mut(scalar);
		assert_eq!(m.project(scalar), m_clone);
	    }
	    let scalar: usize = thread_rng().gen();
	    let m = rand_mvec();
	    let mut m_clone = m.clone();
	    m_clone.project_mut(scalar);
	    assert_eq!(m.project(scalar), m_clone);
	}

	#[test]
	fn reverse_test() {
	    let m = rand_mvec();
	    let mut m_clone = m.clone();
	    m_clone.reverse_mut();
	    assert_eq!(m.reverse(), m_clone);
	}

	#[test]
	fn geometric_product_test() {
	    let u = rand_mvec();
	    let mut u_clone = u.clone();
	    
	    let v = rand_mvec();

	    u_clone.geometric_product_mut(v);
	    assert_eq!(u.geometric_product(v), u_clone);
	}

	#[test]
	fn geometric_division_test() {
	    let u = rand_mvec();
	    let mut u_clone = u.clone();
	    
	    let v = rand_mvec();

	    u_clone.geometric_division_mut(v);
	    assert_eq!(u.geometric_division(v), u_clone);
	}

	#[test]
	fn geometric_division_left_test() {
	    let u = rand_mvec();
	    let mut u_clone = u.clone();
	    
	    let v = rand_mvec();

	    v.geometric_division_left_mut(&mut u_clone);
	    assert_eq!(v.geometric_division_left(u), u_clone);
	}

	#[test]
	fn dual_test() {
	    let u = rand_mvec();
	    let mut u_clone = u.clone();
	    
	    u_clone.dual_mut();
	    assert_eq!(u.dual(), u_clone);
	}

	#[test]
	fn inverse_dual_test() {
	    let u = rand_mvec();
	    let mut u_clone = u.clone();
	    
	    u_clone.inverse_dual_mut();
	    assert_eq!(u.inverse_dual(), u_clone);
	}

	#[test]
	fn outer_product_test() {
	    let u = rand_mvec();
	    let mut u_clone = u.clone();
	    
	    let v = rand_mvec();

	    u_clone.outer_product_mut(v);
	    assert_eq!(u.outer_product(v), u_clone);
	}

	#[test]
	fn inner_product_test() {
	    let u = rand_mvec();
	    let mut u_clone = u.clone();
	    
	    let v = rand_mvec();

	    u_clone.inner_product_mut(v);
	    assert_eq!(u.inner_product(v), u_clone);
	}
    }


    #[cfg(test)]
    mod logic {
        //! Testing for consistent results of operations
        use super::rand_mvec;

        #[test]
        fn reverse_twice() {
            let m = rand_mvec();
            assert_eq!(m, m.reverse().reverse());
        }

        #[test]
        fn outer_product_anticommutes() {
            let u = rand_mvec();
            let v = rand_mvec();
            assert_eq!(u.outer_product(v), -v.outer_product(u));
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
	/*
	#[test]
	fn dual_inverse_dual_is_original() {
	    let u = rand_mvec();
	    assert_eq!(u.dual().inverse_dual(), u);
	}*/
    }
}
