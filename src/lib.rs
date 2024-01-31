//! Library made for utilising geometric algebra with only two basis vectors each squaring to 1 or G(2, 0, 0) in common notation.
//! The main object here is a multivector, however, there is also an exponential function allow the generation of Rotors in 2d if you prefer.
//! This is a personal project of mine but I published it on crates in case someone finds it interesting.
//! If you don't know about geometric algebra, check out this swift introduction: <https://youtu.be/60z_hpEAtD8?si=3R4PWa4RDw0soUXw>
//! Also, vga refers to Vanilla Geometric Algebra

use std::ops::{MulAssign, AddAssign};
use std::ops::{Div, Mul, Neg, Sub, Add};

use std::cmp::PartialOrd;

use num_traits::{identities::One, Float as Flt, Pow, Zero};

use radians::{Angle, Float, Radians};

/// In 2D vga, Multivectors have 4 components
#[derive(Clone, Debug, Copy, PartialEq, PartialOrd, Hash, Eq, Ord)]
pub struct MultiVector<S> {
    pub scalar: S,
    pub e1: S,
    pub e2: S,
    pub e12: S,
}

impl<S:AddAssign> AddAssign for MultiVector<S>{
    fn add_assign(&mut self, other: Self) {
	self.scalar += other.scalar;
	self.e1 += other.e1;
	self.e2 += other.e2;
	self.e12 += other.e12;
    }
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

impl<S> Sub for MultiVector<S>
where
    S: Sub<Output = S>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            scalar: self.scalar - rhs.scalar,
            e1: self.e1 - rhs.e1,
            e2: self.e2 - rhs.e2,
            e12: self.e12 - rhs.e12,
        }
    }
}

/// Most of these functions have an in place mutation counterpart, they have "_mut" at the end of it
impl<S> MultiVector<S>
where
    S: MulAssign
        + Mul<Output = S>
        + Neg<Output = S>
        + Div<Output = S>
        + Copy
        + Zero
        + Sub<Output = S>
        + One
        + Pow<u8, Output = S>
        + Flt
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

    /// Convenient function for initialising an empty multivector
    pub fn zero() -> Self {
	MultiVector::new(S::zero(), S::zero(), S::zero(), S::zero())
    }

    /// Convenient function to access the psuedoscalar in 2d vga
    pub fn ps() -> Self {
        Self {
            scalar: S::zero(),
            e1: S::zero(),
            e2: S::zero(),
            e12: S::one(),
        }
    }

    // todo: Mul implementation using the geometric project

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

    /// Grade projection. This function returns only the elements of the given grade.
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

    /// The reverse only negates the e12 component in 2d vga.
    pub fn reverse(self) -> Self {
        Self {
            scalar: self.scalar,
            e1: self.e1,
            e2: self.e2,
            e12: -self.e12,
        }
    }
    pub fn reverse_mut(&mut self) {
        self.e12 = -self.e12
    }

    /// Hack for normalizing a vector
    pub fn normalize_vector_part(self) -> Self {
	self.project(1).scale(S::one() / (self.project(1).magnitude_squared().sqrt()))
    }

    pub fn normalize_vector_part_mut(&mut self) {
	self.project_mut(1);
	self.scale_mut(S::one() / (self.magnitude_squared().sqrt()));
    }

    /// The Geometric product is distributive across addittion and is defined as: uv = 1/2(u.v + u^v), where u and v are vectors. The Geometric product for vectors commutes when the vectors are parallel and anticommutes when they are perpendicular
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
    /// Convenience function for the geometric product where the right hand side is a mutable reference
    pub fn geometric_product_swapped_mut(self, rhs: &mut Self) {
        let a = self.scalar;
        let b = self.e1;
        let c = self.e2;
        let d = self.e12;
        let e = rhs.scalar;
        let f = rhs.e1;
        let g = rhs.e2;
        let h = rhs.e12;
        *rhs = Self {
            scalar: (a * e) + (b * f) + (c * g) - (d * h),
            e1: (a * f) + (b * e) - (c * h) + (d * g),
            e2: (a * g) + (b * h) + (c * e) - (d * f),
            e12: (a * h) + (b * g) - (c * f) + (d * e),
        }
    }

    /// Gets the maximum grade element above 0. There may be some imprecise floating point errors
    pub fn grade(&self) -> usize {
        let mut grade = 0;
        if self.e1 > S::zero() || self.e2 > S::zero() {
            grade += 1;
        }
        if self.e2 > S::zero() {
            grade += 1;
        }
        grade
    }

    /// Helpful function that gets the scalar and e12 component of the multivector; getting the even subspace entails using the `.project(1)` function, this also makes the multivector a versor.
    pub fn odd_subspace(&self) -> Self {
	self.project(0).project(2)
    }

    pub fn odd_subspace_mut(&mut self) {
	self.project_mut(0);
	self.project_mut(2);
    }
    
    /// The outer product of a k-vector and a j-vector is a j+k-vector
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
    /// The inner product of a k-vector and a j-vector is a |j-k|-vector
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

    pub fn versor_inverse(self) -> Self {
        self.reverse()
            .scale(S::one() / self.geometric_product(self.reverse()).scalar)
    }

    /// A versor is the product of an arbitrary amount of vectors.
    /// There is no arbitrary multivector inverse for 2d vga, however, there is a n inverse for versor in general
    pub fn versor_inverse_mut(&mut self) {
        let self_clone = *self;
        self.reverse_mut();
        self.scale_mut(S::one() / self_clone.geometric_product(self_clone.reverse()).scalar);
    }

    /// a * dual(a) = e12.
    /// The dual requires the multivector to have an inverse and not all 2d vga multivectors have an inverse but versors do, so this function only works for versors. For this function to work correctly, you must pass in a versor
    /// A versors will either come in vector for or scalar + e12 form
    pub fn versor_dual(self) -> Self {
        self.versor_inverse().geometric_product(Self::ps())
    }

    pub fn versor_dual_mut(&mut self) {
        self.versor_inverse_mut();
        self.geometric_product_mut(Self::ps());
    }

    pub fn versor_inverse_dual(self) -> Self {
        Self::ps().geometric_product(self.versor_inverse())
    }

    pub fn versor_inverse_dual_mut(&mut self) {
        self.versor_inverse_mut();
        Self::ps().geometric_product_swapped_mut(self);
    }

    /// Magnitude is sometimes negative in GA so the magintude sqaured is more useful. This function doesn't have a mut version because I think that would be weird.
    pub fn magnitude_squared(self) -> S {
        self.geometric_product(self.reverse()).scalar
    }

    pub fn round_components(self, dp: u8) -> Self {
        let ten = S::one()
            + S::one()
            + S::one()
            + S::one()
            + S::one()
            + S::one()
            + S::one()
            + S::one()
            + S::one()
            + S::one();
        Self {
            scalar: (self.scalar * S::pow(ten, dp)).round() / S::pow(ten, dp),
            e1: (self.e1 * S::pow(ten, dp)).round() / S::pow(ten, dp),
            e2: (self.e2 * S::pow(ten, dp)).round() / S::pow(ten, dp),
            e12: (self.e12 * S::pow(ten, dp)).round() / S::pow(ten, dp),
        }
    }
}

/// The exponential of the 2d psuedoscalar squares to -1 so it is also i. This means we can use eulers formula. Also if you multiply a vector on the left and right by the exponential of angle in 2d, rotates the vector anticlockwise by double the angle.
pub fn exp<F, S>(theta: Angle<F, Radians>) -> MultiVector<S>
where
    F: Float,
    S: MulAssign
        + Mul<Output = S>
        + Neg<Output = S>
        + Div<Output = S>
        + Copy
        + Zero
        + Sub<Output = S>
        + One
        + From<F>
        + Flt
        + Pow<u8, Output = S>
        + PartialOrd,
{
    let (sin, cos) = theta.sin_cos();
    MultiVector::<S>::new(cos.into(), S::zero(), S::zero(), sin.into())
}

#[cfg(test)]
mod tests {
    use super::*;
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

    fn rand_versor() -> MultiVector<f32> {
        rand_mvec().project(1)
    }

    fn fair_versor() -> MultiVector<f32> {
        if rand::random() {
            rand_mvec().project(1)
        } else {
            rand_mvec()
                .project(1)
                .geometric_product(rand_mvec().project(1))
        }
    }

    /// Testing for parity between the mut and non-mut versions of functions
    #[cfg(test)]
    mod parity {
        use super::*;

        #[test]
        fn scale_test() {
            for _ in 0..3 {
                let m = rand_mvec();
                let mut m_clone = m;

                let mut rng = thread_rng();
                let scalar: f32 = rng.gen();

                m_clone.scale_mut(scalar);
                assert_eq!(m.scale(scalar), m_clone)
            }
        }

        #[test]
        fn project_test() {
            for _ in 0..3 {
                for scalar in 0..2 {
                    let m = rand_mvec();
                    let mut m_clone = m;

                    m_clone.project_mut(scalar);
                    assert_eq!(m.project(scalar), m_clone);
                }
                let scalar: usize = thread_rng().gen();
                let m = rand_mvec();
                let mut m_clone = m;
                m_clone.project_mut(scalar);
                assert_eq!(m.project(scalar), m_clone);
            }
        }

        #[test]
        fn reverse_test() {
            for _ in 0..3 {
                let m = rand_mvec();
                let mut m_clone = m;
                m_clone.reverse_mut();
                assert_eq!(m.reverse(), m_clone);
            }
        }

        #[test]
        fn geometric_product_test() {
            for _ in 0..3 {
                let u = rand_mvec();
                let mut u_clone = u;

                let v = rand_mvec();

                u_clone.geometric_product_mut(v);
                assert_eq!(u.geometric_product(v), u_clone);
            }
        }

        #[test]
        fn outer_product_test() {
            for _ in 0..3 {
                let u = rand_mvec();
                let mut u_clone = u;

                let v = rand_mvec();

                u_clone.outer_product_mut(v);
                assert_eq!(u.outer_product(v), u_clone);
            }
        }

        #[test]
        fn inner_product_test() {
            for _ in 0..3 {
                let u = rand_mvec();
                let mut u_clone = u;

                let v = rand_mvec();

                u_clone.inner_product_mut(v);
                assert_eq!(u.inner_product(v), u_clone);
            }
        }

        #[test]
        fn versor_dual_test() {
            for _ in 0..3 {
                let u = fair_versor();

                let mut u_clone = u;
                u_clone.versor_dual_mut();

                assert_eq!(u.versor_dual(), u_clone);
            }
        }

        #[test]
        fn geometric_product_swapped_test() {
            for _ in 0..3 {
                let u = rand_mvec();
                let v = rand_mvec();

                let mut u_clone = u;
                let mut v_clone = v;
                u.geometric_product_swapped_mut(&mut v_clone);
                u_clone.geometric_product_mut(v);

                assert_eq!(u_clone, v_clone);
            }
        }

        #[test]
        fn versor_inverse_dual_test() {
            for _ in 0..3 {
                let u = fair_versor();

                let mut u_clone = u;
                u_clone.versor_inverse_dual_mut();

                assert_eq!(u.versor_inverse_dual(), u_clone);
            }
        }

	#[test]
	fn normalize_vector_part_test() {
	    for _ in 0..3 {
		let u = rand_mvec();
		let mut u_clone = u.clone();
		u_clone.normalize_vector_part_mut();

		assert_eq!(u.normalize_vector_part(), u_clone);
	    }
	}
    }

    #[cfg(test)]
    mod logic {
        //! Testing for consistent results of operations
        use super::*;

        #[test]
        fn reverse_twice() {
            for _ in 0..3 {
                let m = rand_mvec();
                assert_eq!(m, m.reverse().reverse());
            }
        }

        #[test]
        fn outer_product_anticommutes() {
            for _ in 0..3 {
                let u = rand_mvec();
                let v = rand_mvec();
                assert_eq!(u.outer_product(v), -v.outer_product(u));
            }
        }

        #[test]
        fn inner_product_commutes() {
            for _ in 0..3 {
                let u = rand_mvec();
                let v = rand_mvec();
                assert_eq!(u.inner_product(v), v.inner_product(u));
            }
        }

        #[test]
        fn sum_of_projections_is_self() {
            for _ in 0..3 {
                let u = rand_mvec();
                let u_0 = u.project(0);
                let u_1 = u.project(1);
                let u_2 = u.project(2);
                assert_eq!(u_0 + u_1 + u_2, u);
            }
        }

        #[test]
        fn versor_inverse_odd() {
            for _ in 0..3 {
                let u = fair_versor();
                assert_eq!(
                    u.geometric_product(u.versor_inverse())
                        .round_components(5)
                        .scalar,
                    1f32
                );
            }
        }

        #[test]
        fn versor_inverse_even() {
            for _ in 0..3 {
                let prod = fair_versor().geometric_product(rand_versor());
                assert_eq!(
                    prod.geometric_product(prod.versor_inverse())
                        .round_components(5)
                        .scalar,
                    1f32
                );
            }
        }

        #[test]
        fn versor_dual_inverse_dual_is_original() {
            for _ in 0..3 {
                let u = fair_versor();
                assert_eq!(
                    u.versor_dual().versor_inverse_dual().round_components(5),
                    u.round_components(5)
                );
            }
        }
    }
}
