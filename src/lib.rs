//! Transistor — artefact théorique (F17, état matriciel, φ, et tour AES-like).
//! - F17: arithmétique mod 17
//! - FSMState<N>: état concret en F17
//! - φ: extraction (4,6,12,14) en row-major
//! - Digit17 + Mat17<T,N>: version générique
//! - AES-like (SubBytes/ShiftRows/MixColumns/AddRoundKey) sur F17
//! - demo_aes_like(): affiche les coefficients avant/après

use core::fmt;

// ===============================
//  F_17 : arithmétique mod 17
// ===============================
#[derive(Copy, Clone, Eq, PartialEq, Default)]
pub struct F17(u8); // invariant : valeur dans [0,16]

impl F17 {
    pub const MOD: u8 = 17;

    #[inline]
    pub fn new(x: i32) -> Self {
        let mut r = x % (Self::MOD as i32);
        if r < 0 { r += Self::MOD as i32; }
        F17(r as u8)
    }

    #[inline] pub fn zero() -> Self { F17(0) }
    #[inline] pub fn one()  -> Self { F17(1) }
    #[inline] pub fn value(self) -> u8 { self.0 }

    #[inline] pub fn add(self, other: Self) -> Self { F17::new(self.0 as i32 + other.0 as i32) }
    #[inline] pub fn sub(self, other: Self) -> Self { F17::new(self.0 as i32 - other.0 as i32) }
    #[inline] pub fn mul(self, other: Self) -> Self { F17::new((self.0 as i32) * (other.0 as i32)) }

    /// Puissance mod 17 (exponentiation rapide).
    pub fn pow(self, mut e: u32) -> Self {
        let mut base = self;
        let mut acc = F17::one();
        while e > 0 {
            if e & 1 == 1 { acc = acc.mul(base); }
            base = base.mul(base);
            e >>= 1;
        }
        acc
    }
}

impl fmt::Debug for F17 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{} (mod 17)", self.0) }
}
impl fmt::Display for F17 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{}", self.0) }
}

// ============================================
//  État FSM concret : matrice N×N sur F17
// ============================================
#[derive(Clone)]
pub struct FSMState<const N: usize> {
    /// X[i][j] en indices 0-based.
    pub x: [[F17; N]; N],
}

impl<const N: usize> FSMState<N> {
    /// Construit X via un générateur f(i,j), sans exiger Copy/Default.
    pub fn from_fn<F: Fn(usize, usize) -> F17>(f: F) -> Self {
        let x: [[F17; N]; N] =
            core::array::from_fn(|i| core::array::from_fn(|j| f(i, j)));
        Self { x }
    }

    /// Convertit t (1-based, row-major) en (i,j) 0-based ; vérifie les bornes.
    #[inline]
    pub fn t_to_ij(t_1based: usize) -> (usize, usize) {
        assert!(t_1based >= 1 && t_1based <= N * N, "t hors bornes");
        let k = t_1based - 1;
        (k / N, k % N)
    }

    /// Accès à x^{lin}_t (t en 1-based).
    #[inline]
    pub fn x_lin(&self, t_1based: usize) -> F17 {
        let (i, j) = Self::t_to_ij(t_1based);
        self.x[i][j]
    }
}

/// Filtrage φ concret : (x_lin_4, x_lin_6, x_lin_12, x_lin_14).
#[inline]
pub fn phi<const N: usize>(state: &FSMState<N>) -> [F17; 4] {
    assert!(N * N >= 14, "m = N² doit être ≥ 14");
    [ state.x_lin(4), state.x_lin(6), state.x_lin(12), state.x_lin(14) ]
}

// =====================================================
//  Version générique : Digit17 + Mat17<T, N> + φ générique
// =====================================================
pub trait Digit17: Sized + Clone {
    fn add_mod(&self, rhs: &Self) -> Self;
    fn sub_mod(&self, rhs: &Self) -> Self;
    fn mul_mod(&self, rhs: &Self) -> Self;
}
impl Digit17 for F17 {
    #[inline] fn add_mod(&self, rhs: &Self) -> Self { self.add(*rhs) }
    #[inline] fn sub_mod(&self, rhs: &Self) -> Self { self.sub(*rhs) }
    #[inline] fn mul_mod(&self, rhs: &Self) -> Self { self.mul(*rhs) }
}

#[derive(Clone)]
pub struct Mat17<T, const N: usize> {
    pub x: [[T; N]; N],
}

impl<T: Digit17 + Clone, const N: usize> Mat17<T, N> {
    /// Construit la matrice via un générateur f(i,j).
    pub fn from_fn<F: Fn(usize, usize) -> T>(f: F) -> Self {
        let x: [[T; N]; N] =
            core::array::from_fn(|i| core::array::from_fn(|j| f(i, j)));
        Self { x }
    }

    /// t (1-based, row-major) → (i,j) 0-based.
    #[inline]
    pub fn t_to_ij(t_1based: usize) -> (usize, usize) {
        assert!(t_1based >= 1 && t_1based <= N * N, "t hors bornes");
        let k = t_1based - 1;
        (k / N, k % N)
    }

    /// Référence sur l’élément linéarisé.
    #[inline]
    pub fn lin(&self, t_1based: usize) -> &T {
        let (i, j) = Self::t_to_ij(t_1based);
        &self.x[i][j]
    }
}

/// φ générique sur Mat17<T, N>.
#[inline]
pub fn phi_generic<T: Digit17 + Clone, const N: usize>(m: &Mat17<T, N>) -> [T; 4] {
    assert!(N * N >= 14, "m = N² doit être ≥ 14");
    [ m.lin(4).clone(), m.lin(6).clone(), m.lin(12).clone(), m.lin(14).clone() ]
}

// ============================================
//  AES-like round sur F17 (N = 4)
// ============================================

/// S-box : S(x) = 5 * x^{-1} + 7 (mod 17) avec S(0) = 7.
#[inline]
pub fn sbox_f17(x: F17) -> F17 {
    let a = F17::new(5);
    let b = F17::new(7);
    if x.value() == 0 { return b; }
    let inv = x.pow(15); // x^{p-2} mod p pour p=17
    a.mul(inv).add(b)
}

/// SubBytes (cellule à cellule).
pub fn subbytes_f17<const N: usize>(state: &[[F17; N]; N]) -> [[F17; N]; N] {
    core::array::from_fn(|i| core::array::from_fn(|j| sbox_f17(state[i][j])))
}

/// Rotation gauche d'une ligne de longueur N de r positions.
#[inline]
fn rotl_row<const N: usize>(row: &[F17; N], r: usize) -> [F17; N] {
    let mut out = [F17::zero(); N];
    for j in 0..N { out[j] = row[(j + r) % N]; }
    out
}

/// ShiftRows à la AES : la ligne i est décalée de i positions à gauche.
pub fn shiftrows_f17<const N: usize>(state: &[[F17; N]; N]) -> [[F17; N]; N] {
    core::array::from_fn(|i| rotl_row(&state[i], i % N))
}

/// Matrice MixColumns 4x4 (inversible mod 17).
const MC4: [[F17; 4]; 4] = [
    [F17(2), F17(3), F17(1), F17(1)],
    [F17(1), F17(2), F17(3), F17(1)],
    [F17(1), F17(1), F17(2), F17(3)],
    [F17(3), F17(1), F17(1), F17(2)],
];

/// Multiplie une colonne v par MC4 (mod 17).
fn mix_one_col(v: [F17; 4]) -> [F17; 4] {
    let mut out = [F17::zero(); 4];
    for i in 0..4 {
        let mut acc = F17::zero();
        for j in 0..4 { acc = acc.add(MC4[i][j].mul(v[j])); }
        out[i] = acc;
    }
    out
}

/// MixColumns (N=4).
pub fn mixcolumns_f17(state: &[[F17; 4]; 4]) -> [[F17; 4]; 4] {
    let mut out = [[F17::zero(); 4]; 4];
    for col in 0..4 {
        let v = [state[0][col], state[1][col], state[2][col], state[3][col]];
        let w = mix_one_col(v);
        for i in 0..4 { out[i][col] = w[i]; }
    }
    out
}

/// AddRoundKey : addition élément par élément modulo 17.
pub fn add_round_key_f17<const N: usize>(state: &[[F17; N]; N], key: &[[F17; N]; N]) -> [[F17; N]; N] {
    core::array::from_fn(|i| core::array::from_fn(|j| state[i][j].add(key[i][j])))
}

/// Un tour AES-like (N=4) : SubBytes → ShiftRows → MixColumns → AddRoundKey.
pub fn aes_like_round_f17(
    state: &[[F17; 4]; 4],
    round_key: &[[F17; 4]; 4],
) -> [[F17; 4]; 4] {
    let s1 = subbytes_f17(state);
    let s2 = shiftrows_f17(&s1);
    let s3 = mixcolumns_f17(&s2);
    add_round_key_f17(&s3, round_key)
}

/// Démo : applique un tour AES-like sur une matrice 4×4, puis imprime les coefficients.
#[allow(dead_code)]
pub fn demo_aes_like() {
    // État initial : valeurs 0..15 en row-major.
    let x0: [[F17; 4]; 4] = core::array::from_fn(|i| core::array::from_fn(|j| {
        let t = (i * 4 + j) as i32; F17::new(t)
    }));
    // Clé simple : lignes constantes 1,2,3,4.
    let key: [[F17; 4]; 4] = [
        [F17::new(1), F17::new(1), F17::new(1), F17::new(1)],
        [F17::new(2), F17::new(2), F17::new(2), F17::new(2)],
        [F17::new(3), F17::new(3), F17::new(3), F17::new(3)],
        [F17::new(4), F17::new(4), F17::new(4), F17::new(4)],
    ];

    let y = aes_like_round_f17(&x0, &key);
    println!("État initial :");
    for r in 0..4 { println!("{:?}", [x0[r][0].value(), x0[r][1].value(), x0[r][2].value(), x0[r][3].value()]); }
    println!("\nAprès un tour AES-like :");
    for r in 0..4 { println!("{:?}", [y[r][0].value(), y[r][1].value(), y[r][2].value(), y[r][3].value()]); }
}

// ==================
//  Tests unitaires
// ==================
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic_mod_17() {
        let x = F17::new(16);
        let y = F17::new(3);
        assert_eq!(x.add(y).value(), 2);  // 16+3=19 ≡ 2 (mod 17)
        assert_eq!(x.sub(y).value(), 13); // 16-3=13 (mod 17)
        assert_eq!(x.mul(y).value(), 14); // 16*3=48 ≡ 14 (mod 17)
    }

    #[test]
    fn test_phi_on_4x4_concrete() {
        const N: usize = 4;
        let s: FSMState<N> = FSMState::from_fn(|i, j| {
            let t = (i * N + j) + 1; // 1..16
            F17::new(t as i32)
        });
        let [a, b, c, d] = phi(&s);
        assert_eq!(a.value(), 4);
        assert_eq!(b.value(), 6);
        assert_eq!(c.value(), 12);
        assert_eq!(d.value(), 14);
    }

    #[test]
    fn test_phi_generic_on_4x4() {
        const N: usize = 4;
        let m: Mat17<F17, N> = Mat17::from_fn(|i, j| {
            let t = (i * N + j) + 1;
            F17::new(t as i32)
        });
        let [a, b, c, d] = phi_generic(&m);
        assert_eq!(a.value(), 4);
        assert_eq!(b.value(), 6);
        assert_eq!(c.value(), 12);
        assert_eq!(d.value(), 14);
    }

    #[test]
    fn test_aes_like_properties() {
        // On vérifie des propriétés simples (sans figer une matrice attendue).
        let x0: [[F17; 4]; 4] = core::array::from_fn(|i| core::array::from_fn(|j| {
            F17::new((i * 4 + j) as i32)
        }));
        let key: [[F17; 4]; 4] = [
            [F17::new(1), F17::new(1), F17::new(1), F17::new(1)],
            [F17::new(2), F17::new(2), F17::new(2), F17::new(2)],
            [F17::new(3), F17::new(3), F17::new(3), F17::new(3)],
            [F17::new(4), F17::new(4), F17::new(4), F17::new(4)],
        ];
        let y1 = aes_like_round_f17(&x0, &key);
        let y2 = aes_like_round_f17(&x0, &key);
        // Déterminisme
        for i in 0..4 { for j in 0..4 { assert_eq!(y1[i][j].value(), y2[i][j].value()); } }
        // Valeurs bien dans F17
        for i in 0..4 { for j in 0..4 { assert!(y1[i][j].value() <= 16); } }
    }
}
