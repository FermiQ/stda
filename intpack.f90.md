# intpack.f90

## Overview

The `intpack` module provides subroutines for the calculation of one-electron integrals over Cartesian Gaussian basis functions. It implements the Obara-Saika scheme or similar methods involving Gaussian product theorem and recurrence relations to compute integrals like overlap, dipole moment, quadrupole moment, kinetic energy, and angular momentum.

## Key Components

*   **`SUBROUTINE propa0(aname, c, va, nt, ij, kl, efact)` and `SUBROUTINE propa1(aname, c, va, nt, ij, kl, efact)`:** These work in tandem. `propa0` sets up common terms for the product of two Gaussian primitives (from `stdacommon%eta` for basis functions `ij` and `kl`), including the product Gaussian's center, exponent, and prefactors using `divpt`, `rhftce`, and `prod`. `propa1` then uses these precomputed values to call a specific integral function `aname` (e.g., `opab1` for dipole) for various angular momentum combinations.
*   **`SUBROUTINE propa(aname, c, va, nt, ij, kl)`:** A general driver routine for calculating one-electron integrals. It takes an external function `aname` (which defines the type of integral, e.g., overlap, dipole) and indices `ij`, `kl` for the two basis functions. It uses the Gaussian product theorem (`divpt`, `rhftce`, `prod`) and then calls `aname` for the actual integral evaluation over combinations of angular momenta up to f-functions.
*   **`PURE SUBROUTINE divpt(a, alpha, b, beta, ca, cb, e, fgama, ffact)`:** Calculates the center (`e`), exponent (`fgama`), and prefactor (`ffact`) of a single Gaussian resulting from the product of two Gaussians.
*   **`PURE SUBROUTINE rhftce(cfs, a, e, iff)`:** Computes Cartesian prefactors (`cfs`) for a Gaussian centered at `a` when its center is shifted to `e`. `iff` denotes the angular momentum.
*   **`PURE SUBROUTINE prod(c, d, s, iff1, iff2)`:** Calculates the product of two sets of Cartesian prefactors (`c`, `d`) to yield the combined prefactors `s`.
*   **`ELEMENTAL FUNCTION olap(l, m, n, gama)`:** Computes the overlap integral for s-type Gaussian functions (or the s-component from a product of higher angular momentum Gaussians) with exponent `gama` and effective angular momentum components `l,m,n`.
*   **`PURE SUBROUTINE opab1(l, m, n, ga, v, d)`:** Calculates the x, y, z components of the electronic dipole moment integral.
*   **`PURE SUBROUTINE opab4(l, m, n, ga, v, d)`:** Calculates the components of the electronic second moment (quadrupole) integral.
*   **`PURE SUBROUTINE opad1(l, m, n, ga, v, d)`:** Calculates the electronic overlap integral.
*   **`SUBROUTINE opap4(l, m, n, gc, v, rc)`:** Calculates kinetic energy integrals and potentially relativistic corrections.
*   **`SUBROUTINE opam(l, m, n, gama, v, d)`:** Calculates angular momentum integrals, using `bip` as a helper.
*   **`SUBROUTINE bip(la, lb, a, b)`:** A helper subroutine, likely for recurrence relations in angular momentum integral calculations.
*   Other specific integral routines: `opac3` (charge density).
*   Helper routines: `lmnpre`, `olap2`.

## Important Variables/Constants

*   **`eta`:** Array from `stdacommon` module, expected to hold basis function parameters (center coordinates, exponents, angular momentum type index `iff`).
*   **`aname`:** External function passed to `propa` or `propa1`, determining the specific one-electron integral to be computed (e.g., `opad1` for overlap, `opab1` for dipole).
*   **`c`:** Real array, reference point (e.g., origin) for moment integrals.
*   **`va`:** Real array, output array to store the computed integral values.
*   **`nt`:** Integer, number of integral components to return (e.g., 1 for overlap, 3 for dipole).
*   **`ij`, `kl`:** Integers, indices pointing to the basis function parameters in the `eta` array.
*   Common blocks:
    *   **`/abfunc/`:** Stores intermediate values from Gaussian product rule (centers `ra`, `rb`, `e`; exponents `ga`, `gb`, `gama`; angular momentum indices `ia`, `ib`; Cartesian prefactors `aa`, `bb`, `dd`). Used by `propa`, `propa0`, `propa1` and the `op*` routines.
    *   **`/prptyp/ mprp`:** Integer, likely specifies the property type being calculated, controlling behavior in `propa` (e.g., `mprp=16` for magnetic moment).
    *   **`/gf/`:** Used by `opam` and `bip`.

## Usage Examples

```fortran
! This module is typically used by intslvm.f to calculate
! specific one-electron integrals needed by the main program.

! Example: Calculating an overlap integral (conceptual, actual call is from intslvm)
! Assuming 'opad1' is the function for overlap integrals
! and 'eta' array in stdacommon is populated.
! integer, parameter :: nt_overlap = 1
! real(wp) :: center(3), value_array(nt_overlap)
! integer :: basis_idx1, basis_idx2
!
! center = 0.0_wp ! Reference point (not always used by all integral types)
! basis_idx1 = 1  ! Index to first basis function in eta
! basis_idx2 = 2  ! Index to second basis function in eta
!
! call propa(opad1, center, value_array, nt_overlap, basis_idx1, basis_idx2)
! ! value_array(1) now holds the overlap integral S_12

! Inside propa (simplified):
! call divpt(a_coords, a_exp, b_coords, b_exp, ca, cb, prod_center, prod_exp, prefactor)
! call rhftce(aa_coeffs, a_coords, prod_center, a_angmom_idx)
! call rhftce(bb_coeffs, b_coords, prod_center, b_angmom_idx)
! call prod(aa_coeffs, bb_coeffs, dd_coeffs, a_angmom_idx, b_angmom_idx)
!
! ! Loop for different components if aname returns multiple values
! call aname(l_eff, m_eff, n_eff, prod_exp, v_temp, d_vec_to_ref)
! va(j) = prefactor * dd_coeffs(idx) * v_temp(j)
```

## Dependencies and Interactions

*   Modules: `iso_fortran_env` (for `wp => real64`), `stdacommon` (for basis set parameters in `eta` array and common variables like `ia`, `ib`, `etaij4`, `etakl4`, `iff1`, `iff2` used within `propa0`/`propa1`/`propa`).
*   External Functions: The `aname` argument to `propa` and `propa1` is a placeholder for actual integral calculation routines like `opad1` (overlap), `opab1` (dipole), `opab4` (quadrupole), `opap4` (kinetic energy), `opam` (angular momentum), `opac3` (charge density). These `op*` routines are defined within the `intpack` module itself.
*   Common Blocks: Uses `/abfunc/`, `/prptyp/`, and `/gf/` for internal communication between its subroutines.
*   Called by: Likely by routines in `intslvm.f` or `libcint.f` (if this is a fallback) to compute specific matrix elements of one-electron operators.
```
