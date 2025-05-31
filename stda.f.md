# stda.f

## Overview

`stda.f` contains the core subroutine `stda` which implements the simplified Time-Dependent Density Functional Theory (sTDA and sTD-DFT) methods for calculating electronic excited states and response properties. It handles CSF generation, Hamiltonian construction, diagonalization, and property calculation for restricted Kohn-Sham (RKS) references.

## Key Components

*   **`SUBROUTINE stda(...)`:** The main driver for sTDA and sTD-DFT calculations. It orchestrates reading/transforming integrals, CSF selection, building and diagonalizing the Hamiltonian (either TDA or RPA-like), and computing various spectroscopic properties.
*   **`subroutine makel(nao, s, x)`:** Computes the Löwdin orthogonalization matrix `X = S^(-1/2)` from the overlap matrix `S`.
*   **`subroutine shrink(n, x, y)` / `subroutine shrink4(n, x, y)`:** Converts a symmetric square matrix `x` into a packed linear array `y` (double and single precision).
*   **`subroutine lo12pop(i, j, ncent, nao, iaoat, ca, q)`:** Calculates atomic Löwdin populations for a given transition density between MOs `i` and `j`.
*   **`subroutine kshift_to_ediag(de, ek)`:** Applies a K(ia,ia)-dependent energy shift to diagonal Hamiltonian elements, typically used with xTB methods.
*   **`subroutine ptselect(...)`:** Selects additional Configuration State Functions (CSFs) beyond an initial energy threshold using a perturbation theory (PT2) approach.
*   **`subroutine rrpamat(...)`:** Constructs the `A+B` and `A-B` matrices required for sTD-DFT (RPA) calculations for restricted Kohn-Sham (RKS) references.
*   **`subroutine rtdamat(...)`:** Constructs the TDA Hamiltonian matrix for RKS references.
*   **`subroutine setrep(rep)`:** Initializes an array `rep` with atomic hardness parameters (eta values) used in calculating repulsion integrals.
*   **`function lin(i1, i2)` / `function lin8(i1, i2)`:** Calculates a linear index for a pair of orbital indices (i1, i2) for storing elements of a symmetric matrix in packed form.
*   **`subroutine cofc(nat, xyz, coc)`:** Calculates the center of nuclear charge.
*   **`subroutine print_tdadat(...)`:** Writes calculated spectral data (excitation energies, oscillator strengths, rotatory strengths) to an output file (e.g., `tda.dat`) for plotting.

## Important Variables/Constants

*   **`ncent`, `nmo`, `nao`, `nprims`:** Integers for number of centers, MOs, AOs, and primitive GTOs.
*   **`thr`:** Real, energy threshold (eV, converted to a.u.) for selecting primary CSFs.
*   **`thrp`:** Real, threshold for perturbation selection of CSFs.
*   **`ax`:** Real, Fock exchange mixing parameter.
*   **`alphak`, `betaj`:** Real, parameters controlling the screening of K and J integrals, respectively.
*   **`fthr`:** Real, maximum energy range for CSF selection in PT2.
*   **`triplet`:** Logical, `.TRUE.` if triplet states are to be calculated.
*   **`rpachk`:** Logical, `.TRUE.` if sTD-DFT (RPA) is performed instead of sTDA.
*   **`eigvec`:** Logical, `.TRUE.` to print eigenvectors.
*   **`dokshift`:** Logical, `.TRUE.` to apply K(ia,ia)-dependent shifts.
*   **`velcorr`:** Logical, `.TRUE.` to apply velocity correction for rotatory strengths.
*   **`Xcore`:** Logical, `.TRUE.` if core excitations are considered.
*   **`XsTD`:** Logical, `.TRUE.` if using exact-exchange sTD-DFT/sTDA.
*   **`RSH_flag`:** Logical, if `.true.`, a range-separated hybrid functional is used.
*   **`iconf`:** Integer array, stores the occupied and virtual orbital indices for each selected CSF.
*   **`ed`:** Real array, stores the diagonal elements (initial energies) of the CSFs.
*   **`xl, yl, zl`:** Real arrays, store x, y, z components of dipole length integrals in MO basis.
*   **`xv, yv, zv`:** Real arrays, store x, y, z components of dipole velocity integrals in MO basis.
*   **`xm, ym, zm`:** Real arrays, store x, y, z components of magnetic dipole integrals in MO basis.
*   **`gamj, gamk`:** Real arrays, store screened Coulomb (J) and exchange-like (K) repulsion terms between atomic centers.
*   **`hci`:** Real array, the TDA Hamiltonian matrix.
*   **`apb, ambsqr`:** Real arrays, `A+B` and `(A-B)^1/2` (or `A-B` directly if response) matrices for sTD-DFT.
*   **`eci, uci`:** Real arrays, eigenvalues and eigenvectors of the TDA/RPA Hamiltonian.

## Usage Examples

```fortran
! This subroutine is called from main.f, e.g.:
! if (imethod.eq.1) then ! RKS
!   if(.not.rw .and. .not.rw_dual .and. .not.spinflip)then
!      call stda(ncent,nmo,nao,xyz,cc,eps,occ,iaoat,thre,
!  .        thrp,ax,alpha,beta,ptlim,nvec,nprims)
!   endif
! endif

! Inside stda:
! Determine active MOs
do i=1,nmo
   if(occ(i).gt.1.990d0.and.eps(i).gt.othr) moci=moci+1 ! Count occupied
   if(occ(i).lt.0.010d0.and.eps(i).lt.vthr) moci=moci+1 ! Count virtual
enddo
! ...
! Transform integrals (e.g., dipole length)
call onetri(1,help,dum,scr,ca,nao,moci) ! Transform AO integrals to MO basis
call shrink(moci,dum,xl)                ! Pack to linear array

! Calculate repulsion terms
call setrep(eta) ! Load atomic hardness
! ... calculate gamj, gamk ...

! Generate CSFs and diagonal elements
do io=1,no ! occ loop
   do iv=no+1,moci ! virt loop
      de=epsi(iv)-epsi(io) ! Orbital energy difference
      ! ... calculate ej (Coulomb) and ek (Exchange) ...
      de=de-ej+ak*ek       ! Zeroth-order energy
      if(de.le.thr)then
         k=k+1
         iconf(k,1)=io
         iconf(k,2)=iv
         ed(k)=de
      endif
   enddo
enddo
! ... PT selection via ptselect ...

! Build and diagonalize Hamiltonian
if(rpachk) then
  call rrpamat(...) ! Build A+B, A-B
  call srpapack(...) ! Diagonalize RPA equations
else
  call rtdamat(...)  ! Build TDA matrix
  call ssyevr(...) ! Diagonalize TDA matrix
endif

! Calculate properties and print
! ... loop over roots ...
! fl = de*p23*xp ! Oscillator strength
! rl = -235.7220d0*xp ! Rotatory strength
! call print_tdadat(...)
```

## Dependencies and Interactions

*   Modules: `commonlogicals`, `commonresp`, `omp_lib`, `kshiftcommon`.
*   Input Data: Expects molecular orbitals (`c`, `eps`, `occ`), atomic coordinates (`xyz`), AO-to-atom map (`iaoat`), various thresholds (`thr`, `thrp`, `fthr`), and method parameters (`ax`, `alphak`, `betaj`). One-electron AO integrals are read from temporary files (`xlint`, `ylint`, `zlint`, `xmint`, `ymint`, `zmint`, `xvint`, `yvint`, `zvint`, `sint`) created by `main.f` or `intslvm`/`libcint`.
*   Output: Writes spectral data to `tda.dat` (and potentially `tdax.dat`, etc., for anisotropic properties). Prints detailed results to standard output. Can generate `TmPvEcInFo` for eigenvector printing. Can call `print_nto`, `print_nto_rpa`, `print_nto_resp` for NTO analysis.
*   External Libraries: LAPACK (e.g., `dsyev`, `dsyevr`, `dgemm`, `ssymm`, `sgemm`, `sdot`). OpenMP for parallelization.
*   Interacts with other subroutines for specific tasks like `onetri` (from `io.f` presumably, for integral transformation), `srpapack` (for RPA solution), various linear response routines (`lresp1`, `lresp_2PA_SP`, `optrot`, `lresp_ESA`, `tpa_sos`), and potentially `Xstda_mat`, `xstd_rpamat` if `XsTD` is true (likely from `xstd.f90`).
```
