# linear_response.f

## Overview

`linear_response.f` is a crucial component of the `std2` package, responsible for calculating various linear and nonlinear optical properties. It implements the necessary machinery to compute frequency-dependent polarizabilities, first hyperpolarizabilities (for SHG and HRS), two-photon absorption cross-sections, excited-state absorption properties, and optical rotation, primarily using sTD-DFT response theory.

## Key Components

*   **`SUBROUTINE lresp(...)`:** Calculates frequency-dependent dynamic linear polarizabilities (alpha).
*   **`SUBROUTINE lresp1(...)`:** Calculates frequency-dependent first hyperpolarizabilities (beta) for Second Harmonic Generation (SHG), including Hyper-Rayleigh Scattering (HRS) quantities. Utilizes helper routines like `beta_resp_fast`, `List_A`, `List_B`.
*   **`SUBROUTINE PrintBeta(beta, unit)`:** Prints the 3x3x3 beta tensor in a formatted way.
*   **`SUBROUTINE HRS(beta, beta2_ZZZ, beta2_XZZ, beta_HRS, DR)`:** Computes Hyper-Rayleigh Scattering intensity (`beta_HRS`) and depolarization ratio (`DR`) from the beta tensor.
*   **`SUBROUTINE lresp_2PA(...)` / `SUBROUTINE lresp_2PA_SP(...)`:** Calculate the two-photon absorption (2PA) tensor and related isotropic averages. `_SP` might be a single-precision or specialized version.
*   **`SUBROUTINE lresp_ESA(...)` / `SUBROUTINE lresp_ESA_tda(...)`:** Calculate state-to-state transition dipole moments and oscillator strengths for Excited State Absorption (ESA) from a specified initial excited state, for TD-DFT and TDA respectively.
*   **`SUBROUTINE ulresp_ESA(...)` / `SUBROUTINE ulresp_ESA_tda(...)`:** Versions of ESA calculation for unrestricted (UKS) references.
*   **`SUBROUTINE sf_lresp_ESA(...)`:** ESA calculation for spin-flip TD-DFT.
*   **`SUBROUTINE optrot(...)` / `SUBROUTINE optrot_velo(...)`:** Calculate frequency-dependent optical rotation. `optrot_velo` uses the velocity-gauge formalism.
*   **`SUBROUTINE pol_sos(...)`:** Calculates polarizabilities using a sum-over-states (SOS) formulation.
*   **`SUBROUTINE hyperpol_sos(...)`:** Calculates first hyperpolarizabilities using an SOS formulation.
*   **`SUBROUTINE tpa_sos(...)`:** Calculates 2PA cross-sections using an SOS formulation.
*   Helper routines for beta calculation: `beta_resp_fast`, `List_A`, `List_B`, `A_beta_fast`, `B_beta_fast`, `beta_resp`, `A_beta1`, `B_beta1`.
*   Helper routines for 2PA calculation: `TPA_resp_fast`, `TPA_resp_fast_SP`, `TPA_resp`, and their lower-level components like `A_2PA_1_fast_SP`, etc.

## Important Variables/Constants

*   **`nci`:** Integer, number of configuration state functions (CSFs) in the response calculation.
*   **`apb`, `amb`:** Real arrays (packed symmetric), representing the `A+B` and `A-B` matrices from sTD-DFT.
*   **`iconf`:** Integer array, stores the orbital indices for each CSF.
*   **`xl, yl, zl`:** Real arrays, dipole length integrals in MO basis (x,y,z components).
*   **`xm, ym, zm`:** Real arrays, magnetic dipole integrals in MO basis (x,y,z components).
*   **`xv, yv, zv`:** Real arrays, dipole velocity integrals in MO basis (x,y,z components, used in `optrot_velo`).
*   **`moci, no, nv`:** Integers, total number of active MOs, number of occupied MOs, number of virtual MOs.
*   **`eci, Xci, Yci`:** Real arrays, eigenvalues, X (excitation), and Y (de-excitation) amplitudes from the underlying sTD-DFT calculation (passed for 2PA and ESA).
*   **`freq`:** Real array, stores the frequencies (in a.u.) at which properties are calculated. Read from `wavelength` file.
*   **`num_freq`:** Integer, number of frequencies to calculate properties for (from `commonresp` module).
*   **`state2opt`:** Integer, specifies the initial excited state for ESA calculations (from `commonresp` module).
*   **`nto`:** Logical, if true, some routines might save NTO information (from `commonlogicals` module).

## Usage Examples

```fortran
! These subroutines are typically called from stda.f or sutda.f
! after the main sTD-DFT diagonalization, if response properties
! are requested via command-line flags.

! Example: Calling lresp1 for SHG beta calculation (inside stda.f)
! if(resp .and. .not. aresp .and. .not. TPA .and. .not. optrota) then
!   allocate( amb(nci*(nci+1)/2), stat=ierr )
!   open(unit=53,file='amb',form='unformatted',status='old')
!   read(53) amb ! Read A-B matrix previously stored
!   close(53,status='delete')
!   call lresp1(nci,apb,amb,iconf,maxconf,xl,yl,zl,moci,no,nv)
! endif

! Inside lresp1:
! Read frequencies
open(unit=101,file='wavelength',form='formatted',status='old')
! ... read freq ...
close(101)

! Invert (A-B) matrix
inv_amb=amb
call ssptrf(uplo,nci,inv_amb,ipiv,info) ! LU factorization
call ssptri(uplo,nci,inv_amb,ipiv,work,info) ! Matrix inverse

! Loop over frequencies
Do ii=1, num_freq+1
  omega=freq(ii)
  ! Solve response equation for X+Y: ( (A+B) - omega^2 * (A-B)^-1 ) * (X+Y) = mu
  inv_resp = apb - omega**2.0 * inv_amb
  XpY_int(:,1)=mu_x(:) ! Dipole vector
  ! ... solve for XpY_int ...
  call ssptrs(uplo,nci,3,inv_resp,ipiv,XpY_int,nci,info)

  ! Calculate X and Y vectors
  ! XmY = omega * inv_amb * XpY
  ! X = (XpY + XmY)/2
  ! Y = XpY - X

  ! Calculate beta tensor components using X, Y and mu
  ! e.g. call beta_resp_fast(1,1,1,X,Y,...)
  ! beta(ix,iy,iz) = Sum_terms (X * mu * Y)
enddo
```

## Dependencies and Interactions

*   Modules: `commonresp` (for `num_freq`, `state2opt`), `commonlogicals` (for `nto`), `omp_lib`.
*   Input Data: Requires `A+B` and `A-B` matrices, CSF definitions (`iconf`), orbital indices (`no`, `nv`, `moci`), and transformed one-electron property integrals (dipole length `xl,yl,zl`; magnetic dipole `xm,ym,zm`; velocity dipole `xv,yv,zv` for `optrot_velo`) from the calling sTDA/sTD-DFT routine. Eigenvectors (`Xci`, `Yci`) and eigenvalues (`eci`) from sTD-DFT are needed for 2PA and ESA. Reads frequencies from a file named `wavelength`.
*   Output: Prints calculated property tensors (alpha, beta, 2PA sigma) and derived quantities (HRS intensity, depolarization ratio, isotropic 2PA values) to standard output. May write detailed tensor components to `beta_HRS` and `beta_tensor`. For 2PA, writes to `2PA-abs`. For ESA, writes to `s2s.dat`. For NTO analysis with `optrot` or `lresp`, may write response vectors to temporary files.
*   External Libraries: LAPACK (e.g., `ssptrf`, `ssptri`, `ssptrs` for solving linear equations and matrix inversion). OpenMP for parallelization.
*   Called from: `stda.f` and `sutda.f` (for unrestricted calculations, though `sutda.f` is not directly analyzed here, its interaction would be analogous for `ulresp_ESA`).
```
