# intslvm.f

## Overview

`intslvm.f` is responsible for computing and storing one-electron Atomic Orbital (AO) integrals. It iterates over pairs of primitive Gaussian basis functions, calculates various types of integrals (overlap, dipole length, magnetic dipole/angular momentum, dipole velocity) using routines from the `intpack` module, and writes these integrals to separate unformatted binary files for later use in the main computational routines.

## Key Components

*   **`SUBROUTINE intslvm(ncent, nmo, nbf, nprims)`:** The primary routine that calculates a comprehensive set of one-electron AO integrals. These include overlap (S), dipole length (R), magnetic dipole (L), and dipole velocity (P) integrals. It writes each type of integral component (e.g., Sx, Sy, Sz for dipole) to a distinct binary file (e.g., `sint`, `xlint`, `xmint`, `xvint`).
*   **`SUBROUTINE intslvm2(ncent, nmo, nbf, nprims)`:** A specialized version that calculates only the magnetic dipole (angular momentum) integrals (Lx, Ly, Lz) and writes them to `xmint`, `ymint`, and `zmint`. This routine is likely used when `libcint` handles other integrals, but the native calculation is preferred for magnetic moments.

## Important Variables/Constants

*   **`ncent, nmo, nbf, nprims`:** Integers representing the number of atoms, molecular orbitals, contracted basis functions, and primitive Gaussian functions, respectively. These define the scope of the integral calculations.
*   **`stdacommon`:** Module containing common data like atomic coordinates (`co`), primitive GTO parameters (`ipat`, `ipao`, `cxip`, `exip`), and basis function information (`eta` used by `intpack`).
*   **`intpack`:** Module providing the core routines (`propa0`, `propa1`, and operator-specific functions like `opad1`, `opab1`, `opam`) for evaluating integrals over Gaussian primitives.
*   **`thr`:** Real, a threshold (1.d-7) based on the overlap prefactor (`s` from `propa0`) to neglect integral contributions from distant or very diffuse primitive pairs.
*   **`mprp`:** Integer variable in common block `/prptyp/`, used to signal to `intpack` routines (specifically `propa1` via `aname`) which type of integral operator to use (e.g., `mprp=16` for angular momentum via `opam`).
*   **`r0`...`r9`:** Allocatable real arrays used to accumulate the different components of the calculated integrals (e.g., `r0` for S, `r1-r3` for R_x,R_y,R_z, `r4-r6` for L_x,L_y,L_z, `r7-r9` for P_x,P_y,P_z).
*   File units and names: e.g., unit 40 for `sint` (overlap), 31-33 for `xlint,ylint,zlint` (dipole length), 34-36 for `xmint,ymint,zmint` (magnetic dipole), 37-39 for `xvint,yvint,zvint` (dipole velocity).

## Usage Examples

```fortran
! intslvm is typically called from main.f if libcint is not used,
! or intslvm2 might be called for magnetic integrals specifically.
! Example call (conceptual, from main.f):
! IF (cint) THEN
!   ! ... use libcint for most integrals ...
!   CALL intslvm2(ncent,nmo,nbf,nprims) ! For magnetic moment integrals
! ELSE
!   CALL intslvm(ncent,nmo,nbf,nprims)  ! Calculate all integrals natively
! ENDIF

! Inside intslvm:
! ... setup, open files ...
! Loop over primitive pairs (i, j)
DO i=1,nprims
  iai=ipao(i)
  c1=cxip(i)
  DO j=1,i-1
    iaj=ipao(j)
    ! ... calculate index ij for packed storage ...
    cf = c1*cxip(j)*2.0d0 ! Combined contraction coefficient factor

    ! Calculate prefactor 's' for neglect test
    CALL propa0(opad1, point, v_dummy, 1, i, j, s)

    IF (s.GT.thr) THEN
      ! Calculate Overlap S
      mprp = 0 ! Signal overlap type (implicitly for opad1)
      CALL propa1(opad1, point, v_overlap, 1, i, j, s_eff_overlap)
      r0(ij) = r0(ij) + v_overlap(1)*cf

      ! Calculate Dipole Length R
      CALL propa1(opab1, point, v_dipole, 3, i, j, s_eff_dipole)
      r1(ij) = r1(ij) + v_dipole(1)*cf ! Rx
      r2(ij) = r2(ij) + v_dipole(2)*cf ! Ry
      r3(ij) = r3(ij) + v_dipole(3)*cf ! Rz

      ! Calculate Magnetic Dipole L
      mprp = 16 ! Signal angular momentum type
      CALL propa1(opam, point, v_angmom, 3, i, j, s_eff_angmom)
      r4(ij) = r4(ij) + v_angmom(1)*cf ! Lx
      r5(ij) = r5(ij) + v_angmom(2)*cf ! Ly
      r6(ij) = r6(ij) + v_angmom(3)*cf ! Lz

      ! Calculate Dipole Velocity P (velo call seems external or missing in provided snippet)
      ! CALL velo(i,j,v_velo)
      ! r7(ij)=r7(ij)-v_velo(1)*cf ! Px
      ! ...
    ENDIF
  ENDDO
  ! ... handle diagonal elements (i=j) ...
ENDDO
! ... scale off-diagonal elements by 0.5 ...
! ... write r0 to sint, r1 to xlint, etc. ...
```

## Dependencies and Interactions

*   Modules: `stdacommon` (provides basis set information, atomic coordinates), `intpack` (provides `propa0`, `propa1`, and the actual integral evaluation routines like `opad1`, `opab1`, `opam`).
*   Input Data: Relies on data populated in `stdacommon` (e.g., `co`, `ipat`, `ipao`, `cxip`, `exip`).
*   Output Files: Creates several unformatted binary files (`sint`, `xlint`, `ylint`, `zlint`, `xmint`, `ymint`, `zmint`, `xvint`, `yvint`, `zvint`) containing the calculated AO integrals. These files are subsequently read by other parts of the program (e.g., in `stda.f`).
*   Common Blocks: `/prptyp/` (for `mprp`), `/cema/` (for `cen`, `xmolw`), `/amass/` (for `ams`).
*   The `velo` subroutine called in `intslvm` is not defined within this file and is assumed to be an external routine or defined elsewhere.
```
