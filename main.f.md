# main.f

## Overview

`main.f` is the main program of the `std2` package. It orchestrates the overall workflow, including parsing inputs, setting up calculations, calling computational modules, and managing resources. It serves as the primary entry point and control center for various quantum chemistry calculations like sTDA, sTD-DFT, and response function calculations.

## Key Components

*   **`program acis_prog`:** The main executable program. It initializes parameters, processes command-line arguments, reads molecular data, calls appropriate computational routines (like `stda`, `sutda`, `sfstda`, `linear_response` - though these calls are indirect via other modules), and manages overall execution flow.
*   **`subroutine byteout(s, mem)`:** A utility subroutine to print memory size in a human-readable format (KB, MB, GB).
*   **`subroutine mulpopcheck(nbf, nmo, c, occ, ireturn)`:** Performs a Mulliken population analysis to verify the consistency of the input molecular orbital data with the overlap matrix. This helps in detecting incorrect input formats.
*   **`subroutine inputcheck_printout(ires, typ, imethod, nbf, nmo)`:** Prints messages to the user regarding the detected or assumed input file format (e.g., ORCA, TURBOMOLE, MOLPRO).

## Important Variables/Constants

*   **`fname`:** Character string, stores the name of the primary input file (e.g., `molden.input`).
*   **`ax`:** Real, stores the amount of Fock-exchange. Critical for defining the DFT functional.
*   **`thre`:** Real, energy threshold for primary CSF space in sTDA/sTD-DFT calculations.
*   **`thrp`:** Real, PT (perturbation theory) threshold.
*   **`mform`:** Integer, specifies the format/style of the Molden input file.
*   **`rpachk`:** Logical, if `.true.`, sTD-DFT (RPA-like) calculation is performed instead of sTDA.
*   **`triplet`:** Logical, if `.true.`, triplet excited states are calculated.
*   **`xtbinp`:** Logical, if `.true.`, input is read from xTB output.
*   **`molden`:** Logical, if `.true.`, input is read from a Molden file.
*   **`resp`:** Logical, if `.true.`, response function calculations are enabled.
*   **`TPA`:** Logical, if `.true.`, two-photon absorption calculations are enabled.
*   **`ESA`:** Logical, if `.true.`, excited-state absorption calculations are enabled.
*   **`spinflip`:** Logical, if `.true.`, spin-flip calculations are performed.
*   **`XsTD`:** Logical, if `.true.`, XsTD methods are used.
*   **`cint`:** Logical, if `.true.`, the `libcint` library is used for integral calculations.
*   **`RSH_flag`:** Logical, if `.true.`, a range-separated hybrid functional is used (requires `XsTD`).
*   **`imethod`:** Integer, indicates the calculation method type (1 for RKS, 2 for UKS).
*   **`ncent`, `nmo`, `nbf`, `nprims`:** Integers, store the number of atoms (centers), molecular orbitals, basis functions, and primitive Gaussian functions, respectively.

## Usage Examples

```fortran
! Example: Running an sTDA calculation from a Molden file
! (Command line, not directly in code)
! std2 -f molden.input -ax 0.20 -e 7.0

! Inside the code, the program flow is:
! 1. Set default parameters
call date_and_time(VALUES=datetimevals) ! Print timestamp
! ... print banner ...
ptlim=1.7976931348623157d308
thre=7.0
! ... other defaults ...

! 2. Parse command line arguments
do i=1,command_argument_count()
  call getarg(i,dummy)
  if(index(dummy,'-f').ne.0)then
     call getarg(i+1,fname)
     molden=.true.
  endif
  ! ... parse other arguments ...
enddo

! 3. Read input data
if(molden)then
  call readmold0(ncent,nmo,nbf,nprims,fname,inpchk)
  ! ... allocate arrays ...
  call readmold(mform,imethod,ncent,nmo,nbf,nprims,cc,ccspin,icdim,fname)
else if(xtbinp) then
  ! ... read xTB input ...
else
  ! ... read other formats ...
endif

! 4. Precalculate integrals
if(cint)then
  call overlap(ncent,nprims,nbf,overlap_AO)
  call dipole_moment(ncent,nprims,nbf,mu)
  ! ... etc. ...
else
  call intslvm(ncent,nmo,nbf,nprims)
endif

! 5. Perform calculation (simplified example for RKS sTDA)
if (imethod.eq.1) then ! RKS
  if(.not.rw .and. .not.rw_dual .and. .not.spinflip)then
     call stda(ncent,nmo,nao,xyz,cc,eps,occ,iaoat,thre,
 .        thrp,ax,alpha,beta,ptlim,nvec,nprims)
  endif
  ! ... other calculation types ...
else ! UKS
  ! ... call sutda ...
endif
```

## Dependencies and Interactions

*   Modules: `stdacommon`, `kshiftcommon`, `commonlogicals`, `commonresp`. These modules provide shared variables, parameters, and potentially subroutines.
*   Input Files: Reads data from `.STDA` (for parameters), `molden.input` (or user-specified via `-f` or `-fo`), or xTB output.
*   Output Files: Generates various output files for NTOs (`NTOao`, `NTOat`, `NTOvar`, `NTOspin`, `fnorm`), integral files (`sint`, `xlint`, etc. which are temporary). Standard output prints calculation progress, results, and timings.
*   Libraries: Can link against `libcint` for efficient integral calculation. MKL or other BLAS/LAPACK libraries are implicitly used for linear algebra.
*   Core Computational Subroutines: Calls subroutines like `readmold0`, `readmold`, `readxtb0`, `readxtb`, `readbas0a`, `readbasa`, `readbasb` for input processing. Calls `intslvm`, `intslvm2` (or `libcint` routines like `overlap`, `dipole_moment`, `velo_moment`) for integral calculations. The main computational tasks are delegated to routines like `stda`, `sutda`, `sfstda`, and others likely found in `stda.f`, `sutda.f`, `sfstda.f`, `linear_response.f` etc. (though these are not directly called by name in `main.f` but rather through conditional logic that leads to their execution).
*   The program uses command-line arguments extensively to control its behavior.
```
