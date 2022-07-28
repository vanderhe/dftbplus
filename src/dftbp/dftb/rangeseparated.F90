!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Contains range separated related routines.
module dftbp_dftb_rangeseparated

  use dftbp_common_accuracy, only : dp, tolSameDist, MinHubDiff
  use dftbp_common_constants, only : pi
  use dftbp_common_environment, only : TEnvironment, globalTimers
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_nonscc, only : TNonSccDiff
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_dftb_sparse2dense, only : blockSymmetrizeHS, symmetrizeHS, hermitianSquareMatrix
  use dftbp_io_message, only : error
  use dftbp_math_blasroutines, only : gemm
  use dftbp_math_sorting, only : index_heap_sort
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_dftb_periodic, only : TNeighbourList, TSymNeighbourList, getCellTranslations,&
      & getLatticePoints, cart2frac
  use dftbp_dftb_nonscc, only : buildS
  use dftbp_dftb_densitymatrix, only : TDensityMatrix
  use dftbp_math_simplealgebra, only : determinant33
  use dftbp_common_parallel, only : getStartAndEndIndex

#:if WITH_OMP
  use omp_lib, only : OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
#:endif

#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip, mpifx_allreduce
#:endif

  implicit none
  private

  public :: TRangeSepSKTag, TRangeSepFunc, TRangeSepFunc_init
  public :: getGammaPrimeValue, rangeSepTypes, checkSupercellFoldingMatrix


  !> Returns the derivative of long-range gamma.
  interface getGammaPrimeValue
    module procedure getLrGammaPrimeValue_cluster
    module procedure getLrGammaPrimeValue_periodic
  end interface getGammaPrimeValue


  type :: TIntArray1D

    !> onedimensional, integer data storage
    integer, allocatable :: array(:)

  end type TIntArray1D


  type :: TRealArray1D

    !> onedimensional, real data storage
    real(dp), allocatable :: array(:)

  end type TRealArray1D


  type :: TRangeSepTypesEnum

    !> Neighbour based
    integer :: neighbour = 0

    !> Threshold based
    integer :: threshold = 1

    !> Matrix based
    integer :: matrixBased = 2

  end type TRangeSepTypesEnum


  !> Container for enumerated range separation types
  type(TRangeSepTypesEnum), parameter :: rangeSepTypes = TRangeSepTypesEnum()


  !> Slater-Koster file RangeSep tag structure
  type :: TRangeSepSKTag

    !> range-separation parameter
    real(dp) :: omega

    !> CAM alpha parameter
    real(dp) :: camAlpha

    !> CAM beta parameter
    real(dp) :: camBeta

  end type TRangeSepSKTag


  !> Range-Sep module structure
  type :: TRangeSepFunc

    !> Real-space coordinates of atoms (relative units), potentially including periodic images
    real(dp), allocatable :: coords(:,:)

    !> Real-space coordinates of atoms (absolute units), potentially including periodic images
    real(dp), allocatable :: rCoords(:,:)

    !> Evaluated long-range gamma of Atom1 and Atom2 (central cell only)
    real(dp), allocatable :: lrGammaEval0(:,:), lrdGammaEval0(:,:,:)

    !> Evaluated Hartree-Fock gamma of Atom1 and Atom2 (central cell only)
    real(dp), allocatable :: hfGammaEval0(:,:), hfdGammaEval0(:,:,:)

    !> Range-separation parameter
    real(dp) :: omega

    !> CAM alpha parameter
    real(dp) :: camAlpha

    !> CAM beta parameter
    real(dp) :: camBeta

    !> True, for global hybrids
    logical :: tHyb

    !> True, for long-range corrected functionals
    logical :: tLc

    !> True, for general CAM range-separation
    logical :: tCam

    !> Hubbard U values for atoms
    real(dp), allocatable :: hubbu(:)

    !> Full overlap as obtained by buildS() for a symmetric neighbour list
    real(dp), allocatable :: overSym(:)

    ! Hamiltonian Screening

    !> Previous hamiltonian in screening by tolerance
    real(dp), allocatable :: hPrev(:,:)

    !> Previous delta density matrix in screening by tolerance
    real(dp), allocatable :: dRhoPrev(:,:)

    !> Previous hamiltonian in screening by tolerance
    complex(dp), allocatable :: hPrevCplxHS(:,:,:)

    !> Previous delta density matrix in screening by tolerance
    real(dp), allocatable :: dRhoPrevCplxHS(:,:,:,:,:,:)

    !> Is screening initialised
    logical :: tScreeningInited

    !> Threshold for screening by value
    real(dp) :: pScreeningThreshold

    !> Total full-range Hartree-Fock energy
    real(dp) :: hfEnergy

    !> Total long-range energy
    real(dp) :: lrEnergy

    !> Is this spin restricted (F) or unrestricted (T)
    logical :: tSpin

    !> Is this DFTB/SSR formalism
    logical :: tREKS

    !> Algorithm for range separation screening
    integer :: rsAlg

    !> Species of atoms in central cell
    integer, allocatable :: species0(:)

    !> Cutoff for real-space g-summation
    real(dp) :: gSummationCutoff

    !> Cutoff for truncated Gamma
    real(dp) :: gammaCutoff

    !> Damping distance for Gamma truncation
    real(dp) :: gammaDamping

    !> Auxiliary gamma damping/screening parameter
    real(dp) :: auxiliaryScreening

    !> Value, 1st and 2nd derivative of gamma integral at damping distance
    real(dp), allocatable :: lrGammaAtDamping(:,:), lrdGammaAtDamping(:,:), lrddGammaAtDamping(:,:)

    !> Value, 1st and 2nd derivative of gamma integral at damping distance
    real(dp), allocatable :: hfGammaAtDamping(:,:), hfdGammaAtDamping(:,:), hfddGammaAtDamping(:,:)

    !> K-point compatible BvK real-space shifts in relative coordinates (units of latVecs)
    real(dp), allocatable :: bvKShifts(:,:)

    !> Supercell folding coefficients (diagonal elements)
    integer, allocatable :: coeffsDiag(:)

  contains

    procedure :: updateCoords_cluster, updateCoords_gamma

    procedure :: foldToBvK => TRangeSepFunc_foldToBvK
    procedure :: foldToBvKIndex => TRangeSepFunc_foldToBvKIndex

    procedure :: addCamHamiltonian_cluster
    procedure :: addCamHamiltonian_gamma
    procedure :: addCamHamiltonian_kpts

    procedure :: addLrHamiltonianMatrixCmplx

    procedure :: addCamEnergy

    procedure :: addCamGradients_cluster
    procedure :: addCamGradients_gamma

    procedure :: getCentralCellSpecies
    procedure :: getLrGammaCluster
    procedure :: getLrGammaDerivCluster
    procedure :: getHfGammaCluster
    procedure :: getHfGammaDerivCluster

    procedure :: evaluateLrEnergyDirect_cluster

  end type TRangeSepFunc


contains


  !> Intitializes the range-sep module.
  subroutine TRangeSepFunc_init(this, nAtom, species0, hubbu, screen, omega, camAlpha, camBeta,&
      & tSpin, tREKS, rsAlg, gammaCutoff, gSummationCutoff, auxiliaryScreening, coeffsDiag)

    !> Instance
    type(TRangeSepFunc), intent(out) :: this

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> List of all atomic species in the central cell
    integer, intent(in) :: species0(:)

    !> Atomic hubbards
    real(dp), intent(in) :: hubbu(:)

    !> Screening threshold value
    real(dp), intent(in) :: screen

    !> Range separation parameter
    real(dp), intent(in) :: omega

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> CAM beta parameter
    real(dp), intent(in) :: camBeta

    !> Is this spin restricted (F) or unrestricted (T)
    logical, intent(in) :: tSpin

    !> Is this DFTB/SSR formalism
    logical, intent(in) :: tREKS

    !> lr-hamiltonian construction algorithm
    integer, intent(in) :: rsAlg

    !> Cutoff for truncated Gamma
    real(dp), intent(in), optional :: gammaCutoff

    !> Cutoff for real-space g-summation
    real(dp), intent(in), optional :: gSummationCutoff

    !> Auxiliary gamma damping/screening parameter
    real(dp), intent(in), optional :: auxiliaryScreening

    !> Supercell folding coefficients (diagonal elements)
    integer, intent(in), optional :: coeffsDiag(:)

    !! Species indices
    integer :: iSp1, iSp2

    !! Number of unique species in system
    integer :: nUniqueSpecies

    this%tScreeningInited = .false.
    this%pScreeningThreshold = screen

    this%omega = omega
    this%rsAlg = rsAlg
    this%tSpin = tSpin
    this%tREKS = tREKS
    this%camAlpha = camAlpha
    this%camBeta = camBeta
    this%hubbu = hubbu
    this%species0 = species0

    this%lrEnergy = 0.0_dp
    this%hfEnergy = 0.0_dp

    if (present(gammaCutoff)) then
      this%gammaCutoff = gammaCutoff

      ! Is this a sensible choice for the beginning of the damping region?
      this%gammaDamping = 0.95_dp * this%gammaCutoff

      if (this%gammaDamping <= 0.0_dp) then
        call error("Beginning of damped region of electrostatics must be positive.")
      end if

      ! Tabulate truncated Gamma properties for all (symmetric) combinations of species
      call getNumberOfUniqueInt(this%species0, nUniqueSpecies)
      allocate(this%lrGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      allocate(this%lrdGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      allocate(this%lrddGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      do iSp2 = 1, nUniqueSpecies
        do iSp1 = 1, nUniqueSpecies
          this%lrGammaAtDamping(iSp1, iSp2) = getAnalyticalLrGammaValue(this, iSp1, iSp2,&
              & this%omega, this%gammaDamping)
          this%lrdGammaAtDamping(iSp1, iSp2) = getdAnalyticalLrGammaDeriv(this, iSp1, iSp2,&
              & this%gammaDamping)
          this%lrddGammaAtDamping(iSp1, iSp2) = getddNumericalLrGammaDeriv(this, iSp1, iSp2,&
              & this%gammaDamping, 1e-08_dp)
        end do
      end do

      call getNumberOfUniqueInt(this%species0, nUniqueSpecies)
      allocate(this%hfGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      allocate(this%hfdGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      allocate(this%hfddGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      do iSp2 = 1, nUniqueSpecies
        do iSp1 = 1, nUniqueSpecies
          this%hfGammaAtDamping(iSp1, iSp2) = getAnalyticalHfGammaValue(this, iSp1, iSp2,&
              & this%gammaDamping)
          this%hfdGammaAtDamping(iSp1, iSp2) = getdAnalyticalHfGammaDeriv(this, iSp1,&
              & iSp2, this%gammaDamping)
          this%hfddGammaAtDamping(iSp1, iSp2) = getddNumericalHfGammaDeriv(this, iSp1,&
              & iSp2, this%gammaDamping, 1e-08_dp)
        end do
      end do
    end if

    if (present(gSummationCutoff)) this%gSummationCutoff = gSummationCutoff
    if (present(auxiliaryScreening)) this%auxiliaryScreening = auxiliaryScreening

    this%tHyb = .false.
    this%tLc = .false.
    this%tCam = .false.

    ! This test is just for saving time while calculating Hamiltonian and gradient contributions.
    ! In theory the most general CAM case covers everything, but is not needed for pure Hyb/LC.
    if ((abs(this%camAlpha) > 1.0e-16_dp) .and. (abs(this%camBeta) < 1.0e-16_dp)) then
      ! apparently this is a pure global hybrid calculation
      this%tHyb = .true.
    elseif ((abs(this%camAlpha) < 1.0e-16_dp) .and.&
        & (abs(this%camBeta - 1.0_dp) < 1.0e-16_dp)) then
      ! apparently this is a pure LC calculation
      this%tLc = .true.
    else
      this%tCam = .true.
    end if

    if (this%tREKS .and. this%tHyb) then
      call error("Global hybrid functionals not currently implemented for REKS.")
    end if

    if (this%tREKS .and. this%tCam) then
      call error("General CAM functionals not currently implemented for REKS.")
    end if

    allocate(this%coords(3, nAtom))
    this%coords(:,:) = 0.0_dp
    allocate(this%rCoords(3, nAtom))
    this%rCoords(:,:) = 0.0_dp

    allocate(this%lrGammaEval0(nAtom, nAtom))
    this%lrGammaEval0(:,:) = 0.0_dp
    allocate(this%lrdGammaEval0(nAtom, nAtom, 3))
    this%lrdGammaEval0(:,:,:) = 0.0_dp

    allocate(this%hfGammaEval0(nAtom, nAtom))
    this%hfGammaEval0(:,:) = 0.0_dp
    allocate(this%hfdGammaEval0(nAtom, nAtom, 3))
    this%hfdGammaEval0(:,:,:) = 0.0_dp

    ! Check for current restrictions
    if (this%tSpin .and. this%rsAlg == rangeSepTypes%threshold) then
      call error("Spin-unrestricted calculation for thresholded range separation not yet&
          & implemented!")
    end if

    if (this%tREKS .and. this%rsAlg == rangeSepTypes%threshold) then
      call error("REKS calculation with thresholded range separation not yet implemented!")
    end if

    if (.not. any([rangeSepTypes%neighbour, rangeSepTypes%threshold,&
          & rangeSepTypes%matrixBased] == this%rsAlg)) then
      call error("Unknown algorithm for screening the exchange in range separation!")
    end if

    if (present(coeffsDiag)) then
      this%coeffsDiag = coeffsDiag
      call getBvKLatticeShifts(this%coeffsDiag, this%bvKShifts)
    end if


  contains

    !> Returns the number of unique integers of an array.
    pure subroutine getNumberOfUniqueInt(array, nUnique)

      !> Array to investigate
      integer, intent(in) :: array(:)

      !> Number of unique entries
      integer, intent(out) :: nUnique

      !! Auxiliary variables
      integer :: ii, jj, tmp(size(array))

      nUnique = 1
      tmp(1) = array(1)

      outer: do ii = 2, size(array)
        do jj = 1, nUnique
          if (tmp(jj) == array(ii)) then
            cycle outer
          end if
        end do
        nUnique = nUnique + 1
        tmp(nUnique) = array(ii)
      end do outer

    end subroutine getNumberOfUniqueInt


    !> Returns BvK real-space shifts, compatible with k-point mesh.
    subroutine getBvKLatticeShifts(coeffsDiag, bvKShifts)

      !> Supercell folding coefficients (diagonal elements)
      integer, intent(in) :: coeffsDiag(:)

      !> K-point compatible BvK real-space shifts in relative coordinates
      real(dp), intent(out), allocatable :: bvKShifts(:,:)

      !! Number of BvK real-space shifts
      integer :: nBvKShifts

      !! Auxiliary variables
      integer :: ii, jj, kk, ind

      nBvKShifts = coeffsDiag(1) * coeffsDiag(2) * coeffsDiag(3)
      allocate(bvKShifts(3, nBvKShifts))

      ind = 1
      do kk = 0, coeffsDiag(3) - 1
        do jj = 0, coeffsDiag(2) - 1
          do ii = 0, coeffsDiag(1) - 1
            bvKShifts(1, ind) = real(ii, dp)
            bvKShifts(2, ind) = real(jj, dp)
            bvKShifts(3, ind) = real(kk, dp)
            ind = ind + 1
          end do
        end do
      end do

    end subroutine getBvKLatticeShifts

  end subroutine TRangeSepFunc_init


  !> Check if obtained supercell folding matrix meets current requirements
  subroutine checkSupercellFoldingMatrix(supercellFoldingMatrix, supercellFoldingDiagOut)

    !> Coefficients of the lattice vectors in the linear combination for the super lattice vectors
    !! (should be integer values) and shift of the grid along the three small reciprocal lattice
    !! vectors (between 0.0 and 1.0)
    real(dp), intent(in), target :: supercellFoldingMatrix(:,:)

    !> Diagonal elements of supercell folding matrix, if present
    integer, intent(out), optional :: supercellFoldingDiagOut(:)

    !! Supercell folding coefficients and shifts
    real(dp), pointer :: coeffs(:,:), shifts(:)

    !! True, if the supercell folding does not correspond to a MP-like scheme
    logical :: tNotMonkhorstPack

    !! Auxiliary variables
    integer :: ii, jj

    if (present(supercellFoldingDiagOut)) then
      @:ASSERT(size(supercellFoldingDiagOut) == 3)
    end if

    coeffs => supercellFoldingMatrix(:, 1:3)
    shifts => supercellFoldingMatrix(:, 4)

    if (abs(determinant33(coeffs)) - 1.0_dp < -1e-06_dp) then
      call error('Determinant of the supercell matrix must be greater than 1.')
    end if

    if (any(abs(modulo(coeffs + 0.5_dp, 1.0_dp) - 0.5_dp) > 1e-6_dp)) then
      call error('The components of the supercell matrix must be integers.')
    end if

    ! Check if k-point mesh is a Monkhorst-Pack sampling with zero shift
    tNotMonkhorstPack = .false.
    lpOuter: do jj = 1, size(coeffs, dim=2)
      do ii = 1, size(coeffs, dim=1)
        if (ii == jj) cycle
        if (coeffs(ii, jj) > 1e-06_dp) then
          tNotMonkhorstPack = .true.
          exit lpOuter
        end if
      end do
    end do lpOuter
    if (tNotMonkhorstPack) then
      call error('Range-separated calculations with k-points require a Monkhorst-Pack-like&
          & sampling, i.e. a uniform extension of the lattice.')
    end if

    ! Check if shifts are zero
    if (any(abs(shifts) > 1e-06_dp)) then
      call error('Range-separated calculations with k-points require a Monkhorst-Pack-like&
          & sampling with zero shift.')
    end if

    ! All checks have passed, continue...

    ! Get diagonal elements as integers, if requested
    if (present(supercellFoldingDiagOut)) then
      do ii = 1, 3
        supercellFoldingDiagOut(ii) = nint(coeffs(ii, ii))
      end do
    end if

  end subroutine checkSupercellFoldingMatrix


  !> Folds relative real-space vector back to BvK region.
  pure function TRangeSepFunc_foldToBvK(this, vector) result(bvKShift)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> Vector (in relative coordinates) to fold back to BvK cell
    real(dp), intent(in) :: vector(:)

    !> Corresponding BvK vector
    integer :: bvKShift(3)

    bvKShift(:) = modulo(nint(vector), this%coeffsDiag)

  end function TRangeSepFunc_foldToBvK


  !> Folds relative real-space vector back to BvK region and returns indices of density matrix.
  pure function TRangeSepFunc_foldToBvKIndex(this, vector) result(bvKIndex)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> Vector (in relative coordinates) to fold back to BvK cell
    real(dp), intent(in) :: vector(:)

    !> Corresponding BvK indexing
    integer :: bvKIndex(3)

    ! additionally shift by 1, so that indices start at 1 and not at 0
    bvKIndex(:) = modulo(nint(vector) + this%coeffsDiag, this%coeffsDiag) + 1

  end function TRangeSepFunc_foldToBvKIndex


  !> Updates the rangeSep module on coordinate change.
  subroutine updateCoords_cluster(this, rCoords)

    !> class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> Atomic coordinates in absolute units
    real(dp), intent(in) :: rCoords(:,:)

    !! Indices of interacting atoms in central cell, as well as their global species index
    integer :: iAtom1, iAtom2, iSp1, iSp2

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Distance between two interacting atoms
    real(dp) :: dist

    this%rCoords(:,:) = rCoords
    nAtom0 = size(this%species0)

    do iAtom1 = 1, nAtom0
      do iAtom2 = 1, iAtom1
        iSp1 = this%species0(iAtom1)
        iSp2 = this%species0(iAtom2)
        dist = norm2(this%rCoords(:, iAtom1) - this%rCoords(:, iAtom2))
        this%lrGammaEval0(iAtom1, iAtom2) = getAnalyticalLrGammaValue(this, iSp1, iSp2, this%omega,&
            & dist)
        this%lrGammaEval0(iAtom2, iAtom1) = this%lrGammaEval0(iAtom1, iAtom2)
      end do
    end do

    if (this%tScreeningInited) then
      this%hPrev(:,:) = 0.0_dp
      this%dRhoPrev(:,:) = 0.0_dp
      this%lrEnergy = 0.0_dp
      this%hfEnergy = 0.0_dp
    end if

  end subroutine updateCoords_cluster


  !> Updates the rangeSep module on coordinate change.
  subroutine updateCoords_gamma(this, env, symNeighbourList, nNeighbourCamSym, skOverCont, orb,&
      & latVecs, recVecs2p)

    !> Class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Symmetric neighbour list version of nNeighbourCam
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Lattice vectors of (periodic) geometry
    real(dp), intent(in) :: latVecs(:,:)

    !> Reciprocal lattice vectors in units of 2pi
    real(dp), intent(in) :: recVecs2p(:,:)

    !! Translation vectors to lattice cells in units of lattice constants
    real(dp), allocatable :: cellVecsG(:,:)

    !! Translation vectors to lattice cells in absolute units
    real(dp), allocatable :: rCellVecsG(:,:)

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Indices of interacting atoms in central cell, as well as their global species index
    integer :: iAtM, iAtN, iSpM, iSpN

    ! range-separated type is initialized with nAtom0 coordinates, therefore re-allocate for
    ! periodic systems, where images beyond the central cell are accounted for
    if (allocated(this%coords)) deallocate(this%coords)
    if (allocated(this%rCoords)) deallocate(this%rCoords)
    this%rCoords = symNeighbourList%coord
    allocate(this%coords(size(this%rCoords, dim=1), size(this%rCoords, dim=2)))

    ! Calculate neighbour list coordinates in relative units
    this%coords(:,:) = this%rCoords
    call cart2frac(this%coords, latVecs)

    ! re-allocate sparse overlap for symmetric neighbour list
    if (allocated(this%overSym)) deallocate(this%overSym)
    allocate(this%overSym(symNeighbourList%sparseSize))

    ! get number of atoms in central cell
    nAtom0 = size(this%species0)

    ! get all cell translations within given cutoff
    call getCellTranslations(cellVecsG, rCellVecsG, latVecs, recVecs2p, this%gSummationCutoff)

    if (this%tLc .or. this%tCam) then
      do iAtM = 1, nAtom0
        iSpM = this%species0(iAtM)
        do iAtN = 1, nAtom0
          iSpN = this%species0(iAtN)
          this%lrGammaEval0(iAtM, iAtN) = getLrGammaGSum(this, iAtM, iAtN, iSpM, iSpN,&
              & rCellVecsG)
          this%lrGammaEval0(iAtN, iAtM) = getLrGammaGSum(this, iAtN, iAtM, iSpN, iSpM,&
              & rCellVecsG)
          this%lrdGammaEval0(iAtM, iAtN, :) = getLrGammaPrimeGSum(this, iAtM, iAtN, iSpM, iSpN,&
              & rCellVecsG)
          this%lrdGammaEval0(iAtN, iAtM, :) = getLrGammaPrimeGSum(this, iAtN, iAtM, iSpN, iSpM,&
              & rCellVecsG)
        end do
      end do
    end if

    if (this%tHyb .or. this%tCam) then
      do iAtM = 1, nAtom0
        iSpM = this%species0(iAtM)
        do iAtN = 1, nAtom0
          iSpN = this%species0(iAtN)
          this%hfGammaEval0(iAtM, iAtN) = getHfGammaGSum(this, iAtM, iAtN, iSpM, iSpN,&
              & rCellVecsG)
          this%hfGammaEval0(iAtN, iAtM) = getHfGammaGSum(this, iAtN, iAtM, iSpN, iSpM,&
              & rCellVecsG)
          this%hfdGammaEval0(iAtM, iAtN, :) = getHfGammaPrimeGSum(this, iAtM, iAtN, iSpM, iSpN,&
              & rCellVecsG)
          this%hfdGammaEval0(iAtN, iAtM, :) = getHfGammaPrimeGSum(this, iAtN, iAtM, iSpN, iSpM,&
              & rCellVecsG)
        end do
      end do
    end if

    ! build symmetric, sparse overlap
    call buildS(env, this%overSym, skOverCont, this%rCoords, nNeighbourCamSym,&
        & symNeighbourList%neighbourList%iNeighbour, symNeighbourList%species,&
        & symNeighbourList%iPair, orb)

    if (this%tScreeningInited) then
      this%hPrev(:,:) = 0.0_dp
      this%dRhoPrev(:,:) = 0.0_dp
      this%lrEnergy = 0.0_dp
      this%hfEnergy = 0.0_dp
    end if

  end subroutine updateCoords_gamma


  !> Interface routine for adding CAM range-separated contributions to the Hamiltonian
  !! (non-periodic).
  subroutine addCamHamiltonian_cluster(this, env, densSqr, over, iNeighbour, nNeighbourCam,&
      & iSquare, iPair, orb, HH, overlap)

    !> Class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    ! Neighbour based screening

    !> Square (unpacked) density matrix
    real(dp), intent(in), target :: densSqr(:,:)

    !> Sparse (packed) overlap matrix.
    real(dp), intent(in) :: over(:)

    !> Neighbour indices.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCam(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> iPair Position of each (neighbour, atom) pair in the sparse matrix.
    !> Shape: (0:maxNeighbour, nAtom)
    integer, intent(in) :: iPair(0:,:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout), target :: HH(:,:)

    ! Threshold based screening

    !> square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    call env%globalTimer%startTimer(globalTimers%rangeSeparatedH)

    ! Add long-range contribution if needed.
    ! For pure Hyb, camBeta would be zero anyway, but we want to save as much time as possible.
    if (this%tLc .or. this%tCam) then
      call addLrHamiltonian_cluster(this, densSqr, over, iNeighbour, nNeighbourCam, iSquare, iPair,&
          & orb, HH, overlap)
    end if

    ! Add full-range Hartree-Fock contribution if needed.
    ! For pure LC, camAlpha would be zero anyway, but we want to save as much time as possible.
    if (this%tHyb .or. this%tCam) then
      call addHfHamiltonian_cluster(this, densSqr, over, iNeighbour, nNeighbourCam,&
          & iSquare, iPair, orb, HH, overlap)
    end if

    call env%globalTimer%stopTimer(globalTimers%rangeSeparatedH)

  end subroutine addCamHamiltonian_cluster


  !> Interface routine for adding CAM range-separated contributions to the Hamiltonian
  !! (periodic systems at the gamma point only).
  subroutine addCamHamiltonian_gamma(this, env, densSqr, symNeighbourList, nNeighbourCamSym,&
      & iSquare, orb, iKS, nKS, HH)

    !> Class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:)

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom (symmetric version)
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Spin/k-point composite index of current diagonalization
    integer, intent(in) :: iKS

    !> Number of Spin/k-point combinations
    integer, intent(in) :: nKS

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout) :: HH(:,:)

    call env%globalTimer%startTimer(globalTimers%rangeSeparatedH)

    ! Add long-range contribution if needed.
    ! For pure Hyb, camBeta would be zero anyway, but we want to save as much time as possible.
    if (this%tLc .or. this%tCam) then
      call addLrHamiltonian_gamma(this, env, densSqr, symNeighbourList, nNeighbourCamSym, iSquare,&
          & orb, iKS, nKS, HH)
    end if

    ! Add full-range Hartree-Fock contribution if needed.
    ! For pure LC, camAlpha would be zero anyway, but we want to save as much time as possible.
    ! if (this%tHyb .or. this%tCam) then
    !   call addHfHamiltonian_gamma()
    ! end if

    call env%globalTimer%stopTimer(globalTimers%rangeSeparatedH)

  end subroutine addCamHamiltonian_gamma


  !> Updates the Hamiltonian with the range separated contribution (k-points version).
  subroutine addCamHamiltonian_kpts(this, env, deltaRhoSqr, deltaRhoSqrCplx, symNeighbourList,&
      & nNeighbourCamSym, iCellVec, rCellVecs, cellVecs, latVecs, recVecs2p, iSquare, orb, kPoint,&
      & kWeight, iKS, iCurSpin, nKS, HSqr)

    !> Class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) delta spin-density matrix at BvK real-space shifts
    real(dp), intent(in), pointer :: deltaRhoSqr(:,:,:,:,:,:)

    !> Square (unpacked) delta spin-density matrix in k-space
    complex(dp), intent(in), pointer :: deltaRhoSqrCplx(:,:)

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Shift vector index for every interacting atom, including periodic images
    integer, intent(in) :: iCellVec(:)

    !> Vectors to neighboring unit cells in absolute units
    real(dp), intent(in) :: rCellVecs(:,:)

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVecs(:,:)

    !> Lattice vectors of (periodic) geometry
    real(dp), intent(in) :: latVecs(:,:)

    !> Reciprocal lattice vectors in units of 2pi
    real(dp), intent(in) :: recVecs2p(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> K-point in relative coordinates to calculate delta H(k) for
    real(dp), intent(in) :: kPoint(:)

    !> K-point weight (for energy contribution)
    real(dp), intent(in) :: kWeight

    !> Spin/k-point composite index of current diagonalization
    integer, intent(in) :: iKS

    !> Spin index of current diagonalization
    integer, intent(in) :: iCurSpin

    !> Number of Spin/k-point combinations
    integer, intent(in) :: nKS

    !> Square (unpacked) Hamiltonian of a certain k-point/spin composite index to be updated
    complex(dp), intent(inout) :: HSqr(:,:)

    call env%globalTimer%startTimer(globalTimers%rangeSeparatedH)

    ! Add long-range contribution if needed.
    ! For pure Hyb, camBeta would be zero anyway, but we want to save as much time as possible.
    if (this%tLc .or. this%tCam) then
      call addLrHamiltonian_kpts(this, env, deltaRhoSqr, deltaRhoSqrCplx, symNeighbourList,&
          & nNeighbourCamSym, iCellVec, rCellVecs, cellVecs, latVecs, recVecs2p, iSquare, orb,&
          & kPoint, kWeight, iKS, iCurSpin, nKS, HSqr)
    end if

    ! Add full-range Hartree-Fock contribution if needed.
    ! For pure LC, camAlpha would be zero anyway, but we want to save as much time as possible.
    ! if (this%tHyb .or. this%tCam) then
    !   call addHfHamiltonian_kpts()
    ! end if

    call env%globalTimer%stopTimer(globalTimers%rangeSeparatedH)

  end subroutine addCamHamiltonian_kpts


  !> Interface routine for adding LC range-separated contributions to the Hamiltonian.
  subroutine addLrHamiltonian_cluster(this, densSqr, over, iNeighbour, nNeighbourCam, iSquare,&
      & iPair, orb, HH, overlap)

    !> Instance
    type(TRangeSepFunc), intent(inout) :: this

    ! Neighbour based screening

    !> Square (unpacked) density matrix
    real(dp), intent(in), target :: densSqr(:,:)

    !> Sparse (packed) overlap matrix.
    real(dp), intent(in) :: over(:)

    !> Neighbour indices.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCam(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> iPair Position of each (neighbour, atom) pair in the sparse matrix.
    !> Shape: (0:maxNeighbour, nAtom)
    integer, intent(in) :: iPair(0:,:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout), target :: HH(:,:)

    ! Threshold based screening

    !> square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    select case(this%rsAlg)
    case (rangeSepTypes%threshold)
      call addLrHamiltonianThreshold_cluster(this, overlap, densSqr, iNeighbour, nNeighbourCam,&
          & iSquare, HH, orb)
    case (rangeSepTypes%neighbour)
      call addLrHamiltonianNeighbour_cluster(this, densSqr, over, iNeighbour, nNeighbourCam,&
          & iSquare, iPair, orb, HH)
    case (rangeSepTypes%matrixBased)
      call addLrHamiltonianMatrix_cluster(this, iSquare, overlap, densSqr, HH)
    end select

  end subroutine addLrHamiltonian_cluster


  !> Interface routine for adding LC range-separated contributions to the Hamiltonian.
  subroutine addLrHamiltonian_gamma(this, env, densSqr, symNeighbourList, nNeighbourCamSym,&
      & iSquare, orb, iKS, nKS, HH)

    !> Instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) density matrix
    real(dp), intent(in), target :: densSqr(:,:)

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom (symmetric version)
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Spin/k-point composite index of current diagonalization
    integer, intent(in) :: iKS

    !> Number of Spin/k-point combinations
    integer, intent(in) :: nKS

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout) :: HH(:,:)

    select case(this%rsAlg)
    case (rangeSepTypes%threshold)
      call error('Thresholded algorithm not implemented for periodic systems.')
    case (rangeSepTypes%neighbour)
      call addLrHamiltonianNeighbour_gamma(this, env, densSqr, symNeighbourList, nNeighbourCamSym,&
          & iSquare, orb, iKS, nKS, HH)
    case (rangeSepTypes%matrixBased)
      call error('Matrix based algorithm not implemented for periodic systems.')
    end select

  end subroutine addLrHamiltonian_gamma


  !> Interface routine for adding LC range-separated contributions to the Hamiltonian.
  subroutine addLrHamiltonian_kpts(this, env, deltaRhoSqr, deltaRhoSqrCplx, symNeighbourList,&
      & nNeighbourCamSym, iCellVec, rCellVecs, cellVecs, latVecs, recVecs2p, iSquare, orb, kPoint,&
      & kWeight, iKS, iCurSpin, nKS, HSqr)

    !> Instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) delta spin-density matrix at BvK real-space shifts
    real(dp), intent(in), pointer :: deltaRhoSqr(:,:,:,:,:,:)

    !> Square (unpacked) delta spin-density matrix in k-space
    complex(dp), intent(in), pointer :: deltaRhoSqrCplx(:,:)

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Shift vector index for every interacting atom, including periodic images
    integer, intent(in) :: iCellVec(:)

    !> Vectors to neighboring unit cells in absolute units
    real(dp), intent(in) :: rCellVecs(:,:)

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVecs(:,:)

    !> Lattice vectors of (periodic) geometry
    real(dp), intent(in) :: latVecs(:,:)

    !> Reciprocal lattice vectors in units of 2pi
    real(dp), intent(in) :: recVecs2p(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> K-point in relative coordinates to calculate delta H(k) for
    real(dp), intent(in) :: kPoint(:)

    !> K-point weight (for energy contribution)
    real(dp), intent(in) :: kWeight

    !> Spin/k-point composite index of current diagonalization
    integer, intent(in) :: iKS

    !> Spin index of current diagonalization
    integer, intent(in) :: iCurSpin

    !> Number of Spin/k-point combinations
    integer, intent(in) :: nKS

    !> Square (unpacked) Hamiltonian of a certain k-point/spin composite index to be updated
    complex(dp), intent(inout) :: HSqr(:,:)

    select case(this%rsAlg)
    case (rangeSepTypes%threshold)
      call error('Thresholded algorithm not implemented for periodic systems.')
    case (rangeSepTypes%neighbour)
      call addLrHamiltonianNeighbour_kpts(this, env, deltaRhoSqr, deltaRhoSqrCplx,&
          & symNeighbourList, nNeighbourCamSym, iCellVec, rCellVecs, cellVecs, latVecs, recVecs2p,&
          & iSquare, orb, kPoint, kWeight, iKS, iCurSpin, nKS, HSqr)
    case (rangeSepTypes%matrixBased)
      call error('Matrix based algorithm not implemented for periodic systems.')
    end select

  end subroutine addLrHamiltonian_kpts


  !> Adds the LR-exchange contribution to hamiltonian using the thresholding algorithm.
  subroutine addLrHamiltonianThreshold_cluster(this, overlap, deltaRho, iNeighbour, nNeighbourCam,&
      & iSquare, hamiltonian, orb)

    !> Class instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    !> Square density matrix (deltaRho in DFTB terms)
    real(dp), intent(in) :: deltaRho(:,:)

    !> Neighbour indices.
    integer, dimension(0:,:), intent(in) :: iNeighbour

    !> Nr. of neighbours for each atom.
    integer, dimension(:), intent(in) :: nNeighbourCam

    !> Mapping atom_number -> number of the first basis function of the atomic block atom_number
    integer, intent(in) :: iSquare(:)

    !> Current hamiltonian
    real(dp), intent(inout) :: hamiltonian(:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    real(dp), allocatable :: tmpOvr(:,:), tmpDRho(:,:), testOvr(:,:), tmpDDRho(:,:), tmpDHam(:,:)
    integer, allocatable :: ovrInd(:,:)
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    call allocateAndInit()
    call evaluateHamiltonian(tmpDHam)
    this%hprev(:,:) = this%hprev + tmpDHam
    hamiltonian(:,:) = hamiltonian + this%camBeta * this%hprev
    this%lrEnergy = this%lrEnergy + evaluateEnergy_real(this%hprev, tmpDRho)

  contains

    !> allocate and initialise some necessary arrays
    subroutine allocateAndInit()

      integer :: matrixSize, nAtom0
      real(dp) :: tmp
      integer :: iAtMu, iAtNu, iNeigh

      matrixSize = size(hamiltonian, dim=1)
      nAtom0 = size(this%species0)

      allocate(tmpOvr(matrixSize, matrixSize))
      tmpOvr(:,:) = overlap
      call blockSymmetrizeHS(tmpOvr, iSquare)

      allocate(tmpDHam(matrixSize, matrixSize))
      tmpDHam(:,:) = 0.0_dp

      allocate(tmpDRho(matrixSize, matrixSize))
      tmpDRho(:,:) = deltaRho
      call symmetrizeHS(tmpDRho)
      call checkAndInitScreening(this, matrixSize, tmpDRho)

      allocate(tmpDDRho(matrixSize, matrixSize))
      tmpDDRho(:,:) = tmpDRho - this%dRhoPrev
      this%dRhoPrev(:,:) = tmpDRho
      allocate(testOvr(nAtom0, nAtom0))
      testOvr(:,:) = 0.0_dp
      allocate(ovrInd(nAtom0, nAtom0))

      do iAtMu = 1, nAtom0
        do iNeigh = 0, nNeighbourCam(iAtMu)
          iAtNu = iNeighbour(iNeigh, iAtMu)
          tmp = maxval(abs(tmpOvr(iSquare(iAtMu) : iSquare(iAtMu + 1) - 1,&
              & iSquare(iAtNu) : iSquare(iAtNu + 1) - 1)))
          testOvr(iAtMu, iAtNu) = tmp
          testOvr(iAtNu, iAtMu) = tmp
        end do
      end do
      do iAtMu = 1, nAtom0
        call index_heap_sort(ovrInd(iAtMu,:), testOvr(iAtMu,:))
      end do

    end subroutine allocateAndInit


    !> Evaluate the update to Hamiltonian due to change in the DM.
    pure subroutine evaluateHamiltonian(tmpDHam)

      !> Update for the old hamiltonian on exit
      real(dp), intent(out) :: tmpDHam(:,:)

      integer :: nAtom0
      real(dp) :: pbound, prb
      real(dp) :: tmpvec1(orb%mOrb), tmpvec2(orb%mOrb)
      real(dp) :: tmp, tstbound, gammabatch, gammabatchtmp
      integer :: iAtMu, iAtNu, iAt1, iAt2, iSp1, iSp2, nOrb1, nOrb2
      integer :: kk, ll, jj, ii, mu, nu
      integer, dimension(descLen) :: desc1, desc2, descM, descN

      nAtom0 = size(this%species0)

      pbound = maxval(abs(tmpDDRho))
      tmpDHam = 0.0_dp
      loopMu: do iAtMu = 1, nAtom0
        descM = getDescriptor(iAtMu, iSquare)
        loopKK: do kk = 1, nAtom0
          iAt1 = ovrInd(iAtMu, nAtom0 + 1 - kk)
          desc1 = getDescriptor(iAt1, iSquare)
          iSp1 = this%species0(iAt1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          prb = pbound * testOvr(iAt1, iAtMu)
          if (abs(prb) < this%pScreeningThreshold) then
            exit loopKK
          end if
          loopNu: do iAtNu = 1, iAtMu
            descN = getDescriptor(iAtNu, iSquare)
            gammabatchtmp = this%lrGammaEval0(iAtMu, iAtNu) + this%lrGammaEval0(iAt1, iAtNu)
            loopLL: do ll = 1, nAtom0
              iAt2 = ovrInd(iAtNu, nAtom0 + 1 - ll)
              iSp2 = this%species0(iAt2)
              nOrb2 = orb%nOrbSpecies(iSp2)
              tstbound = prb * testOvr(iAt2, iAtNu)
              if (abs(tstbound) < this%pScreeningThreshold) then
                exit loopLL
              end if
              desc2 = getDescriptor(iAt2, iSquare)
              gammabatch = (this%lrGammaEval0(iAtMu, iAt2) + this%lrGammaEval0(iAt1, iAt2)&
                  & + gammabatchtmp)
              gammabatch = -0.125_dp * gammabatch
              ! calculate the Q_AB
              do nu = descN(iStart), descN(iEnd)
                jj = 0
                tmpvec2(1:nOrb2) = tmpOvr(desc2(iStart):desc2(iEnd), nu)
                do ii = desc1(iStart), desc1(iEnd)
                  jj = jj + 1
                  tmpvec1(jj) = sum(tmpvec2(1:nOrb2) * tmpDDRho(ii, desc2(iStart):desc2(iEnd)))
                end do
                tmp = 0.0_dp
                do mu = descM(iStart), descM(iEnd)
                  tmp = sum(tmpOvr(desc1(iStart):desc1(iEnd), mu) * tmpvec1(1:nOrb1))
                  tmpDHam(mu, nu) = tmpDHam(mu, nu) + gammabatch * tmp
                end do
              end do
            end do loopLL
          end do loopNu
        end do loopKK
      end do loopMu

    end subroutine evaluateHamiltonian


    !> Initialise the screening matrices.
    subroutine checkAndInitScreening(this, matrixSize, tmpDRho)

      !> Instance
      type(TRangeSepFunc), intent(inout) :: this

      !> linear dimension of matrix
      integer, intent(in) :: matrixSize

      !> Delta rho from iteration
      real(dp), intent(in), allocatable :: tmpDRho(:,:)

      if (.not. this%tScreeningInited) then
        allocate(this%hprev(matrixSize, matrixSize))
        allocate(this%dRhoPrev(matrixSize, matrixSize))
        this%hprev(:,:) = 0.0_dp
        this%dRhoPrev(:,:) = tmpDRho
        this%tScreeningInited = .true.
      end if

    end subroutine checkAndInitScreening

  end subroutine addLrHamiltonianThreshold_cluster


  !> Updates the Hamiltonian with the range separated contribution.
  subroutine addLrHamiltonianNeighbour_cluster(this, densSqr, over, iNeighbour, nNeighbourCam,&
      & iSquare, iPair, orb, HH)

    !> instance of object
    type(TRangeSepFunc), intent(inout) :: this

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:)

    !> Sparse (packed) overlap matrix.
    real(dp), intent(in) :: over(:)

    !> Neighbour indices.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCam(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom0)
    integer, intent(in) :: iSquare(:)

    !> Position of each (neighbour, atom) pair in the sparse matrix. Shape: (0:maxNeighbour, nAtom0)
    integer, intent(in) :: iPair(0:,:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout) :: HH(:,:)

    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3
    real(dp), dimension(orb%mOrb**2) :: Sma, Sam, Snb, Sbn
    real(dp), dimension(orb%mOrb**2) :: Pab, Pmb, Pan, Pmn
    real(dp), dimension(:,:), pointer :: pSma, pSam, pSnb, pSbn
    real(dp), dimension(:,:), pointer :: pPab, pPmb, pPan, pPmn
    real(dp) :: gamma1, gamma2, gammaTot
    integer :: nAtom0
    integer :: iAtM, iAtN, iAtA, iAtB, iNeighN, iNeighA
    integer, dimension(descLen) :: descA, descB, descM, descN
    real(dp), dimension(:,:), allocatable :: tmpDRho
    real(dp), dimension(:,:), allocatable, target :: tmpHH

    call allocateAndInit(tmpHH, tmpDRho)
    call evaluateHamiltonian()

    if (this%tSpin .or. this%tREKS) then
      tmpHH(:,:) = 0.25_dp * tmpHH
    else
      tmpHH(:,:) = 0.125_dp * tmpHH
    end if

    call symmetrizeHS(tmpHH)

    HH(:,:) = HH + this%camBeta * tmpHH
    this%lrEnergy = this%lrEnergy + evaluateEnergy_real(tmpHH, tmpDRho)

  contains

    !> Allocate storage for mapping 1D<->2D array sections
    subroutine allocateAndInit(tmpHH, tmpDRho)

      !> density matrix case
      real(dp), dimension(:,:), allocatable, target, intent(inout) :: tmpDRho

      !> hamiltonian matrix case
      real(dp), dimension(:,:), allocatable, target, intent(inout) :: tmpHH

      allocate(tmpHH(size(HH, dim=1), size(HH, dim=2)))
      tmpHH(:,:) = 0.0_dp
      allocate(tmpDRho(size(densSqr, dim=1), size(densSqr, dim=1)))
      tmpDRho(:,:) = densSqr
      call symmetrizeHS(tmpDRho)

    end subroutine allocateAndInit


    !> actually evaluate the neighbour based cut-off hamiltonian
    subroutine evaluateHamiltonian()

      nAtom0 = size(this%species0)

      loopN: do iAtN = 1, nAtom0
        descN = getDescriptor(iAtN, iSquare)
        loopB: do iNeighN = 0, nNeighbourCam(iAtN)
          iAtB = iNeighbour(iNeighN, iAtN)
          descB = getDescriptor(iAtB, iSquare)
          call copyOverlapBlock(iAtN, iNeighN, descN(iNOrb), descB(iNOrb), Sbn, pSbn)
          call transposeBlock(pSbn, Snb, pSnb)
          loopA: do iAtA = 1, nAtom0
            descA = getDescriptor(iAtA, iSquare)
            call copyDensityBlock(descA, descB, Pab, pPab)
            call copyDensityBlock(descA, descN, Pan, pPan)
            gamma1 = this%lrGammaEval0(iAtA, iAtN) + this%lrGammaEval0(iAtA, iAtB)
            loopM: do iNeighA = 0, nNeighbourCam(iAtA)
              iAtM = iNeighbour(iNeighA, iAtA)
              descM = getDescriptor(iAtM, iSquare)
              call copyOverlapBlock(iAtA, iNeighA, descA(iNOrb), descM(iNOrb), Sma, pSma)
              call transposeBlock(pSma, Sam, pSam)
              gamma2 = this%lrGammaEval0(iAtM, iAtN) + this%lrGammaEval0(iAtM, iAtB)
              gammaTot = gamma1 + gamma2

              if (iAtM >= iAtN) then
                call updateHamiltonianBlock(descM, descN, pSma, pSbn, pPab)
              end if
              if (iAtA >= iAtN .and. iAtM /= iAtA) then
                call copyDensityBlock(descM, descB, Pmb, pPmb)
                call updateHamiltonianBlock(descA, descN, pSam, pSbn, pPmb)
              end if
              if (iAtM >= iAtB .and. iAtN /= iAtB) then
                call updateHamiltonianBlock(descM, descB, pSma, pSnb, pPan)
              end if
              if (iAtA >= iAtB .and. iAtM /= iAtA .and. iAtN /= iAtB) then
                call copyDensityBlock(descM, descN, Pmn, pPmn)
                call updateHamiltonianBlock(descA, descB, pSam, pSnb, pPmn)
              end if
            end do loopM
          end do loopA
        end do loopB
      end do loopN

    end subroutine evaluateHamiltonian


    !> copy atom block from sparse matrix
    pure subroutine copyOverlapBlock(iAt, iNeigh, iNOrbAt, iNOrbNeigh, localBlock, pLocalBlock)

      !> Atom for which this is a neighbour
      integer, intent(in) :: iAt

      !> Number of neighbour for this block
      integer, intent(in) :: iNeigh

      !> Number of orbitals on iAt
      integer, intent(in) :: iNOrbAt

      !> Number of orbitals on neighbour atom
      integer, intent(in) :: iNOrbNeigh

      !> local block
      real(dp), dimension(:), target, intent(inout) :: localBlock

      !> Pointer to local block
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      integer :: ind

      ind = iPair(iNeigh, iAt) + 1
      localBlock(1:iNOrbNeigh*iNOrbAt) = over(ind:ind+iNOrbNeigh*iNOrbAt-1)
      pLocalBlock(1:iNOrbNeigh, 1:iNOrbAt) => localBlock(1:iNOrbNeigh*iNOrbAt)

    end subroutine copyOverlapBlock


    !> copy a density matrix block from sparse matrix
    pure subroutine copyDensityBlock(desc1, desc2, localBlock, pLocalBlock)

      !> start, end and range of first block
      integer, dimension(descLen), intent(in) :: desc1

      !> start, end and range of second block
      integer, dimension(descLen), intent(in) :: desc2

      !> local block in 1D format
      real(dp), dimension(:), target, intent(inout) :: localBlock

      !> Pointer to local block
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      pLocalBlock(1:desc1(iNOrb), 1:desc2(iNOrb)) => localBlock(1:desc1(iNOrb) * desc2(iNOrb))
      pLocalBlock(:,:) = tmpDRho(desc1(iStart):desc1(iEnd), desc2(iStart):desc2(iEnd))

    end subroutine copyDensityBlock


    !> Transpose a block
    pure subroutine transposeBlock(orig, localBlock, pLocalBlock)

      !> Original matrix block
      real(dp), dimension(:,:), intent(in) :: orig

      !> local copy in 1D
      real(dp), dimension(:), target, intent(out) :: localBlock

      !> pointer to local copy
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      pLocalBlock(1:size(orig, dim=2), 1:size(orig, dim=1)) => localBlock(1:size(orig))
      pLocalBlock = transpose(orig)


    end subroutine transposeBlock


    !> Add a contribution to a Hamiltonian block
    subroutine updateHamiltonianBlock(descM, descN, pSma, pSbN, pPab)

      !> start, end and range of row
      integer, dimension(descLen), intent(in) :: descM

      !> start, end and range of column
      integer, dimension(descLen), intent(in) :: descN

      !> First overlap block
      real(dp), dimension(:,:), pointer, intent(in) :: pSma

      !> Second overlap block
      real(dp), dimension(:,:), pointer, intent(in) :: pSbN

      !> density matrix block
      real(dp), dimension(:,:), pointer, intent(in) :: pPab

      real(dp), dimension(:,:), pointer :: pHmn

      pHmn => tmpHH(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd))
      pHmn(:,:) = pHmn - gammaTot * matmul(matmul(pSma, pPab), pSbn)

    end subroutine updateHamiltonianBlock

  end subroutine addLrHamiltonianNeighbour_cluster


  !> Update Hamiltonian with long-range contribution using matrix-matrix multiplications
  !>
  !> The routine provides a matrix-matrix multiplication based implementation of
  !> the 3rd term in Eq. 26 in https://doi.org/10.1063/1.4935095
  !>
  subroutine addLrHamiltonianMatrixCmplx(this, iSquare, overlap, densSqr, HH)

    !> Class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, dimension(:), intent(in) :: iSquare

    !> Square (unpacked) overlap matrix.
    complex(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    complex(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated.
    complex(dp), intent(inout) :: HH(:,:)

    complex(dp), allocatable :: Smat(:,:)
    complex(dp), allocatable :: Dmat(:,:)
    real(dp), allocatable :: lrGammaAO(:,:)
    complex(dp), allocatable :: gammaCmplx(:,:)
    complex(dp), allocatable :: Hlr(:,:)

    integer :: nOrb

    nOrb = size(overlap,dim=1)

    allocate(Smat(nOrb, nOrb))
    allocate(Dmat(nOrb, nOrb))
    allocate(lrGammaAO(nOrb, nOrb))
    allocate(gammaCmplx(nOrb, nOrb))
    allocate(Hlr(nOrb, nOrb))

    call allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, lrGammaAO, gammaCmplx)

    call evaluateHamiltonian(this, Smat, Dmat, gammaCmplx, Hlr)

    HH(:,:) = HH + this%camBeta * Hlr

    this%lrenergy = this%lrenergy + 0.5_dp * real(sum(Dmat * Hlr), dp)

  contains

    subroutine allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, lrGammaAO,&
        & gammaCmplx)

      !> instance
      type(TRangeSepFunc), intent(inout) :: this

      !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
      integer, dimension(:), intent(in) :: iSquare

      !> Square (unpacked) overlap matrix.
      complex(dp), intent(in) :: overlap(:,:)

      !> Square (unpacked) density matrix
      complex(dp), intent(in) :: densSqr(:,:)

      !> Square (unpacked) Hamiltonian to be updated.
      complex(dp), intent(inout) :: HH(:,:)

      !> Symmetrized square overlap matrix
      complex(dp), intent(out) :: Smat(:,:)

      !> Symmetrized square density matrix
      complex(dp), intent(out) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      real(dp), intent(out) :: lrGammaAO(:,:)

      !> Symmetrized long-range gamma matrix
      complex(dp), intent(out) :: gammaCmplx(:,:)

      integer :: nAtom, iAt, jAt

      nAtom = size(this%lrGammaEval0, dim=1)

      !! Symmetrize Hamiltonian, overlap, density matrices
      call hermitianSquareMatrix(HH)
      Smat(:,:) = overlap
      call hermitianSquareMatrix(Smat)
      Dmat(:,:) = densSqr
      call hermitianSquareMatrix(Dmat)

      ! Get long-range gamma variable
      lrGammaAO(:,:) = 0.0_dp
      do iAt = 1, nAtom
        do jAt = 1, nAtom
          lrGammaAO(iSquare(jAt):iSquare(jAt+1)-1,iSquare(iAt):iSquare(iAt+1)-1) =&
              & this%lrGammaEval0(jAt, iAt)
        end do
      end do
      gammaCmplx(:,:) = lrGammaAO

    end subroutine allocateAndInit


    subroutine evaluateHamiltonian(this, Smat, Dmat, gammaCmplx, Hlr)

      !> instance
      type(TRangeSepFunc), intent(inout) :: this

      !> Symmetrized square overlap matrix
      complex(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      complex(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      complex(dp), intent(in) :: gammaCmplx(:,:)

      !> Symmetrized long-range Hamiltonian matrix
      complex(dp), intent(out) :: Hlr(:,:)

      complex(dp), allocatable :: Hmat(:,:)
      complex(dp), allocatable :: tmpMat(:,:)

      integer :: nOrb

      nOrb = size(Smat,dim=1)

      allocate(Hmat(nOrb,nOrb))
      allocate(tmpMat(nOrb,nOrb))

      Hlr(:,:) = cmplx(0.0_dp,0.0_dp,dp)

      call gemm(tmpMat, Smat, Dmat)
      call gemm(Hlr, tmpMat, Smat)
      Hlr(:,:) = Hlr * gammaCmplx

      tmpMat(:,:) = tmpMat * gammaCmplx
      call gemm(Hlr, tmpMat, Smat, alpha=(1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))

      Hmat(:,:) = Dmat * gammaCmplx
      call gemm(tmpMat, Smat, Hmat)
      call gemm(Hlr, tmpMat, Smat, alpha=(1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))

      call gemm(tmpMat, Dmat, Smat)
      tmpMat(:,:) = tmpMat * gammaCmplx
      call gemm(Hlr, Smat, tmpMat, alpha=(1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))

      if (this%tSpin) then
        Hlr(:,:) = -0.25_dp * Hlr
      else
        Hlr(:,:) = -0.125_dp * Hlr
      end if

    end subroutine evaluateHamiltonian

  end subroutine addLrHamiltonianMatrixCmplx


  !> Update Hamiltonian with long-range contribution using matrix-matrix multiplications
  !>
  !> The routine provides a matrix-matrix multiplication based implementation of
  !> the 3rd term in Eq. 26 in https://doi.org/10.1063/1.4935095
  !>
  subroutine addLrHamiltonianMatrix_cluster(this, iSquare, overlap, densSqr, HH)

    !> Class instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, dimension(:), intent(in) :: iSquare

    !> Square (unpacked) overlap matrix.
    real(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout) :: HH(:,:)

    real(dp), allocatable :: Smat(:,:)
    real(dp), allocatable :: Dmat(:,:)
    real(dp), allocatable :: lrGammaAO(:,:)
    real(dp), allocatable :: Hlr(:,:)

    integer :: nOrb

    nOrb = size(overlap,dim=1)

    allocate(Smat(nOrb, nOrb))
    allocate(Dmat(nOrb, nOrb))
    allocate(lrGammaAO(nOrb, nOrb))
    allocate(Hlr(nOrb, nOrb))

    call allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, lrGammaAO)
    call evaluateHamiltonian(this, Smat, Dmat, lrGammaAO, Hlr)
    HH(:,:) = HH + this%camBeta * Hlr
    this%lrEnergy = this%lrEnergy + 0.5_dp * sum(Dmat * Hlr)

  contains

    !> Set up storage and get orbital-by-orbital gamma matrix
    subroutine allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, lrGammaAO)

      !> Class instance
      type(TRangeSepFunc), intent(inout) :: this

      !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
      integer, dimension(:), intent(in) :: iSquare

      !> Square (unpacked) overlap matrix.
      real(dp), intent(in) :: overlap(:,:)

      !> Square (unpacked) density matrix
      real(dp), intent(in) :: densSqr(:,:)

      !> Square (unpacked) Hamiltonian to be updated.
      real(dp), intent(inout) :: HH(:,:)

      !> Symmetrized square overlap matrix
      real(dp), intent(out) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(out) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      real(dp), intent(out) :: lrGammaAO(:,:)

      integer :: nAtom, iAt, jAt

      nAtom = size(this%lrGammaEval0, dim=1)

      ! Symmetrize Hamiltonian, overlap, density matrices
      call symmetrizeHS(HH)
      Smat(:,:) = overlap
      call symmetrizeHS(Smat)
      Dmat(:,:) = densSqr
      call symmetrizeHS(Dmat)

      ! Get long-range gamma variable
      lrGammaAO(:,:) = 0.0_dp
      do iAt = 1, nAtom
        do jAt = 1, nAtom
          lrGammaAO(iSquare(jAt):iSquare(jAt+1)-1,iSquare(iAt):iSquare(iAt+1)-1) =&
              & this%lrGammaEval0(jAt, iAt)
        end do
      end do

    end subroutine allocateAndInit


    !> Evaluate the hamiltonian using GEMM operations
    subroutine evaluateHamiltonian(this, Smat, Dmat, lrGammaAO, Hlr)

      !> Class instance
      type(TRangeSepFunc), intent(inout) :: this

      !> Symmetrized square overlap matrix
      real(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      real(dp), intent(in) :: lrGammaAO(:,:)

      !> Symmetrized long-range Hamiltonian matrix
      real(dp), intent(out) :: Hlr(:,:)

      real(dp), allocatable :: Hmat(:,:)
      real(dp), allocatable :: tmpMat(:,:)

      integer :: nOrb

      nOrb = size(Smat, dim=1)

      allocate(Hmat(nOrb, nOrb))
      allocate(tmpMat(nOrb, nOrb))

      Hlr(:,:) = 0.0_dp

      call gemm(tmpMat, Smat, Dmat)
      call gemm(Hlr, tmpMat, Smat)
      Hlr(:,:) = Hlr * lrGammaAO

      tmpMat(:,:) = tmpMat * lrGammaAO
      call gemm(Hlr, tmpMat, Smat, alpha=1.0_dp, beta=1.0_dp)

      Hmat(:,:) = Dmat * lrGammaAO
      call gemm(tmpMat, Smat, Hmat)
      call gemm(Hlr, tmpMat, Smat, alpha=1.0_dp, beta=1.0_dp)

      call gemm(tmpMat, Dmat, Smat)
      tmpMat(:,:) = tmpMat * lrGammaAO
      call gemm(Hlr, Smat, tmpMat, alpha=1.0_dp, beta=1.0_dp)

      if (this%tSpin .or. this%tREKS) then
        Hlr(:,:) = -0.25_dp * Hlr
      else
        Hlr(:,:) = -0.125_dp * Hlr
      end if

    end subroutine evaluateHamiltonian

  end subroutine addLrHamiltonianMatrix_cluster


  !> Updates the Hamiltonian with the range separated contribution.
  subroutine addLrHamiltonianNeighbour_gamma(this, env, deltaRhoSqr, symNeighbourList,&
      & nNeighbourCamSym, iSquare, orb, iKS, nKS, HSqr)

    !> Instance
    type(TRangeSepFunc), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) delta density matrix
    real(dp), intent(in) :: deltaRhoSqr(:,:)

    !> list of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Spin/k-point composite index of current diagonalization
    integer, intent(in) :: iKS

    !> Number of Spin/k-point combinations
    integer, intent(in) :: nKS

    !> Square (unpacked) Hamiltonian to be updated
    real(dp), intent(inout) :: HSqr(:,:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Atom blocks from sparse, real-space overlap matrices S_{\alpha\mu}, S_{\beta\nu}
    real(dp), pointer :: pSam(:,:), pSbn(:,:)
    real(dp), allocatable :: pSamT(:,:)

    !! Density matrix block from sparse matrix
    real(dp), allocatable :: Pab(:,:)

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descA(descLen), descB(descLen), descM(descLen), descN(descLen)

    !! Temporary storages
    real(dp), allocatable :: tmpDeltaRhoSqr(:,:), tmpDeltaDeltaRhoSqr(:,:), tmpHSqr(:,:)
    type(TRealArray1D), allocatable :: testSquareOver(:)
    type(TIntArray1D), allocatable :: overlapIndices(:)

    !! Temporary arrays for matrix-matrix operations
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab, pSamT_Pab_pSbn, Pab_Sbn
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab_gammaAB, tot

    !! \tilde{\gamma}_{\mu\nu}, \tilde{\gamma}_{\mu\beta},
    !! \tilde{\gamma}_{\alpha\nu}, \tilde{\gamma}_{\alpha\beta}
    real(dp) :: gammaMNMB, gammaTot

    !! Overlap matrix elements
    real(dp) :: Sam, Sbn

    !! Density matrix elements
    real(dp) :: dPmn, dPab

    !! Orbital indices
    integer :: mu, nu, alpha, beta

    !! Product Sam * dPab * Sbn * gammaTot
    real(dp) :: SamdPabSbnGammaTot

    !! Atom indices (central cell)
    integer :: iAtM, iAtN

    !! Neighbour indices (+corresponding atom indices)
    integer :: iNeighN, iNeighM, iAtA, iAtB

    !! Folded (to central cell) atom indices
    integer :: iAtAfold, iAtBfold

    !! Sorted (according to max overlap estimates) neighbour indices
    integer :: iNeighMsort, iNeighNsort

    !! Auxiliary variables for setting up 2D pointer to sparse overlap
    integer :: ind, nOrbAt, nOrbNeigh

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Size of square matrices (e.g. Hamiltonian)
    integer :: squareSize

    !! Max estimate for difference of square delta rho to previous SCC iteration and products with
    !! max overlap estimates
    real(dp) :: pMax, pMaxSbn, pMaxSbnSam, maxEstimate, pSbnMax, pMaxpSbnMax

    !! Start and end index for MPI parallelization, if applicable
    integer :: iParallelStart, iParallelEnd

    !! Composite index iAtM/iAtN
    integer :: ii, iAtMN(2, size(this%species0)**2)

    squareSize = size(HSqr, dim=1)
    nAtom0 = size(this%species0)

    ! Build up composite index iAtMN for collapsing iAtM and iAtN
    ind = 1
    loopM: do iAtM = 1, nAtom0
      loopN: do iAtN = 1, nAtom0
        iAtMN(1, ind) = iAtM
        iAtMN(2, ind) = iAtN
        ind = ind + 1
      end do loopN
    end do loopM

  #:if WITH_MPI
    call getStartAndEndIndex(nAtom0**2, env%mpi%groupComm%size, env%mpi%groupComm%rank,&
        & iParallelStart, iParallelEnd)
  #:else
    iParallelStart = 1
    iParallelEnd = nAtom0**2
  #:endif

    ! allocate delta Hamiltonian
    allocate(tmpHSqr(squareSize, squareSize))
    tmpHSqr(:,:) = 0.0_dp

    tmpDeltaRhoSqr = deltaRhoSqr
    call symmetrizeHS(tmpDeltaRhoSqr)

    ! check and initialize screening
    if (.not. this%tScreeningInited) then
      if (iKS == 1) then
        allocate(this%hPrev(squareSize, squareSize))
        this%hPrev(:,:) = 0.0_dp
        this%dRhoPrev = tmpDeltaRhoSqr
      end if
      if (iKS == nKS) this%tScreeningInited = .true.
      ! there is no previous delta density matrix, therefore just copy over
      tmpDeltaDeltaRhoSqr = tmpDeltaRhoSqr
    else
      ! allocate and initialize difference of delta rho to previous SCC iteration
      tmpDeltaDeltaRhoSqr = tmpDeltaRhoSqr - this%dRhoPrev
    end if

    pMax = maxval(abs(tmpDeltaDeltaRhoSqr))
    ! store delta density matrix only once per spin and SCC iteration
    if (iKS == nKS) then
      this%dRhoPrev(:,:) = tmpDeltaRhoSqr
    end if
    ! skip whole procedure if delta density matrix is close to zero, e.g. in the first SCC iteration
    if (pMax < 1e-16_dp) return

    ! allocate max estimates of square overlap blocks and index array for sorting
    allocate(testSquareOver(nAtom0))
    allocate(overlapIndices(nAtom0))

    do iAtN = 1, nAtom0
      descN = getDescriptor(iAtN, iSquare)
      allocate(testSquareOver(iAtN)%array(nNeighbourCamSym(iAtN) + 1))
      do iNeighN = 0, nNeighbourCamSym(iAtN)
        iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighN, iAtN)
        iAtBfold = symNeighbourList%img2CentCell(iAtB)
        descB = getDescriptor(iAtBfold, iSquare)
        ! get 2D pointer to Sbn overlap block
        ind = symNeighbourList%iPair(iNeighN, iAtN) + 1
        nOrbAt = descN(iNOrb)
        nOrbNeigh = descB(iNOrb)
        pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
        testSquareOver(iAtN)%array(iNeighN + 1) = maxval(abs(pSbn))
      end do
    end do

    ! sort max square overlap estimates (descending)
    ! this way we can exit the whole loop a.s.a. threshold has been undershot for the first time
    do iAtN = 1, nAtom0
      allocate(overlapIndices(iAtN)%array(nNeighbourCamSym(iAtN) + 1))
      call index_heap_sort(overlapIndices(iAtN)%array, testSquareOver(iAtN)%array)
      ! switch from ascending to descending
      overlapIndices(iAtN)%array(:)&
          & = overlapIndices(iAtN)%array(size(overlapIndices(iAtN)%array):1:-1)
    end do

    loopMN: do ii = iParallelStart, iParallelEnd
      iAtM = iAtMN(1, ii)
      iAtN = iAtMN(2, ii)
      descM = getDescriptor(iAtM, iSquare)
      descN = getDescriptor(iAtN, iSquare)
      loopB: do iNeighN = 0, nNeighbourCamSym(iAtN)
        iNeighNsort = overlapIndices(iAtN)%array(iNeighN + 1) - 1
        pSbnMax = testSquareOver(iAtN)%array(iNeighNsort + 1)
        pMaxpSbnMax = pMax * pSbnMax
        if (pMaxpSbnMax < this%pScreeningThreshold) exit loopB
        iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighNsort, iAtN)
        iAtBfold = symNeighbourList%img2CentCell(iAtB)
        descB = getDescriptor(iAtBfold, iSquare)
        ! get 2D pointer to S_{\beta\nu}(\vec{l}) overlap block
        ind = symNeighbourList%iPair(iNeighNsort, iAtN) + 1
        nOrbAt = descN(iNOrb)
        nOrbNeigh = descB(iNOrb)
        pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
        ! \tilde{\gamma}_{\mu\nu} + \tilde{\gamma}_{\mu\beta}
        gammaMNMB = this%lrGammaEval0(iAtM, iAtN) + this%lrGammaEval0(iAtM, iAtBfold)
        loopA: do iNeighM = 0, nNeighbourCamSym(iAtM)
          iNeighMsort = overlapIndices(iAtM)%array(iNeighM + 1) - 1
          maxEstimate = pMaxpSbnMax * testSquareOver(iAtM)%array(iNeighMsort + 1)
          if (maxEstimate < this%pScreeningThreshold) exit loopA
          iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighMsort, iAtM)
          iAtAfold = symNeighbourList%img2CentCell(iAtA)
          descA = getDescriptor(iAtAfold, iSquare)
          ! get continuous 2D copy of Pab density matrix block
          Pab = tmpDeltaDeltaRhoSqr(descA(iStart):descA(iEnd), descB(iStart):descB(iEnd))
          ! get 2D pointer to S_{\alpha\mu}(\vec{h}) overlap block
          ind = symNeighbourList%iPair(iNeighMsort, iAtM) + 1
          nOrbAt = descM(iNOrb)
          nOrbNeigh = descA(iNOrb)
          pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

          gammaTot = gammaMNMB + this%lrGammaEval0(iAtAfold, iAtN)&
              & + this%lrGammaEval0(iAtAfold, iAtBfold)

          do mu = 1, descM(iNOrb)
            do nu = 1, descN(iNOrb)
              do beta = 1, descB(iNOrb)
                Sbn = pSbn(beta, nu)
                do alpha = 1, descA(iNOrb)
                  Sam = pSam(alpha, mu)
                  dPab = Pab(alpha, beta)
                  SamdPabSbnGammaTot = Sam * dPab * Sbn * gammaTot

                  tmpHSqr(descM(iStart) + mu - 1, descN(iStart) + nu - 1)&
                      & = tmpHSqr(descM(iStart) + mu - 1, descN(iStart) + nu - 1)&
                      & + SamdPabSbnGammaTot
                end do
              end do
            end do
          end do

          ! \transpose{S_{\alpha\mu}(\vec{-h})}
          ! pSamT = transpose(pSam(1:nOrbNeigh, 1:nOrbAt))

          ! pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) = matmul(pSamT, Pab)
          ! Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) = matmul(Pab, pSbn)

          ! ! term #1
          ! pSamT_Pab_pSbn(1:descM(iNOrb), 1:descN(iNOrb))&
          !     & = matmul(pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)), pSbn)
          ! tot(1:descM(iNOrb), 1:descN(iNOrb)) = pSamT_Pab_pSbn(1:descM(iNOrb), 1:descN(iNOrb))&
          !     & * gammaMN

          ! ! term #2
          ! tot(1:descM(iNOrb), 1:descN(iNOrb)) = tot(1:descM(iNOrb), 1:descN(iNOrb))&
          !     & + matmul(pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) * gammaMB, pSbn)

          ! ! term #3
          ! tot(1:descM(iNOrb), 1:descN(iNOrb)) = tot(1:descM(iNOrb), 1:descN(iNOrb))&
          !     & + matmul(pSamT, Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) * gammaAN)

          ! ! term #4
          ! pSamT_Pab_gammaAB(1:descM(iNorb), 1:descB(iNorb)) = matmul(pSamT, Pab * gammaAB)
          ! tot(1:descM(iNOrb), 1:descN(iNOrb)) = tot(1:descM(iNOrb), 1:descN(iNOrb))&
          !     & + matmul(pSamT_Pab_gammaAB(1:descM(iNOrb), 1:descB(iNOrb)), pSbn)

          ! tmpHSqr(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd))&
          !     & = tmpHSqr(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd))&
          !     & + tot(1:descM(iNOrb), 1:descN(iNOrb))

        end do loopA
      end do loopB
    end do loopMN

    ! tmpHSqr(:,:) = 0.0_dp
    if (this%tSpin .or. this%tREKS) then
      tmpHSqr(:,:) = -0.25_dp * tmpHSqr
    else
      tmpHSqr(:,:) = -0.125_dp * tmpHSqr
    end if

  #:if WITH_MPI
    ! Sum up contributions of current MPI group
    call mpifx_allreduceip(env%mpi%groupComm, tmpHSqr, MPI_SUM)
  #:endif

    this%hprev(:,:) = this%hprev + tmpHSqr
    HSqr(:,:) = HSqr + this%camBeta * this%hprev

    ! Add energy contribution but divide by the number of processes
  #:if WITH_MPI
    this%lrEnergy = this%lrEnergy + evaluateEnergy_real(this%hprev, tmpDeltaRhoSqr)&
        & / real(env%mpi%groupComm%size, dp)
  #:else
    this%lrEnergy = this%lrEnergy + evaluateEnergy_real(this%hprev, tmpDeltaRhoSqr)
  #:endif

  end subroutine addLrHamiltonianNeighbour_gamma


  !> Updates the Hamiltonian with the range separated contribution (k-points version).
  subroutine addLrHamiltonianNeighbour_kpts(this, env, deltaRhoSqr, deltaRhoOutSqrCplx,&
      & symNeighbourList, nNeighbourCamSym, iCellVec, rCellVecs, cellVecs, latVecs, recVecs2p,&
      & iSquare, orb, kPoint, kWeight, iKS, iCurSpin, nKS, HSqr)

    !> Instance
    type(TRangeSepFunc), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) delta spin-density matrix at BvK real-space shifts
    real(dp), intent(in), pointer :: deltaRhoSqr(:,:,:,:,:,:)

    !> Square (unpacked) delta spin-density matrix in k-space
    complex(dp), intent(in), pointer :: deltaRhoOutSqrCplx(:,:)

    !> list of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Shift vector index for every interacting atom, including periodic images
    integer, intent(in) :: iCellVec(:)

    !> Vectors to neighboring unit cells in absolute units
    real(dp), intent(in) :: rCellVecs(:,:)

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVecs(:,:)

    !> Lattice vectors of (periodic) geometry
    real(dp), intent(in) :: latVecs(:,:)

    !> Reciprocal lattice vectors in units of 2pi
    real(dp), intent(in) :: recVecs2p(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> K-point in relative coordinates to calculate delta H(k) for
    real(dp), intent(in) :: kPoint(:)

    !> K-point weight (for energy contribution)
    real(dp), intent(in) :: kWeight

    !> Spin/k-point composite index of current diagonalization
    integer, intent(in) :: iKS

    !> Spin index of current diagonalization
    integer, intent(in) :: iCurSpin

    !> Number of Spin/k-point combinations
    integer, intent(in) :: nKS

    !> Square (unpacked) Hamiltonian of a certain k-point/spin composite index to be updated
    complex(dp), intent(inout) :: HSqr(:,:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Atom blocks from sparse, real-space overlap matrices S_{\alpha\mu}, S_{\beta\nu}
    real(dp), pointer :: pSam(:,:), pSbn(:,:)
    real(dp), allocatable :: pSamT(:,:)

    !! Density matrix block \alpha\beta
    real(dp), allocatable :: Pab(:,:,:,:,:)

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descA(descLen), descB(descLen), descM(descLen), descN(descLen)

    !! Temporary storages
    real(dp), allocatable :: deltaDeltaRhoSqr(:,:,:,:,:)
    complex(dp), allocatable :: tmpHSqr(:,:)
    type(TRealArray1D), allocatable :: testSquareOver(:)
    type(TIntArray1D), allocatable :: overlapIndices(:)

    !! Index arrays for descending sorting of Gamma arrays
    integer, allocatable :: gammaMNsortIdx(:), gammaMBsortIdx(:)
    integer, allocatable :: gammaANsortIdx(:), gammaABsortIdx(:)

    !! Number of non-vanishing elements in descending Gamma arrays
    integer :: nNonZeroGammaMN, nNonZeroGammaMB, nNonZeroGammaAN, nNonZeroGammaAB

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Real-space \vec{h} and \vec{l} vectors in absolute and relative coordinates
    real(dp) :: rVecH(3), rVecL(3), vecH(3), vecL(3)

    !! Translation vectors to lattice cells in units of lattice constants
    real(dp), allocatable :: cellVecsG(:,:)

    !! Vectors to unit cells in absolute units
    real(dp), allocatable :: rCellVecsG(:,:)

    !! Temporary arrays for gemm operations
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab, pSamT_Pab_pSbn, Pab_Sbn
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab_gammaAB, tot

    !! \gamma_{\mu\nu}(\vec{g}), \gamma_{\mu\beta}(\vec{g} - \vec{l}),
    !! \gamma_{\alpha\nu}(\vec{g} - \vec{h}), \gamma_{\alpha\beta}(\vec{g} - \vec{l} - \vec{h})
    real(dp), allocatable :: gammaMN(:), gammaMB(:), gammaAN(:), gammaAB(:)

    !! Species of atom where orbitals \alpha, \beta, \mu and \nu are located
    integer :: iSpA, iSpB, iSpM, iSpN

    !! Atom indices (central cell)
    integer :: iAtM, iAtN

    !! Neighbour indices (+corresponding atom indices)
    integer :: iNeighN, iNeighM, iAtA, iAtB

    !! Folded (to central cell) atom indices
    integer :: iAtAfold, iAtBfold

    !! Sorted (according to max overlap estimates) neighbour indices
    integer :: iNeighMsort, iNeighNsort

    !! Auxiliary variables for setting up 2D pointer to sparse overlap
    integer :: ind, nOrbAt, nOrbNeigh

    !! Size of square matrices (e.g. Hamiltonian)
    integer :: squareSize

    !! Max estimate for difference of square delta rho to previous SCC iteration and products with
    !! max overlap estimates
    real(dp) :: pMax, pMaxSbn, pMaxSbnSam, maxEstimate, pSbnMax, pMaxpSbnMax

    !! Max. number of neighbours in neighbour list
    integer :: nNeigh

    !! Integer BvK real-space shift in relative coordinates
    integer :: bvKShift(3), bvKIndex(3)

    !! Phase factor
    complex(dp) :: phase

    !! Iterates over all BvK real-space vectors
    integer :: iG, iGMN, iGMB, iGAN, iGAB

    !! Dummy array with zeros
    real(dp) :: zeros(3)

    !! Start and end index for MPI parallelization, if applicable
    integer :: iParallelStart, iParallelEnd

    !! Composite index iAtM/iAtN
    integer :: ii, iAtMN(2, size(this%species0)**2)

    zeros(:) = 0.0_dp

    ! this initialization is also valid for OMP, since tot is firstprivate
    tot(:,:) = 0.0_dp

    squareSize = size(HSqr, dim=1)
    nAtom0 = size(this%species0)

    ! Build up composite index iAtMN for collapsing iAtM and iAtN
    ind = 1
    loopM: do iAtM = 1, nAtom0
      loopN: do iAtN = 1, nAtom0
        iAtMN(1, ind) = iAtM
        iAtMN(2, ind) = iAtN
        ind = ind + 1
      end do loopN
    end do loopM

  #:if WITH_MPI
    call getStartAndEndIndex(nAtom0**2, env%mpi%groupComm%size, env%mpi%groupComm%rank,&
        & iParallelStart, iParallelEnd)
  #:else
    iParallelStart = 1
    iParallelEnd = nAtom0**2
  #:endif

    ! allocate delta Hamiltonian
    allocate(tmpHSqr(squareSize, squareSize))
    tmpHSqr(:,:) = (0.0_dp, 0.0_dp)

    ! check and initialize screening
    if (.not. this%tScreeningInited) then
      if (iKS == 1) then
        allocate(this%hprevCplxHS(squareSize, squareSize, nKS))
        this%hprevCplxHS(:,:,:) = cmplx(0, 0, dp)
        this%dRhoPrevCplxHS = deltaRhoSqr
      end if
      if (iKS == nKS .or. iCurSpin == 2) this%tScreeningInited = .true.
      ! there is no previous delta density matrix, therefore just copy over
      deltaDeltaRhoSqr = deltaRhoSqr(:,:,:,:,:, iCurSpin)
    else
      ! allocate and initialize difference of delta rho to previous SCC iteration
      deltaDeltaRhoSqr = deltaRhoSqr(:,:,:,:,:, iCurSpin) - this%dRhoPrevCplxHS(:,:,:,:,:, iCurSpin)
    end if

    pMax = maxval(abs(deltaDeltaRhoSqr))
    ! store delta density matrix only once per spin and SCC iteration
    if (iKS == nKS .or. iCurSpin == 2) then
      this%dRhoPrevCplxHS(:,:,:,:,:, iCurSpin) = deltaRhoSqr(:,:,:,:,:, iCurSpin)
    end if
    ! skip whole procedure if delta density matrix is close to zero, e.g. in the first SCC iteration
    if (pMax < 1e-16_dp) return

    ! allocate max estimates of square overlap blocks and index array for sorting
    allocate(testSquareOver(nAtom0))
    allocate(overlapIndices(nAtom0))

    do iAtN = 1, nAtom0
      descN = getDescriptor(iAtN, iSquare)
      allocate(testSquareOver(iAtN)%array(nNeighbourCamSym(iAtN) + 1))
      do iNeighN = 0, nNeighbourCamSym(iAtN)
        iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighN, iAtN)
        iAtBfold = symNeighbourList%img2CentCell(iAtB)
        descB = getDescriptor(iAtBfold, iSquare)
        ! get 2D pointer to Sbn overlap block
        ind = symNeighbourList%iPair(iNeighN, iAtN) + 1
        nOrbAt = descN(iNOrb)
        nOrbNeigh = descB(iNOrb)
        pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
        testSquareOver(iAtN)%array(iNeighN + 1) = maxval(abs(pSbn))
      end do
    end do

    ! sort max square overlap estimates (descending)
    ! this way we can exit the whole loop a.s.a. threshold has been undershot for the first time
    do iAtN = 1, nAtom0
      allocate(overlapIndices(iAtN)%array(nNeighbourCamSym(iAtN) + 1))
      call index_heap_sort(overlapIndices(iAtN)%array, testSquareOver(iAtN)%array)
      ! switch from ascending to descending
      overlapIndices(iAtN)%array(:)&
          & = overlapIndices(iAtN)%array(size(overlapIndices(iAtN)%array):1:-1)
    end do

    ! get all cell translations within given cutoff
    call getCellTranslations(cellVecsG, rCellVecsG, latVecs, recVecs2p, this%gSummationCutoff)

    allocate(gammaMN(size(cellVecsG, dim=2)))
    allocate(gammaMB(size(cellVecsG, dim=2)))
    allocate(gammaAN(size(cellVecsG, dim=2)))
    allocate(gammaAB(size(cellVecsG, dim=2)))

    ! re-indexing for descending sorting
    allocate(gammaMNsortIdx(size(cellVecsG, dim=2)))
    allocate(gammaMBsortIdx(size(cellVecsG, dim=2)))
    allocate(gammaANsortIdx(size(cellVecsG, dim=2)))
    allocate(gammaABsortIdx(size(cellVecsG, dim=2)))

    loopMN: do ii = 1, nAtom0**2
      iAtM = iAtMN(1, ii)
      iAtN = iAtMN(2, ii)
      iSpM = this%species0(iAtM)
      descM = getDescriptor(iAtM, iSquare)
      iSpN = this%species0(iAtN)
      descN = getDescriptor(iAtN, iSquare)
      ! pre-tabulate g-resolved \gamma_{\mu\nu}(\vec{g})
      gammaMN(:) = getGammaGResolved(this, iAtM, iAtN, iSpM, iSpN, this%coords, this%rCoords,&
          & cellVecsG, rCellVecsG, zeros, zeros)
      call index_heap_sort(gammaMNsortIdx, gammaMN)
      gammaMNsortIdx(:) = gammaMNsortIdx(size(gammaMNsortIdx):1:-1)
      nNonZeroGammaMN = getNumberOfNonZeroElements(gammaMN(gammaMNsortIdx),&
          & this%pScreeningThreshold)
      loopB: do iNeighN = 0, nNeighbourCamSym(iAtN)
        iNeighNsort = overlapIndices(iAtN)%array(iNeighN + 1) - 1
        pSbnMax = testSquareOver(iAtN)%array(iNeighNsort + 1)
        pMaxpSbnMax = pMax * pSbnMax
        if (pMaxpSbnMax < this%pScreeningThreshold) exit loopB
        iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighNsort, iAtN)
        iAtBfold = symNeighbourList%img2CentCell(iAtB)
        iSpB = this%species0(iAtBfold)
        descB = getDescriptor(iAtBfold, iSquare)
        ! get real-space \vec{l} for gamma arguments
        rVecL(:) = -rCellVecs(:, symNeighbourList%iCellVec(iAtB))
        vecL(:) = -cellVecs(:, symNeighbourList%iCellVec(iAtB))
        ! get 2D pointer to Sbn overlap block
        ind = symNeighbourList%iPair(iNeighNsort, iAtN) + 1
        nOrbAt = descN(iNOrb)
        nOrbNeigh = descB(iNOrb)
        pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
        ! pre-tabulate g-resolved \gamma_{\mu\beta}(\vec{g}-\vec{l})
        gammaMB(:) = getGammaGResolved(this, iAtM, iAtBfold, iSpM, iSpB, this%coords, this%rCoords,&
            & cellVecsG, rCellVecsG, vecL, rVecL)
        call index_heap_sort(gammaMBsortIdx, gammaMB)
        gammaMBsortIdx(:) = gammaMBsortIdx(size(gammaMBsortIdx):1:-1)
        nNonZeroGammaMB = getNumberOfNonZeroElements(gammaMB(gammaMBsortIdx),&
            & this%pScreeningThreshold)
        loopA: do iNeighM = 0, nNeighbourCamSym(iAtM)
          iNeighMsort = overlapIndices(iAtM)%array(iNeighM + 1) - 1
          maxEstimate = pMaxpSbnMax * testSquareOver(iAtM)%array(iNeighMsort + 1)
          if (maxEstimate < this%pScreeningThreshold) exit loopA
          iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighMsort, iAtM)
          iAtAfold = symNeighbourList%img2CentCell(iAtA)
          iSpA = this%species0(iAtAfold)
          descA = getDescriptor(iAtAfold, iSquare)
          ! get continuous 2D copy of Pab density matrix block
          Pab = deltaDeltaRhoSqr(descA(iStart):descA(iEnd), descB(iStart):descB(iEnd), :,:,:)
          ! get real-space \vec{h} for gamma arguments
          rVecH(:) = rCellVecs(:, symNeighbourList%iCellVec(iAtA))
          vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
          ! pre-tabulate g-resolved \gamma_{\alpha\nu}(\vec{g}+\vec{h})
          gammaAN(:) = getGammaGResolved(this, iAtAfold, iAtN, iSpA, iSpN, this%coords,&
              & this%rCoords, cellVecsG, rCellVecsG, vecH, rVecH)
          call index_heap_sort(gammaANsortIdx, gammaAN)
          gammaANsortIdx(:) = gammaANsortIdx(size(gammaANsortIdx):1:-1)
          nNonZeroGammaAN = getNumberOfNonZeroElements(gammaAN(gammaANsortIdx),&
              & this%pScreeningThreshold)
          ! pre-tabulate g-resolved \gamma_{\alpha\beta}(\vec{g}+\vec{h}-\vec{l})
          gammaAB(:) = getGammaGResolved(this, iAtAfold, iAtBfold, iSpA, iSpB, this%coords,&
              & this%rCoords, cellVecsG, rCellVecsG, vecH + vecL, rVecH + rVecL)
          call index_heap_sort(gammaABsortIdx, gammaAB)
          gammaABsortIdx(:) = gammaABsortIdx(size(gammaABsortIdx):1:-1)
          nNonZeroGammaAB = getNumberOfNonZeroElements(gammaAB(gammaABsortIdx),&
              & this%pScreeningThreshold)
          ! get 2D pointer to Sam overlap block
          ind = symNeighbourList%iPair(iNeighMsort, iAtM) + 1
          nOrbAt = descM(iNOrb)
          nOrbNeigh = descA(iNOrb)
          ! S_{\alpha\mu}(h)
          pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          ! S_{\mu\alpha}(-h)
          pSamT = transpose(pSam(1:nOrbNeigh, 1:nOrbAt))

          tot(1:descM(iNOrb), 1:descN(iNOrb)) = 0.0_dp

          loopGMN: do iG = 1, nNonZeroGammaMN
            iGMN = gammaMNsortIdx(iG)
            bvKIndex(:) = this%foldToBvKIndex(cellVecsG(:, iGMN) + vecH + vecL)
            phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoint, cellVecsG(:, iGMN)))

            pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) = matmul(pSamT,&
                & Pab(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3)))
            Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) = matmul(Pab(:,:, bvKIndex(1), bvKIndex(2),&
                & bvKIndex(3)), pSbn)

            ! term #1
            pSamT_Pab_pSbn(1:descM(iNOrb), 1:descN(iNOrb)) = matmul(pSamT_Pab(1:descM(iNOrb),&
                & 1:descB(iNOrb)), pSbn)

            tot(1:descM(iNOrb), 1:descN(iNOrb)) = tot(1:descM(iNOrb), 1:descN(iNOrb))&
                & + pSamT_Pab_pSbn(1:descM(iNOrb), 1:descN(iNOrb)) * gammaMN(iGMN) * phase
          end do loopGMN

          loopGMB: do iG = 1, nNonZeroGammaMB
            iGMB = gammaMBsortIdx(iG)
            bvKIndex(:) = this%foldToBvKIndex(cellVecsG(:, iGMB) + vecH + vecL)
            phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoint, cellVecsG(:, iGMB)))

            pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) = matmul(pSamT,&
                & Pab(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3)))
            Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) = matmul(Pab(:,:, bvKIndex(1), bvKIndex(2),&
                & bvKIndex(3)), pSbn)

            ! term #2
            tot(1:descM(iNOrb), 1:descN(iNOrb)) = tot(1:descM(iNOrb), 1:descN(iNOrb))&
                & + matmul(pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) * gammaMB(iGMB), pSbn) * phase
          end do loopGMB

          loopGAN: do iG = 1, nNonZeroGammaAN
            iGAN = gammaANsortIdx(iG)
            bvKIndex(:) = this%foldToBvKIndex(cellVecsG(:, iGAN) + vecH + vecL)
            phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoint, cellVecsG(:, iGAN)))

            pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) = matmul(pSamT,&
                & Pab(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3)))
            Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) = matmul(Pab(:,:, bvKIndex(1), bvKIndex(2),&
                & bvKIndex(3)), pSbn)

            ! term #3
            tot(1:descM(iNOrb), 1:descN(iNOrb)) = tot(1:descM(iNOrb), 1:descN(iNOrb))&
                & + matmul(pSamT, Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) * gammaAN(iGAN)) * phase
          end do loopGAN

          loopGAB: do iG = 1, nNonZeroGammaAB
            iGAB = gammaABsortIdx(iG)
            bvKIndex(:) = this%foldToBvKIndex(cellVecsG(:, iGAB) + vecH + vecL)
            phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoint, cellVecsG(:, iGAB)))

            pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) = matmul(pSamT,&
                & Pab(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3)))
            Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) = matmul(Pab(:,:, bvKIndex(1), bvKIndex(2),&
                & bvKIndex(3)), pSbn)

            ! term #4
            pSamT_Pab_gammaAB(1:descM(iNorb), 1:descB(iNorb)) = matmul(pSamT, Pab(:,:,&
                & bvKIndex(1), bvKIndex(2), bvKIndex(3)) * gammaAB(iGAB))
            tot(1:descM(iNOrb), 1:descN(iNOrb)) = tot(1:descM(iNOrb), 1:descN(iNOrb))&
                & + matmul(pSamT_Pab_gammaAB(1:descM(iNOrb), 1:descB(iNOrb)), pSbn) * phase
          end do loopGAB

          tmpHSqr(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd))&
              & = tmpHSqr(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd))&
              & + tot(1:descM(iNOrb), 1:descN(iNOrb))

        end do loopA
      end do loopB
    end do loopMN

    if (this%tSpin .or. this%tREKS) then
      tmpHSqr(:,:) = -0.25_dp * tmpHSqr
    else
      tmpHSqr(:,:) = -0.125_dp * tmpHSqr
    end if

  #:if WITH_MPI
    ! Sum up contributions of current MPI group
    call mpifx_allreduceip(env%mpi%groupComm, tmpHSqr, MPI_SUM)
  #:endif

    ! HSqr(:,:) = HSqr + this%camBeta * tmpHSqr
    ! this%lrEnergy = this%lrEnergy + evaluateEnergy_cplx_kptrho(tmpHSqr, kWeight, deltaRhoOutSqrCplx)

    this%hprevCplxHS(:,:, iKS) = this%hprevCplxHS(:,:, iKS) + tmpHSqr
    HSqr(:,:) = HSqr + this%camBeta * this%hprevCplxHS(:,:, iKS)

    ! Add energy contribution but divide by the number of processes working on iKS
  #:if WITH_MPI
    this%lrEnergy = this%lrEnergy + evaluateEnergy_cplx_kptrho(this%hprevCplxHS(:,:, iKS), kWeight,&
        & deltaRhoOutSqrCplx) / real(env%mpi%groupComm%size, dp)
  #:else
    this%lrEnergy = this%lrEnergy + evaluateEnergy_cplx_kptrho(this%hprevCplxHS(:,:, iKS), kWeight,&
        & deltaRhoOutSqrCplx)
  #:endif

  contains

    !> Returns the number of non-zero elements in a descending array of reals.
    function getNumberOfNonZeroElements(array, threshold) result(nNonZeroEntries)

      !> Descending one-dimensional, real-valued array to search
      real(dp), intent(in) :: array(:)

      !> Screening threshold value
      real(dp), intent(in) :: threshold

      !> Number of non-zero entries
      real(dp) :: nNonZeroEntries

      !! iterates over all array elements
      integer :: ii

      nNonZeroEntries = 0

      do ii = 1, size(array)
        if (array(ii) < 1e-16_dp) return
        ! if (array(ii) < threshold) return
        nNonZeroEntries = ii
      end do

    end function getNumberOfNonZeroElements

  end subroutine addLrHamiltonianNeighbour_kpts


  !> Add the CAM-energy contribution to the total energy.
  subroutine addCamEnergy(this, env, energy)

    !> Class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Total energy
    real(dp), intent(inout) :: energy

    call addLrEnergy(this, env, energy)
    call addHfEnergy(this, env, energy)

  end subroutine addCamEnergy


  !> Add the LR-energy contribution to the total energy.
  subroutine addLrEnergy(this, env, energy)

    !> Instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Total energy
    real(dp), intent(inout) :: energy

  #:if WITH_MPI
    call mpifx_allreduceip(env%mpi%globalComm, this%lrEnergy, MPI_SUM)
  #:endif

    energy = energy + this%camBeta * this%lrEnergy

    ! hack for spin unrestricted calculation
    this%lrEnergy = 0.0_dp

  end subroutine addLrEnergy


  !> Finds location of relevant atomic block indices in a dense matrix.
  pure function getDescriptor(iAt, iSquare) result(desc)

    !> relevant atom
    integer, intent(in) :: iAt

    !> indexing array for start of atom orbitals
    integer, intent(in) :: iSquare(:)

    !> resulting location ranges
    integer :: desc(3)

    desc(:) = [iSquare(iAt), iSquare(iAt + 1) - 1, iSquare(iAt + 1) - iSquare(iAt)]

  end function getDescriptor


  !> Evaluates energy from triangles of the Hamiltonian and density matrix (real version).
  pure function evaluateEnergy_real(hamiltonian, densityMat) result(energy)

    !> Hamiltonian matrix
    real(dp), intent(in) :: hamiltonian(:,:)

    !> Density matrix
    real(dp), intent(in) :: densityMat(:,:)

    !> Resulting energy due to CAM contribution
    real(dp) :: energy

    integer :: mu

    energy = 0.0_dp

    do mu = 1, size(hamiltonian, dim=2)
      energy = energy + hamiltonian(mu, mu) * densityMat(mu, mu)&
          & + 2.0_dp * sum(hamiltonian(mu + 1 :, mu) * densityMat(mu + 1 :, mu))
    end do

    energy = 0.5_dp * energy

  end function evaluateEnergy_real


  !> Evaluates energy from triangles of the Hamiltonian and density matrix (complex version).
  pure function evaluateEnergy_cplx(hamiltonian, densityMat, bvKShifts, coeffsDiag, kPoint)&
      & result(energy)

    !> Hamiltonian matrix
    complex(dp), intent(in) :: hamiltonian(:,:)

    !> Density matrix
    real(dp), intent(in), pointer :: densityMat(:,:,:,:,:)

    !> K-point compatible BvK real-space shifts in relative coordinates
    real(dp), intent(in) :: bvKShifts(:,:)

    !> Supercell folding coefficients (diagonal elements)
    integer, intent(in) :: coeffsDiag(:)

    !> K-point in relative coordinates
    real(dp), intent(in) :: kPoint(:)

    !> Resulting energy due to CAM contribution
    real(dp) :: energy

    !! BvK summed k-space density matrix
    complex(dp) :: recDensityMat(size(hamiltonian, dim=1), size(hamiltonian, dim=1))

    !! \mu orbital index
    integer :: mu

    !! Integer BvK real-space shift translated to density matrix indices
    integer :: bvKIndex(3)

    !! Phase factor
    complex(dp) :: phase

    !! Iterates over all BvK real-space vectors
    integer :: iG

    recDensityMat(:,:) = 0.0_dp

    do iG = 1, size(bvKShifts, dim=2)
      phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoint, bvKShifts(:, iG)))
      bvKIndex(:) = nint(bvKShifts(:, iG)) + 1
      recDensityMat(:,:) = recDensityMat&
          & + densityMat(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3)) * phase
    end do

    energy = 0.5_dp * real(sum(hamiltonian * recDensityMat), dp)

  end function evaluateEnergy_cplx


  !> Evaluates energy from triangles of the Hamiltonian and density matrix (complex version).
  pure function evaluateEnergy_cplx_kptrho(hamiltonian, kWeight, densityMat) result(energy)

    !> Hamiltonian matrix
    complex(dp), intent(in) :: hamiltonian(:,:)

    !> K-point weight (for energy contribution)
    real(dp), intent(in) :: kWeight

    !> Density matrix in k-space
    complex(dp), intent(in), pointer :: densityMat(:,:)

    !> Resulting energy due to CAM contribution
    real(dp) :: energy

    energy = 0.5_dp * real(sum(hamiltonian * densityMat), dp) * kWeight

  end function evaluateEnergy_cplx_kptrho


  !> Returns the value of a polynomial of 5th degree at x (or its derivative).
  !! The polynomial is created with the following boundary conditions:
  !! Its value, its 1st and 2nd derivatives are zero at x = rCut and agree with the provided values
  !! at x = rDamp, i.e. x = rCut - delta.
  !! WARNING: To avoid additional branches, there are no consistency checks, e.g. that rDamp < rCut
  pure function poly5zero(y0, y0p, y0pp, xx, rDamp, rCut, tDerivative) result(yy)

    !> Value of the polynom at x = rDamp
    real(dp), intent(in) :: y0

    !> Value of the 1st derivative at x = rDamp
    real(dp), intent(in) :: y0p

    !> Value of the 2nd derivative at x = rDamp
    real(dp), intent(in) :: y0pp

    !> Point where the polynomial should be calculated
    real(dp), intent(in) :: xx

    !> Point, where the polynomials value and first two derivatives should take the provided values
    real(dp), intent(in) :: rDamp

    !> Point, where the polynomial (and its 1st/2nd derivative) is supposed to be zero
    real(dp), intent(in) :: rCut

    !> True, if the derivative at xx is desired, otherwise the function value it returned
    logical, intent(in), optional :: tDerivative

    !! True, if the derivative at xx is desired, otherwise the function value it returned
    !! Default: .false.
    logical :: tPrime

    !! Value of the polynomial at xx (in general should satisfy rDamp <= x <= rCut)
    real(dp) :: yy

    !! Polynomial coefficients
    real(dp) :: aa, bb, cc, dd, ee, ff

    !! rCut - rDamp
    real(dp) :: delta

    if (present(tDerivative)) then
      tPrime = tDerivative
    else
      tPrime = .false.
    end if

    ! 5th order polynomial definition
    ! f(x) = ax^5 + bx^4 + cx^3 + dx^2 + ex + f
    ! f'(x) = 5ax^4 + 4bx^3 + 3cx^2 + 2dx + e
    ! f''(x) = 20ax^3 + 12bx^2 + 6cx + 2d

    ! boundary conditions:
    ! f(rCut) = 0, f'(rCut) = 0, f''(rCut) = 0
    ! f(rDamp) = y0, f'(rDamp) = y0p, f''(rDamp) = y0pp

    delta = rCut - rDamp

    aa = -(delta**2 * y0pp + 6.0_dp * delta * y0p + 12.0_dp * y0) / (2.0_dp * delta**5)

    bb = (rCut * (5.0_dp * delta**2 * y0pp + 30.0_dp * delta * y0p + 60.0_dp * y0)&
        & - 2.0_dp * delta**3 * y0pp - 14.0_dp * delta**2 * y0p - 30.0_dp * delta * y0)&
        & / (2.0_dp * delta**5)

    cc = -(rCut * (-8.0_dp * delta**3 * y0pp - 56.0_dp * delta**2 * y0p - 120.0_dp * delta * y0)&
        & + rCut**2 * (10.0_dp * delta**2 * y0pp + 60.0_dp * delta * y0p + 120.0_dp * y0)&
        & + delta**4 * y0pp + 8.0_dp * delta**3 * y0p + 20.0_dp * delta**2 * y0)&
        & / (2.0_dp * delta**5)

    dd = (rCut * (3.0_dp * delta**4 * y0pp + 24.0_dp * delta**3 * y0p + 60.0_dp * delta**2 * y0)&
        & + rCut**2 * (-12.0_dp * delta**3 * y0pp - 84.0_dp * delta**2 * y0p - 180.0_dp * delta&
        & * y0) + rCut**3 * (10.0_dp * delta**2 * y0pp + 60.0_dp * delta * y0p + 120.0_dp * y0))&
        & / (2.0_dp * delta**5)

    ee = -(rCut**2 * (3.0_dp * delta**4 * y0pp + 24.0_dp * delta**3 * y0p + 60.0_dp * delta**2&
        & * y0) + rCut**3*(-8.0_dp * delta**3 * y0pp - 56.0_dp * delta**2 * y0p - 120.0_dp * delta&
        & * y0) + rCut**4 * (5.0_dp * delta**2 * y0pp + 30.0_dp * delta * y0p + 60.0_dp * y0))&
        & / (2.0_dp * delta**5)

    ff = (rCut**3 * (delta**4 * y0pp + 8.0_dp * delta**3 * y0p + 20.0_dp * delta**2 * y0)&
        & + rCut**4 * (-2.0_dp * delta**3 * y0pp - 14.0_dp * delta**2 * y0p - 30.0_dp * delta * y0)&
        & + rCut**5 * (delta**2 * y0pp + 6.0_dp * delta * y0p + 12.0_dp * y0)) / (2.0_dp * delta**5)

    if (tPrime) then
      yy = 5.0_dp * aa * xx**4 + 4.0_dp * bb * xx**3 + 3.0_dp * cc * xx**2 + 2.0_dp * dd * xx + ee
    else
      yy = aa * xx**5 + bb * xx**4 + cc * xx**3 + dd * xx**2 + ee * xx + ff
    end if

  end function poly5zero


  !> Returns the numerical second derivative of gamma, by means of a central finite difference.
  function getddNumericalLrGammaDeriv(this, iSp1, iSp2, dist, delta) result(ddGamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Delta for finite differences
    real(dp), intent(in) :: delta

    !> Numerical gamma derivative
    real(dp) :: ddGamma

    ddGamma = (getdAnalyticalLrGammaDeriv(this, iSp1, iSp2, dist + delta)&
        & - getdAnalyticalLrGammaDeriv(this, iSp1, iSp2, dist - delta)) / (2.0_dp * delta)

  end function getddNumericalLrGammaDeriv


  !> Returns the numerical second derivative of gamma, by means of a central finite difference.
  function getddNumericalHfGammaDeriv(this, iSp1, iSp2, dist, delta) result(ddGamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Delta for finite differences
    real(dp), intent(in) :: delta

    !> Numerical gamma derivative
    real(dp) :: ddGamma

    ddGamma = (getdAnalyticalHfGammaDeriv(this, iSp1, iSp2, dist + delta)&
        & - getdAnalyticalHfGammaDeriv(this, iSp1, iSp2, dist - delta)) / (2.0_dp * delta)

  end function getddNumericalHfGammaDeriv


  !> Calculates analytical, screened Coulomb, long-range gamma.
  function getLrScreenedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    gamma = getAnalyticalLrGammaValue(this, iSp1, iSp2, this%omega + this%auxiliaryScreening, dist)&
        & - getAnalyticalLrGammaValue(this, iSp1, iSp2, this%auxiliaryScreening, dist)

  end function getLrScreenedGammaValue


  !> Calculates analytical, truncated Coulomb, long-range gamma.
  function getLrTruncatedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    if (dist > this%gammaDamping .and. dist < this%gammaCutoff) then
      gamma = poly5zero(this%lrGammaAtDamping(iSp1, iSp2), this%lrdGammaAtDamping(iSp1, iSp2),&
          & this%lrddGammaAtDamping(iSp1, iSp2), dist, this%gammaDamping, this%gammaCutoff,&
          & tDerivative=.false.)
    elseif (dist >= this%gammaCutoff) then
      gamma = 0.0_dp
    else
      gamma = getAnalyticalLrGammaValue(this, iSp1, iSp2, this%omega, dist)
    end if

    ! if (dist >= this%gammaCutoff) then
    !   gamma = 0.0_dp
    ! else
    !   gamma = getAnalyticalLrGammaValue(this, iSp1, iSp2, this%omega, dist)
    ! end if

  end function getLrTruncatedGammaValue


  !> Calculates analytical, truncated Coulomb, full-range Hartree-Fock gamma.
  function getHfTruncatedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    if (dist > this%gammaDamping .and. dist < this%gammaCutoff) then
      gamma = poly5zero(this%hfGammaAtDamping(iSp1, iSp2), this%hfdGammaAtDamping(iSp1, iSp2),&
          & this%hfddGammaAtDamping(iSp1, iSp2), dist, this%gammaDamping, this%gammaCutoff,&
          & tDerivative=.false.)
    elseif (dist >= this%gammaCutoff) then
      gamma = 0.0_dp
    else
      gamma = getAnalyticalHfGammaValue(this, iSp1, iSp2, dist)
    end if

    ! if (dist >= this%gammaCutoff) then
    !   gamma = 0.0_dp
    ! else
    !   gamma = getAnalyticalHfGammaValue(this, iSp1, iSp2, dist)
    ! end if

  end function getHfTruncatedGammaValue


  function getLrGammaGSum(this, iAt1, iAt2, iSp1, iSp2, rCellVecsG) result(gamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> Resulting Gamma, summed up for g-vectors
    real(dp) :: gamma

    !! Distance between the two atoms
    real(dp) :: dist

    !! Index of real-space \vec{g} summation
    integer :: iG

    gamma = 0.0_dp

    loopG: do iG = 1, size(rCellVecsG, dim=2)
      dist = norm2(this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG)))
      gamma = gamma + getLrTruncatedGammaValue(this, iSp1, iSp2, dist)
      ! gamma = gamma + getLrScreenedGammaValue(this, iSp1, iSp2, dist)
    end do loopG

  end function getLrGammaGSum


  function getHfGammaGSum(this, iAt1, iAt2, iSp1, iSp2, rCellVecsG) result(gamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> Resulting Gamma, summed up for g-vectors
    real(dp) :: gamma

    !! Distance between the two atoms
    real(dp) :: dist

    !! Index of real-space \vec{g} summation
    integer :: iG

    gamma = 0.0_dp

    loopG: do iG = 1, size(rCellVecsG, dim=2)
      dist = norm2(this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG)))
      gamma = gamma + getHfTruncatedGammaValue(this, iSp1, iSp2, dist)
      ! gamma = gamma + getHfScreenedGammaValue(this, iSp1, iSp2, dist)
    end do loopG

  end function getHfGammaGSum


  function getGammaGResolved(this, iAt1, iAt2, iSp1, iSp2, coords, rCoords, cellVecsG, rCellVecsG,&
    & shift, rShift) result(gammas)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Relative coordinates, including periodic images
    real(dp), intent(in) :: coords(:,:)

    !> Absolute coordinates, including periodic images
    real(dp), intent(in) :: rCoords(:,:)

    !> Vectors to unit cells in relative units
    real(dp), intent(in) :: cellVecsG(:,:)

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> Relative shift of atom position (Gamma arguments)
    real(dp), intent(in) :: shift(:)

    !> Absolute shift of atom position (Gamma arguments)
    real(dp), intent(in) :: rShift(:)

    !> Resulting Gammas, g-vector resolved
    real(dp) :: gammas(size(rCellVecsG, dim=2))

    !! Total shift of atom position (Gamma arguments) in relative coordinates
    real(dp) :: totshift(3)

    !! Total shift of atom position (Gamma arguments) in absolute coordinates
    real(dp) :: rTotshift(3)

    !! Distance between the two atoms
    real(dp) :: dist

    !> Distance vector (relative coordinates) between atoms
    real(dp) :: distVec(3)

    !! Index of real-space \vec{g} summation
    integer :: iG

    loopG: do iG = 1, size(rCellVecsG, dim=2)
      totshift(:) = cellVecsG(:, iG) + shift
      rTotshift(:) = rCellVecsG(:, iG) + rShift
      dist = norm2(rCoords(:, iAt1) - (rCoords(:, iAt2) + rTotshift))
      distVec(:) = coords(:, iAt1) - (coords(:, iAt2) + totshift)
      gammas(iG) = getLrTruncatedGammaValue(this, iSp1, iSp2, dist)
      ! gammas(iG) = getLrScreenedGammaValue(this, iSp1, iSp2, dist)
    end do loopG

  end function getGammaGResolved


  !> Calculates analytical, truncated Coulomb, long-range gamma derivative.
  function getTruncatedLrGammaPrimeValue(this, iSp1, iSp2, dist) result(dgamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative (1st)
    real(dp) :: dgamma

    if (dist > this%gammaDamping .and. dist < this%gammaCutoff) then
      dgamma = poly5zero(this%lrGammaAtDamping(iSp1, iSp2), this%lrdGammaAtDamping(iSp1, iSp2),&
          & this%lrddGammaAtDamping(iSp1, iSp2), dist, this%gammaDamping, this%gammaCutoff,&
          & tDerivative=.true.)
    elseif (dist >= this%gammaCutoff) then
      dgamma = 0.0_dp
    else
      dgamma = getdAnalyticalLrGammaDeriv(this, iSp1, iSp2, dist)
    end if

  end function getTruncatedLrGammaPrimeValue


  !> Calculates analytical, truncated Coulomb, long-range gamma derivative.
  function getTruncatedHfGammaPrimeValue(this, iSp1, iSp2, dist) result(dgamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative (1st)
    real(dp) :: dgamma

    if (dist > this%gammaDamping .and. dist < this%gammaCutoff) then
      dgamma = poly5zero(this%hfGammaAtDamping(iSp1, iSp2), this%hfdGammaAtDamping(iSp1, iSp2),&
          & this%hfddGammaAtDamping(iSp1, iSp2), dist, this%gammaDamping, this%gammaCutoff,&
          & tDerivative=.true.)
    elseif (dist >= this%gammaCutoff) then
      dgamma = 0.0_dp
    else
      dgamma = getdAnalyticalHfGammaDeriv(this, iSp1, iSp2, dist)
    end if

  end function getTruncatedHfGammaPrimeValue


  function getLrGammaPrimeGSum(this, iAt1, iAt2, iSp1, iSp2, rCellVecsG) result(dgamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> Resulting Gamma derivative (1st), summed up for g-vectors
    real(dp) :: dgamma(3)

    !! Temporary distance vector
    real(dp) :: distVect(3)

    !! Distance between the two atoms
    real(dp) :: dist

    !! Index of real-space \vec{g} summation
    integer :: iG

    dgamma(:) = 0.0_dp

    loopG: do iG = 1, size(rCellVecsG, dim=2)
      distVect(:) = this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG))
      dist = norm2(distVect)
      distVect(:) = distVect / dist
      dgamma(:) = dgamma + distVect * getTruncatedLrGammaPrimeValue(this, iSp1, iSp2, dist)
    end do loopG

  end function getLrGammaPrimeGSum


  function getHfGammaPrimeGSum(this, iAt1, iAt2, iSp1, iSp2, rCellVecsG) result(dgamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> Resulting Gamma derivative (1st), summed up for g-vectors
    real(dp) :: dgamma(3)

    !! Temporary distance vector
    real(dp) :: distVect(3)

    !! Distance between the two atoms
    real(dp) :: dist

    !! Index of real-space \vec{g} summation
    integer :: iG

    dgamma(:) = 0.0_dp

    loopG: do iG = 1, size(rCellVecsG, dim=2)
      distVect(:) = this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG))
      dist = norm2(distVect)
      distVect(:) = distVect / dist
      dgamma(:) = dgamma + distVect * getTruncatedHfGammaPrimeValue(this, iSp1, iSp2, dist)
    end do loopG

  end function getHfGammaPrimeGSum


  !> Calculates analytical long-range gamma.
  function getAnalyticalLrGammaValue(this, Sp1, Sp2, omega, dist) result(gamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: Sp1

    !> Second species
    integer, intent(in) :: Sp2

    !> Screening parameter
    real(dp), intent(in) :: omega

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma
    real(dp) :: gamma

    real(dp) :: tauA, tauB
    real(dp) :: prefac, tmp, tmp2, tau

    tauA = 3.2_dp * this%hubbu(Sp1)
    tauB = 3.2_dp * this%hubbu(Sp2)

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff) then
        tau = 0.5_dp * (tauA + tauB)
        tmp = 5.0_dp * tau**6 + 15.0_dp * tau**4 * omega**2 - 5.0_dp * tau**2 * omega**4 + omega**6
        tmp = tmp * 0.0625_dp / tau**5 - omega
        tmp = tmp * tau**8 / (tau**2 - omega**2)**4
        gamma = tau * 0.3125_dp - tmp
      else
        call error("Error(RangeSep): R = 0, Ua != Ub")
      end if
    else
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff ) then
        tauA = 0.5_dp * (tauA + tauB)
        tmp2 = ((dist * tauA)**3 / 48.0_dp + 0.1875_dp * (dist * tauA)**2 +&
            & 0.6875_dp * (dist * tauA) + 1.0_dp) * exp(-tauA * dist) / dist
        tmp = -tauA**8 / (tauA**2 - omega**2)**4 * (tmp2 + exp(-tauA*dist) * &
            & (dist**2 * (3.0_dp * tauA**4 * omega**4 - 3.0_dp * tauA**6 * omega**2 - &
            & tauA**2 * omega**6) + dist * (15.0_dp * tauA**3 * omega**4 - &
            & 21.0_dp * tauA**5 * omega**2 - 3.0_dp * tauA * omega**6) + &
            & (15.0_dp * tauA**2 * omega**4 - 45.0_dp * tauA**4 * omega**2 - &
            & 3.0_dp * omega**6)) / (48.0_dp * tauA**5))
        gamma = 1.0_dp/dist - tmp2 - (tauA**8 / (tauA**2 - omega**2)**4 *&
            & exp(-omega * dist) / dist + tmp)
      else
        ! off-site, Ua != Ub
        prefac = tauA**4 / (tauA * tauA - omega * omega)**2
        prefac = prefac * tauB**4 / (tauB * tauB - omega * omega)**2
        prefac = prefac * exp(-omega * dist) / dist
        tmp = prefac&
            & - getYGammaSubPart(tauA, tauB, dist, omega)&
            & - getYGammaSubPart(tauB, tauA, dist, omega)
        tmp = 1.0_dp / dist - tmp
        tmp = tmp&
            & - getYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getYGammaSubPart(tauB, tauA, dist, 0.0_dp)
        gamma = tmp
      end if
    end if

  end function getAnalyticalLrGammaValue


  !> Returns the subexpression for the evaluation of the off-site Y-Gamma-integral.
  pure function getYGammaSubPart(tauA, tauB, R, omega) result(yGamma)

    !> decay constant site A
    real(dp), intent(in) :: tauA

    !> decay constant site B
    real(dp), intent(in) :: tauB

    !> separation of the sites A and B
    real(dp), intent(in) :: R

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> resulting off-site Y-Gamma-integral
    real(dp) :: yGamma

    !!
    real(dp) :: prefac, tmp

    tmp = (tauA - omega)
    tmp = tmp * (tauA + omega)
    prefac = tauA * tauA / tmp
    tmp = (tauB**6 - 3.0_dp * tauA * tauA * tauB**4 + 2.0_dp * omega * omega * tauB**4) / R
    tmp = tmp * prefac * prefac / (tauA * tauA - tauB * tauB)**3
    tmp = tauA * tauB**4 * 0.5_dp * prefac / (tauB * tauB - tauA * tauA)**2 - tmp
    yGamma = tmp * exp(-tauA * R)

  end function getYGammaSubPart


  !> Derivative of analytical long-range Gamma
  function getdAnalyticalLrGammaDeriv(this, iSp1, iSp2, dist) result(dAnalyticalGammaDeriv)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Species index of first and second atom
    integer, intent(in) :: iSp1, iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting d gamma / d dist
    real(dp) :: dAnalyticalGammaDeriv

    real(dp) :: tauA, tauB, omega
    real(dp) :: prefac, tmp, tmp2, dTmp, dTmp2

    tauA = 3.2_dp * this%hubbu(iSp1)
    tauB = 3.2_dp * this%hubbu(iSp2)
    omega = this%omega

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff) then
        dAnalyticalGammaDeriv = 0.0_dp
      else
        call error("Error(RangeSep1): R = 0, Ua != Ub")
      end if
    else
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff ) then
        tauA = 0.5_dp * (tauA + tauB)

        tmp = dist**2 * (3.0_dp*tauA**4*omega**4 - 3.0_dp * tauA**6 * omega**2 - tauA**2*omega**6)&
            & + dist * (15.0_dp*tauA**3*omega**4 - 21.0_dp*tauA**5*omega**2 - 3.0_dp*tauA*omega**6)&
            & + (15.0_dp * tauA**2 * omega**4 - 45.0_dp * tauA**4 * omega**2 - 3.0_dp * omega**6)

        dTmp = 2.0_dp*dist*(3.0_dp*tauA**4*omega**4 - 3.0_dp*tauA**6*omega**2 - tauA**2*omega**6)&
            & + (15.0_dp*tauA**3*omega**4 - 21.0_dp*tauA**5*omega**2 - 3.0_dp*tauA*omega**6)

        dtmp = (dtmp*exp(-tauA*dist) -tmp*tauA*exp(-tauA*dist))/ (48.0_dp * tauA**5)

        tmp2 = ( dist**2 * tauA**3 / 48.0_dp + 0.1875_dp * dist * tauA**2 + 0.6875_dp * tauA&
            & + 1.0_dp / dist ) * exp(-tauA * dist)

        dTmp2 = &
            & (2.0_dp*dist*tauA**3/48.0_dp + 0.1875_dp*tauA**2 -1.0_dp/dist**2) * exp(-tauA*dist)&
            & -(dist**2*tauA**3/48.0_dp + 0.1875_dp*dist*tauA**2 + 0.6875_dp*tauA +1.0_dp/dist)&
            & * tauA * exp(-tauA * dist)

        dAnalyticalGammaDeriv = -1.0_dp/dist**2 -dtmp2&
            & + (tauA**8 / (tauA**2 - omega**2)**4) * (dtmp + dtmp2 + omega*exp(-omega * dist)/dist&
            & +exp(-omega * dist) / dist**2)

      else
        ! off-site, Ua != Ub
        prefac = tauA**4 / (tauA * tauA - omega * omega )**2
        prefac = prefac * tauB**4 / (tauB * tauB - omega * omega )**2
        prefac = prefac * (-omega * exp(-omega * dist) / dist - exp(-omega * dist) / dist**2)
        dAnalyticalGammaDeriv = -1.0_dp / (dist**2) - prefac&
            & + getdYGammaSubPart(tauA, tauB, dist, omega)&
            & + getdYGammaSubPart(tauB, tauA, dist, omega)&
            & - getdYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getdYGammaSubPart(tauB, tauA, dist, 0.0_dp)
      end if
    end if

  end function getdAnalyticalLrGammaDeriv


  !> Returns the derivative of the subexpression for the evaluation of the off-site
  !! Y-Gamma-integral. Note that tauA /= tauB.
  pure function getdYGammaSubPart(tauA, tauB, R, omega) result(dYGammaSubPart)

    !> decay constant site A
    real(dp), intent(in) :: tauA

    !> decay constant site B
    real(dp), intent(in) :: tauB

    !> separation of the sites A and B
    real(dp), intent(in) :: R

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> resulting derivative of the subexpression
    real(dp) :: dYGammaSubPart

    !! auxiliary variables
    real(dp) :: prefac, tmp, tmp2, dtmp

    tmp = tauA**2 - omega**2
    prefac = tauA * tauA / tmp
    tmp = prefac * prefac / (tauA * tauA - tauB * tauB)**3
    dtmp = tmp * (tauB**6 - 3.0_dp * tauA * tauA * tauB**4 + 2.0_dp * omega * omega * tauB**4)/R**2
    tmp = tmp * (tauB**6 - 3.0_dp * tauA * tauA * tauB**4 + 2.0_dp * omega * omega * tauB**4) / R
    tmp2 = tauA * tauB**4 * 0.5_dp * prefac / (tauB * tauB - tauA * tauA )**2 - tmp

    dYGammaSubPart = (dtmp - tmp2 * tauA) * exp(-tauA * R)

  end function getdYGammaSubPart


  !> Returns the derivative of long-range gamma for iAtom1, iAtom2 (non-periodic systems only).
  subroutine getLrGammaPrimeValue_cluster(this, grad, iAtom1, iAtom2, coords0, species0)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Gradient of gamma between atoms
    real(dp), intent(out) :: grad(3)

    !> First atom
    integer, intent(in) :: iAtom1

    !> Second atom
    integer, intent(in) :: iAtom2

    !> Coordinates of atoms in central cell
    real(dp), intent(in) :: coords0(:,:)

    !> List of all atomic species in central cell
    integer, intent(in) :: species0(:)

    !! Species index of first and second atom
    integer :: iSp1, iSp2

    !! Distance(-vector) of the two atoms
    real(dp) :: vect(3), dist

    iSp1 = species0(iAtom1)
    iSp2 = species0(iAtom2)

    ! analytical derivatives
    vect(:) = coords0(:, iAtom1) - coords0(:, iAtom2)
    dist = sqrt(sum(vect**2))
    vect(:) = vect / dist
    grad(:) = vect * getdAnalyticalLrGammaDeriv(this, iSp1, iSp2, dist)

  end subroutine getLrGammaPrimeValue_cluster


  !> Returns the derivative of long-range gamma for iAtom1, iAtom2 (periodic systems).
  subroutine getLrGammaPrimeValue_periodic(this, grad, iAtom1, iAtom2, coords, species0,&
      & img2CentCell)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Gradient of gamma between atoms
    real(dp), intent(out) :: grad(3)

    !> First atom
    integer, intent(in) :: iAtom1

    !> Second atom
    integer, intent(in) :: iAtom2

    !> Coordinates of atoms, including periodic images
    real(dp), intent(in) :: coords(:,:)

    !> List of all atomic species in central cell
    integer, intent(in) :: species0(:)

    !> Map images of atoms to the central cell
    integer, intent(in) :: img2CentCell(:)

    !! Species index of first and second atom
    integer :: iSp1, iSp2

    !! Distance(-vector) of the two atoms
    real(dp) :: vect(3), dist

    iSp1 = species0(img2CentCell(iAtom1))
    iSp2 = species0(img2CentCell(iAtom2))

    ! analytical derivatives
    vect(:) = coords(:, iAtom1) - coords(:, iAtom2)
    dist = sqrt(sum(vect**2))
    vect(:) = vect / dist
    grad(:) = vect * getdAnalyticalLrGammaDeriv(this, iSp1, iSp2, dist)

  end subroutine getLrGammaPrimeValue_periodic


  !> Adds CAM gradients due to full-/long-range HF-contributions (non-periodic version).
  subroutine addCamGradients_cluster(this, gradients, derivator, deltaRho, skOverCont, coords0,&
      & species0, orb, iSquare, ovrlapMat, iNeighbour, nNeighbourSK)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> Density matrix difference from reference q0
    real(dp), intent(in) :: deltaRho(:,:,:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> Atomic coordinates in central cell
    real(dp), intent(in) :: coords0(:,:)

    !> Species of atoms in central cell
    integer, intent(in) :: species0(:)

    !> Orbital information for system
    type(TOrbitals), intent(in) :: orb

    !> Index for dense arrays
    integer, intent(in) :: iSquare(:)

    !> Overlap matrix
    real(dp), intent(in) :: ovrlapMat(:,:)

    !> Neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of atoms neighbouring each site where the overlap is non-zero
    integer, intent(in) :: nNeighbourSK(:)

    ! Add long-range contribution if needed.
    ! For pure Hyb, camBeta would be zero anyway, but we want to save as much time as possible.
    if (this%tLc .or. this%tCam) then
      call addLrGradients_cluster(this, gradients, derivator, deltaRho, skOverCont, coords0,&
          & species0, orb, iSquare, ovrlapMat, iNeighbour, nNeighbourSK)
    end if

    ! Add full-range Hartree-Fock contribution if needed.
    ! For pure LC, camAlpha would be zero anyway, but we want to save as much time as possible.
    if (this%tHyb .or. this%tCam) then
      call addHfGradients_cluster(this, gradients, derivator, deltaRho, skOverCont,&
          & coords0, species0, orb, iSquare, ovrlapMat, iNeighbour, nNeighbourSK)
    end if

  end subroutine addCamGradients_cluster


  !> Adds CAM gradients due to full-/long-range HF-contributions (Gamma-point version).
  subroutine addCamGradients_gamma(this, deltaRhoSqr, skOverCont, symNeighbourList,&
      & nNeighbourCamSym, iSquare, orb, derivator, gradients)

    !> Class instance
    class(TRangeSepFunc), intent(in), target :: this

    !> Square (unpacked) delta density matrix
    real(dp), intent(in) :: deltaRhoSqr(:,:,:)

    !> SK overlap container
    type(TSlakoCont), intent(in) :: skOverCont

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    ! Add long-range contribution if needed.
    ! For pure Hyb, camBeta would be zero anyway, but we want to save as much time as possible.
    if (this%tLc .or. this%tCam) then
      call addLrGradients_gamma(this, deltaRhoSqr, skOverCont, symNeighbourList, nNeighbourCamSym,&
          & iSquare, orb, derivator, gradients)
    end if

    ! Add full-range Hartree-Fock contribution if needed.
    ! For pure LC, camAlpha would be zero anyway, but we want to save as much time as possible.
    ! if (this%tHyb .or. this%tCam) then
    !   call addHfGradients_cluster()
    ! end if

  end subroutine addCamGradients_gamma


  !> Adds gradients due to long-range HF-contribution (non-periodic version).
  subroutine addLrGradients_cluster(this, gradients, derivator, deltaRho, skOverCont, coords0,&
      & species0, orb, iSquare, ovrlapMat, iNeighbour, nNeighbourSK)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> Density matrix difference from reference q0
    real(dp), intent(in) :: deltaRho(:,:,:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> Atomic coordinates in central cell
    real(dp), intent(in) :: coords0(:,:)

    !> Species of atoms in central cell
    integer, intent(in) :: species0(:)

    !> Orbital information for system
    type(TOrbitals), intent(in) :: orb

    !> Index for dense arrays
    integer, intent(in) :: iSquare(:)

    !> Overlap matrix
    real(dp), intent(in) :: ovrlapMat(:,:)

    !> Neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of atoms neighbouring each site where the overlap is non-zero
    integer, intent(in) :: nNeighbourSK(:)

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    integer :: nAtom0, iAtK, iNeighK, iAtB, iNeighB, iAtC, iAtA, kpa
    real(dp) :: tmpgamma1, tmpgamma2
    real(dp) :: tmpforce(3), tmpforce_r(3), tmpforce2, tmpmultvar1
    integer :: nSpin, iSpin, mu, alpha, beta, ccc, kkk
    real(dp) :: sPrimeTmp(orb%mOrb, orb%mOrb, 3)
    real(dp) :: sPrimeTmp2(orb%mOrb, orb%mOrb, 3)
    real(dp), allocatable :: gammaPrimeTmp(:,:,:), tmpOvr(:,:), tmpRho(:,:,:), tmpderiv(:,:)

    nSpin = size(deltaRho, dim=3)
    call allocateAndInit(tmpOvr, tmpRho, gammaPrimeTmp, tmpderiv)
    nAtom0 = size(this%species0)
    tmpderiv(:,:) = 0.0_dp
    ! sum K
    loopK: do iAtK = 1, nAtom0
      ! C >= K
      loopC: do iNeighK = 0, nNeighbourSK(iAtK)
        iAtC = iNeighbour(iNeighK, iAtK)
        ! evaluate the ovr_prime
        sPrimeTmp2(:,:,:) = 0.0_dp
        sPrimeTmp(:,:,:) = 0.0_dp
        if (iAtK /= iAtC) then
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords0, species0, iAtK, iAtC, orb)
          call derivator%getFirstDeriv(sPrimeTmp2, skOverCont, coords0, species0, iAtC, iAtK, orb)
        end if
        loopB: do iAtB = 1, nAtom0
          ! A > B
          loopA: do iNeighB = 0, nNeighbourSK(iAtB)
            iAtA = iNeighbour(iNeighB, iAtB)
            tmpgamma1 = this%lrGammaEval0(iAtK, iAtB) + this%lrGammaEval0(iAtC, iAtB)
            tmpgamma2 = tmpgamma1 + this%lrGammaEval0(iAtK, iAtA) + this%lrGammaEval0(iAtC, iAtA)
            tmpforce(:) = 0.0_dp
            tmpforce_r(:) = 0.0_dp
            tmpforce2 = 0.0_dp
            ccc = 0
            do mu = iSquare(iAtC), iSquare(iAtC + 1) - 1
              ccc = ccc + 1
              kkk = 0
              do kpa = iSquare(iAtK), iSquare(iAtK + 1) - 1
                kkk = kkk + 1
                tmpmultvar1 = 0.0_dp
                do iSpin = 1, nSpin
                  do alpha = iSquare(iAtA), iSquare(iAtA + 1) - 1
                    do beta = iSquare(iAtB), iSquare(iAtB + 1) - 1
                      tmpmultvar1 = tmpmultvar1 + tmpOvr(beta, alpha)&
                          & * (tmpRho(beta, kpa, iSpin) * tmpRho(alpha, mu, iSpin)&
                          & + tmpRho(alpha, kpa, iSpin) * tmpRho(beta, mu, iSpin))
                    end do
                  end do
                end do
                tmpforce(:) = tmpforce + tmpmultvar1 * sPrimeTmp(ccc, kkk, :)
                tmpforce_r(:) = tmpforce_r + tmpmultvar1 * sPrimeTmp2(kkk, ccc, :)
                tmpforce2 = tmpforce2 + tmpmultvar1 * tmpOvr(kpa, mu)
              end do
            end do

            ! C /= K
            if (iAtK /= iAtC) then
              if (iAtB /= iAtA) then
                tmpforce(:) = tmpforce * tmpgamma2
                tmpforce_r(:) = tmpforce_r * tmpgamma2
                tmpforce(:) = tmpforce + tmpforce2 * (gammaPrimeTmp(:, iAtK, iAtA)&
                    & + gammaPrimeTmp(:, iAtK, iAtB))
                tmpforce_r(:) = tmpforce_r + tmpforce2 * (gammaPrimeTmp(:, iAtC, iAtA)&
                    & + gammaPrimeTmp(:, iAtC, iAtB))
              else
                tmpforce(:) = tmpforce * tmpgamma1
                tmpforce_r(:) = tmpforce_r * tmpgamma1
                tmpforce(:) = tmpforce + tmpforce2 * gammaPrimeTmp(:, iAtK, iAtA)
                tmpforce_r(:) = tmpforce_r + tmpforce2 * gammaPrimeTmp(:, iAtC, iAtA)
              end if
            else
              if (iAtB /= iAtA) then
                tmpforce(:) = tmpforce + tmpforce2 * (gammaPrimeTmp(:, iAtK, iAtA)&
                    & + gammaPrimeTmp(:, iAtK, iAtB))
              else
                tmpforce(:) = tmpforce + tmpforce2 * (gammaPrimeTmp(:, iAtK, iAtA))
              end if
            end if
            tmpderiv(:, iAtK) = tmpderiv(:, iAtK) + tmpforce
            tmpderiv(:, iAtC) = tmpderiv(:, iAtC) + tmpforce_r
          end do loopA
        end do loopB
      end do loopC
    end do loopK

    if (this%tREKS) then
      gradients(:,:) = gradients - 0.5_dp * this%camBeta * tmpderiv
    else
      gradients(:,:) = gradients - 0.25_dp * this%camBeta * nSpin * tmpderiv
    end if


  contains

    !> Initialise the
    subroutine allocateAndInit(tmpOvr, tmpRho, gammaPrimeTmp, tmpderiv)

      !> Storage for the overlap
      real(dp), allocatable, intent(inout) :: tmpOvr(:,:)

      !> storage for density matrix
      real(dp), allocatable, intent(inout) :: tmpRho(:,:,:)

      !> storage for derivative of gamma interaction, shape: [nCoords (x,y,z), nAtom0, nAtom0]
      real(dp), allocatable, intent(inout) :: gammaPrimeTmp(:,:,:)

      !> workspace for the derivatives
      real(dp), allocatable, intent(inout) :: tmpderiv(:,:)

      !! Holds long-range gamma derivatives of a single interaction
      real(dp) :: tmp(3)

      !!
      integer :: iSpin, iAt1, iAt2, nAtom0

      nAtom0 = size(this%species0)

      allocate(tmpOvr(size(ovrlapMat, dim=1), size(ovrlapMat, dim=1)))
      allocate(tmpRho(size(deltaRho, dim=1), size(deltaRho, dim=1), size(deltaRho, dim=3)))
      allocate(gammaPrimeTmp(3, nAtom0, nAtom0))
      allocate(tmpderiv(3, size(gradients, dim=2)))
      tmpOvr = ovrlapMat
      tmpRho = deltaRho

      call symmetrizeHS(tmpOvr)

      do iSpin = 1, size(deltaRho, dim=3)
        call symmetrizeHS(tmpRho(:,:, iSpin))
      end do

      ! precompute the gamma derivatives
      gammaPrimeTmp(:,:,:) = 0.0_dp
      do iAt1 = 1, nAtom0
        do iAt2 = 1, nAtom0
          if (iAt1 /= iAt2) then
            call getLrGammaPrimeValue_cluster(this, tmp, iAt1, iAt2, coords0, species0)
            gammaPrimeTmp(:, iAt1, iAt2) = tmp
          end if
        end do
      end do
    end subroutine allocateAndInit

  end subroutine addLrGradients_cluster


  !> Adds gradients due to long-range HF-contribution (Gamma-point version).
  subroutine addLrGradients_gamma(this, deltaRhoSqr, skOverCont, symNeighbourList,&
      & nNeighbourCamSym, iSquare, orb, derivator, gradients)

    !> Instance
    type(TRangeSepFunc), intent(in), target :: this

    !> Square (unpacked) delta density matrix
    real(dp), intent(in) :: deltaRhoSqr(:,:,:)

    !> Sparse overlap container
    type(TSlakoCont), intent(in) :: skOverCont

    !> list of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Atom blocks from sparse, real-space overlap matrices S_{\alpha\mu}, S_{\beta\nu}
    real(dp), pointer :: pSam(:,:), pSbn(:,:)

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descA(descLen), descB(descLen), descM(descLen), descN(descLen)

    !! Temporary storages
    real(dp), allocatable :: tmpDeltaRhoSqr(:,:,:), tmpGradients(:,:)

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Overlap matrix elements
    real(dp) :: Sam, Sbn

    !! Density matrix elements
    real(dp) :: dPmn, dPab

    !! Overlap derivatives
    real(dp), dimension(orb%mOrb, orb%mOrb, 3) :: SbnPrimeKequalsB, SbnPrimeKequalsN
    real(dp), dimension(orb%mOrb, orb%mOrb, 3) :: SamPrimeKequalsA, SamPrimeKequalsM

    !! \tilde{\gamma}_{\mu\nu}, \tilde{\gamma}_{\mu\beta},
    !! \tilde{\gamma}_{\alpha\nu}, \tilde{\gamma}_{\alpha\beta}
    real(dp) :: gammaMN, gammaAN, gammaAB, gammaMNMB, gammaTot

    !! 1st derivatives of
    !! \tilde{\gamma}_{\mu\nu}, \tilde{\gamma}_{\mu\beta},
    !! \tilde{\gamma}_{\alpha\nu}, \tilde{\gamma}_{\alpha\beta}
    real(dp), dimension(3) :: dGammaMN, dGammaMB, dGammaAN, dGammaAB, dGammaMNMB, dGammaTot

    !! Atom to calculate energy gradient components for
    integer :: iAtK

    !! Atom indices (central cell)
    integer :: iAtM, iAtN

    !! Neighbour indices (+corresponding atom indices)
    integer :: iNeighN, iNeighM, iAtA, iAtB

    !! Folded (to central cell) atom indices
    integer :: iAtAfold, iAtBfold

    !! Auxiliary variables for setting up 2D pointer to sparse overlap
    integer :: ind, nOrbAt, nOrbNeigh

    !! Spin index and total number of spin channels
    integer :: iSpin, nSpin

    !! Orbital indices
    integer :: mu, nu, alpha, beta

    !! Product dPmn * Sam
    real(dp) :: dPmnSam

    !! Product dPmn * gammaTot
    real(dp) :: dPmnGammaTot

    !! Product Sam * Sbn * dPab * dPmn
    real(dp) :: dPabdPmnSamSbn

    !! Product dPab * dPmn * gammaTot
    real(dp) :: dPabdPmnGammaTot

    !! Composite index iAtM/iAtN
    integer :: ii, iAtMN(2, size(this%species0)**2)

    nAtom0 = size(this%species0)
    nSpin = size(deltaRhoSqr, dim=3)

    ! Build up composite index iAtMN for collapsing iAtM and iAtN
    ind = 1
    loopiM: do iAtM = 1, nAtom0
      loopiN: do iAtN = 1, nAtom0
        iAtMN(1, ind) = iAtM
        iAtMN(2, ind) = iAtN
        ind = ind + 1
      end do loopiN
    end do loopiM

    ! allocate gradient contribution
    allocate(tmpGradients(3, size(gradients, dim=2)))
    tmpGradients(:,:) = 0.0_dp

    tmpDeltaRhoSqr = deltaRhoSqr
    do iSpin = 1, nSpin
      call symmetrizeHS(tmpDeltaRhoSqr(:,:, iSpin))
    end do

    loopK: do iAtK = 1, nAtom0
      loopMN: do ii = 1, nAtom0**2
        iAtM = iAtMN(1, ii)
        iAtN = iAtMN(2, ii)
        descM = getDescriptor(iAtM, iSquare)
        descN = getDescriptor(iAtN, iSquare)
        ! \tilde{\gamma}_{\mu\nu}(\vec{g})
        gammaMN = this%lrGammaEval0(iAtM, iAtN)
        dGammaMN(:) = 0.0_dp
        if (iAtK == iAtM .and. iAtM /= iAtN) then
          dGammaMN(:) = this%lrdGammaEval0(iAtM, iAtN, :)
        elseif (iAtK == iAtN .and. iAtM /= iAtN) then
          dGammaMN(:) = -this%lrdGammaEval0(iAtM, iAtN, :)
        end if
        loopB: do iNeighN = 0, nNeighbourCamSym(iAtN)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighN, iAtN)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          descB = getDescriptor(iAtBfold, iSquare)
          ! get 2D pointer to Sbn overlap block
          ind = symNeighbourList%iPair(iNeighN, iAtN) + 1
          nOrbAt = descN(iNOrb)
          nOrbNeigh = descB(iNOrb)
          pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          ! \tilde{\gamma}_{\mu\beta}
          gammaMNMB = gammaMN + this%lrGammaEval0(iAtM, iAtBfold)
          dGammaMB(:) = 0.0_dp
          if (iAtK == iAtM .and. iAtM /= iAtBfold) then
            dGammaMB(:) = this%lrdGammaEval0(iAtM, iAtBfold, :)
          elseif (iAtK == iAtBfold .and. iAtM /= iAtBfold) then
            dGammaMB(:) = -this%lrdGammaEval0(iAtM, iAtBfold, :)
          end if
          dGammaMNMB(:) = dGammaMN + dGammaMB
          SbnPrimeKequalsB(:,:,:) = 0.0_dp
          SbnPrimeKequalsN(:,:,:) = 0.0_dp
          if (iAtK == iAtBfold .and. iAtN /= iAtBfold) then
            call derivator%getFirstDeriv(SbnPrimeKequalsB, skOverCont, this%rCoords,&
                & symNeighbourList%species, iAtB, iAtN, orb)
          elseif (iAtK == iAtN .and. iAtN /= iAtBfold) then
            call derivator%getFirstDeriv(SbnPrimeKequalsN, skOverCont, this%rCoords,&
                & symNeighbourList%species, iAtN, iAtB, orb)
          end if
          loopA: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            descA = getDescriptor(iAtAfold, iSquare)
            ! \tilde{\gamma}_{\alpha\nu}
            gammaAN = this%lrGammaEval0(iAtAfold, iAtN)
            dGammaAN(:) = 0.0_dp
            if (iAtK == iAtAfold .and. iAtAfold /= iAtN) then
              dGammaAN(:) = this%lrdGammaEval0(iAtAfold, iAtN, :)
            elseif (iAtK == iAtN .and. iAtAfold /= iAtN) then
              dGammaAN(:) = -this%lrdGammaEval0(iAtAfold, iAtN, :)
            end if
            ! \tilde{\gamma}_{\alpha\beta}
            gammaAB = this%lrGammaEval0(iAtAfold, iAtBfold)
            dGammaAB(:) = 0.0_dp
            if (iAtK == iAtAfold .and. iAtAfold /= iAtBfold) then
              dGammaAB(:) = this%lrdGammaEval0(iAtAfold, iAtBfold, :)
            elseif (iAtK == iAtBfold .and. iAtAfold /= iAtBfold) then
              dGammaAB(:) = -this%lrdGammaEval0(iAtAfold, iAtBfold, :)
            end if
            ! get 2D pointer to S_{\alpha\mu}(\vec{h}) overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            gammaTot = gammaMNMB + gammaAN + gammaAB
            dGammaTot(:) = dGammaMNMB + dGammaAN + dGammaAB

            SamPrimeKequalsA(:,:,:) = 0.0_dp
            SamPrimeKequalsM(:,:,:) = 0.0_dp
            if (iAtK == iAtAfold .and. iAtAfold /= iAtM) then
              call derivator%getFirstDeriv(SamPrimeKequalsA, skOverCont, this%rCoords,&
                  & symNeighbourList%species, iAtA, iAtM, orb)
            elseif (iAtK == iAtM .and. iAtAfold /= iAtM) then
              call derivator%getFirstDeriv(SamPrimeKequalsM, skOverCont, this%rCoords,&
                  & symNeighbourList%species, iAtM, iAtA, orb)
            end if

            do mu = 1, descM(iNOrb)
              do nu = 1, descN(iNOrb)
                do iSpin = 1, nSpin
                  dPmn = tmpDeltaRhoSqr(descM(iStart) + mu - 1, descN(iStart) + nu - 1, iSpin)
                  dPmnGammaTot = dPmn * gammaTot
                  do alpha = 1, descA(iNOrb)
                    Sam = pSam(alpha, mu)
                    dPmnSam = dPmn * Sam
                    do beta = 1, descB(iNOrb)
                      Sbn = pSbn(beta, nu)
                      dPab = tmpDeltaRhoSqr(descA(iStart) + alpha - 1,&
                          & descB(iStart) + beta - 1, iSpin)
                      dPabdPmnSamSbn = dPmnSam * Sbn * dPab
                      dPabdPmnGammaTot = dPab * dPmnGammaTot

                      tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                          & + dPabdPmnGammaTot * Sam * SbnPrimeKequalsB(nu, beta, :)
                      tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                          & + dPabdPmnGammaTot * Sam * SbnPrimeKequalsN(beta, nu, :)

                      tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                          & + dPabdPmnGammaTot * Sbn * SamPrimeKequalsA(mu, alpha, :)
                      tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                          & + dPabdPmnGammaTot * Sbn * SamPrimeKequalsM(alpha, mu, :)

                      tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                          & + dPabdPmnSamSbn * dGammaTot
                    end do
                  end do
                end do
              end do
            end do

          end do loopA
        end do loopB
      end do loopMN
    end do loopK

    if (this%tREKS) then
      gradients(:,:) = gradients - 0.125_dp * this%camBeta * tmpGradients
    else
      gradients(:,:) = gradients - 0.0625_dp * this%camBeta * nSpin * tmpGradients
    end if

  end subroutine addLrGradients_gamma


  !> evaluate the LR-Energy contribution directly. Very slow, use addLrEnergy instead.
  function evaluateLrEnergyDirect_cluster(this, env, deltaRho, overlap, iSquare) result(energy)

    !> instance of LR
    class(TRangeSepFunc), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> square density matrix
    real(dp), intent(in) :: deltaRho(:,:)

    !> square overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    !> Dense matrix atom indexing
    integer, intent(in) :: iSquare(:)

    !> resulting energy
    real(dp) :: energy

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Atom indices, where orbitals \mu and \nu are located
    integer :: iAt1, iAt2, mu, nu, alpha, beta

    !!
    real(dp), allocatable :: tmpOvr(:,:), tmpDRho(:,:)
    real(dp) :: tmp

    call env%globalTimer%startTimer(globalTimers%energyEval)

    nAtom0 = size(this%species0)

    allocate(tmpOvr(size(overlap, dim=1), size(overlap, dim=1)))
    allocate(tmpDRho(size(deltaRho, dim=1), size(deltaRho, dim=1)))
    tmpOvr(:,:) = overlap
    tmpDRho(:,:) = deltaRho

    call symmetrizeHS(tmpOvr)
    call symmetrizeHS(tmpDRho)

    energy = 0.0_dp
    do iAt1 = 1, nAtom0
      do iAt2 = 1, nAtom0
        tmp = 0.0_dp
        do mu = iSquare(iAt1), iSquare(iAt1 + 1) - 1
          do nu = iSquare(iAt2), iSquare(iAt2 + 1) - 1
            do alpha = 1, size(tmpOvr, dim=1)
              do beta = 1, size(tmpOvr, dim=1)
                tmp = tmp + (&
                    & tmpDRho(alpha, beta) * tmpDRho(mu, nu)&
                    & + tmpDRho(mu,beta) * tmpDRho(alpha,nu)) * tmpOvr(mu,alpha) * tmpOvr(nu,beta)
              end do
            end do
          end do
        end do
        energy = energy + tmp * this%lrGammaEval0(iAt1, iAt2)
      end do
    end do
    energy = -energy / 8.0_dp

    call env%globalTimer%stopTimer(globalTimers%energyEval)

  end function evaluateLrEnergyDirect_cluster


  !> Returns the array of atomic species in central cell.
  subroutine getCentralCellSpecies(this, species0)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> 1D array for output, will be allocated
    integer, intent(out), allocatable :: species0(:)

    species0 = this%species0

  end subroutine getCentralCellSpecies


  !> Returns tabulated long-range gamma integrals (non-periodic systems only).
  subroutine getLrGammaCluster(this, lrGamma0)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> Long-range gamma integrals in AO basis
    real(dp), intent(out) :: lrGamma0(:,:)

    lrGamma0(:,:) = this%lrGammaEval0

  end subroutine getLrGammaCluster


  !> Calculates long-range gamma derivative integrals (non-periodic).
  subroutine getLrGammaDerivCluster(this, coords0, species0, lrGammaDeriv0)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> Atomic coordinates in central cell
    real(dp), intent(in) :: coords0(:,:)

    !> Species of all atoms in central cell
    integer, intent(in) :: species0(:)

    !> Long-range gamma derivative integrals
    real(dp), intent(out) :: lrGammaDeriv0(:,:,:)

    !! Holds long-range gamma derivatives of a single interaction
    real(dp) :: tmp(3)

    !! Number of atoms and indices of interacting atoms
    integer :: nAtom0, iAt1, iAt2

    nAtom0 = size(lrGammaDeriv0, dim=1)

    do iAt1 = 1, nAtom0
      do iAt2 = 1, nAtom0
        if (iAt1 /= iAt2) then
          call getLrGammaPrimeValue_cluster(this, tmp, iAt1, iAt2, coords0, species0)
          lrGammaDeriv0(iAt2, iAt1, :) = tmp
        end if
      end do
    end do

  end subroutine getLrGammaDerivCluster


  !!
  !! Full-range Hartree-Fock related Routines.
  !!

  !> Interface routine for the full-range Hartree-Fock contribution to the Hamiltonian.
  subroutine addHfHamiltonian_cluster(this, densSqr, over, iNeighbour, nNeighbourCam,&
      & iSquare, iPair, orb, HH, overlap)

    !> Instance
    type(TRangeSepFunc), intent(inout) :: this

    ! Neighbour based screening

    !> Square (unpacked) density matrix
    real(dp), intent(in), target :: densSqr(:,:)

    !> Sparse (packed) overlap matrix.
    real(dp), intent(in) :: over(:)

    !> Neighbour indices
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCam(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Position of each (neighbour, atom) pair in the sparse matrix.
    !> Shape: (0:maxNeighbour, nAtom)
    integer, intent(in) :: iPair(0:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) Hamiltonian to be updated
    real(dp), intent(inout), target :: HH(:,:)

    ! Threshold based screening

    !> Square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    select case(this%rsAlg)
    case (rangeSepTypes%threshold)
      call error('Thresholded algorithm not yet implemented for HF non-periodic systems.')
    case (rangeSepTypes%neighbour)
      call error('Neighbour based algorithm not yet implemented for HF non-periodic systems.')
    case (rangeSepTypes%matrixBased)
      call addHfHamiltonianMatrix_cluster(this, iSquare, overlap, densSqr, HH)
    end select

  end subroutine addHfHamiltonian_cluster


  !> Updates Hamiltonian with Hartree-Fock contribution using matrix-matrix multiplications.
  subroutine addHfHamiltonianMatrix_cluster(this, iSquare, overlap, densSqr, HH)

    !> Instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Square (unpacked) overlap matrix.
    real(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout) :: HH(:,:)

    real(dp), allocatable :: Smat(:,:)
    real(dp), allocatable :: Dmat(:,:)
    real(dp), allocatable :: hfGammaAO(:,:)
    real(dp), allocatable :: Hhf(:,:)

    integer :: nOrb

    nOrb = size(overlap, dim=1)

    allocate(Smat(nOrb, nOrb))
    allocate(Dmat(nOrb, nOrb))
    allocate(hfGammaAO(nOrb, nOrb))
    allocate(Hhf(nOrb, nOrb))

    call allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, hfGammaAO)
    call evaluateHamiltonian(this, Smat, Dmat, hfGammaAO, Hhf)
    HH(:,:) = HH + this%camAlpha * Hhf
    this%hfEnergy = this%hfEnergy + 0.5_dp * sum(Dmat * Hhf)

  contains

    !> Set up storage and get orbital-by-orbital gamma matrix
    subroutine allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, hfGammaAO)

      !> Instance
      type(TRangeSepFunc), intent(in) :: this

      !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
      integer, intent(in) :: iSquare(:)

      !> Square (unpacked) overlap matrix.
      real(dp), intent(in) :: overlap(:,:)

      !> Square (unpacked) density matrix
      real(dp), intent(in) :: densSqr(:,:)

      !> Square (unpacked) Hamiltonian to be updated.
      real(dp), intent(inout) :: HH(:,:)

      !> Symmetrized square overlap matrix
      real(dp), intent(out) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(out) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      real(dp), intent(out) :: hfGammaAO(:,:)

      integer :: nAtom0, iAt, jAt

      nAtom0 = size(this%hfGammaEval0, dim=1)

      ! Symmetrize Hamiltonian, overlap, density matrices
      call symmetrizeHS(HH)
      Smat(:,:) = overlap
      call symmetrizeHS(Smat)
      Dmat(:,:) = densSqr
      call symmetrizeHS(Dmat)

      ! Get long-range gamma variable
      hfGammaAO(:,:) = 0.0_dp
      do iAt = 1, nAtom0
        do jAt = 1, nAtom0
          hfGammaAO(iSquare(jAt):iSquare(jAt + 1) - 1, iSquare(iAt):iSquare(iAt + 1) - 1) =&
              & this%hfGammaEval0(jAt, iAt)
        end do
      end do

    end subroutine allocateAndInit


    !> Evaluate the Hamiltonian using GEMM operations.
    subroutine evaluateHamiltonian(this, Smat, Dmat, hfGammaAO, Hhf)

      !> Instance
      type(TRangeSepFunc), intent(in) :: this

      !> Symmetrized square overlap matrix
      real(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      real(dp), intent(in) :: hfGammaAO(:,:)

      !> Symmetrized full-range Hartree-Fock Hamiltonian matrix
      real(dp), intent(out) :: Hhf(:,:)

      real(dp), allocatable :: Hmat(:,:)
      real(dp), allocatable :: tmpMat(:,:)

      integer :: nOrb

      nOrb = size(Smat,dim=1)

      allocate(Hmat(nOrb,nOrb))
      allocate(tmpMat(nOrb,nOrb))

      Hhf(:,:) = 0.0_dp

      call gemm(tmpMat, Smat, Dmat)
      call gemm(Hhf, tmpMat, Smat)
      Hhf(:,:) = Hhf * hfGammaAO

      tmpMat(:,:) = tmpMat * hfGammaAO
      call gemm(Hhf, tmpMat, Smat, alpha=1.0_dp, beta=1.0_dp)

      Hmat(:,:) = Dmat * hfGammaAO
      call gemm(tmpMat, Smat, Hmat)
      call gemm(Hhf, tmpMat, Smat, alpha=1.0_dp, beta=1.0_dp)

      call gemm(tmpMat, Dmat, Smat)
      tmpMat(:,:) = tmpMat * hfGammaAO
      call gemm(Hhf, Smat, tmpMat, alpha=1.0_dp, beta=1.0_dp)

      if (this%tSpin .or. this%tREKS) then
        Hhf(:,:) = -0.25_dp * Hhf
      else
        Hhf(:,:) = -0.125_dp * Hhf
      end if

    end subroutine evaluateHamiltonian

  end subroutine addHfHamiltonianMatrix_cluster


  !> Add the full-range Hartree-Fock Energy contribution to the total energy.
  subroutine addHfEnergy(this, env, energy)

    !> Instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Total energy
    real(dp), intent(inout) :: energy

    !! Total full-range Hartree-Fock energy of all MPI processes
    real(dp) :: hfEnergy

  #:if WITH_MPI
    call mpifx_allreduce(env%mpi%globalComm, this%hfEnergy, hfEnergy, MPI_SUM)
  #:else
    hfEnergy = this%hfEnergy
  #:endif

    energy = energy + this%camAlpha * hfEnergy

    ! hack for spin unrestricted calculation
    this%hfEnergy = 0.0_dp

  end subroutine addHfEnergy


  !> Calculates analytical full-range Hartree-Fock gamma.
  function getAnalyticalHfGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma
    real(dp) :: gamma

    real(dp) :: tauA, tauB
    real(dp) :: tmp, tau

    tauA = 3.2_dp * this%hubbu(iSp1)
    tauB = 3.2_dp * this%hubbu(iSp2)

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff) then
        tau = 0.5_dp * (tauA + tauB)
        gamma = tau * 0.3125_dp
      else
        call error("Error(RangeSep): R = 0, Ua != Ub")
      end if
    else
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff ) then
        tauA = 0.5_dp * (tauA + tauB)
        tmp = ((dist * tauA)**3 / 48.0_dp + 0.1875_dp * (dist * tauA)**2 +&
            & 0.6875_dp * (dist * tauA) + 1.0_dp) * exp(-tauA * dist) / dist
        gamma = 1.0_dp / dist - tmp
      ! off-site, Ua != Ub
      else
        gamma = 1.0_dp / dist&
            & - getYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getYGammaSubPart(tauB, tauA, dist, 0.0_dp)
      end if
    end if

  end function getAnalyticalHfGammaValue


  !> Derivative of analytical full-range Hartree-Fock gamma.
  function getdAnalyticalHfGammaDeriv(this, Sp1, Sp2, dist) result(dGamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: Sp1

    !> Second species
    integer, intent(in) :: Sp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting d gamma / d dist
    real(dp) :: dGamma

    real(dp) :: tauA, tauB
    real(dp) :: dTmp

    tauA = 3.2_dp * this%hubbu(Sp1)
    tauB = 3.2_dp * this%hubbu(Sp2)

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff) then
        dGamma = 0.0_dp
      else
        call error("Error(RangeSep): R = 0, Ua != Ub")
      end if
    else
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff) then
        tauA = 0.5_dp * (tauA + tauB)

        dTmp = &
            & (2.0_dp * dist * tauA**3 / 48.0_dp + 0.1875_dp * tauA**2 - 1.0_dp / dist**2)&
            & * exp(-tauA * dist) - (dist**2 * tauA**3 / 48.0_dp + 0.1875_dp * dist * tauA**2&
            & + 0.6875_dp * tauA + 1.0_dp / dist) * tauA * exp(-tauA * dist)

        dGamma = -1.0_dp / dist**2 - dtmp

      ! off-site, Ua != Ub
      else
        dGamma = -1.0_dp / (dist**2)&
            & - getdYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getdYGammaSubPart(tauB, tauA, dist, 0.0_dp)
      end if
    end if

  end function getdAnalyticalHfGammaDeriv


  !> !> Returns the derivative of full-range Hartree-Fock gamma for iAtom1, iAtom2
  !! (non-periodic systems only).
  subroutine getHfGammaPrimeValue_cluster(this, grad, iAtom1, iAtom2, coords0, species0)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Gradient of gamma between atoms
    real(dp), intent(out) :: grad(3)

    !> First atom
    integer, intent(in) :: iAtom1

    !> Second atom
    integer, intent(in) :: iAtom2

    !> Coordinates of atoms in central cell
    real(dp), intent(in) :: coords0(:,:)

    !> List of all atomic species in central cell
    integer, intent(in) :: species0(:)

    !! Species index of first and second atom
    integer :: iSp1, iSp2

    !! Distance(-vector) of the two atoms
    real(dp) :: vect(3), dist

    iSp1 = species0(iAtom1)
    iSp2 = species0(iAtom2)

    ! analytical derivatives
    vect(:) = coords0(:, iAtom1) - coords0(:, iAtom2)
    dist = sqrt(sum(vect**2))
    vect(:) = vect / dist
    grad(:) = vect * getdAnalyticalHfGammaDeriv(this, iSp1, iSp2, dist)

  end subroutine getHfGammaPrimeValue_cluster


  !> !> Returns the derivative of full-range Hartree-Fock gamma for iAtom1, iAtom2
  !! (periodic systems).
  subroutine getHfGammaPrimeValue_periodic(this, grad, iAtom1, iAtom2, coords, species0,&
      & img2CentCell)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Gradient of gamma between atoms
    real(dp), intent(out) :: grad(3)

    !> First atom
    integer, intent(in) :: iAtom1

    !> Second atom
    integer, intent(in) :: iAtom2

    !> Coordinates of atoms, including periodic images
    real(dp), intent(in) :: coords(:,:)

    !> List of all atomic species in central cell
    integer, intent(in) :: species0(:)

    !> Map images of atoms to the central cell
    integer, intent(in) :: img2CentCell(:)

    !! Species index of first and second atom
    integer :: iSp1, iSp2

    !! Distance(-vector) of the two atoms
    real(dp) :: vect(3), dist

    iSp1 = species0(img2CentCell(iAtom1))
    iSp2 = species0(img2CentCell(iAtom2))

    ! analytical derivatives
    vect(:) = coords(:, iAtom1) - coords(:, iAtom2)
    dist = sqrt(sum(vect**2))
    vect(:) = vect / dist
    grad(:) = vect * getdAnalyticalHfGammaDeriv(this, iSp1, iSp2, dist)

  end subroutine getHfGammaPrimeValue_periodic


  !> Returns full-range HartreeFock gamma integrals (non-periodic).
  pure subroutine getHfGammaCluster(this, hfGamma0)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> Full-range HartreeFock gamma integrals in AO basis
    real(dp), intent(out) :: hfGamma0(:,:)

    hfGamma0(:,:) = this%hfGammaEval0

  end subroutine getHfGammaCluster


  !> Calculate long-range gamma derivative integrals.
  subroutine getHfGammaDerivCluster(this, coords0, species0, hfGammaDeriv0)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> Atomic coordinates in central cell
    real(dp), intent(in) :: coords0(:,:)

    !> Species of all atoms in central cell
    integer, intent(in) :: species0(:)

    !> Long-range gamma derivative integrals
    real(dp), intent(out) :: hfGammaDeriv0(:,:,:)

    !! Holds long-range gamma derivatives of a single interaction
    real(dp) :: tmp(3)

    !! Number of atoms and indices of interacting atoms
    integer :: nAtom0, iAt1, iAt2

    nAtom0 = size(hfGammaDeriv0, dim=1)

    do iAt1 = 1, nAtom0
      do iAt2 = 1, nAtom0
        if (iAt1 /= iAt2) then
          call getHfGammaPrimeValue_cluster(this, tmp, iAt1, iAt2, coords0, species0)
          hfGammaDeriv0(iAt2, iAt1, :) = tmp
        end if
      end do
    end do

  end subroutine getHfGammaDerivCluster


  !> Adds gradients due to full-range Hartree-Fock contribution.
  subroutine addHfGradients_cluster(this, gradients, derivator, deltaRho, skOverCont,&
      & coords, species, orb, iSquare, overlapMat, iNeighbour, nNeighbourSK)

    !> Class instance
    type(TRangeSepFunc), intent(in) :: this

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> Density matrix difference from reference q0
    real(dp), intent(in) :: deltaRho(:,:,:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Chemical species of atoms
    integer, intent(in) :: species(:)

    !> Orbital information for system
    type(TOrbitals), intent(in) :: orb

    !> Index for dense arrays
    integer, intent(in) :: iSquare(:)

    !> Overlap matrix
    real(dp), intent(in) :: overlapMat(:,:)

    !> Neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of atoms neighbouring each site where the overlap is non-zero
    integer, intent(in) :: nNeighbourSK(:)

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    integer :: nAtom0, iAtK, iNeighK, iAtB, iNeighB, iAtC, iAtA, kpa
    real(dp) :: tmpgamma1, tmpgamma2
    real(dp) :: tmpforce(3), tmpforce_r(3), tmpforce2, tmpmultvar1
    integer :: nSpin, iSpin, mu, alpha, beta, ccc, kkk
    real(dp) :: sPrimeTmp(orb%mOrb, orb%mOrb, 3)
    real(dp) :: sPrimeTmp2(orb%mOrb, orb%mOrb, 3)
    real(dp), allocatable :: gammaPrimeTmp(:,:,:), tmpOvr(:,:), tmpRho(:,:,:), tmpderiv(:,:)

    nAtom0 = size(this%species0)
    nSpin = size(deltaRho, dim=3)
    call allocateAndInit(tmpOvr, tmpRho, gammaPrimeTmp, tmpderiv)

    ! sum K
    loopK: do iAtK = 1, nAtom0
      ! C >= K
      loopC: do iNeighK = 0, nNeighbourSK(iAtK)
        iAtC = iNeighbour(iNeighK, iAtK)
        ! evaluate the ovr_prime
        sPrimeTmp2(:,:,:) = 0.0_dp
        sPrimeTmp(:,:,:) = 0.0_dp
        if (iAtK /= iAtC) then
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species, iAtK, iAtC, orb)
          call derivator%getFirstDeriv(sPrimeTmp2, skOverCont, coords, species, iAtC, iAtK, orb)
        end if
        loopB: do iAtB = 1, nAtom0
          ! A > B
          loopA: do iNeighB = 0, nNeighbourSK(iAtB)
            iAtA = iNeighbour(iNeighB, iAtB)
            tmpgamma1 = this%hfGammaEval0(iAtK, iAtB) + this%hfGammaEval0(iAtC, iAtB)
            tmpgamma2 = tmpgamma1 + this%hfGammaEval0(iAtK, iAtA) + this%hfGammaEval0(iAtC, iAtA)
            tmpforce(:) = 0.0_dp
            tmpforce_r(:) = 0.0_dp
            tmpforce2 = 0.0_dp
            ccc = 0
            do mu = iSquare(iAtC), iSquare(iAtC + 1) - 1
              ccc = ccc + 1
              kkk = 0
              do kpa = iSquare(iAtK), iSquare(iAtK + 1) - 1
                kkk = kkk + 1
                tmpmultvar1 = 0.0_dp
                do iSpin = 1, nSpin
                  do alpha = iSquare(iAtA), iSquare(iAtA + 1) - 1
                    do beta = iSquare(iAtB), iSquare(iAtB + 1) - 1
                      tmpmultvar1 = tmpmultvar1&
                          & + tmpOvr(beta, alpha) * (tmpRho(beta, kpa, iSpin)&
                          & * tmpRho(alpha, mu, iSpin) + tmpRho(alpha, kpa, iSpin)&
                          & * tmpRho(beta, mu, iSpin))
                    end do
                  end do
                end do
                tmpforce(:) = tmpforce + tmpmultvar1 * sPrimeTmp(ccc, kkk, :)
                tmpforce_r(:) = tmpforce_r + tmpmultvar1 * sPrimeTmp2(kkk, ccc, :)
                tmpforce2 = tmpforce2 + tmpmultvar1 * tmpOvr(kpa, mu)
              end do
            end do

            ! C /= K
            if (iAtK /= iAtC) then
              if (iAtB /= iAtA) then
                tmpforce(:) = tmpforce * tmpgamma2
                tmpforce_r(:) = tmpforce_r * tmpgamma2
                tmpforce(:) = tmpforce + tmpforce2 * (gammaPrimeTmp(:, iAtK, iAtA)&
                    & + gammaPrimeTmp(:, iAtK, iAtB))
                tmpforce_r(:) = tmpforce_r + tmpforce2 * (gammaPrimeTmp(:, iAtC, iAtA)&
                    & + gammaPrimeTmp(:, iAtC, iAtB))
              else
                tmpforce(:) = tmpforce * tmpgamma1
                tmpforce_r(:) = tmpforce_r * tmpgamma1
                tmpforce(:) = tmpforce + tmpforce2 * gammaPrimeTmp(:, iAtK, iAtA)
                tmpforce_r(:) = tmpforce_r + tmpforce2 * gammaPrimeTmp(:, iAtC, iAtA)
              end if
            else
              if (iAtB /= iAtA) then
                tmpforce(:) = tmpforce + tmpforce2 * (gammaPrimeTmp(:, iAtK, iAtA)&
                    & + gammaPrimeTmp(:, iAtK, iAtB))
              else
                tmpforce(:) = tmpforce + tmpforce2 * gammaPrimeTmp(:, iAtK, iAtA)
              end if
            end if
            tmpderiv(:, iAtK) = tmpderiv(:, iAtK) + tmpforce
            tmpderiv(:, iAtC) = tmpderiv(:, iAtC) + tmpforce_r
          end do loopA
        end do loopB
      end do loopC
    end do loopK

    if (this%tREKS) then
      gradients(:,:) = gradients - 0.5_dp * this%camAlpha * tmpderiv
    else
      gradients(:,:) = gradients - 0.25_dp * this%camAlpha * nSpin * tmpderiv
    end if

  contains

    !> Initialise.
    subroutine allocateAndInit(tmpOvr, tmpRho, gammaPrimeTmp, tmpderiv)

      !> Storage for the overlap
      real(dp), intent(inout), allocatable :: tmpOvr(:,:)

      !> Storage for density matrix
      real(dp), intent(inout), allocatable :: tmpRho(:,:,:)

      !> Storage for derivative of gamma interaction
      real(dp), intent(inout), allocatable :: gammaPrimeTmp(:,:,:)

      !> Workspace for the derivatives
      real(dp), intent(inout), allocatable :: tmpderiv(:,:)

      !! Holds gamma derivative for a single iAt1-iAt2 interaction
      real(dp) :: tmp(3)

      !! Spin index (up/down)
      integer :: iSpin

      !! Indices of interacting atoms
      integer :: iAt1, iAt2

      !! Number of atoms in central cell
      integer :: nAtom0

      nAtom0 = size(this%species0)
      allocate(tmpOvr(size(overlapMat, dim=1), size(overlapMat, dim=1)))
      allocate(tmpRho(size(deltaRho, dim=1), size(deltaRho, dim=1), size(deltaRho, dim=3)))
      allocate(gammaPrimeTmp(3, nAtom0, nAtom0))
      allocate(tmpderiv(3, size(gradients, dim=2)))

      tmpOvr(:,:) = overlapMat
      tmpRho(:,:,:) = deltaRho
      gammaPrimeTmp(:,:,:) = 0.0_dp
      tmpderiv(:,:) = 0.0_dp

      call symmetrizeHS(tmpOvr)
      do iSpin = 1, size(deltaRho, dim=3)
        call symmetrizeHS(tmpRho(:,:, iSpin))
      end do

      ! precompute the gamma derivatives
      do iAt1 = 1, nAtom0
        do iAt2 = 1, nAtom0
          if (iAt1 /= iAt2) then
            call getHfGammaPrimeValue_cluster(this, tmp, iAt1, iAt2, coords, species)
            gammaPrimeTmp(:, iAt1, iAt2) = tmp
          end if
        end do
      end do

    end subroutine allocateAndInit

  end subroutine addHfGradients_cluster

end module dftbp_dftb_rangeseparated
