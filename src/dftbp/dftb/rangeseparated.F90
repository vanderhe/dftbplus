!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Contains hybrid xc-functional related routines.
module dftbp_dftb_hybridxc

  use dftbp_common_accuracy, only : dp, tolSameDist, MinHubDiff
  use dftbp_common_constants, only : pi, Bohr__AA
  use dftbp_common_environment, only : TEnvironment, globalTimers
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_nonscc, only : TNonSccDiff
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_dftb_sparse2dense, only : blockSymmetrizeHS, symmetrizeHS, hermitianSquareMatrix,&
      & lowerTriangleSquareMatrix
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
  use dftbp_math_wigner, only : generateWignerSeitzGrid

#:if WITH_OMP
  use omp_lib, only : OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
#:endif

#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip, mpifx_allreduce, mpifx_bcast
#:endif

  implicit none
  private

  public :: THybridXcSKTag, THybridXcFunc, THybridXcFunc_init
  public :: getDirectedCamGammaPrimeValue
  public :: hybridXcFunc, hybridXcAlgo, hybridXcGammaTypes, checkSupercellFoldingMatrix


  !> Returns the (directed) derivative of long-range + full-range Hartree-Fock gamma.
  interface getDirectedCamGammaPrimeValue
    module procedure getDirectedCamGammaPrimeValue_cluster
    module procedure getDirectedCamGammaPrimeValue_periodic
  end interface getDirectedCamGammaPrimeValue


  type :: TIntArray1D

    !> Onedimensional, integer data storage
    integer, allocatable :: array(:)

  end type TIntArray1D


  type :: TRealArray1D

    !> Onedimensional, real data storage
    real(dp), allocatable :: array(:)

  end type TRealArray1D


  type :: TRealArray2D

    !> Twodimensional, real data storage
    real(dp), allocatable :: array(:,:)

  end type TRealArray2D


  !> Enumerator for type of hybrid functional used.
  !! (global hybrid, purely long-range corrected, general CAM range-separated form)
  type :: THybridXcFuncEnum

    !> Global hybrid
    integer :: hyb = 0

    !> Long-range corrected
    integer :: lc = 1

    !> General Coulomb-attenuated method
    integer :: cam = 2

  end type THybridXcFuncEnum


  !> Enumerator for algorithms to build up the range-separated Hamiltonian.
  type :: THybridXcAlgoEnum

    !> Neighbour based
    integer :: neighbour = 0

    !> Threshold based
    integer :: threshold = 1

    !> Matrix based
    integer :: matrixBased = 2

  end type THybridXcAlgoEnum


  !> Enumerator for gamma function types.
  type :: THybridXcGammaTypesEnum

    !> Full, unaltered gamma function (full)
    integer :: full = 0

    !> Truncated gamma function with hard cutoff (truncated)
    integer :: truncated = 1

    !> Truncated and poly5zero damped gamma function (truncated+damping)
    integer :: truncatedAndDamped = 2

    !> Artificially screened gamma function (screened)
    integer :: screened = 3

    !> Artificially screened and poly5zero damped gamma function (screened+damping)
    integer :: screenedAndDamped = 4

    !> Minimum image convention, based on density matrix criteria (mic)
    integer :: mic = 5

  end type THybridXcGammaTypesEnum


  !> Container for enumerated types of hybrid functionals.
  type(THybridXcFuncEnum), parameter :: hybridXcFunc = THybridXcFuncEnum()

  !> Container for enumerated range separation algorithms
  type(THybridXcAlgoEnum), parameter :: hybridXcAlgo = THybridXcAlgoEnum()

  !> Container for enumerated range separation gamma function types
  type(THybridXcGammaTypesEnum), parameter :: hybridXcGammaTypes = THybridXcGammaTypesEnum()


  !> Slater-Koster file HybridXc tag structure.
  type :: THybridXcSKTag

    !> range-separation parameter
    real(dp) :: omega

    !> CAM alpha parameter
    real(dp) :: camAlpha

    !> CAM beta parameter
    real(dp) :: camBeta

  end type THybridXcSKTag


  !> Abstract base type for range-separated functionality, some procedures need to be overridden.
  type, abstract :: THybridXcFunc

    !> Real-space coordinates of atoms (relative units), potentially including periodic images
    real(dp), allocatable :: coords(:,:)

    !> Real-space coordinates of atoms (absolute units), potentially including periodic images
    real(dp), allocatable :: rCoords(:,:)

    !> Evaluated (long-range + Hartree-Fock) gamma of Atom1 and Atom2 (central cell only)
    real(dp), allocatable :: camGammaEval0(:,:), camdGammaEval0(:,:,:)

    !> Evaluated (long-range + Hartree-Fock) gamma in the general k-point case
    type(TRealArray1D), allocatable :: camGammaEvalG(:,:)
    type(TRealArray2D), allocatable :: camdGammaEvalG(:,:)

    !> Number of numerically non-zero gamma's
    type(TIntArray1D), allocatable :: nNonZeroGammaG(:,:)

    !> Range-separation parameter
    real(dp) :: omega

    !> CAM alpha parameter
    real(dp) :: camAlpha

    !> CAM beta parameter
    real(dp) :: camBeta

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

    !> Total (long-range + full-range Hartree-Fock) CAM energy
    real(dp) :: camEnergy

    !> Is this spin restricted (F) or unrestricted (T)
    logical :: tSpin

    !> Is this DFTB/SSR formalism
    logical :: tREKS

    !> Algorithm for range separation screening
    integer :: hybridXcAlg

    !> Hybrid xc-functional type, as extracted from SK-file(s)
    integer :: hybridXcType

    !> Species of atoms in central cell
    integer, allocatable :: species0(:)

    !> Cutoff for real-space g-summation
    real(dp) :: gSummationCutoff

    !> Number of unitcells along each supercell folding direction to substract from MIC Wigner-Seitz
    !! cell construction
    integer :: wignerSeitzReduction

    !> Cutoff for truncated Gamma
    real(dp) :: gammaCutoff

    !> Damping distance for Gamma truncation
    real(dp) :: gammaDamping

    !> Auxiliary gamma damping/screening parameter
    real(dp) :: auxiliaryScreening

    !> Value, 1st and 2nd derivative of gamma integral at damping distance
    real(dp), allocatable :: lrGammaAtDamping(:,:), lrdGammaAtDamping(:,:), lrddGammaAtDamping(:,:)

    !> Value, 1st and 2nd derivative of screened gamma integral at damping distance
    real(dp), allocatable :: lrScreenedGammaAtDamping(:,:), lrdScreenedGammaAtDamping(:,:)
    real(dp), allocatable :: lrddScreenedGammaAtDamping(:,:)

    !> Value, 1st and 2nd derivative of gamma integral at damping distance
    real(dp), allocatable :: hfGammaAtDamping(:,:), hfdGammaAtDamping(:,:), hfddGammaAtDamping(:,:)

    !> Overlap estimates
    type(TRealArray1D), allocatable :: testSquareOver(:)

    !> Descending neighbour indices in terms of overlap estimates
    type(TIntArray1D), allocatable :: overlapIndices(:)

    !> K-point compatible BvK real-space shifts in relative coordinates (units of latVecs)
    real(dp), allocatable :: bvKShifts(:,:)

    !> Supercell folding coefficients (diagonal elements)
    integer, allocatable :: coeffsDiag(:)

    !> Gamma function type (mostly for periodic cases)
    integer :: gammaType

    !> Wigner-Seitz grid points in units of lattice vectors
    integer, allocatable :: wsVectors(:,:)

    !> Translation vectors to lattice cells in units of lattice constants
    real(dp), allocatable :: cellVecsG(:,:)

    !> Vectors to unit cells in absolute units
    real(dp), allocatable :: rCellVecsG(:,:)

  contains

    procedure :: updateCoords_cluster, updateCoords_gamma, updateCoords_kpts

    procedure :: foldToBvK => THybridXcFunc_foldToBvK
    procedure :: foldToBvKIndex => THybridXcFunc_foldToBvKIndex

    procedure :: addCamHamiltonian_cluster
    procedure :: addCamHamiltonian_gamma
    procedure :: addCamHamiltonian_kpts

    procedure :: addCamHamiltonianMatrix_cluster_cmplx

    procedure :: addCamEnergy
    procedure :: addCamEnergy_kpts

    procedure :: addCamGradients_cluster
    procedure :: addCamGradients_gamma
    procedure :: addCamGradients_kpts_ct

    procedure :: getCentralCellSpecies
    procedure :: getCamGammaCluster
    procedure :: getCamGammaDerivCluster

    procedure :: evaluateLrEnergyDirect_cluster

    procedure(gammaFunc), deferred :: getLrGammaValue
    procedure(gammaFunc), deferred :: getLrGammaPrimeValue

    procedure(gammaFunc), deferred :: getHfGammaValue
    procedure(gammaFunc), deferred :: getHfGammaPrimeValue

  end type THybridXcFunc


  abstract interface
    !> Calculates analytical long-range or full-range HF gamma (derivative) of given type.
    function gammaFunc(this, iSp1, iSp2, dist) result(gamma)

      import :: THybridXcFunc, dp

      !> Class instance
      class(THybridXcFunc), intent(in) :: this

      !> First species
      integer, intent(in) :: iSp1

      !> Second species
      integer, intent(in) :: iSp2

      !> Distance between atoms
      real(dp), intent(in) :: dist

      !> Resulting truncated gamma
      real(dp) :: gamma

    end function gammaFunc
  end interface


  !> Base type extension for unaltered analytical gamma functions.
  type, extends(THybridXcFunc) :: THybridXcFunc_full

  contains

    procedure :: getLrGammaValue => getLrAnalyticalGammaValue
    procedure :: getLrGammaPrimeValue => getdLrAnalyticalGammaValue

    procedure :: getHfGammaValue => getHfAnalyticalGammaValue
    procedure :: getHfGammaPrimeValue => getdHfAnalyticalGammaValue

  end type THybridXcFunc_full


  !> Base type extension for truncated analytical gamma functions.
  type, extends(THybridXcFunc) :: THybridXcFunc_truncated

  contains

    procedure :: getLrGammaValue => getLrTruncatedGammaValue
    procedure :: getLrGammaPrimeValue => getdLrTruncatedGammaValue

    procedure :: getHfGammaValue => getHfTruncatedGammaValue
    procedure :: getHfGammaPrimeValue => getdHfTruncatedGammaValue

  end type THybridXcFunc_truncated


  !> Base type extension for truncated and poly5zero damped analytical gamma functions.
  type, extends(THybridXcFunc) :: THybridXcFunc_truncatedAndDamped

  contains

    procedure :: getLrGammaValue => getLrTruncatedAndDampedGammaValue
    procedure :: getLrGammaPrimeValue => getdLrTruncatedAndDampedGammaValue

    procedure :: getHfGammaValue => getHfTruncatedAndDampedGammaValue
    procedure :: getHfGammaPrimeValue => getdHfTruncatedAndDampedGammaValue

  end type THybridXcFunc_truncatedAndDamped


  !> Base type extension for artificially screened analytical gamma functions.
  type, extends(THybridXcFunc) :: THybridXcFunc_screened

  contains

    procedure :: getLrGammaValue => getLrScreenedGammaValue
    procedure :: getLrGammaPrimeValue => getdLrScreenedGammaValue

    procedure :: getHfGammaValue => getHfScreenedGammaValue
    procedure :: getHfGammaPrimeValue => getdHfScreenedGammaValue

  end type THybridXcFunc_screened


  !> Base type extension for artificially screened and poly5zero damped analytical gamma functions.
  type, extends(THybridXcFunc) :: THybridXcFunc_screenedAndDamped

  contains

    procedure :: getLrGammaValue => getLrScreenedAndDampedGammaValue
    procedure :: getLrGammaPrimeValue => getdLrScreenedAndDampedGammaValue

    procedure :: getHfGammaValue => getHfScreenedAndDampedGammaValue
    procedure :: getHfGammaPrimeValue => getdHfScreenedAndDampedGammaValue

  end type THybridXcFunc_screenedAndDamped


  !> Base type extension for minimum image convention, based on density matrix criteria.
  type, extends(THybridXcFunc) :: THybridXcFunc_mic

  contains

    procedure :: getLrGammaValue => getLrMicGammaValue
    procedure :: getLrGammaPrimeValue => getdLrMicGammaValue

    procedure :: getHfGammaValue => getHfMicGammaValue
    procedure :: getHfGammaPrimeValue => getdHfMicGammaValue

  end type THybridXcFunc_mic


contains

  !> Intitializes the range-separated module.
  subroutine THybridXcFunc_init(this, nAtom, species0, hubbu, screen, omega, camAlpha, camBeta,&
      & tSpin, tREKS, hybridXcAlg, hybridXcType, gammaType, tPeriodic, tRealHS, gammaCutoff,&
      & gSummationCutoff, wignerSeitzReduction, auxiliaryScreening, coeffsDiag, latVecs)

    !> Class instance
    class(THybridXcFunc), intent(out), allocatable :: this

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
    integer, intent(in) :: hybridXcAlg

    !> Hybrid xc-functional type, as extracted from SK-file(s)
    integer, intent(in) :: hybridXcType

    !> Gamma function type (mostly for periodic cases)
    integer, intent(in) :: gammaType

    !> True, if system is periodic (i.e. Gamma-only or k-points)
    logical, intent(in) :: tPeriodic

    !> True, if overlap and Hamiltonian are real-valued
    logical, intent(in) :: tRealHS

    !> Cutoff for truncated Gamma
    real(dp), intent(in), optional :: gammaCutoff

    !> Cutoff for real-space g-summation
    real(dp), intent(in), optional :: gSummationCutoff

    !> Number of unitcells along each supercell folding direction to substract from MIC Wigner-Seitz
    !! cell construction
    integer, intent(in), optional :: wignerSeitzReduction

    !> Auxiliary gamma damping/screening parameter
    real(dp), intent(in), optional :: auxiliaryScreening

    !> Supercell folding coefficients (diagonal elements)
    integer, intent(in), optional :: coeffsDiag(:)

    !> Lattice vectors of (periodic) geometry
    real(dp), intent(in), optional :: latVecs(:,:)

    !! Species indices
    integer :: iSp1, iSp2

    !! Number of unique species in system
    integer :: nUniqueSpecies

    real(dp), allocatable :: wsDistances(:)
    integer, allocatable :: wsDegeneracy(:)

    ! integer :: ii, fd
    ! real(dp) :: dist, gamma, dgamma

    ! Perform basic consistency checks for optional arguments
    if (tPeriodic .and. (.not. present(gSummationCutoff))) then
      call error('Range-separated Module: Periodic systems require g-summation cutoff, which is not&
          & present.')
    end if
    if ((.not. tRealHS) .and. (.not. present(coeffsDiag))) then
      call error('Range-separated Module: General k-point case requires supercell folding&
          & coefficients, which are not present.')
    end if

    ! Allocate selected gamma function types
    select case(gammaType)
    case (hybridXcGammaTypes%full)
      allocate(THybridXcFunc_full:: this)
    case (hybridXcGammaTypes%mic)
      allocate(THybridXcFunc_mic:: this)
    case (hybridXcGammaTypes%truncated)
      if (.not. present(gammaCutoff)) then
        call error('Range-separated Module: Coulomb truncation requires cutoff, which is not&
            & present.')
      end if
      allocate(THybridXcFunc_truncated:: this)
    case (hybridXcGammaTypes%truncatedAndDamped)
      if (.not. present(gammaCutoff)) then
        call error('Range-separated Module: Coulomb truncation requires cutoff, which is not&
            & present.')
      end if
      allocate(THybridXcFunc_truncatedAndDamped:: this)
    case (hybridXcGammaTypes%screened)
      if (.not. present(auxiliaryScreening)) then
        call error('Range-separated Module: Coulomb screening requires additional range-separation&
            & parameter, which is not present.')
      end if
      allocate(THybridXcFunc_screened:: this)
    case (hybridXcGammaTypes%screenedAndDamped)
      if (.not. present(gammaCutoff)) then
        call error('Range-separated Module: Coulomb truncation requires cutoff, which is not&
            & present.')
      end if
      if (.not. present(auxiliaryScreening)) then
        call error('Range-separated Module: Coulomb screening requires additional range-separation&
            & parameter, which is not present.')
      end if
      allocate(THybridXcFunc_screenedAndDamped:: this)
    case default
      call error('Range-separated Module: Invalid gamma function type obtained.')
    end select

    this%gammaType = gammaType

    this%tScreeningInited = .false.
    this%pScreeningThreshold = screen

    this%omega = omega
    this%hybridXcAlg = hybridXcAlg
    this%hybridXcType = hybridXcType
    this%tSpin = tSpin
    this%tREKS = tREKS
    this%camAlpha = camAlpha
    this%camBeta = camBeta
    this%hubbu = hubbu
    this%species0 = species0

    this%camEnergy = 0.0_dp

    if (present(gSummationCutoff)) this%gSummationCutoff = gSummationCutoff
    if (present(wignerSeitzReduction)) this%wignerSeitzReduction = wignerSeitzReduction
    if (present(auxiliaryScreening)) this%auxiliaryScreening = auxiliaryScreening

    if (present(gammaCutoff)) then
      this%gammaCutoff = gammaCutoff

      ! Start beginning of the damping region and 95% of the gamma cutoff.
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
          this%lrGammaAtDamping(iSp1, iSp2) = getLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1),&
              & this%hubbu(iSp2), this%omega, this%gammaDamping)
          this%lrdGammaAtDamping(iSp1, iSp2)&
              & = getdLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2),&
              & this%omega, this%gammaDamping)
          this%lrddGammaAtDamping(iSp1, iSp2)&
              & = getddLrNumericalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2),&
              & this%omega, this%gammaDamping, 1e-08_dp)
        end do
      end do

      if (present(auxiliaryScreening)) then
        allocate(this%lrScreenedGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
        allocate(this%lrdScreenedGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
        allocate(this%lrddScreenedGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
        do iSp2 = 1, nUniqueSpecies
          do iSp1 = 1, nUniqueSpecies
            this%lrScreenedGammaAtDamping(iSp1, iSp2)&
                & = getLrScreenedGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2),&
                & this%omega, this%auxiliaryScreening, this%gammaDamping)
            this%lrdScreenedGammaAtDamping(iSp1, iSp2)&
                & = getdLrScreenedGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2),&
                & this%omega, this%auxiliaryScreening, this%gammaDamping)
            this%lrddScreenedGammaAtDamping(iSp1, iSp2)&
                & = getddLrNumericalScreenedGammaValue_workhorse(this%hubbu(iSp1),&
                & this%hubbu(iSp2), this%omega, this%auxiliaryScreening, this%gammaDamping,&
                & 1e-08_dp)
          end do
        end do
      end if

      call getNumberOfUniqueInt(this%species0, nUniqueSpecies)
      allocate(this%hfGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      allocate(this%hfdGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      allocate(this%hfddGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      do iSp2 = 1, nUniqueSpecies
        do iSp1 = 1, nUniqueSpecies
          this%hfGammaAtDamping(iSp1, iSp2) = getHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1),&
              & this%hubbu(iSp2), this%gammaDamping)
          this%hfdGammaAtDamping(iSp1, iSp2)&
              & = getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2),&
              & this%gammaDamping)
          this%hfddGammaAtDamping(iSp1, iSp2) = getddHfNumericalGammaDeriv(this, iSp1, iSp2,&
              & this%gammaDamping, 1e-08_dp)
        end do
      end do
    end if

    if (this%tREKS .and. this%hybridXcType == hybridXcFunc%hyb) then
      call error("Global hybrid functionals not currently implemented for REKS.")
    end if

    if (this%tREKS .and. this%hybridXcType == hybridXcFunc%cam) then
      call error("General CAM functionals not currently implemented for REKS.")
    end if

    allocate(this%coords(3, nAtom))
    this%coords(:,:) = 0.0_dp
    allocate(this%rCoords(3, nAtom))
    this%rCoords(:,:) = 0.0_dp

    allocate(this%camGammaEval0(nAtom, nAtom))
    this%camGammaEval0(:,:) = 0.0_dp
    allocate(this%camdGammaEval0(nAtom, nAtom, 3))
    this%camdGammaEval0(:,:,:) = 0.0_dp

    ! Check for current restrictions
    if (this%tSpin .and. this%hybridXcAlg == hybridXcAlgo%threshold) then
      call error("Spin-unrestricted calculation for thresholded range separation not yet&
          & implemented!")
    end if

    if (this%tREKS .and. this%hybridXcAlg == hybridXcAlgo%threshold) then
      call error("REKS calculation with thresholded range separation not yet implemented!")
    end if

    if (.not. any([hybridXcAlgo%neighbour, hybridXcAlgo%threshold,&
          & hybridXcAlgo%matrixBased] == this%hybridXcAlg)) then
      call error("Unknown algorithm for screening the exchange in range separation!")
    end if

    if (tPeriodic .and. (.not. present(latVecs))) then
      call error("Range-separated Module: Periodic structure, but no lattice vectors handed over.")
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

      nBvKShifts = product(coeffsDiag)
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

  end subroutine THybridXcFunc_init


  ! !> Brute force search of phase periodicity for a given k-point.
  ! function findPhasePeriodicity(kPoint, uBound) result(periodicity)

  !   !> Single (relative) k-point coordinates for phase factor
  !   real(dp), intent(in) :: kPoint(:)

  !   !> Optional upper bound as limit for brute-force search
  !   integer, intent(in), optional :: uBound

  !   !! G-vector in relative coordinates
  !   real(dp) :: gVec(3)

  !   !! Periodicity of phase factor
  !   integer :: periodicity(3)

  !   !! Upper bound as limit for brute-force search
  !   integer :: uBound_

  !   !! True, if the periodicity has not been found yet
  !   logical :: tSearchIncomplete

  !   !! Phase factor
  !   complex(dp) :: phase

  !   !! Auxiliary variables
  !   integer :: iBound, ii, jj, kk

  !   complex(dp), allocatable :: phaseMesh(:,:,:)
  !   logical, allocatable :: tMatch(:,:,:)

  !   tSearchIncomplete = .true.

  !   if (present(uBound)) then
  !     uBound_ = uBound
  !   else
  !     uBound_ = 10000
  !   end if

  !   do iBound = 1, uBound_
  !     gVec = real(ii, dp)

  !     if (allocated(phaseMesh)) deallocate(phaseMesh)
  !     allocate(phaseMesh(iBound + 1, iBound + 1, iBound + 1))
  !     do kk = 0, iBound
  !       gVec(3) = real(kk, dp)
  !       do jj = 0, iBound
  !         gVec(2) = real(jj, dp)
  !         do ii = 0, iBound
  !           gVec(1) = real(ii, dp)
  !           phaseMesh(ii, jj, kk) = exp(cmplx(0, -1, dp) * dot_product(2.0_dp * pi * kPoint, gVec))
  !         end do
  !       end do
  !     end do

  !     tMatch = abs(phaseMesh - (1.0_dp, 0.0_dp)) > 1e-12_dp
  !     print *, tMatch
  !     stop
  !   end do

  !   periodicity(:) = 0

  ! end function findPhasePeriodicity1D


  !> Checks if obtained supercell folding matrix meets current requirements.
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
  pure function THybridXcFunc_foldToBvK(this, vector) result(bvKShift)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Vector (in relative coordinates) to fold back to BvK cell
    real(dp), intent(in) :: vector(:)

    !> Corresponding BvK vector
    integer :: bvKShift(3)

    bvKShift(:) = modulo(nint(vector), this%coeffsDiag)

  end function THybridXcFunc_foldToBvK


  !> Folds relative real-space vector back to BvK region and returns indices of density matrix.
  pure function THybridXcFunc_foldToBvKIndex(this, vector) result(bvKIndex)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Vector (in relative coordinates) to fold back to BvK cell
    real(dp), intent(in) :: vector(:)

    !> Corresponding BvK indexing
    integer :: bvKIndex(3)

    ! additionally shift by 1, so that indices start at 1 and not at 0
    bvKIndex(:) = modulo(nint(vector), this%coeffsDiag) + 1

  end function THybridXcFunc_foldToBvKIndex


  !> Updates the range-separated module on coordinate change (non-periodic version).
  subroutine updateCoords_cluster(this, rCoords)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

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
        this%camGammaEval0(iAtom1, iAtom2) = getCamAnalyticalGammaValue_workhorse(this%hubbu(iSp1),&
            & this%hubbu(iSp2), this%omega, this%camAlpha, this%camBeta, dist)
        this%camGammaEval0(iAtom2, iAtom1) = this%camGammaEval0(iAtom1, iAtom2)
      end do
    end do

    if (this%tScreeningInited) then
      this%hPrev(:,:) = 0.0_dp
      this%dRhoPrev(:,:) = 0.0_dp
      this%camEnergy = 0.0_dp
    end if

  end subroutine updateCoords_cluster


  !> Updates the range-separated module on coordinate change (Gamma-only version).
  subroutine updateCoords_gamma(this, env, symNeighbourList, nNeighbourCamSym, skOverCont, orb,&
      & latVecs, recVecs2p, iSquare)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

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

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descN(descLen), descB(descLen)

    !! Atom blocks from sparse, real-space overlap matrices S_{\alpha\mu}, S_{\beta\nu}
    real(dp), pointer :: pSbn(:,:)

    !! Neighbour indices (+corresponding atom indices)
    integer :: iNeighN, iAtB

    !! Folded (to central cell) atom index
    integer :: iAtBfold

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Indices of interacting atoms in central cell, as well as their global species index
    integer :: iAtM, iAtN, iSpM, iSpN

    !! Auxiliary variables for setting up 2D pointer to sparse overlap
    integer :: ind, nOrbAt, nOrbNeigh

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
    call getCellTranslations(this%cellVecsG, this%rCellVecsG, latVecs, recVecs2p,&
        & this%gSummationCutoff)

    ! build symmetric, sparse overlap
    call buildS(env, this%overSym, skOverCont, this%rCoords, nNeighbourCamSym,&
        & symNeighbourList%neighbourList%iNeighbour, symNeighbourList%species,&
        & symNeighbourList%iPair, orb)

    if (this%tScreeningInited) then
      this%hPrev(:,:) = 0.0_dp
      this%dRhoPrev(:,:) = 0.0_dp
      this%camEnergy = 0.0_dp
    end if

    ! ##################### Beginning of \gamma(\bm{g}) pre-tabulation #############################

    ! allocate max estimates of square overlap blocks and index array for sorting
    if (allocated(this%testSquareOver)) deallocate(this%testSquareOver)
    allocate(this%testSquareOver(nAtom0))
    if (allocated(this%overlapIndices)) deallocate(this%overlapIndices)
    allocate(this%overlapIndices(nAtom0))

    do iAtN = 1, nAtom0
      descN = getDescriptor(iAtN, iSquare)
      allocate(this%testSquareOver(iAtN)%array(nNeighbourCamSym(iAtN) + 1))
      do iNeighN = 0, nNeighbourCamSym(iAtN)
        iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighN, iAtN)
        iAtBfold = symNeighbourList%img2CentCell(iAtB)
        descB = getDescriptor(iAtBfold, iSquare)
        ! get 2D pointer to Sbn overlap block
        ind = symNeighbourList%iPair(iNeighN, iAtN) + 1
        nOrbAt = descN(iNOrb)
        nOrbNeigh = descB(iNOrb)
        pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
        this%testSquareOver(iAtN)%array(iNeighN + 1) = maxval(abs(pSbn))
      end do
    end do

    ! sort max square overlap estimates (descending)
    ! this way we can exit the whole loop a.s.a. threshold has been undershot for the first time
    do iAtN = 1, nAtom0
      allocate(this%overlapIndices(iAtN)%array(nNeighbourCamSym(iAtN) + 1))
      call index_heap_sort(this%overlapIndices(iAtN)%array, this%testSquareOver(iAtN)%array)
      ! switch from ascending to descending
      this%overlapIndices(iAtN)%array(:)&
          & = this%overlapIndices(iAtN)%array(size(this%overlapIndices(iAtN)%array):1:-1)
    end do

    do iAtM = 1, nAtom0
      iSpM = this%species0(iAtM)
      do iAtN = 1, nAtom0
        iSpN = this%species0(iAtN)
        this%camGammaEval0(iAtM, iAtN) = getCamGammaGSum(this, iAtM, iAtN, iSpM, iSpN,&
            & this%rCellVecsG)
        this%camGammaEval0(iAtN, iAtM) = getCamGammaGSum(this, iAtN, iAtM, iSpN, iSpM,&
            & this%rCellVecsG)
        this%camdGammaEval0(iAtM, iAtN, :) = getCamGammaPrimeGSum(this, iAtM, iAtN, iSpM, iSpN,&
            & this%rCellVecsG)
        this%camdGammaEval0(iAtN, iAtM, :) = getCamGammaPrimeGSum(this, iAtN, iAtM, iSpN, iSpM,&
            & this%rCellVecsG)
      end do
    end do

  end subroutine updateCoords_gamma


  !> Updates the range-separated module on coordinate change (k-point version).
  subroutine updateCoords_kpts(this, env, symNeighbourList, nNeighbourCamSym, skOverCont, orb,&
      & latVecs, recVecs2p, cellVecs, iSquare)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

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

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVecs(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descN(descLen), descB(descLen)

    !! Index array for descending sorting of Gamma arrays
    integer, allocatable :: gammasortIdx(:)

    !! Temporary storage for g-resolved gamma values (+ directed derivatives)
    real(dp), allocatable :: gammaEvalGTmp(:), dGammaEvalGTmp(:,:)

    !! Atom blocks from sparse, real-space overlap matrices S_{\alpha\mu}, S_{\beta\nu}
    real(dp), pointer :: pSbn(:,:)

    !! Neighbour indices (+corresponding atom indices)
    integer :: iNeighN, iNeighNsort, iAtB, iNeighM, iNeighMsort, iAtA

    !! Auxiliary variables for setting up 2D pointer to sparse overlap
    integer :: ind, nOrbAt, nOrbNeigh

    !! Folded (to central cell) atom indices
    integer :: iAtBfold

    !! Number of non-zero entries of (sorted) array
    integer :: nNonZeroEntries

    !! Indices of interacting atoms in central cell, as well as their global species index
    integer :: iAtM, iAtN, iSpM, iSpN

    !! Number of atoms in the central cell
    integer :: nAtom0

    !! Number of g-summation summands
    integer :: nGShifts

    !! Composite index iAtM/iAtN
    integer :: ii, iAtMN(2, size(this%species0)**2)

    !! Dummy array with zeros
    real(dp) :: zeros(3)

    !! Minimum and maximum cell translations of neighbour list in each direction
    integer :: minVecX, maxVecX, minVecY, maxVecY, minVecZ, maxVecZ

    !! Real-space shift components of neighbour list atoms
    integer :: vecLx, vecLy, vecLz, vecHx, vecHy, vecHz

    !! Real-space \vec{h} and \vec{l} vectors in relative coordinates
    integer :: vecH(3), vecL(3)

    !! Real-space \vec{h} and \vec{l} vectors in absolute coordinates
    real(dp) :: rVecH(3), rVecL(3)

    !! Real-space shift of delta density matrix argument
    integer :: dPabArg(3)

    !! Iterates over g-vectors
    integer :: iG

    ! Max estimate for overlap and products
    real(dp) :: maxEstimate, pSbnMax

    zeros(:) = 0.0_dp

    nAtom0 = size(this%species0)

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

    ! build symmetric, sparse overlap
    call buildS(env, this%overSym, skOverCont, this%rCoords, nNeighbourCamSym,&
        & symNeighbourList%neighbourList%iNeighbour, symNeighbourList%species,&
        & symNeighbourList%iPair, orb)

    if (this%tScreeningInited) then
      this%hPrevCplxHS(:,:,:) = 0.0_dp
      this%dRhoPrevCplxHS(:,:,:,:,:,:) = 0.0_dp
      this%camEnergy = 0.0_dp
    end if

    ! ##################### Beginning of overlap screening pre-tabulation ##########################

    ! allocate max estimates of square overlap blocks and index array for sorting
    if (allocated(this%testSquareOver)) deallocate(this%testSquareOver)
    allocate(this%testSquareOver(nAtom0))
    if (allocated(this%overlapIndices)) deallocate(this%overlapIndices)
    allocate(this%overlapIndices(nAtom0))

    do iAtN = 1, nAtom0
      descN = getDescriptor(iAtN, iSquare)
      allocate(this%testSquareOver(iAtN)%array(nNeighbourCamSym(iAtN) + 1))
      do iNeighN = 0, nNeighbourCamSym(iAtN)
        iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighN, iAtN)
        iAtBfold = symNeighbourList%img2CentCell(iAtB)
        descB = getDescriptor(iAtBfold, iSquare)
        ! get 2D pointer to Sbn overlap block
        ind = symNeighbourList%iPair(iNeighN, iAtN) + 1
        nOrbAt = descN(iNOrb)
        nOrbNeigh = descB(iNOrb)
        pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
        this%testSquareOver(iAtN)%array(iNeighN + 1) = maxval(abs(pSbn))
      end do
    end do

    ! sort max square overlap estimates (descending)
    ! this way we can exit the whole loop a.s.a. threshold has been undershot for the first time
    do iAtN = 1, nAtom0
      allocate(this%overlapIndices(iAtN)%array(nNeighbourCamSym(iAtN) + 1))
      call index_heap_sort(this%overlapIndices(iAtN)%array, this%testSquareOver(iAtN)%array)
      ! switch from ascending to descending
      this%overlapIndices(iAtN)%array(:)&
          & = this%overlapIndices(iAtN)%array(size(this%overlapIndices(iAtN)%array):1:-1)
    end do

    ! ##################### Beginning of \gamma(\bm{g}) pre-tabulation #############################

    if (this%gammaType == hybridXcGammaTypes%mic) then
      if (allocated(this%wsVectors)) deallocate(this%wsVectors)
      ! Generate "save" Wigner-Seitz vectors for density matrix arguments
      call generateWignerSeitzGrid(max(this%coeffsDiag - this%wignerSeitzReduction, 1), latVecs,&
          & this%wsVectors, tExcludeSurface=.false.)

      ! The Wigner-Seitz grid actually defines all the relevant g-vectors, therefore copy over
      this%cellVecsG = this%wsVectors
      if (allocated(this%rCellVecsG)) deallocate(this%rCellVecsG)
      allocate(this%rCellVecsG(3, size(this%cellVecsG, dim=2)))
      do iG = 1, size(this%rCellVecsG, dim=2)
        this%rCellVecsG(:, iG) = matmul(latVecs, this%cellVecsG(:, iG))
      end do

    else

      ! Build up composite index iAtMN for collapsing iAtM and iAtN
      ind = 1
      loopM: do iAtM = 1, nAtom0
        loopN: do iAtN = 1, nAtom0
          iAtMN(1, ind) = iAtM
          iAtMN(2, ind) = iAtN
          ind = ind + 1
        end do loopN
      end do loopM

      ! get all cell translations within given cutoff
      call getCellTranslations(this%cellVecsG, this%rCellVecsG, latVecs, recVecs2p,&
          & this%gSummationCutoff)

      if (allocated(this%nNonZeroGammaG)) deallocate(this%nNonZeroGammaG)
      allocate(this%nNonZeroGammaG(nAtom0, nAtom0))

      if (allocated(this%camGammaEvalG)) deallocate(this%camGammaEvalG)
      allocate(this%camGammaEvalG(nAtom0, nAtom0))
      if (allocated(this%camdGammaEvalG)) deallocate(this%camdGammaEvalG)
      allocate(this%camdGammaEvalG(nAtom0, nAtom0))

      nGShifts = size(this%cellVecsG, dim=2)
      allocate(gammaEvalGTmp(nGShifts))
      allocate(dGammaEvalGTmp(3, nGShifts))
      allocate(gammaSortIdx(nGShifts))

      ! Pre-tabulate CAM gamma integrals (+ directed derivatives)
      do ii = 1, nAtom0**2
        iAtM = iAtMN(1, ii)
        iAtN = iAtMN(2, ii)
        iSpM = this%species0(iAtM)
        iSpN = this%species0(iAtN)
        ! pre-tabulate g-resolved \gamma_{\mu\nu}(\vec{g})
        gammaEvalGTmp(:) = getCamGammaGResolved(this, iAtM, iAtN, iSpM, iSpN, this%rCellVecsG,&
            & zeros)
        call index_heap_sort(gammaSortIdx, gammaEvalGTmp)
        gammaSortIdx(:) = gammaSortIdx(size(gammaSortIdx):1:-1)
        gammaEvalGTmp(:) = gammaEvalGTmp(gammaSortIdx)
        nNonZeroEntries = getNumberOfNonZeroElements(gammaEvalGTmp)
        this%nNonZeroGammaG(iAtM, iAtN)%array = gammaSortIdx(1:nNonZeroEntries)
        this%camGammaEvalG(iAtM, iAtN)%array = gammaEvalGTmp(1:nNonZeroEntries)
        ! pre-tabulate g-resolved \partial\gamma_{\mu\nu}(\vec{g})
        dGammaEvalGTmp(:,:) = getCamGammaPrimeGResolved(this, iAtM, iAtN, iSpM, iSpN,&
            & this%rCellVecsG(:, gammaSortIdx))
        this%camdGammaEvalG(iAtM, iAtN)%array = dGammaEvalGTmp(:, 1:nNonZeroEntries)
      end do

    end if


  contains

    !> Returns the number of non-zero elements in a descending array of reals.
    function getNumberOfNonZeroElements(array) result(nNonZeroEntries)

      !> Descending one-dimensional, real-valued array to search
      real(dp), intent(in) :: array(:)

      !> Number of non-zero entries
      real(dp) :: nNonZeroEntries

      !! iterates over all array elements
      integer :: ii

      nNonZeroEntries = 0

      do ii = 1, size(array)
        if (array(ii) < 1e-16_dp) return
        nNonZeroEntries = ii
      end do

    end function getNumberOfNonZeroElements

  end subroutine updateCoords_kpts


  !> Builds MPI rank dependent composite index for nested HFX loops.
  subroutine getFourLoopCompositeIndex(this, env, nNeighbourCamSym, pMax, compositeIndex)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !! Max estimate for delta-delta density matrix
    real(dp), intent(in) :: pMax

    !> Composite index for four nested loops
    integer, intent(out), allocatable :: compositeIndex(:,:)

    !! Auxiliary variables
    integer :: ii, indGlobal, indLocal

    !! Atom indices (central cell)
    integer :: iAtM, iAtN

    !! Neighbour indices
    integer :: iNeighN, iNeighM

    !! Sorted (according to max overlap estimates) neighbour indices
    integer :: iNeighMsort, iNeighNsort

    !! Max estimate for products of delta-delta density matrix and max overlap estimates
    real(dp) :: maxEstimate, pSbnPabMax

    !! Start and end index for MPI parallelization, if applicable
    integer :: iParallelStart, iParallelEnd

    !! Loop counter
    integer :: ind

    !! Number of atoms in central cell
    integer :: nAtom0

    nAtom0 = size(this%species0)

    ! Lead process counts composite index entries for MPI parallelization
    ! (includes integral screening)
    if (env%tGlobalLead) then
      ! First run over loops to count non-vanishing cycles
      ind = 0
      loopM1: do iAtM = 1, nAtom0
        loopN1: do iAtN = 1, nAtom0
          loopB1: do iNeighN = 0, nNeighbourCamSym(iAtN)
            iNeighNsort = this%overlapIndices(iAtN)%array(iNeighN + 1) - 1
            pSbnPabMax = pMax * this%testSquareOver(iAtN)%array(iNeighNsort + 1)
            if (pSbnPabMax < this%pScreeningThreshold) exit loopB1
            loopA1: do iNeighM = 0, nNeighbourCamSym(iAtM)
              iNeighMsort = this%overlapIndices(iAtM)%array(iNeighM + 1) - 1
              maxEstimate = pSbnPabMax * this%testSquareOver(iAtM)%array(iNeighMsort + 1)
              if (maxEstimate < this%pScreeningThreshold) exit loopA1

              ind = ind + 1

            end do loopA1
          end do loopB1
        end do loopN1
      end do loopM1
    end if

  #:if WITH_MPI
    call mpifx_bcast(env%mpi%globalComm, ind)
    call getStartAndEndIndex(ind, env%mpi%globalComm%size, env%mpi%globalComm%rank, iParallelStart,&
        & iParallelEnd)
  #:else
    iParallelStart = 1
    iParallelEnd = ind
  #:endif

    ! Composite index needs to be allocated on all ranks
    allocate(compositeIndex(4, iParallelEnd - iParallelStart + 1))

    ! Flatten four nested loops into a single composite index for MPI parallelization
    indGlobal = 1
    indLocal = 1
    loopM2: do iAtM = 1, nAtom0
      loopN2: do iAtN = 1, nAtom0
        loopB2: do iNeighN = 0, nNeighbourCamSym(iAtN)
          iNeighNsort = this%overlapIndices(iAtN)%array(iNeighN + 1) - 1
          pSbnPabMax = pMax * this%testSquareOver(iAtN)%array(iNeighNsort + 1)
          if (pSbnPabMax < this%pScreeningThreshold) exit loopB2
          loopA2: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iNeighMsort = this%overlapIndices(iAtM)%array(iNeighM + 1) - 1
            maxEstimate = pSbnPabMax * this%testSquareOver(iAtM)%array(iNeighMsort + 1)
            if (maxEstimate < this%pScreeningThreshold) exit loopA2
            if (indGlobal >= iParallelStart .and. indGlobal <= iParallelEnd) then
              compositeIndex(1, indLocal) = iAtM
              compositeIndex(2, indLocal) = iAtN
              compositeIndex(3, indLocal) = iNeighNsort
              compositeIndex(4, indLocal) = iNeighMsort
              indLocal = indLocal + 1
            end if
            indGlobal = indGlobal + 1
          end do loopA2
        end do loopB2
      end do loopN2
    end do loopM2

  end subroutine getFourLoopCompositeIndex


  !> Interface routine for adding CAM range-separated contributions to the Hamiltonian.
  !! (non-periodic version)
  subroutine addCamHamiltonian_cluster(this, env, densSqr, over, iNeighbour, nNeighbourCam,&
      & iSquare, iPair, orb, HH, overlap)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

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

    call env%globalTimer%startTimer(globalTimers%hybridXcH)

    select case(this%hybridXcAlg)
    case (hybridXcAlgo%threshold)
      call addCamHamiltonianThreshold_cluster(this, overlap, densSqr, iNeighbour, nNeighbourCam,&
          & iSquare, HH, orb)
    case (hybridXcAlgo%neighbour)
      call addCamHamiltonianNeighbour_cluster(this, densSqr, over, iNeighbour, nNeighbourCam,&
          & iSquare, iPair, orb, HH)
    case (hybridXcAlgo%matrixBased)
      call addCamHamiltonianMatrix_cluster_real(this, iSquare, overlap, densSqr, HH)
    end select

    call env%globalTimer%stopTimer(globalTimers%hybridXcH)

  end subroutine addCamHamiltonian_cluster


  !> Interface routine for adding CAM range-separated contributions to the Hamiltonian.
  !! (Gamma-only version)
  subroutine addCamHamiltonian_gamma(this, env, densSqr, overSqr, symNeighbourList,&
      & nNeighbourCamSym, iSquare, orb, iKS, nKS, HH)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

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

    call env%globalTimer%startTimer(globalTimers%hybridXcH)

    select case(this%hybridXcAlg)
    case (hybridXcAlgo%threshold)
      call error('Thresholded algorithm not implemented for periodic systems.')
    case (hybridXcAlgo%neighbour)
      call addCamHamiltonianNeighbour_gamma(this, env, densSqr, symNeighbourList, nNeighbourCamSym,&
          & iSquare, orb, iKS, nKS, HH)
    case (hybridXcAlgo%matrixBased)
      call addCamHamiltonianMatrix_cluster_real(this, iSquare, overSqr, densSqr, HH)
    end select

    call env%globalTimer%stopTimer(globalTimers%hybridXcH)

  end subroutine addCamHamiltonian_gamma


  !> Interface routine for adding CAM range-separated contributions to the Hamiltonian.
  !! (k-point version)
  subroutine addCamHamiltonian_kpts(this, env, deltaRhoSqr, symNeighbourList, nNeighbourCamSym,&
      & rCellVecs, cellVecs, iSquare, orb, kPoints, kWeights, iKiSToiGlobalKS, HSqr)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) delta spin-density matrix at BvK real-space shifts
    real(dp), intent(in) :: deltaRhoSqr(:,:,:,:,:,:)

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Vectors to neighboring unit cells in absolute units
    real(dp), intent(in) :: rCellVecs(:,:)

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVecs(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> K-points in relative coordinates to calculate delta H(k) for
    real(dp), intent(in) :: kPoints(:,:)

    !> K-point weights (for energy contribution)
    real(dp), intent(in) :: kWeights(:)

    !> Composite index for mapping iK/iS --> iGlobalKS for arrays present at every MPI rank
    integer, intent(in) :: iKiSToiGlobalKS(:,:)

    !> Square (unpacked) Hamiltonian for all k-point/spin composite indices to be updated
    complex(dp), intent(inout) :: HSqr(:,:,:)

    call env%globalTimer%startTimer(globalTimers%hybridXcH)

    select case(this%hybridXcAlg)
    case (hybridXcAlgo%threshold)
      call error('Thresholded algorithm not implemented for periodic systems.')
    case (hybridXcAlgo%neighbour)
      if (this%gammaType == hybridXcGammaTypes%mic) then
        call addCamHamiltonianNeighbour_kpts_mic(this, env, deltaRhoSqr, symNeighbourList,&
            & nNeighbourCamSym, rCellVecs, cellVecs, iSquare, orb, kPoints, kWeights,&
            & iKiSToiGlobalKS, HSqr)
    else
      call addCamHamiltonianNeighbour_kpts_ct(this, env, deltaRhoSqr, symNeighbourList,&
          & nNeighbourCamSym, cellVecs, iSquare, orb, kPoints, kWeights, iKiSToiGlobalKS, HSqr)
      end if
    case (hybridXcAlgo%matrixBased)
      call error('Matrix based algorithm not implemented for periodic systems.')
    end select

    call env%globalTimer%stopTimer(globalTimers%hybridXcH)

  end subroutine addCamHamiltonian_kpts


  !> Adds CAM range-separated contributions to Hamiltonian, using the thresholding algorithm.
  !! (non-periodic version)
  subroutine addCamHamiltonianThreshold_cluster(this, overlap, deltaRho, iNeighbour, nNeighbourCam,&
      & iSquare, hamiltonian, orb)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

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
    hamiltonian(:,:) = hamiltonian + this%hprev
    this%camEnergy = this%camEnergy + evaluateEnergy_real(this%hprev, tmpDRho)

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
            gammabatchtmp = this%camGammaEval0(iAtMu, iAtNu) + this%camGammaEval0(iAt1, iAtNu)
            loopLL: do ll = 1, nAtom0
              iAt2 = ovrInd(iAtNu, nAtom0 + 1 - ll)
              iSp2 = this%species0(iAt2)
              nOrb2 = orb%nOrbSpecies(iSp2)
              tstbound = prb * testOvr(iAt2, iAtNu)
              if (abs(tstbound) < this%pScreeningThreshold) then
                exit loopLL
              end if
              desc2 = getDescriptor(iAt2, iSquare)
              gammabatch = (this%camGammaEval0(iAtMu, iAt2) + this%camGammaEval0(iAt1, iAt2)&
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

      !> Class instance
      class(THybridXcFunc), intent(inout) :: this

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

  end subroutine addCamHamiltonianThreshold_cluster


  !> Adds CAM range-separated contributions to Hamiltonian, using neighbour-list based algorithm.
  !! (non-periodic version)
  subroutine addCamHamiltonianNeighbour_cluster(this, densSqr, over, iNeighbour, nNeighbourCam,&
      & iSquare, iPair, orb, HH)

    !> instance of object
    class(THybridXcFunc), intent(inout) :: this

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

    HH(:,:) = HH + tmpHH
    this%camEnergy = this%camEnergy + evaluateEnergy_real(tmpHH, tmpDRho)

  contains

    !> Allocates storage for mapping 1D<->2D array sections.
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


    !> Actually evaluates the neighbour based cut-off hamiltonian.
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
            gamma1 = this%camGammaEval0(iAtA, iAtN) + this%camGammaEval0(iAtA, iAtB)
            loopM: do iNeighA = 0, nNeighbourCam(iAtA)
              iAtM = iNeighbour(iNeighA, iAtA)
              descM = getDescriptor(iAtM, iSquare)
              call copyOverlapBlock(iAtA, iNeighA, descA(iNOrb), descM(iNOrb), Sma, pSma)
              call transposeBlock(pSma, Sam, pSam)
              gamma2 = this%camGammaEval0(iAtM, iAtN) + this%camGammaEval0(iAtM, iAtB)
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


    !> Copies atom block from sparse matrix.
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


    !> Copies density matrix block from sparse matrix.
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


    !> Transposes a block.
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


    !> Adds a contribution to a Hamiltonian block.
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

  end subroutine addCamHamiltonianNeighbour_cluster


  !> Update Hamiltonian with CAM range-separated contributions, using matrix-based algorithm.
  !! (non-periodic, real version)
  !!
  !! The routine provides a matrix-matrix multiplication based implementation of the 3rd term in
  !! Eq. 26 in https://doi.org/10.1063/1.4935095
  subroutine addCamHamiltonianMatrix_cluster_real(this, iSquare, overlap, densSqr, HH)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

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
    real(dp), allocatable :: camGammaAO(:,:)
    real(dp), allocatable :: Hcam(:,:)

    integer :: nOrb

    nOrb = size(overlap,dim=1)

    allocate(Smat(nOrb, nOrb))
    allocate(Dmat(nOrb, nOrb))
    allocate(camGammaAO(nOrb, nOrb))
    allocate(Hcam(nOrb, nOrb))

    call allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, camGammaAO)
    call evaluateHamiltonian(this, Smat, Dmat, camGammaAO, Hcam)
    HH(:,:) = HH + Hcam
    this%camEnergy = this%camEnergy + 0.5_dp * sum(Dmat * Hcam)

  contains

    !> Set up storage and get orbital-by-orbital gamma matrix
    subroutine allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, camGammaAO)

      !> Class instance
      class(THybridXcFunc), intent(inout) :: this

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

      !> Symmetrized CAM gamma matrix
      real(dp), intent(out) :: camGammaAO(:,:)

      integer :: nAtom, iAt, jAt

      nAtom = size(this%camGammaEval0, dim=1)

      ! Symmetrize Hamiltonian, overlap, density matrices
      call symmetrizeHS(HH)
      Smat(:,:) = overlap
      call symmetrizeHS(Smat)
      Dmat(:,:) = densSqr
      call symmetrizeHS(Dmat)

      ! Get CAM gamma variable
      camGammaAO(:,:) = 0.0_dp
      do iAt = 1, nAtom
        do jAt = 1, nAtom
          camGammaAO(iSquare(jAt):iSquare(jAt+1)-1, iSquare(iAt):iSquare(iAt+1)-1)&
              & = this%camGammaEval0(jAt, iAt)
        end do
      end do

    end subroutine allocateAndInit


    !> Evaluates the Hamiltonian, using GEMM operations.
    subroutine evaluateHamiltonian(this, Smat, Dmat, camGammaAO, Hcam)

      !> Class instance
      class(THybridXcFunc), intent(inout) :: this

      !> Symmetrized square overlap matrix
      real(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized CAM gamma matrix
      real(dp), intent(in) :: camGammaAO(:,:)

      !> Symmetrized CAM Hamiltonian matrix
      real(dp), intent(out) :: Hcam(:,:)

      real(dp), allocatable :: Hmat(:,:)
      real(dp), allocatable :: tmpMat(:,:)

      integer :: nOrb

      nOrb = size(Smat, dim=1)

      allocate(Hmat(nOrb, nOrb))
      allocate(tmpMat(nOrb, nOrb))

      Hcam(:,:) = 0.0_dp

      call gemm(tmpMat, Smat, Dmat)
      call gemm(Hcam, tmpMat, Smat)
      Hcam(:,:) = Hcam * camGammaAO

      tmpMat(:,:) = tmpMat * camGammaAO
      call gemm(Hcam, tmpMat, Smat, alpha=1.0_dp, beta=1.0_dp)

      Hmat(:,:) = Dmat * camGammaAO
      call gemm(tmpMat, Smat, Hmat)
      call gemm(Hcam, tmpMat, Smat, alpha=1.0_dp, beta=1.0_dp)

      call gemm(tmpMat, Dmat, Smat)
      tmpMat(:,:) = tmpMat * camGammaAO
      call gemm(Hcam, Smat, tmpMat, alpha=1.0_dp, beta=1.0_dp)

      if (this%tSpin .or. this%tREKS) then
        Hcam(:,:) = -0.25_dp * Hcam
      else
        Hcam(:,:) = -0.125_dp * Hcam
      end if

    end subroutine evaluateHamiltonian

  end subroutine addCamHamiltonianMatrix_cluster_real


  !> Update Hamiltonian with CAM range-separated contributions, using matrix-based algorithm.
  !! (non-periodic, complex version)
  !!
  !! The routine provides a matrix-matrix multiplication based implementation of the 3rd term in
  !! Eq. 26 in https://doi.org/10.1063/1.4935095
  subroutine addCamHamiltonianMatrix_cluster_cmplx(this, iSquare, overlap, densSqr, HH)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Square (unpacked) overlap matrix.
    complex(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    complex(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated.
    complex(dp), intent(inout) :: HH(:,:)

    complex(dp), allocatable :: Smat(:,:)
    complex(dp), allocatable :: Dmat(:,:)
    real(dp), allocatable :: camGammaAO(:,:)
    complex(dp), allocatable :: gammaCmplx(:,:)
    complex(dp), allocatable :: Hcam(:,:)

    integer :: nOrb

    nOrb = size(overlap,dim=1)

    allocate(Smat(nOrb, nOrb))
    allocate(Dmat(nOrb, nOrb))
    allocate(camGammaAO(nOrb, nOrb))
    allocate(gammaCmplx(nOrb, nOrb))
    allocate(Hcam(nOrb, nOrb))

    call allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, camGammaAO, gammaCmplx)

    call evaluateHamiltonian(this, Smat, Dmat, gammaCmplx, Hcam)

    HH(:,:) = HH + Hcam

    this%camEnergy = this%camEnergy + 0.5_dp * real(sum(Dmat * Hcam), dp)

  contains

    subroutine allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, camGammaAO,&
        & gammaCmplx)

      !> instance
      class(THybridXcFunc), intent(inout) :: this

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

      !> Symmetrized CAM gamma matrix
      real(dp), intent(out) :: camGammaAO(:,:)

      !> Symmetrized CAM gamma matrix
      complex(dp), intent(out) :: gammaCmplx(:,:)

      integer :: nAtom, iAt, jAt

      nAtom = size(this%camGammaEval0, dim=1)

      !! Symmetrize Hamiltonian, overlap, density matrices
      call hermitianSquareMatrix(HH)
      Smat(:,:) = overlap
      call hermitianSquareMatrix(Smat)
      Dmat(:,:) = densSqr
      call hermitianSquareMatrix(Dmat)

      ! Get CAM gamma variable
      camGammaAO(:,:) = 0.0_dp
      do iAt = 1, nAtom
        do jAt = 1, nAtom
          camGammaAO(iSquare(jAt):iSquare(jAt+1)-1,iSquare(iAt):iSquare(iAt+1)-1) =&
              & this%camGammaEval0(jAt, iAt)
        end do
      end do
      gammaCmplx(:,:) = camGammaAO

    end subroutine allocateAndInit


    subroutine evaluateHamiltonian(this, Smat, Dmat, gammaCmplx, Hcam)

      !> instance
      class(THybridXcFunc), intent(inout) :: this

      !> Symmetrized square overlap matrix
      complex(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      complex(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized CAM gamma matrix
      complex(dp), intent(in) :: gammaCmplx(:,:)

      !> Symmetrized CAM Hamiltonian matrix
      complex(dp), intent(out) :: Hcam(:,:)

      complex(dp), allocatable :: Hmat(:,:)
      complex(dp), allocatable :: tmpMat(:,:)

      integer :: nOrb

      nOrb = size(Smat,dim=1)

      allocate(Hmat(nOrb,nOrb))
      allocate(tmpMat(nOrb,nOrb))

      Hcam(:,:) = cmplx(0.0_dp,0.0_dp,dp)

      call gemm(tmpMat, Smat, Dmat)
      call gemm(Hcam, tmpMat, Smat)
      Hcam(:,:) = Hcam * gammaCmplx

      tmpMat(:,:) = tmpMat * gammaCmplx
      call gemm(Hcam, tmpMat, Smat, alpha=(1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))

      Hmat(:,:) = Dmat * gammaCmplx
      call gemm(tmpMat, Smat, Hmat)
      call gemm(Hcam, tmpMat, Smat, alpha=(1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))

      call gemm(tmpMat, Dmat, Smat)
      tmpMat(:,:) = tmpMat * gammaCmplx
      call gemm(Hcam, Smat, tmpMat, alpha=(1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))

      if (this%tSpin) then
        Hcam(:,:) = -0.25_dp * Hcam
      else
        Hcam(:,:) = -0.125_dp * Hcam
      end if

    end subroutine evaluateHamiltonian

  end subroutine addCamHamiltonianMatrix_cluster_cmplx


  !> Adds CAM range-separated contributions to Hamiltonian, using neighbour-list based algorithm.
  !! (Gamma-only version)
  subroutine addCamHamiltonianNeighbour_gamma(this, env, deltaRhoSqr, symNeighbourList,&
      & nNeighbourCamSym, iSquare, orb, iKS, nKS, HSqr)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

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

    !! Density matrix block from sparse matrix
    real(dp), allocatable :: Pab(:,:)

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descA(descLen), descB(descLen), descM(descLen), descN(descLen)

    !! Temporary storages
    real(dp), allocatable :: tmpDeltaRhoSqr(:,:), tmpDeltaDeltaRhoSqr(:,:), tmpHSqr(:,:)

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
    real(dp) :: pMax, maxEstimate, pSbnMax

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

    loopMN: do ii = iParallelStart, iParallelEnd
      iAtM = iAtMN(1, ii)
      iAtN = iAtMN(2, ii)
      descM = getDescriptor(iAtM, iSquare)
      descN = getDescriptor(iAtN, iSquare)
      loopB: do iNeighN = 0, nNeighbourCamSym(iAtN)
        iNeighNsort = this%overlapIndices(iAtN)%array(iNeighN + 1) - 1
        pSbnMax = this%testSquareOver(iAtN)%array(iNeighNsort + 1)
        if (pSbnMax < this%pScreeningThreshold) exit loopB
        iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighNsort, iAtN)
        iAtBfold = symNeighbourList%img2CentCell(iAtB)
        descB = getDescriptor(iAtBfold, iSquare)
        ! get 2D pointer to S_{\beta\nu}(\vec{l}) overlap block
        ind = symNeighbourList%iPair(iNeighNsort, iAtN) + 1
        nOrbAt = descN(iNOrb)
        nOrbNeigh = descB(iNOrb)
        pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
        ! \tilde{\gamma}_{\mu\nu} + \tilde{\gamma}_{\mu\beta}
        gammaMNMB = this%camGammaEval0(iAtM, iAtN) + this%camGammaEval0(iAtM, iAtBfold)
        loopA: do iNeighM = 0, nNeighbourCamSym(iAtM)
          iNeighMsort = this%overlapIndices(iAtM)%array(iNeighM + 1) - 1
          maxEstimate = pSbnMax * this%testSquareOver(iAtM)%array(iNeighMsort + 1)
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

          gammaTot = gammaMNMB + this%camGammaEval0(iAtAfold, iAtN)&
              & + this%camGammaEval0(iAtAfold, iAtBfold)

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
    HSqr(:,:) = HSqr + this%hprev

    ! HSqr(:,:) = HSqr + tmpHSqr

    ! Add energy contribution but divide by the number of processes
  #:if WITH_MPI
    this%camEnergy = this%camEnergy + evaluateEnergy_real(this%hprev, tmpDeltaRhoSqr)&
        & / real(env%mpi%groupComm%size, dp)
    ! this%camEnergy = this%camEnergy + evaluateEnergy_real(tmpHSqr, tmpDeltaRhoSqr)&
    !     & / real(env%mpi%groupComm%size, dp)
  #:else
    this%camEnergy = this%camEnergy + evaluateEnergy_real(this%hprev, tmpDeltaRhoSqr)
    ! this%camEnergy = this%camEnergy + evaluateEnergy_real(tmpHSqr, tmpDeltaRhoSqr)
  #:endif

  end subroutine addCamHamiltonianNeighbour_gamma


  !> Adds CAM range-separated contributions to Hamiltonian, using neighbour-list based algorithm.
  !! (k-point version)
  subroutine addCamHamiltonianNeighbour_kpts_mic(this, env, deltaRhoSqr, symNeighbourList,&
      & nNeighbourCamSym, rCellVecs, cellVecs, iSquare, orb, kPoints, kWeights, iKiSToiGlobalKS,&
      & HSqr)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) delta spin-density matrix at BvK real-space shifts
    real(dp), intent(in) :: deltaRhoSqr(:,:,:,:,:,:)

    !> list of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Vectors to neighboring unit cells in absolute units
    real(dp), intent(in) :: rCellVecs(:,:)

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVecs(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> K-points in relative coordinates to calculate delta H(k) for
    real(dp), intent(in) :: kPoints(:,:)

    !> K-point weights (for energy contribution)
    real(dp), intent(in) :: kWeights(:)

    !> Composite index for mapping iK/iS --> iGlobalKS for arrays present at every MPI rank
    integer, intent(in) :: iKiSToiGlobalKS(:,:)

    !> Square (unpacked) Hamiltonian for all k-point/spin composite indices to be updated
    complex(dp), intent(inout) :: HSqr(:,:,:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Atom blocks from sparse, real-space overlap matrices S_{\alpha\mu}, S_{\beta\nu}
    real(dp), pointer :: pSam(:,:), pSbn(:,:)
    real(dp), allocatable :: pSamT(:,:)

    !! \tilde{\gamma}_{\mu\nu}, \tilde{\gamma}_{\mu\beta},
    !! \tilde{\gamma}_{\alpha\nu}, \tilde{\gamma}_{\alpha\beta}
    real(dp), allocatable :: gammaMN(:), gammaMB(:), gammaAN(:), gammaAB(:)

    !! Density matrix block \alpha\beta
    real(dp), allocatable :: Pab(:,:,:,:,:)

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descA(descLen), descB(descLen), descM(descLen), descN(descLen)

    !! Temporary storages
    real(dp), allocatable :: deltaDeltaRhoSqr(:,:,:,:,:,:)
    complex(dp), allocatable :: tmpHSqr(:,:,:)

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Real-space \vec{h} and \vec{l} vectors in relative coordinates
    real(dp) :: vecH(3), vecL(3)

    !! Real-space \vec{h} and \vec{l} vectors in absolute coordinates
    real(dp) :: rVecH(3), rVecL(3)

    !! Temporary arrays for gemm operations
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab, pSamT_Pab_pSbn, Pab_Sbn
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab_gammaAB, pSamT_Pab_gammaMB_pSbn
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab_Sbn_gammaAN, pSamT_Pab_gammaAB_pSbn
    complex(dp) :: tot(orb%mOrb, orb%mOrb, size(kPoints, dim=2) * size(deltaRhoSqr, dim=6))

    !! K-point-spin compound index
    integer :: iGlobalKS, iLocalKS

    !! K-point index
    integer :: iK

    !! Spin index
    integer :: iS

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
    real(dp) :: pMax, maxEstimate, pSbnPabMax

    !! Integer BvK index
    integer :: bvKIndex(3)

    !! Phase factor
    complex(dp) :: phase

    !! Species indices
    integer :: iSpM, iSpN, iSpA, iSpB

    !! Iterates over all BvK real-space vectors
    integer :: iG

    !! Start and end index for MPI parallelization, if applicable
    integer :: iParallelStart, iParallelEnd

    !! Composite index for four nested loops
    integer, allocatable :: compositeIndex(:,:)
    integer :: ii, indGlobal, indLocal

    !! Dummy array with zeros
    real(dp) :: zeros(3)

    !! Number of k-points and spins
    integer :: nK, nS

    !! Number of k-point-spin compound indices
    integer :: nKS

    !! Pre-factor
    real(dp) :: fac

    zeros(:) = 0.0_dp

    tot(:,:,:) = (0.0_dp, 0.0_dp)

    squareSize = size(HSqr, dim=1)
    nAtom0 = size(this%species0)
    nS = size(deltaRhoSqr, dim=6)
    nK = size(kPoints, dim=2)
    nKS = nK * nS

    ! check and initialize screening
    if (.not. this%tScreeningInited) then
      allocate(this%hprevCplxHS(squareSize, squareSize, nKS))
      this%hprevCplxHS(:,:,:) = (0.0_dp, 0.0_dp)
      this%dRhoPrevCplxHS = deltaRhoSqr
      ! there is no previous delta density matrix, therefore just copy over
      deltaDeltaRhoSqr = deltaRhoSqr
      this%tScreeningInited = .true.
    else
      ! allocate and initialize difference of delta rho to previous SCC iteration
      deltaDeltaRhoSqr = deltaRhoSqr - this%dRhoPrevCplxHS
    end if

    pMax = maxval(abs(deltaDeltaRhoSqr))
    ! store delta density matrix
    this%dRhoPrevCplxHS = deltaRhoSqr

    ! skip whole procedure if delta density matrix is close to zero, e.g. in the first SCC iteration
    if (pMax < 1e-16_dp) return

    call getFourLoopCompositeIndex(this, env, nNeighbourCamSym, pMax, compositeIndex)

    ! allocate delta Hamiltonian
    allocate(tmpHSqr(squareSize, squareSize, nKS))
    tmpHSqr(:,:,:) = (0.0_dp, 0.0_dp)

    allocate(gammaMN(size(this%cellVecsG, dim=2)))
    allocate(gammaMB(size(this%cellVecsG, dim=2)))
    allocate(gammaAN(size(this%cellVecsG, dim=2)))
    allocate(gammaAB(size(this%cellVecsG, dim=2)))

    loopMNBA: do ii = 1, size(compositeIndex, dim=2)
      ! Recover indices from composite
      iAtM = compositeIndex(1, ii)
      iAtN = compositeIndex(2, ii)
      iNeighNsort = compositeIndex(3, ii)
      iNeighMsort = compositeIndex(4, ii)

      ! Former loopMN
      iSpM = this%species0(iAtM)
      iSpN = this%species0(iAtN)
      descM = getDescriptor(iAtM, iSquare)
      descN = getDescriptor(iAtN, iSquare)

      ! Former loopB
      iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighNsort, iAtN)
      iAtBfold = symNeighbourList%img2CentCell(iAtB)
      descB = getDescriptor(iAtBfold, iSquare)
      iSpB = this%species0(iAtBfold)
      ! get real-space \vec{l} for gamma arguments
      rVecL(:) = rCellVecs(:, symNeighbourList%iCellVec(iAtB))
      vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
      ! get 2D pointer to Sbn overlap block
      ind = symNeighbourList%iPair(iNeighNsort, iAtN) + 1
      nOrbAt = descN(iNOrb)
      nOrbNeigh = descB(iNOrb)
      pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

      ! Former loopA
      iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighMsort, iAtM)
      iAtAfold = symNeighbourList%img2CentCell(iAtA)
      descA = getDescriptor(iAtAfold, iSquare)
      iSpA = this%species0(iAtAfold)
      ! get real-space \vec{h} for gamma arguments
      rVecH(:) = rCellVecs(:, symNeighbourList%iCellVec(iAtA))
      vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
      ! get 2D pointer to Sam overlap block
      ind = symNeighbourList%iPair(iNeighMsort, iAtM) + 1
      nOrbAt = descM(iNOrb)
      nOrbNeigh = descA(iNOrb)
      ! S_{\alpha\mu}
      pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
      ! S_{\mu\alpha}
      pSamT = transpose(pSam(1:nOrbNeigh, 1:nOrbAt))

      gammaMN(:) = getCamGammaGResolved(this, iAtM, iAtN, iSpM, iSpN, this%rCellVecsG,&
          & rVecH - rVecL)
      gammaMB(:) = getCamGammaGResolved(this, iAtM, iAtBfold, iSpM, iSpB, this%rCellVecsG, rVecH)
      gammaAN(:) = getCamGammaGResolved(this, iAtAfold, iAtN, iSpA, iSpN, this%rCellVecsG, -rVecL)
      gammaAB(:) = getCamGammaGResolved(this, iAtAfold, iAtBfold, iSpA, iSpB, this%rCellVecsG,&
          & zeros)

      tot(:,:,:) = (0.0_dp, 0.0_dp)

      loopS: do iS = 1, nS

        ! get continuous 2D copy of Pab density matrix block
        Pab = deltaDeltaRhoSqr(descA(iStart):descA(iEnd), descB(iStart):descB(iEnd), :,:,:, iS)

        loopG: do iG = 1, size(this%cellVecsG, dim=2)
          bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iG))

          ! term #1/2
          pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) = matmul(pSamT, Pab(:,:, bvKIndex(1),&
              & bvKIndex(2), bvKIndex(3)))

          ! term #1
          pSamT_Pab_pSbn(1:descM(iNOrb), 1:descN(iNOrb)) = matmul(pSamT_Pab(1:descM(iNOrb),&
              & 1:descB(iNOrb)), pSbn)

          ! term #2
          pSamT_Pab_gammaMB_pSbn(1:descM(iNOrb), 1:descN(iNOrb))&
              & = matmul(pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) * gammaMB(iG), pSbn)

          ! term #3
          Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) = matmul(Pab(:,:, bvKIndex(1), bvKIndex(2),&
              & bvKIndex(3)), pSbn)
          pSamT_Pab_Sbn_gammaAN(1:descM(iNOrb), 1:descN(iNOrb))&
              & = matmul(pSamT, Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) * gammaAN(iG))

          ! term #4
          pSamT_Pab_gammaAB(1:descM(iNorb), 1:descB(iNorb)) = matmul(pSamT, Pab(:,:,&
              & bvKIndex(1), bvKIndex(2), bvKIndex(3)) * gammaAB(iG))
          pSamT_Pab_gammaAB_pSbn(1:descM(iNOrb), 1:descN(iNOrb))&
              & = matmul(pSamT_Pab_gammaAB(1:descM(iNOrb), 1:descB(iNOrb)), pSbn)

          loopK: do iK = 1, nK

            iGlobalKS = iKiSToiGlobalKS(iK, iS)

            phase = exp(cmplx(0, -1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                & this%cellVecsG(:, iG) - vecL + vecH))

            ! term #1
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_pSbn(1:descM(iNOrb), 1:descN(iNOrb)) * gammaMN(iG), 0, dp)&
                & * phase

            ! term #2
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_gammaMB_pSbn(1:descM(iNOrb), 1:descN(iNOrb)), 0, dp) * phase

            ! term #3
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_Sbn_gammaAN(1:descM(iNOrb), 1:descN(iNOrb)), 0, dp) * phase

            ! term #4
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_gammaAB_pSbn(1:descM(iNOrb), 1:descN(iNOrb)), 0, dp) * phase

          end do loopK

        end do loopG

      end do loopS

      tmpHSqr(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd), :)&
          & = tmpHSqr(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd), :)&
          & + tot(1:descM(iNOrb), 1:descN(iNOrb), :)

    end do loopMNBA

    if (this%tSpin .or. this%tREKS) then
      tmpHSqr(:,:,:) = -0.25_dp * tmpHSqr
    else
      tmpHSqr(:,:,:) = -0.125_dp * tmpHSqr
    end if

  #:if WITH_MPI
    ! Sum up contributions of current MPI group
    call mpifx_allreduceip(env%mpi%globalComm, tmpHSqr, MPI_SUM)
  #:endif

    if (env%tGlobalLead) then
      this%hprevCplxHS(:,:,:) = this%hprevCplxHS + tmpHSqr
      HSqr(:,:,:) = HSqr + this%hprevCplxHS
    end if

  #:if WITH_MPI
    call mpifx_bcast(env%mpi%globalComm, this%hprevCplxHS)
    call mpifx_bcast(env%mpi%globalComm, HSqr)
  #:endif

  end subroutine addCamHamiltonianNeighbour_kpts_mic


  !> Adds range-separated contributions to Hamiltonian, using neighbour-list based algorithm.
  !! (k-point version)
  subroutine addCamHamiltonianNeighbour_kpts_ct(this, env, deltaRhoSqr, symNeighbourList,&
      & nNeighbourCamSym, cellVecs, iSquare, orb, kPoints, kWeights, iKiSToiGlobalKS, HSqr)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) delta spin-density matrix at BvK real-space shifts
    real(dp), intent(in) :: deltaRhoSqr(:,:,:,:,:,:)

    !> list of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVecs(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> K-points in relative coordinates to calculate delta H(k) for
    real(dp), intent(in) :: kPoints(:,:)

    !> K-point weights (for energy contribution)
    real(dp), intent(in) :: kWeights(:)

    !> Composite index for mapping iK/iS --> iGlobalKS for arrays present at every MPI rank
    integer, intent(in) :: iKiSToiGlobalKS(:,:)

    !> Square (unpacked) Hamiltonian for all k-point/spin composite indices to be updated
    complex(dp), intent(inout) :: HSqr(:,:,:)

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
    real(dp), allocatable :: deltaDeltaRhoSqr(:,:,:,:,:,:)
    complex(dp), allocatable :: tmpHSqr(:,:,:)

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Real-space \vec{h} and \vec{l} vectors in relative coordinates
    real(dp) :: vecH(3), vecL(3)

    !! Temporary arrays for gemm operations
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab, pSamT_Pab_pSbn, Pab_Sbn
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab_gammaAB, pSamT_Pab_gammaMB_pSbn
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab_Sbn_gammaAN, pSamT_Pab_gammaAB_pSbn
    complex(dp) :: tot(orb%mOrb, orb%mOrb, size(kPoints, dim=2) * size(deltaRhoSqr, dim=6))

    !! K-point-spin compound index
    integer :: iGlobalKS

    !! K-point index
    integer :: iK

    !! Spin index
    integer :: iS

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
    real(dp) :: pMax, maxEstimate, pSbnMax, pSbnPabMax, PabMax

    !! Integer BvK index
    integer :: bvKIndex(3)

    !! Phase factor
    complex(dp) :: phase

    !! Iterates over all BvK real-space vectors
    integer :: iG, iGMN, iGMB, iGAN, iGAB

    !! Start and end index for MPI parallelization, if applicable
    integer :: iParallelStart, iParallelEnd

    !! Composite index for four nested loops
    integer, allocatable :: compositeIndex(:,:)
    integer :: ii, indGlobal, indLocal

    !! Number of k-points and spins
    integer :: nK, nS

    !! Number of k-point-spin compound indices
    integer :: nKS

    !! Pre-factor
    real(dp) :: fac

    tot(:,:,:) = (0.0_dp, 0.0_dp)

    squareSize = size(HSqr, dim=1)
    nAtom0 = size(this%species0)
    nS = size(deltaRhoSqr, dim=6)
    nK = size(kPoints, dim=2)
    nKS = nK * nS

    ! check and initialize screening
    if (.not. this%tScreeningInited) then
      allocate(this%hprevCplxHS(squareSize, squareSize, nKS))
      this%hprevCplxHS(:,:,:) = (0.0_dp, 0.0_dp)
      this%dRhoPrevCplxHS = deltaRhoSqr
      ! there is no previous delta density matrix, therefore just copy over
      deltaDeltaRhoSqr = deltaRhoSqr
      this%tScreeningInited = .true.
    else
      ! allocate and initialize difference of delta rho to previous SCC iteration
      deltaDeltaRhoSqr = deltaRhoSqr - this%dRhoPrevCplxHS
    end if

    pMax = maxval(abs(deltaDeltaRhoSqr))
    ! store delta density matrix
    this%dRhoPrevCplxHS = deltaRhoSqr

    ! skip whole procedure if delta density matrix is close to zero, e.g. in the first SCC iteration
    if (pMax < 1e-16_dp) return

    call getFourLoopCompositeIndex(this, env, nNeighbourCamSym, pMax, compositeIndex)

    ! allocate delta Hamiltonian
    allocate(tmpHSqr(squareSize, squareSize, nKS))
    tmpHSqr(:,:,:) = (0.0_dp, 0.0_dp)

    loopMNBA: do ii = 1, size(compositeIndex, dim=2)
      ! Recover indices from composite
      iAtM = compositeIndex(1, ii)
      iAtN = compositeIndex(2, ii)
      iNeighNsort = compositeIndex(3, ii)
      iNeighMsort = compositeIndex(4, ii)

      ! Former loopMN
      descM = getDescriptor(iAtM, iSquare)
      descN = getDescriptor(iAtN, iSquare)

      ! Former loopB
      iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighNsort, iAtN)
      iAtBfold = symNeighbourList%img2CentCell(iAtB)
      descB = getDescriptor(iAtBfold, iSquare)
      ! get real-space \vec{l} for gamma arguments
      vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
      ! get 2D pointer to Sbn overlap block
      ind = symNeighbourList%iPair(iNeighNsort, iAtN) + 1
      nOrbAt = descN(iNOrb)
      nOrbNeigh = descB(iNOrb)
      pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

      ! Former loopA
      iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighMsort, iAtM)
      iAtAfold = symNeighbourList%img2CentCell(iAtA)
      descA = getDescriptor(iAtAfold, iSquare)
      ! get real-space \vec{h} for gamma arguments
      vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
      ! get 2D pointer to Sam overlap block
      ind = symNeighbourList%iPair(iNeighMsort, iAtM) + 1
      nOrbAt = descM(iNOrb)
      nOrbNeigh = descA(iNOrb)
      ! S_{\alpha\mu}
      pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
      pSamT = transpose(pSam(1:nOrbNeigh, 1:nOrbAt))

      tot(:,:,:) = (0.0_dp, 0.0_dp)

      loopS: do iS = 1, nS

        ! get continuous 2D copy of Pab density matrix block
        Pab = deltaDeltaRhoSqr(descA(iStart):descA(iEnd), descB(iStart):descB(iEnd), :,:,:, iS)

        loopGMN: do iG = 1, size(this%nNonZeroGammaG(iAtM, iAtN)%array)
          iGMN = this%nNonZeroGammaG(iAtM, iAtN)%array(iG)
          bvKIndex(:) = this%foldToBvKIndex(vecH - vecL - this%cellVecsG(:, iGMN))

          pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) = matmul(pSamT,&
              & Pab(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3)))
          pSamT_Pab_pSbn(1:descM(iNOrb), 1:descN(iNOrb)) = matmul(pSamT_Pab(1:descM(iNOrb),&
              & 1:descB(iNOrb)), pSbn)

          loopKMN: do iK = 1, nK
            iGlobalKS = iKiSToiGlobalKS(iK, iS)
            phase = exp(cmplx(0, -1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                & this%cellVecsG(:, iGMN)))

            ! term #1
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_pSbn(1:descM(iNOrb), 1:descN(iNOrb))&
                & * this%camGammaEvalG(iAtM, iAtN)%array(iG), 0, dp) * phase
          end do loopKMN
        end do loopGMN

        loopGMB: do iG = 1, size(this%nNonZeroGammaG(iAtM, iAtBfold)%array)
          iGMB = this%nNonZeroGammaG(iAtM, iAtBfold)%array(iG)
          bvKIndex(:) = this%foldToBvKIndex(vecH - this%cellVecsG(:, iGMB))

          pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) = matmul(pSamT,&
              & Pab(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3)))
          pSamT_Pab_gammaMB_pSbn(1:descM(iNOrb), 1:descN(iNOrb))&
              & = matmul(pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb))&
              & * this%camGammaEvalG(iAtM, iAtBfold)%array(iG), pSbn)

          loopKMB: do iK = 1, nK
            iGlobalKS = iKiSToiGlobalKS(iK, iS)
            phase = exp(cmplx(0, -1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                & this%cellVecsG(:, iGMB) - vecL))

            ! term #2
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_gammaMB_pSbn(1:descM(iNOrb), 1:descN(iNOrb)), 0, dp) * phase
          end do loopKMB
        end do loopGMB

        loopGAN: do iG = 1, size(this%nNonZeroGammaG(iAtAfold, iAtN)%array)
          iGAN = this%nNonZeroGammaG(iAtAfold, iAtN)%array(iG)
          bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iGAN) - vecL)

          Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) = matmul(Pab(:,:, bvKIndex(1), bvKIndex(2),&
              & bvKIndex(3)), pSbn)
          pSamT_Pab_Sbn_gammaAN(1:descM(iNOrb), 1:descN(iNOrb))&
              & = matmul(pSamT, Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb))&
              & * this%camGammaEvalG(iAtAfold, iAtN)%array(iG))

          loopKAN: do iK = 1, nK
            iGlobalKS = iKiSToiGlobalKS(iK, iS)
            phase = exp(cmplx(0, -1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                & this%cellVecsG(:, iGAN) + vecH))

            ! term #3
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_Sbn_gammaAN(1:descM(iNOrb), 1:descN(iNOrb)), 0, dp) * phase
          end do loopKAN
        end do loopGAN

        loopGAB: do iG = 1, size(this%nNonZeroGammaG(iAtAfold, iAtBfold)%array)
          iGAB = this%nNonZeroGammaG(iAtAfold, iAtBfold)%array(iG)
          bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iGAB))

          pSamT_Pab_gammaAB(1:descM(iNorb), 1:descB(iNorb)) = matmul(pSamT, Pab(:,:,&
              & bvKIndex(1), bvKIndex(2), bvKIndex(3))&
              & * this%camGammaEvalG(iAtAfold, iAtBfold)%array(iG))
          pSamT_Pab_gammaAB_pSbn(1:descM(iNOrb), 1:descN(iNOrb))&
              & = matmul(pSamT_Pab_gammaAB(1:descM(iNOrb), 1:descB(iNOrb)), pSbn)

          loopKAB: do iK = 1, nK
            iGlobalKS = iKiSToiGlobalKS(iK, iS)
            phase = exp(cmplx(0, -1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                & this%cellVecsG(:, iGAB) - vecL + vecH))

            ! term #4
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_gammaAB_pSbn(1:descM(iNOrb), 1:descN(iNOrb)), 0, dp) * phase
          end do loopKAB
        end do loopGAB

      end do loopS

      tmpHSqr(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd), :)&
          & = tmpHSqr(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd), :)&
          & + tot(1:descM(iNOrb), 1:descN(iNOrb), :)

    end do loopMNBA

    if (this%tSpin .or. this%tREKS) then
      tmpHSqr(:,:,:) = -0.25_dp * tmpHSqr
    else
      tmpHSqr(:,:,:) = -0.125_dp * tmpHSqr
    end if

  #:if WITH_MPI
    ! Sum up contributions of current MPI group
    call mpifx_allreduceip(env%mpi%globalComm, tmpHSqr, MPI_SUM)
  #:endif

    if (env%tGlobalLead) then
      this%hprevCplxHS(:,:,:) = this%hprevCplxHS + tmpHSqr
      HSqr(:,:,:) = HSqr + this%hprevCplxHS
    end if

  #:if WITH_MPI
    call mpifx_bcast(env%mpi%globalComm, this%hprevCplxHS)
    call mpifx_bcast(env%mpi%globalComm, HSqr)
  #:endif

  end subroutine addCamHamiltonianNeighbour_kpts_ct


  !> Adds the CAM-energy contribution to the total energy.
  subroutine addCamEnergy(this, env, energy)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Total energy
    real(dp), intent(inout) :: energy

  #:if WITH_MPI
    call mpifx_allreduceip(env%mpi%globalComm, this%camEnergy, MPI_SUM)
  #:endif

    energy = energy + this%camEnergy

    ! hack for spin unrestricted calculation
    this%camEnergy = 0.0_dp

  end subroutine addCamEnergy


  !> Adds the CAM-energy contribution to the total energy.
  subroutine addCamEnergy_kpts(this, env, iKiSToiGlobalKS, kWeights, deltaRhoOutSqrCplx, energy)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Composite index for mapping iK/iS --> iGlobalKS for arrays present at every MPI rank
    integer, intent(in) :: iKiSToiGlobalKS(:,:)

    !> K-point weights
    real(dp), intent(in) :: kWeights(:)

    !> Complex, dense, square k-space delta density matrix of all spins/k-points
    complex(dp), intent(in) :: deltaRhoOutSqrCplx(:,:,:)

    !> Total energy
    real(dp), intent(inout) :: energy

    !! Spin and k-point indices
    integer :: iS, iK

    !! Global iKS for arrays present at every MPI rank
    integer :: iGlobalKS

    if (env%tGlobalLead) then
      do iK = 1, size(kWeights)
        do iS = 1, nint(real(size(deltaRhoOutSqrCplx, dim=3), dp) / real(size(kWeights), dp))
          iGlobalKS = iKiSToiGlobalKS(iK, iS)
          this%camEnergy = this%camEnergy&
              & + evaluateEnergy_cplx(this%hprevCplxHS(:,:, iGlobalKS), kWeights(iK),&
              & transpose(deltaRhoOutSqrCplx(:,:, iGlobalKS)))
        end do
      end do
    end if

  #:if WITH_MPI
    call mpifx_bcast(env%mpi%globalComm, this%camEnergy)
  #:endif

    energy = energy + this%camEnergy

    ! hack for spin unrestricted calculation
    this%camEnergy = 0.0_dp

  end subroutine addCamEnergy_kpts


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
  pure function evaluateEnergy_cplx(hamiltonian, kWeight, densityMat) result(energy)

    !> Hamiltonian matrix
    complex(dp), intent(in) :: hamiltonian(:,:)

    !> K-point weight (for energy contribution)
    real(dp), intent(in) :: kWeight

    !> Density matrix in k-space
    complex(dp), intent(in) :: densityMat(:,:)

    !> Resulting energy due to CAM contribution
    real(dp) :: energy

    energy = 0.5_dp * real(sum(hamiltonian * densityMat), dp) * kWeight

  end function evaluateEnergy_cplx


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


  !> Returns the numerical second derivative of long-range gamma, by means of a central finite
  !! difference.
  function getddLrNumericalGammaValue(this, iSp1, iSp2, dist, delta) result(ddGamma)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

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

    ddGamma = getddLrNumericalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
        & dist, delta)

  end function getddLrNumericalGammaValue


  !> Workhorse routine for getddLrNumericalGammaValue.
  function getddLrNumericalGammaValue_workhorse(hubbu1, hubbu2, omega, dist, delta) result(ddGamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Delta for finite differences
    real(dp), intent(in) :: delta

    !> Numerical gamma derivative
    real(dp) :: ddGamma

    ddGamma = (getdLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega, dist + delta)&
        & - getdLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega, dist - delta))&
        & / (2.0_dp * delta)

  end function getddLrNumericalGammaValue_workhorse


  !> Returns the numerical second derivative of full-range Hartree-Fock gamma, by means of a central
  !! finite difference.
  function getddHfNumericalGammaDeriv(this, iSp1, iSp2, dist, delta) result(ddGamma)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

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

    ddGamma = (&
        & getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist + delta)&
        & - getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist - delta))&
        & / (2.0_dp * delta)

  end function getddHfNumericalGammaDeriv


  !> Calculates analytical, screened Coulomb, long-range gamma.
  function getLrScreenedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_screened), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting screened gamma
    real(dp) :: gamma

    gamma = getLrScreenedGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
        & this%auxiliaryScreening, dist)

  end function getLrScreenedGammaValue


  !> Calculates analytical, screened Coulomb, long-range gamma.
  function getLrScreenedGammaValue_workhorse(hubbu1, hubbu2, omega, auxiliaryScreening, dist)&
      & result(gamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Auxiliary screening parameter
    real(dp), intent(in) :: auxiliaryScreening

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting screened gamma
    real(dp) :: gamma

    gamma = getLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega + auxiliaryScreening, dist)&
        & - getLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, auxiliaryScreening, dist)

  end function getLrScreenedGammaValue_workhorse


  !> Calculates analytical, screened and poly5zero damped Coulomb, long-range gamma.
  function getLrScreenedAndDampedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_screenedAndDamped), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting screened gamma
    real(dp) :: gamma

    gamma = getLrScreenedAndDampedGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2),&
        & this%omega, this%auxiliaryScreening, dist, this%gammaCutoff, this%gammaDamping,&
        & this%lrScreenedGammaAtDamping(iSp1, iSp2), this%lrdScreenedGammaAtDamping(iSp1, iSp2),&
        & this%lrddScreenedGammaAtDamping(iSp1, iSp2))

  end function getLrScreenedAndDampedGammaValue


  !> Calculates analytical, screened and poly5zero damped Coulomb, long-range gamma.
  function getLrScreenedAndDampedGammaValue_workhorse(hubbu1, hubbu2, omega, auxiliaryScreening,&
      & dist, gammaCutoff, gammaDamping, lrScreenedGammaAtDamping, lrdScreenedGammaAtDamping,&
      & lrddScreenedGammaAtDamping) result(gamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Auxiliary screening parameter
    real(dp), intent(in) :: auxiliaryScreening

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Cutoff for truncated Gamma
    real(dp), intent(in) :: gammaCutoff

    !> Damping distance for Gamma truncation
    real(dp), intent(in) :: gammaDamping

    !> Value, 1st and 2nd derivative of screened gamma integral at damping distance
    real(dp), intent(in) :: lrScreenedGammaAtDamping
    real(dp), intent(in) :: lrdScreenedGammaAtDamping
    real(dp), intent(in) :: lrddScreenedGammaAtDamping

    !> Resulting screened gamma
    real(dp) :: gamma

    if (dist > gammaDamping .and. dist < gammaCutoff) then
      gamma = poly5zero(lrScreenedGammaAtDamping, lrdScreenedGammaAtDamping,&
          & lrddScreenedGammaAtDamping, dist, gammaDamping, gammaCutoff, tDerivative=.false.)
    elseif (dist >= gammaCutoff) then
      gamma = 0.0_dp
    else
      gamma = getLrScreenedGammaValue_workhorse(hubbu1, hubbu2, omega, auxiliaryScreening, dist)
    end if

  end function getLrScreenedAndDampedGammaValue_workhorse


  !> Calculates analytical, screened Coulomb, full-range Hartree-Fock gamma.
  function getHfScreenedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_screened), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting screened gamma
    real(dp) :: gamma

    gamma = 0.0_dp

  end function getHfScreenedGammaValue


  !> Calculates analytical, screened and poly5zero damped Coulomb, full-range Hartree-Fock gamma.
  function getHfScreenedAndDampedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_screenedAndDamped), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting screened gamma
    real(dp) :: gamma

    gamma = 0.0_dp

  end function getHfScreenedAndDampedGammaValue


  !> Calculates analytical, truncated Coulomb, long-range gamma.
  function getLrTruncatedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_truncated), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    if (dist >= this%gammaCutoff) then
      gamma = 0.0_dp
    else
      gamma = getLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
          & dist)
    end if

  end function getLrTruncatedGammaValue


  !> Calculates analytical, long-range gamma for MIC algorithm.
  function getLrMicGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_mic), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    gamma = getLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
        & dist)

  end function getLrMicGammaValue


  !> Calculates analytical, truncated and poly5zero damped Coulomb, long-range gamma.
  function getLrTruncatedAndDampedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_truncatedAndDamped), intent(in) :: this

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
      gamma = getLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
          & dist)
    end if

  end function getLrTruncatedAndDampedGammaValue


  !> Calculates analytical, truncated Coulomb, full-range Hartree-Fock gamma.
  function getHfTruncatedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_truncated), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    if (dist >= this%gammaCutoff) then
      gamma = 0.0_dp
    else
      gamma = getHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)
    end if

  end function getHfTruncatedGammaValue


  !> Calculates analytical, full-range Hartree-Fock gamma for MIC algorithm.
  function getHfMicGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_mic), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    gamma = getHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)

  end function getHfMicGammaValue


  !> Calculates analytical, truncated and poly5zero damped Coulomb, full-range Hartree-Fock gamma.
  function getHfTruncatedAndDampedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_truncatedAndDamped), intent(in) :: this

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
      gamma = getHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)
    end if

  end function getHfTruncatedAndDampedGammaValue


  !> Returns g-sum of long-range and full-range Hartree-Fock gamma's.
  function getCamGammaGSum(this, iAt1, iAt2, iSp1, iSp2, rCellVecsG) result(gamma)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

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

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGLr: do iG = 1, size(rCellVecsG, dim=2)
        dist = norm2(this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG)))
        gamma = gamma + this%camBeta * this%getLrGammaValue(iSp1, iSp2, dist)
      end do loopGLr
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGHf: do iG = 1, size(rCellVecsG, dim=2)
        dist = norm2(this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG)))
        gamma = gamma + this%camAlpha * this%getHfGammaValue(iSp1, iSp2, dist)
      end do loopGHf
    end if

  end function getCamGammaGSum


  !> Returns g-resolved, long-range + full-range Hartree-Fock gamma.
  function getCamGammaGResolved(this, iAt1, iAt2, iSp1, iSp2, rCellVecsG, rShift)&
      & result(gammas)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> Absolute shift of atom position (Gamma arguments)
    real(dp), intent(in) :: rShift(:)

    !> Resulting Gammas, g-vector resolved
    real(dp) :: gammas(size(rCellVecsG, dim=2))

    !! Total shift of atom position (Gamma arguments) in absolute coordinates
    real(dp) :: rTotshift(3)

    !! Distance between the two atoms
    real(dp) :: dist

    !! Index of real-space \vec{g} summation
    integer :: iG

    gammas(:) = 0.0_dp

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGLr: do iG = 1, size(rCellVecsG, dim=2)
        rTotshift(:) = rCellVecsG(:, iG) + rShift
        dist = norm2(this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rTotshift))
        gammas(iG) = gammas(iG) + this%camBeta * this%getLrGammaValue(iSp1, iSp2, dist)
      end do loopGLr
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGHf: do iG = 1, size(rCellVecsG, dim=2)
        rTotshift(:) = rCellVecsG(:, iG) + rShift
        dist = norm2(this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rTotshift))
        gammas(iG) = gammas(iG) + this%camAlpha * this%getHfGammaValue(iSp1, iSp2, dist)
      end do loopGHf
    end if

  end function getCamGammaGResolved


  !> Calculates analytical, unaltered long-range gamma derivative.
  function getdLrAnalyticalGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_full), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma derivative
    real(dp) :: dGamma

    dGamma = getdLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
        & dist)

  end function getdLrAnalyticalGammaValue


  !> Calculates analytical, unaltered full-range Hartree-Fock gamma derivative.
  function getdHfAnalyticalGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_full), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma derivative
    real(dp) :: dGamma

    dGamma = getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)

  end function getdHfAnalyticalGammaValue


  !> Calculates analytical, truncated Coulomb, long-range gamma derivative.
  function getdLrTruncatedGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_truncated), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative
    real(dp) :: dGamma

    if (dist >= this%gammaCutoff) then
      dGamma = 0.0_dp
    else
      dGamma = getdLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
          & dist)
    end if

  end function getdLrTruncatedGammaValue


  !> Calculates analytical, long-range gamma derivative for MIC algorithm.
  function getdLrMicGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_mic), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative
    real(dp) :: dGamma

    dGamma = getdLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
        & dist)

  end function getdLrMicGammaValue


  !> Calculates analytical, truncated Coulomb, full-range Hartree-Fock gamma derivative.
  function getdHfTruncatedGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_truncated), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative
    real(dp) :: dGamma

    if (dist >= this%gammaCutoff) then
      dGamma = 0.0_dp
    else
      dGamma = getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)
    end if

  end function getdHfTruncatedGammaValue


  !> Calculates analytical, full-range Hartree-Fock gamma derivative for MIC algorithm.
  function getdHfMicGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_mic), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative
    real(dp) :: dGamma

    dGamma = getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)

  end function getdHfMicGammaValue


  !> Calculates analytical, truncated and poly5zero damped Coulomb, long-range gamma derivative.
  function getdLrTruncatedAndDampedGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_truncatedAndDamped), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative
    real(dp) :: dGamma

    if (dist > this%gammaDamping .and. dist < this%gammaCutoff) then
      dGamma = poly5zero(this%lrGammaAtDamping(iSp1, iSp2), this%lrdGammaAtDamping(iSp1, iSp2),&
          & this%lrddGammaAtDamping(iSp1, iSp2), dist, this%gammaDamping, this%gammaCutoff,&
          & tDerivative=.true.)
    elseif (dist >= this%gammaCutoff) then
      dGamma = 0.0_dp
    else
      dGamma = getdLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
          & dist)
    end if

  end function getdLrTruncatedAndDampedGammaValue


  !> Calculates analytical, truncated and poly5zero damped Coulomb, full-range Hartree-Fock gamma
  !! derivative.
  function getdHfTruncatedAndDampedGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_truncatedAndDamped), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative (1st)
    real(dp) :: dGamma

    if (dist > this%gammaDamping .and. dist < this%gammaCutoff) then
      dGamma = poly5zero(this%lrGammaAtDamping(iSp1, iSp2), this%lrdGammaAtDamping(iSp1, iSp2),&
          & this%lrddGammaAtDamping(iSp1, iSp2), dist, this%gammaDamping, this%gammaCutoff,&
          & tDerivative=.true.)
    elseif (dist >= this%gammaCutoff) then
      dGamma = 0.0_dp
    else
      dGamma = getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)
    end if

  end function getdHfTruncatedAndDampedGammaValue


  !> Calculates analytical, screened Coulomb, long-range gamma derivative.
  function getdLrScreenedGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_screened), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting screened gamma derivative
    real(dp) :: dGamma

    dGamma = getdLrScreenedGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2),&
        & this%omega, this%auxiliaryScreening, dist)

  end function getdLrScreenedGammaValue


  !> Calculates analytical, screened Coulomb, long-range gamma derivative.
  function getdLrScreenedGammaValue_workhorse(hubbu1, hubbu2, omega, auxiliaryScreening, dist)&
      & result(dGamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Auxiliary screening parameter
    real(dp), intent(in) :: auxiliaryScreening

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting screened gamma derivative
    real(dp) :: dGamma

    dGamma = getdLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega + auxiliaryScreening, dist)&
        & - getdLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, auxiliaryScreening, dist)

  end function getdLrScreenedGammaValue_workhorse


  !> Calculates analytical, screened Coulomb, long-range gamma (2nd) derivative.
  function getddLrNumericalScreenedGammaValue_workhorse(hubbu1, hubbu2, omega, auxiliaryScreening,&
      & dist, delta) result(ddGamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Auxiliary screening parameter
    real(dp), intent(in) :: auxiliaryScreening

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Delta for finite differences
    real(dp), intent(in) :: delta

    !> Resulting screened gamma (2nd) derivative
    real(dp) :: ddGamma

    ddGamma = getddLrNumericalGammaValue_workhorse(hubbu1, hubbu2, omega + auxiliaryScreening,&
        & dist, delta)&
        & - getddLrNumericalGammaValue_workhorse(hubbu1, hubbu2, auxiliaryScreening, dist, delta)

  end function getddLrNumericalScreenedGammaValue_workhorse


  !> Calculates analytical, screened and poly5zero damped Coulomb, long-range gamma derivative.
  function getdLrScreenedAndDampedGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_screenedAndDamped), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting screened gamma derivative
    real(dp) :: dGamma

    dGamma = getdLrScreenedAndDampedGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2),&
        & this%omega, this%auxiliaryScreening, dist, this%gammaCutoff, this%gammaDamping,&
        & this%lrScreenedGammaAtDamping(iSp1, iSp2), this%lrdScreenedGammaAtDamping(iSp1, iSp2),&
        & this%lrddScreenedGammaAtDamping(iSp1, iSp2))

  end function getdLrScreenedAndDampedGammaValue


  !> Calculates analytical, screened and poly5zero damped Coulomb, long-range gamma derivative.
  function getdLrScreenedAndDampedGammaValue_workhorse(hubbu1, hubbu2, omega, auxiliaryScreening,&
      & dist, gammaCutoff, gammaDamping, lrScreenedGammaAtDamping, lrdScreenedGammaAtDamping,&
      & lrddScreenedGammaAtDamping) result(dGamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Auxiliary screening parameter
    real(dp), intent(in) :: auxiliaryScreening

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Cutoff for truncated Gamma
    real(dp), intent(in) :: gammaCutoff

    !> Damping distance for Gamma truncation
    real(dp), intent(in) :: gammaDamping

    !> Value, 1st and 2nd derivative of screened gamma integral at damping distance
    real(dp), intent(in) :: lrScreenedGammaAtDamping
    real(dp), intent(in) :: lrdScreenedGammaAtDamping
    real(dp), intent(in) :: lrddScreenedGammaAtDamping

    !> Resulting screened gamma derivative
    real(dp) :: dGamma

    if (dist > gammaDamping .and. dist < gammaCutoff) then
      dGamma = poly5zero(lrScreenedGammaAtDamping, lrdScreenedGammaAtDamping,&
          & lrddScreenedGammaAtDamping, dist, gammaDamping, gammaCutoff, tDerivative=.true.)
    elseif (dist >= gammaCutoff) then
      dGamma = 0.0_dp
    else
      dGamma = getdLrScreenedGammaValue_workhorse(hubbu1, hubbu2, omega, auxiliaryScreening, dist)
    end if

  end function getdLrScreenedAndDampedGammaValue_workhorse


  !> Calculates analytical, screened Coulomb, full-range Hartree-Fock gamma derivative.
  function getdHfScreenedGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_screened), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting screened gamma derivative
    real(dp) :: dGamma

    dGamma = 0.0_dp

  end function getdHfScreenedGammaValue


  !> Calculates analytical, screened and poly5zero damped Coulomb, full-range Hartree-Fock gamma
  !! derivative.
  function getdHfScreenedAndDampedGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_screenedAndDamped), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting screened gamma derivative
    real(dp) :: dGamma

    dGamma = 0.0_dp

  end function getdHfScreenedAndDampedGammaValue


  !> Returns g-sum of long-range and full-range Hartree-Fock gamma derivatives.
  function getCamGammaPrimeGSum(this, iAt1, iAt2, iSp1, iSp2, rCellVecsG) result(dGamma)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> Resulting Gamma derivative (1st), summed up for g-vectors
    real(dp) :: dGamma(3)

    !! Temporary distance vector
    real(dp) :: distVect(3)

    !! Distance between the two atoms
    real(dp) :: dist

    !! Index of real-space \vec{g} summation
    integer :: iG

    dGamma(:) = 0.0_dp

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGLr: do iG = 1, size(rCellVecsG, dim=2)
        distVect(:) = this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG))
        dist = norm2(distVect)
        distVect(:) = distVect / dist
        dGamma(:) = dGamma + distVect * this%camBeta * this%getLrGammaPrimeValue(iSp1, iSp2, dist)
      end do loopGLr
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGHf: do iG = 1, size(rCellVecsG, dim=2)
        distVect(:) = this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG))
        dist = norm2(distVect)
        distVect(:) = distVect / dist
        dGamma(:) = dGamma + distVect * this%camAlpha * this%getHfGammaPrimeValue(iSp1, iSp2, dist)
      end do loopGHf
    end if

  end function getCamGammaPrimeGSum


  !> Returns g-resolved long-range and full-range Hartree-Fock gamma derivatives.
  function getCamGammaPrimeGResolved(this, iAt1, iAt2, iSp1, iSp2, rCellVecsG) result(dGammas)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> Resulting Gamma derivative (1st), for each g-vector
    real(dp) :: dGammas(3, size(rCellVecsG, dim=2))

    !! Temporary distance vector
    real(dp) :: distVect(3)

    !! Distance between the two atoms
    real(dp) :: dist

    !! Index of real-space \vec{g} summation
    integer :: iG

    dGammas(:,:) = 0.0_dp

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGLr: do iG = 1, size(rCellVecsG, dim=2)
        distVect(:) = this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG))
        dist = norm2(distVect)
        distVect(:) = distVect / dist
        dGammas(:, iG) = dGammas(:, iG)&
            & + distVect * this%camBeta * this%getLrGammaPrimeValue(iSp1, iSp2, dist)
      end do loopGLr
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGHf: do iG = 1, size(rCellVecsG, dim=2)
        distVect(:) = this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG))
        dist = norm2(distVect)
        distVect(:) = distVect / dist
        dGammas(:, iG) = dGammas(:, iG)&
            & + distVect * this%camAlpha * this%getHfGammaPrimeValue(iSp1, iSp2, dist)
      end do loopGHf
    end if

  end function getCamGammaPrimeGResolved


  !> Workhorse for calculating CAM gamma integral.
  function getCamAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega, camAlpha, camBeta, dist)&
      & result(gamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> CAM alpha and beta parameter
    real(dp), intent(in) :: camAlpha, camBeta

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting CAM gamma
    real(dp) :: gamma

    gamma =&
        & camAlpha * getHfAnalyticalGammaValue_workhorse(hubbu1, hubbu2, dist)&
        & + camBeta * getLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega, dist)

  end function getCamAnalyticalGammaValue_workhorse


  !> Calculates analytical long-range gamma.
  function getLrAnalyticalGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_full), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma
    real(dp) :: gamma

    gamma = getLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
        & dist)

  end function getLrAnalyticalGammaValue


  !> Workhorse for getLrAnalyticalGammaValue wrapper.
  function getLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega, dist) result(gamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma
    real(dp) :: gamma

    real(dp) :: tauA, tauB
    real(dp) :: prefac, tmp, tmp2, tau

    tauA = 3.2_dp * hubbu1
    tauB = 3.2_dp * hubbu2

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff) then
        tau = 0.5_dp * (tauA + tauB)
        tmp = 5.0_dp * tau**6 + 15.0_dp * tau**4 * omega**2 - 5.0_dp * tau**2 * omega**4 + omega**6
        tmp = tmp * 0.0625_dp / tau**5 - omega
        tmp = tmp * tau**8 / (tau**2 - omega**2)**4
        gamma = tau * 0.3125_dp - tmp
      else
        call error("Error(HybridXc): R = 0, Ua != Ub")
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

  end function getLrAnalyticalGammaValue_workhorse


  !> Returns the subexpression for the evaluation of the off-site Y-Gamma-integral.
  pure function getYGammaSubPart(tauA, tauB, R, omega) result(yGamma)

    !> Decay constant site A
    real(dp), intent(in) :: tauA

    !> Decay constant site B
    real(dp), intent(in) :: tauB

    !> Separation of the sites A and B
    real(dp), intent(in) :: R

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Resulting off-site Y-Gamma-integral
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


  !> Returns 1st analytical derivative of long-range gamma.
  function getdLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega, dist) result(dGamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting d gamma / d dist
    real(dp) :: dGamma

    real(dp) :: tauA, tauB
    real(dp) :: prefac, tmp, tmp2, dTmp, dTmp2

    tauA = 3.2_dp * hubbu1
    tauB = 3.2_dp * hubbu2

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff) then
        dGamma = 0.0_dp
      else
        call error("Error(HybridXc1): R = 0, Ua != Ub")
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

        dGamma = -1.0_dp/dist**2 -dtmp2&
            & + (tauA**8 / (tauA**2 - omega**2)**4) * (dtmp + dtmp2 + omega*exp(-omega * dist)/dist&
            & +exp(-omega * dist) / dist**2)

      else
        ! off-site, Ua != Ub
        prefac = tauA**4 / (tauA * tauA - omega * omega )**2
        prefac = prefac * tauB**4 / (tauB * tauB - omega * omega )**2
        prefac = prefac * (-omega * exp(-omega * dist) / dist - exp(-omega * dist) / dist**2)
        dGamma = -1.0_dp / (dist**2) - prefac&
            & + getdYGammaSubPart(tauA, tauB, dist, omega)&
            & + getdYGammaSubPart(tauB, tauA, dist, omega)&
            & - getdYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getdYGammaSubPart(tauB, tauA, dist, 0.0_dp)
      end if
    end if

  end function getdLrAnalyticalGammaValue_workhorse


  !> Returns the derivative of the subexpression for the evaluation of the off-site
  !! Y-Gamma-integral. Note that tauA /= tauB.
  pure function getdYGammaSubPart(tauA, tauB, R, omega) result(dYGammaSubPart)

    !> Decay constant site A
    real(dp), intent(in) :: tauA

    !> Decay constant site B
    real(dp), intent(in) :: tauB

    !> Separation of the sites A and B
    real(dp), intent(in) :: R

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Resulting derivative of the subexpression
    real(dp) :: dYGammaSubPart

    !! Auxiliary variables
    real(dp) :: prefac, tmp, tmp2, dtmp

    tmp = tauA**2 - omega**2
    prefac = tauA * tauA / tmp
    tmp = prefac * prefac / (tauA * tauA - tauB * tauB)**3
    dtmp = tmp * (tauB**6 - 3.0_dp * tauA * tauA * tauB**4 + 2.0_dp * omega * omega * tauB**4)/R**2
    tmp = tmp * (tauB**6 - 3.0_dp * tauA * tauA * tauB**4 + 2.0_dp * omega * omega * tauB**4) / R
    tmp2 = tauA * tauB**4 * 0.5_dp * prefac / (tauB * tauB - tauA * tauA )**2 - tmp

    dYGammaSubPart = (dtmp - tmp2 * tauA) * exp(-tauA * R)

  end function getdYGammaSubPart


  !> Returns the derivative of CAM range-separated gamma for iAtom1, iAtom2 (non-periodic version).
  subroutine getDirectedCamGammaPrimeValue_cluster(this, grad, iAtom1, iAtom2)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Gradient of gamma between atoms
    real(dp), intent(out) :: grad(3)

    !> First atom
    integer, intent(in) :: iAtom1

    !> Second atom
    integer, intent(in) :: iAtom2

    !! Species index of first and second atom
    integer :: iSp1, iSp2

    !! Distance(-vector) of the two atoms
    real(dp) :: vect(3), dist

    grad(:) = 0.0_dp

    iSp1 = this%species0(iAtom1)
    iSp2 = this%species0(iAtom2)

    ! analytical derivatives
    vect(:) = this%rCoords(:, iAtom1) - this%rCoords(:, iAtom2)
    dist = sqrt(sum(vect**2))
    vect(:) = vect / dist

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      grad(:) = grad + vect * this%camBeta * this%getLrGammaPrimeValue(iSp1, iSp2, dist)
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      grad(:) = grad + vect * this%camAlpha * this%getHfGammaPrimeValue(iSp1, iSp2, dist)
    end if

  end subroutine getDirectedCamGammaPrimeValue_cluster


  !> Returns the derivative of CAM range-separated gamma for iAtom1, iAtom2 (periodic version).
  subroutine getDirectedCamGammaPrimeValue_periodic(this, grad, iAtom1, iAtom2, img2CentCell)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Gradient of gamma between atoms
    real(dp), intent(out) :: grad(3)

    !> First atom
    integer, intent(in) :: iAtom1

    !> Second atom
    integer, intent(in) :: iAtom2

    !> Map images of atoms to the central cell
    integer, intent(in) :: img2CentCell(:)

    !! Species index of first and second atom
    integer :: iSp1, iSp2

    !! Distance(-vector) of the two atoms
    real(dp) :: vect(3), dist

    iSp1 = this%species0(img2CentCell(iAtom1))
    iSp2 = this%species0(img2CentCell(iAtom2))

    ! analytical derivatives
    vect(:) = this%rCoords(:, iAtom1) - this%rCoords(:, iAtom2)
    dist = sqrt(sum(vect**2))
    vect(:) = vect / dist

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      grad(:) = grad + vect * this%camBeta * this%getLrGammaPrimeValue(iSp1, iSp2, dist)
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      grad(:) = grad + vect * this%camAlpha * this%getHfGammaPrimeValue(iSp1, iSp2, dist)
    end if

  end subroutine getDirectedCamGammaPrimeValue_periodic


  !> Adds gradients due to CAM range-separated contributions (non-periodic version).
  subroutine addCamGradients_cluster(this, gradients, derivator, deltaRho, skOverCont, orb,&
      & iSquare, ovrlapMat, iNeighbour, nNeighbourSK)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> Density matrix difference from reference q0
    real(dp), intent(in) :: deltaRho(:,:,:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

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
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, this%rCoords, this%species0, iAtK,&
              & iAtC, orb)
          call derivator%getFirstDeriv(sPrimeTmp2, skOverCont, this%rCoords, this%species0, iAtC,&
              & iAtK, orb)
        end if
        loopB: do iAtB = 1, nAtom0
          ! A > B
          loopA: do iNeighB = 0, nNeighbourSK(iAtB)
            iAtA = iNeighbour(iNeighB, iAtB)
            tmpgamma1 = this%camGammaEval0(iAtK, iAtB) + this%camGammaEval0(iAtC, iAtB)
            tmpgamma2 = tmpgamma1 + this%camGammaEval0(iAtK, iAtA) + this%camGammaEval0(iAtC, iAtA)
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
      gradients(:,:) = gradients - 0.5_dp * tmpderiv
    else
      gradients(:,:) = gradients - 0.25_dp * nSpin * tmpderiv
    end if


  contains

    !> Initialise.
    subroutine allocateAndInit(tmpOvr, tmpRho, gammaPrimeTmp, tmpderiv)

      !> Storage for the overlap
      real(dp), allocatable, intent(inout) :: tmpOvr(:,:)

      !> Storage for density matrix
      real(dp), allocatable, intent(inout) :: tmpRho(:,:,:)

      !> Storage for derivative of gamma interaction, shape: [nCoords (x,y,z), nAtom0, nAtom0]
      real(dp), allocatable, intent(inout) :: gammaPrimeTmp(:,:,:)

      !> Workspace for the derivatives
      real(dp), allocatable, intent(inout) :: tmpderiv(:,:)

      !! Holds long-range gamma derivatives of a single interaction
      real(dp) :: tmp(3)

      !! Atom indices
      integer :: iAt1, iAt2

      !! Spin channel index
      integer :: iSpin

      !! Number of atoms
      integer :: nAtom0

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

      ! Precompute the gamma derivatives
      gammaPrimeTmp(:,:,:) = 0.0_dp
      do iAt1 = 1, nAtom0
        do iAt2 = 1, nAtom0
          if (iAt1 /= iAt2) then
            call getDirectedCamGammaPrimeValue_cluster(this, tmp, iAt1, iAt2)
            gammaPrimeTmp(:, iAt1, iAt2) = tmp
          end if
        end do
      end do
    end subroutine allocateAndInit

  end subroutine addCamGradients_cluster


  !> Adds CAM gradients due to CAM range-separated contributions (Gamma-only version).
  subroutine addCamGradients_gamma(this, deltaRhoSqr, skOverCont, symNeighbourList,&
      & nNeighbourCamSym, iSquare, orb, derivator, gradients)

    !> Class instance
    class(THybridXcFunc), intent(in), target :: this

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

    !! Atom blocks from sparse, real-space overlap matrices
    real(dp), pointer :: pSam(:,:), pSak(:,:), pSbn(:,:), pSbm(:,:), pSbk(:,:)

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descA(descLen), descB(descLen), descM(descLen), descN(descLen), descK(descLen)

    !! Temporary storages
    real(dp), allocatable :: tmpDeltaRhoSqr(:,:,:), tmpGradients(:,:)

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Overlap matrix elements
    real(dp) :: Sam, Sak, Sbm, Sbn, Sbk

    !! Density matrix elements
    real(dp) :: dPmn, dPab, dPba, dPmk, dPkm

    !! Overlap derivatives
    real(dp), dimension(orb%mOrb, orb%mOrb, 3) :: SbnPrimeKequalsN

    !! \tilde{\gamma}_{\mu\kappa}, \tilde{\gamma}_{\alpha\kappa},
    !! \tilde{\gamma}_{\alpha\beta}
    real(dp) :: gammaMK, gammaAK, gammaAB, gammaMKMB, gammaTot

    !! 1st derivatives of
    !! \tilde{\gamma}_{\mu\nu}, \tilde{\gamma}_{\mu\beta},
    !! \tilde{\gamma}_{\alpha\nu}, \tilde{\gamma}_{\alpha\beta}
    real(dp), dimension(3) :: dGammaMK, dGammaAK, dGammaKB, dGammaAB

    !! Atom to calculate energy gradient components for
    integer :: iAtK, iNeighK

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
    integer :: mu, nu, kk, alpha, beta

    !! Product dPmn * Sam
    real(dp) :: dPmnSam

    !! Product dPkm * Sak
    real(dp) :: dPkmSak

    !! Product dPmn * gammaTot
    real(dp) :: dPmnGammaTot

    !! Product dPmk * Sam
    real(dp) :: dPmkSam

    !! Product dPmk * gammaTot
    real(dp) :: dPmkGammaTot

    !! Product dPmk * gammaTot * Sam
    real(dp) :: dPmkGammaTotSam

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

    loopK1: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopM1: do iAtM = 1, nAtom0
        descM = getDescriptor(iAtM, iSquare)
        gammaMK = this%camGammaEval0(iAtM, iAtK)
        loopB1: do iNeighK = 0, nNeighbourCamSym(iAtK)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          if (iAtK == iAtBfold) cycle
          descB = getDescriptor(iAtBfold, iSquare)
          ! \tilde{\gamma}_{\mu\beta}
          gammaMKMB = gammaMK + this%camGammaEval0(iAtM, iAtBfold)
          SbnPrimeKequalsN(:,:,:) = 0.0_dp
          call derivator%getFirstDeriv(SbnPrimeKequalsN, skOverCont, this%rCoords,&
              & symNeighbourList%species, iAtK, iAtB, orb)
          loopA1: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            descA = getDescriptor(iAtAfold, iSquare)
            ! \tilde{\gamma}_{\alpha\nu}
            gammaAK = this%camGammaEval0(iAtAfold, iAtK)
            ! \tilde{\gamma}_{\alpha\beta}
            gammaAB = this%camGammaEval0(iAtAfold, iAtBfold)
            ! get 2D pointer to Sam overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            gammaTot = gammaMKMB + gammaAK + gammaAB

            do mu = 1, descM(iNOrb)
              do kk = 1, descK(iNOrb)
                do iSpin = 1, nSpin
                  dPmk = tmpDeltaRhoSqr(descM(iStart) + mu - 1, descK(iStart) + kk - 1, iSpin)
                  dPmkGammaTot = dPmk * gammaTot
                  do alpha = 1, descA(iNOrb)
                    Sam = pSam(alpha, mu)
                    dPmkGammaTotSam = dPmkGammaTot * Sam
                    do beta = 1, descB(iNOrb)
                      dPab = tmpDeltaRhoSqr(descA(iStart) + alpha - 1,&
                          & descB(iStart) + beta - 1, iSpin)

                      tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                          & + dPmkGammaTotSam * dPab * SbnPrimeKequalsN(beta, kk, :)
                    end do
                  end do
                end do
              end do
            end do

          end do loopA1
        end do loopB1
      end do loopM1
    end do loopK1

    loopK2: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopM2: do iAtM = 1, nAtom0
        if (iAtK == iAtM) cycle
        descM = getDescriptor(iAtM, iSquare)
        dGammaMK(:) = -this%camdGammaEval0(iAtM, iAtK, :)
        loopB2: do iNeighK = 0, nNeighbourCamSym(iAtK)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          descB = getDescriptor(iAtBfold, iSquare)
          ! get 2D pointer to Sbk overlap block
          ind = symNeighbourList%iPair(iNeighK, iAtK) + 1
          nOrbAt = descK(iNOrb)
          nOrbNeigh = descB(iNOrb)
          pSbk(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          loopA2: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            descA = getDescriptor(iAtAfold, iSquare)
            ! get 2D pointer to Sam overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            do mu = 1, descM(iNOrb)
              do kk = 1, descK(iNOrb)
                do iSpin = 1, nSpin
                  dPmk = tmpDeltaRhoSqr(descM(iStart) + mu - 1, descK(iStart) + kk - 1, iSpin)
                  do alpha = 1, descA(iNOrb)
                    Sam = pSam(alpha, mu)
                    dPmkSam = dPmk * Sam
                    do beta = 1, descB(iNOrb)
                      Sbk = pSbk(beta, kk)
                      dPab = tmpDeltaRhoSqr(descA(iStart) + alpha - 1,&
                          & descB(iStart) + beta - 1, iSpin)

                      tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                          & + dPmkSam * dPab * Sbk * dGammaMK
                    end do
                  end do
                end do
              end do
            end do

          end do loopA2
        end do loopB2
      end do loopM2
    end do loopK2

    loopK3: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopM3: do iAtM = 1, nAtom0
        descM = getDescriptor(iAtM, iSquare)
        loopA3: do iNeighM = 0, nNeighbourCamSym(iAtM)
          iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
          iAtAfold = symNeighbourList%img2CentCell(iAtA)
          if (iAtK == iAtAfold) cycle
          dGammaAK(:) = -this%camdGammaEval0(iAtAfold, iAtK, :)
          descA = getDescriptor(iAtAfold, iSquare)
          ! get 2D pointer to Sam overlap block
          ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
          nOrbAt = descM(iNOrb)
          nOrbNeigh = descA(iNOrb)
          pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          loopB3: do iNeighK = 0, nNeighbourCamSym(iAtK)
            iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
            iAtBfold = symNeighbourList%img2CentCell(iAtB)
            descB = getDescriptor(iAtBfold, iSquare)
            ! get 2D pointer to Sbk overlap block
            ind = symNeighbourList%iPair(iNeighK, iAtK) + 1
            nOrbAt = descK(iNOrb)
            nOrbNeigh = descB(iNOrb)
            pSbk(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            do mu = 1, descM(iNOrb)
              do kk = 1, descK(iNOrb)
                do iSpin = 1, nSpin
                  dPmk = tmpDeltaRhoSqr(descM(iStart) + mu - 1, descK(iStart) + kk - 1, iSpin)
                  do alpha = 1, descA(iNOrb)
                    Sam = pSam(alpha, mu)
                    dPmkSam = dPmk * Sam
                    do beta = 1, descB(iNOrb)
                      Sbk = pSbk(beta, kk)
                      dPab = tmpDeltaRhoSqr(descA(iStart) + alpha - 1,&
                          & descB(iStart) + beta - 1, iSpin)

                      tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                          & + dPmkSam * dPab * Sbk * dGammaAK
                    end do
                  end do
                end do
              end do
            end do

          end do loopB3
        end do loopA3
      end do loopM3
    end do loopK3

    if (this%tREKS) then
      gradients(:,:) = gradients - 0.5_dp * tmpGradients
    else
      gradients(:,:) = gradients - 0.25_dp * nSpin * tmpGradients
    end if

  end subroutine addCamGradients_gamma


  !> Adds CAM gradients due to CAM range-separated contributions (k-point version).
  subroutine addCamGradients_kpts_ct(this, deltaRhoSqr, deltaRhoOutSqrCplx, symNeighbourList,&
      & nNeighbourCamSym, cellVecs, iSquare, orb, kPoints, kWeights, skOverCont, derivator,&
      & gradients)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> Square (unpacked) delta spin-density matrix at BvK real-space shifts
    real(dp), intent(in) :: deltaRhoSqr(:,:,:,:,:,:)

    !> Square (unpacked) delta spin-density matrix in k-space
    complex(dp), intent(in) :: deltaRhoOutSqrCplx(:,:,:)

    !> list of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVecs(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> K-points in relative coordinates to calculate delta H(k) for
    real(dp), intent(in) :: kPoints(:,:)

    !> K-point weights (for energy contribution)
    real(dp), intent(in) :: kWeights(:)

    !> Sparse overlap container
    type(TSlakoCont), intent(in) :: skOverCont

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Atom blocks from sparse, real-space overlap matrices
    real(dp), pointer :: pSak(:,:), pSam(:,:), pSbk(:,:), pSbn(:,:)

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descA(descLen), descB(descLen), descM(descLen), descN(descLen), descK(descLen)

    !! Temporary storages
    real(dp), allocatable :: tmpGradients(:,:)

    !! Overlap derivatives
    real(dp), dimension(orb%mOrb, orb%mOrb, 3) :: SbnPrimeKequalsN

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Real-space \vec{h} and \vec{l} vectors in relative coordinates
    real(dp) :: vecH(3), vecL(3)

    !! K-point index
    integer :: iK

    !! Spin index
    integer :: iS

    !! Atom indices (central cell)
    integer :: iAtM, iAtN, iAtK

    !! Neighbour indices (+corresponding atom indices)
    integer :: iNeighK, iNeighM, iNeighN, iAtA, iAtB

    !! Folded (to central cell) atom indices
    integer :: iAtAfold, iAtBfold

    !! Auxiliary variables for setting up 2D pointer to sparse overlap
    integer :: ind, nOrbAt, nOrbNeigh

    !! Integer BvK index
    integer :: bvKIndex(3)

    !! Phase factor
    complex(dp) :: phase

    !! Iterates over all BvK real-space vectors
    integer :: iG, iGKB, iGMK, iGMB, iGAK, iGAB

    !! Orbital indices
    integer :: mu, nu, kk, alpha, beta

    !! Number of k-points and spins
    integer :: nK, nS

    !! Overlap matrix elements
    real(dp) :: Sam, Sak, Sbk, Sbn

    !! Density matrix elements
    real(dp) :: dPab
    complex(dp) :: dPkm, dPnk, dPnm

    !! Products phase * gammaMK, phase * gammaMB, phase * gammaAK, phase * gammaAB
    complex(dp) :: phaseGammaMK, phaseGammaMB, phaseGammaAK, phaseGammaAB

    !! Product dPkm * phase * gammaMK
    complex(dp) :: dPkmPhaseGammaMK

    !! Product dPkm * phase * gammaMB
    complex(dp) :: dPkmPhaseGammaMB

    !! Product dPkm * phase * gammaAK
    complex(dp) :: dPkmPhaseGammaAK

    !! Product dPkm * phase * gammaAB
    complex(dp) :: dPkmPhaseGammaAB

    !! Product dPkm * phase * gammaMK * Sam
    complex(dp) :: dPkmPhaseGammaMKSam

    !! Product dPkm * phase * gammaMB * Sam
    complex(dp) :: dPkmPhaseGammaMBSam

    !! Product dPkm * phase * gammaAK * Sam
    complex(dp) :: dPkmPhaseGammaAKSam

    !! Product dPkm * phase * gammaAB * Sam
    complex(dp) :: dPkmPhaseGammaABSam

    !! Directed derivatives
    real(dp), allocatable :: dGammaMK(:,:), dGammaAK(:,:), dGammaKB(:,:)

    !! Product phase * dGammaMK(:, iG)
    complex(dp) :: phasedGammaMK(3), phasedGammaKB(3), phasedGammaAK(3)

    complex(dp) :: dPkmPhasedGammaMK(3), dPkmPhasedGammaMKSam(3), dPnkPhasedGammaKB(3)
    complex(dp) :: dPnkPhasedGammaKBSak(3), dPkmPhasedGammaAK(3), dPkmPhasedGammaAKSam(3)
    complex(dp) :: dPnmPhasedGammaAK(3), dPnmPhasedGammaAKSam(3)

    !! Composite index for mapping iK/iS --> iGlobalKS
    integer, allocatable :: iKiSToiGlobalKS(:,:)
    integer :: iGlobalKS

    nAtom0 = size(this%species0)
    nS = size(deltaRhoSqr, dim=6)
    nK = size(kPoints, dim=2)

    ! Build spin/k-point composite index for all spins and k-points (global)
    iGlobalKS = 1
    allocate(iKiSToiGlobalKS(nK, nS))
    do iS = 1, nS
      do iK = 1, nK
        iKiSToiGlobalKS(iK, iS) = iGlobalKS
        iGlobalKS = iGlobalKS + 1
      end do
    end do

    ! allocate gradient contribution
    allocate(tmpGradients(3, size(gradients, dim=2)))
    tmpGradients(:,:) = 0.0_dp

    ! First terms of (48)
    loopK1: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopM1: do iAtM = 1, nAtom0
        descM = getDescriptor(iAtM, iSquare)
        loopB1: do iNeighK = 0, nNeighbourCamSym(iAtK)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          if (iAtBfold == iAtK) cycle
          descB = getDescriptor(iAtBfold, iSquare)
          vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
          SbnPrimeKequalsN(:,:,:) = 0.0_dp
          call derivator%getFirstDeriv(SbnPrimeKequalsN, skOverCont, this%rCoords,&
              & symNeighbourList%species, iAtK, iAtB, orb)
          loopA1: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            descA = getDescriptor(iAtAfold, iSquare)
            vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
            ! get 2D pointer to Sam overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            loopGMK1: do iG = 1, size(this%nNonZeroGammaG(iAtM, iAtK)%array)
              iGMK = this%nNonZeroGammaG(iAtM, iAtK)%array(iG)
              bvKIndex(:) = this%foldToBvKIndex(vecH - vecL - this%cellVecsG(:, iGMK))

              loopKptsMK1: do iK = 1, nK
                phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & -this%cellVecsG(:, iGMK)))
                phaseGammaMK = phase * cmplx(this%camGammaEvalG(iAtM, iAtK)%array(iG)&
                    & * kWeights(iK), 0, dp)

                ! term #1
                do mu = 1, descM(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPkm = deltaRhoOutSqrCplx(descK(iStart) + kk - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPkmPhaseGammaMK = dPkm * phaseGammaMK
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPkmPhaseGammaMKSam = dPkmPhaseGammaMK * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(2.0_dp * dPkmPhaseGammaMKSam&
                              & * cmplx(dPab * SbnPrimeKequalsN(beta, kk, :), 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsMK1
            end do loopGMK1

            loopGMB: do iG = 1, size(this%nNonZeroGammaG(iAtM, iAtBfold)%array)
              iGMB = this%nNonZeroGammaG(iAtM, iAtBfold)%array(iG)
              bvKIndex(:) = this%foldToBvKIndex(vecH - this%cellVecsG(:, iGMB))

              loopKptsMB: do iK = 1, nK
                phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & vecL - this%cellVecsG(:, iGMB)))
                phaseGammaMB = phase * cmplx(this%camGammaEvalG(iAtM, iAtBfold)%array(iG)&
                    & * kWeights(iK), 0, dp)

                ! term #2
                do mu = 1, descM(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPkm = deltaRhoOutSqrCplx(descK(iStart) + kk - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPkmPhaseGammaMB = dPkm * phaseGammaMB
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPkmPhaseGammaMBSam = dPkmPhaseGammaMB * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(2.0_dp * dPkmPhaseGammaMBSam&
                              & * cmplx(dPab * SbnPrimeKequalsN(beta, kk, :), 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsMB
            end do loopGMB

            loopGAK1: do iG = 1, size(this%nNonZeroGammaG(iAtAfold, iAtK)%array)
              iGAK = this%nNonZeroGammaG(iAtAfold, iAtK)%array(iG)
              bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iGAK) - vecL)

              loopKptsAK1: do iK = 1, nK
                phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & -this%cellVecsG(:, iGAK) - vecH))
                phaseGammaAK = phase * cmplx(this%camGammaEvalG(iAtAfold, iAtK)%array(iG)&
                    & * kWeights(iK), 0, dp)

                ! term #3
                do mu = 1, descM(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPkm = deltaRhoOutSqrCplx(descK(iStart) + kk - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPkmPhaseGammaAK = dPkm * phaseGammaAK
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPkmPhaseGammaAKSam = dPkmPhaseGammaAK * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(2.0_dp * dPkmPhaseGammaAKSam&
                              & * cmplx(dPab * SbnPrimeKequalsN(beta, kk, :), 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsAK1
            end do loopGAK1

            loopGAB1: do iG = 1, size(this%nNonZeroGammaG(iAtAfold, iAtBfold)%array)
              iGAB = this%nNonZeroGammaG(iAtAfold, iAtBfold)%array(iG)
              bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iGAB))

              loopKptsAB1: do iK = 1, nK
                phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & -this%cellVecsG(:, iGAB) + vecL - vecH))
                phaseGammaAB = phase * cmplx(this%camGammaEvalG(iAtAfold, iAtBfold)%array(iG)&
                    & * kWeights(iK), 0, dp)

                ! term #4
                do mu = 1, descM(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPkm = deltaRhoOutSqrCplx(descK(iStart) + kk - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPkmPhaseGammaAB = dPkm * phaseGammaAB
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPkmPhaseGammaABSam = dPkmPhaseGammaAB * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(2.0_dp * dPkmPhaseGammaABSam&
                              & * cmplx(dPab * SbnPrimeKequalsN(beta, kk, :), 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsAB1
            end do loopGAB1

          end do loopA1
        end do loopB1
      end do loopM1
    end do loopK1

    ! term #1 with prefactor 1/8 of eq.(48)
    loopK2: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopM2: do iAtM = 1, nAtom0
        if (iAtK == iAtM) cycle
        descM = getDescriptor(iAtM, iSquare)
        dGammaMK = -this%camdGammaEvalG(iAtM, iAtK)%array
        loopB2: do iNeighK = 0, nNeighbourCamSym(iAtK)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          descB = getDescriptor(iAtBfold, iSquare)
          vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
          ! get 2D pointer to Sbk overlap block
          ind = symNeighbourList%iPair(iNeighK, iAtK) + 1
          nOrbAt = descK(iNOrb)
          nOrbNeigh = descB(iNOrb)
          pSbk(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          loopA2: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            descA = getDescriptor(iAtAfold, iSquare)
            vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
            ! get 2D pointer to Sam overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            loopGMK: do iG = 1, size(this%nNonZeroGammaG(iAtM, iAtK)%array)
              iGMK = this%nNonZeroGammaG(iAtM, iAtK)%array(iG)
              bvKIndex(:) = this%foldToBvKIndex(vecH - vecL - this%cellVecsG(:, iGMK))

              loopKptsMK: do iK = 1, nK
                phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & -this%cellVecsG(:, iGMK)))
                phasedGammaMK(:) = phase * cmplx(dGammaMK(:, iG) * kWeights(iK), 0, dp)

                do mu = 1, descM(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPkm = deltaRhoOutSqrCplx(descK(iStart) + kk - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPkmPhasedGammaMK(:) = dPkm * phasedGammaMK
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPkmPhasedGammaMKSam(:) = dPkmPhasedGammaMK * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          Sbk = pSbk(beta, kk)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(dPkmPhasedGammaMKSam * cmplx(Sbk * dPab, 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsMK
            end do loopGMK

          end do loopA2
        end do loopB2
      end do loopM2
    end do loopK2

    ! term #2 with prefactor 1/8 of eq.(48)
    loopK3: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopN1: do iAtN = 1, nAtom0
        descN = getDescriptor(iAtN, iSquare)
        loopB3: do iNeighN = 0, nNeighbourCamSym(iAtN)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighN, iAtN)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          if (iAtK == iAtBfold) cycle
          descB = getDescriptor(iAtBfold, iSquare)
          vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
          dGammaKB = this%camdGammaEvalG(iAtK, iAtBfold)%array
          ! get 2D pointer to Sbn overlap block
          ind = symNeighbourList%iPair(iNeighN, iAtN) + 1
          nOrbAt = descN(iNOrb)
          nOrbNeigh = descB(iNOrb)
          pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          loopA3: do iNeighK = 0, nNeighbourCamSym(iAtK)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            descA = getDescriptor(iAtAfold, iSquare)
            vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
            ! get 2D pointer to Sak overlap block
            ind = symNeighbourList%iPair(iNeighK, iAtK) + 1
            nOrbAt = descK(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSak(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            loopGKB: do iG = 1, size(this%nNonZeroGammaG(iAtK, iAtBfold)%array)
              iGKB = this%nNonZeroGammaG(iAtK, iAtBfold)%array(iG)
              bvKIndex(:) = this%foldToBvKIndex(vecH - this%cellVecsG(:, iGKB))

              loopKptsKB: do iK = 1, nK
                phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & vecL - this%cellVecsG(:, iGKB)))
                phasedGammaKB(:) = phase * cmplx(dGammaKB(:, iG) * kWeights(iK), 0, dp)

                do nu = 1, descN(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPnk = deltaRhoOutSqrCplx(descN(iStart) + nu - 1, descK(iStart) + kk - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPnkPhasedGammaKB(:) = dPnk * phasedGammaKB
                      do alpha = 1, descA(iNOrb)
                        Sak = pSak(alpha, kk)
                        dPnkPhasedGammaKBSak(:) = dPnkPhasedGammaKB * cmplx(Sak, 0, dp)
                        do beta = 1, descB(iNOrb)
                          Sbn = pSbn(beta, nu)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(dPnkPhasedGammaKBSak * cmplx(Sbn * dPab, 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsKB
            end do loopGKB

          end do loopA3
        end do loopB3
      end do loopN1
    end do loopK3

    ! term #3 with prefactor 1/8 of eq.(48)
    loopK4: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopM3: do iAtM = 1, nAtom0
        descM = getDescriptor(iAtM, iSquare)
        loopB4: do iNeighK = 0, nNeighbourCamSym(iAtK)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          descB = getDescriptor(iAtBfold, iSquare)
          vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
          ! get 2D pointer to Sbk overlap block
          ind = symNeighbourList%iPair(iNeighK, iAtK) + 1
          nOrbAt = descK(iNOrb)
          nOrbNeigh = descB(iNOrb)
          pSbk(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          loopA4: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            if (iAtK == iAtAfold) cycle
            descA = getDescriptor(iAtAfold, iSquare)
            vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
            dGammaAK = -this%camdGammaEvalG(iAtAfold, iAtK)%array
            ! get 2D pointer to Sam overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            loopGAK2: do iG = 1, size(this%nNonZeroGammaG(iAtAfold, iAtK)%array)
              iGAK = this%nNonZeroGammaG(iAtAfold, iAtK)%array(iG)
              bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iGAK) - vecL)

              loopKptsAK2: do iK = 1, nK
                phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & -this%cellVecsG(:, iGAK) - vecH))
                phasedGammaAK(:) = phase * cmplx(dGammaAK(:, iG) * kWeights(iK), 0, dp)

                do mu = 1, descM(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPkm = deltaRhoOutSqrCplx(descK(iStart) + kk - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPkmPhasedGammaAK(:) = dPkm * phasedGammaAK
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPkmPhasedGammaAKSam(:) = dPkmPhasedGammaAK * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          Sbk = pSbk(beta, kk)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(dPkmPhasedGammaAKSam * cmplx(Sbk * dPab, 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsAK2
            end do loopGAK2

          end do loopA4
        end do loopB4
      end do loopM3
    end do loopK4

    ! term #4 with prefactor 1/8 of eq.(48)
    loopM4: do iAtM = 1, nAtom0
      descM = getDescriptor(iAtM, iSquare)
      loopN3: do iAtN = 1, nAtom0
        descN = getDescriptor(iAtN, iSquare)
        loopB5: do iNeighN = 0, nNeighbourCamSym(iAtN)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighN, iAtN)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          iAtK = iAtBfold
          descK = getDescriptor(iAtK, iSquare)
          descB = getDescriptor(iAtBfold, iSquare)
          vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
          ! get 2D pointer to Sbn overlap block
          ind = symNeighbourList%iPair(iNeighN, iAtN) + 1
          nOrbAt = descN(iNOrb)
          nOrbNeigh = descB(iNOrb)
          pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          loopA5: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            if (iAtK == iAtAfold) cycle
            descA = getDescriptor(iAtAfold, iSquare)
            vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
            dGammaAK = -this%camdGammaEvalG(iAtAfold, iAtK)%array
            ! get 2D pointer to Sam overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            loopGAK3: do iG = 1, size(this%nNonZeroGammaG(iAtAfold, iAtK)%array)
              iGAK = this%nNonZeroGammaG(iAtAfold, iAtK)%array(iG)
              bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iGAK))

              loopKptsAK3: do iK = 1, nK
                phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & vecL - this%cellVecsG(:, iGAK) - vecH))
                phasedGammaAK(:) = phase * cmplx(dGammaAK(:, iG) * kWeights(iK), 0, dp)

                do mu = 1, descM(iNOrb)
                  do nu = 1, descN(iNOrb)
                    do iS = 1, nS
                      dPnm = deltaRhoOutSqrCplx(descN(iStart) + nu - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPnmPhasedGammaAK(:) = dPnm * phasedGammaAK
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPnmPhasedGammaAKSam(:) = dPnmPhasedGammaAK * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          Sbn = pSbn(beta, nu)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(dPnmPhasedGammaAKSam * cmplx(Sbn * dPab, 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsAK3
            end do loopGAK3

          end do loopA5
        end do loopB5
      end do loopN3
    end do loopM4

    if (this%tREKS) then
      tmpGradients(:,:) = -0.25_dp * tmpGradients
    else
      tmpGradients(:,:) = -0.125_dp * nS * tmpGradients
    end if

    gradients(:,:) = gradients + tmpGradients

  ! #:if WITH_MPI
  !   ! Sum up contributions of current MPI group
  !   call mpifx_allreduceip(env%mpi%globalComm, tmpGradients, MPI_SUM)
  ! #:endif

  !   if (env%tGlobalLead) then
  !     gradients(:,:) = gradients + tmpGradients
  !   end if

  ! #:if WITH_MPI
  !   call mpifx_bcast(env%mpi%globalComm, gradients)
  ! #:endif

  end subroutine addCamGradients_kpts_ct


  !> Explicitly evaluates the LR-Energy contribution (very slow, use addLrEnergy instead).
  function evaluateLrEnergyDirect_cluster(this, env, deltaRho, overlap, iSquare) result(energy)

    !> instance of LR
    class(THybridXcFunc), intent(in) :: this

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
        energy = energy + tmp * this%camGammaEval0(iAt1, iAt2)
      end do
    end do
    energy = -energy / 8.0_dp

    call env%globalTimer%stopTimer(globalTimers%energyEval)

  end function evaluateLrEnergyDirect_cluster


  !> Returns the array of atomic species in central cell.
  subroutine getCentralCellSpecies(this, species0)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> 1D array for output, will be allocated
    integer, intent(out), allocatable :: species0(:)

    species0 = this%species0

  end subroutine getCentralCellSpecies


  !> Returns tabulated (long-range + full-range Hartree-Fock) gamma integrals.
  !! (non-periodic systems only)
  subroutine getCamGammaCluster(this, camGamma0)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Long-range + full-range Hartree-Fock gamma integrals in AO basis
    real(dp), intent(out) :: camGamma0(:,:)

    camGamma0(:,:) = this%camGammaEval0

  end subroutine getCamGammaCluster


  !> Calculates (long-range + full-range Hartree-Fock) gamma derivative integrals.
  !! (non-periodic systems only)
  subroutine getCamGammaDerivCluster(this, camGammaDeriv0)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Long-range + full-range Hartree-Fock gamma derivative integrals
    real(dp), intent(out) :: camGammaDeriv0(:,:,:)

    !! Holds long-range + full-range Hartree-Fock gamma derivatives of a single interaction
    real(dp) :: tmp(3)

    !! Number of atoms and indices of interacting atoms
    integer :: nAtom0, iAt1, iAt2

    nAtom0 = size(camGammaDeriv0, dim=1)

    do iAt1 = 1, nAtom0
      do iAt2 = 1, nAtom0
        if (iAt1 /= iAt2) then
          call getDirectedCamGammaPrimeValue_cluster(this, tmp, iAt1, iAt2)
          camGammaDeriv0(iAt2, iAt1, :) = tmp
        end if
      end do
    end do

  end subroutine getCamGammaDerivCluster


  !> Calculates analytical full-range Hartree-Fock gamma.
  function getHfAnalyticalGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_full), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma
    real(dp) :: gamma

    gamma = getHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)

  end function getHfAnalyticalGammaValue


  !> Workhorse for getHfAnalyticalGammaValue wrapper.
  function getHfAnalyticalGammaValue_workhorse(hubbu1, hubbu2, dist) result(gamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma
    real(dp) :: gamma

    real(dp) :: tauA, tauB
    real(dp) :: tmp, tau

    tauA = 3.2_dp * hubbu1
    tauB = 3.2_dp * hubbu2

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff) then
        tau = 0.5_dp * (tauA + tauB)
        gamma = tau * 0.3125_dp
      else
        call error("Error(HybridXc): R = 0, Ua != Ub")
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

  end function getHfAnalyticalGammaValue_workhorse


  !> Returns analytical derivative of full-range Hartree-Fock gamma.
  function getdHfAnalyticalGammaValue_workhorse(hubbu1, hubbu2, dist) result(dGamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting d gamma / d dist
    real(dp) :: dGamma

    real(dp) :: tauA, tauB
    real(dp) :: dTmp

    tauA = 3.2_dp * hubbu1
    tauB = 3.2_dp * hubbu2

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff) then
        dGamma = 0.0_dp
      else
        call error("Error(HybridXc): R = 0, Ua != Ub")
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

  end function getdHfAnalyticalGammaValue_workhorse

end module dftbp_dftb_hybridxc
