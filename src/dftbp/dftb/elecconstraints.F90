!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module to impose constraints on the electronic ground state.
module dftbp_dftb_elecconstraints
  use dftbp_common_accuracy, only : dp
  use dftbp_math_angmomentum, only : getLOperators
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_typegeometry, only : TGeometry
  use dftbp_extlibs_xmlf90, only : fnode, string, char, getNodeName, getLength, getItem1, fnodeList
  use dftbp_type_wrappedintr, only : TWrappedInt1, TWrappedReal1, TWrappedReal2
  use dftbp_geoopt_package, only : TOptimizerInput
  use dftbp_io_message, only : error
  use dftbp_io_hsdutils, only : getChildValue, getChildren, detailedError, getSelectedAtomIndices
  use dftbp_io_hsdutils2, only : renameChildren, convertUnitHsd
  use dftbp_common_unitconversion, only : timeUnits
  use dftbp_dftbplus_input_geoopt, only : readOptimizerInput
  use dftbp_geoopt_geoopt, only : TGeoOpt, geoOptTypes, next, reset, init
  use dftbp_geoopt_steepdesc, only : TSteepDesc
  use dftbp_geoopt_conjgrad, only : TConjGrad
  use dftbp_geoopt_fire, only : TFire, TFire_init
  use dftbp_geoopt_gdiis, only : TDIIS
  use dftbp_geoopt_lbfgs, only : TLbfgs, TLbfgs_init
  use dftbp_common_globalenv, only : stdOut
  implicit none

  private
  public :: TElecConstraint, TElecConstraint_init, TElecConstraintInput
  public :: readElecConstraintInput


  type TElecConstraintInput

    !> Group of atoms in a constraint
    type(TWrappedInt1), allocatable :: atomGrp(:)

    !> Constraint targets for atom groups
    real(dp), allocatable :: atomNc(:)

    !> Direction of constraint in (q,m) space
    type(TWrappedReal1), allocatable :: atomSpinDir(:)

    !> Optimisation algorithm
    integer :: iConstrOpt

    !> Derivative tolerance for constraint
    real(dp) :: constrTol

    !> Number of iterations for enforcing constraint
    integer :: nConstrIter

    !> Alpha factor to introduce variation into the space
    real(dp) :: diisAlpha

    !> Generations in the history
    integer :: diisGens

    !> Number of stored steps
    integer :: lbfgsMemory

    !> Timestep
    real(dp) :: fireTimeStep

  end type TElecConstraintInput


  type TElecConstraint

    !> Value of the constraint
    real(dp), allocatable :: Nc(:)

    !> Potential
    real(dp), allocatable :: Vc(:)

    ! Weighting function for constrain

    !> Atom(s) involved in each constrain
    type(TWrappedInt1), allocatable :: wAt(:)

    !> Atomic orbital(s) involved in each constrain
    type(TWrappedInt1), allocatable :: wAtOrb(:)

    !> Atomic orbital charge/spin quaternion involved in each constrain
    type(TWrappedReal2), allocatable :: wAtSpin(:)

    !> General optimiser
    type(TGeoOpt), allocatable :: pVcOpt

    !> Number of iterations for enforcing constraint
    integer :: nConstrIter

  contains

    procedure constrain
    procedure getMaxIter

  end type TElecConstraint


contains


  !> General entry point to read constraint on the electronic ground state.
  subroutine readElecConstraintInput(node, geo, input, tSpinPol)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Control structure to be filled
    type(TElecConstraintInput), intent(out) :: input

    !> True, if this is a spin polarized calculation
    logical, intent(in) :: tSpinPol

    type(fnode), pointer :: val, child, child2, child3, value1, field
    type(fnodeList), pointer :: children
    type(string) :: modifier, buffer
    integer :: iConstr, nConstr

    call getChildValue(node, "Driver", child, "lbfgs", child=child2)
    call getNodeName(child, buffer)

    select case(char(buffer))
    case("steepestdescent")
      input%iConstrOpt = geoOptTypes%steepestDesc
    case("conjugategradient")
      input%iConstrOpt = geoOptTypes%conjugateGrad
    case("gdiis")
      input%iConstrOpt = geoOptTypes%diis
      call getChildValue(child, "Alpha", input%diisAlpha, 1.0e-01_dp)
      call getChildValue(child, "Generations", input%diisGens, 8)
    case("lbfgs")
      input%iConstrOpt = geoOptTypes%lbfgs
      call getChildValue(child, "Memory", input%lbfgsMemory, 20)
    case("fire")
      input%iConstrOpt = geoOptTypes%fire
      call getChildValue(child, "TimeStep", input%fireTimeStep, 1.0_dp, modifier=modifier,&
          & child=field)
      call convertUnitHsd(char(modifier), timeUnits, field, input%fireTimeStep)
    case default
      call detailedError(child, "Unknown driver '" // char(buffer) // "'")
    end select

    call getChildValue(node, "Tolerance", input%constrTol, 1.0e-08_dp)
    call getChildValue(node, "Iterations", input%nConstrIter, 200)

    call getChildValue(node, "Regions", val, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! Read specification for regions of atoms
    call getChildren(child, "Atoms", children)
    nConstr = getLength(children)

    allocate(input%atomGrp(nConstr))
    allocate(input%atomNc(nConstr))
    allocate(input%atomSpinDir(nConstr))

    do iConstr = 1, nConstr
      call getItem1(children, iConstr, child2)
      call getChildValue(child2, "Domain", buffer, child=child3, multiple=.true.)
      call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species,&
          & input%atomGrp(iConstr)%data)
      call getChildValue(child2, "Population", input%atomNc(iConstr))
      ! Functionality currently restricted to charges
      if (tSpinPol) then
        allocate(input%atomSpinDir(iConstr)%data(2))
        input%atomSpinDir(iConstr)%data(1) = 1.0_dp
      else
        allocate(input%atomSpinDir(iConstr)%data(1))
        input%atomSpinDir(iConstr)%data(1) = 1.0_dp
      end if
    end do

  end subroutine readElecConstraintInput


  !> Initialises the constraints structure.
  subroutine TElecConstraint_init(this, input, orb)

    !> Constrain structure instance
    type(TElecConstraint), intent(out) :: this

    !> Input data structure
    type(TElecConstraintInput), intent(in) :: input

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    ! Constrain optimiser algorithms

    ! !> Steepest descent driver
    ! type(TSteepDesc), allocatable :: pSteepDesc

    !> Conjugate gradient driver
    type(TConjGrad), allocatable :: pConjGrad

    !> gradient DIIS driver
    type(TDIIS), allocatable :: pDIIS

    !> lBFGS driver for geometry optimisation
    type(TLbfgs), allocatable :: pLbfgs

    !> FIRE driver for geometry optimisation
    type(TFire), allocatable :: pFire

    integer :: iConstr, nConstr, ii, jj, iAt, iOrb, nOrb, nSpin

    nConstr = size(input%atomGrp)

    allocate(this%Vc(nConstr))
    ! should enable optional initialization of Vc from input
    this%Vc(:) = 0.0_dp

    allocate(this%pVcOpt)
    select case(input%iConstrOpt)
    case(geoOptTypes%steepestDesc)
      call error("Steepest descent not yet implemented for electronic constrains.")
    case(geoOptTypes%conjugateGrad)
      allocate(pConjGrad)
      call init(pConjGrad, nConstr, input%constrTol, 5.0e-01_dp)
      call init(this%pVcOpt, pConjGrad)
    case(geoOptTypes%diis)
      allocate(pDIIS)
      call init(pDIIS, nConstr, input%constrTol, input%diisAlpha, input%diisGens)
      call init(this%pVcOpt, pDIIS)
    case(geoOptTypes%lbfgs)
      allocate(pLbfgs)
      call TLbfgs_init(pLbfgs, nConstr, input%constrTol, 1.0e-08_dp, 1.0e-01_dp, input%lbfgsMemory,&
          & .true., .false., .true.)
      call init(this%pVcOpt, pLbfgs)
    case (geoOptTypes%fire)
      allocate(pFire)
      call TFire_init(pFire, nConstr, input%constrTol, input%fireTimeStep)
      call init(this%pVcOpt, pFire)
    end select
    call reset(this%pVcOpt, this%Vc)

    this%nConstrIter = input%nConstrIter
    this%Nc = input%AtomNc

    allocate(this%wAt(nConstr))
    allocate(this%wAtOrb(nConstr))
    allocate(this%wAtSpin(nConstr))

    ! Allocate + initialize arrays and build index mappings
    do iConstr = 1, nConstr
      ! Count orbitals subject to constraints
      nOrb = 0
      do ii = 1, size(input%atomGrp(iConstr)%data)
        iAt = input%atomGrp(iConstr)%data(ii)
        nOrb = nOrb + orb%nOrbAtom(iAt)
      end do
      allocate(this%wAt(iConstr)%data(nOrb))
      allocate(this%wAtOrb(iConstr)%data(nOrb))
      nSpin = size(input%atomSpinDir(iConstr)%data)
      allocate(this%wAtSpin(iConstr)%data(nOrb,nSpin))
      this%wAt(iConstr)%data(:) = 0
      this%wAtOrb(iConstr)%data(:) = 0
      nOrb = 0
      do ii = 1, size(input%atomGrp(iConstr)%data)
        iAt = input%atomGrp(iConstr)%data(ii)
        do jj = 1, orb%nOrbAtom(iAt)
          this%wAt(iConstr)%data(nOrb+jj) = iAt
          this%wAtOrb(iConstr)%data(nOrb+jj) = jj
          this%wAtSpin(iConstr)%data(nOrb+jj,:) = input%atomSpinDir(iConstr)%data
        end do
        nOrb = nOrb + orb%nOrbAtom(iAt)
      end do
    end do

  end subroutine TElecConstraint_init


  !> Returns maximum number of iterations for constraint driver.
  pure function getMaxIter(this) result(maxIter)

    !> Class instance
    class(TElecConstraint), intent(in) :: this

    !> Obtained maximum number of iterations
    integer :: maxIter

    maxIter = this%nConstrIter

  end function getMaxIter


  !> Applies electronic constraints to system.
  subroutine constrain(this, qq, energy, shift, tConverged, dWdVcMax)

    !> Class instance
    class(TElecConstraint), intent(inout) :: this

    !> Mulliken populations
    real(dp), intent(in) :: qq(:,:,:)

    !> Energy
    real(dp), intent(in) :: energy

    !> Potential from constraints
    real(dp), intent(out) :: shift(:,:,:,:)

    !> Gradient convergence achieved
    logical, intent(out) :: tConverged

    !> Maximum derivative of energy functional with respect to Vc
    real(dp), intent(out) :: dWdVcMax

    !> Contribution to free energy functional from constraint(s)
    real(dp) :: deltaW

    !! Derivative of energy functional with respect to Vc
    real(dp), allocatable :: dWdVc(:)

    !! Iterates over constraints
    integer :: iConstr

    !! Number of constraints requested by the user
    integer :: nConstr

    shift(:,:,:,:) = 0.0_dp

    nConstr = size(this%wAt)
    allocate(dWdVc(nConstr))
    dWdVc(:) = 0.0_dp
    deltaW = 0.0_dp

    do iConstr = 1, nConstr

      call constrainQ(shift, deltaW, dWdVc(iConstr), this%Vc(iConstr), this%Nc(iConstr),&
          & this%wAt(iConstr)%data, this%wAtOrb(iConstr)%data, this%wAtSpin(iConstr)%data,&
          & qq)

    end do

    ! Optimizers set up to minimize, therefore sign change in total energy and gradients
    call next(this%pVcOpt, -(energy + deltaW), -dWdVc, this%Vc, tConverged)

    dWdVcMax = maxval(dWdVc)

  end subroutine constrain


  !> Constraint on atomic charge.
  subroutine constrainQ(shift, deltaW, dWdV, Vc, Nc, wAt, wOrb, wSp, qq)

    !> Shift to which contribution is appended
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> Free energy contribution from current contraint
    real(dp), intent(inout) :: deltaW

    !> Derivative of free energy with respect to potential
    real(dp), intent(out) :: dWdV

    !> Potential / Lagrange multiplier
    real(dp), intent(in) :: Vc

    !> Target value of current constrain
    real(dp), intent(in) :: Nc

    !> Atom(s) involved in current constrain
    integer, intent(in) :: wAt(:)

    !> Orbital(s) involved in current constrain
    integer, intent(in) :: wOrb(:)

    !> Spin(s) involved in current constrain
    real(dp), intent(in) :: wSp(:,:)

    !> Mulliken populations
    real(dp), intent(in) :: qq(:,:,:)

    !! Number of spin channels to be constrained
    integer :: nSpin

    !! Index of spin channel
    integer :: iSpin

    !! Index of atomic orbital
    integer :: iW

    !!
    real(dp) :: wn

    nSpin = size(wSp, dim=2)
    wn = 0.0_dp
    do iSpin = 1, nSpin
      do iW = 1, size(wAt)
        wn = wn + wSp(iW, iSpin) * qq(wOrb(iW), wAt(iW), iSpin)
      end do
    end do

    dWdV = wn - Nc
    deltaW = deltaW + Vc * dWdV

    do iSpin = 1, nSpin
      do iW = 1, size(wAt)
        shift(wOrb(iW), wOrb(iW), wAt(iW), iSpin) = shift(wOrb(iW), wOrb(iW), wAt(iW), iSpin)&
            & + Vc * wSp(iW, iSpin)
      end do
    end do

  end subroutine constrainQ

end module dftbp_dftb_elecconstraints
