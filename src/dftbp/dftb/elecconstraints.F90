!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module to impose constraints on the electronic ground state.
module dftbp_dftb_elecconstraints

  use dftbp_common_accuracy, only : dp, sc
  use dftbp_common_globalenv, only : stdOut
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_io_message, only : error
  use dftbp_type_typegeometry, only : TGeometry
  use dftbp_extlibs_xmlf90, only : fnode, string, char, getLength, getItem1, fnodeList
  use dftbp_type_wrappedintr, only : TWrappedInt1, TWrappedReal1, TWrappedReal2
  use dftbp_geoopt_package, only : TOptimizer, TOptimizerInput, createOptimizer
  use dftbp_io_hsdutils, only : detailedError, getChildValue, getChildren, getSelectedAtomIndices
  use dftbp_dftbplus_input_geoopt, only : readOptimizerInput
  use dftbp_extlibs_xmlf90, only : destroyNodeList
  use dftbp_math_blasroutines, only : hemv
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_orbitals, only : orbitalNames, getShellNames
  use dftbp_math_bisect, only : bisection
  use dftbp_io_charmanip, only : i2c
  use dftbp_type_linkedlist, only : len, TListInt, init, asArray, destruct
  use dftbp_dftb_spin, only : qm2ud, ud2qm
  implicit none

  private
  public :: TElecConstraint, TElecConstraint_init, TElecConstraintInput
  public :: readElecConstraintInput
  public :: printHighestAO


  !> Contains input data for electronic constraints.
  type TElecConstraintInput

    !> Optimiser input choice
    class(TOptimizerInput), allocatable :: optimiser

    !> Group of atoms in a constraint
    type(TWrappedInt1), allocatable :: atomGrp(:)

    !> Constraint targets for atom groups
    real(dp), allocatable :: atomNc(:)

    !> Constraint spin for atom groups
    integer, allocatable :: atomSpin(:)

    !> Direction of constraint in (q,m) space
    type(TWrappedReal1), allocatable :: atomSpinDir(:)

    !> Derivative tolerance for constraint
    real(dp) :: constrTol

    !> Number of iterations for enforcing constraint
    integer :: nConstrIter

    !> True, if converged micro-iterations are required
    logical :: isConstrConvRequired

  end type TElecConstraintInput


  !> Represents electronic contraints.
  type TElecConstraint

    !> Value of the constraint
    real(dp), allocatable :: Nc(:)

    !> Spin of constraints
    integer, allocatable :: spin(:)

    !> Potential
    real(dp), allocatable :: Vc(:)

    !> Contribution to free energy functional from constraints
    real(dp), allocatable :: deltaW(:)

    !> Derivative of energy functional with respect to Vc
    real(dp), allocatable :: dWdVc(:)

    ! Weighting function for constrain

    !> Atom(s) involved in each constrain
    type(TWrappedInt1), allocatable :: wAt(:)

    !> Atomic orbital(s) involved in each constrain
    type(TWrappedInt1), allocatable :: wAtOrb(:)

    !> General optimiser
    class(TOptimizer), allocatable :: potOpt

    !> Derivative tolerance for constraint
    real(dp) :: constrTol

    !> Number of iterations for enforcing constraint
    integer :: nConstrIter

    !> True, if converged micro-iterations are required
    logical :: isConstrConvRequired

  contains

    procedure :: getConstrainShift
    procedure :: propagateConstraints
    procedure :: getMaxIter
    procedure :: getFreeEnergy
    procedure :: getMaxEnergyDerivWrtVc

  end type TElecConstraint


contains


  !> Prints atom index, angular and magnetic quantum number of occupied AO with highest energy.
  subroutine printHighestAO(HSqrReal, eigvecs, SSqr, orb, denseDesc, species0, speciesName, filling)

    !> Square (unpacked) Hamiltonian of certain spin channel
    real(dp), intent(in) :: HSqrReal(:,:)

    !> Eigenvectors
    real(dp), intent(in) :: eigvecs(:,:,:)

    !> Square overlap
    real(dp), intent(in) :: SSqr(:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Species of atoms
    integer, intent(in) :: species0(:)

    !> Name of each species
    character(*), intent(in) :: speciesName(:)

    !> Fillings of the states
    real(dp), intent(in) :: filling(:,:,:)

    !! Index of highest occupied kanonical orbital
    integer :: iHoko

    character(sc), allocatable :: shellNamesTmp(:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Stores start/end index and number of orbitals of square matrices
    integer :: desc(descLen)

    integer :: ll, iAtMax, nOrb, iOrbStart, iOrbMax, iSpMax, iOrbOnAtom, mm
    integer :: iSh, iBlockStart, iBlockEnd, degeneracy, iDiag, iState
    integer, allocatable :: magQN(:)
    real(dp) :: probMax, diagMax
    real(dp), allocatable :: hamBlock(:,:)
    character(len=:), allocatable :: fmt
    real(dp), allocatable :: prob(:)

    ! determine HOKO
    lpStates: do iState = size(filling, dim=1), 1, -1
      if (filling(iState, 1, 1) > 1.0e-06_dp) then
        iHoko = iState
        exit lpStates
      end if
    end do lpStates

    write(stdOut, "(/,A)") '################AO-Analysis################'

    ! Check if we are dealing with distributed matrices
    @:ASSERT(size(HSqrReal, dim=1) == size(HSqrReal, dim=2))

    nOrb = size(eigvecs, dim=1)

    allocate(prob(nOrb), source=0.0_dp)
    call hemv(prob, SSqr, eigvecs(:, iHoko, 1))

    ! |<\Psi_iHoko|\Phi_\alpha>|^2
    prob(:) = prob * eigvecs(:, iHoko, 1)
    ! write(stdOut, "(F10.6)") prob

    iOrbMax = maxloc(prob, dim=1)
    probMax = prob(iOrbMax)
    diagMax = HSqrReal(iOrbMax, iOrbMax)

    ! get atom index from orbital index
    call bisection(iAtMax, denseDesc%iAtomStart, iOrbMax)

    ! get species index from atom index
    iSpMax = species0(iAtMax)
    call getShellNames(iSpMax, orb, shellNamesTmp)

    ! first orbital in global indexing on atom
    iOrbStart = denseDesc%iAtomStart(iAtMax)

    ! orbital index of related to the associated atom
    iOrbOnAtom = iOrbMax - iOrbStart + 1

    ll = orb%iShellOrb(iOrbOnAtom, iSpMax) - 1
    iSh = ll + 1
    mm = iOrbOnAtom - orb%posShell(ll + 1, iSpMax) - ll

    desc = getDescriptor(iAtMax, denseDesc%iAtomStart)
    iBlockStart = desc(iStart) + orb%posShell(iSh, iSpMax) - 1
    iBlockEnd = desc(iStart) + orb%posShell(iSh + 1, iSpMax) - 2
    hamBlock = HSqrReal(iBlockStart:iBlockEnd, iBlockStart:iBlockEnd)

    degeneracy = 0
    allocate(magQN(size(hamBlock, dim=1)))
    do iDiag = 1, size(hamBlock, dim=1)
      if (abs(diagMax - hamBlock(iDiag, iDiag)) < 1e-08_dp) then
        degeneracy = degeneracy + 1
        magQN(degeneracy) = iDiag - 1 - ll
      end if
    end do
    magQN = magQN(1:degeneracy)

    write(stdOut, "(A,I0)") "Analysis performed for HOKO: ", iHoko
    write(stdOut, "(A,E20.12,/)") "Highest projection found: ", probMax
    write(stdOut, "(A,I0,5A)") "Global orbital index/type: ", iOrbMax,&
        & " (", trim(shellNamesTmp(ll + 1)), "_", trim(orbitalNames(mm, ll)), ")"
    write(stdOut, "(3A,I0,A)") "Atom type/ global index: ", trim(speciesName(iSpMax)),&
        & " (", iAtMax, ")"
    write(stdOut, "(A)") "Hamiltonian block of full " // trim(shellNamesTmp(ll + 1)) // "-shell:"
    fmt = "(" // i2c(size(hamBlock, dim=1)) // "F16.10)"
    write(stdOut, fmt) transpose(hamBlock)
    fmt = "(A," // i2c(size(magQN, dim=1)) // "I5)"
    write(stdOut, fmt) "Magnetic QNs: ", magQN
    write(stdOut, "(A,I0)") "degeneracy: ", degeneracy
    write(stdOut, "(A,1E24.14,A)") "energy: ", diagMax, " Ha"

    write(stdOut, "(A,/)") '####################END####################'

  contains

    !> Finds location of relevant atomic block indices in a dense matrix.
    pure function getDescriptor(iAt, iSquare) result(desc)

      !> Relevant atom
      integer, intent(in) :: iAt

      !> Indexing array for start of atom orbitals
      integer, intent(in) :: iSquare(:)

      !> Resulting location ranges
      integer :: desc(3)

      desc(:) = [iSquare(iAt), iSquare(iAt + 1) - 1, iSquare(iAt + 1) - iSquare(iAt)]

    end function getDescriptor

  end subroutine printHighestAO


  !> General entry point to read constraint on the electronic ground state.
  subroutine readElecConstraintInput(node, geo, input, isSpinPol)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Control structure to be filled
    type(TElecConstraintInput), intent(out) :: input

    !> True, if this is a spin polarized calculation
    logical, intent(in) :: isSpinPol

    type(fnode), pointer :: val, child1, child2, child3
    type(fnodeList), pointer :: children
    type(string) :: buffer
    integer :: iConstr, nConstr

    call getChildValue(node, "Optimiser", child1, "FIRE")
    call readOptimizerInput(child1, input%optimiser)

    call getChildValue(node, "ConstrTolerance", input%constrTol, 1.0e-08_dp)
    call getChildValue(node, "MaxConstrIterations", input%nConstrIter, 100)
    call getChildValue(node, "ConvergentConstrOnly", input%isConstrConvRequired, .true.)

    call getChildValue(node, "Regions", val, "", child=child1, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! Read specification for regions of atoms
    call getChildren(child1, "Atom", children)
    nConstr = getLength(children)

    if (nConstr > 1) then
      call error("This modified version of DFTB+ only supports a single constraint.")
    end if

    allocate(input%atomGrp(nConstr))
    allocate(input%atomNc(nConstr))
    allocate(input%atomSpin(nConstr))
    allocate(input%atomSpinDir(nConstr))

    do iConstr = 1, nConstr
      call getItem1(children, iConstr, child2)
      call getChildValue(child2, "Index", buffer, child=child3, multiple=.true.)
      call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species,&
          & input%atomGrp(iConstr)%data)
      if (size(input%atomGrp(iConstr)%data) /= 1) then
        call detailedError(child3, "This modified version of DFTB+ only supports constraints on a&
            & single atom.")
      end if
      call getChildValue(child2, "Population", input%atomNc(iConstr))

      if (isSpinPol) then
        call getChildValue(child2, "iSpin", input%atomSpin(iConstr))
      else
        input%atomSpin(iConstr) = 1
      end if

      ! Functionality currently restricted to charges
      if (isSpinPol) then
        input%atomSpinDir(iConstr)%data = [1.0_dp, 0.0_dp]
      else
        input%atomSpinDir(iConstr)%data = [1.0_dp]
      end if
    end do

    call destroyNodeList(children)

  end subroutine readElecConstraintInput


  !> Initialises the constraints structure.
  subroutine TElecConstraint_init(this, input, orb, species0)

    !> Constrain structure instance
    type(TElecConstraint), intent(out) :: this

    !> Input data structure
    type(TElecConstraintInput), intent(inout) :: input

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Chemical species of atoms
    integer, intent(in) :: species0(:)

    integer :: iConstr, nConstr, ii, jj, iAt, nOrb, nSpin

    nConstr = size(input%atomGrp)

    ! We should enable optional initialization of Vc from input at some point.
    allocate(this%Vc(nConstr), source=0.0_dp)

    allocate(this%dWdVc(nConstr), source=0.0_dp)
    allocate(this%deltaW(nConstr), source=0.0_dp)

    call createOptimizer(input%optimiser, nConstr, this%potOpt)

    this%nConstrIter = input%nConstrIter
    this%isConstrConvRequired = input%isConstrConvRequired
    this%constrTol = input%constrTol
    this%Nc = input%atomNc

    ! consistency checks
    if (input%atomSpin(1) /= 1 .and. input%atomSpin(1) /= 2) then
      call error("Invalid spin channel specified.")
    end if

    this%spin = input%atomSpin

    allocate(this%wAt(nConstr))
    allocate(this%wAtOrb(nConstr))

    ! Allocate + initialize arrays and build index mappings
    do iConstr = 1, nConstr
      ! Count orbitals subject to constraints
      nOrb = 0
      do ii = 1, size(input%atomGrp(iConstr)%data)
        iAt = input%atomGrp(iConstr)%data(ii)
        nOrb = nOrb + orb%nOrbAtom(iAt)
      end do
      allocate(this%wAt(iConstr)%data(nOrb), source=0)
      allocate(this%wAtOrb(iConstr)%data(nOrb), source=0)
      nSpin = size(input%atomSpinDir(iConstr)%data)
      nOrb = 0
      do ii = 1, size(input%atomGrp(iConstr)%data)
        iAt = input%atomGrp(iConstr)%data(ii)
        do jj = 1, orb%nOrbAtom(iAt)
          this%wAt(iConstr)%data(nOrb+jj) = iAt
          this%wAtOrb(iConstr)%data(nOrb+jj) = jj
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


  !> Returns total contribution to free energy functional from constraints.
  pure function getFreeEnergy(this) result(deltaWTotal)

    !> Class instance
    class(TElecConstraint), intent(in) :: this

    !> Summed up contribution to free energy functional from constraints
    real(dp) :: deltaWTotal

    deltaWTotal = sum(this%deltaW)

  end function getFreeEnergy


  !> Returns maximum derivative of energy functional with respect to Vc.
  pure function getMaxEnergyDerivWrtVc(this) result(dWdVcMax)

    !> Class instance
    class(TElecConstraint), intent(in) :: this

    !> Maximum derivative of energy functional with respect to Vc
    real(dp) :: dWdVcMax

    dWdVcMax = maxval(abs(this%dWdVc))

  end function getMaxEnergyDerivWrtVc


  !> Applies electronic constraints to system.
  subroutine propagateConstraints(this, qq, energy, tConverged)

    !> Class instance
    class(TElecConstraint), intent(inout) :: this

    !> Mulliken populations
    real(dp), intent(in) :: qq(:,:,:)

    !> Energy
    real(dp), intent(in) :: energy

    !> Gradient convergence achieved
    logical, intent(out) :: tConverged

    ! Summed up contribution to free energy functional from constraints
    real(dp) :: deltaWTotal

    ! Maximum derivative of energy functional with respect to Vc
    real(dp) :: dWdVcMax

    ! Potential displacement proposed by optimizer
    real(dp) :: potDisplace(size(this%Vc))

    integer :: iConstr, nConstr

    nConstr = size(this%wAt)

    do iConstr = 1, nConstr
      call getConstrainEnergyAndPotQ(this%Vc(iConstr), this%Nc(iConstr),&
          & this%wAt(iConstr)%data, this%wAtOrb(iConstr)%data, this%spin(1), qq,&
          & this%deltaW(iConstr), this%dWdVc(iConstr))
    end do

    ! Sum up all free energy contributions
    deltaWTotal = this%getFreeEnergy()

    ! Get maximum derivative of energy functional with respect to Vc
    dWdVcMax = this%getMaxEnergyDerivWrtVc()

    call this%potOpt%step(energy + deltaWTotal, -this%dWdVc, potDisplace)
    this%Vc(:) = this%Vc + potDisplace

    ! In this case dWdVc is equivalent to the condition itself,
    ! so we can use it to measure convergence.
    tConverged = dWdVcMax < this%constrTol

  end subroutine propagateConstraints


  !> Calculate artificial potential to realize constraint on atomic charge.
  subroutine getConstrainEnergyAndPotQ(Vc, Nc, wAt, wOrb, iSpin, qq, deltaW, dWdV)

    !> Potential / Lagrange multiplier
    real(dp), intent(in) :: Vc

    !> Target population
    real(dp), intent(in) :: Nc

    !> Atom(s) involved in current constrain
    integer, intent(in) :: wAt(:)

    !> Orbital(s) involved in current constrain
    integer, intent(in) :: wOrb(:)

    !> Spin channel to be constrained
    integer, intent(in) :: iSpin

    !> Mulliken populations [q, m]
    real(dp), intent(in) :: qq(:,:,:)

    !> Free energy contribution from current contraint
    real(dp), intent(out) :: deltaW

    !> Derivative of free energy with respect to potential
    real(dp), intent(out) :: dWdV

    !> Mulliken populations [up, dn]
    real(dp), allocatable :: qqUpDown(:,:,:)

    ! Present population
    real(dp) :: wn

    ! Index of atomic orbital
    integer :: iW

    ! copy over populations
    qqUpDown = qq
    call qm2ud(qqUpDown)

    wn = 0.0_dp

    do iW = 1, size(wAt)
      wn = wn + qqUpDown(wOrb(iW), wAt(iW), iSpin)
    end do

    dWdV = wn - Nc
    deltaW = Vc * dWdV

  end subroutine getConstrainEnergyAndPotQ


  !> Get total shift of all constraints.
  subroutine getConstrainShift(this, shift)

    !> Class instance
    class(TElecConstraint), intent(inout) :: this

    !> Total shift of all constraints
    real(dp), intent(out) :: shift(:,:,:,:)

    integer :: iConstr, nConstr

    shift(:,:,:,:) = 0.0_dp
    nConstr = size(this%wAt)

    do iConstr = 1, nConstr
      call getConstrainShiftQ(shift, this%Vc(iConstr), this%wAt(iConstr)%data,&
          & this%wAtOrb(iConstr)%data, this%spin(1))
    end do

  end subroutine getConstrainShift


  !> Get shift for atomic charge constraint.
  subroutine getConstrainShiftQ(shift, Vc, wAt, wOrb, iSpin)

    !> Shift to which contribution is appended
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> Potential / Lagrange multiplier
    real(dp), intent(in) :: Vc

    !> Atom(s) involved in current constrain
    integer, intent(in) :: wAt(:)

    !> Orbital(s) involved in current constrain
    integer, intent(in) :: wOrb(:)

    !> Spin channel to be constrained
    integer, intent(in) :: iSpin

    ! Index of atomic orbital
    integer :: iW

    do iW = 1, size(wAt)
      shift(wOrb(iW), wOrb(iW), wAt(iW), iSpin) = shift(wOrb(iW), wOrb(iW), wAt(iW), iSpin) + Vc
    end do

    ! switch to [q, m] representation
    call ud2qm(shift)

  end subroutine getConstrainShiftQ

end module dftbp_dftb_elecconstraints
