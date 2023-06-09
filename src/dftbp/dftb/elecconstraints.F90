!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Constraints on the electronic ground state
module dftbp_dftb_elecconstraints
  use dftbp_common_accuracy, only : dp, sc
  use dftbp_math_angmomentum, only : getLOperators
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_orbitals, only : orbitalNames, getShellNames
  use dftbp_math_bisect, only : bisection
  use dftbp_io_charmanip, only : i2c
  implicit none

  private

  public :: constrainQ, constrainS, constrainL, constrainJ, constrainMj
  public :: printHighestAO


contains


  !> Prints atom index, angular and magnetic quantum number of occupied AO with highest energy.
  subroutine printHighestAO(HSqrReal, orb, denseDesc, iAtInCentralRegion, species, speciesName, qq)

    !> Square (unpacked) Hamiltonian of certain spin channel
    real(dp), intent(in) :: HSqrReal(:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> List of atoms in the central cell (or device region if transport)
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Species of atoms
    integer, intent(in) :: species(:)

    !> Name of each species
    character(*), intent(in) :: speciesName(:)

    !> Mulliken orbital charges on output (mOrb, nAtom, nSpin)
    real(dp), intent(in) :: qq(:,:,:)

    character(sc), allocatable :: shellNamesTmp(:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Stores start/end index and number of orbitals of square matrices
    integer :: desc(descLen)

    !! Number of orbitals in square matrices
    integer :: nOrb

    integer :: ll, iAt
    integer :: ii, iSp, iSh, iBlockStart, iBlockEnd, iDiag
    real(dp) :: occ, energy
    real(dp), allocatable :: hamBlock(:,:)
    character(len=:), allocatable :: fmt

    print "(/,A)", '############### AO Analysis ###############'

    ! Check if we are dealing with distributed matrices
    @:ASSERT(size(HSqrReal, dim=1) == size(HSqrReal, dim=2))

    nOrb = size(HSqrReal, dim=1)

    do ii = 1, size(iAtInCentralRegion)
      iAt = iAtInCentralRegion(ii)
      desc = getDescriptor(iAt, denseDesc%iAtomStart)
      iSp = species(iAt)
      call getShellNames(iSp, orb, shellNamesTmp)
      do iSh = 1, orb%nShell(iSp)
        ll = orb%angShell(iSh, iSp)
        occ = sum(qq(orb%posShell(iSh, iSp):orb%posShell(iSh + 1, iSp) - 1, iAt, 1))
        if (occ > 1.0_dp) then
          iBlockStart = desc(iStart) + orb%posShell(iSh, iSp) - 1
          iBlockEnd = desc(iStart) + orb%posShell(iSh + 1, iSp) - 2
          hamBlock = HSqrReal(iBlockStart:iBlockEnd, iBlockStart:iBlockEnd)
          energy = 0.0_dp
          do iDiag = 1, size(hamBlock, dim=1)
            energy = energy + hamBlock(iDiag, iDiag)
          end do
          energy = energy / size(hamBlock, dim=1)
          print "(A,I0,4A,A,2E26.16)", "iAtom, species, shell, occupation, energy: ", iAt, ", ",&
              & trim(speciesName(iSp)), ", ", trim(shellNamesTmp(ll + 1)), ", ", occ, energy
          fmt = "(" // i2c(size(hamBlock, dim=1)) // "F16.10)"
          ! print fmt, transpose(hamBlock)
        end if
      end do
      deallocate(shellNamesTmp)
    end do

    ! diagMax = minval(HSqrReal)

    ! do iOrb = 1, nOrb
    !   call bisection(iAt, denseDesc%iAtomStart, iOrb)
    !   iOrbStart = denseDesc%iAtomStart(iAt)
    !   iOrbOnAtom = iOrb - iOrbStart + 1

    !   diag = HSqrReal(iOrb, iOrb)
    !   ! occ = qq(orb%posShell(iSh, iSp):orb%posShell(iSh + 1, iSp) - 1, iAt, 1)
    !   occ = 1.0_dp

    !   if (occ > 1.0e-08_dp) then
    !     if (diag > diagMax) then
    !       diagMax = diag
    !       iOrbMax = iOrb
    !     end if
    !   end if
    ! end do

    ! ! get atom index from orbital index
    ! call bisection(iAtMax, denseDesc%iAtomStart, iOrbMax)
    ! print *, iAtMax, iOrbMax

    ! ! get species index from atom index
    ! iSpMax = species(iAtMax)
    ! print *, orb%iShellOrb(:, iSpMax)
    ! call getShellNames(iSpMax, orb, shellNamesTmp)

    ! ! first orbital in global indexing on atom
    ! iOrbStart = denseDesc%iAtomStart(iAtMax)

    ! ! orbital index of related to the associated atom
    ! iOrbOnAtom = iOrbMax - iOrbStart + 1

    ! ll = orb%iShellOrb(iOrbOnAtom, iSpMax) - 1
    ! mm = iOrbOnAtom - orb%posShell(ll + 1, iSpMax) - ll

    ! print "(A,E20.12,A,/)", "Found highest (occupied) AO at: ", diagMax, " Ha"
    ! print "(A,I0,5A)", "Orbital-index/-type: ", iOrbMax, " (", trim(shellNamesTmp(ll + 1)),&
    !     & "_", trim(orbitalNames(mm, ll)), ")"
    ! print "(3A,I0,A)", "Atom-type/-index: ", trim(speciesName(iSpMax)), " (", iAtMax, ")"
    ! print "(A,E20.12)", "Occupation: ", qq(iOrbOnAtom, iAtMax, 1)

    print "(A,/)", '###############     END     ###############'

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


  !> Quadratic constraint on atomic charge
  subroutine constrainQ(shift, qIn, orb, species, conAt, conSh, Qtarget, V)

    !> shift to append contribution
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> charges
    real(dp), intent(in) :: qIn(:,:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Chemical species of atoms
    integer, intent(in) :: species(:)

    !> atom to be constrained
    integer, intent(in) :: conAt

    !> shell of atom to be constrained
    integer, intent(in) :: conSh

    !> target value
    real(dp), intent(in) :: Qtarget

    !> weight of the constraint
    real(dp), intent(in) :: V

    integer :: iOrb
    real(dp) :: Qshell

    Qshell = sum(qIn(orb%posShell(conSh,species(conAt)):orb%posShell(conSh+1,species(conAt))-1, &
        & conAt,1))

    ! Push q towards required value
    do iOrb = orb%posShell(conSh,species(conAt)), orb%posShell(conSh+1,species(conAt))-1
      shift(iOrb,iOrb,conAt,1) = shift(iOrb,iOrb,conAt,1) + V * 0.5_dp*(Qshell - Qtarget)
    end do

  end subroutine constrainQ


  !> Quadratic constraint on local spin (non-colinear)
  subroutine constrainS(shift, qIn, orb, species, conAt, conSh, Starget, V, vec)

    !> shift to append contribution
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> charges
    real(dp), intent(in) :: qIn(:,:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Chemical species of atoms
    integer, intent(in) :: species(:)

    !> atom to be constrained
    integer, intent(in) :: conAt

    !> shell of atom to be constrained
    integer, intent(in) :: conSh

    !> target value
    real(dp), intent(in) :: Starget

    !> weight of the constraint
    real(dp), intent(in) :: V

    !> direction of spin
    real(dp), intent(in) :: vec(3)

    integer :: iOrb, nSpin, iSpin
    real(dp) :: Sshell(3), W, vecNorm(3)

    nSpin = size(shift,dim=4)

    vecNorm = vec / sqrt(sum(vec**2))

    Sshell = sum(qIn(orb%posShell(conSh,species(conAt)):orb%posShell(conSh+1,species(conAt))-1, &
        & conAt,2:4),dim=1)

    if (sqrt(sum(Sshell**2)) < 1.0E-8_dp) then
      Sshell = Sshell + 1.0E-8_dp*(/1,1,1/)
    end if

    vecNorm = Sshell  / sqrt(sum(Sshell**2))

    ! Push S towards required value
    w = V * 0.5_dp*(dot_product(Sshell,vecNorm) - Starget)

    do iSpin = 2, nSpin
      do iOrb = orb%posShell(conSh,species(conAt)),orb%posShell(conSh+1,species(conAt))-1
        shift(iOrb,iOrb,conAt,iSpin) = shift(iOrb,iOrb,conAt,iSpin) + w * vecNorm(iSpin-1)
      end do
    end do

  end subroutine constrainS


  !> Quadratic constraint on orbital angular momentum
  subroutine constrainL(iShift,qBlockSkew, orb, species, conAt, conSh, Ltarget, V, vec)

    !> shift block shift
    real(dp), intent(inout) :: iShift(:,:,:,:)

    !> Antisymmetric Mulliken block populations for imaginary coefficients of
  !> Pauli matrics
    real(dp), intent(in) :: qBlockSkew(:,:,:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Atom for constraint
    integer, intent(in) :: conAt

    !> Shell for constraint
    integer, intent(in) :: conSh

    !> value of L
    real(dp), intent(in) :: Ltarget

    !> strength of constraint
    real(dp), intent(in) :: V

    !> direction of constrain
    real(dp), intent(in) :: vec(3)

    integer :: ii, iSp, iSh, iOrb, iStart, iEnd
    real(dp), allocatable :: SpeciesL(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: tmpBlock(:,:)
    real(dp) :: Lshell(3), W, vecNorm(3)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    vecNorm = vec / sqrt(sum(vec**2))

    allocate(SpeciesL(orb%mOrb,orb%mOrb,3))
    SpeciesL = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))

    iSp = species(conAt)
    Lz = 0.0_dp
    Lplus = 0.0_dp
    iSh = orb%angShell(conSh,iSp)
    call getLoperators(iSh, Lplus(1:2*iSh+1,1:2*iSh+1),Lz(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,1) = aimag(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,2) = -real(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,3) = aimag(Lz(1:2*iSh+1,1:2*iSh+1))

    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    Lshell = 0.0_dp

    iOrb = orb%nOrbSpecies(iSp)
    tmpBlock(:,:) = 0.0_dp
    ! identity part
    tmpBlock(1:iOrb,1:iOrb) = qBlockSkew(1:iOrb,1:iOrb,conAt,1)
    iStart = orb%posShell(conSh,iSp)
    iEnd = orb%posShell(conSh+1,iSp)-1
    do ii = 1, 3
      Lshell(ii) = &
          & - sum(SpeciesL(iStart:iEnd,iStart:iEnd,ii) &
          &  * transpose(tmpBlock(iStart:iEnd,iStart:iEnd)))
    end do

    if (sqrt(sum(Lshell**2)) < 1.0E-8_dp) then
      Lshell = Lshell + 1.0E-8_dp*(/1,1,1/)
    end if

    vecNorm = Lshell / sqrt(sum(Lshell**2))

    ! Push L towards required value
    w = V * 0.5_dp*(dot_product(lshell,vecNorm)-Ltarget)

    do ii = 1, 3
      iShift(iStart:iEnd,iStart:iEnd,conAt,1) = &
          & iShift(iStart:iEnd,iStart:iEnd,conAt,1) &
          & + w * vecNorm(ii) * SpeciesL(iStart:iEnd,iStart:iEnd,ii)
    end do

  end subroutine constrainL


  !> Quadratic constraint on total angular momentum
  subroutine constrainJ(shift, qIn, iShift, qBlockSkew, orb, species, conAt, conSh, Jtarget, V, &
      & vec)
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> charges
    real(dp), intent(in) :: qIn(:,:,:)

    !> Imaginary block shift
    real(dp), intent(inout) :: iShift(:,:,:,:)

    !> Antisymmetric Mulliken block populations for imaginary coefficients of Pauli matrics
    real(dp), intent(in) :: qBlockSkew(:,:,:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Atom for constraint
    integer, intent(in) :: conAt

    !> Shell for constraint
    integer, intent(in) :: conSh

    !> value of J
    real(dp), intent(in) :: Jtarget

    !> strength of constraint
    real(dp), intent(in) :: V

    !> direction of constrain
    real(dp), intent(in) :: vec(3)

    integer :: ii, iSp, iSh, iOrb, iStart, iEnd, nSpin, iSpin
    real(dp), allocatable :: SpeciesL(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: tmpBlock(:,:)
    real(dp) :: Lshell(3), Sshell(3), W, vecNorm(3)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    nSpin = size(shift,dim=4)

    vecNorm = vec / sqrt(sum(vec**2))

    allocate(SpeciesL(orb%mOrb,orb%mOrb,3))
    SpeciesL = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))

    iSp = species(conAt)
    Lz = 0.0_dp
    Lplus = 0.0_dp
    iSh = orb%angShell(conSh,iSp)
    call getLoperators(iSh, Lplus(1:2*iSh+1,1:2*iSh+1),Lz(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,1) = aimag(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,2) = -real(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,3) = aimag(Lz(1:2*iSh+1,1:2*iSh+1))

    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    Lshell = 0.0_dp

    iOrb = orb%nOrbSpecies(iSp)
    tmpBlock(:,:) = 0.0_dp

    ! identity part
    tmpBlock(1:iOrb,1:iOrb) = qBlockSkew(1:iOrb,1:iOrb,conAt,1)
    iStart = orb%posShell(conSh,iSp)
    iEnd = orb%posShell(conSh+1,iSp)-1
    do ii = 1, 3
      Lshell(ii) = -sum(SpeciesL(iStart:iEnd,iStart:iEnd,ii) &
          & * transpose(tmpBlock(iStart:iEnd,iStart:iEnd)))
    end do

    Sshell = sum(qIn(orb%posShell(conSh,species(conAt)): &
        & orb%posShell(conSh+1,species(conAt))-1,conAt,2:4),dim=1)

    if ( sqrt(sum((lshell + 0.5_dp*Sshell)**2)) < 1.0E-8_dp) then
      Sshell = Sshell + 1.0E-8_dp*(/1,1,1/)
    end if

    vecNorm = (lshell + 0.5_dp*Sshell) / sqrt(sum((lshell + 0.5_dp*Sshell)**2))

    ! Push J towards required value
    w = V * 0.5_dp*(dot_product(lshell,vecNorm)+ &
        & 0.5_dp*dot_product(Sshell,vecNorm) -Jtarget)

    do ii = 1, 3
      iShift(iStart:iEnd,iStart:iEnd,conAt,1) = &
          & iShift(iStart:iEnd,iStart:iEnd,conAt,1) &
          & + w * vecNorm(ii) * SpeciesL(iStart:iEnd,iStart:iEnd,ii)
    end do

    do iSpin = 2, nSpin
      do iOrb = orb%posShell(conSh,species(conAt)), &
          & orb%posShell(conSh+1,species(conAt))-1
        shift(iOrb,iOrb,conAt,iSpin) = shift(iOrb,iOrb,conAt,iSpin) &
            & + w * vecNorm(iSpin-1)
      end do
    end do

  end subroutine constrainJ


  !> Quadratic constraint on projection of angular momentum
  subroutine constrainMj(shift, qIn, iShift, qBlockSkew, orb, species, conAt, conSh, MjTarget, V, &
      & vec)

    !> block shift
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> charges
    real(dp), intent(in) :: qIn(:,:,:)

    !> Imaginary block shift
    real(dp), intent(inout) :: iShift(:,:,:,:)

    !> Antisymmetric Mulliken block populations for imaginary coefficients of Pauli matrics
    real(dp), intent(in) :: qBlockSkew(:,:,:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Atom for constraint
    integer, intent(in) :: conAt

    !> Shell for constraint
    integer, intent(in) :: conSh

    !> value of Mj
    real(dp), intent(in) :: MjTarget

    !> strength of constraint
    real(dp), intent(in) :: V

    !> direction of constrain
    real(dp), intent(in) :: vec(3)

    integer :: ii, iSp, iSh, iOrb, iStart, iEnd, nSpin, iSpin
    real(dp), allocatable :: SpeciesL(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: tmpBlock(:,:)
    real(dp) :: Lshell(3), Sshell(3), W, vecNorm(3)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    nSpin = size(shift,dim=4)

    vecNorm = vec / sqrt(sum(vec**2))

    allocate(SpeciesL(orb%mOrb,orb%mOrb,3))
    SpeciesL = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))

    iSp = species(conAt)
    Lz = 0.0_dp
    Lplus = 0.0_dp
    iSh = orb%angShell(conSh,iSp)
    call getLOperators(iSh, Lplus(1:2*iSh+1,1:2*iSh+1),Lz(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,1) = aimag(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,2) = -real(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,3) = aimag(Lz(1:2*iSh+1,1:2*iSh+1))

    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    Lshell = 0.0_dp

    iOrb = orb%nOrbSpecies(iSp)
    tmpBlock(:,:) = 0.0_dp

    ! identity part
    tmpBlock(1:iOrb,1:iOrb) = qBlockSkew(1:iOrb,1:iOrb,conAt,1)
    iStart = orb%posShell(conSh,iSp)
    iEnd = orb%posShell(conSh+1,iSp)-1
    do ii = 1, 3
      Lshell(ii) = &
          & - sum(SpeciesL(iStart:iEnd,iStart:iEnd,ii) &
          &  * transpose(tmpBlock(iStart:iEnd,iStart:iEnd)))
    end do

    Sshell = sum(qIn(orb%posShell(conSh,species(conAt)): &
        & orb%posShell(conSh+1,species(conAt))-1,conAt,2:4),dim=1)

    ! Push Mj towards required value
    w = V * 0.5_dp*(dot_product(lshell,vecNorm)+ &
        & 0.5_dp*dot_product(Sshell,vecNorm)-MjTarget)

    do ii = 1, 3
      iShift(iStart:iEnd,iStart:iEnd,conAt,1) = &
          & iShift(iStart:iEnd,iStart:iEnd,conAt,1) &
          & + w * vecNorm(ii) * SpeciesL(iStart:iEnd,iStart:iEnd,ii)
    end do

    do iSpin = 2, nSpin
      do iOrb = orb%posShell(conSh,species(conAt)), &
          & orb%posShell(conSh+1,species(conAt))-1
        shift(iOrb,iOrb,conAt,iSpin) = shift(iOrb,iOrb,conAt,iSpin) &
            & + w * vecNorm(iSpin-1)
      end do
    end do

  end subroutine constrainMj

end module dftbp_dftb_elecconstraints
