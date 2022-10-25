!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Calculate either the whole single particle density matrix and energy
!> weighted density matrix or only the elements dictated by a neighbour map.
!> Calculation of the whole matrix scales as O(N**3), the sparse form as
!> O(N**2) but with a larger pre-factor.
!> Note: Dense code based on suggestions from Thomas Heine
!> Caveat: The routines create the transposed and complex conjugated of the density matrices! (cc*
!> instead of the conventional c*c)
module dftbp_dftb_densitymatrix

  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_math_blasroutines, only : herk
  use dftbp_math_sorting, only : unique, heap_sort
  use dftbp_type_commontypes, only : TOrbitals, TParallelKS
  use dftbp_common_globalenv, only : stdOut
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : blacsgrid, blocklist, size, pblasfx_pgemm, pblasfx_ptran,&
      & pblasfx_ptranc
#:endif

  implicit none
  private

  public :: makeDensityMatrix, TDensityMatrix
  public :: transformDualSpaceToBvKRealSpace, transformBvKRealSpaceToDualSpace

#:if WITH_SCALAPACK
  public :: makeDensityMtxRealBlacs, makeDensityMtxCplxBlacs
#:endif


  !> Provides an interface to calculate the two types of dm - regular and weighted and put them into
  !! packed storage.
  interface makeDensityMatrix
    module procedure fullDensityMatrix_real
    module procedure fullDensityMatrix_cmplx
    module procedure fullEnergyDensityMatrix_real
    module procedure fullEnergyDensityMatrix_cmplx
  end interface makeDensityMatrix


  !> Holds real and complex delta density matrices and pointers.
  type :: TDensityMatrix

    !> DeltaRho input for calculation of range separated Hamiltonian
    real(dp), allocatable :: deltaRhoIn(:)

    !> DeltaRho output from calculation of range separated Hamiltonian
    real(dp), allocatable :: deltaRhoOut(:)

    !> DeltaRho output from calculation of range separated Hamiltonian
    complex(dp), allocatable :: deltaRhoOutCplx(:)

    !> Holds change in deltaRho between SCC steps for range separation
    real(dp), allocatable :: deltaRhoDiff(:)

    !> DeltaRho input for range separation in matrix form
    real(dp), pointer :: deltaRhoInSqr(:,:,:) => null()

    !> DeltaRho output from range separation in matrix form
    real(dp), pointer :: deltaRhoOutSqr(:,:,:) => null()

    !> Complex, square dual-space deltaRho output from range separation
    complex(dp), pointer :: deltaRhoOutSqrCplx(:,:,:) => null()

    !> Real-space, square deltaRho input for range separation
    real(dp), pointer :: deltaRhoInSqrCplxHS(:,:,:,:,:,:) => null()

    !> Real-space, square deltaRho output for range separation
    real(dp), pointer :: deltaRhoOutSqrCplxHS(:,:,:,:,:,:) => null()

  end type TDensityMatrix


  real(dp), parameter :: arbitraryConstant = 0.1_dp

contains


  !> Transforms dense, square density matrix for all spins/k-points to real-space (BvK cell).
  subroutine transformDualSpaceToBvKRealSpace(rhoSqrDual, parallelKS, kPoint, kWeight, bvKShifts,&
      & coeffsDiag, rhoSqrBvK)

    !> Complex, dense, square dual-space rho of all spins/k-points
    complex(dp), intent(in), pointer :: rhoSqrDual(:,:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> k-points in relative units
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights of k-points
    real(dp), intent(in) :: kWeight(:)

    !> K-point compatible BvK real-space shifts in relative coordinates (units of latVecs)
    real(dp), intent(in) :: bvKShifts(:,:)

    !> Supercell folding coefficients (diagonal elements)
    integer, intent(in) :: coeffsDiag(:)

    !> Real-space, dense, square rho for BvK cell
    real(dp), intent(in), pointer :: rhoSqrBvK(:,:,:,:,:,:)

    !> Real-space, dense, square rho for BvK cell (complex version)
    ! complex(dp) :: rhoSqrBvKCplx(size(rhoSqrBvK, dim=1), size(rhoSqrBvK, dim=2),&
    !     & size(rhoSqrBvK, dim=3), size(rhoSqrBvK, dim=4), size(rhoSqrBvK, dim=5),&
    !     & size(rhoSqrBvK, dim=6))
    complex(dp), allocatable :: rhoSqrBvKCplx(:,:,:,:,:,:)

    !! K-point-spin composite index and k-point/spin index
    integer :: iKS, iK, iSpin

    !! Iterates over all BvK real-space vectors
    integer :: iG

    !! Phase factor
    complex(dp) :: phase

    !! Integer BvK real-space shift translated to density matrix indices
    integer :: bvKIndex(3)

    allocate(rhoSqrBvKCplx(size(rhoSqrBvK, dim=1), size(rhoSqrBvK, dim=2),&
        & size(rhoSqrBvK, dim=3), size(rhoSqrBvK, dim=4), size(rhoSqrBvK, dim=5),&
        & size(rhoSqrBvK, dim=6)))

    rhoSqrBvKCplx(:,:,:,:,:,:) = (0.0_dp, 0.0_dp)

    do iG = 1, size(bvKShifts, dim=2)
      bvKIndex(:) = nint(bvKShifts(:, iG)) + 1
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iSpin = parallelKS%localKS(2, iKS)
        phase = exp(cmplx(0, -1, dp) * dot_product(2.0_dp * pi * kPoint(:, iK), bvKShifts(:, iG)))
        rhoSqrBvKCplx(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3), iSpin)&
            & = rhoSqrBvKCplx(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3), iSpin)&
            & + kWeight(iK) * rhoSqrDual(:,:, iKS) * phase
      end do
    end do

    rhoSqrBvK(:,:,:,:,:,:) = rhoSqrBvK + real(rhoSqrBvKCplx, dp)

  end subroutine transformDualSpaceToBvKRealSpace


  !> Transforms dense, square real-space density matrix (in BvK cell) to dual-space.
  subroutine transformBvKRealSpaceToDualSpace(rhoSqrBvK, parallelKS, kPoint, kWeight, bvKShifts,&
      & coeffsDiag, rhoSqrDual)

    !> Real-space, dense, square rho for BvK cell
    real(dp), intent(in), pointer :: rhoSqrBvK(:,:,:,:,:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> k-points in relative units
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights of k-points
    real(dp), intent(in) :: kWeight(:)

    !> K-point compatible BvK real-space shifts in relative coordinates (units of latVecs)
    real(dp), intent(in) :: bvKShifts(:,:)

    !> Supercell folding coefficients (diagonal elements)
    integer, intent(in) :: coeffsDiag(:)

    !> Complex, dense, square dual-space rho of all spins/k-points
    complex(dp), intent(in), pointer :: rhoSqrDual(:,:,:)

    !! K-point-spin composite index and k-point/spin index
    integer :: iKS, iK, iSpin

    !! Iterates over all BvK real-space vectors
    integer :: iG

    !! Phase factor
    complex(dp) :: phase

    !! Integer BvK real-space shift translated to density matrix indices
    integer :: bvKIndex(3)

    rhoSqrDual(:,:,:) = 0.0_dp

    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iSpin = parallelKS%localKS(2, iKS)
      do iG = 1, size(bvKShifts, dim=2)
        phase = exp(cmplx(0, 1, dp) * dot_product(2.0_dp * pi * kPoint(:, iK), bvKShifts(:, iG)))
        bvKIndex(:) = nint(bvKShifts(:, iG)) + 1
        rhoSqrDual(:,:, iKS) = rhoSqrDual(:,:, iKS)&
            & + rhoSqrBvK(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3), iSpin) * phase
      end do
    end do

  end subroutine transformBvKRealSpaceToDualSpace


  !> Make a regular density matrix for the real wave-function case
  !> Note: In order to save memory, the eigenvectors (which should be intent(in) parameters) are
  !> overwritten and then restored again
  subroutine fullDensityMatrix_real(dm, eigenvecs, filling)

    !> the resulting nOrb*nOrb density matrix
    real(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    integer :: ii, nLevels
    real(dp) :: shift

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))

    dm(:,:) = 0.0_dp
    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do
    shift = minval(filling(1:nLevels))
    if (shift > epsilon(1.0_dp)) then
      ! all fillings are definitely positive

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling operations are positive
      call herk(dm, eigenvecs(:,1:nLevels))
      shift = shift - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)-shift) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), beta=shift)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if
  end subroutine fullDensityMatrix_real


  !> Make a regular density matrix for the complex wave-function case
  !> Note: in order to save memory, the eigenvectors (which should be intent(in) parameters) are
  !> overwritten and then restored again
  subroutine fullDensityMatrix_cmplx(dm, eigenvecs, filling)

    !> the resulting nOrb*nOrb density matrix
    complex(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    complex(dp), intent(inout) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    integer :: ii, nLevels
    real(dp) :: shift

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))

    dm(:,:) = cmplx(0.0_dp,0.0_dp,dp)

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do
    shift = minval(filling(1:nLevels))
    if (shift > epsilon(1.0_dp)) then
      ! all fillings are definitely positive
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii))*eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii))
      end do
      !$OMP  END PARALLEL DO

    else
      ! shift matrix so that filling operations are positive
      call herk(dm, eigenvecs(:,1:nLevels))
      shift = shift - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)-shift)*eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), beta=shift)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if
  end subroutine fullDensityMatrix_cmplx


  !> Make an energy weighted density matrix for the real wave-function case
  !> Note: in order to save memory, the eigenvectors (which should be intent(in) parameters) are
  !> overwritten and then restored again
  subroutine fullEnergyDensityMatrix_real(dm, eigenvecs, filling, eigen)

    !> the resulting nOrb*nOrb density matrix
    real(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> eigenvalues of the system
    real(dp), intent(in) :: eigen(:)

    integer :: ii, nLevels
    real(dp) :: shift
    real(dp) :: fillProduct(size(filling))

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    dm(:,:) = 0.0_dp
    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do
    fillProduct(1:nLevels) = filling(1:nLevels) * eigen(1:nLevels)
    if ((minval(fillProduct(1:nLevels)) < 0.0_dp&
        & .eqv. maxval(fillProduct(1:nLevels)) < 0.0_dp)&
        & .and. abs(minval(fillProduct(1:nLevels))) > epsilon(1.0_dp)&
        & .and. abs(maxval(fillProduct(1:nLevels))) > epsilon(1.0_dp)) then
      ! all fillings the same sign, and fairly large

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(abs(fillProduct(ii)))*eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), alpha=sign(1.0_dp, maxval(fillProduct(1:nLevels))))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(abs(fillProduct(ii)))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling operations are positive
      call herk(dm, eigenvecs(:,1:nLevels))
      shift = minval(fillProduct(1:nLevels)) - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(fillProduct(ii)-shift) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), beta=shift)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(fillProduct(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if
  end subroutine fullEnergyDensityMatrix_real


  !> Make an energy weighted density matrix for the complex wave-function case
  !> Note: in order to save memory, the eigenvectors (which should be intent(in) parameters) are
  !> overwritten and then restored again
  subroutine fullEnergyDensityMatrix_cmplx(dm, eigenvecs, filling, eigen)

    !> the resulting nOrb*nOrb density matrix
    complex(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    complex(dp), intent(inout) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> eigenvalues of the system
    real(dp), intent(in) :: eigen(:)

    integer :: ii, nLevels
    real(dp) :: shift
    real(dp) :: fillProduct(size(filling))

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    dm(:,:) = cmplx(0.0_dp,0.0_dp,dp)

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    fillProduct(1:nLevels) = filling(1:nLevels) * eigen(1:nLevels)
    if ((minval(fillProduct(1:nLevels)) < 0.0_dp&
        & .eqv. maxval(fillProduct(1:nLevels)) < 0.0_dp)&
        & .and. abs(minval(fillProduct(1:nLevels))) > epsilon(1.0_dp)&
        & .and. abs(maxval(fillProduct(1:nLevels))) > epsilon(1.0_dp)) then
      ! all fillings the same sign, and fairly large
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(abs(fillProduct(ii))) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), alpha=sign(1.0_dp, maxval(fillProduct(1:nLevels))))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(abs(fillProduct(ii)))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling operations are positive
      call herk(dm, eigenvecs(:,1:nLevels))
      shift = minval(fillProduct(1:nLevels)) - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(fillProduct(ii)-shift)*eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), beta=shift)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(fillProduct(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if
  end subroutine fullEnergyDensityMatrix_cmplx


#:if WITH_SCALAPACK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Scalapack routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create density or energy weighted density matrix (real) for both triangles.
  subroutine makeDensityMtxRealBlacs(myBlacs, desc, filling, eigenVecs, densityMtx, eigenVals)

    !> BLACS grid information
    type(blacsgrid), intent(in) :: myBlacs

    !> Matrix descriptor
    integer, intent(in) :: desc(:)

    !> Occupations of levels
    real(dp), target, intent(in) :: filling(:)

    !> Eigenvectors of system
    real(dp), intent(inout) :: eigenVecs(:,:)

    !> Resulting (symmetric) density matrix
    real(dp), intent(out) :: densityMtx(:,:)

    !> Eigenvalues, if energy weighted density matrix required
    real(dp), intent(in), optional :: eigenVals(:)

    integer  :: ii, jj, iGlob, iLoc, blockSize
    type(blocklist) :: blocks
    real(dp), allocatable :: work(:,:)

    densityMtx(:, :) = 0.0_dp
    work = densityMtx

    ! Scale a copy of the eigenvectors
    call blocks%init(myBlacs, desc, "c")
    if (present(eigenVals)) then
      do ii = 1, size(blocks)
        call blocks%getblock(ii, iGlob, iLoc, blockSize)
        do jj = 0, blockSize - 1
          work(:, iLoc + jj) = eigenVecs(:, iLoc + jj) * eigenVals(iGlob + jj) * filling(iGlob + jj)
        end do
      end do
    else
      do ii = 1, size(blocks)
        call blocks%getblock(ii, iGlob, iLoc, blockSize)
        do jj = 0, blockSize - 1
          work(:, iLoc + jj) = eigenVecs(:, iLoc + jj) * filling(iGlob + jj)
        end do
      end do
    end if

    ! Create matrix
    call pblasfx_pgemm(eigenVecs, desc, work, desc, densityMtx, desc, transb="T")
    ! symmetrize
    work(:,:) = densityMtx
    call pblasfx_ptran(work, desc, densityMtx, desc, alpha=0.5_dp, beta=0.5_dp)

  end subroutine makeDensityMtxRealBlacs


  !> Create density or energy weighted density matrix (complex) for both triangles.
  subroutine makeDensityMtxCplxBlacs(myBlacs, desc, filling, eigenVecs, densityMtx, eigenVals)

    !> BLACS grid information
    type(blacsgrid), intent(in) :: myBlacs

    !> Matrix descriptor
    integer, intent(in) :: desc(:)

    !> Occupations of levels
    real(dp), target, intent(in) :: filling(:)

    !> Eigenvectors of system
    complex(dp), intent(inout) :: eigenVecs(:,:)

    !> Resulting (symmetric) density matrix
    complex(dp), intent(out) :: densityMtx(:,:)

    !> Eigenvalues, if energy weighted density matrix required
    real(dp), intent(in), optional :: eigenVals(:)

    integer  :: ii, jj, iGlob, iLoc, blockSize
    type(blocklist) :: blocks
    complex(dp), allocatable :: work(:,:)

    densityMtx(:, :) = cmplx(0,0,dp)
    work = densityMtx

    ! Scale a copy of the eigenvectors
    call blocks%init(myBlacs, desc, "c")
    if (present(eigenVals)) then
      do ii = 1, size(blocks)
        call blocks%getblock(ii, iGlob, iLoc, blockSize)
        do jj = 0, blockSize - 1
          work(:, iLoc + jj) = eigenVecs(:, iLoc + jj) * eigenVals(iGlob + jj) * filling(iGlob + jj)
        end do
      end do
    else
      do ii = 1, size(blocks)
        call blocks%getblock(ii, iGlob, iLoc, blockSize)
        do jj = 0, blockSize - 1
          work(:, iLoc + jj) = eigenVecs(:, iLoc + jj) * filling(iGlob + jj)
        end do
      end do
    end if

    ! Create matrix
    call pblasfx_pgemm(eigenVecs, desc, work, desc, densityMtx, desc, transb="C")
    ! hermitian symmetrize
    work(:,:) = densityMtx
    call pblasfx_ptranc(work, desc, densityMtx, desc, alpha=(0.5_dp,0.0_dp), beta=(0.5_dp,0.0_dp))

  end subroutine makeDensityMtxCplxBlacs

#:endif

end module dftbp_dftb_densitymatrix
