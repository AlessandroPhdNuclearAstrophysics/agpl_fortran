!> \file integration.f90
!! \brief Numerical integration utilities for variational scattering codes.
!! \defgroup integration Numerical Integration
!! \ingroup math
!!
!! This module provides routines for block-adaptive integration, Gauss-Laguerre quadrature,
!! and exponentially growing grids, with OpenMP support for parallelism.
!!
!! \author Alessandro
!! \date 2025
!!
!! \note This module must be compiled with OpenMP support enabled.
MODULE INTEGRATION_MOD
  USE OMP_LIB ! OpenMP library for parallel processing
CONTAINS

  !> \brief Block-adaptive numerical integration over up to 3 blocks.
  !! \ingroup integration
  !! \param[in] JB    Number of blocks (1-3)
  !! \param[in] M1,M2,M3 Number of steps in each block
  !! \param[in] H1,H2,H3 Step sizes for each block
  !! \param[in] A     Data array to integrate
  !! \param[in] IAS   Starting index
  !! \return   I5     Integrated result
  FUNCTION B5(JB, M1, M2, M3, H1, H2, H3, A, IAS) RESULT(I5)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: JB          ! Number of blocks (1 to 3)
    INTEGER, INTENT(IN) :: M1, M2, M3  ! Steps in blocks 1,2,3
    INTEGER, INTENT(IN) :: IAS         ! Starting index into A
    DOUBLE PRECISION, INTENT(IN) :: H1, H2, H3 ! Step sizes for each block
    DOUBLE PRECISION, INTENT(IN) :: A(*)       ! Data array (size >= M1+M2+M3+3)
    DOUBLE PRECISION :: I5                     ! Integrated result

    INTEGER :: BLOCK, N, NB, N1, L0, I, L_base
    DOUBLE PRECISION :: H, HL5, S2, S3, S4, AX
    INTEGER :: IDX

    ! Validate block count
    IF (JB < 1 .OR. JB > 3) THEN
      WRITE(*,*) "Error: JB must be 1, 2, or 3."
      STOP
    END IF

    I5 = 0.D0

    ! Loop over blocks serially
    DO BLOCK = 1, JB

      ! Select block parameters
      SELECT CASE (BLOCK)
      CASE (1)
        N   = M1
        H   = H1
        L0  = IAS
      CASE (2)
        N   = M2
        H   = H2
        L0  = IAS + M1
      CASE (3)
        N   = M3
        H   = H3
        L0  = IAS + M1 + M2
      END SELECT

      ! Require at least 4 intervals per block
      IF (N < 4) THEN
        WRITE(*,*) "Error: Block ", BLOCK, " must have at least 4 intervals."
        STOP
      END IF

      ! Precompute
      HL5 = H / 32.D0
      NB  = N / 4
      N1  = N - 4*NB + 1
      L_base = L0 ! base index for leftover calc

      ! Initialize sums
      S2 = 0.D0
      S3 = 0.D0
      S4 = 0.D0

      ! Parallel accumulate group sums
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) REDUCTION(+:S2,S3,S4)
        !$OMP DO SCHEDULE(GUIDED,16)
        DO I = 1, NB
          idx = L0 + 4*(I-1)
          S2 = S2 + A(idx+2)
          S3 = S3 + A(idx+1) + A(idx+3)
          S4 = S4 + A(idx)   + A(idx+4)
        END DO
        !$OMP END DO
      !$OMP END PARALLEL

      ! Compute starting index for leftover points
      L_base = L0 + 4*NB

      ! Correction for leftover points
      SELECT CASE (N1)
      CASE (1)
        AX = 0.D0
      CASE (2)
        AX = HL5 * (-19.D0*A(L_base-3) +106.D0*A(L_base-2) -264.D0*A(L_base-1) + &
                    646.D0*A(L_base)   +251.D0*A(L_base+1))
      CASE (3)
        AX = HL5 * (-8.D0*A(L_base-2)  + 32.D0*A(L_base-1) +192.D0*A(L_base) + &
                    992.D0*A(L_base+1)+232.D0*A(L_base+2))
      CASE (4)
        AX = HL5 * (-27.D0*A(L_base-1) +378.D0*A(L_base)   +648.D0*A(L_base+1) + &
                    918.D0*A(L_base+2)+243.D0*A(L_base+3))
      END SELECT

      ! Accumulate block contribution
      I5 = I5 + AX + H*(7.D0*S4 + 32.D0*S3 + 12.D0*S2)

    END DO

  END FUNCTION B5

  !> \brief Numerical integration over a single block using a high-order rule.
  !! \ingroup integration
  !! \param[in]  M    Number of steps in the block (>=4)
  !! \param[in]  H    Step size for the block
  !! \param[in]  A    Data array to integrate
  !! \param[in]  IAS  Starting index
  !! \return     I5   Integrated result
  FUNCTION B5_SINGLE(M, H, A, IAS) RESULT(I5)
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: M           ! Number of steps in the block (>=4)
    INTEGER, INTENT(IN)          :: IAS         ! Starting index into A
    DOUBLE PRECISION, INTENT(IN) :: H           ! Step size for the block
    DOUBLE PRECISION, INTENT(IN) :: A(*)        ! Data array (size >= M+3)
    DOUBLE PRECISION             :: I5          ! Integrated result

    INTEGER :: NB, N1, I, idx, L0, L_base
    DOUBLE PRECISION :: HL5, S2, S3, S4, AX

    ! Require at least 4 intervals
    IF (M < 4) THEN
      WRITE(*,*) "Error: Number of intervals must be >= 4."
      STOP
    END IF

    ! Setup
    L0   = IAS
    HL5  = H / 32.D0
    NB   = M / 4
    N1   = M - 4*NB + 1
    L_base = L0 + 4*NB

    ! Initialize sums
    S2 = 0.D0; S3 = 0.D0; S4 = 0.D0
    I5 = 0.D0

    ! Parallel inner loop
    ! Parallel inner loop with static schedule and minimal overhead
    ! Choose parallel or serial based on workload
    IF (NB .GT. 1000) THEN
      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idx) REDUCTION(+:S2,S3,S4) SCHEDULE(STATIC)
      DO I = 1, NB
        idx = L0 + 4*(I-1)
        S2 = S2 + A(idx+2)
        S3 = S3 + A(idx+1) + A(idx+3)
        S4 = S4 + A(idx)   + A(idx+4)
      END DO
      !!$OMP END PARALLEL DO
    ELSE
      DO I = 1, NB
        idx = L0 + 4*(I-1)
        S2 = S2 + A(idx+2)
        S3 = S3 + A(idx+1) + A(idx+3)
        S4 = S4 + A(idx)   + A(idx+4)
      END DO
    END IF

    ! Leftover correction
    SELECT CASE (N1)
    CASE (1)
      AX = 0.D0
    CASE (2)
      AX = HL5 * (-19.D0*A(L_base-3) +106.D0*A(L_base-2) -264.D0*A(L_base-1) + &
                  646.D0*A(L_base)   +251.D0*A(L_base+1))
    CASE (3)
      AX = HL5 * (-8.D0*A(L_base-2)  + 32.D0*A(L_base-1) +192.D0*A(L_base) + &
                  992.D0*A(L_base+1)+232.D0*A(L_base+2))
    CASE (4)
      AX = HL5 * (-27.D0*A(L_base-1) +378.D0*A(L_base)   +648.D0*A(L_base+1) + &
                  918.D0*A(L_base+2)+243.D0*A(L_base+3))
    END SELECT

    ! Final result
    I5 = AX + H*(7.D0*S4 + 32.D0*S3 + 12.D0*S2)

  END FUNCTION B5_SINGLE

  !> \brief Computes Gauss-Laguerre quadrature points and weights.
  !! \ingroup integration
  !! \param[in]  N        Number of quadrature points (1 ≤ N ≤ 600)
  !! \param[out] XPNT     Array(N) of quadrature points (roots)
  !! \param[out] PWEIGHT  Array(N) of quadrature weights
  SUBROUTINE GAULAG(N,XPNT,PWEIGHT)
    ! IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    IMPLICIT NONE
    INTEGER, PARAMETER :: NNR=600, MAXIT=10
    DOUBLE PRECISION, PARAMETER :: EPS=1.D-12
    INTEGER :: I, J, ITS
    DOUBLE PRECISION :: AI, Z, Z1, P1, P2, P3, PP, DN, ALF
    DOUBLE PRECISION :: X(NNR)
    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION, INTENT(OUT) ::  XPNT(NNR),PWEIGHT(NNR)


    IF(N.LT.1)    STOP 'GAULAG.F: BAD ARGUMENT N '
    IF(N.GT.NNR)  STOP 'GAULAG.F: N > NNR TOO LARGE'

    ALF=0.D0
    Z = 0.D0
    DO I=1,N
  ! APPROXIMATE THE ITH ROOT
      IF(I.EQ.1)THEN
        Z=(1.+ALF)*(3.+.92*ALF)/(1.+2.4*N+1.8*ALF)
      ELSEIF(I.EQ.2)THEN
        Z=Z+(15.+6.25*ALF)/(1.+.9*ALF+2.5*N)
      ELSE
        AI=I-2
        Z=Z+((1.+2.55*AI)/(1.9*AI)+1.26*AI*ALF/(1.+3.5*AI))*(Z-X(I-2))/(1.+.3*ALF)
      ENDIF

      DO ITS=1,MAXIT
        P1=1.D0
        P2=0.D0
        DO J=1,N ! RECURRENCE RELATION FOR LAGUERRE POLYNOMIAL IN Z
          P3=P2
          P2=P1
          P1=((2*J-1+ALF-Z)*P2-(J-1+ALF)*P3)/J
        ENDDO
        PP=(N*P1-(N+ALF)*P2)/Z ! DERIVATIVE OF LAGUERRE POLYNOMIAL
        Z1=Z
        Z=Z1-P1/PP ! NEWTON'S METHOD TO REFINE ROOT
        IF(DABS(Z-Z1).LE.EPS) EXIT
      ENDDO

      IF ( ITS > MAXIT ) THEN
        WRITE(*,'("N = ", I3," I = ",I3)') N,I
        STOP 'GAULAG.F: TOO MANY ITERATIONS'
      ENDIF
      X(I)=Z
      XPNT(I)=Z
      DN=REAL(N, KIND=KIND(1.0D0))
      PWEIGHT(I)=-1.D0/(PP*N*P2) ! WEIGHT
    ENDDO
  END SUBROUTINE GAULAG

  !> \brief Generate an exponentially growing radial grid and weights.
  !! \ingroup integration
  !! \param[in]  H     Initial step size
  !! \param[in]  AF    Exponential growth factor (>1)
  !! \param[inout] RANGE  Maximum range (may be updated)
  !! \param[out] R     Allocatable array of grid points
  !! \param[out] W     Allocatable array of weights
  !! \param[out] N     Number of grid points
  SUBROUTINE EXPONENTIALLY_GROWING_GRID(H, AF, RANGE, R, W, N)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: H, AF
    DOUBLE PRECISION, INTENT(INOUT) :: RANGE
    DOUBLE PRECISION, ALLOCATABLE, INTENT(OUT):: R(:), W(:)
    INTEGER, INTENT(OUT) :: N

    INTEGER :: NP, I


    N = INT(1 + DLOG(1.D0 + RANGE*(AF-1.D0)/H)/DLOG(AF))

    IF (.NOT.ALLOCATED(R)) THEN
      ALLOCATE(R(N))
    ELSE
      IF (SIZE(R) .NE. N) THEN
        DEALLOCATE(R)
        ALLOCATE(R(N))
      ENDIF
    ENDIF
    IF (.NOT.ALLOCATED(W)) THEN
      ALLOCATE(W(N))
    ELSE
      IF (SIZE(W) .NE. N) THEN
        DEALLOCATE(W)
        ALLOCATE(W(N))
      ENDIF
    ENDIF

    NP = N+1
    RANGE = H*(AF**NP - 1.D0)/(AF - 1.D0)
    DO I=1, N
      R(I) = H*(AF**I-1.D0)/(AF-1.D0)
      W(I) = R(I)**2*AF**I*DLOG(AF)/(AF-1.D0)
    ENDDO
  END SUBROUTINE EXPONENTIALLY_GROWING_GRID

END MODULE INTEGRATION_MOD