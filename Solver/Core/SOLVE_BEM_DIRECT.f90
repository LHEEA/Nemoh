!--------------------------------------------------------------------------------------
!
!   Copyright 2014 Ecole Centrale de Nantes, 1 rue de la Noë, 44300 Nantes, France
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
!
!   Contributors list:
!   - G. Delhommeau
!   - J. Singh
!   - P. Guével
!   - J.C. Daubisse
!
!--------------------------------------------------------------------------------------
MODULE SOLVE_BEM_DIRECT
  ! Resolution of the boundary elements problem

  USE Constants
  USE MMesh,              ONLY: TMesh
  USE MEnvironment,       ONLY: TEnvironment

  ! Green functions
  USE INITIALIZE_GREEN_2, ONLY: LISC
  USE GREEN_1,            ONLY: VAV
  USE GREEN_2,            ONLY: VNSINFD, VNSFD

  ! Solver for linear problem
  USE M_SOLVER,           ONLY: GAUSSZ

  IMPLICIT NONE

  PUBLIC :: SOLVE_POTENTIAL_DIRECT

  PRIVATE

  ! Those variables will be conserved between calls of the subroutine.
  REAL :: Omega_previous = -1.0
  COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: S
  COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: Vinv

CONTAINS
  !--------------------------------------------------------------------------!
  SUBROUTINE SOLVE_POTENTIAL_DIRECT  &
    ( Mesh, Env, omega, wavenumber,  &
      NVel, ZIGB, ZIGS, Potential    &
      )
    ! In this subroutine the linear system Ax=b is constructed and solved

    ! Input/output
    TYPE(TMesh),                      INTENT(IN)  :: Mesh
    TYPE(TEnvironment),               INTENT(IN)  :: Env
    COMPLEX, DIMENSION(*),            INTENT(IN)  :: NVel
    REAL,                             INTENT(IN)  :: omega, wavenumber
    COMPLEX, DIMENSION(Mesh%Npanels), INTENT(OUT) :: ZIGB, ZIGS ! Source distribution
    COMPLEX, DIMENSION(Mesh%Npanels*2**Mesh%Isym), INTENT(OUT) :: Potential

    INTEGER :: I, J

    ! Return of GREEN_1 module
    REAL :: FSP, FSM
    REAL, DIMENSION(3) :: VSXP, VSXM

    ! Return of GREEN_2 module
    COMPLEX :: SP, SM
    COMPLEX, DIMENSION(3) :: VSP, VSM

    COMPLEX, DIMENSION(Mesh%NPanels, Mesh%NPanels, 2**Mesh%ISym) :: V
    COMPLEX, DIMENSION(Mesh%NPanels, 2**Mesh%ISym)               :: ZOL

    ! Allocation of the matrices
    IF (.NOT. ALLOCATED(Vinv)) THEN
      ALLOCATE(S   (Mesh%NPanels, Mesh%NPanels, 2**Mesh%ISym))
      ALLOCATE(Vinv(Mesh%NPanels, Mesh%NPanels, 2**Mesh%ISym))
    END IF

    ! IF (ABS(Omega) <= 1.E-4) THEN
    !   WRITE(*,*)'ABS(Omega)  = ',ABS(Omega),' < 1.E-4'
    !   STOP
    ! END IF

    ! TODO: Clean up those warning messages...
    ! IF (AMH-AKH <= 1e-3) THEN
    !   PRINT*, 'Depth is almost infinite: the infinite depth solver might be more efficient'
    ! END IF
    ! IF (omega**2*depth/g <= 0.1) THEN
    !   PRINT*, 'Depth is too low for the given wavelength'
    ! END IF

    !=====================================
    ! Construction of the influence matrix
    !=====================================

    IF (omega /= omega_previous) THEN
      ! Do not recompute if the same frequency is studied twice
      omega_previous = omega

      ! Initialization of the Green function computations
      IF (.NOT. Env%depth == INFINITE_DEPTH) THEN
        CALL LISC(omega**2*Env%depth/Env%g, wavenumber*Env%depth)
      END IF

      DO I = 1, Mesh%NPanels
        DO J = 1, Mesh%NPanels

          ! First part of the Green function
          CALL VAV                                &
          ( I, Mesh%XM(:, I), J, Mesh, Env%depth, &
            FSP, FSM, VSXP, VSXM                  &
            )
          ! These output are actually independent of omega and could be computed only once.

          ! Second part of the Green function
          IF (Env%depth == INFINITE_DEPTH) THEN
            CALL VNSINFD                          &
            ( wavenumber, Mesh%XM(:, I), J, Mesh, &
              SP, SM, VSP, VSM                    &
              )
          ELSE
            CALL VNSFD                                       &
            ( wavenumber, Mesh%XM(:, I), J, Mesh, Env%depth, &
              SP, SM, VSP, VSM                               &
              )
          END IF

          ! Store into influence matrix
          S(I, J, 1) = FSP + SP                              ! Green function
          V(I, J, 1) = DOT_PRODUCT(Mesh%N(:, I), VSXP + VSP) ! Gradient of the Green function

          IF (Mesh%ISym == Y_SYMMETRY) THEN
            S(I, J, 2) = FSM + SM
            V(I, J, 2) = DOT_PRODUCT(Mesh%N(:, I), VSXM + VSM)
          ENDIF

        END DO
      END DO

      ! Invert matrix V
      CALL GAUSSZ(V(:, :, 1), Mesh%NPanels, Vinv(:, :, 1))
      IF (Mesh%ISym == Y_SYMMETRY) THEN
        CALL GAUSSZ(V(:, :, 2), Mesh%NPanels, Vinv(:, :, 2))
      END IF
    END IF

    !=========================
    ! Solve the linear problem
    !=========================
     
    IF (Mesh%ISym == NO_Y_SYMMETRY) THEN
      ZIGB(:) = MATMUL(Vinv(:, :, 1), NVEL(1:Mesh%NPanels))

    ELSE IF (Mesh%ISym == Y_SYMMETRY) THEN
      ZOL(:, 1) = MATMUL(                                              &
        Vinv(:, :, 1),                                                 &
        (NVEL(1:Mesh%NPanels) + NVEL(Mesh%NPanels+1:2*Mesh%NPanels))/2 &
        )

      ZOL(:, 2) = MATMUL(                                              &
        Vinv(:, :, 2),                                                 &
        (NVEL(1:Mesh%NPanels) - NVEL(Mesh%NPanels+1:2*Mesh%NPanels))/2 &
        )

      ZIGB(:) = ZOL(:, 1) + ZOL(:, 2)
      ZIGS(:) = ZOL(:, 1) - ZOL(:, 2)
    END IF

    !=============================================================
    ! Computation of potential phi = S*source on the floating body
    !=============================================================

    IF (Mesh%ISym == NO_Y_SYMMETRY) THEN
      Potential(:) = MATMUL(S(:, :, 1), ZIGB(:))

    ELSE IF (Mesh%ISym == Y_SYMMETRY) THEN
      Potential(1:Mesh%NPanels) =                    &
        MATMUL(S(:, :, 1) + S(:, :, 2), ZIGB(:))/2   &
        + MATMUL(S(:, :, 1) - S(:, :, 2), ZIGS(:))/2
      Potential(Mesh%NPanels+1:) =                   &
        MATMUL(S(:, :, 1) - S(:, :, 2), ZIGB(:))/2   &
        + MATMUL(S(:, :, 1) + S(:, :, 2), ZIGS(:))/2
    END IF

    RETURN
  END SUBROUTINE
END MODULE
