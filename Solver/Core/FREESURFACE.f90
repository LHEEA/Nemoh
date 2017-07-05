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
!   - P. Guével
!   - J.C. Daubisse
!   - J. Singh  
!
!--------------------------------------------------------------------------------------
MODULE FREESURFACE
  ! Computation of the free surface elevation from the computed source distribution.

  USE Constants
  USE MMesh,                     ONLY: TMesh, CreateTMesh
  USE MEnvironment,              ONLY: TEnvironment

  ! Green functions
  USE GREEN_1,                   ONLY: VAV
  USE GREEN_2,                   ONLY: VNSFD, VNSINFD

  USE OUTPUT

  IMPLICIT NONE

  PUBLIC  :: COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION
  PRIVATE :: READ_FREESURFACE_PARAMETERS, COMPUTE_POTENTIAL_AT_POINT, WRITE_FS

  PRIVATE
  TYPE(TMesh) :: MeshFS ! Mesh of the free surface

CONTAINS

  !-------------------------------------------

  SUBROUTINE COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION  &
    ! Main subroutine of the module. Called in NEMOH.f90.
    ( parameters_file,                                 &
      Mesh, Env, omega, wavenumber, ZIGB, ZIGS,        &
      output_file                                      &
      )

    CHARACTER(LEN=*),                 INTENT(IN) :: parameters_file, output_file
    TYPE(TMesh),                      INTENT(IN) :: Mesh
    TYPE(TEnvironment),               INTENT(IN) :: Env
    REAL,                             INTENT(IN) :: omega, wavenumber
    COMPLEX, DIMENSION(Mesh%NPanels), INTENT(IN) :: ZIGB, ZIGS ! Sources

    ! Local variables
    INTEGER :: j
    COMPLEX :: PHI
    COMPLEX, DIMENSION(:), ALLOCATABLE :: ETA

    CALL READ_FREESURFACE_PARAMETERS(parameters_file)

    ALLOCATE(ETA(MeshFS%NPoints))

    DO j = 1, MeshFS%Npoints

      CALL COMPUTE_POTENTIAL_AT_POINT             &
      !==============================
      ( Mesh, Env, omega, wavenumber, ZIGB, ZIGS, &
        MeshFS%X(1, j), MeshFS%X(2, j), 0.0,      &
        PHI                                       &
        )

      ! Get elevation ETA from potential PHI
      ETA(j) = II*omega/Env%G*PHI
    END DO

    CALL WRITE_FS(output_file, ETA)

    DEALLOCATE(ETA)

  END SUBROUTINE COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION

  !-------------------------------------------
  !-------------------------------------------
  !-------------------------------------------

  SUBROUTINE READ_FREESURFACE_PARAMETERS(filename)
    ! Load mesh of the free surface

    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: j, u

    IF (.NOT. ALLOCATED(MeshFS%X)) THEN
      OPEN(NEWUNIT=u, FILE=filename, STATUS='OLD', ACTION='READ')
      READ(u, *) MeshFS%Npoints, MeshFS%Npanels
      IF (MeshFS%Npoints > 0) THEN
        CALL CreateTMesh(MeshFS, MeshFS%Npoints, MeshFS%Npanels, 1)
        DO j = 1, MeshFS%Npoints
          READ(u, *) MeshFS%X(1,j), MeshFS%X(2,j)
        END DO
        DO j = 1, MeshFS%Npanels
          READ(u, *) MeshFS%P(1,j), MeshFS%P(2,j), MeshFS%P(3,j), MeshFS%P(4,j)
        END DO
      END IF
      CLOSE(u)
    ELSE
      ! A file has already been loaded.
      ! We assume only one file should be used for each run of the program.
      ! Thus nothing happens.
    END IF

  END SUBROUTINE READ_FREESURFACE_PARAMETERS

  !-------------------------------------------
  !-------------------------------------------
  !-------------------------------------------

  SUBROUTINE COMPUTE_POTENTIAL_AT_POINT         &
    ( Mesh, Env, omega, wavenumber, ZIGB, ZIGS, &
      XC, YC, ZC,                               &
      PHI                                       &
      )
    ! Compute the potential PHI at the given coordinate (XC, YC, ZC) using the
    ! source distributions ZIGB abd ZIGS.

    ! Inputs
    TYPE(TMesh),                      INTENT(IN) :: Mesh ! Mesh of the floating body
    TYPE(TEnvironment),               INTENT(IN) :: Env
    REAL,                             INTENT(IN) :: omega, wavenumber
    COMPLEX, DIMENSION(Mesh%NPanels), INTENT(IN) :: ZIGB, ZIGS
    REAL,                             INTENT(IN) :: XC, YC, ZC

    ! Output
    COMPLEX,                          INTENT(OUT) :: PHI

    ! Local variables
    INTEGER               :: J
    REAL                  :: FSP, FSM
    REAL, DIMENSION(3)    :: VSXP, VSXM
    COMPLEX               :: SP, SM
    COMPLEX, DIMENSION(3) :: VSP, VSM

    PHI = 0.0

    DO J = 1, Mesh%NPanels

      ! Compute the Greem function coefficients between the point and the face of index J.
      ! Compute also its gradient although it is not used.
      ! Values could be stored for more efficiency as in SOLVE_BEM.

      CALL VAV(0, (/XC, YC, ZC/), J, Mesh, Env%depth, FSP, FSM, VSXP, VSXM)
      IF ((Env%depth == INFINITE_DEPTH) .OR. (wavenumber*Env%depth >= 20)) THEN
        CALL VNSINFD(wavenumber, (/XC, YC, ZC/), J, Mesh, SP, SM, VSP, VSM)
      ELSE
        CALL VNSFD(wavenumber, (/XC, YC, ZC/), J, Mesh, Env%depth, SP, SM, VSP, VSM)
      ENDIF

      ! Compute potential from sources and Green function.
      IF (Mesh%Isym == NO_Y_SYMMETRY) THEN
        PHI = PHI + (SP+FSP)*ZIGB(J)
      ELSE IF (Mesh%Isym == Y_SYMMETRY) THEN
        PHI = PHI + (FSP+SP+FSM+SM)/2*ZIGB(J) + (FSP+SP-FSM-SM)/2*ZIGS(J)
      ENDIF
    END DO

    RETURN
  END SUBROUTINE COMPUTE_POTENTIAL_AT_POINT

  !-------------------------------------------
  !-------------------------------------------
  !-------------------------------------------

  SUBROUTINE WRITE_FS(filename, ETA)
    ! Write the field ETA into a text file.

    CHARACTER(LEN=*),      INTENT(IN) :: filename
    COMPLEX, DIMENSION(*), INTENT(IN) :: ETA

    ! Local variables
    INTEGER          :: i, u

    OPEN(NEWUNIT=u, FILE=filename, ACTION='WRITE')

    IF (output_format == TECPLOT_OUTPUT) THEN
      WRITE(u, *) 'VARIABLES="X" "Y" "abs(eta) (m)" "angle(phi) (rad)" "PRE1" "PRE2"'
      WRITE(u, '(A,I7,A,I7,A)') 'ZONE N=', MeshFS%Npoints, ' , E = ', MeshFS%Npanels, ' , F=FEPOINT,ET=QUADRILATERAL'
    END IF

    DO i = 1, MeshFS%Npoints
      WRITE(u, '(6(X, E14.7))') MeshFS%X(1, i), MeshFS%X(2, i), ABS(eta(i)), ATAN2(IMAG(eta(i)), REAL(eta(I))), REAL(eta(i)), IMAG(eta(i))
    END DO

    IF (output_format == TECPLOT_OUTPUT) THEN
      DO i=1,MeshFS%Npanels
        WRITE(u, *) MeshFS%P(1, i), MeshFS%P(2, i), MeshFS%P(3, i), MeshFS%P(4, i)
      END DO
    END IF

    CLOSE(u)

  END SUBROUTINE WRITE_FS

  !-------------------------------------------

END MODULE FREESURFACE
