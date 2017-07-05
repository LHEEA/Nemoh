!--------------------------------------------------------------------------------------
! NEMOH Solver
! See license and contributors list in the main directory.
!--------------------------------------------------------------------------------------
MODULE FORCES
  ! Deduce from the potential the forces on the floating body (for the diffraction problems)
  ! and the added mass and damping (for the radiation problem).

  USE Constants
  USE MMesh,     ONLY: TMesh

  IMPLICIT NONE

  PUBLIC  :: COMPUTE_AND_WRITE_FORCES
  PRIVATE :: INITIALIZE_FORCE_PARAMETERS, APPEND_FORCE_TO_FILE

  PRIVATE
  INTEGER :: Nintegration
  REAL, DIMENSION(:, :), ALLOCATABLE :: NDS

CONTAINS

  !-------------------------------------------

  SUBROUTINE COMPUTE_AND_WRITE_FORCES  &
      ! Main subroutine. Called from NEMOH.f90.
    ( parameters_file,                 &
      Mesh, rho, omega, Potential,     &
      switch_type, output_file         &
      )

    CHARACTER(LEN=*),                               INTENT(IN) :: parameters_file, output_file
    TYPE(TMesh),                                    INTENT(IN) :: Mesh
    REAL,                                           INTENT(IN) :: omega, rho
    INTEGER,                                        INTENT(IN) :: switch_type
    COMPLEX, DIMENSION(Mesh%NPanels*(Mesh%Isym+1)), INTENT(IN) :: Potential

    ! Local variables
    INTEGER :: i, j
    COMPLEX, DIMENSION(:), ALLOCATABLE :: Force

    CALL INITIALIZE_FORCE_PARAMETERS(parameters_file, mesh, output_file)

    ALLOCATE(Force(Nintegration))
    Force(:) = 0.0

    ! Compute force coefficients
    DO i = 1, Nintegration
      DO j = 1, Mesh%NPanels*(Mesh%Isym+1)
        Force(i) = Force(i) - RHO * potential(j)*NDS(i, j)
      END DO
    END DO

    CALL APPEND_FORCE_TO_FILE(output_file, switch_type, omega, force)

    DEALLOCATE(Force)

  END SUBROUTINE COMPUTE_AND_WRITE_FORCES

  !-------------------------------------------
  !-------------------------------------------
  !-------------------------------------------

  SUBROUTINE INITIALIZE_FORCE_PARAMETERS(parameters_file, mesh, output_file)
    ! Load NDS array from the given file.

    CHARACTER(LEN=*), INTENT(IN) :: parameters_file, output_file
    TYPE(TMesh),      INTENT(IN) :: Mesh

    INTEGER :: i, j, u

    IF (.NOT. ALLOCATED(NDS)) THEN
      OPEN(NEWUNIT=u, FILE=parameters_file, STATUS='OLD', ACTION='READ')
      READ(u, *) Nintegration
      ALLOCATE(NDS(Nintegration, Mesh%NPanels*2**Mesh%Isym))
      DO i = 1, Nintegration
        READ(u, *) (NDS(i,j), j=1,Mesh%NPanels*2**Mesh%Isym)
      END DO
      CLOSE(u)

      ! Initialize output file...
      OPEN(NEWUNIT=u, FILE=output_file, ACTION='WRITE')
      DO i = 1, Nintegration
        WRITE(u, '(A)') ''
      END DO
      CLOSE(u)
    ELSE
      ! A file has already been loaded.
      ! We assume only one file should be used for each run of the program.
      ! Thus nothing happens.
    END IF

  END SUBROUTINE INITIALIZE_FORCE_PARAMETERS

  !-------------------------------------------
  !-------------------------------------------
  !-------------------------------------------

  SUBROUTINE APPEND_FORCE_TO_FILE(filename, switch_type, omega, force)
    ! Append the newly computed forces at the end of each line.
    ! By reading each line and then rewriting it.

    ! Done this way for backward compatibility with Nemoh 2.0.
    ! (Might be more efficient to add each new case on a new line.)

    CHARACTER(LEN=*),                 INTENT(IN) :: filename
    INTEGER,                          INTENT(IN) :: switch_type
    REAL,                             INTENT(IN) :: Omega
    COMPLEX, DIMENSION(Nintegration), INTENT(IN) :: force

    ! Local variables
    INTEGER :: i, u
    CHARACTER(LEN=10000), DIMENSION(Nintegration) :: old_data

    ! Read existing data.
    OPEN(NEWUNIT=u, FILE=filename, ACTION='READ')
    DO i = 1, Nintegration ! Loop on all lines
      READ(u, '(A)') old_data(i)
    END DO
    CLOSE(u)

    ! Rewrite file with new data appended at the end of the line.
    OPEN(NEWUNIT=u, FILE=filename, ACTION='WRITE')
    DO i = 1, Nintegration ! Loop on all lines

      IF (switch_type == DIFFRACTION_PROBLEM) THEN
        WRITE(u, *) TRIM(old_data(i)), ' ', ABS(omega*Force(i)), ATAN2(AIMAG(Force(i)),REAL(Force(i))) + PI/2

      ELSE IF (switch_type == RADIATION_PROBLEM) THEN
        WRITE(u, *) TRIM(old_data(i)), ' ', REAL(Force(i)), omega*AIMAG(Force(i)) ! Added mass and damping

      END IF
    END DO
    CLOSE(u)

  END SUBROUTINE APPEND_FORCE_TO_FILE

  !-------------------------------------------

END MODULE FORCES
