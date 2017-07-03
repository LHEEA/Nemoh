!--------------------------------------------------------------------------------------
!
!   NEMOH V1.0 - BVP solver - January 2014
!
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
!   - A. Babarit
!
!--------------------------------------------------------------------------------------

PROGRAM Main

  USE Constants
  USE MMesh,                ONLY: TMesh,           ReadTMesh
  USE MEnvironment,         ONLY: TEnvironment,    ReadTEnvironment
  USE MBodyConditions,      ONLY: TBodyConditions, ReadTBodyConditions

  ! Preprocessing and initialization
  USE INITIALIZE_GREEN_2,   ONLY: INITIALIZE_GREEN
  USE Elementary_functions, ONLY: X0

  ! Resolution
  USE SOLVE_BEM_DIRECT,     ONLY: SOLVE_POTENTIAL_DIRECT

  ! Post processing and output
  USE OUTPUT,               ONLY: WRITE_DATA_ON_MESH
  USE FORCES,               ONLY: COMPUTE_AND_WRITE_FORCES
  USE KOCHIN,               ONLY: COMPUTE_AND_WRITE_KOCHIN
  USE FREESURFACE,          ONLY: COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION

  IMPLICIT NONE

  CHARACTER(LEN=1000)   :: wd             ! Working directory path (max length: 1000 characters, increase if necessary)
  TYPE(TMesh)           :: Mesh           ! Mesh of the floating body
  TYPE(TBodyConditions) :: BodyConditions ! Physical conditions on the floating body
  TYPE(TEnvironment)    :: Env            ! Physical conditions of the environment

  INTEGER                            :: i_problem          ! Index of the current problem
  REAL                               :: omega, wavenumber  ! Wave frequency and wavenumber
  COMPLEX, DIMENSION(:), ALLOCATABLE :: ZIGB, ZIGS         ! Computed source distribution
  COMPLEX, DIMENSION(:), ALLOCATABLE :: Potential          ! Computed potential

  ! Initialization ---------------------------------------------------------------------

  WRITE(*,*) ' '
  WRITE(*,'(A,$)') '  -> Initialisation '

  ! Get working directory from command line argument
  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, wd)
  ELSE
    wd = "."
  END IF

  CALL ReadTMesh(Mesh, TRIM(wd)//'/mesh/')

  ALLOCATE(ZIGB(Mesh%NPanels), ZIGS(Mesh%NPanels))
  ALLOCATE(Potential(Mesh%NPanels*2**Mesh%Isym))

  CALL ReadTBodyConditions            &
  ( BodyConditions,                   &
    Mesh%Npanels*2**Mesh%Isym,        &
    TRIM(wd)//'/Normalvelocities.dat' &
    )

  CALL ReadTEnvironment(Env, file=TRIM(wd)//'/Nemoh.cal')

  call INITIALIZE_GREEN()

  WRITE(*, *) '. Done !'
  WRITE(*, *) ' '

  ! Solve BVPs and calculate forces ----------------------------------------------------

  WRITE(*, *) ' -> Solve BVPs and calculate forces '
  WRITE(*, *) ' '

  DO i_problem = 1, BodyConditions%Nproblems

    WRITE(*, '(A,I5,A,I5,A,$)') ' Problem ', i_problem, ' / ', BodyConditions%Nproblems, ' .'

    omega = BodyConditions%omega(i_problem) ! Wave frequency

    ! Compute wave number
    IF ((Env%depth == INFINITE_DEPTH) .OR. (omega**2*Env%depth/Env%g >= 20)) THEN
      wavenumber = omega**2/Env%g
    ELSE
      wavenumber = X0(omega**2*Env%depth/Env%g)/Env%depth
      ! X0(y) returns the solution of y = x * tanh(x)
    END IF

    !===============
    ! BEM Resolution
    !===============
    
    CALL SOLVE_POTENTIAL_DIRECT                                              &
    !==========================
    ( Mesh, Env, omega, wavenumber,                                          &
      BodyConditions%NormalVelocity(1:Mesh%Npanels*2**Mesh%Isym, i_problem), &
      ZIGB, ZIGS,                                                            &
      Potential(:))


    !===========================
    ! Post processing and output
    !===========================

    CALL COMPUTE_AND_WRITE_FORCES            &
    !============================
    ( TRIM(wd)//'/mesh/Integration.dat',     &
      Mesh, Env%rho, omega, Potential,       &
      Bodyconditions%Switch_type(i_problem), &
      TRIM(wd)//'/results/Forces.dat'        &
      )

    IF (BodyConditions%Switch_Potential(i_problem) == 1) THEN
      ! Write pressure field on the floating body in file
      CALL WRITE_DATA_ON_MESH                                      &
      !=======================
      ( Mesh,                                                      &
        Env%rho*II*omega*Potential(:),                             &
        TRIM(wd)//'/results/pressure.'//string(i_problem)//'.dat'  &
        )
    END IF

    IF (BodyConditions%Switch_Kochin(i_problem) == 1) THEN
      CALL COMPUTE_AND_WRITE_KOCHIN                             &
      !============================
      ( TRIM(wd)//'/mesh/Kochin.dat',                           &
        Mesh, Env, wavenumber, ZIGB, ZIGS,                      &
        TRIM(wd)//'/results/Kochin.'//string(i_problem)//'.dat' &
        )
    END IF

    IF (BodyConditions%Switch_FreeSurface(i_problem) == 1) THEN
      CALL COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION                  &
      !============================================
      ( TRIM(wd)//'/mesh/Freesurface.dat',                           &
        Mesh, Env, omega, wavenumber, ZIGB, ZIGS,                    &
        TRIM(wd)//'/results/freesurface.'//string(i_problem)//'.dat' &
        )
    END IF

    WRITE(*,*) '. Done !'
  END DO
  WRITE(*,*) ' '

  ! Finalize ---------------------------------------------------------------------------

  DEALLOCATE(ZIGB, ZIGS, Potential)

CONTAINS

  FUNCTION string (i) result (s)
    ! For example 5 -> "00005"
    INTEGER :: i
    CHARACTER(LEN=5) :: s
    WRITE(s, '(I0.5)') i
  END FUNCTION

END PROGRAM Main

