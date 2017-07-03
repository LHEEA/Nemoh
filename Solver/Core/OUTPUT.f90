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
!   - A. Babarit 
!
!--------------------------------------------------------------------------------------
MODULE OUTPUT

  USE Constants
  USE MMesh

  IMPLICIT NONE

  PUBLIC :: WRITE_DATA_ON_MESH

  PUBLIC 
  INTEGER, PARAMETER :: RAW_OUTPUT=0, TECPLOT_OUTPUT=1
  INTEGER :: output_format = TECPLOT_OUTPUT

CONTAINS

  SUBROUTINE WRITE_DATA_ON_MESH(Mesh, cdata, filename)
    ! Save in a file a field of complex values on a mesh.

    TYPE(TMesh),                                   INTENT(IN) :: Mesh
    COMPLEX, DIMENSION(Mesh%NPanels*2**Mesh%ISym), INTENT(IN) :: cdata
    CHARACTER(LEN=*),                              INTENT(IN) :: filename ! Output file
    
    ! Local variables
    INTEGER            :: u, i, j
    REAL, DIMENSION(3) :: x

    OPEN(NEWUNIT=u, FILE=filename, ACTION='WRITE')

    IF (output_format == TECPLOT_OUTPUT) THEN
      WRITE(u,*) 'VARIABLES="X" "Y" "Z" "abs(p) (Pa)" "angle(p) (rad)"'
      WRITE(u,'(A,I7,A,I7,A)') 'ZONE N=', Mesh%Npoints*2**Mesh%ISym,',E = ', Mesh%Npanels*2**Mesh%ISym,', F=FEPOINT, ET=QUADRILATERAL'
    END IF 

    DO i = 1, Mesh%NPanels
      DO j = 1, 4
        WRITE(u, *) Mesh%X(:, Mesh%P(j, i)), ABS(cdata(i)), ATAN2(AIMAG(cdata(i)),REAL(cdata(i)))
      END DO
    END DO

    IF (Mesh%ISym == Y_SYMMETRY) THEN
      DO i = Mesh%NPanels+1, 2*Mesh%NPanels
        DO j = 1, 4
          x(:) = Mesh%X(:, Mesh%P(j, i-Mesh%NPanels))
          x(2) = -x(2)
          WRITE(u, *) x(:), ABS(cdata(i)), ATAN2(AIMAG(cdata(i)),REAL(cdata(i)))
        END DO
      END DO
    END IF

    IF (output_format == TECPLOT_OUTPUT) THEN
      DO i=1,Mesh%Npanels
        WRITE(u, *) 1+(i-1)*4, 2+(i-1)*4, 3+(i-1)*4, 4+(i-1)*4
      END DO

      IF (Mesh%ISym == Y_SYMMETRY) THEN
        DO i=1,Mesh%Npanels
          WRITE(u, *) 1+(i-1)*4+4*Mesh%NPanels, 2+(i-1)*4+4*Mesh%NPanels, 3+(i-1)*4+4*Mesh%NPanels, 4+(i-1)*4+4*Mesh%NPanels
        END DO
      END IF
    END IF

    CLOSE(u)

  END SUBROUTINE WRITE_DATA_ON_MESH

END MODULE OUTPUT
