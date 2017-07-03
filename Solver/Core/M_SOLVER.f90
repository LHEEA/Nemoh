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
MODULE M_SOLVER

  IMPLICIT NONE

  PUBLIC :: GAUSSZ

CONTAINS
  !---------------------------------------------------------------------------!

  SUBROUTINE GAUSSZ(A, N, Ainv)

    ! Input
    COMPLEX, DIMENSION(N, N), INTENT(IN)  :: A
    INTEGER,                  INTENT(IN)  :: N

    ! Output
    COMPLEX, DIMENSION(N, N), INTENT(OUT) :: Ainv

    ! Local variables
    INTEGER :: I, J, K, L, IL
    COMPLEX :: C, P
    COMPLEX, DIMENSION(N, 2*N) :: WS ! Working matrix
    REAL, PARAMETER :: EPS=1E-20

    ! Initialize working matrix
    WS(:, 1:N)     = A
    WS(:, N+1:2*N) = CMPLX(0., 0.)
    DO I = 1, N
      WS(I, N+I) = CMPLX(1., 0.)
    END DO

    ! Gauss pivot inversion
    DO J = 1, N-1
      K = J
      DO I = J+1, N
        IF ((ABS(WS(K, J))-ABS(WS(I, J))) <= 0.0) K = I
      END DO
      IF (K <= J) THEN 
        DO L = J, 2*N
          C = WS(J, L)
          WS(J, L) = WS(K, L)
          WS(K, L) = C
        END DO
      ELSE 
        IF (ABS(WS(J, J)) <= EPS) THEN
          WRITE(*, '(A,E16.6)') 'PIVOT INFERIEUR A ', EPS
          STOP
        END IF
      END IF
      DO K = J+1, 2*N
        P = WS(J, K)/WS(J, J)
        DO I = J+1, N
          WS(I, K) = WS(I, K) - WS(I, J)*P
        END DO
      END DO
    END DO
    IF (ABS(WS(N, N)) < EPS) THEN 
      WRITE(*, '(A,E16.6)') 'PIVOT INFERIEUR A ', EPS
      STOP
    END IF
    DO IL = N+1, 2*N
      DO J = N, 1, -1
        WS(J, IL) = WS(J, IL)/WS(J, J)
        DO I = 1, J-1
          WS(I, IL) = WS(I, IL) - WS(I, J)*WS(J, IL)
        END DO
      END DO
    END DO

    ! Extract output
    Ainv(:, :) = WS(:, N+1:2*N)

    RETURN

  END SUBROUTINE

END MODULE 
