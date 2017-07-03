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
MODULE MIdentification

  ! Definition of TYPE TID
  TYPE TID
    CHARACTER(LEN=80) :: ID
    INTEGER :: lID
  END TYPE TID

CONTAINS

  SUBROUTINE ReadTID(ID)
    ! Read ID data from commad line argument

    IMPLICIT NONE
    TYPE(TID) :: ID

    IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
      CALL GET_COMMAND_ARGUMENT(1, ID%ID)
    ELSE
      ID%ID = "."
    END IF
    ID%lID = LEN(TRIM(ID%ID)) ! for legacy...

  END SUBROUTINE

END MODULE MIdentification

