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
MODULE Integration

IMPLICIT NONE

CONTAINS
!-- SUBROUTINE ComputeNDS
    SUBROUTINE ComputeNDS(Mesh,c,iCase,Direction,Axis,NDS)  
    USE MMesh
    IMPLICIT NONE
    TYPE(TMesh) :: Mesh
    INTEGER :: c,iCase
    REAL,DIMENSION(3) :: Direction,Axis
    REAL,DIMENSION(*) :: NDS
    REAL,DIMENSION(3) :: VEL
    INTEGER :: i
    SELECT CASE (iCase)
    CASE (1)
        DO i=1,Mesh%Npanels
            IF (Mesh%cPanel(i).EQ.c) THEN
                VEL(1)=Direction(1)
                VEL(2)=Direction(2)
                VEL(3)=Direction(3)
                NDS(i)=(Mesh%N(1,i)*VEL(1)+Mesh%N(2,i)*VEL(2)+Mesh%N(3,i)*VEL(3))*Mesh%A(i)
            ELSE
                NDS(i)=0.
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                IF (Mesh%cPanel(i).EQ.c) THEN
                    VEL(1)=Direction(1)
                    VEL(2)=Direction(2)
                    VEL(3)=Direction(3)
                    NDS(i+Mesh%Npanels)=(Mesh%N(1,i)*VEL(1)-Mesh%N(2,i)*VEL(2)+Mesh%N(3,i)*VEL(3))*Mesh%A(i)
                 ELSE
                    NDS(i+Mesh%Npanels)=0.
                 END IF          
            END IF
        END DO
    CASE (2)
        DO i=1,Mesh%Npanels
            IF (Mesh%cPanel(i).EQ.c) THEN
                VEL(1)=Direction(2)*(Mesh%XM(3,i)-Axis(3))-Direction(3)*(Mesh%XM(2,i)-Axis(2))
                VEL(2)=Direction(3)*(Mesh%XM(1,i)-Axis(1))-Direction(1)*(Mesh%XM(3,i)-Axis(3))
                VEL(3)=Direction(1)*(Mesh%XM(2,i)-Axis(2))-Direction(2)*(Mesh%XM(1,i)-Axis(1))  
                NDS(i)=(Mesh%N(1,i)*VEL(1)+Mesh%N(2,i)*VEL(2)+Mesh%N(3,i)*VEL(3))*Mesh%A(i)
            ELSE
                NDS(i)=0.
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                IF (Mesh%cPanel(i).EQ.c) THEN
                    VEL(1)=Direction(2)*(Mesh%XM(3,i)-Axis(3))-Direction(3)*(-Mesh%XM(2,i)-Axis(2))
                    VEL(2)=Direction(3)*(Mesh%XM(1,i)-Axis(1))-Direction(1)*(Mesh%XM(3,i)-Axis(3))
                    VEL(3)=Direction(1)*(-Mesh%XM(2,i)-Axis(2))-Direction(2)*(Mesh%XM(1,i)-Axis(1))  
                    NDS(i+Mesh%Npanels)=(Mesh%N(1,i)*VEL(1)-Mesh%N(2,i)*VEL(2)+Mesh%N(3,i)*VEL(3))*Mesh%A(i)
                ELSE
                    NDS(i+Mesh%Npanels)=0.
                 END IF 
            END IF
        END DO
    CASE (3)
        WRITE(*,*) 'Error: force case 3 not implemented yet'
        STOP
    CASE DEFAULT
        WRITE(*,*) 'Error: unknown radiation case'
        STOP
    END SELECT
    END SUBROUTINE
END MODULE Integration