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
    MODULE MBodyConditions
!
    TYPE TBodyConditions
        INTEGER :: Nproblems
        INTEGER :: Npanels
        REAL,DIMENSION(:),ALLOCATABLE :: Omega
        COMPLEX,DIMENSION(:,:),ALLOCATABLE :: NormalVelocity
        INTEGER,DIMENSION(:),ALLOCATABLE :: Switch_Potential
        INTEGER,DIMENSION(:),ALLOCATABLE :: Switch_FreeSurface
        INTEGER,DIMENSION(:),ALLOCATABLE :: Switch_Kochin
        REAL,DIMENSION(:),ALLOCATABLE :: Switch_Type
    END TYPE TBodyConditions
!
    CONTAINS
!
!       Operators for creation, copy, initialisation and destruction
!
        SUBROUTINE CreateTBodyConditions(BodyConditions,Nproblems,Npanels)
        IMPLICIT NONE
        TYPE(TBodyConditions) :: BodyConditions
        INTEGER :: Nproblems,Npanels
        BodyConditions%Nproblems=Nproblems
        BodyConditions%Npanels=Npanels
        ALLOCATE(BodyConditions%NormalVelocity(Npanels,Nproblems))
        ALLOCATE(BodyConditions%Omega(Nproblems),BodyConditions%Switch_Potential(Nproblems),BodyConditions%Switch_FreeSurface(Nproblems),BodyConditions%Switch_Kochin(Nproblems),BodyConditions%Switch_Type(Nproblems))
        END SUBROUTINE CreateTBodyConditions
!       --- 
        SUBROUTINE CopyTBodyConditions(BodyConditionsTarget,BodyConditionsSource)
        IMPLICIT NONE
        INTEGER :: i,j,k
        TYPE(TBodyConditions) :: BodyConditionsTarget,BodyConditionsSource      
        CALL CreateTBodyConditions(BodyConditionsTarget,BodyConditionsSource%Nproblems,BodyConditionsSource%Npanels)
        DO i=1,BodyConditionsTarget%Nproblems
            BodyConditionsTarget%Omega(i)=BodyConditionsSource%Omega(i)
            BodyConditionsTarget%Switch_Potential(i)=BodyConditionsSource%Switch_Potential(i)
            BodyConditionsTarget%Switch_FreeSurface(i)=BodyConditionsSource%Switch_FreeSurface(i)
            BodyConditionsTarget%Switch_Kochin(i)=BodyConditionsSource%Switch_Kochin(i)
            BodyConditionsTarget%Switch_Type(i)=BodyConditionsSource%Switch_Type(i)
            DO k=1,BodyConditionsTarget%Npanels
                BodyConditionsTarget%NormalVelocity(k,i)=BodyConditionsSource%NormalVelocity(k,i)
            END DO
        END DO
        END SUBROUTINE CopyTBodyConditions
!       ---
        SUBROUTINE ReadTBodyConditions(BodyConditions,Npanels,namefile)
        IMPLICIT NONE	    
        TYPE(TBodyConditions) :: BodyConditions
        CHARACTER*(*) :: namefile
        INTEGER :: Nproblems,Npanels
        INTEGER :: i,j,k
        REAL,DIMENSION(:),ALLOCATABLE :: RBC,IBC
        OPEN(10,FILE=namefile)
        READ(10,*) Nproblems
        CALL CreateTBodyConditions(BodyConditions,Nproblems,Npanels)
        ALLOCATE(RBC(Nproblems),IBC(Nproblems))
        READ(10,*) (BodyConditions%Omega(i),i=1,Nproblems)
        READ(10,*) (BodyConditions%Switch_Type(i),i=1,Nproblems)
        READ(10,*) (BodyConditions%Switch_Potential(i),i=1,Nproblems)
        READ(10,*) (BodyConditions%Switch_FreeSurface(i),i=1,Nproblems)
        READ(10,*) (BodyConditions%Switch_Kochin(i),i=1,Nproblems)
        DO k=1,Npanels
            READ(10,*) (RBC(i),IBC(i),i=1,Nproblems)
            DO i=1,Nproblems
                BodyConditions%NormalVelocity(k,i)=CMPLX(RBC(i),IBC(i))
            END DO
        END DO
        CLOSE(10)
        DEALLOCATE(RBC,IBC)
        END SUBROUTINE ReadTBodyConditions
!       --- 
        SUBROUTINE DeleteTBodyConditions(BodyConditions)
        IMPLICIT NONE
        TYPE(TBodyConditions) :: BodyConditions
        DEALLOCATE(BodyConditions%Omega,BodyConditions%Switch_potential,BodyConditions%Switch_FreeSurface,BodyConditions%Switch_Kochin,BodyConditions%Switch_Type)
        DEALLOCATE(BodyConditions%NormalVelocity)
        END SUBROUTINE DeleteTBodyConditions  
!       ---
END MODULE