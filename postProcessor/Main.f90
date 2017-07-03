!--------------------------------------------------------------------------------------
!
!   NEMOH V1.0 - postProcessor - January 2014
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
!   - A. Babarit
!
!--------------------------------------------------------------------------------------
!
    PROGRAM Main
!
    USE MIdentification
    USE MEnvironment
    USE MResults
    USE MIRF
#ifndef GNUFORT
    USE iflport
#endif
!
    IMPLICIT NONE
!   ID
    TYPE(TID)               :: ID
!   Environment
    TYPE(TEnvironment) :: Environment
!   Hydrodynamic coefficients cases
    TYPE(TResults) :: Results
!   IRFs
    TYPE(TIRF) :: IRF
!   RAOs
    COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: RAOs
!   Plot Wave elevation
    REAL :: Switch_Plot_WaveElevation

!
!   --- Initialisation -------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    WRITE(*,*) ' '
    WRITE(*,'(A,$)') '  -> Initialisation '
!   Read case ID
    CALL ReadTID(ID)
    WRITE(*,'(A,$)') '.'
!   Read environement
    CALL ReadTEnvironment(Environment,TRIM(ID%ID)//'/Nemoh.cal')
!   Read results
    CALL ReadTResults(Results,TRIM(ID%ID)//'/results/Forces.dat',TRIM(ID%ID)//'/results/index.dat',TRIM(ID%ID)//'/results/FKForce.tec')
    CALL SaveTResults(Results,TRIM(ID%ID)//'/results')
    WRITE(*,*) '. Done !'
    WRITE(*,*) ' '
!
!   --- Compute IRFs -------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    CALL Initialize_IRF(IRF,Results,TRIM(ID%ID)//'/Nemoh.cal')
    IF (IRF%Switch.EQ.1) THEN
        CALL Compute_IRF(IRF,Results)
        CALL Save_IRF(IRF,TRIM(ID%ID)//'/results/IRF.tec')
    END IF
!
!   --- Compute RAOs -------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    ALLOCATE(RAOs(Results%Nintegration,Results%Nw,Results%Nbeta))
    CALL Compute_RAOs(RAOs,Results)
!
!   --- Save results -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    WRITE(*,*) ' -> Save results '
    WRITE(*,*) ' '
    CALL Initialize_Plot_WaveElevation(Switch_Plot_WaveElevation,TRIM(ID%ID)//'/Nemoh.cal')
    IF (Switch_Plot_WaveElevation.GT.1) THEN
        CALL Plot_WaveElevation(ID,Environment,1,1,RAOs,Results)
    END IF

!
!   --- Finalize -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    CALL DeleteTResults(Results)
    DEALLOCATE(RAOs)
!
    END PROGRAM Main
!
