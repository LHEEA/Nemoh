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
!   Locals
    REAL                    :: PI
    INTEGER			        :: c,d,M,l,i,j,k
    REAL                    :: Discard
    COMPLEX,PARAMETER       :: II=CMPLX(0.,1.)
      
!
!   --- Initialisation -------------------------------------------------------------------------------------------------------------------------------------------------------------------
!     
    PI=4.*ATAN(1.)
    WRITE(*,*) ' '
    WRITE(*,'(A,$)') '  -> Initialisation ' 
!   Read case ID
    CALL ReadTID(ID,'ID.dat')    
    WRITE(*,'(A,$)') '.'
!   Read environement
    CALL ReadTEnvironment(Environment,ID%ID(1:ID%lID)//'/Nemoh.cal') 
!   Read results
    CALL ReadTResults(Results,ID%ID(1:ID%lID)//'/results/Forces.dat',ID%ID(1:ID%lID)//'/results/index.dat',ID%ID(1:ID%lID)//'/results/FKForce.tec')
    CALL SaveTResults(Results,ID%ID(1:ID%lID)//'/results')
    WRITE(*,*) '. Done !'
    WRITE(*,*) ' '
!
!   --- Compute IRFs -------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    CALL Initialize_IRF(IRF,Results,ID%ID(1:ID%lID)//'/Nemoh.cal')  
    IF (IRF%Switch.EQ.1) THEN
        CALL Compute_IRF(IRF,Results)
        CALL Save_IRF(IRF,ID%ID(1:ID%lID)//'/results/IRF.tec')
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
    CALL Initialize_Plot_WaveElevation(Switch_Plot_WaveElevation,ID%ID(1:ID%lID)//'/Nemoh.cal')   
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
