! ----------------------------------------------
!
!   postProcessing
!
! ----------------------------------------------
!
    PROGRAM Main
!   
    USE MIdentification
    USE MEnvironment
    USE MResults
    USE MIRF
    USE iflport
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
    CALL ReadTEnvironment(Environment,ID%ID(1:ID%lID)//'/aquaplus.cal') 
!   Read results
    CALL ReadTResults(Results,ID%ID(1:ID%lID)//'/results/forces.dat',ID%ID(1:ID%lID)//'/results/index.dat',ID%ID(1:ID%lID)//'/results/FKForce.tec')
    CALL SaveTResults(Results,ID%ID(1:ID%lID)//'/results')
!
!   --- Compute IRFs -------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    CALL Initialize_IRF(IRF,Results,ID%ID(1:ID%lID)//'/aquaplus.cal')  
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
    CALL Plot_WaveElevation(ID,Environment,1,1,RAOs,Results)
     
!
!   --- Finalize -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    CALL DeleteTResults(Results)
    DEALLOCATE(RAOs)
!
    END PROGRAM Main
!      