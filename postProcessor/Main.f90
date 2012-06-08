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
    USE iflport
!    
    IMPLICIT NONE
!   ID
    TYPE(TID)               :: ID
!   Environment
    TYPE(TEnvironment) :: Environment
!   Hydrodynamic coefficients cases
    TYPE(TResults) :: RadiationResults
    TYPE(TResults) :: DiffractionResults
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
!   Read Environement
    CALL ReadTEnvironment(Environment,ID%ID(1:ID%lID)//'/aquaplus.cal') 
!   Read radiation results
    CALL ReadTResults(RadiationResults,ID%ID(1:ID%lID)//'/results/Radiation.dat')
!   Read radiation results
    CALL ReadTResults(DiffractionResults,ID%ID(1:ID%lID)//'/results/Diffraction.dat')
!
!   --- Compute RAOs -------------------------------------------------------------------------------------------------------------------------------------------------------------------
!    
    ALLOCATE(RAOs(RadiationResults%Ncase,DiffractionResults%Nperiod,DiffractionResults%Ncase))
    CALL Compute_RAOs(RAOs,RadiationResults,DiffractionResults)
!
!   --- Save results -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    WRITE(*,*) ' -> Save results ' 
    WRITE(*,*) ' '
    CALL Plot_WaveElevation(ID,Environment,1,1,RAOs,RadiationResults,DiffractionResults)
     
!
!   --- Finalize -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    CALL DeleteTResults(RadiationResults)
    CALL DeleteTResults(DiffractionResults)
    DEALLOCATE(RAOs)
!
    END PROGRAM Main
!      