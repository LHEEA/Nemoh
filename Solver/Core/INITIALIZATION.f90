MODULE INITIALIZATION

IMPLICIT NONE

CONTAINS
!---------------------------------------------------------------------------
    SUBROUTINE INITIALIZE(ID,NF,NSYM,XF,YF,Mesh)
!
    USE MIDENTIFICATION
    USE COM_VAR
    USE PREPARE_MESH
    USE MMesh
!
    IMPLICIT NONE
!   ID
    TYPE(TID) :: ID
!   Geometry
    INTEGER :: NF,NSYM
    REAL :: XF,YF
    TYPE(TMesh) :: Mesh    
!
!   Read input file and geometry
    OPEN(10,file=ID%ID(1:ID%lID)//'/nemoh.cal',form='formatted',status='old')
    READ(10,*)
    READ(10,*) RHO
    READ(10,*) G
    READ(10,*) DEPTH
    READ(10,*) XF,YF
    CLOSE(10)
    OPEN(10,file=ID%ID(1:ID%lID)//'/input.txt',form='formatted',status='old')
    READ(10,*) 
    READ(10,*) Indiq_solver
    READ(10,*) IRES
    READ(10,*) TOL_GMRES
    READ(10,*) MAXIT
    CLOSE(10)   
    XEFF=XF
    YEFF=YF
    NFA=Mesh%Npanels
    NP=Mesh%Npoints
    NSYMY=Mesh%Isym
    IF (NSYMY.NE.1) NSYMY=0
    IF (NSYMY.EQ.1) YEFF=0
    NF=NFA
    NSYM=NSYMY
!   Initialise Nemoh
    CALL ALLOCATE_DATA
    CALL PRE_PROC_MESH(Mesh)
    CALL CREK(251)
!
    END SUBROUTINE INITIALIZE  

END MODULE INITIALIZATION