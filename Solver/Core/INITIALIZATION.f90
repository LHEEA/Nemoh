MODULE INITIALIZATION

IMPLICIT NONE

CONTAINS
!---------------------------------------------------------------------------
    SUBROUTINE INITIALIZE(ID,NF,NSYM)
!
    USE MIDENTIFICATION
    USE COM_VAR
    USE PREPARE_MESH
!
    IMPLICIT NONE
!   ID
    TYPE(TID) :: ID
!   Geometry
    INTEGER :: NF,NSYM
!
!   Read input file and geometry
    OPEN(10,file=ID%ID(1:ID%lID)//'/aquaplus.cal',form='formatted',status='old')
    READ(10,*)
    READ(10,*) RHO
    READ(10,*) G
    READ(10,*) DEPTH
    READ(10,*) XEFF,YEFF
    CLOSE(10)
    OPEN(10,file=ID%ID(1:ID%lID)//'/input.txt',form='formatted',status='old')
    READ(10,*) 
    READ(10,*) Indiq_solver
    READ(10,*) IRES
    READ(10,*) TOL_GMRES
    READ(10,*) MAXIT
    READ(10,*) Sav_potential
    CLOSE(10)
    MESHFILE=ID%ID(1:ID%lID)//'/Mesh/L12.dat'
    LFILE=ID%lID+13
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Mesh/L10.dat')
    READ(10,*)
    READ(10,*) NSYMY,NP,NFA
    CLOSE(10)
    IF (NSYMY.NE.1) NSYMY=0
    IF (NSYMY.EQ.1) YEFF=0
    NF=NFA
    NSYM=NSYMY
!   Initialise Aquaplus
    CALL ALLOCATE_DATA
    CALL PRE_PROC_MESH
    CALL CREK(251)
!
    END SUBROUTINE INITIALIZE  

END MODULE INITIALIZATION