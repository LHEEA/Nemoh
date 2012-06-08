MODULE OUTPUT

TYPE TMeshFS
    INTEGER :: Npoints,Npanels
    REAL,DIMENSION(:),ALLOCATABLE :: XC,YC,ZC
    INTEGER,DIMENSION(:,:),ALLOCATABLE :: P
END TYPE TMeshFS

CONTAINS
!
    SUBROUTINE CreateTMeshFS(MeshFS,Npoints,Npanels)
    IMPLICIT NONE
    TYPE(TMeshFS):: MeshFS
    INTEGER :: Npoints,Npanels
    MeshFS%NPoints=Npoints
    MeshFS%Npanels=Npanels
    ALLOCATE(MeshFS%XC(Npoints),MeshFS%YC(Npoints),MeshFS%ZC(Npoints))
    ALLOCATE(MeshFS%P(4,Npanels))
    END SUBROUTINE
!
    SUBROUTINE DeleteTMeshFS(MeshFS)
    IMPLICIT NONE
    TYPE(TMeshFS):: MeshFS
    DEALLOCATE(MeshFS%XC,MeshFS%YC,MeshFS%ZC,MeshFS%P)
    END SUBROUTINE
!
    SUBROUTINE WRITE_FS(ID,MeshFS,PHI)
    USE MIDENTIFICATION
    USE COM_VAR
    TYPE(TID) :: ID
    TYPE(TMeshFS) :: MeshFS
    INTEGER:: i,j
    COMPLEX,DIMENSION(*) :: PHI
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/FS.tec')
    WRITE(11,'(A)') 'VARIABLES="X" "Y" "Z" "Module" "Phase" "Real" "Imag" "Re(NVEL)" "Im(NVEL)"'
    WRITE(11,'(A,I7,A,I7,A)') 'ZONE t="Free surface",N=',MeshFS%Npoints,', E=',MeshFS%Npanels,' , F=FEPOINT,ET=QUADRILATERAL'  
    DO i=1,MeshFS%Npoints
        WRITE(11,'(9(2X,E14.7))') MeshFS%XC(i),MeshFS%YC(i),MeshFS%ZC(i),ABS(PHI(i)),ATAN2(IMAG(PHI(i)),REAL(PHI(i))),REAL(PHI(i)),IMAG(PHI(i)),0.,0.
    END DO
    DO i=1,MeshFS%Npanels
        WRITE(11,'(4(X,I7))') (MeshFS%P(j,i),j=1,4)
    END DO
    CLOSE(11)
    END SUBROUTINE WRITE_FS 
!
    SUBROUTINE WRITE_POTENTIAL(ID,NVEL)
    USE MIDENTIFICATION
    USE COM_VAR
    TYPE(TID) :: ID
    INTEGER:: i,NPT_temp,NFA_temp
    COMPLEX:: ZPBT(NP),ZPST(NP)
    REAL:: amp,phase,REZ,IMZ
    COMPLEX,DIMENSION(NP) :: NVELTEMPB,NVELTEMPS
!   Body condition    
    COMPLEX,DIMENSION(*) :: NVEL  
    DO I=1,NFA         
	    ZPBT(M1(I))=ZPB(I)
	    ZPBT(M2(I))=ZPB(I)
	    ZPBT(M3(I))=ZPB(I)
	    ZPBT(M4(I))=ZPB(I)
        ZPST(M1(I))=ZPS(I)
	    ZPST(M2(I))=ZPS(I)
	    ZPST(M3(I))=ZPS(I)
	    ZPST(M4(I))=ZPS(I) 
        NVELTEMPB(M1(I))=NVEL(I)
        NVELTEMPB(M2(I))=NVEL(I)
        NVELTEMPB(M3(I))=NVEL(I)
        NVELTEMPB(M4(I))=NVEL(I)
        NVELTEMPS(M1(I))=NVEL(I+NFA)
        NVELTEMPS(M2(I))=NVEL(I+NFA)
        NVELTEMPS(M3(I))=NVEL(I+NFA)
        NVELTEMPS(M4(I))=NVEL(I+NFA)
    END DO  
    IF (NSYMY .eq. 1) THEN
        NPT_temp=2*NP
        NFA_temp=2*NFA
    ELSE
        NPT_temp=NP
        NFA_temp=NFA
    END IF
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/potential.tec')
    WRITE(11,'(A)') 'VARIABLES="X" "Y" "Z" "Module" "Phase" "Real" "Imag" "Re(NVEL)" "Im(NVEL)"'
    WRITE(11,*) 'ZONE t="Body",N=',NPT_temp,', E=',NFA_temp,' , F=FEPOINT,ET=QUADRILATERAL'   
    DO i=1,NP
        amp=ABS(ZPBT(i))
        REZ=REAL(ZPBT(i))
        IMZ=AIMAG(ZPBT(i))
        phase=ATAN2(IMZ,REZ)
        WRITE(11,'(8(2X,E14.7))') X(i),Y(i),Z(i),amp,phase,REZ,IMZ,REAL(NVELTEMPB(I)),IMAG(NVELTEMPB(I))
    END DO
    IF (NSYMY .eq. 1) THEN
        DO i=1,NP        
	        REZ=REAL(ZPST(i))
	        IMZ=AIMAG(ZPST(i))
            phase=ATAN2(IMZ,REZ)
            amp=ABS(ZPST(i))
	        WRITE(11,'(8(2X,E14.7))') X(i),-Y(i),Z(i),amp,phase,REZ,IMZ,REAL(NVELTEMPS(I)),IMAG(NVELTEMPS(I))
        END DO
    END IF
    DO i=1,NFA
        WRITE(11,*)M1(I),M2(I),M3(I),M4(I)
    END DO
    IF (NSYMY .eq. 1) THEN
        DO i=1,NFA
	        WRITE(11,*)M1(I)+NP,M2(I)+NP,M3(I)+NP,M4(I)+NP
        END DO
    END IF
    CLOSE(11)    
    END SUBROUTINE
  
 END MODULE