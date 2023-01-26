!-------------------------------------------------------------------------------------------------
!
!  NEMOH v1.0 - Mesh generator - January 14
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

	PROGRAM Mesh

    USE MIdentification
#ifndef GNUFORT
    USE iflport
#endif

	IMPLICIT NONE
	
!	Working directory and name of file describing the geometry
	TYPE(TID) :: ID,DSCRPT
	LOGICAL existence
!   Nombre de symétries du flotteur
	INTEGER :: Nsym
	INTEGER :: nXsym,nYsym
!   Position du maillage dans le plan (x0y)
	REAL :: tX,tY
!   Description du flotteur en faces
	INTEGER,PARAMETER :: nFacemx=3000
!   Description de tout le flotteur
	INTEGER :: nFace						! Nombre reel de faces
	REAL,DIMENSION(4,3,nFacemx) :: Coin		! Coins des faces
!   Description de la partie immergee
	INTEGER :: nFacem						! Nombre reel de faces
	REAL,DIMENSION(4,3,nFacemx) :: Coinm	! Coins des faces
!   Parametres de reglage du maillage
	INTEGER  nmaille,nmaillage
	REAL maille
	REAL Tcol
	REAL :: lambda                                  ! Facteur d echelle
!   Maillage proprement dit
	INTEGER,PARAMETER :: NFMX=20000					! Nombre de facettes max
	INTEGER,PARAMETER :: NPMX=20000					! Nombre de points max
	INTEGER :: Nmailmx		! Nombre de facettes du maillage std max
!   Maillage du corps
	INTEGER :: NF,NP
	INTEGER,DIMENSION(4,NFMX) :: Facette
	REAL,DIMENSION(NPMX) :: X,Y,Z
!   Partie immergee du maillage
	INTEGER :: NFm,NPm
	INTEGER,DIMENSION(4,NFMX) :: Facettem
	REAL,DIMENSION(NPMX) :: Xm,Ym,Zm
!   Calcul hydrostatique
	REAL DEPLACEMENT,XF,YF,ZF,SF
	REAL,DIMENSION(6,6) :: KH
	REAL :: xG,yG,zG
	REAL :: RHO,G
!   Calcul coque
	REAL,DIMENSION(3,3) :: Icoque
	REAL,DIMENSION(3) :: Gcoque,CDG
!
	INTEGER i,j,indx,k,c,cont
	INTEGER,DIMENSION(4) :: p
	CHARACTER*2 :: Num
	INTEGER :: cNum

	WRITE(*,*)
	WRITE(*,*) ' -> Read input data '
	WRITE(*,*) ' '
	CALL ReadTID(ID,'ID.dat') 
	OPEN(10,FILE='Mesh.cal')
	READ(10,*) DSCRPT%ID
	DSCRPT%lID=LNBLNK(DSCRPT%ID)
	READ(10,*) Nsym
	READ(10,*) tX,tY
	READ(10,*) xG,yG,zG
	READ(10,*) Nmailmx
	CLOSE(10)
	OPEN(10,file=ID%ID(1:ID%lID)//'/mesh/'//DSCRPT%ID(1:DSCRPT%lID))
	READ(10,*) Np
	READ(10,*) nFace
	DO i=1,Np
	    READ(10,*) X(i),Y(i),Z(i)
	END DO
	DO j=1,nFace
	    READ(10,*) p(1),P(2),P(3),p(4)
		DO i=1,4
			Coin(i,1,j)=1.0*X(p(i))
			Coin(i,2,j)=1.0*Y(p(i))
			Coin(i,3,j)=1.0*Z(p(i))
		END DO
	END DO
	CLOSE(10)
!   Prise en compte d'une collerette ?
	OPEN(10,FILE='Mesh.cal')
	DO i=1,6
	    READ(10,*)
	END DO
	READ(10,*) Tcol
	READ(10,*) lambda
	READ(10,*) RHO
	READ(10,*) G
	CLOSE(10)
!   Mise a l echelle
    IF ((lambda).NE.(1.)) THEN
        WRITE(*,'(A,F7.3)') '   - Mesh is scaled by ', lambda
        xG=xG*lambda
        yG=yG*lambda
        zG=zG*lambda
        tX=tX*lambda
        tY=tY*lambda
        DO i=1,nFace
		    DO j=1,4
			    Coin(j,1,i)=Coin(j,1,i)*lambda
			    Coin(j,2,i)=Coin(j,2,i)*lambda
			    Coin(j,3,i)=Coin(j,3,i)*lambda
			END DO
		END DO
    END IF
!   Translate
    xG=xG+tX
    yG=yG+tY
    DO i=1,nFace
		    DO j=1,4
			    Coin(j,1,i)=Coin(j,1,i)+tX
			    Coin(j,2,i)=Coin(j,2,i)+tY
			END DO
		END DO
    WRITE(*,*) ' -> Calculate hydrostatics '
    WRITE(*,*) ' '     
	Np=0
	Nf=nFace
	DO i=1,nFace
	    DO j=1,4
		    X(Np+j)=Coin(j,1,i)-xG
			Y(Np+j)=Coin(j,2,i)-yG
			Z(Np+j)=Coin(j,3,i)+Tcol
			Facette(j,i)=Np+j
		END DO
		Np=Np+4
	END DO	
	CALL HYDRO(X,Y,Z,NP,FACETTE,NF,DEPLACEMENT,XF,YF,ZF,SF,KH,Xm,Ym,Zm,NPm,FACETTEm,NFm,RHO,G)
	IF (Tcol.GT.0) THEN
		CALL calCol(NFm,Xm,Ym,Zm,Facettem,Tcol,nFacemx)
    END IF
	nFacem=NFm
	DO i=1,NFm
		DO j=1,4
			Coinm(j,1,i)=Xm(Facettem(j,i))+xG
			Coinm(j,2,i)=Ym(Facettem(j,i))+yG
			Coinm(j,3,i)=Zm(Facettem(j,i))
		END DO
	END DO
!	Prise en compte de la symétrie
	IF (Nsym.EQ.1) THEN
		DEPLACEMENT=2.0*DEPLACEMENT
		YF=0.
		SF=2.0*SF
	END IF
	WRITE(*,'(A,I3)') '   - Coordinates of buoyancy centre '
	WRITE(*,'(A,F7.3,A)') '     XB = ',XF+xG,'m'
	WRITE(*,'(A,F7.3,A)') '     YB = ',YF+yG,'m'
	WRITE(*,'(A,F7.3,A)') '     ZB = ',ZF,'m'
	WRITE(*,'(A,E14.7,A)') '    - Displacement  = ',DEPLACEMENT,' m^3'
	WRITE(*,'(A,E14.7,A)') '    - Waterplane area  = ',SF, ' m^2'
	WRITE(*,*) ' '
	IF ((ABS(XF).GT.1E-02).OR.(ABS(YF).GT.1E-02)) THEN
	    WRITE(*,*) ' '
	    WRITE(*,*) ' !!! WARNING !!! '
		WRITE(*,*) ' '
		WRITE(*,'(A,I3)') ' Buoyancy center and gravity center are not vertically aligned. '
		WRITE(*,*) ' This is not an equilibrium position.'
		WRITE(*,'(A,F7.3,1X,A,F7.3)') ' XF = ',XF+xG,' XG = ',xG
		WRITE(*,'(A,F7.3,1X,A,F7.3)') ' YF = ',YF+yG,' YG = ',yG
	END IF
	OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/Hydrostatics.dat')
	WRITE(10,'(A,F7.3,A,F7.3)') ' XF = ',XF+xG,' - XG = ',xG
	WRITE(10,'(A,F7.3,A,F7.3)') ' YF = ',YF+yG,' - YG = ',yG
	WRITE(10,'(A,F7.3,A,F7.3)') ' ZF = ',ZF,' - ZG = ',zG
	WRITE(10,'(A,E14.7)') ' Displacement = ',DEPLACEMENT
	WRITE(10,'(A,E14.7)') ' Waterplane area = ',SF
	CLOSE(10)
!	Sauvegarde de la description du maillage et de la partie immergee
	OPEN(10,file=ID%ID(1:ID%lID)//'/mesh/Description_Full.tec')
	WRITE(10,*) 'ZONE N=',4*nFace,',E=',nFace,',F=FEPOINT,ET=QUADRILATERAL'
	DO i=1,nFace
		DO j=1,4
			WRITE(10,'(6(2X,E14.7))') Coin(j,1,i),Coin(j,2,i),Coin(j,3,i),0.,0.,0.
		END DO
	END DO
	DO i=1,nFace
		WRITE(10,'(I4,3(2X,I4))') 1+(i-1)*4,2+(i-1)*4,3+(i-1)*4,4+(i-1)*4
	END DO
	CLOSE(10)
	OPEN(10,file=ID%ID(1:ID%lID)//'/mesh/Description_Wetted.tec')
	WRITE(10,*) 'ZONE N=',4*nFacem,',E=',nFacem,',F=FEPOINT,ET=QUADRILATERAL'
	do i=1,nFacem
		do j=1,4
			WRITE(10,'(6(2X,E14.7))') Coinm(j,1,i),Coinm(j,2,i),Coinm(j,3,i),0.,0.,0.
		end do
	end do
	do i=1,nFacem
		WRITE(10,'(I4,3(2X,I4))') 1+(i-1)*4,2+(i-1)*4,3+(i-1)*4,4+(i-1)*4
	end do
	CLOSE(10)
	write(*,*) ' -> Make mesh '
	WRITE(*,*) ' '
	maille=5.0			
	nf=2*Nmailmx
	nmaillage=0
	do while (((nf.GT.Nmailmx).OR.(nf.LT.(0.9*Nmailmx))).AND.(nmaillage.LT.100))
		call Maillage(nFacem,Coinm,maille,X,Y,Z,NP,facette,NF)
		if (nf.GT.Nmailmx) then
			maille=maille*sqrt(1.0*nf/Nmailmx)
		else 
			maille=maille*sqrt(1.0*nf/(0.9*Nmailmx))
		end if
		nmaillage=nmaillage+1
	end do
	WRITE(*,'(A,I5)') '    Number of panels in final mesh : ',NF
	call ExMaillage(ID,DSCRPT,X,Y,Z,NP,NPMX,facette,NF,NFMX,NSYM)
!	Calcul et sauvegarde de la matrice de raideur hydrostatique
	DO j=1,NP
		X(j)=X(j)-xG
		Y(j)=Y(j)-yG
	END DO
	CALL HYDRO(X,Y,Z,NP,FACETTE,NF,DEPLACEMENT,XF,YF,ZF,SF,KH,Xm,Ym,Zm,NPm,FACETTEm,NFm,RHO,G)
	DO j=1,NP
		X(j)=X(j)+xG
		Y(j)=Y(j)+yG
	END DO
	IF (Nsym.EQ.1) THEN
		DEPLACEMENT=2.0*DEPLACEMENT
		YF=0.
		SF=2.0*SF
		KH(3,3)=2.*KH(3,3)
		KH(3,4)=0.
		KH(4,3)=0.
		KH(3,5)=2.*KH(3,5)
		KH(5,3)=KH(3,5)
		KH(4,4)=2.*KH(4,4)
		KH(4,5)=0.
		KH(5,4)=0.
		KH(5,5)=2.*KH(5,5)
	END IF
	KH(4,4)=KH(4,4)+deplacement*RHO*G*(ZF-ZG)
	KH(5,5)=KH(5,5)+deplacement*RHO*G*(ZF-ZG)
	OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/KH.dat')
	DO i=1,6
		WRITE(10,'(6(1X,E14.7))') (KH(i,j),j=1,6)
	END DO
	CLOSE(10)
	write(*,*) ' -> Calculate hull mass and inertia '
	WRITE(*,*) ' '
	CDG(1)=xG
	CDG(2)=yG
	CDG(3)=zG	
    CALL coque(X,Y,Z,NP,facette,NF,Deplacement,Icoque,Gcoque,CDG,Nsym,rho)
	OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/GC_hull.dat')
	WRITE(10,'(3(1X,E14.7))') Gcoque(1),Gcoque(2),Gcoque(3)
	CLOSE(10)
	OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/Inertia_hull.dat')
	DO i=1,3
		WRITE(10,'(3(1X,E14.7))') (Icoque(i,j),j=1,3)
	END DO
	CLOSE(10)

	end program Mesh

