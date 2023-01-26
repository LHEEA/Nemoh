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
	SUBROUTINE ExMaillage(ID,DSCRPT,X,Y,Z,NP,NPMX,facette,NF,NFMX,NSYM)
	
	USE MIdentification

	implicit none

	! Entrees
	TYPE(TID) :: ID,DSCRPT
	INTEGER NP,NF,NPMX,NFMX
	REAL,DIMENSION(*) :: X,Y,Z
	INTEGER,DIMENSION(4,*) :: facette
	INTEGER :: NSYM
	! Locales
	REAL,DIMENSION(NFMX) :: Xn,Yn,Zn	! Normale a la facette
	REAL,DIMENSION(NFMX) :: Xg,Yg,Zg	! CdG des facettes
	REAL,DIMENSION(NFMX) :: Aire		! Aire des facettes
	REAL,DIMENSION(NFMX) :: Dist		! Plus grande longueur des facettes
	INTEGER i,j,k
	REAL T 								! Tirant d eau (valeur adimensionnalisation ACHIL3D)
	REAL,DIMENSION(3) :: u,v,w
	REAL norme,norme1,norme2,norme3
	
	! Recherche du tirant d eau
	T=0.0
	do i=1,NP
		IF (ABS(Z(i)).GT.T) T=ABS(Z(i))
	end do
	! Calcul des normales
	do i=1,NF
		u(1)=X(facette(2,i))-X(facette(1,i))
		u(2)=Y(facette(2,i))-Y(facette(1,i))
		u(3)=Z(facette(2,i))-Z(facette(1,i))
		norme1=u(1)**2+u(2)**2+u(3)**2
		v(1)=X(facette(3,i))-X(facette(1,i))
		v(2)=Y(facette(3,i))-Y(facette(1,i))
		v(3)=Z(facette(3,i))-Z(facette(1,i))
		norme2=v(1)**2+v(2)**2+v(3)**2
		call prodvect(u,v,w)
		Xn(i)=w(1)
		Yn(i)=w(2)
		Zn(i)=w(3)
		u(1)=X(facette(3,i))-X(facette(1,i))
		u(2)=Y(facette(3,i))-Y(facette(1,i))
		u(3)=Z(facette(3,i))-Z(facette(1,i))
		v(1)=X(facette(4,i))-X(facette(1,i))
		v(2)=Y(facette(4,i))-Y(facette(1,i))
		v(3)=Z(facette(4,i))-Z(facette(1,i))
		norme3=v(1)**2+v(2)**2+v(3)**2
		call prodvect(u,v,w)
		Xn(i)=0.5*(w(1)+Xn(i))
		Yn(i)=0.5*(w(2)+Yn(i))
		Zn(i)=0.5*(w(3)+Zn(i))
		Dist(i)=max(norme1,norme2,norme3)
		Aire(i)=sqrt(Xn(i)**2+Yn(i)**2+Zn(i)**2)
		norme=Aire(i)
		Xn(i)=Xn(i)/norme
		Yn(i)=Yn(i)/norme
		Zn(i)=Zn(i)/norme
		Xg(i)=0.25*(X(facette(1,i))+X(facette(2,i))+X(facette(3,i))+X(facette(4,i)))
		Yg(i)=0.25*(Y(facette(1,i))+Y(facette(2,i))+Y(facette(3,i))+Y(facette(4,i)))
		Zg(i)=0.25*(Z(facette(1,i))+Z(facette(2,i))+Z(facette(3,i))+Z(facette(4,i)))
	end do	
	! Maillage AQUAPLUS
	! Fichier Tecplot
	open(10,file=ID%ID(1:ID%lID)//'/mesh/'//DSCRPT%ID(1:DSCRPT%lID)//'.tec')
	write(10,*) 'ZONE N=',np,', E=',nf,' , F=FEPOINT,ET=QUADRILATERAL'
	do i=1,np
		write(10,'(6(2X,E14.7))') X(i),Y(i),Z(i),0.0,0.0,0.0
	end do
	do i=1,nf
		write(10,'(I6,3(2X,I6))') facette(1,i),facette(2,i),facette(3,i),facette(4,i)
	end do
	write(10,*) 'ZONE t="normales", F=POINT, I=',nf
	do i=1,nf
		write(10,'(6(2X,E14.7))') Xg(i),Yg(i),Zg(i),Xn(i),Yn(i),Zn(i) 
	end do
	close(10)
	! Fichier de maillage
	open(10,file=ID%ID(1:ID%lID)//'/mesh/'//DSCRPT%ID(1:DSCRPT%lID)//'.dat')
	write(10,'(20X,I1,10X,I1)') 2,NSYM
	do i=1,np
		write(10,'(10X,I4,3(10X,F14.7))') i,X(i),Y(i),Z(i)
	end do
	write(10,'(10X,I4,3(10X,F4.2))') 0,0.,0.,0.
	do i=1,nf
		write(10,'(4(10X,I6))') facette(1,i),facette(2,i),facette(3,i),facette(4,i)
	end do
	write(10,'(4(10X,I1))') 0,0,0,0
	close(10)
	! Fichiers de maillage
	open(10,file=ID%ID(1:ID%lID)//'/mesh/'//DSCRPT%ID(1:DSCRPT%lID)//'_info.dat')
	write(10,'(I7,X,I7,A)') np,nf, ' Number of points and number of panels'
	close(10)
	
	RETURN

	end 

!*******************************************************************
!
!	Calcul du produit vectoriel w = u ^ v
!
!******************************************************************

	SUBROUTINE prodvect(u,v,w)

	IMPLICIT NONE

	REAL,DIMENSION(3) :: u,v,w

	w(1)=u(2)*v(3)-u(3)*v(2)
	w(2)=u(3)*v(1)-u(1)*v(3)
	w(3)=u(1)*v(2)-u(2)*v(1)

	RETURN

	END