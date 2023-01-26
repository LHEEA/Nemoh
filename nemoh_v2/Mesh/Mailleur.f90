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
	SUBROUTINE MAILLAGE(nFace,Coin,maille,X,Y,Z,NP,facette,NF)

	implicit none

	! Entrees
	REAL maille					!   Longueur maximale de maille
	INTEGER nFace
	REAL,DIMENSION(4,3,*) :: Coin
	! Sorties
	INTEGER NP,NF
	REAL,DIMENSION(*) :: X,Y,Z
	INTEGER,DIMENSION(4,*) :: facette
	! Locales
	INTEGER i,j,k
	INTEGER,PARAMETER :: ndpmx=10000,ndfmx=10000
	INTEGER :: ndp,ndf			!	ndp : nombre de points
								!	nf : nombre de facettes
	REAL,DIMENSION(3,ndpmx) :: noeud	!	noeuds : du maillage
	INTEGER,DIMENSION(4,ndfmx) :: face	!	faces : du maillage
	INTEGER,DIMENSION(ndfmx) :: indx	!	Index de maillage
	REAL,DIMENSION(3) :: X1,X2,X3,X4
	REAL norme1,norme2
	
	! Analyse des faces. Si une face est trop petite par rapport a la taille de la maille
	! elle ne sera pas maillee�et constituera une facette a�part entiere.
	do i=1,Nface
		do j=1,3
			X1(j)=Coin(1,j,i)
			X2(j)=Coin(2,j,i)
			X3(j)=Coin(3,j,i)
			X4(j)=Coin(4,j,i)
		end do
		call nprodvect(X1,X2,X3,norme1)
		call nprodvect(X1,X3,X4,norme2)	
!		if (0.5*(norme1+norme2).LT.maille) then
!			indx(i)=0
!		else
			indx(i)=1
!		end if
	end do
	! Calcul du maillage
	np=0
	nf=0
	! Maillage des faces 
	do i=1,Nface		
		if (indx(i).EQ.1) then
			do j=1,3
				X1(j)=Coin(1,j,i)
				X2(j)=Coin(2,j,i)
				X3(j)=Coin(3,j,i)
				X4(j)=Coin(4,j,i)
			end do
			call plan(X2,X3,X4,X1,maille,noeud,ndp,face,ndf)
			do j=1,ndp
				X(j+np)=noeud(1,j)
				Y(j+np)=noeud(2,j)
				Z(j+np)=noeud(3,j)
			end do
			do j=1,ndf
				do k=1,4
					facette(k,j+nf)=face(k,j)+np
				end do
			end do
			np=np+ndp
			nf=nf+ndf
		else
			do j=1,4
				X(j+np)=Coin(j,1,i)
				Y(j+np)=Coin(j,2,i)
				Z(j+np)=Coin(j,3,i)
			end do
			do k=1,4
				facette(k,1+nf)=np+k
			end do
			np=np+4
			nf=nf+1
		end if
	end do

	RETURN

	end 
!*************************************************************
!
!	Cette subroutine realise le maillage d'une face 
!	sommmets X1,X2,X3,X4 avec des mailles de dimension max maille.
!		X3 peut etre confondu avec X4
!		Les noeuds sont stockes dans noeud.
!		Les facettes dans facette.
!
!**************************************************************

	SUBROUTINE plan(X1,X2,X3,X4,maille,noeud,n,facette,nf)

	INTEGER i,j,k
	REAL X1,X2,X3,X4,Pi,Pf
	DIMENSION X1(3),X2(3),X3(3),X4(3),Pi(3),Pf(3)
	REAL noeud,maille
	INTEGER,DIMENSION(4,*) :: facette
	DIMENSION noeud(3,*)
	REAL norme1,norme2,norme3,norme
	REAL angle1,angle2,angle3
	REAL A1(3),A2(3),A3(3)
!	REAL,PARAMETER :: njmn=3		! Nombre minimum de points de discretisation
	INTEGER :: njmn

	OPEN(10,FILE='Mesh.cal')
	DO i=1,5
	    READ(10,*)
	END DO
	READ(10,*) njmn
	CLOSE(10)
	call nprodvect(X1,X2,X3,norme1)
	call nprodvect(X1,X2,X4,norme2)
	call nprodvect(X1,X3,X4,norme3)
	if ((norme1.EQ.0.0).AND.(norme2.EQ.0.0).AND.(norme3.EQ.0.0)) then
		nf=0
		n=0
	else
		if ((X3(1).EQ.X4(1)).AND.(X3(2).EQ.X4(2)).AND.(X3(3).EQ.X4(3))) then
			norme1=((X2(1)-X1(1))**2+(X2(2)-X1(2))**2+(X2(3)-X1(3))**2)**0.5
			norme2=((X3(1)-X2(1))**2+(X3(2)-X2(2))**2+(X3(3)-X2(3))**2)**0.5
			norme3=((X1(1)-X3(1))**2+(X1(2)-X3(2))**2+(X1(3)-X3(3))**2)**0.5
			Angle1=((X2(1)-X1(1))*(X3(3)-X1(3))+(X2(3)-X1(3))*(X3(3)-X1(3))+(X2(3)-X1(3))*(X3(3)-X1(3)))/norme1/norme3
			Angle2=((X3(1)-X2(1))*(X1(3)-X2(3))+(X3(3)-X2(3))*(X1(3)-X2(3))+(X3(3)-X2(3))*(X1(3)-X2(3)))/norme2/norme1
			Angle3=((X1(1)-X3(1))*(X2(3)-X3(3))+(X1(3)-X3(3))*(X2(3)-X3(3))+(X1(3)-X3(3))*(X2(3)-X3(3)))/norme3/norme2
			if ((angle1.GT.angle2).AND.(angle1.GT.angle3)) then
				do i=1,3
					A1(i)=X1(i)
					A2(i)=X2(i)
					A3(i)=X3(i)
				end do
			else 
				if (angle2.GT.angle3) then
					do i=1,3
						A1(i)=X2(i)
						A2(i)=X3(i)
						A3(i)=X1(i)
					end do
				else
					do i=1,3
						A1(i)=X3(i)
						A2(i)=X1(i)
						A3(i)=X2(i)
					end do
				end if
			end if
			norme1=((A2(1)-A1(1))**2+(A2(2)-A1(2))**2+(A2(3)-A1(3))**2)**0.5
			norme2=((A3(1)-A2(1))**2+(A3(2)-A2(2))**2+(A3(3)-A2(3))**2)**0.5
			norme3=((A1(1)-A3(1))**2+(A1(2)-A3(2))**2+(A1(3)-A3(3))**2)**0.5
			if (norme1.LT.norme3) then
				nj=int(norme1/maille)+njmn
				n=0
				do j=1,nj
					do k=1,3
						Pi(k)=A1(k)+(A2(k)-A1(k))*(j-1)/(nj-1)
					end do
					norme=((Pi(1)-A1(1))**2+(Pi(2)-A1(2))**2+(Pi(3)-A1(3))**2)**0.5
					do k=1,3				
						Pf(k)=A1(k)+(A3(k)-A1(k))*norme/norme1
					end do
					do i=1,j
						n=n+1
						do k=1,3
							if (j.GT.1) then
								noeud(k,n)=Pi(k)+(Pf(k)-Pi(k))*(i-1)/(j-1)
							else
								noeud(k,n)=Pi(k)
							end if
						end do
					end do
				end do
				nf=0	
				do j=1,nj-1
					do i=1,j-1
						nf=nf+1
						facette(1,nf)=nf
						facette(2,nf)=nf+j
						facette(3,nf)=nf+j+1
						facette(4,nf)=nf+1
					end do
					nf=nf+1
					facette(1,nf)=nf
					facette(2,nf)=nf+j
					facette(3,nf)=nf+j+1
					facette(4,nf)=nf
				end do			
			else
				nj=int(norme3/maille)+njmn
				n=0
				do j=1,nj
					do k=1,3
						Pi(k)=A1(k)+(A3(k)-A1(k))*(j-1)/(nj-1)
					end do
					norme=((Pi(1)-A1(1))**2+(Pi(2)-A1(2))**2+(Pi(3)-A1(3))**2)**0.5
					do k=1,3				
						Pf(k)=A1(k)+(A2(k)-A1(k))*norme/norme3
					end do
					do i=1,j
						n=n+1
						do k=1,3
							if (j.GT.1) then
								noeud(k,n)=Pi(k)+(Pf(k)-Pi(k))*(i-1)/(j-1)
							else
								noeud(k,n)=Pi(k)
							end if
						end do
					end do
				end do
				nf=0	
				do j=1,nj-1
					do i=1,j-1
						nf=nf+1
						facette(1,nf)=nf
						facette(4,nf)=nf+j
						facette(3,nf)=nf+j+1
						facette(2,nf)=nf+1
					end do
					nf=nf+1
					facette(1,nf)=nf
					facette(3,nf)=nf+j
					facette(2,nf)=nf+j+1
					facette(4,nf)=nf
				end do	
			end if
		else
			ni=int(((X2(1)-X1(1))**2+(X2(2)-X1(2))**2+(X2(3)-X1(3))**2)**0.5/maille)+njmn
			nj=int(((X4(1)-X1(1))**2+(X4(2)-X1(2))**2+(X4(3)-X1(3))**2)**0.5/maille)+njmn
			n=0
			do j=1,nj
				do k=1,3
					Pi(k)=X1(k)+(X4(k)-X1(k))*(j-1)/(nj-1)
					Pf(k)=X2(k)+(X3(k)-X2(k))*(j-1)/(nj-1)
				end do		
				do i=1,ni
					n=n+1
					do k=1,3
						noeud(k,n)=Pi(k)+(Pf(k)-Pi(k))*(i-1)/(ni-1)
					end do
				end do
			end do
			nf=0	
			do j=1,nj-1
				do i=1,ni-1
					nf=nf+1
					facette(1,nf)=i+(j-1)*ni
					facette(2,nf)=i+1+(j-1)*ni
					facette(4,nf)=i+j*ni
					facette(3,nf)=i+1+j*ni
				end do
			end do
		end if
	end if

	END 

!*******************************************************************
!
!	Calcule la norme du  produit vectoriel des vecteurs X1X2 et X1X3
!
!******************************************************************

	SUBROUTINE nprodvect(X1,X2,X3,norme)

	IMPLICIT NONE

	REAL X1(3),X2(3),X3(3)
	REAL N(3)
	REAL norme

	N(1)=(X2(2)-X1(2))*(X3(3)-X1(3))-(X2(3)-X1(3))*(X3(2)-X1(2))
	N(2)=(X2(3)-X1(3))*(X3(1)-X1(1))-(X2(1)-X1(1))*(X3(3)-X1(3))
	N(3)=(X2(1)-X1(1))*(X3(2)-X1(2))-(X2(2)-X1(2))*(X3(1)-X1(1))
	norme=sqrt(N(1)**2+N(2)**2+N(3)**2)

	RETURN

	END