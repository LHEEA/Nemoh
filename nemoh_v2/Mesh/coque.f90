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
	SUBROUTINE coque(X,Y,Z,NP,facettes,NF,Deplacement,Icoque,Gcoque,CG,Nsym,rho)

	IMPLICIT NONE

!	Globale
	REAL,PARAMETER :: PI=3.141592653589
!	Description du maillage
	INTEGER Np,Nf
	REAL,DIMENSION(*) :: X,Y,Z
	INTEGER,DIMENSION(4,*) :: Facettes
	REAL Deplacement
	INTEGER Nsym
	REAL :: rho
!	Caracteristiques de la coque
	REAL,DIMENSION(3) :: Gcoque,CG
	REAL mcoque
	REAL,DIMENSION(3,3) :: Icoque
	REAL decoque
!	Locales
	REAL,DIMENSION(3) :: U,V,W
	REAL,DIMENSION(3,Nf) :: CdG
	REAL :: N1,N2
	REAL,DIMENSION(Nf) :: Aire
!	Indices
	INTEGER i,j

	mcoque=0.
	Gcoque(1)=0.
	Gcoque(2)=0.
	Gcoque(3)=0.
	DO i=1,nF
		CdG(1,i)=0.25*(X(Facettes(1,i))+X(Facettes(2,i))+X(Facettes(3,i))+X(Facettes(4,i)))
		CdG(2,i)=0.25*(Y(Facettes(1,i))+Y(Facettes(2,i))+Y(Facettes(3,i))+Y(Facettes(4,i)))
		CdG(3,i)=0.25*(Z(Facettes(1,i))+Z(Facettes(2,i))+Z(Facettes(3,i))+Z(Facettes(4,i)))
		U(1)=X(Facettes(2,i))-X(Facettes(1,i))
		U(2)=Y(Facettes(2,i))-Y(Facettes(1,i))
		U(3)=Z(Facettes(2,i))-Z(Facettes(1,i))
		V(1)=X(Facettes(3,i))-X(Facettes(1,i))
		V(2)=Y(Facettes(3,i))-Y(Facettes(1,i))
		V(3)=Z(Facettes(3,i))-Z(Facettes(1,i))
		CALL prdvct(U,V,W)
		N1=0.5*SQRT(W(1)*W(1)+W(2)*W(2)+W(3)*W(3))
		U(1)=X(Facettes(3,i))-X(Facettes(1,i))
		U(2)=Y(Facettes(3,i))-Y(Facettes(1,i))
		U(3)=Z(Facettes(3,i))-Z(Facettes(1,i))
		V(1)=X(Facettes(4,i))-X(Facettes(1,i))
		V(2)=Y(Facettes(4,i))-Y(Facettes(1,i))
		V(3)=Z(Facettes(4,i))-Z(Facettes(1,i))
		CALL prdvct(U,V,W)
		N2=0.5*SQRT(W(1)*W(1)+W(2)*W(2)+W(3)*W(3))
		Aire(i)=N1+N2
		mcoque=mcoque+Aire(i)
		Gcoque(1)=Gcoque(1)+Aire(i)*CdG(1,i)
		Gcoque(2)=Gcoque(2)+Aire(i)*CdG(2,i)
		Gcoque(3)=Gcoque(3)+Aire(i)*CdG(3,i)
	END DO
	decoque=Deplacement*rho/mcoque
	Gcoque(1)=Gcoque(1)/mcoque
	IF (Nsym.EQ.1) THEN
	    Gcoque(2)=0.
	ELSE
	    Gcoque(2)=Gcoque(2)/mcoque
	END IF
	Gcoque(3)=Gcoque(3)/mcoque
	DO i=1,3
		DO j=1,3
			Icoque(i,j)=0.
		END DO
	END DO
	DO i=1,nF
!		Icoque(1,1)=Icoque(1,1)+Aire(i)*decoque*((CdG(2,i)-Gcoque(2))**2+(CdG(3,i)-Gcoque(3))**2)
!		Icoque(1,2)=Icoque(1,2)-Aire(i)*decoque*((CdG(1,i)-Gcoque(1))*(CdG(2,i)-Gcoque(2)))
!		Icoque(1,3)=Icoque(1,3)-Aire(i)*decoque*((CdG(1,i)-Gcoque(1))*(CdG(3,i)-Gcoque(3)))
!		Icoque(2,2)=Icoque(2,2)+Aire(i)*decoque*((CdG(1,i)-Gcoque(1))**2+(CdG(3,i)-Gcoque(3))**2)
!		Icoque(2,3)=Icoque(2,3)-Aire(i)*decoque*((CdG(2,i)-Gcoque(2))*(CdG(3,i)-Gcoque(3)))
!		Icoque(3,3)=Icoque(3,3)+Aire(i)*decoque*((CdG(1,i)-Gcoque(1))**2+(CdG(2,i)-Gcoque(2))**2)
		Icoque(1,1)=Icoque(1,1)+Aire(i)*decoque*((CdG(2,i)-Gcoque(2))**2+(CdG(3,i)-CG(3))**2)
		Icoque(1,2)=Icoque(1,2)-Aire(i)*decoque*((CdG(1,i)-Gcoque(1))*(CdG(2,i)-CG(2)))
		Icoque(1,3)=Icoque(1,3)-Aire(i)*decoque*((CdG(1,i)-Gcoque(1))*(CdG(3,i)-CG(3)))
		Icoque(2,2)=Icoque(2,2)+Aire(i)*decoque*((CdG(1,i)-Gcoque(1))**2+(CdG(3,i)-CG(3))**2)
		Icoque(2,3)=Icoque(2,3)-Aire(i)*decoque*((CdG(2,i)-Gcoque(2))*(CdG(3,i)-CG(3)))
		Icoque(3,3)=Icoque(3,3)+Aire(i)*decoque*((CdG(1,i)-Gcoque(1))**2+(CdG(2,i)-CG(2))**2)
	END DO
	IF (Nsym.EQ.1) THEN
	    Icoque(1,2)=0.
	    Icoque(2,3)=0.    
	END IF
	!
	Icoque(2,1)=Icoque(1,2)
	Icoque(3,1)=Icoque(1,3)
	Icoque(3,2)=Icoque(2,3)

	RETURN

	END SUBROUTINE
	
