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
    SUBROUTINE calCol(NFm,Xm,Ym,Zm,Facettem,Tcol,nFacemx)
!    
    IMPLICIT NONE   
!   Maillage de la partie mouille
    INTEGER NFm,nFacemx
    REAL,DIMENSION(*) :: Xm,Ym,Zm
    INTEGER,DIMENSION(4,*) :: Facettem
!   Hauteur de col
    REAL Tcol
!   Locales
    INTEGER :: nFacem2,Np
	REAL,DIMENSION(4,3,nFacemx) :: Coinm2	
    INTEGER :: i,j,k


   	nFacem2=NFm
	DO i=1,NFm
		DO j=1,4
			Coinm2(j,1,i)=Xm(Facettem(j,i))
			Coinm2(j,2,i)=Ym(Facettem(j,i))
			Coinm2(j,3,i)=Zm(Facettem(j,i))-Tcol
		END DO
	END DO
	DO i=1,Nfm
		IF (((Coinm2(1,3,i)+Tcol).GT.-1.0E-03).AND.((Coinm2(2,3,i)+Tcol).GT.-1.0E-03)) THEN
			Nfacem2=Nfacem2+1
			Coinm2(1,1,Nfacem2)=Coinm2(2,1,i)
			Coinm2(1,2,Nfacem2)=Coinm2(2,2,i)
			Coinm2(1,3,Nfacem2)=Coinm2(2,3,i)
			Coinm2(2,1,Nfacem2)=Coinm2(1,1,i)
			Coinm2(2,2,Nfacem2)=Coinm2(1,2,i)
			Coinm2(2,3,Nfacem2)=Coinm2(1,3,i)
			Coinm2(3,1,Nfacem2)=Coinm2(1,1,i)
			Coinm2(3,2,Nfacem2)=Coinm2(1,2,i)
			Coinm2(3,3,Nfacem2)=0.
			Coinm2(4,1,Nfacem2)=Coinm2(2,1,i)
			Coinm2(4,2,Nfacem2)=Coinm2(2,2,i)
			Coinm2(4,3,Nfacem2)=0.
		ELSE
			IF (((Coinm2(2,3,i)+Tcol).GT.-1.0E-03).AND.((Coinm2(3,3,i)+Tcol).GT.-1.0E-03)) THEN
				Nfacem2=Nfacem2+1
	    		Coinm2(1,1,Nfacem2)=Coinm2(3,1,i)
				Coinm2(1,2,Nfacem2)=Coinm2(3,2,i)
				Coinm2(1,3,Nfacem2)=Coinm2(3,3,i)
				Coinm2(2,1,Nfacem2)=Coinm2(2,1,i)
				Coinm2(2,2,Nfacem2)=Coinm2(2,2,i)
				Coinm2(2,3,Nfacem2)=Coinm2(2,3,i)
				Coinm2(3,1,Nfacem2)=Coinm2(2,1,i)
				Coinm2(3,2,Nfacem2)=Coinm2(2,2,i)
				Coinm2(3,3,Nfacem2)=0.
				Coinm2(4,1,Nfacem2)=Coinm2(3,1,i)
				Coinm2(4,2,Nfacem2)=Coinm2(3,2,i)
				Coinm2(4,3,Nfacem2)=0.
			ELSE
				IF (((Coinm2(3,3,i)+Tcol).GT.-1.0E-03).AND.((Coinm2(4,3,i)+Tcol).GT.-1.0E-03)) THEN
					Nfacem2=Nfacem2+1
					Coinm2(1,1,Nfacem2)=Coinm2(4,1,i)
					Coinm2(1,2,Nfacem2)=Coinm2(4,2,i)
					Coinm2(1,3,Nfacem2)=Coinm2(4,3,i)
					Coinm2(2,1,Nfacem2)=Coinm2(3,1,i)
					Coinm2(2,2,Nfacem2)=Coinm2(3,2,i)
					Coinm2(2,3,Nfacem2)=Coinm2(3,3,i)
					Coinm2(3,1,Nfacem2)=Coinm2(3,1,i)
					Coinm2(3,2,Nfacem2)=Coinm2(3,2,i)
					Coinm2(3,3,Nfacem2)=0.
					Coinm2(4,1,Nfacem2)=Coinm2(4,1,i)
					Coinm2(4,2,Nfacem2)=Coinm2(4,2,i)
					Coinm2(4,3,Nfacem2)=0.
				ELSE
					IF (((Coinm2(4,3,i)+Tcol).GT.-1.0E-03).AND.((Coinm2(1,3,i)+Tcol).GT.-1.0E-03)) THEN
						Nfacem2=Nfacem2+1
						Coinm2(1,1,Nfacem2)=Coinm2(1,1,i)
						Coinm2(1,2,Nfacem2)=Coinm2(1,2,i)
						Coinm2(1,3,Nfacem2)=Coinm2(1,3,i)
						Coinm2(2,1,Nfacem2)=Coinm2(4,1,i)
						Coinm2(2,2,Nfacem2)=Coinm2(4,2,i)
						Coinm2(2,3,Nfacem2)=Coinm2(4,3,i)
						Coinm2(3,1,Nfacem2)=Coinm2(4,1,i)
						Coinm2(3,2,Nfacem2)=Coinm2(4,2,i)
						Coinm2(3,3,Nfacem2)=0.
						Coinm2(4,1,Nfacem2)=Coinm2(1,1,i)
						Coinm2(4,2,Nfacem2)=Coinm2(1,2,i)
						Coinm2(4,3,Nfacem2)=0.
					END IF
				END IF
			END IF
		END IF
	END DO
	Np=0
	NFm=Nfacem2
	DO i=1,nFacem2
		DO j=1,4
			Xm(Np+j)=Coinm2(j,1,i)
			Ym(Np+j)=Coinm2(j,2,i)
			Zm(Np+j)=Coinm2(j,3,i)
			Facettem(j,i)=Np+j
		END DO
		Np=Np+4
	END DO
    
    END SUBROUTINE