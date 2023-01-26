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
	SUBROUTINE HYDRO(X,Y,Z,NP,FACETTES,NF,DEPLACEMENT,XF,YF,ZF,SF,KH,Xmouillee,Ymouillee,Zmouillee&
,Npmouillee,Facettesmouillee,Nfmouillee,RHO,G)

	IMPLICIT NONE

!	DESCRIPTION SUR LE MAILLAGE
	INTEGER NP,NF
	REAL,DIMENSION(*) :: X,Y,Z
	INTEGER,DIMENSION(4,*) :: FACETTES
!	Resultats hydrostatique
	REAL :: DEPLACEMENT
	REAL ::XF,YF,ZF
	REAL :: SF
    REAL :: RHO,G
!	Maillage de la surface mouille
	INTEGER Npmouillee,Nfmouillee
	REAL,DIMENSION(*) :: Xmouillee,Ymouillee,Zmouillee
	INTEGER,DIMENSION(4,*) :: Facettesmouillee
!	LOCALES
	REAL,DIMENSION(4,3) :: P
	REAL,DIMENSION(3) :: P1,P2,P3,P4,PGF,PL,PK
	REAL VLF,ZMN,C,ZMX,SEF
	REAL,DIMENSION(6,6) :: KHe,KH
	INTEGER I,J,K,L,IMN
	INTEGER,PARAMETER :: ITEC=0

	DEPLACEMENT=0.0
	XF=0.0
	YF=0.0
	ZF=0.0
	SF=0.0
	KH=0.0
	Npmouillee=0
	Nfmouillee=0
	IF (ITEC.EQ.1) THEN
		OPEN(10,FILE='SW.TEC')
	END IF
	DO I=1,NF
		P(1,1)=X(FACETTES(1,I))
		P(1,2)=Y(FACETTES(1,I))
		P(1,3)=Z(FACETTES(1,I))
		P(2,1)=X(FACETTES(2,I))
		P(2,2)=Y(FACETTES(2,I))
		P(2,3)=Z(FACETTES(2,I))
		P(3,1)=X(FACETTES(3,I))
		P(3,2)=Y(FACETTES(3,I))
		P(3,3)=Z(FACETTES(3,I))
		P(4,1)=X(FACETTES(4,I))
		P(4,2)=Y(FACETTES(4,I))
		P(4,3)=Z(FACETTES(4,I))
		IMN=1
		ZMN=P(IMN,3)
		ZMX=P(IMN,3)
		DO J=IMN+1,4
			IF (P(J,3).LT.ZMN) THEN
				ZMN=P(J,3)
				IMN=J
			END IF
			IF (P(J,3).GT.ZMX) THEN
				ZMX=P(J,3)
			END IF
		END DO
		K=IMN
		DO J=1,4			
			IF (K.GT.4) THEN
				K=K-4
			END IF
			P(J,1)=X(FACETTES(K,I))
			P(J,2)=Y(FACETTES(K,I))
			P(J,3)=Z(FACETTES(K,I))
			K=K+1
		END DO
		IF (ZMN.LT.0.0) THEN
			! CAS OU Z2 ET Z4 SONT POSITIFS
			IF ((P(2,3).GT.0.0).AND.(P(4,3).GT.0.0)) THEN
				J=2
				C=P(1,3)/(P(1,3)-P(J,3))
				DO K=1,3
					P(J,K)=P(1,K)+C*(P(J,K)-P(1,K))
				END DO
				J=4
				C=P(1,3)/(P(1,3)-P(J,3))
				DO K=1,3
					P(J,K)=P(1,K)+C*(P(J,K)-P(1,K))
				END DO
				DO K=1,3
					P1(K)=P(1,K)
					P2(K)=P(2,K)
					P3(K)=P(4,K)
					P4(K)=P(1,K)
				END DO
			ELSE 
				IF ((P(2,3).GT.0.0).AND.(P(3,3).GT.0.0)) THEN
					! CAS OU Z2 ET Z3 SONT POSITIFS, Z4 NEG
					J=2
					C=P(1,3)/(P(1,3)-P(J,3))
					DO K=1,3
						P(J,K)=P(1,K)+C*(P(J,K)-P(1,K))
					END DO
					J=3
					C=P(4,3)/(P(4,3)-P(J,3))
					DO K=1,3
						P(J,K)=P(4,K)+C*(P(J,K)-P(4,K))
					END DO
					DO K=1,3
						P1(K)=P(1,K)
						P2(K)=P(2,K)
						P3(K)=P(3,K)
						P4(K)=P(4,K)
					END DO
				ELSE 
					IF ((P(4,3).GT.0.0).AND.(P(3,3).GT.0.0)) THEN
						! CAS OU Z3 ET Z4 SONT POSITIFS, Z2 NEG
						J=4
						C=P(1,3)/(P(1,3)-P(J,3))
						DO K=1,3
							P(J,K)=P(1,K)+C*(P(J,K)-P(1,K))
						END DO
						J=3
						C=P(2,3)/(P(2,3)-P(J,3))
						DO K=1,3
							P(J,K)=P(2,K)+C*(P(J,K)-P(2,K))
						END DO
						DO K=1,3
							P1(K)=P(1,K)
							P2(K)=P(2,K)
							P3(K)=P(3,K)
							P4(K)=P(4,K)
						END DO
					ELSE
						IF (P(3,3).GT.0.0) THEN
						! CAS OU Z3 EST POSITIF, Z2 ET Z4 NEG
							J=3
							C=P(2,3)/(P(2,3)-P(J,3))
							DO K=1,3
								PL(K)=P(2,K)+C*(P(J,K)-P(2,K))
							END DO
							DO K=1,3
								P1(K)=P(1,K)
								P2(K)=P(2,K)
								P3(K)=PL(K)
								P4(K)=P(1,K)
							END DO	
							IF (ITEC.EQ.1) THEN
								WRITE(10,*) (P1(K),K=1,3)
								WRITE(10,*) (P2(K),K=1,3)
								WRITE(10,*) (P3(K),K=1,3)
								WRITE(10,*) (P4(K),K=1,3)
							END IF
							Xmouillee(Npmouillee+1)=P1(1)
							Xmouillee(Npmouillee+2)=P2(1)
							Xmouillee(Npmouillee+3)=P3(1)
							Xmouillee(Npmouillee+4)=P4(1)
							Ymouillee(Npmouillee+1)=P1(2)
							Ymouillee(Npmouillee+2)=P2(2)
							Ymouillee(Npmouillee+3)=P3(2)
							Ymouillee(Npmouillee+4)=P4(2)
							Zmouillee(Npmouillee+1)=P1(3)
							Zmouillee(Npmouillee+2)=P2(3)
							Zmouillee(Npmouillee+3)=P3(3)
							Zmouillee(Npmouillee+4)=P4(3)
							DO k=1,4
								Facettesmouillee(k,Nfmouillee+1)=Npmouillee+k
							END DO
							Npmouillee=4+Npmouillee
							Nfmouillee=1+Nfmouillee	
							CALL VOLELMT(P1,P2,P3,P4,VLF,PGF,SEF,KHE,RHO,G)
							DEPLACEMENT=DEPLACEMENT+VLF
							XF=XF+VLF*PGF(1)
							YF=YF+VLF*PGF(2)
							ZF=ZF+VLF*PGF(3)
							SF=SF+SEF
							KH=KH+KHE
							J=3
							C=P(4,3)/(P(4,3)-P(J,3))
							DO K=1,3
								PK(K)=P(4,K)+C*(P(J,K)-P(4,K))
							END DO
							DO K=1,3
								P1(K)=P(4,K)
								P2(K)=P(1,K)
								P3(K)=PL(K)
								P4(K)=PK(K)
							END DO
						ELSE
							IF (P(4,3).GT.0.0) THEN
								! CAS OU Z4 EST POSITIF, Z2 ET Z3 NEG
								J=4
								C=P(3,3)/(P(3,3)-P(J,3))
								DO K=1,3
									PL(K)=P(3,K)+C*(P(J,K)-P(3,K))
								END DO
								DO K=1,3
									P1(K)=P(2,K)
									P2(K)=P(3,K)
									P3(K)=PL(K)
									P4(K)=P(2,K)
								END DO
								IF (ITEC.EQ.1) THEN
									WRITE(10,*) (P1(K),K=1,3)
									WRITE(10,*) (P2(K),K=1,3)
									WRITE(10,*) (P3(K),K=1,3)
									WRITE(10,*) (P4(K),K=1,3)
								END IF
								Xmouillee(Npmouillee+1)=P1(1)
								Xmouillee(Npmouillee+2)=P2(1)
								Xmouillee(Npmouillee+3)=P3(1)
								Xmouillee(Npmouillee+4)=P4(1)
								Ymouillee(Npmouillee+1)=P1(2)
								Ymouillee(Npmouillee+2)=P2(2)
								Ymouillee(Npmouillee+3)=P3(2)
								Ymouillee(Npmouillee+4)=P4(2)
								Zmouillee(Npmouillee+1)=P1(3)
								Zmouillee(Npmouillee+2)=P2(3)
								Zmouillee(Npmouillee+3)=P3(3)
								Zmouillee(Npmouillee+4)=P4(3)
								DO k=1,4
									Facettesmouillee(k,Nfmouillee+1)=Npmouillee+k
								END DO
								Npmouillee=4+Npmouillee
								Nfmouillee=1+Nfmouillee	
								CALL VOLELMT(P1,P2,P3,P4,VLF,PGF,SEF,KHE,RHO,G)
								DEPLACEMENT=DEPLACEMENT+VLF
								XF=XF+VLF*PGF(1)
								YF=YF+VLF*PGF(2)
								ZF=ZF+VLF*PGF(3)
								SF=SF+SEF
								KH=KH+KHE
								J=4
								C=P(1,3)/(P(1,3)-P(J,3))
								DO K=1,3
									PK(K)=P(1,K)+C*(P(J,K)-P(1,K))
								END DO
								DO K=1,3
									P1(K)=P(1,K)
									P2(K)=P(2,K)
									P3(K)=PL(K)
									P4(K)=PK(K)
								END DO								
							ELSE
								IF (P(4,3).GT.0.0) THEN
									! CAS OU Z4 EST POSITIF, Z2 ET Z3 NEG
									J=4
									C=P(3,3)/(P(3,3)-P(J,3))
									DO K=1,3
										PL(K)=P(3,K)+C*(P(J,K)-P(3,K))
									END DO
									DO K=1,3
										P1(K)=P(2,K)
										P2(K)=P(3,K)
										P3(K)=PL(K)
										P4(K)=P(2,K)
									END DO
									IF (ITEC.EQ.1) THEN
										WRITE(10,*) (P1(K),K=1,3)
										WRITE(10,*) (P2(K),K=1,3)
										WRITE(10,*) (P3(K),K=1,3)
										WRITE(10,*) (P4(K),K=1,3)
									END IF
									Xmouillee(Npmouillee+1)=P1(1)
									Xmouillee(Npmouillee+2)=P2(1)
									Xmouillee(Npmouillee+3)=P3(1)
									Xmouillee(Npmouillee+4)=P4(1)
									Ymouillee(Npmouillee+1)=P1(2)
									Ymouillee(Npmouillee+2)=P2(2)
									Ymouillee(Npmouillee+3)=P3(2)
									Ymouillee(Npmouillee+4)=P4(2)
									Zmouillee(Npmouillee+1)=P1(3)
									Zmouillee(Npmouillee+2)=P2(3)
									Zmouillee(Npmouillee+3)=P3(3)
									Zmouillee(Npmouillee+4)=P4(3)
									DO k=1,4
										Facettesmouillee(k,Nfmouillee+1)=Npmouillee+k
									END DO
									Npmouillee=4+Npmouillee
									Nfmouillee=1+Nfmouillee	
									CALL VOLELMT(P1,P2,P3,P4,VLF,PGF,SEF,KHe,RHO,G)
									DEPLACEMENT=DEPLACEMENT+VLF
									XF=XF+VLF*PGF(1)
									YF=YF+VLF*PGF(2)
									ZF=ZF+VLF*PGF(3)
									SF=SF+SEF
									KH=KH+KHe
									J=4
									C=P(1,3)/(P(1,3)-P(J,3))
									DO K=1,3
										PK(K)=P(1,K)+C*(P(J,K)-P(1,K))
									END DO
									DO K=1,3
										P1(K)=P(1,K)
										P2(K)=P(2,K)
										P3(K)=PL(K)
										P4(K)=PK(K)
									END DO								
								ELSE
									IF (P(2,3).GT.0.0) THEN
										! CAS OU Z2 EST POSITIF, Z3 ET Z4 NEG
										J=2
										C=P(1,3)/(P(1,3)-P(J,3))
										DO K=1,3
											PL(K)=P(1,K)+C*(P(J,K)-P(1,K))
										END DO
										DO K=1,3
											P1(K)=P(4,K)
											P2(K)=P(1,K)
											P3(K)=PL(K)
											P4(K)=P(4,K)
										END DO
										IF (ITEC.EQ.1) THEN
											WRITE(10,*) (P1(K),K=1,3)
											WRITE(10,*) (P2(K),K=1,3)
											WRITE(10,*) (P3(K),K=1,3)
											WRITE(10,*) (P4(K),K=1,3)
										END IF
										Xmouillee(Npmouillee+1)=P1(1)
										Xmouillee(Npmouillee+2)=P2(1)
										Xmouillee(Npmouillee+3)=P3(1)
										Xmouillee(Npmouillee+4)=P4(1)
										Ymouillee(Npmouillee+1)=P1(2)
										Ymouillee(Npmouillee+2)=P2(2)
										Ymouillee(Npmouillee+3)=P3(2)
										Ymouillee(Npmouillee+4)=P4(2)
										Zmouillee(Npmouillee+1)=P1(3)
										Zmouillee(Npmouillee+2)=P2(3)
										Zmouillee(Npmouillee+3)=P3(3)
										Zmouillee(Npmouillee+4)=P4(3)
										DO k=1,4
											Facettesmouillee(k,Nfmouillee+1)=Npmouillee+k
										END DO
										Npmouillee=4+Npmouillee
										Nfmouillee=1+Nfmouillee	
										CALL VOLELMT(P1,P2,P3,P4,VLF,PGF,SEF,KHe,RHO,G)
										DEPLACEMENT=DEPLACEMENT+VLF
										XF=XF+VLF*PGF(1)
										YF=YF+VLF*PGF(2)
										ZF=ZF+VLF*PGF(3)
										SF=SF+SEF
										KH=KH+KHe
										J=2
										C=P(3,3)/(P(3,3)-P(J,3))
										DO K=1,3
											PK(K)=P(3,K)+C*(P(J,K)-P(3,K))
										END DO
										DO K=1,3
											P1(K)=P(3,K)
											P2(K)=P(4,K)
											P3(K)=PL(K)
											P4(K)=PK(K)
										END DO								
									ELSE
										IF ((P(2,3).GT.0.0).AND.(P(4,3).EQ.0.0)) THEN
											J=2
											C=P(1,3)/(P(1,3)-P(J,3))
											DO K=1,3
												PL(K)=P(1,K)+C*(P(J,K)-P(1,K))
											END DO
											DO K=1,3
												P1(K)=P(1,K)
												P2(K)=PL(K)
												P3(K)=P(3,K)
												P4(K)=P(4,K)
											END DO	
										ELSE
											DO K=1,3
												P1(K)=P(1,K)
												P2(K)=P(2,K)
												P3(K)=P(3,K)
												P4(K)=P(4,K)
											END DO
										END IF
									END IF
								END IF
							END IF
						END IF
					END IF					
				END IF
			END IF
			Xmouillee(Npmouillee+1)=P1(1)
			Xmouillee(Npmouillee+2)=P2(1)
			Xmouillee(Npmouillee+3)=P3(1)
			Xmouillee(Npmouillee+4)=P4(1)
			Ymouillee(Npmouillee+1)=P1(2)
			Ymouillee(Npmouillee+2)=P2(2)
			Ymouillee(Npmouillee+3)=P3(2)
			Ymouillee(Npmouillee+4)=P4(2)
			Zmouillee(Npmouillee+1)=P1(3)
			Zmouillee(Npmouillee+2)=P2(3)
			Zmouillee(Npmouillee+3)=P3(3)
			Zmouillee(Npmouillee+4)=P4(3)
			DO k=1,4
				Facettesmouillee(k,Nfmouillee+1)=Npmouillee+k
			END DO
			Npmouillee=4+Npmouillee
			Nfmouillee=1+Nfmouillee
			CALL VOLELMT(P1,P2,P3,P4,VLF,PGF,SEF,KHe,RHO,G)
			IF (ITEC.EQ.1) THEN
				WRITE(10,*) (P1(K),K=1,3)
				WRITE(10,*) (P2(K),K=1,3)
				WRITE(10,*) (P3(K),K=1,3)
				WRITE(10,*) (P4(K),K=1,3)
			END IF
			DEPLACEMENT=DEPLACEMENT+VLF
			XF=XF+VLF*PGF(1)
			YF=YF+VLF*PGF(2)
			ZF=ZF+VLF*PGF(3)
			SF=SF+SEF
			KH=KH+KHe
		END IF
	END DO
	IF (DEPLACEMENT.GT.0.0) THEN
		XF=XF/DEPLACEMENT
		YF=YF/DEPLACEMENT
		ZF=ZF/DEPLACEMENT
	END IF
	IF (ITEC.EQ.1) THEN
		CLOSE(10)
	END IF

	RETURN

	END SUBROUTINE HYDRO

!	CALCUL DU VOLUME COMPRIS ENTRE UNE FACETTE DU MAILLAGE
!	ET LE PLAN Z=0

	SUBROUTINE VOLELMT(P1,P2,P3,P4,VOLUME,PG,SF,KH,RHO,G)

	IMPLICIT NONE

	REAL,DIMENSION(3) :: P1,P2,P3,P4,P01,P02,P03,P04
	REAL,DIMENSION(3) :: PG,P0G,PG1,PG2,PG3,PG4,PG5,PG6
	REAL,DIMENSION(3) :: AB,AC,AD,U,V,W
	REAL,DIMENSION(6,6) :: KH
	REAL VOLUME,V1,V2,V3,V4,V5,V6,SF,RHO,G,S1,S2
	INTEGER I

	DO I=1,2
		P01(I)=P1(I)
		P02(I)=P2(I)
		P03(I)=P3(I)
		P04(I)=P4(I)
	END DO
	I=3
	P01(I)=0.0
	P02(I)=0.0
	P03(I)=0.0
	P04(I)=0.0
!	CALCUL DU VOLUME PAR DECOMPOSITION EN 6 TETRAEDRES
	CALL VOLTETRA(P1,P2,P4,P01,V1,PG1)	
	CALL VOLTETRA(P2,P4,P01,P04,V2,PG2)
	CALL VOLTETRA(P01,P02,P04,P2,V3,PG3)
	CALL VOLTETRA(P2,P3,P4,P04,V4,PG4)
	CALL VOLTETRA(P02,P03,P04,P2,V5,PG5)
	CALL VOLTETRA(P2,P3,P04,P03,V6,PG6)
	VOLUME=V1+V2+V3+V4+V5+V6
	DO I=1,3
		PG(I)=V1*PG1(I)+V2*PG2(I)+V3*PG3(I)+V4*PG4(I)+V5*PG5(I)+V6*PG6(I)
		IF (ABS(VOLUME).GT.0.0) THEN
			PG(I)=PG(I)/VOLUME
		ELSE
			PG(I)=0.0
		END IF
	END DO
!	CALCUL DE LA SURFACE DE FLOTTAISON
	AB=P02-P01
	AC=P03-P01
	AD=P04-P01
	CALL PRDVCT(AB,AC,U)
	U=0.5*U
	S1=SQRT(U(1)*U(1)+U(2)*U(2)+U(3)*U(3))
	CALL PRDVCT(AC,AD,V)
	V=0.5*V
	S2=SQRT(V(1)*V(1)+V(2)*V(2)+V(3)*V(3))
	W=(U+V)
	SF=SQRT(W(1)*W(1)+W(2)*W(2)+W(3)*W(3))
!	SIGNE DU VOLUME ET DE LA SURFACE EN FONCTION DE L ORIENTATION DE LA FACETTE
	AB=P2-P1
	AC=P3-P1
	AD=P4-P1
	CALL PRDVCT(AB,AC,U)
	CALL PRDVCT(AC,AD,V)
	W=U+V
	IF (W(3).GT.0.0) THEN
		VOLUME=-1.0*VOLUME
		S1=-1.*S1
		S2=-1.*S2
		SF=-1.0*SF
	END IF
!	Calcul des contributions elementaires a la matrice de raideur hydrostatique
    IF (SF.EQ.0) THEN 
        P0G=0.
    ELSE
        P0G=1./3.*(S1*(P01+P02+P03)+S2*(P01+P03+P04))/SF
    END IF
!    IF (SF.EQ.0) THEN 
!        WRITE(*,'(A,E14.7,A)') 'WARNING: At least one panel with 0 area has been found!'
!    END IF
    
	! Ajout d'un test pour prendre en compte les facettes triangulaires	
!     IF ( P01(1)==P02(1) ) THEN
!	     IF (P01(2)==P02(2)) THEN
!		 	IF (P01(3)==P02(3)) THEN
 !       	   P0G=(1./3.)*(P02+P03+P04)
!		 	ELSE
!			   P0G=0.25*(P01+P02+P03+P04)
!			END IF
!		 ELSE
!		    P0G=0.25*(P01+P02+P03+P04)
!		 END IF
!		 
 !     ELSEIF ( P02(1)==P03(1) ) THEN
!	     IF ( P02(2)==P03(2) ) THEN
!		 	IF ( P02(3)==P03(3) ) THEN     
!			   P0G=(1./3.)*(P01+P03+P04)
!		 	ELSE
!			   P0G=0.25*(P01+P02+P03+P04)
!		    END IF
!	     ELSE
!	        P0G=0.25*(P01+P02+P03+P04)
!		 END IF 
!		 	  
 !     ELSEIF ( P03(1)==P04(1) ) THEN
!	     IF ( P03(2)==P04(2) ) THEN
!		 	IF ( P03(3)==P04(3) ) THEN             	  
!			   P0G=(1./3.)*(P01+P02+P04)
!			ELSE
!			   P0G=0.25*(P01+P02+P03+P04)
!		 	END IF
!		 ELSE
!			P0G=0.25*(P01+P02+P03+P04)
!		 END IF 
!		 	  
 !     ELSEIF ( P04(1)==P01(1) ) THEN
!	     IF (P04(2)==P01(2)) THEN
!		 	IF (P04(3)==P01(3)) THEN   
!			   P0G=(1./3.)*(P01+P02+P03)
!			ELSE
!			   P0G=0.25*(P01+P02+P03+P04)
!		 	END IF
!		 ELSE
!		    P0G=0.25*(P01+P02+P03+P04)
!		 END IF 
		 	  
!      ELSE
!    	 P0G=0.25*(P01+P02+P03+P04)
!      END IF
	! 
	
	KH=0.
	KH(3,3)=RHO*G*SF
	KH(4,4)=RHO*G*SF*P0G(2)**2
	KH(5,5)=RHO*G*SF*P0G(1)**2
	KH(3,4)=RHO*G*SF*P0G(2)
	KH(3,5)=-RHO*G*SF*P0G(1)
	KH(4,5)=-RHO*G*SF*P0G(2)*P0G(1)
	KH(4,3)=KH(3,4)
	KH(5,3)=KH(3,5)
	KH(4,5)=KH(5,4)
	
	RETURN

	END SUBROUTINE VOLELMT

!	CALCUL DU VOLUME ET DU CENTRE DE GRAVITE D UN TETRAEDRE.

	SUBROUTINE VOLTETRA(PA,PB,PC,PD,VOLUME,PG)

	IMPLICIT NONE

	REAL,DIMENSION(3) :: PA,PB,PC,PD,PG,PK
	REAL,DIMENSION(3) :: AB,AC,AD,W
	REAL :: VOLUME
	REAL,PARAMETER :: S3=1.0/3.0
	INTEGER I

	DO I=1,3
		AB(I)=PB(I)-PA(I)
		AC(I)=PC(I)-PA(I)
		AD(I)=PD(I)-PA(I)
		PK(I)=PA(I)+s3*(AC(I)+AB(I))
		PG(I)=PD(I)+3.0*0.25*(PK(I)-PD(I))
	END DO
	CALL PRDVCT(AB,AC,W)
	VOLUME=0.5*s3*ABS((AD(1)*W(1)+AD(2)*W(2)+AD(3)*W(3)))
!	VOLUME=0.5*s3*(AD(1)*W(1)+AD(2)*W(2)+AD(3)*W(3))
	
	RETURN

	END SUBROUTINE VOLTETRA

!	PRODUIT VECTORIEL W=UxV

	SUBROUTINE PRDVCT(U,V,W)

	IMPLICIT NONE
	
	REAL,DIMENSION(3) :: U,V,W
	
	W(1)=U(2)*V(3)-U(3)*V(2)
	W(2)=U(3)*V(1)-U(1)*V(3)
	W(3)=U(1)*V(2)-U(2)*V(1)
	
	END
