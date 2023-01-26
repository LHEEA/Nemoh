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
!   - G. Delhommeau
!   - P. Guével
!   - J.C. Daubisse
!   - J. Singh 
!   - A. Babarit 
!
!--------------------------------------------------------------------------------------
    SUBROUTINE COMPUTE_KOCHIN(kwave,Theta,HKochin)
!
    USE COM_VAR
    USE ELEMENTARY_FNS
!    
    IMPLICIT NONE
!
!   Inputs / Outputs
    REAL :: kwave
    REAL :: THeta
    COMPLEX :: HKochin,HKochin1    
!   Locals
    COMPLEX,PARAMETER :: II=CMPLX(0.,1.)
    REAL :: PI
    INTEGER :: i,k,l
    REAL :: wbar
    REAL :: CTE,CTM
    COMPLEX,DIMENSION(NFA*2**NSYMY) :: ZS
    COMPLEX,DIMENSION(NP) :: CEP,CEM,ZJ
    INTEGER,DIMENSION(5) :: KK
    COMPLEX :: ZMS,ZNS,ZCS,ZC,ZHS,ZHD
!
    PI=4.*ATAN(1.)
!   Compute Kochin coefficients (integration using Gauss - 1 point)
    HKochin1=0.
!    DO i=1,NFA
!        wbar=(XG(i)-XEFF)*COS(Theta)+(YG(i)-YEFF)*SIN(Theta)
!        HKochin1=HKochin1+ZIGB(i)*CIH(kwave,ZG(i),Depth)*CEXP(-II*kwave*wbar)*AIRE(i)
!        IF (NSYMY.EQ.1) THEN
!            wbar=(XG(i)-XEFF)*COS(Theta)+(-YG(i)-YEFF)*SIN(Theta)
!            HKochin1=HKochin1+ZIGS(i)*CIH(kwave,ZG(i),Depth)*CEXP(-II*kwave*wbar)*AIRE(i)
!        END IF
!    END DO  
!   Compute Kochin coefficients (analytical integration)   
    DO k=1,NP
        wbar=(X(k)-XEFF)*COS(Theta)+(Y(k)-YEFF)*SIN(Theta)
        ZJ(k)=kwave*(Z(k)-II*wbar)
        IF (REAL(ZJ(k)).LT.-18.) THEN
            CEP(k)=(0.,0.)
        ELSE
            CEP(k)=CEXP(ZJ(k))
        END IF
        IF ((REAL(-ZJ(k)-2.*kwave*Depth).LT.-18.).OR.(-kwave*Depth.GE.0.)) THEN
            CEM(k)=(0.,0.)
        ELSE
            CEM(k)=CEXP(-ZJ(k)-2.*kwave*Depth)
        END IF
    END DO    
    DO i=1,NFA
        KK(1)=M1(i)
		KK(2)=M2(i)
		KK(3)=M3(i)
		KK(4)=M4(i)
		KK(5)=KK(1)
		ZMS=YN(i)-II*ZN(i)*SIN(Theta)
		ZNS=XN(i)-II*ZN(i)*COS(Theta)
		ZS(i)=(0.,0.)
		DO L=1,4
	        ZCS=ZMS*(X(KK(l+1))-X(KK(l)))-ZNS*(Y(KK(l+1))-Y(KK(l)))		
			ZC=ZJ(KK(L+1))-ZJ(KK(L))
			IF(ABS(AIMAG(ZC)).LT.1.E-04.AND.ABS(REAL(ZC)).LT.1.E-04)THEN
				ZHS=0.5*(CEP(KK(L+1))+CEP(KK(L)))
				ZHD=-0.5*(CEM(KK(L+1))+CEM(KK(L)))
			ELSE
				ZHS=(CEP(KK(L+1))-CEP(KK(L)))/ZC
				ZHD=(CEM(KK(L+1))-CEM(KK(L)))/ZC
			ENDIF
			ZS(i)=ZS(i)+ZCS*ZHS+CONJG(ZCS*ZHD)
	    END DO
	    IF ((kwave*Depth.GT.18.).OR.(kwave*Depth.LE.0.)) THEN
	        ZS(i)=ZS(i)/kwave    
	    ELSE
	        ZS(i)=0.5*ZS(i)*EXP(kwave*Depth)/SINH(kwave*Depth)*1./kwave
	    END IF
	 END DO
     IF (NSYMY.EQ.1) THEN
         DO k=1,NP
            wbar=(X(k)-XEFF)*COS(Theta)+(-Y(k)-YEFF)*SIN(Theta)
            ZJ(k)=kwave*(Z(k)-II*wbar)
            IF (REAL(ZJ(k)).LT.-18.) THEN
                CEP(k)=(0.,0.)
            ELSE
                CEP(k)=CEXP(ZJ(k))
            END IF
            IF ((REAL(-ZJ(k)-2.*kwave*Depth).LT.-18.).OR.(-kwave*Depth.GE.0.)) THEN
                CEM(k)=(0.,0.)
            ELSE
                CEM(k)=CEXP(-ZJ(k)-2.*kwave*Depth)
            END IF
        END DO    
        DO i=1,NFA
            KK(1)=M1(i)
		    KK(2)=M2(i)
		    KK(3)=M3(i)
		    KK(4)=M4(i)
		    KK(5)=KK(1)
		    ZMS=-YN(i)-II*ZN(i)*SIN(Theta)
		    ZNS=XN(i)-II*ZN(i)*COS(Theta)
    !		Initialisation de l integrale sur la facette J
		    ZS(i+NFA)=(0.,0.)
		    DO L=1,4
	            ZCS=ZMS*(X(KK(l+1))-X(KK(l)))+ZNS*(Y(KK(l+1))-Y(KK(l)))		
			    ZC=ZJ(KK(L+1))-ZJ(KK(L))
			    IF(ABS(AIMAG(ZC)).LT.1.E-04.AND.ABS(REAL(ZC)).LT.1.E-04)THEN
				    ZHS=0.5*(CEP(KK(L+1))+CEP(KK(L)))
				    ZHD=-0.5*(CEM(KK(L+1))+CEM(KK(L)))
			    ELSE
				    ZHS=(CEP(KK(L+1))-CEP(KK(L)))/ZC
				    ZHD=(CEM(KK(L+1))-CEM(KK(L)))/ZC
			    ENDIF
			    ZS(i+NFA)=ZS(i+NFA)+ZCS*ZHS+CONJG(ZCS*ZHD)
	        END DO
	        IF ((kwave*Depth.GT.18.).OR.(kwave*Depth.LE.0.)) THEN
	            ZS(i+NFA)=-ZS(i+NFA)/kwave    
	        ELSE
	            ZS(i+NFA)=-0.5*ZS(i+NFA)*EXP(kwave*Depth)/SINH(kwave*Depth)*1./kwave
	        END IF
	     END DO
     END IF
     HKochin=0.
     DO i=1,NFA
        HKochin=HKochin+ZIGB(i)*ZS(i)
        IF (NSYMY.EQ.1) THEN
            HKochin=HKochin+ZIGS(i)*ZS(i+NFA)
        END IF 
     END DO
     HKochin=0.25/PI*HKochin
!  
    END SUBROUTINE COMPUTE_KOCHIN