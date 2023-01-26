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
!   - J. Singh
!   - P. Guével
!   - J.C. Daubisse 
!
!--------------------------------------------------------------------------------------
MODULE SOLVE_BEM_FD_DIRECT

  USE COM_VAR
  USE COMPUTE_GREEN_FD
  USE ELEMENTARY_FNS
  USE M_SOLVER

  CONTAINS
!--------------------------------------------------------------------------!
  SUBROUTINE SOLVE_POTENTIAL_FD_DIRECT(NVEL,AMH,NEXP)
!In this subroutine the linear system Ax=b is constructed

!    LOGICAL:: RHS
    COMPLEX,DIMENSION(*) :: NVEL
    INTEGER:: BX,ISYM,NJJ,N
    INTEGER:: I,J,ISP,IFP
    INTEGER:: NP1,I1,JJ,NEXP
    REAL:: W,BETA,BJ,GM,DIJ,AKK 
    REAL:: CB,SB,CP,CR,COEFB,COEFS,PCOS,PSIN
    REAL:: AKAD,AM0,AKH,AMH,SD1B,SD1S,SD2B,SD2S,PI,DPI,ZERO
    COMPLEX:: ZOL(IMX,2),B(IMX)
    


      NJJ=NSYMY+1
      PI=4.*ATAN(1.)
      DPI=2.*PI
      W=DPI/T
      ZERO=0.0
      AKH=W**2*Depth/G                                                 
!       AMH=X0(AKH)
      AKH=AMH*TANH(AMH)                                                                                                    
      AM0=AMH/Depth 


      IF(ABS(W).LT.1.E-4)THEN
	WRITE(*,*)'ABS(WR)  = ',ABS(W),' < 1.E-4'
	STOP     
      ENDIF
      IF(AMH-AKH-1.E-03)7106,7106,7101                                          
 7106 CONTINUE !WRITE(*,7102)                                                            
 7102 FORMAT(/5X,'PROFONDEUR QUASI-INFINIE'/5X, &                                
     & 'LE PROGRAMME EN PROFONDEUR INFINIE SERAIT PLUS ADAPTE')                  
      GOTO 7104                                                                 
 7101 IF(AKH-0.1)7103,7103,7104                                                 
 7103 CONTINUE !WRITE(*,7105)                                                            
 7105 FORMAT(/5X,'PROFONDEUR TROP FAIBLE POUR LA LONGUEUR D''ONDE')             
 7104 CONTINUE 
!--------------Initilizations---------------
      CQ=0.0
      QQ=0.0
      AMBDA=0.0
      AR=0.0
      VSXP=0.
      VSXM=0.
      VSYP=0.
      VSYM=0.
      VSZP=0.
      VSZM=0.
      ZOL=CMPLX(0.,0.)
      ZIJ=CMPLX(0.,0.)
!--------------------------------------------                                                                                                                          
    GM=0.
    NP1=NP-1
    DO 7100 I=1,NP1
	    I1=I+1
	    DO 7200 JJ=1,NJJ
	        BJ=(-1.)**(JJ+1)
	        DO 7200 J=I1,NP
	            DIJ=SQRT((X(J)-X(I))**2+(Y(I)-Y(J)*BJ)**2)
	            GM=AMAX1(DIJ,GM)
        7200 CONTINUE
    7100 CONTINUE                                              
    AKK=AM0*GM     
    CALL CINT_FD(AKK,N)
    NQ=N
    CALL LISC(AKH,AMH,NEXP)

!   Construction of the influence matrix 
    IF (w.NE.w_previous) THEN
        w_previous=w
        DO ISYM=1,NJJ
            BX=(-1)**(ISYM+1)
	        DO 1 ISP=1,IMX
    	        IF(ISYM.EQ.1)THEN
    	            DO 31 IFP=1,IMX
    	                call VAVFD(2,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP)  !1/r+1/r1           
    	                call VNSFD(AM0,AMH,NEXP,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP))  
    	                PCOS=VSXP1*XN(ISP)+VSYP1*YN(ISP)+VSZP1*ZN(ISP)
    	                PSIN=VSXP2*XN(ISP)+VSYP2*YN(ISP)+VSZP2*ZN(ISP)
    	                ZIJ(ISP,IFP)=CMPLX(PCOS,PSIN)
    	        31 CONTINUE
           	    ELSE
    	            DO 32 IFP=1,IMX
    	                call VAVFD(2,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP)
                       call VNSFD(AM0,AMH,NEXP,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP))  
    	                PCOS=VSXM1*XN(ISP)+VSYM1*YN(ISP)+VSZM1*ZN(ISP)
    	                PSIN=VSXM2*XN(ISP)+VSYM2*YN(ISP)+VSZM2*ZN(ISP)
    	                ZIJ(ISP,IFP)=CMPLX(PCOS,PSIN)           
	            32 CONTINUE
    	        ENDIF
           1 CONTINUE
    	    DO I=1,IMX
    	        DO J=1,IMX
    	            ZIJ(I,J+IMX)=CMPLX(0.,0.)
    	        END DO
    	        ZIJ(I,I+IMX)=CMPLX(1.,0.)
    	    END DO
!------------------------------------------------!
           CALL GAUSSZ(ZIJ,NFA,IMX,2*IMX)
!------------------------------------------------!
    	    DO I=1,IMX
    	        DO J=1,IMX
    	            AInv(I,J,(ISYM-1)+1)=ZIJ(I,J+IMX)
    	        END DO
    	    END DO
        END DO
    END IF
    
!   Solution of the linear problem A*ZOL=B and calculation of source distribution (ZI)
    DO ISYM=1,NJJ
        BX=(-1)**(ISYM+1)
        DO I=1,IMX
    	    IF (NSYMY.EQ.1) THEN
    	        B(I)=(NVEL(I)+BX*NVEL(I+NFA))*0.5
    	    ELSE
    	        B(I)=NVEL(I)
    	    END IF
    	END DO
    	DO I=1,IMX
    	    ZOL(I,(ISYM-1)+1)=(0.,0.)
    	    DO K=1,IMX
    	        ZOL(I,(ISYM-1)+1)=ZOL(I,(ISYM-1)+1)+AInv(I,K,(ISYM-1)+1)*B(K)    
    	    END DO
    	END DO
    END DO
    ZIGB=(0.,0.)
	ZIGS=(0.,0.)
	DO I=1,IMX
	    IF(NSYMY .EQ. 0)THEN
	        ZIGB(I)=ZOL(I,1)    ! Factor 2 is removed in comparison with previous version of Nemoh because
	                            ! Normal velocity is halved in Nemoh (because of symmetry)
	        ZIGS(I)=0.0
	    ELSE
	        ZIGB(I)=(ZOL(I,1)+ZOL(I,2))
	        ZIGS(I)=(ZOL(I,1)-ZOL(I,2))
        ENDIF
    END DO  

!   Computation of potential phi=S*sigma on the boundary
    ZPB=(0.,0.)
    ZPS=(0.,0.)
    DO I=1,IMX
	    DO J=1,IMX
	        call VAVFD(2,XG(I),YG(I),ZG(I),I,J)
            call VNSFD(AM0,AMH,NEXP,I,J,XG(I),YG(I),ZG(I))  
	        ZPB(I)=ZPB(I)+0.5*(ZIGB(J)*CMPLX(SP1+SM1,SP2+SM2)+ZIGS(J)*CMPLX(SP1-SM1,SP2-SM2))
	        ZPS(I)=ZPS(I)+0.5*(ZIGS(J)*CMPLX(SP1+SM1,SP2+SM2)+ZIGB(J)*CMPLX(SP1-SM1,SP2-SM2))
        END DO
    END DO
    
END SUBROUTINE
!-------------------------------------------------
END MODULE