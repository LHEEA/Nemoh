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
!
!--------------------------------------------------------------------------------------
MODULE PREPARE_MESH

USE COM_VAR
USE FIC_COM

IMPLICIT NONE

CONTAINS
!-------------------------------------------------------------------------------------!
    SUBROUTINE PRE_PROC_MESH(Mesh)
    
    USE MMesh
     
    IMPLICIT NONE

    INTEGER:: I,J,K,L,NC1 
    REAL :: GL,BL 
    REAL,DIMENSION(3) :: N13,N24  
    REAL :: DPOINT 
    REAL :: XL,YL,ZL
    TYPE(TMesh) :: Mesh                    

!-------------------------------------------------------------------------------------!
!   Copy mesh in Nemoh format
    DO i=1,Mesh%Npoints
        X(i)=Mesh%X(1,i)
        Y(i)=Mesh%X(2,i)
        Z(i)=Mesh%X(3,i)
    END DO
    DO i=1,Mesh%Npanels
        M1(i)=Mesh%P(1,i)
        M2(i)=Mesh%P(2,i)
        M3(i)=Mesh%P(3,i)
        M4(i)=Mesh%P(4,i)
    END DO	
!-------------------------------------------------------------------------------------!
    IMX=NFA
    GL=0.
    BL=0.
    DO I=1,NP
        DO J=I+1,NP                                                          
	        GL=AMAX1(ABS(X(J)-X(I)),GL)     
	        IF(NSYMY.EQ.0) BL=AMAX1(ABS(Y(J)-Y(I)),BL)             
	        IF(NSYMY.EQ.1) BL=AMAX1(ABS(Y(J)),ABS(Y(I)),BL)  
	    END DO
	END DO                          
    IF(NSYMY.EQ.1)BL=2.*BL                                                   
    GL=AMAX1(GL,BL)
    ZER=-1.E-5*GL
!-------------------------------------------------------------------------------------!
    IXX=IMX
!-------------------------------------------------------------------------------------!
    CALL GOMG

    RETURN
    
    END SUBROUTINE
!----------------------------------------------------------------------------

    SUBROUTINE GOMG

    USE COM_VAR
    USE FIC_COM
    IMPLICIT NONE

    INTEGER :: I,J,K,L,M,N,IJJ,NJJ,IMETH
    REAL :: XAVER,YAVER,ZAVER
    REAL:: EPPM,D1,T1,T2,T1X,T1Y,T1Z,T2X,T2Y,T2Z,XNQ,XNX,XNY,XNZ
    REAL:: ETA1,ETA2,ETA3,ETA4,ETAN1,ETAN2,ETAN3,ETAN4,ETAO
    REAL:: T1UNX,T1UNY,T1UNZ,T2UNX,T2UNY,T2UNZ
    REAL:: PSCA,XI1,XI2,XI3,XI4,XIN1,XIN2,XIN3,XIN4,XIO
    REAL:: DPOINT,TI1,TT1,TT2,TT3,TT4
    REAL:: XL1,XL2,XL3,XL4,X41,X32,X34,X21,Y41,Y32,Y34,Y21
    REAL:: A,B,C,D,AA,BB,CC,DD,XG1,YG1
    REAL:: YL1,YL2,YL3,YL4
    Real:: at1,at2,at3,at4,at5,at6,at7,at8,at9
    INTEGER :: NSYM
    LOGICAL FICHEX
!~ C
!~ C     PREPARATION DU CALCUL DES COEFFICIENTS D'INFLUENCE:
!~ C     DETERMINATION DE LA POSITION DES POINTS DE GAUSS
!~ C---------------------GAUSS 1 4 9 OU 16--------------------------------
!~ C     ET DU JACOBIEN POUR CHAQUE FACETTE
!~ C      NG =1,4,9 OU 16
      REAL :: XX(16,4),YY(16,4),WG(16,4)

      DATA((XX(I,L),I=1,16),L=1,4) /&
       16*0., &
!     1  -0.5,15*0.,
       .57735027, .57735027,-.57735027,-.57735027,12*0.,&
       .77459667, .77459667, .77459667,3*0.,-.77459667,&
      -.77459667,-.77459667,7*0.,&
       .86113631, .86113631, .86113631, .86113631,&
       .33998104, .33998104, .33998104, .33998104,&
      -.33998104,-.33998104,-.33998104,-.33998104,&
      -.86113631,-.86113631,-.86113631,-.86113631/
      DATA((YY(I,L),I=1,16),L=1,4)/ &
       16*0., &
!     1  -0.5,15*0.,
      -.57735027,.57735027,-.57735027,.57735027,12*0.,&
      -.77459667,0.,.77459667,-.77459667,0.,.77459667,&
      -.77459667,0.,.77459667,7*0.,&
!     3 -.75459667,0.2,.79459667,-.75459667,0.2,.79459667,
!     3 -.75459667,0.2,.79459667,7*0.,
      -.86113363,-.33998104,.33998104,.86113363,&
      -.86113363,-.33998104,.33998104,.86113363,&
      -.86113363,-.33998104,.33998104,.86113363,&
      -.86113363,-.33998104,.33998104,.86113363/
      DATA((WG(I,L),I=1,16),L=1,4)/ &
       1,15*0., &
      .25,.25,.25,.25,12*0.,&
      .07716049,.12345679,.07716049,.12345679,.19753086,&
      .12345679,.07716049,.12345679,.07716049,7*0.,&
      .30250748E-1,.56712963E-1,.56712963E-1,.30250748E-1,&
      .56712963E-1,.10632333,.10632333,.56712963E-1,&
      .56712963E-1,.10632333,.10632333,.56712963E-1,&
      .30250748E-1,.56712963E-1,.56712963E-1,.30250748E-1/

      NG=4

      IF(NG.EQ.1)THEN
		IMETH=1
      ELSE
      IF(NG.EQ.4)THEN
		IMETH=2
      ELSE
      IF(NG.EQ.9) THEN
		IMETH=3
      ELSE
      IF(NG.EQ.16) THEN
		IMETH=4
      ELSE
		WRITE(*,*)'ARGUMENT INCORRECT DANS SOMGO'
      STOP
      ENDIF
      ENDIF
      ENDIF
      ENDIF

    NSYM=NSYMY+1
    EPPM=1.E-25
    DO I=1,IMX                                                         
	    K=M1(I)
	    L=M2(I)                                                                   
	    M=M3(I)                                                                   
	    N=M4(I)
	    T1X=X(M)-X(K)                                                             
	    T1Y=Y(M)-Y(K)
	    T1Z=Z(M)-Z(K)
	    T2X=X(N)-X(L)                                                             
	    T2Y=Y(N)-Y(L)                                                             
	    T2Z=Z(N)-Z(L)                                                             
	    XNX=T1Y*T2Z-T2Y*T1Z
	    XNY=T1Z*T2X-T2Z*T1X                                                       !
	    XNZ=T1X*T2Y-T2X*T1Y
	    XNQ=SQRT(XNX**2+XNY**2+XNZ**2)
	    IF(XNQ.LE.EPPM)THEN                                                         !
	        WRITE(*,*) 
            WRITE(*,'(A,I7,A)') 'Error: area of panel ',i,' is too small'
            WRITE(*,'(A,E14.7)') '       - Area = ',0.5*XNQ
            WRITE(*,'(A,E14.7)') '       - XG   = ',0.25*(X(M1(I))+X(M2(I))+X(M3(I))+X(M4(I)))
            WRITE(*,'(A,E14.7)') '       - YG   = ',0.25*(Y(M1(I))+Y(M2(I))+Y(M3(I))+Y(M4(I)))
            WRITE(*,'(A,E14.7)') '       - ZG   = ',0.25*(Z(M1(I))+Z(M2(I))+Z(M3(I))+Z(M4(I)))
            STOP
        ENDIF
        XN(I)=XNX/XNQ
        YN(I)=XNY/XNQ
        ZN(I)=XNZ/XNQ
        XAVER=(X(K)+X(M)+X(L)+X(N))*0.25
        YAVER=(Y(K)+Y(M)+Y(N)+Y(L))*0.25
        ZAVER=(Z(K)+Z(M)+Z(N)+Z(L))*0.25
        D1=ABS(XN(I)*(XAVER-X(K))+YN(I)*(YAVER-Y(K))+ZN(I)*(ZAVER-Z(K)))
        T1=SQRT(T1X**2+T1Y**2+T1Z**2)
      T1UNX=T1X/T1
      T1UNY=T1Y/T1
      T1UNZ=T1Z/T1
      T2UNX=YN(I)*T1UNZ-ZN(I)*T1UNY
      T2UNY=ZN(I)*T1UNX-XN(I)*T1UNZ                                             
      T2UNZ=XN(I)*T1UNY-YN(I)*T1UNX
      AT3=-XN(I)
      AT6=-YN(I)
      AT9=-ZN(I)
      AT1=T1UNX
      AT4=T1UNY
      AT7=T1UNZ
      AT2=-T2UNX
      AT5=-T2UNY
      AT8=-T2UNZ
        PSCA=-YN(I)*YAVER
        IF(PSCA.GT.0.)THEN  
            WRITE(*,*)
            WRITE(*,'(A,I7,A)') 'Warning: normal vector of panel ',i,' points towards the x axis'                                                     
        ENDIF                                                                                         
      XI1=T1UNX*(X(K)-XAVER)+T1UNY*(Y(K)-YAVER)+T1UNZ*(Z(K)-ZAVER)
      XI2=T1UNX*(X(L)-XAVER)+T1UNY*(Y(L)-YAVER)+T1UNZ*(Z(L)-ZAVER)
      XI3=T1UNX*(X(M)-XAVER)+T1UNY*(Y(M)-YAVER)+T1UNZ*(Z(M)-ZAVER)
      XI4=T1UNX*(X(N)-XAVER)+T1UNY*(Y(N)-YAVER)+T1UNZ*(Z(N)-ZAVER)
      ETA1=T2UNX*(X(K)-XAVER)+T2UNY*(Y(K)-YAVER)+T2UNZ*(Z(K)-ZAVER)
      ETA2=T2UNX*(X(L)-XAVER)+T2UNY*(Y(L)-YAVER)+T2UNZ*(Z(L)-ZAVER)
      ETA3=T2UNX*(X(M)-XAVER)+T2UNY*(Y(M)-YAVER)+T2UNZ*(Z(M)-ZAVER)             
      ETA4=T2UNX*(X(N)-XAVER)+T2UNY*(Y(N)-YAVER)+T2UNZ*(Z(N)-ZAVER)             
      XIO=(XI4*(ETA1-ETA2)+XI2*(ETA4-ETA1))/(3.*(ETA2-ETA4))
      ETAO=-ETA1*0.33333333                                                  
        ETAN1=ETA1-ETAO
        ETAN2=ETA2-ETAO
        ETAN3=ETA3-ETAO
        ETAN4=ETA4-ETAO                                                           
        XIN1=XI1-XIO
        XIN2=XI2-XIO                                                              
        XIN3=XI3-XIO                                                              
        XIN4=XI4-XIO                                                              
        XG(I)=XAVER+T1UNX*XIO+T2UNX*ETAO
        YG(I)=YAVER+T1UNY*XIO+T2UNY*ETAO                                          
        ZG(I)=ZAVER+T1UNZ*XIO+T2UNZ*ETAO                                         
        DIST(I)=0.                                                                
        DO IJJ=1,NSYM
            NJJ=(-1)**(IJJ+1)                                                         
            DO K=1,NP                                                            
                DPOINT=AMAX1(ABS(XG(I)-X(K)),ABS(YG(I)-Y(K)*NJJ))                         
                DIST(I)=AMAX1(DIST(I),DPOINT)                                             
            END DO
        END DO                                                                
       T1=(XIN3-XIN1)*(XIN3-XIN1)                                                
       T2=(XIN4-XIN2)*(XIN4-XIN2)+(ETAN4-ETAN2)*(ETAN4-ETAN2)                    
       TI1=SQRT(MAX(T1,T2))                                                     
       AIRE(I)=(XIN3-XIN1)*(ETAN4-ETAN2)*0.5                                        
       TT1=SQRT(XIN1*XIN1+ETAN1*ETAN1)                                         
       TT2=SQRT(XIN2*XIN2+ETAN2*ETAN2)                                         
       TT3=SQRT(XIN3*XIN3+ETAN3*ETAN3)                                         
       TT4=SQRT(XIN4*XIN4+ETAN4*ETAN4)                                         
       TDIS(I)=MAX(TT1,TT2,TT3,TT4)
      XL1=AT1*(X(M1(I))-XG(I))+AT4*(Y(M1(I))-YG(I))+AT7*(Z(M1(I))-ZG(I))
      YL1=AT2*(X(M1(I))-XG(I))+AT5*(Y(M1(I))-YG(I))+AT8*(Z(M1(I))-ZG(I))
      XL2=AT1*(X(M2(I))-XG(I))+AT4*(Y(M2(I))-YG(I))+AT7*(Z(M2(I))-ZG(I))
      YL2=AT2*(X(M2(I))-XG(I))+AT5*(Y(M2(I))-YG(I))+AT8*(Z(M2(I))-ZG(I))
      XL3=AT1*(X(M3(I))-XG(I))+AT4*(Y(M3(I))-YG(I))+AT7*(Z(M3(I))-ZG(I))
      YL3=AT2*(X(M3(I))-XG(I))+AT5*(Y(M3(I))-YG(I))+AT8*(Z(M3(I))-ZG(I))
      XL4=AT1*(X(M4(I))-XG(I))+AT4*(Y(M4(I))-YG(I))+AT7*(Z(M4(I))-ZG(I))
      YL4=AT2*(X(M4(I))-XG(I))+AT5*(Y(M4(I))-YG(I))+AT8*(Z(M4(I))-ZG(I))
      X41=XL4-XL1
      X32=XL3-XL2
      Y34=YL3-YL4
      Y21=YL2-YL1
      Y41=YL4-YL1
      Y32=YL3-YL2
      X34=XL3-XL4
      X21=XL2-XL1
      DO L=1,NG
      A=(1.-YY(L,IMETH))*X41+(1+YY(L,IMETH))*X32
      B=(1.-XX(L,IMETH))*Y21+(1+XX(L,IMETH))*Y34
      C=(1.-XX(L,IMETH))*X21+(1+XX(L,IMETH))*X34
      D=(1.-YY(L,IMETH))*Y41+(1+YY(L,IMETH))*Y32
      XJAC(L,I)=(ABS(A*B-C*D)*WG(L,IMETH)*.25)/AIRE(I)
      ENDDO
!~ C
!~ C    COORDONNES ABSOLUES DES POINTS DE GAUSS
!~ C
      DO L=1,NG
      AA=.25*(1-XX(L,IMETH))*(1-YY(L,IMETH))
      BB=.25*(1-XX(L,IMETH))*(1+YY(L,IMETH))
      CC=.25*(1+XX(L,IMETH))*(1+YY(L,IMETH))
      DD=.25*(1+XX(L,IMETH))*(1-YY(L,IMETH))
      XG1=AA*XL1+BB*XL2+CC*XL3+DD*XL4
      YG1=AA*YL1+BB*YL2+CC*YL3+DD*YL4
      XGA(L,I)=XG(I)+AT1*XG1+AT2*YG1 
      YGA(L,I)=YG(I)+AT4*XG1+AT5*YG1 
      ZGA(L,I)=ZG(I)+AT7*XG1+AT8*YG1
      ENDDO
!      IF(NG.EQ.1)THEN
!		  XGA(1,I)=XG(I)
!		  YGA(1,I)=YG(I)
!		  ZGA(1,I)=ZG(I)
!      ENDIF
    END DO
      PRINT *,'NG = ',NG
!~ C
!~ C    CREATION OU LECTURE DU FICHIER DES FONCTIONs DE GREEN
!~ C
    INQUIRE(FILE='GRIN.QAT',EXIST=FICHEX)
       print *,'grin.qat ',fichex
    IF(.NOT.FICHEX)THEN
    CALL CREK
     OPEN(UNIT=44,FILE='GRIN.QAT',FORM='UNFORMATTED',STATUS='NEW')
      WRITE(44)IIR,JJZ,(Xr(I),I=1,IIR),(xZ(J),J=1,JJZ)
      DO J=1,JJZ
      WRITE(44)(APD1X(I,J),I=1,IIR),(APD1Z(I,J),I=1,IIR),&
     (APD2X(I,J),I=1,IIR),(APD2Z(I,J),I=1,IIR)

      ENDDO
      ELSE
     OPEN(UNIT=44,FILE='GRIN.QAT',FORM='UNFORMATTED',STATUS='OLD')
      READ(44)IIR,JJZ,(Xr(I),I=1,IIR),(xZ(J),J=1,JJZ)
      DO J=1,JJZ
      READ(44)(APD1X(I,J),I=1,IIR),(APD1Z(I,J),I=1,IIR),&
     (APD2X(I,J),I=1,IIR),(APD2Z(I,J),I=1,IIR)
      ENDDO
      ENDIF

    RETURN

    END SUBROUTINE
    SUBROUTINE CREK
    USE COM_VAR
    USE FIC_COM
    IMPLICIT NONE
    INTEGER:: I,J,KJ,KI
    REAL:: AKZ,AKR,ZKJ
      JJZ=124
      IIR=676
    DO 8006 J=1,JJZ
      AKZ=10**(J/10.-10)
     XZ(J)=-AKZ
      KJ=10*(ALOG10(-XZ(J))+10.)
      ZKJ=10*(ALOG10(-XZ(J))+10.)
      KJ=MAX(KJ,2)
      KJ=MIN(KJ,123)
8006 CONTINUE
      XR(1)=1.e-8
    DO 8007 I=2,IIR
      IF(I.LT.81)THEN
      AKR=AMIN1(10**((I-1.)/5-6),4./3.+ABS((I-32.)/3.))
      AKR=10**((I-1.)/10-8)
      ELSE
      AKR=1+ABS((I-81.)/6.)
      ENDIF
      XR(I)=AKR
      IF(AKR.LT.1.)THEN
      KI=10*(ALOG10(AKR+1.E-10)+8)+1
      ELSE
      KI=6*AKR+75
      ENDIF
      KI=MAX(KI,2)
      KI=MIN(KI,674)
8007 CONTINUE

    DO 8009 J=1,JJZ
      DO 8008 I=1,IIR
	CALL VNS(XZ(J),XR(I),I,J)
8008 CONTINUE
8009 CONTINUE

    RETURN

    END SUBROUTINE
!---------------------------------------------------------                                                                    

      SUBROUTINE VNS(AKZ,AKR,I,J) 
      USE COM_VAR
      USE ELEMENTARY_FNS
    IMPLICIT NONE         
      INTEGER:: I,J,IT
      REAL:: AKZ,AKR,PI,CT,ST,CQIT,TETA
      REAL:: QQT(NPINTE),CQT(NPINTE)
      REAL:: FD1JX,FD1JZ,FD2JX,FD2JZ
      COMPLEX:: IM,C1,C2,ZIK,GZ,CEX            
      IM=(0.,1.)                                                                
      PI=4.*ATAN(1.)
      CALL COFINT(CQT,QQT)         
      FD1JX=0.                                                              
      FD1JZ=0.                                                              
      FD2JX=0.                                                              
      FD2JZ=0.                                                              
      DO 30 IT=1,NPINTE                                           
	TETA=QQT(IT)                                                           
	CQIT=CQT(IT)                                                   
	CT=COS(TETA)                                                              
	ST=SIN(TETA)                                                              
	ZIK=AKZ+IM*AKR*CT    
	IF(REAL(ZIK)+30.)2,2,1
      2 CEX=(0.,0.)
	GOTO 3 
      1 CEX=CEXP(ZIK)                                               
      3 GZ=GG(ZIK,CEX)                      
	C1=CQIT*(GZ-1./ZIK)
	C2=CQIT*CEX
	FD1JX=FD1JX+CT*AIMAG(C1)                                 
	FD1JZ=FD1JZ+REAL(C1)                                     
	FD2JX=FD2JX+CT*AIMAG(C2)                                 
	FD2JZ=FD2JZ+REAL(C2)                                     
   30 CONTINUE                                                                
      APD1X(I,J)=FD1JX                                               
      APD1Z(I,J)=FD1JZ                                               
      APD2X(I,J)=FD2JX                                               
      APD2Z(I,J)=FD2JZ                                              
      RETURN 
                                                                   
      END SUBROUTINE  
!-------------------------------------------------------------------

      SUBROUTINE COFINT(CQT,QQT)  
      USE COM_VAR
    IMPLICIT NONE 
      INTEGER :: J
      REAL:: PI,QQT(NPINTE),CQT(NPINTE)

      PI=4.*ATAN(1.)
      DO 160 J=1,NPINTE   
      QQT(J)=-PI/2.+(J-1.)/(NPINTE-1.)*PI
      IF(J-1)161,161,162
  161 CQT(J)=PI/(3.*(NPINTE-1.))
      GOTO 160
  162 IF(J-NPINTE)163,161,161
  163 IF(MOD(J,2))164,165,164
  164 CQT(J)=2./(3.*(NPINTE-1.))*PI
      GOTO 160
  165 CQT(J)=4./(3.*(NPINTE-1.))*PI
  160 CONTINUE

      RETURN                                                                    
      END SUBROUTINE   
!------------------------------------------------------------------

END MODULE
