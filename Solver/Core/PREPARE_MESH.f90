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
    
    IMPLICIT NONE

    INTEGER :: I,K,L,M,N,IJJ,NJJ
    REAL :: XAVER,YAVER,ZAVER
    REAL:: EPPM,D1,T1,T1X,T1Y,T1Z,T2X,T2Y,T2Z,XNQ,XNX,XNY,XNZ
    REAL:: ETA1,ETA2,ETA3,ETA4,ETAN1,ETAN2,ETAN3,ETAN4,ETAO
    REAL:: T1UNX,T1UNY,T1UNZ,T2UNX,T2UNY,T2UNZ
    REAL:: PSCA,XI1,XI2,XI3,XI4,XIN1,XIN2,XIN3,XIN4,XIO
    REAL:: DPOINT,T2,TI1,TT1,TT2,TT3,TT4,BA,BB,BC,BD,DNUL,DST,U
    INTEGER :: NSYM                                             
       
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
        ETAO=-ETA1/3.                                                
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
    END DO
    
      NG=4
      CALL SOMGO(I)
      IF(NG.EQ.1)THEN
		  XGA(1,I)=XG(I)
		  YGA(1,I)=YG(I)
		  ZGA(1,I)=ZG(I)
      ENDIF
      DO L=1,NG
		XJAC(L,I)=XJAC(L,I)/AIRE(I)
      END DO
      PRINT *,'NG = ',NG
      
    RETURN
    
    END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE SOMGO(J)

    USE COM_VAR
    
    IMPLICIT NONE

    INTEGER:: I,J,L,IMETH
    REAL:: XAVER,YAVER,ZAVER
    REAL:: T1X,T1Y,T1Z,T2X,T2Y,T2Z,XNQ,XNX,XNY,XNZ
    REAL:: AT1,AT2,AT3,AT4,AT5,AT6,AT7,AT8,AT9
    REAL:: XL1,XL2,XL3,XL4,X41,X32,X34,X21,Y41,Y32,Y34,Y21
    REAL:: A,B,C,D,AA,BB,CC,DD,EPSN,XG1,YG1
    REAL:: TT,YL1,YL2,YL3,YL4
!~     REAL, DIMENSION(:)::XJAC,XGA,YGA,ZGA
!~     INTEGER:: M1,M2,M3,M4

!~ C
!~ C     PREPARATION DU CALCUL DES COEFFICIENTS D'INFLUENCE:
!~ C     DETERMINATION DE LA POSITION DES POINTS DE GAUSS
!~ C---------------------GAUSS 1 4 9 OU 16--------------------------------
!~ C     ET DU JACOBIEN POUR CHAQUE FACETTE
!~ C      NG =1,4,9 OU 16

      REAL :: XX(16,4),YY(16,4),WG(16,4)

      DATA((XX(I,L),I=1,16),L=1,4)/ &
       16*0., &
      .57735027, .57735027,-.57735027,-.57735027,12*0., &
      .77459667, .77459667, .77459667,3*0.,-.77459667, &
      -.77459667,-.77459667,7*0., &
      .86113631, .86113631, .86113631, .86113631, &
      .33998104, .33998104, .33998104, .33998104, &
      -.33998104,-.33998104,-.33998104,-.33998104, &
      -.86113631,-.86113631,-.86113631,-.86113631/
      DATA((YY(I,L),I=1,16),L=1,4)/ &
       16*0., &
      -.57735027,.57735027,-.57735027,.57735027,12*0., &
      -.77459667,0.,.77459667,-.77459667,0.,.77459667, &
      -.77459667,0.,.77459667,7*0., &
      -.86113363,-.33998104,.33998104,.86113363, &
      -.86113363,-.33998104,.33998104,.86113363, &
      -.86113363,-.33998104,.33998104,.86113363, &
      -.86113363,-.33998104,.33998104,.86113363/
      DATA((WG(I,L),I=1,16),L=1,4)/ &
       1,15*0., &
      .25,.25,.25,.25,12*0., &
      .07716049,.12345679,.07716049,.12345679,.19753086, &
      .12345679,.07716049,.12345679,.07716049,7*0., &
      .30250748E-1,.56712963E-1,.56712963E-1,.30250748E-1, &
      .56712963E-1,.10632333,.10632333,.56712963E-1, &
      .56712963E-1,.10632333,.10632333,.56712963E-1, &
      .30250748E-1,.56712963E-1,.56712963E-1,.30250748E-1/


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
      
      T1X=X(M3(J))-X(M1(J))
      T1Y=Y(M3(J))-Y(M1(J))
      T1Z=Z(M3(J))-Z(M1(J))
      T2X=X(M4(J))-X(M2(J))
      T2Y=Y(M4(J))-Y(M2(J))
      T2Z=Z(M4(J))-Z(M2(J))
      XNX=T2Y*T1Z-T1Y*T2Z
      XNY=T1X*T2Z-T2X*T1Z
      XNZ=T2X*T1Y-T1X*T2Y
      XNQ=SQRT(XNX**2+XNY**2+XNZ**2)
      AT3=XNX/XNQ
      AT6=XNY/XNQ
      AT9=XNZ/XNQ
      TT=SQRT(T1X**2+T1Y**2+T1Z**2)
      AT1=T1X/TT
      AT4=T1Y/TT
      AT7=T1Z/TT
      AT2=AT6*AT7-AT9*AT4
      AT5=AT9*AT1-AT3*AT7
      AT8=AT3*AT4-AT6*AT1
!~       XG=0.25*(X(M1(J))+X(M2(J))+X(M3(J))+X(M4(J)))
!~       YG=0.25*(Y(M1(J))+Y(M2(J))+Y(M3(J))+Y(M4(J)))
!~       ZG=0.25*(Z(M1(J))+Z(M2(J))+Z(M3(J))+Z(M4(J)))
      XL1=AT1*(X(M1(J))-XG(J))+AT4*(Y(M1(J))-YG(J))+AT7*(Z(M1(J))-ZG(J))
      YL1=AT2*(X(M1(J))-XG(J))+AT5*(Y(M1(J))-YG(J))+AT8*(Z(M1(J))-ZG(J))
      XL2=AT1*(X(M2(J))-XG(J))+AT4*(Y(M2(J))-YG(J))+AT7*(Z(M2(J))-ZG(J))
      YL2=AT2*(X(M2(J))-XG(J))+AT5*(Y(M2(J))-YG(J))+AT8*(Z(M2(J))-ZG(J))
      XL3=AT1*(X(M3(J))-XG(J))+AT4*(Y(M3(J))-YG(J))+AT7*(Z(M3(J))-ZG(J))
      YL3=AT2*(X(M3(J))-XG(J))+AT5*(Y(M3(J))-YG(J))+AT8*(Z(M3(J))-ZG(J))
      XL4=AT1*(X(M4(J))-XG(J))+AT4*(Y(M4(J))-YG(J))+AT7*(Z(M4(J))-ZG(J))
      YL4=AT2*(X(M4(J))-XG(J))+AT5*(Y(M4(J))-YG(J))+AT8*(Z(M4(J))-ZG(J))
!~ C
!~ C    DETERMINATION DU JACOBIEN EN CHACUN DES POINTS DE GAUSS
!~ C
      X41=XL4-XL1
      X32=XL3-XL2
      Y34=YL3-YL4
      Y21=YL2-YL1
      Y41=YL4-YL1
      Y32=YL3-YL2
      X34=XL3-XL4
      X21=XL2-XL1
      DO 15 L=1,NG
      A=(1.-YY(L,IMETH))*X41+(1+YY(L,IMETH))*X32
      B=(1.-XX(L,IMETH))*Y21+(1+XX(L,IMETH))*Y34
      C=(1.-XX(L,IMETH))*X21+(1+XX(L,IMETH))*X34
      D=(1.-YY(L,IMETH))*Y41+(1+YY(L,IMETH))*Y32
      XJAC(L,J)=ABS(A*B-C*D)*WG(L,IMETH)*.25
   15 CONTINUE
!~ C
!~ C    COORDONNES ABSOLUES DES POINTS DE GAUSS
!~ C
      EPSN=0.
      DO 20 L=1,NG
      AA=.25*(1-XX(L,IMETH))*(1-YY(L,IMETH))
      BB=.25*(1-XX(L,IMETH))*(1+YY(L,IMETH))
      CC=.25*(1+XX(L,IMETH))*(1+YY(L,IMETH))
      DD=.25*(1+XX(L,IMETH))*(1-YY(L,IMETH))
      XG1=AA*XL1+BB*XL2+CC*XL3+DD*XL4
      YG1=AA*YL1+BB*YL2+CC*YL3+DD*YL4
      XGA(L,J)=XG(J)+AT1*XG1+AT2*YG1+AT3*EPSN
      YGA(L,J)=YG(J)+AT4*XG1+AT5*YG1+AT6*EPSN
      ZGA(L,J)=ZG(J)+AT7*XG1+AT8*YG1+AT9*EPSN
  20  CONTINUE

      RETURN
      END SUBROUTINE

END MODULE
