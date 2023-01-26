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
    RETURN
    END SUBROUTINE

END MODULE