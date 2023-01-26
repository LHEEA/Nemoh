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
!
!--------------------------------------------------------------------------------------
MODULE COMPUTE_GREEN_INFD

  USE COM_VAR
  USE ELEMENTARY_FNS

  IMPLICIT NONE

  CONTAINS   
!-------------------------------------------------------------------------------!
      SUBROUTINE VAVINFD(KKK,XGI,YGI,ZGI,ISP,IFP) 
                                                  
      INTEGER:: ISP,IFP
      INTEGER:: KKK,I,J,IMXX,MK,NJJ,JJ,L,MH,MY,MZ,MJJ
      INTEGER:: KK(5)
      REAL:: DH,XOI,YOI,ZOI,XGI,YGI,ZGI                      
      REAL:: RR(5),DRX(5),DRY(5),DRZ(5)                    
      REAL:: PI,PI4,DPI,QPI
      REAL:: TXN(5),TYN(5),TZN(5),AIJS(4),VXS(4),VYS(4),VZS(4)
      REAL:: A3J,A6J,A9J,ALDEN,ANL,ANLX,ANLY,ANLZ,ANTX,ANTY,ANTZ
      REAL:: ARG,ASRO,AT,ATX,ATY,ATZ,DAT,DDK,DEN,DENL,DENT,DK,DLOGG
      REAL:: ANT,DNL,DNT,DNTX,DNTY,DNTZ,DR,DS,GY,GYX,GYZ,GZ,PJ,QJ,RJ,RO,SGN,W
      REAL:: GYY,XOJ,YOJ,ZOJ

      PI4=ATAN(1.)
      PI=4.*PI4                                                             
      DPI=2.*PI                                                                 
      QPI=4.*PI
      NJJ=2*(NSYMY+1)                                                           
      DH=0.
      MK=(-1)**(KKK+1)
      IF(KKK.EQ.1)IMXX=IMX
      IF(KKK.EQ.2)IMXX=IXX

       I=ISP
       J=IFP                                                  
       XOI=XGI                                                             
       YOI=YGI
       ZOI=ZGI
       IF(KKK.EQ.1)THEN 
	 IF(I.LE.IMX)THEN
	   IF(ZGI.GT.ZER)THEN
	     ZOI=ZER
	   ELSE
	     ZOI=ZGI
	   ENDIF
	 ENDIF
       ELSE
	 IF(I.LE.IMX)THEN
	   IF(ZGI.GT.ZER)THEN
	     ZOI=20*ZER
	   ELSE
	     ZOI=ZGI
	   ENDIF
	 ENDIF
       ENDIF
                                                         
	DO 25 JJ=1,NJJ
	MJJ=(-1)**(JJ+1)                                                          
	MY=(-1)**(JJ/3+2)
	MZ=(-1)**(JJ/2+2)
	MH=(1-(-1)**(JJ/2+2))/2
	XOJ=XG(J)                                                               
	YOJ=YG(J)*MY                                                            
	ZOJ=ZG(J)*MZ-DH*MH                                               
	A3J=XN(J)                                                               
	A6J=YN(J)*MY                                                            
	A9J=ZN(J)*MZ                                                           
	RO=SQRT((XOI-XOJ)**2+(YOI-YOJ)**2+(ZOI-ZOJ)**2)
	IF(RO.GT.7.*TDIS(J))THEN                                               
	  AIJS(JJ)=AIRE(J)/RO                                                     
	  ASRO=AIJS(JJ)/RO**2                                                       
	  VXS(JJ)=-(XOI-XOJ)*ASRO                                              
	  VYS(JJ)=-(YOI-YOJ)*ASRO                                              
	  VZS(JJ)=-(ZOI-ZOJ)*ASRO                                              
	ELSE
	  AIJS(JJ)=0.                                                               
	  VXS(JJ)=0.                                                                
	  VYS(JJ)=0.                                                                
	  VZS(JJ)=0.                                                                
	  KK(1)=M1(J)                                                               
	  KK(2)=M2(J)                                                               
	  KK(3)=M3(J)                                                               
	  KK(4)=M4(J)                                                               
	  KK(5)=KK(1)                                                               
	  DO 211 L=1,4                                                              
	    TXN(L)=X(KK(L))                                                            
	    TYN(L)=Y(KK(L))*MY                                                         
	    TZN(L)=Z(KK(L))*MZ-DH*MH                                              
      211 CONTINUE                                                                  
	  TXN(5)=TXN(1)                                                               
	  TYN(5)=TYN(1)                                                               
	  TZN(5)=TZN(1)                                                               
	  DO 212 L=1,4                                                              
	    RR(L)=SQRT((XOI-TXN(L))**2+(YOI-TYN(L))**2+(ZOI-TZN(L))**2)          
	    DRX(L)=(XOI-TXN(L))/RR(L)                                                
	    DRY(L)=(YOI-TYN(L))/RR(L)                                                
	    DRZ(L)=(ZOI-TZN(L))/RR(L)                                           
      212 CONTINUE             
	  RR(5)=RR(1)
	  DRX(5)=DRX(1)                                                           
	  DRY(5)=DRY(1)                                                           
	  DRZ(5)=DRZ(1)                                                           
	  GZ=(XOI-XOJ)*A3J+(YOI-YOJ)*A6J+(ZOI-ZOJ)*A9J
	  DO 29 L=1,4                                                               
	  DK=SQRT((TXN(L+1)-TXN(L))**2+(TYN(L+1)-TYN(L))**2+(TZN(L+1)-TZN(L))**2)
	    IF(DK.GE.1.E-3*TDIS(J))THEN                                             
	    PJ=(TXN(L+1)-TXN(L))/DK                                                     
	    QJ=(TYN(L+1)-TYN(L))/DK                                                     
	    RJ=(TZN(L+1)-TZN(L))/DK                                                     
	    GYX=A6J*RJ-A9J*QJ                                                     
	    GYY=A9J*PJ-A3J*RJ                                                     
	    GYZ=A3J*QJ-A6J*PJ                                                     
	    GY=(XOI-TXN(L))*GYX+(YOI-TYN(L))*GYY+(ZOI-TZN(L))*GYZ
	    SGN=SIGN(1.,GZ)                                                           
	    DDK=2.*DK                                                                 
	    ANT=GY*DDK                                                                
	    DNT=(RR(L+1)+RR(L))**2-DK*DK+2.*ABS(GZ)*(RR(L+1)+RR(L))          
	    ARG=ANT/DNT                                                               
	    ANL=RR(L+1)+RR(L)+DK                                                      
	    DNL=RR(L+1)+RR(L)-DK                                                      
	    DEN=ANL/DNL                                                               
	    ALDEN=ALOG(DEN)                                                           
	      IF(ABS(GZ).GE.1.E-4*TDIS(J))THEN                                       
	      AT=ATAN(ARG)
	      ELSE
	      AT=0.                                                                     
	      ENDIF                                                                  
	    AIJS(JJ)=AIJS(JJ)+GY*ALDEN-2.*ABS(GZ)*AT                                  
	    DAT=2.*AT*SGN                                                             
	    ANTX=GYX*DDK                                                              
	    ANTY=GYY*DDK                                                              
	    ANTZ=GYZ*DDK                                                              
	    ANLX=DRX(L+1)+DRX(L)                                                      
	    ANLY=DRY(L+1)+DRY(L)                                                      
	    ANLZ=DRZ(L+1)+DRZ(L)                                                      
	    DR=2.*(RR(L+1)+RR(L)+ABS(GZ))                                             
	    DS=2.*(RR(L+1)+RR(L))*SGN                                                 
	    DNTX=DR*ANLX+A3J*DS                                                     
	    DNTY=DR*ANLY+A6J*DS                                                     
	    DNTZ=DR*ANLZ+A9J*DS                                                     
	    DENL=ANL*DNL                                                              
	    DENT=ANT*ANT+DNT*DNT                                                      
	    ATX=(ANTX*DNT-DNTX*ANT)/DENT                                              
	    ATY=(ANTY*DNT-DNTY*ANT)/DENT                                              
	    ATZ=(ANTZ*DNT-DNTZ*ANT)/DENT                                              
	    DLOGG=(DNL-ANL)/DENL                                                       
	    VXS(JJ)=VXS(JJ)+GYX*ALDEN+GY*ANLX*DLOGG-2.*ABS(GZ)*ATX-DAT*A3J
	    VYS(JJ)=VYS(JJ)+GYY*ALDEN+GY*ANLY*DLOGG-2.*ABS(GZ)*ATY-DAT*A6J             
	    VZS(JJ)=VZS(JJ)+GYZ*ALDEN+GY*ANLZ*DLOGG-2.*ABS(GZ)*ATZ-DAT*A9J            
	    ENDIF
      29 CONTINUE                                                                  
	  IF(I.EQ.J.AND.JJ.EQ.1)THEN                                               
	  VXS(1)=VXS(1)-DPI*A3J                                                   
	  VYS(1)=VYS(1)-DPI*A6J                                                   
	  VZS(1)=VZS(1)-DPI*A9J                                                   
	  ELSE                                             
	  AIJS(JJ)=AIJS(JJ)*MJJ                                                     
	  VXS(JJ)=VXS(JJ)*MJJ                                                       
	  VYS(JJ)=VYS(JJ)*MJJ                                                       
	  VZS(JJ)=VZS(JJ)*MJJ                                                       
	  ENDIF
	ENDIF
    25 CONTINUE                                                                  
	  IF(NSYMY.EQ.1)THEN                                                       
	    W=AIJS(1)-MK*(AIJS(2)+AIJS(3))+AIJS(4)                                   
	    FSP=-W/QPI                                                            
	    W=AIJS(1)-MK*(AIJS(2)-AIJS(3))-AIJS(4)                                  
	    FSM=-W/QPI                                                            
	    W=VXS(1)-MK*(VXS(2)+VXS(3))+VXS(4)                                 
	    VSXP=-W/QPI                                                           
	    W=VYS(1)-MK*(VYS(2)+VYS(3))+VYS(4)                               
	    VSYP=-W/QPI                                                            
	    W=VZS(1)-MK*(VZS(2)+VZS(3))+VZS(4)                              
	    VSZP=-W/QPI                                                            
	    W=VXS(1)-MK*(VXS(2)-VXS(3))-VXS(4)                               
	    VSXM=-W/QPI                                                            
	    W=VYS(1)-MK*(VYS(2)-VYS(3))-VYS(4)                             
	    VSYM=-W/QPI                                                            
	    W=VZS(1)-MK*(VZS(2)-VZS(3))-VZS(4)                             
	    VSZM=-W/QPI                                                          
	  ELSE                                                                      
	    W=AIJS(1)-MK*AIJS(2)                                                     
	    FSP=-W/QPI                                                            
	    FSM=FSP                                                           
	    W=VXS(1)-MK*VXS(2)                                                    
	    VSXP=-W/QPI                                                           
	    VSXM=VSXP                                                          
	    W=VYS(1)-MK*VYS(2)                                                     
	    VSYP=-W/QPI                                                            
	    VSYM=VSYP                                                           
	    W=VZS(1)-MK*VZS(2)                                                     
	    VSZP=-W/QPI                                                            
	    VSZM=VSZP                                                                                                                           
	  ENDIF                                                                     

       RETURN                                                                    
      END SUBROUTINE 
!--------------------------------------------------------------------------!

  SUBROUTINE VNSINFD(KKK,ISP,IFP,XGI,YGI,ZGI)  

      INTEGER:: ISP,IFP
      INTEGER:: I,J,L,JJ,KK(5),NJJ,IJUMP,BX,KI,KJ,IT,KKK,IMXX
      REAL:: XGI,YGI,ZGI                                                    
      REAL:: FS1(NFA,2),FS2(NFA,2)                        
      REAL:: VSX1(NFA,2),VSY1(NFA,2),VSZ1(NFA,2)                                
      REAL:: VSX2(NFA,2),VSY2(NFA,2),VSZ2(NFA,2)                                        
      REAL:: PI4,PI,DPI,QPI,DPI2
      REAL:: WH,WR,AK0
      REAL:: EPS,ADPI,ADPI2,AKAIR,AKDPI,AKDPI2,AKP4
      REAL:: AKR,AKZ,DD,PSURR,QJJJ,RRR,ZZZ,YMJJJ,ZMIII,CVX,CVY,VR1,VR2
      REAL:: CSK,DSK,EPZ,F1,F2,F3,CT,ST,TETA
      REAL:: PD1X,PD2X,PD1Z,PD2Z,SIK,SQ,VZ1,VZ2,XL1,XL2,XL3,ZL1,ZL2,ZL3                   
      COMPLEX*8 ZIJ(5),CEX(5),GZ(5)
      COMPLEX*8 ZI,C1,C2,ASD,BSD,CSD,ZA,ZB,ZVS

      PI4=ATAN(1.)
      PI=4.*PI4                                                             
      DPI=2.*PI                                                                 
      QPI=4.*PI
      ZI=(0.,1.)
      DPI2=2.*PI**2                                                             
      EPS=0.0001
                                                             
      IF(KKK.EQ.1)IMXX=IMX
      IF(KKK.EQ.2)IMXX=IXX
      WH=DPI/T  
      WR=WH                                           
      IF(ABS(WR).LT.1.E-4)THEN
       WRITE(*,*)'ABS(WR)  = ',ABS(WR),' < 1.E-4'
       STOP
      ENDIF

      AK0=WR**2/G
      NJJ=NSYMY+1                                     
      IJUMP=0                                                           

      I=ISP  !source point
      J=IFP  !field point
                                               
      IF(I.LE.IMX)THEN !1B  
	ZMIII=ZGI
	IF(ZGI.GT.ZER)ZMIII=ZER
      ENDIF!1F

 	DO 22 JJ=1,NJJ                                                            
 	 BX=(-1)**(JJ+1)                                                           
 	    IF(ZGI.LT.ZER.AND.ZG(J).LT.ZER)THEN !2B
	      YMJJJ=YG(J)*BX                                                           
	      QJJJ=YN(J)*BX                                                              
	      RRR=SQRT((XGI-XG(J))**2+(YGI-YMJJJ)**2)
	      AKR=AK0*RRR               
	      ZZZ=ZMIII+ZG(J)
	      AKZ=AK0*ZZZ
	      DD=SQRT(RRR**2+ZZZ**2)
	      IF(DD.GT.EPS)THEN
		PSURR=PI/(AK0*DD)**3
	      ELSE
		PSURR=0.
	      ENDIF
	      IF(AKZ.GT.-1.5E-6)THEN !3B
		IF(IJUMP.NE.1)THEN
		  WRITE(*,*)'AKZ < -1.5 E-6'
		  IJUMP=1
		ENDIF
	      ELSE !3E
		IF(AKZ.GT.-16.)THEN !4B
		  IF(AKR.LT.99.7)THEN !5B
		    IF(AKZ.LT.-1.E-2)THEN
		      KJ=8*(ALOG10(-AKZ)+4.5)
		    ELSE
		      KJ=5*(ALOG10(-AKZ)+6)
		    ENDIF
		    KJ=MAX(KJ,2)
		    KJ=MIN(KJ,45)
		    IF(AKR.LT.1.)THEN
		      KI=5*(ALOG10(AKR+1.E-20)+6)+1
		    ELSE
		      KI=3*AKR+28
		    ENDIF        
		    KI=MAX(KI,2)
		    KI=MIN(KI,327)
		    XL1=PL2(XR(KI),XR(KI+1),XR(KI-1),AKR)
		    XL2=PL2(XR(KI+1),XR(KI-1),XR(KI),AKR)
		    XL3=PL2(XR(KI-1),XR(KI),XR(KI+1),AKR)
		    ZL1=PL2(XZ(KJ),XZ(KJ+1),XZ(KJ-1),AKZ)
		    ZL2=PL2(XZ(KJ+1),XZ(KJ-1),XZ(KJ),AKZ)
		    ZL3=PL2(XZ(KJ-1),XZ(KJ),XZ(KJ+1),AKZ)
		    F1=XL1*APD1Z(KI-1,KJ-1)+XL2*APD1Z(KI,KJ-1)+XL3*APD1Z(KI+1,KJ-1)
		    F2=XL1*APD1Z(KI-1,KJ)+XL2*APD1Z(KI,KJ)+XL3*APD1Z(KI+1,KJ)
		    F3=XL1*APD1Z(KI-1,KJ+1)+XL2*APD1Z(KI,KJ+1)+XL3*APD1Z(KI+1,KJ+1) 
		    PD1Z=ZL1*F1+ZL2*F2+ZL3*F3
		    F1=XL1*APD2Z(KI-1,KJ-1)+XL2*APD2Z(KI,KJ-1)+XL3*APD2Z(KI+1,KJ-1)
		    F2=XL1*APD2Z(KI-1,KJ)+XL2*APD2Z(KI,KJ)+XL3*APD2Z(KI+1,KJ)
		    F3=XL1*APD2Z(KI-1,KJ+1)+XL2*APD2Z(KI,KJ+1)+XL3*APD2Z(KI+1,KJ+1) 
		    PD2Z=ZL1*F1+ZL2*F2+ZL3*F3
		  ELSE !5E
		    EPZ=EXP(AKZ)
		    AKP4=AKR-PI4                                        
		    SQ=SQRT(DPI/AKR)
		    CSK=COS(AKP4)
		    SIK=SIN(AKP4)
		    PD1Z=PSURR*AKZ-PI*EPZ*SQ*SIK
		    PD2Z=EPZ*SQ*CSK
		  ENDIF !5F
		    VZ1=PD1Z-PSURR*AKZ
		    VZ2=PD2Z
		ELSE !4E
		  PD1Z=PSURR*AKZ
		  PD2Z=0.
		  VZ1=0.
		  VZ2=0.           
		ENDIF !4F
 	      ENDIF !3F
 	      FS1(J,JJ)=PD1Z                                        
 	      FS2(J,JJ)=PD2Z                                    
! 	      IF(I.LE.IMX)THEN !6B
		IF(RRR.GT.EPS)THEN !7B
		  IF(AKZ.LE.-1.5E-6)THEN !8B
		    IF(AKZ.GT.-16.)THEN !9B
		      IF(AKR.LT.99.7)THEN !10B
			F1=XL1*APD1X(KI-1,KJ-1)+XL2*APD1X(KI,KJ-1)+XL3*APD1X(KI+1,KJ-1)
			F2=XL1*APD1X(KI-1,KJ)+XL2*APD1X(KI,KJ)+XL3*APD1X(KI+1,KJ)
			F3=XL1*APD1X(KI-1,KJ+1)+XL2*APD1X(KI,KJ+1)+XL3*APD1X(KI+1,KJ+1) 
			PD1X=ZL1*F1+ZL2*F2+ZL3*F3
			F1=XL1*APD2X(KI-1,KJ-1)+XL2*APD2X(KI,KJ-1)+XL3*APD2X(KI+1,KJ-1)
			F2=XL1*APD2X(KI-1,KJ)+XL2*APD2X(KI,KJ)+XL3*APD2X(KI+1,KJ)
			F3=XL1*APD2X(KI-1,KJ+1)+XL2*APD2X(KI,KJ+1)+XL3*APD2X(KI+1,KJ+1) 
			PD2X=ZL1*F1+ZL2*F2+ZL3*F3
		      ELSE !10E
			DSK=0.5/AKR
!                       PD1X=-PSURR*AKR-PI*EPZ*SQ*(CSK-DSK*SIK) !coorection par GD le 17/09/2010
                        PD1X=-PSURR*AKR+PI*EPZ*SQ*(CSK-DSK*SIK)
			PD2X=EPZ*SQ*(SIK+DSK*CSK)
		      ENDIF !10F
		      VR1=-PD1X-PSURR*AKR
		      VR2=-PD2X
		    ELSE !9E
		      PD1X=-PSURR*AKR
		      PD2X=0.
		      VR1=0.
		      VR2=0.
		    ENDIF !9F
		  ENDIF !8F
		    CVX=(XGI-XG(J))/RRR
		    CVY=(YGI-YMJJJ)/RRR
		    VSX1(J,JJ)=VR1*CVX                                   
		    VSX2(J,JJ)=VR2*CVX                                    
		    VSY1(J,JJ)=VR1*CVY                                 
		    VSY2(J,JJ)=VR2*CVY                                    
		    VSZ1(J,JJ)=VZ1                                        
		    VSZ2(J,JJ)=VZ2
		ELSE !7E  
		    VSX1(J,JJ)=0.
		    VSX2(J,JJ)=0.                                         
		    VSY1(J,JJ)=0.                                      
		    VSY2(J,JJ)=0.                                         
		    VSZ1(J,JJ)=VZ1                                        
		    VSZ2(J,JJ)=VZ2
		ENDIF !7F
! 	      ENDIF !6F	  
	   ELSE !2E            
	      FS1(J,JJ)=0.                                        
	      FS2(J,JJ)=0.                                    
	      VSX1(J,JJ)=0.
	      VSX2(J,JJ)=0.                                         
	      VSY1(J,JJ)=0.                                      
	      VSY2(J,JJ)=0.                                         
	      VSZ1(J,JJ)=0.                                       
	      VSZ2(J,JJ)=0.
	      KK(1)=M1(J)
	      KK(2)=M2(J)
	      KK(3)=M3(J)
	      KK(4)=M4(J)
	      KK(5)=KK(1)
	   DO 30 IT=1,NQ
	      TETA=QQ(IT)
	      CT=COS(TETA)
	      ST=SIN(TETA)
	      DO 20 L=1,4
		ZIJ(L)=AK0*(ZMIII+Z(KK(L))+ZI*((XGI-X(KK(L)))*CT+(YGI-Y(KK(L))*BX)*ST))
		IF(REAL(ZIJ(L)).GT.-25.)THEN
		  CEX(L)=CEXP(ZIJ(L))
		ELSE
		  CEX(L)=(0.,0.)
		ENDIF
		GZ(L)=GG(ZIJ(L),CEX(L))
	    20 CONTINUE
	      ZIJ(5)=ZIJ(1)
	      CEX(5)=CEX(1)
	      GZ(5)=GZ(1)
	      ZA=(0.,0.)
	      ZB=(0.,0.)
	      ZVS=(0.,0.)
	      DO 37 L=1,4
		C1=ZIJ(L+1)-ZIJ(L)
		IF(ABS(AIMAG(C1)).LT.EPS.AND.ABS(REAL(C1)).LT.EPS)THEN
		  ASD=(GZ(L+1)+GZ(L))*0.5
		  BSD=(CEX(L+1)+CEX(L))*0.5
		  CSD=ASD-(0.5/ZIJ(L+1)+0.5/ZIJ(L))
		ELSE
		  CSD=(GZ(L+1)-GZ(L))/C1
		  ASD=(GZ(L+1)-GZ(L)+CLOG(ZIJ(L+1)/ZIJ(L)))/C1
		  BSD=(CEX(L+1)-CEX(L))/C1
		ENDIF
		C2=(YN(J)*BX-ZI*ZN(J)*ST)*(X(KK(L+1))-X(KK(L)))-&
	      & (XN(J)-ZI*ZN(J)*CT)*(Y(KK(L+1))-Y(KK(L)))*BX
		ZA=ZA+ASD*C2
		ZB=ZB+BSD*C2
		ZVS=ZVS+CSD*C2
 	  37 CONTINUE
	      FS1(J,JJ)=FS1(J,JJ)+CQ(IT)*REAL(ZA)
	      FS2(J,JJ)=FS2(J,JJ)+CQ(IT)*REAL(ZB)
	      VSX1(J,JJ)=VSX1(J,JJ)-CQ(IT)*CT*AIMAG(ZVS)
	      VSX2(J,JJ)=VSX2(J,JJ)-CQ(IT)*CT*AIMAG(ZB)
	      VSY1(J,JJ)=VSY1(J,JJ)-CQ(IT)*ST*AIMAG(ZVS)
	      VSY2(J,JJ)=VSY2(J,JJ)-CQ(IT)*ST*AIMAG(ZB)
	      VSZ1(J,JJ)=VSZ1(J,JJ)+CQ(IT)*REAL(ZVS)
	      VSZ2(J,JJ)=VSZ2(J,JJ)+CQ(IT)*REAL(ZB)
 	30 CONTINUE                 
	      AKAIR=AK0*AIRE(J)            
	      FS1(J,JJ)=FS1(J,JJ)*BX/AKAIR
	      FS2(J,JJ)=FS2(J,JJ)*BX/AKAIR
	      VSX1(J,JJ)=VSX1(J,JJ)*BX/AKAIR
	      VSX2(J,JJ)=VSX2(J,JJ)*BX/AKAIR
	      VSY1(J,JJ)=VSY1(J,JJ)*BX/AKAIR
	      VSY2(J,JJ)=VSY2(J,JJ)*BX/AKAIR
	      VSZ1(J,JJ)=VSZ1(J,JJ)*BX/AKAIR
	      VSZ2(J,JJ)=VSZ2(J,JJ)*BX/AKAIR
 	   ENDIF !2F                                                             
    22 CONTINUE    
                                                              
      IF(NSYMY.EQ.1)THEN !11B        
	  AKAIR=AK0*AIRE(J)                  
	  ADPI2=AKAIR/DPI2
	  ADPI=AKAIR/DPI
	  SM1=FSM-(FS1(J,1)-FS1(J,2))*ADPI2                             
	  SP1=FSP-(FS1(J,1)+FS1(J,2))*ADPI2
	  SM2=-(FS2(J,1)-FS2(J,2))*ADPI                                   
	  SP2=-(FS2(J,1)+FS2(J,2))*ADPI
	  AKDPI2=ADPI2*AK0
	  AKDPI=ADPI*AK0
	  VSXP1=VSXP-(VSX1(J,1)+VSX1(J,2))*AKDPI2                             
	  VSXM1=VSXM-(VSX1(J,1)-VSX1(J,2))*AKDPI2                            
	  VSYP1=VSYP-(VSY1(J,1)+VSY1(J,2))*AKDPI2                           
	  VSYM1=VSYM-(VSY1(J,1)-VSY1(J,2))*AKDPI2                            
	  VSZP1=VSZP-(VSZ1(J,1)+VSZ1(J,2))*AKDPI2                            
	  VSZM1=VSZM-(VSZ1(J,1)-VSZ1(J,2))*AKDPI2
	  VSXP2=-(VSX2(J,1)+VSX2(J,2))*AKDPI                                    
	  VSXM2=-(VSX2(J,1)-VSX2(J,2))*AKDPI                                    
	  VSYP2=-(VSY2(J,1)+VSY2(J,2))*AKDPI                                    
	  VSYM2=-(VSY2(J,1)-VSY2(J,2))*AKDPI                                    
	  VSZP2=-(VSZ2(J,1)+VSZ2(J,2))*AKDPI                                    
	  VSZM2=-(VSZ2(J,1)-VSZ2(J,2))*AKDPI                                 
      ELSE   !11E
	  AKAIR=AK0*AIRE(J)                  
	  ADPI2=AKAIR/DPI2
	  ADPI=AKAIR/DPI
	  SP1=FSP-FS1(J,1)*ADPI2                                              
	  SM1=SP1                                                             
	  SP2=-FS2(J,1)*ADPI                                                    
	  SM2=SP2
	  AKDPI2=ADPI2*AK0
	  AKDPI=ADPI*AK0
	  VSXP1=VSXP-VSX1(J,1)*AKDPI2                                        
	  VSXM1=VSXP1                                                        
	  VSYP1=VSYP-VSY1(J,1)*AKDPI2                                        
	  VSYM1=VSYP1                                                        
	  VSZP1=VSZP-VSZ1(J,1)*AKDPI2                                       
	  VSZM1=VSZP1                                                        
	  VSXP2=-VSX2(J,1)*AKDPI                                            
	  VSXM2=VSXP2                                                        
	  VSYP2=-VSY2(J,1)*AKDPI                                             
	  VSYM2=VSYP2                                                      
	  VSZP2=-VSZ2(J,1)*AKDPI                                             
	  VSZM2=VSZP2                                                         
      ENDIF !11F
                            
      RETURN                                                                    
      END SUBROUTINE
!-------------------------------------------------------------------------------!
      SUBROUTINE CINT_INFD(AKK,N)

      INTEGER::I,J,N,NPIN
      REAL::AKK,PI                      
      REAL:: Q8(8),CQ8(8),Q12(12),CQ12(12),Q16(16),CQ16(16)
      REAL:: Q24(24),CQ24(24),Q32(32),CQ32(32)
      DATA Q8/.4801449,.3983332,.2627662,.09171732,4*0./
      DATA CQ8/.05061427,.1111905,.1568533,.1813419,4*0./
      DATA Q12/.4907803,.4520586,.3849513,.2936589,.1839157,.06261670,&
     & 6*0./
      DATA CQ12/.2358766E-1,.5346966E-1,.8003916E-1,.1015837,&
     &.1167462,.1245735,6*0./
      DATA Q16/.4947004,.4722875,.4328156,.3777022,.3089381,.2290084,&
     &.1408017,.04750625,8*0./
      DATA CQ16/.01357622,.03112676,.04757925,.06231448,.07479799,&
     &.08457826,.09130170,.09472530,8*0./
      DATA Q24/.4975936,.4873642,.469137,.4432077,.4100009,.3700621,&
     &.3240468,.2727107,.2168967,.1575213,.09555943,.032028446,12*0./
      DATA CQ24/.6170615E-2,.1426569E-1,.2213872E-1,.2964929E-1,&
     &.366732E-1,.4309508E-1,.4880932E-1,.5372213E-1,.5775283E-1,&
     &.6083523E-1,.6291873E-1,.6396909E-1,12*0./
      DATA Q32/.4986319,.4928057,.4823811,.4674530,.4481605,.4246838,&
     &.3972418,.3660910,.3315221,.2938578,.2534499,.2106756,.1659343,&
     &.1196436,.07223598,.02415383,16*0./
      DATA CQ32/.350930E-2,.8137197E-2,.1269603E-1,.1713693E-1,&
     & .2141794E-1,.2549903E-1,.2934204E-1,.3291111E-1,&
     & .3617289E-1,.3909694E-1,.4165596E-1,.4382604E-1,&
     & .4558693E-1,.4692219E-1,.4781936E-1,.4827004E-1,16*0./

      NPIN=101
      PI=4.*ATAN(1.)

      IF(AKK-.4)801,801,810
  810 IF(AKK-2.5)802,802,811
  811 IF(AKK-4.)803,803,812
  812 IF(AKK-8.)804,804,813
  813 IF(AKK-25.)805,805,806

  801 N=8
      DO 700 I=1,4
	Q8(I)=Q8(I)
	Q8(9-I)=-Q8(I)
	CQ8(I)=CQ8(I)
	CQ8(9-I)=CQ8(I)
  700 CONTINUE
      DO 110 J=1,N
	QQ(J)=Q8(J)*PI
	CQ(J)=CQ8(J)*PI
  110 CONTINUE
      GOTO 815

  802 N=12
      DO 701 I=1,6
	Q12(I)=Q12(I)
	Q12(13-I)=-Q12(I)
	CQ12(I)=CQ12(I)
	CQ12(13-I)=CQ12(I)
  701 CONTINUE
      DO 120 J=1,N
	QQ(J)=Q12(J)*PI
	CQ(J)=CQ12(J)*PI
  120 CONTINUE
      GOTO 815

  803 N=16
      DO 702 I=1,8
	Q16(I)=Q16(I)
	Q16(17-I)=-Q16(I)
	CQ16(I)=CQ16(I)
	CQ16(17-I)=CQ16(I)
  702 CONTINUE
      DO 130 J=1,N
	QQ(J)=Q16(J)*PI
	CQ(J)=CQ16(J)*PI
  130 CONTINUE
      GOTO 815

  804 N=24
      DO 703 I=1,12
	Q24(I)=Q24(I)
	Q24(25-I)=-Q24(I)
	CQ24(I)=CQ24(I)
	CQ24(25-I)=CQ24(I)
  703 CONTINUE
      DO 140 J=1,N
	QQ(J)=Q24(J)*PI
	CQ(J)=CQ24(J)*PI
  140 CONTINUE
      GOTO 815

  805 N=32
      DO 704 I=1,16
	Q32(I)=Q32(I)
	Q32(33-I)=-Q32(I)
	CQ32(I)=CQ32(I)
	CQ32(33-I)=CQ32(I)
  704 CONTINUE
      DO 150 J=1,N
	QQ(J)=Q32(J)*PI
	CQ(J)=CQ32(J)*PI
  150 CONTINUE

  815 CONTINUE
 !     WRITE(*,3001)N,AKK
 3001 FORMAT( 2X,'INTEGRATION PAR QG',I2 &
     & /2X,'K0 ADIMENSIONNEL PAR RAPPORT LA PLUS GRANDE DIMENSION = '&
     & ,1PE16.6)
      GOTO 807
  806 CONTINUE

      N=51
      IF(AKK.GT.40.)THEN
 !     WRITE(*,822)
  822 FORMAT(10X,'L''INTEGRATION RISQUE D''ETRE PEU PRECISE')
      N=NPIN
      ENDIF
      DO 160 J=1,N
      QQ(J)=-PI/2.+(J-1.)/(N-1.)*PI
      IF(J-1)161,161,162
  161 CQ(J)=PI/(3.*(N-1.))
      GOTO 160
  162 IF(J-N)163,161,161
  163 IF(MOD(J,2))164,165,164
  164 CQ(J)=2./(3.*(N-1.))*PI
      GOTO 160
  165 CQ(J)=4./(3.*(N-1.))*PI
  160 CONTINUE
!      WRITE(*,3002)N,AKK
 3002 FORMAT(2X,'SIMPSON ',I3,' POINTS'&
     & /2X,'K0 ADIMENSIONNEL PAR RAPPORT A LA PLUS GRANDE DIMENSION = '&
     &,1PE16.6)
  807 CONTINUE
      RETURN
      END SUBROUTINE
!------------------------------------------------------------------------

      REAL FUNCTION PL2(U1,U2,U3,XU)
	REAL::U1,U2,U3,XU
	PL2=((XU-U1)*(XU-U2))/((U3-U1)*(U3-U2))
	RETURN
      END FUNCTION

END MODULE