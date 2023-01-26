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
MODULE COMPUTE_GREEN_FD

  USE COM_VAR
  USE ELEMENTARY_FNS
  IMPLICIT NONE

  CONTAINS
!-------------------------------------------------------------------------------!
      SUBROUTINE VAVFD(KKK,XGI,YGI,ZGI,ISP,IFP)                                                     
                                                  
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
      DH=2*Depth
      MK=(-1)**(KKK+1)
      IF(KKK.EQ.1)IMXX=IMX
      IF(KKK.EQ.2)IMXX=IXX

       I=ISP
       J=IFP                                                  
       XOI=XGI                                                             
       YOI=YGI
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
!-------------------------------------------------------------------------------!
                                                                      
      SUBROUTINE VNSFD(AM0,AMH,NEXP,ISP,IFP,XGI,YGI,ZGI)                     
  
      INTEGER::ISP,IFP
      REAL:: AM0,AMH,XGI,YGI,ZGI                                                   
      REAL:: FS1(NFA,2),FS2(NFA,2)    
      INTEGER::I,J,JJ,IJUMP,NJJ,NEXP,NEXP1
      INTEGER::KK(5),BX,BKL,IT,KE,KI,KJ1,KJ2,KJ3,KJ4,KL,L
      REAL::H,A,ADPI,ADPI2,AKH,COE1,COE2,COE3,COE4,EPS
      REAL::PI,PI4,DPI,QPI,WH                
      REAL:: VSX1(NFA,2),VSY1(NFA,2),VSZ1(NFA,2)                                
      REAL:: VSX2(NFA,2),VSY2(NFA,2),VSZ2(NFA,2)
      REAL:: ZMIII,ACT,AKP4,AKR,AKZ1,AKZ2,AKZ3,AKZ4,AQT,ASRO1,ASRO2
      REAL:: ASRO3,ASRO4,C1V3,C2V3,COF1,COF2,COF3,COF4
      REAL:: CSK,CT,ST,CVX,CVY,DD1,DD2,DD3,DD4,DSK,DXL,DYL
      REAL:: EPZ1,EPZ2,EPZ3,EPZ4,F1,F2,F3,FFS1,FFS2,FFS3,FFS4,FTS1,FTS2,FTS3,FTS4              
      REAL:: OM,PD1X1,PD1X2,PD1X3,PD1X4,PD1Z1,PD1Z2, PD1Z3,PD1Z4,PD2X1,PD2X2
      REAL:: PD2X3,PD2X4,PD2Z1,PD2Z2,PD2Z3,PD2Z4
      REAL:: PSK,PSR1,PSR2,PSR3,PSR4,PSURR1,PSURR2,PSURR3,PSURR4,QJJJ,QTQQ
      REAL:: RO1,RO2,RO3,RO4,RRR,RR1,RR2,RR3,RR4,SCDS,SSDS,STSS
      REAL:: SQ,SIK,SCK,TETA,VR21,VR22,VR23,VR24
      REAL:: VX1,VX2,VX3,VX4,VY1,VY2,VY3,VY4,VZ1,VZ2,VZ3,VZ4
      REAL:: VXS1,VXS2,VXS3,VXS4,VZ11,VZ12,VZ13,VZ14,VZ21,VZ22,VZ23,VZ24
      REAL:: VYS1,VYS2,VYS3,VYS4,VZS1,VZS2,VZS3,VZS4
      REAL:: XL1,XL2,XL3,XPG,YPG,YMJJJ,ZL11,ZL12,ZL13,ZL14
      REAL:: ZL21,ZL22,ZL23,ZL24,ZL31,ZL32,ZL33,ZL34,ZPG1,ZPG2,ZPG3,ZPG4
      REAL:: ZZZ1,ZZZ2,ZZZ3,ZZZ4                   
      COMPLEX ZIJ(4,5),CEX(4,5),GZ(4,5),CL(4,5)
      COMPLEX:: S1,ZV1,S2,ZV2,ZV3,ZV4,Z1,Z0,CL1,CL0,G1,G0
      COMPLEX:: CEX1,CEX0,AUX,ZAM,ZI,ZIRS,ZIRC
      REAL:: ZP(4),XFT(5),YFT(5),ZFT(5)
      PI4=ATAN(1.)
      PI=4.*ATAN(1.)                                                             
      DPI=2.*PI
      QPI=4.*PI
      H=Depth
      ZI=(0.,1.)
      WH=DPI/T                                                                                                                                                                           
      AKH=AMH*TANH(AMH)                                                                                
      A=(AMH+AKH)**2/(H*(AMH**2-AKH**2+AKH))
      NEXP1=NEXP+1                                                              
      AMBDA(NEXP1)=0.                                                           
      AR(NEXP1)=2. 
      ADPI2=-A/(8.*PI**2)
      ADPI=-A/(8*PI)                                                              
      COE1=ADPI2/AM0                                                       
      COE2=ADPI/AM0                                                           
      COE3=ADPI2                                                     
      COE4=ADPI                                                           
      ijump=1
      EPS=0.0001 
      NJJ=NSYMY+1

      I=ISP  !source point N
      J=IFP  ! observation point

	IF(I.LE.IMX)THEN
	  ZMIII=ZGI
	  IF(ZGI.GT.ZER)ZMIII=20*ZER
	ELSE                    
	  ZMIII=20*ZER
	ENDIF
             
	DO 7122 JJ=1,NJJ 
	  BX=(-1)**(JJ+1)                                                                                                            
	    IF(ZGI.LT.ZER.AND.ZG(J).LT.ZER)THEN   !0000B
	      QJJJ=BX*YN(J)                                                              
	      YMJJJ=BX*YG(J)                                                          
	      COF1=COE3*AIRE(J)
	      COF2=COE4*AIRE(J)
	      COF3=AM0*COF1                                             
	      COF4=AM0*COF2                                             
	      RRR=SQRT((XGI-XG(J))**2+(YGI-YMJJJ)**2)
	      AKR=AM0*RRR
	      ZZZ1=ZMIII+ZG(J)
	      AKZ1=AM0*ZZZ1
	      DD1=SQRT(RRR**2+ZZZ1**2)
	      IF(DD1.GT.EPS)THEN
		RR1=AM0*DD1
		PSR1=PI/RR1
		PSURR1=PI/RR1**3
	      ELSE
		PSR1=0.
		PSURR1=0.
	      ENDIF

	      IF(AKZ1.GT.-1.5E-6)THEN             !0001B
		IF(IJUMP.NE.1)THEN
		  WRITE(*,*)'AKZ < -1.5 E-6'
		  IJUMP=1
		ENDIF
	      ELSE                               !0001E
		IF(AKZ1.GT.-16.)THEN             !0002B
		  IF(AKR.LT.99.7)THEN            !0003B
		    IF(AKZ1.LT.-1.E-2)THEN
		      KJ1=8*(ALOG10(-AKZ1)+4.5)
		    ELSE
		      KJ1=5*(ALOG10(-AKZ1)+6)
		    ENDIF
		    KJ1=MAX(KJ1,2)
		    KJ1=MIN(KJ1,45)
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
		    ZL11=PL2(XZ(KJ1),XZ(KJ1+1),XZ(KJ1-1),AKZ1)
		    ZL21=PL2(XZ(KJ1+1),XZ(KJ1-1),XZ(KJ1),AKZ1)
		    ZL31=PL2(XZ(KJ1-1),XZ(KJ1),XZ(KJ1+1),AKZ1)
		    F1=XL1*APD1Z(KI-1,KJ1-1)+XL2*APD1Z(KI,KJ1-1)+XL3*APD1Z(KI+1,KJ1-1)
		    F2=XL1*APD1Z(KI-1,KJ1)+XL2*APD1Z(KI,KJ1)+XL3*APD1Z(KI+1,KJ1)
		    F3=XL1*APD1Z(KI-1,KJ1+1)+XL2*APD1Z(KI,KJ1+1)+XL3*APD1Z(KI+1,KJ1+1) 
		    PD1Z1=ZL11*F1+ZL21*F2+ZL31*F3
		    F1=XL1*APD2Z(KI-1,KJ1-1)+XL2*APD2Z(KI,KJ1-1)+XL3*APD2Z(KI+1,KJ1-1)
		    F2=XL1*APD2Z(KI-1,KJ1)+XL2*APD2Z(KI,KJ1)+XL3*APD2Z(KI+1,KJ1)
		    F3=XL1*APD2Z(KI-1,KJ1+1)+XL2*APD2Z(KI,KJ1+1)+XL3*APD2Z(KI+1,KJ1+1) 
		    PD2Z1=ZL11*F1+ZL21*F2+ZL31*F3
		  ELSE                           !0003E
		    EPZ1=EXP(AKZ1)
		    AKP4=AKR-PI4                                        
		    SQ=SQRT(DPI/AKR)
		    CSK=COS(AKP4)
		    SIK=SIN(AKP4)
		    PSK=PI*SQ*SIK
		    SCK=SQ*CSK  
		    PD1Z1=PSURR1*AKZ1-PSK*EPZ1
		    PD2Z1=EPZ1*SCK
		  ENDIF                         !0003F
		  VZ11=PD1Z1-PSURR1*AKZ1
		  VZ21=PD2Z1
		ELSE                            !0002E
		  PD1Z1=PSURR1*AKZ1
		  PD2Z1=0.
		  VZ11=0.
		  VZ21=0.
		ENDIF                           !0002F
	      ENDIF                             !0001F
	      ZZZ2=ZG(J)-ZMIII-2*H
	      AKZ2=AM0*ZZZ2
	      DD2=SQRT(RRR**2+ZZZ2**2)
	      IF(DD2.GT.EPS)THEN
		RR2=AM0*DD2
		PSR2=PI/RR2
		PSURR2=PI/RR2**3
	      ELSE
		PSR2=0.
		PSURR2=0.
	      ENDIF

	      IF(AKZ2.GT.-1.5E-6)THEN            !0004B
		IF(IJUMP.NE.1)THEN
		  WRITE(*,*)'AKZ < -1.5 E-6'
		  IJUMP=1
		ENDIF
	      ELSE                               !0004E
		IF(AKZ2.GT.-16.)THEN             !0005B
		  IF(AKR.LT.99.7)THEN            !0006B
		    IF(AKZ2.LT.-1.E-2)THEN
		      KJ2=8*(ALOG10(-AKZ2)+4.5)
		    ELSE
		      KJ2=5*(ALOG10(-AKZ2)+6)
		    ENDIF
		    KJ2=MAX(KJ2,2)
		    KJ2=MIN(KJ2,45)
		    ZL12=PL2(XZ(KJ2),XZ(KJ2+1),XZ(KJ2-1),AKZ2)
		    ZL22=PL2(XZ(KJ2+1),XZ(KJ2-1),XZ(KJ2),AKZ2)
		    ZL32=PL2(XZ(KJ2-1),XZ(KJ2),XZ(KJ2+1),AKZ2)
		    F1=XL1*APD1Z(KI-1,KJ2-1)+XL2*APD1Z(KI,KJ2-1)+XL3*APD1Z(KI+1,KJ2-1)
		    F2=XL1*APD1Z(KI-1,KJ2)+XL2*APD1Z(KI,KJ2)+XL3*APD1Z(KI+1,KJ2)
		    F3=XL1*APD1Z(KI-1,KJ2+1)+XL2*APD1Z(KI,KJ2+1)+XL3*APD1Z(KI+1,KJ2+1) 
		    PD1Z2=ZL12*F1+ZL22*F2+ZL32*F3
		    F1=XL1*APD2Z(KI-1,KJ2-1)+XL2*APD2Z(KI,KJ2-1)+XL3*APD2Z(KI+1,KJ2-1)
		    F2=XL1*APD2Z(KI-1,KJ2)+XL2*APD2Z(KI,KJ2)+XL3*APD2Z(KI+1,KJ2)
		    F3=XL1*APD2Z(KI-1,KJ2+1)+XL2*APD2Z(KI,KJ2+1)+XL3*APD2Z(KI+1,KJ2+1) 
		    PD2Z2=ZL12*F1+ZL22*F2+ZL32*F3
		  ELSE                          !0006E
		    EPZ2=EXP(AKZ2)
		    PD1Z2=PSURR2*AKZ2-PSK*EPZ2
		    PD2Z2=EPZ2*SCK
		  ENDIF                         !0006F
		  VZ12=PD1Z2-PSURR2*AKZ2
		  VZ22=PD2Z2
		ELSE                            !0005E
		  PD1Z2=PSURR2*AKZ2
		  PD2Z2=0.
		  VZ12=0.
		  VZ22=0.
		ENDIF                            !0005F
	      ENDIF                              !0004F
	      ZZZ3=ZMIII-ZG(J)-2*H
	      AKZ3=AM0*ZZZ3
	      DD3=SQRT(RRR**2+ZZZ3**2)
	      IF(DD3.GT.EPS)THEN
		RR3=AM0*DD3
		PSR3=PI/RR3
		PSURR3=PI/RR3**3
	      ELSE
		PSR3=0.
		PSURR3=0.
	      ENDIF

	      IF(AKZ3.GT.-1.5E-6)THEN   !0007B
		IF(IJUMP.NE.1)THEN
		  WRITE(*,*)'AKZ < -1.5 E-6'
		  IJUMP=1
		ENDIF
	      ELSE                    !0007E
		IF(AKZ3.GT.-16.)THEN    !0008B
		  IF(AKR.LT.99.7)THEN     !0009B
		    IF(AKZ3.LT.-1.E-2)THEN
		      KJ3=8*(ALOG10(-AKZ3)+4.5)
		    ELSE
		      KJ3=5*(ALOG10(-AKZ3)+6)
		    ENDIF
		    KJ3=MAX(KJ3,2)
		    KJ3=MIN(KJ3,45)
		    ZL13=PL2(XZ(KJ3),XZ(KJ3+1),XZ(KJ3-1),AKZ3)
		    ZL23=PL2(XZ(KJ3+1),XZ(KJ3-1),XZ(KJ3),AKZ3)
		    ZL33=PL2(XZ(KJ3-1),XZ(KJ3),XZ(KJ3+1),AKZ3)
		    F1=XL1*APD1Z(KI-1,KJ3-1)+XL2*APD1Z(KI,KJ3-1)+XL3*APD1Z(KI+1,KJ3-1)
		    F2=XL1*APD1Z(KI-1,KJ3)+XL2*APD1Z(KI,KJ3)+XL3*APD1Z(KI+1,KJ3)
		    F3=XL1*APD1Z(KI-1,KJ3+1)+XL2*APD1Z(KI,KJ3+1)+XL3*APD1Z(KI+1,KJ3+1) 
		    PD1Z3=ZL13*F1+ZL23*F2+ZL33*F3
		    F1=XL1*APD2Z(KI-1,KJ3-1)+XL2*APD2Z(KI,KJ3-1)+XL3*APD2Z(KI+1,KJ3-1)
		    F2=XL1*APD2Z(KI-1,KJ3)+XL2*APD2Z(KI,KJ3)+XL3*APD2Z(KI+1,KJ3)
		    F3=XL1*APD2Z(KI-1,KJ3+1)+XL2*APD2Z(KI,KJ3+1)+XL3*APD2Z(KI+1,KJ3+1) 
		    PD2Z3=ZL13*F1+ZL23*F2+ZL33*F3
		  ELSE                    !0009E
		    EPZ3=EXP(AKZ3)
		    PD1Z3=PSURR3*AKZ3-PSK*EPZ3
		    PD2Z3=EPZ3*SCK
		  ENDIF                   !0009F
		  VZ13=PD1Z3-PSURR3*AKZ3
		  VZ23=PD2Z3                
		ELSE                    !0008E
		  PD1Z3=PSURR3*AKZ3
		  PD2Z3=0.
		  VZ13=0.
		  VZ23=0.
		ENDIF                    !0008F
	      ENDIF                    !0007F
	      ZZZ4=-ZG(J)-ZMIII-4*H
	      AKZ4=AM0*ZZZ4
	      DD4=SQRT(RRR**2+ZZZ4**2)
	      IF(DD4.GT.EPS)THEN
		RR4=AM0*DD4
		PSR4=PI/RR4
		PSURR4=PI/RR4**3
	      ELSE
		PSR4=0.
		PSURR4=0.
	      ENDIF

	      IF(AKZ4.GT.-1.5E-6)THEN       !0010B
		IF(IJUMP.NE.1)THEN
		  WRITE(*,*)'AKZ < -1.5 E-6'
		  IJUMP=1
		ENDIF
	      ELSE                         !0010E
		IF(AKZ4.GT.-16.)THEN       !0011B
		  IF(AKR.LT.99.7)THEN      !0012B
		    IF(AKZ4.LT.-1.E-2)THEN
		      KJ4=8*(ALOG10(-AKZ4)+4.5)
		    ELSE
		      KJ4=5*(ALOG10(-AKZ4)+6)
		    ENDIF
		    KJ4=MAX(KJ4,2)
		    KJ4=MIN(KJ4,45)
		    ZL14=PL2(XZ(KJ4),XZ(KJ4+1),XZ(KJ4-1),AKZ4)
		    ZL24=PL2(XZ(KJ4+1),XZ(KJ4-1),XZ(KJ4),AKZ4)
		    ZL34=PL2(XZ(KJ4-1),XZ(KJ4),XZ(KJ4+1),AKZ4)
		    F1=XL1*APD1Z(KI-1,KJ4-1)+XL2*APD1Z(KI,KJ4-1)+XL3*APD1Z(KI+1,KJ4-1)
		    F2=XL1*APD1Z(KI-1,KJ4)+XL2*APD1Z(KI,KJ4)+XL3*APD1Z(KI+1,KJ4)
		    F3=XL1*APD1Z(KI-1,KJ4+1)+XL2*APD1Z(KI,KJ4+1)+XL3*APD1Z(KI+1,KJ4+1) 
		    PD1Z4=ZL14*F1+ZL24*F2+ZL34*F3
		    F1=XL1*APD2Z(KI-1,KJ4-1)+XL2*APD2Z(KI,KJ4-1)+XL3*APD2Z(KI+1,KJ4-1)
		    F2=XL1*APD2Z(KI-1,KJ4)+XL2*APD2Z(KI,KJ4)+XL3*APD2Z(KI+1,KJ4)
		    F3=XL1*APD2Z(KI-1,KJ4+1)+XL2*APD2Z(KI,KJ4+1)+XL3*APD2Z(KI+1,KJ4+1) 
		    PD2Z4=ZL14*F1+ZL24*F2+ZL34*F3
		  ELSE                     !0012E
		    EPZ4=EXP(AKZ4)
		    PD1Z4=PSURR4*AKZ4-PSK*EPZ4
		    PD2Z4=EPZ4*SCK
		  ENDIF                    !0012F
		  VZ14=PD1Z4-PSURR4*AKZ4
		  VZ24=PD2Z4
		ELSE                       !0011E
		  PD1Z4=PSURR4*AKZ4
		  PD2Z4=0.
		  VZ14=0.
		  VZ24=0.
		ENDIF                      !0011F
	      ENDIF                      !0010F
	      QTQQ=PD1Z1+PD1Z2+PD1Z3+PD1Z4
	      FS1(J,JJ)=COF1*(QTQQ-PSR1-PSR2-PSR3-PSR4)              
	      STSS=PD2Z1+PD2Z2+PD2Z3+PD2Z4
	      FS2(J,JJ)=COF2*STSS
                                
	      IF(I.LE.IMX)THEN   !501 B
		IF(RRR.GT.EPS)THEN !601B
		  IF(AKZ1.LE.-1.5E-6)THEN !701B
		    IF(AKZ1.GT.-16.)THEN  !801B
		      IF(AKR.LT.99.7)THEN
			F1=XL1*APD1X(KI-1,KJ1-1)+XL2*APD1X(KI,KJ1-1)+XL3*APD1X(KI+1,KJ1-1)
			F2=XL1*APD1X(KI-1,KJ1)+XL2*APD1X(KI,KJ1)+XL3*APD1X(KI+1,KJ1)
			F3=XL1*APD1X(KI-1,KJ1+1)+XL2*APD1X(KI,KJ1+1)+XL3*APD1X(KI+1,KJ1+1) 
			PD1X1=ZL11*F1+ZL21*F2+ZL31*F3
			F1=XL1*APD2X(KI-1,KJ1-1)+XL2*APD2X(KI,KJ1-1)+XL3*APD2X(KI+1,KJ1-1)
			F2=XL1*APD2X(KI-1,KJ1)+XL2*APD2X(KI,KJ1)+XL3*APD2X(KI+1,KJ1)
			F3=XL1*APD2X(KI-1,KJ1+1)+XL2*APD2X(KI,KJ1+1)+XL3*APD2X(KI+1,KJ1+1) 
			PD2X1=ZL11*F1+ZL21*F2+ZL31*F3
		      ELSE
			DSK=0.5/AKR
			SCDS=PI*SQ*(CSK-DSK*SIK)
			SSDS=SQ*(SIK+DSK*CSK)
			PD1X1=-PSURR1*AKR-EPZ1*SCDS
			PD2X1=EPZ1*SSDS
		      ENDIF
		      VR21=-PD2X1
		    ELSE !801E
		      PD1X1=-PSURR1*AKR
		      PD2X1=0.
		      VR21=0.
		    ENDIF   !801F
		  ENDIF !701F
		  IF(AKZ2.LE.-1.5E-6)THEN
		    IF(AKZ2.GT.-16.)THEN
		      IF(AKR.LT.99.7)THEN
			F1=XL1*APD1X(KI-1,KJ2-1)+XL2*APD1X(KI,KJ2-1)+XL3*APD1X(KI+1,KJ2-1)
			F2=XL1*APD1X(KI-1,KJ2)+XL2*APD1X(KI,KJ2)+XL3*APD1X(KI+1,KJ2)
			F3=XL1*APD1X(KI-1,KJ2+1)+XL2*APD1X(KI,KJ2+1)+XL3*APD1X(KI+1,KJ2+1) 
			PD1X2=ZL12*F1+ZL22*F2+ZL32*F3
			F1=XL1*APD2X(KI-1,KJ2-1)+XL2*APD2X(KI,KJ2-1)+XL3*APD2X(KI+1,KJ2-1)
			F2=XL1*APD2X(KI-1,KJ2)+XL2*APD2X(KI,KJ2)+XL3*APD2X(KI+1,KJ2)
			F3=XL1*APD2X(KI-1,KJ2+1)+XL2*APD2X(KI,KJ2+1)+XL3*APD2X(KI+1,KJ2+1) 
			PD2X2=ZL12*F1+ZL22*F2+ZL32*F3
		      ELSE
			PD1X2=-PSURR2*AKR-EPZ2*SCDS
			PD2X2=EPZ2*SSDS
		      ENDIF
		      VR22=-PD2X2
		    ELSE
		      PD1X2=-PSURR2*AKR
		      PD2X2=0.
		      VR22=0.
		    ENDIF
		  ENDIF
		  IF(AKZ3.LE.-1.5E-6)THEN
		    IF(AKZ3.GT.-16.)THEN
		      IF(AKR.LT.99.7)THEN
			F1=XL1*APD1X(KI-1,KJ3-1)+XL2*APD1X(KI,KJ3-1)+XL3*APD1X(KI+1,KJ3-1)
			F2=XL1*APD1X(KI-1,KJ3)+XL2*APD1X(KI,KJ3)+XL3*APD1X(KI+1,KJ3)
			F3=XL1*APD1X(KI-1,KJ3+1)+XL2*APD1X(KI,KJ3+1)+XL3*APD1X(KI+1,KJ3+1) 
			PD1X3=ZL13*F1+ZL23*F2+ZL33*F3
			F1=XL1*APD2X(KI-1,KJ3-1)+XL2*APD2X(KI,KJ3-1)+XL3*APD2X(KI+1,KJ3-1)
			F2=XL1*APD2X(KI-1,KJ3)+XL2*APD2X(KI,KJ3)+XL3*APD2X(KI+1,KJ3)
			F3=XL1*APD2X(KI-1,KJ3+1)+XL2*APD2X(KI,KJ3+1)+XL3*APD2X(KI+1,KJ3+1) 
			PD2X3=ZL13*F1+ZL23*F2+ZL33*F3
		      ELSE
			PD1X3=-PSURR3*AKR-EPZ3*SCDS
			PD2X3=EPZ3*SSDS
		      ENDIF
		      VR23=-PD2X3
		    ELSE
		      PD1X3=-PSURR3*AKR
		      PD2X3=0.
		      VR23=0.
		    ENDIF
		  ENDIF
		  IF(AKZ4.LE.-1.5E-6)THEN
		    IF(AKZ4.GT.-16.)THEN
		      IF(AKR.LT.99.7)THEN
			F1=XL1*APD1X(KI-1,KJ4-1)+XL2*APD1X(KI,KJ4-1)+XL3*APD1X(KI+1,KJ4-1)
			F2=XL1*APD1X(KI-1,KJ4)+XL2*APD1X(KI,KJ4)+XL3*APD1X(KI+1,KJ4)
			F3=XL1*APD1X(KI-1,KJ4+1)+XL2*APD1X(KI,KJ4+1)+XL3*APD1X(KI+1,KJ4+1) 
			PD1X4=ZL14*F1+ZL24*F2+ZL34*F3
			F1=XL1*APD2X(KI-1,KJ4-1)+XL2*APD2X(KI,KJ4-1)+XL3*APD2X(KI+1,KJ4-1)
			F2=XL1*APD2X(KI-1,KJ4)+XL2*APD2X(KI,KJ4)+XL3*APD2X(KI+1,KJ4)
			F3=XL1*APD2X(KI-1,KJ4+1)+XL2*APD2X(KI,KJ4+1)+XL3*APD2X(KI+1,KJ4+1) 
			PD2X4=ZL14*F1+ZL24*F2+ZL34*F3
		      ELSE
			PD1X4=-PSURR4*AKR-EPZ4*SCDS
			PD2X4=EPZ4*SSDS
		      ENDIF
		      VR24=-PD2X4
		    ELSE
		      PD1X4=-PSURR4*AKR
		      PD2X4=0.
		      VR24=0.
		    ENDIF  
		  ENDIF 
		  C1V3=-COF3*(PD1X1+PD1X2+PD1X3+PD1X4)
		  C2V3=COF4*(VR21+VR22+VR23+VR24)
		  CVX=(XGI-XG(J))/RRR
		  CVY=(YGI-YMJJJ)/RRR
		  VSX1(J,JJ)=C1V3*CVX                                 
		  VSX2(J,JJ)=C2V3*CVX                                    
		  VSY1(J,JJ)=C1V3*CVY                                 
		  VSY2(J,JJ)=C2V3*CVY                                    
		  ELSE                                 !601E
		    VSX1(J,JJ)=0.                                 
		    VSX2(J,JJ)=0.                                    
		    VSY1(J,JJ)=0.                                 
		    VSY2(J,JJ)=0.                                    
		  ENDIF                                 !601F
		    VSZ1(J,JJ)=COF3*(PD1Z1-PD1Z2+PD1Z3-PD1Z4)
		    VSZ2(J,JJ)=COF4*(VZ21-VZ22+VZ23-VZ24)  
		ENDIF                                   !501F
	
		XPG=XGI-XG(J)
		YPG=YGI-YMJJJ
		ACT=-0.5*AIRE(J)/QPI
		DO 7234 KE=1,NEXP1           
		  AQT=ACT*AR(KE)
		  ZPG1=ZMIII-2.*H+H*AMBDA(KE)-ZG(J)
		  ZPG2=-ZMIII-H*AMBDA(KE)-ZG(J)
		  ZPG3=-ZMIII-4.*H+H*AMBDA(KE)-ZG(J)
		  ZPG4=ZMIII+2.*H-H*AMBDA(KE)-ZG(J)
		  RR1=RRR**2+ZPG1**2
		  RO1=SQRT(RR1)
		  IF(RO1.GT.EPS)THEN
		    FTS1=AQT/RO1                                                     
		    ASRO1=FTS1/RR1                                              
		  ELSE
		    FTS1=0.
		    ASRO1=0.
		  ENDIF
		  IF(I.LE.IMX)THEN
		    VXS1=-XPG*ASRO1                                                         
		    VYS1=-YPG*ASRO1                                                     
		    VZS1=-ZPG1*ASRO1                                                         
		  ENDIF
		    RR2=RRR**2+ZPG2**2
		    RO2=SQRT(RR2)
		  IF(RO2.GT.EPS)THEN
		    FTS2=AQT/RO2                                                     
		    ASRO2=FTS2/RR2                                              
		  ELSE
		    FTS2=0.
		    ASRO2=0.
		  ENDIF
		  IF(I.LE.IMX)THEN
		    VXS2=-XPG*ASRO2                                                         
		    VYS2=-YPG*ASRO2                                                     
		    VZS2=-ZPG2*ASRO2                                                         
		  ENDIF
		  RR3=RRR**2+ZPG3**2
		  RO3=SQRT(RR3)
		  IF(RO3.GT.EPS)THEN
		    FTS3=AQT/RO3                                                     
		    ASRO3=FTS3/RR3                                              
		  ELSE
		    FTS3=0.
		    ASRO3=0.
		  ENDIF
		  IF(I.LE.IMX)THEN
		    VXS3=-XPG*ASRO3                                                         
		    VYS3=-YPG*ASRO3                                                     
		    VZS3=-ZPG3*ASRO3                                                         
		  ENDIF 
		    RR4=RRR**2+ZPG4**2
		    RO4=SQRT(RR4)
		  IF(RO4.GT.EPS)THEN
		    FTS4=AQT/RO4                                                     
		    ASRO4=FTS4/RR4                                              
		  ELSE
		    FTS4=0.
		    ASRO4=0.
		  ENDIF
		  IF(I.LE.IMX)THEN
		    VXS4=-XPG*ASRO4                                                         
		    VYS4=-YPG*ASRO4                                                     
		    VZS4=-ZPG4*ASRO4                                                         
		  ENDIF
		  FS1(J,JJ)=FS1(J,JJ)+FTS1+FTS2+FTS3+FTS4    
		  IF(I.LE.IMX)THEN
		    VSX1(J,JJ)=VSX1(J,JJ)+(VXS1+VXS2+VXS3+VXS4)                             
		    VSY1(J,JJ)=VSY1(J,JJ)+(VYS1+VYS2+VYS3+VYS4)
		    VSZ1(J,JJ)=VSZ1(J,JJ)+(VZS1-VZS2-VZS3+VZS4)                             
		  ENDIF
	  7234 CONTINUE
	    ELSE !0000E

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
		  OM=(XGI-X(KK(L)))*CT+(YGI-BX*Y(KK(L)))*ST
		  ZIJ(1,L)=AM0*(Z(KK(L))+ZMIII+ZI*OM)
		  ZIJ(2,L)=AM0*(Z(KK(L))-ZMIII-2*H+ZI*OM)
		  ZIJ(3,L)=AM0*(ZMIII-Z(KK(L))-2*H+ZI*OM)
		  ZIJ(4,L)=AM0*(-Z(KK(L))-ZMIII-4*H+ZI*OM)
		  DO 23 KL=1,4
		  IF(REAL(ZIJ(KL,L)).GT.-25.)THEN
		  CEX(KL,L)=CEXP(ZIJ(KL,L))
		  ELSE
		  CEX(KL,L)=(0.,0.)
		  ENDIF
		  GZ(KL,L)=GG(ZIJ(KL,L),CEX(KL,L))
		  CL(KL,L)=CLOG(-ZIJ(KL,L))
		  23 CONTINUE
		  20 CONTINUE
		  DO 24 KL=1,4
		  ZIJ(KL,5)=ZIJ(KL,1)
		  CEX(KL,5)=CEX(KL,1)
		  GZ(KL,5)=GZ(KL,1)
		  CL(KL,5)=CL(KL,1)
		  24 CONTINUE
		  S1=(0.,0.)
		  S2=(0.,0.)
		  ZV1=(0.,0.)
		  ZV2=(0.,0.)
		  ZV3=(0.,0.)
		  ZV4=(0.,0.)
		  ZIRS=ZI*ZN(J)*ST
		  ZIRC=ZI*ZN(J)*CT
		  DO 40 L=1,4
		    DXL=(X(KK(L+1))-X(KK(L)))
		    DYL=(Y(KK(L+1))-Y(KK(L)))*BX
		    DO 50 KL=1,4
		      BKL=(-1)**(KL+1)
		      IF(KL.LT.3)THEN
			AUX=DXL*(YN(J)*BX-ZIRS)-DYL*(XN(J)-ZIRC)
		      ELSE
			AUX=DXL*(-YN(J)*BX-ZIRS)-DYL*(-XN(J)-ZIRC)
		      ENDIF
		      Z1=ZIJ(KL,L+1)
		      Z0=ZIJ(KL,L)
		      CL1=CL(KL,L+1)
		      CL0=CL(KL,L)
		      G1=GZ(KL,L+1)
		      G0=GZ(KL,L)
		      CEX1=CEX(KL,L+1)
		      CEX0=CEX(KL,L)
		      ZAM=Z1-Z0
		      IF(ABS(AIMAG(ZAM)).LT.EPS.AND.ABS(REAL(ZAM)).LT.EPS)THEN
			S1=S1+AUX*(G1+G0+CL1+CL0)*0.5
			ZV1=ZV1+AUX*BKL*(G1+G0)*0.5
			ZV3=ZV3+AUX*(G1+G0)*0.5
			S2=S2+AUX*(CEX1+CEX0)*0.5
			ZV2=ZV2+AUX*BKL*(CEX1+CEX0)*0.5
			ZV4=ZV4+AUX*(CEX1+CEX0)*0.5
		      ELSE
			S1=S1+AUX*(G1-G0+CL1-CL0+Z1*CL1-Z0*CL0-ZAM)/ZAM
			ZV1=ZV1+AUX*BKL*(G1-G0+CL1-CL0)/ZAM
			ZV3=ZV3+AUX*(G1-G0+CL1-CL0)/ZAM
			S2=S2+AUX*(CEX1-CEX0)/ZAM
			ZV2=ZV2+AUX*BKL*(CEX1-CEX0)/ZAM
			ZV4=ZV4+AUX*(CEX1-CEX0)/ZAM
		      ENDIF
		50 CONTINUE
	      40 CONTINUE
		  FS1(J,JJ)=FS1(J,JJ)+CQ(IT)*REAL(S1)*COE1*BX
		  FS2(J,JJ)=FS2(J,JJ)+CQ(IT)*REAL(S2)*COE2*BX
		  VSX1(J,JJ)=VSX1(J,JJ)-CQ(IT)*CT*AIMAG(ZV3)*COE3*BX
		  VSX2(J,JJ)=VSX2(J,JJ)-CQ(IT)*CT*AIMAG(ZV4)*COE4*BX
		  VSY1(J,JJ)=VSY1(J,JJ)-CQ(IT)*ST*AIMAG(ZV3)*COE3*BX
		  VSY2(J,JJ)=VSY2(J,JJ)-CQ(IT)*ST*AIMAG(ZV4)*COE4*BX
		  VSZ1(J,JJ)=VSZ1(J,JJ)+CQ(IT)*REAL(ZV1)*COE3*BX
		  VSZ2(J,JJ)=VSZ2(J,JJ)+CQ(IT)*REAL(ZV2)*COE4*BX
	    30 CONTINUE
		ZP(1)=-ZMIII
		ZP(2)=ZMIII+2*H
		ZP(3)=ZMIII-2.*H
		ZP(4)=-ZMIII-4*H
		DO 48 L=1,5
		  XFT(L)=X(KK(L))
		  YFT(L)=Y(KK(L))
		  ZFT(L)=Z(KK(L))
	    48 CONTINUE
		DO 51 KE=1,NEXP1
		  CALL VSD(XFT,YFT,ZFT,JJ,XN(J),YN(J),ZN(J),AIRE(J),TDIS(J),XG(J),&
		  & YG(J),ZG(J),XGI,YGI,ZP(1)-H*AMBDA(KE),FFS1,VX1,VY1,VZ1)
		  CALL VSD(XFT,YFT,ZFT,JJ,XN(J),YN(J),ZN(J),AIRE(J),TDIS(J),XG(J),&
		  & YG(J),ZG(J),XGI,YGI,ZP(2)-H*AMBDA(KE),FFS2,VX2,VY2,VZ2)
		  CALL VSD(XFT,YFT,ZFT,JJ,XN(J),YN(J),ZN(J),AIRE(J),TDIS(J),XG(J),&
		  & YG(J),ZG(J),XGI,YGI,ZP(3)+H*AMBDA(KE),FFS3,VX3,VY3,VZ3)
		  CALL VSD(XFT,YFT,ZFT,JJ,XN(J),YN(J),ZN(J),AIRE(J),TDIS(J),XG(J),&
		  & YG(J),ZG(J),XGI,YGI,ZP(4)+H*AMBDA(KE),FFS4,VX4,VY4,VZ4)
		  FS1(J,JJ)=FS1(J,JJ)+(FFS1+FFS2+FFS3+FFS4)*AR(KE)
		  IF(I.LE.IMX)THEN
		    VSX1(J,JJ)=VSX1(J,JJ)+(VX1+VX2+VX3+VX4)*AR(KE)
		    VSY1(J,JJ)=VSY1(J,JJ)+(VY1+VY2+VY3+VY4)*AR(KE)
		    VSZ1(J,JJ)=VSZ1(J,JJ)+(-VZ1+VZ2+VZ3-VZ4)*AR(KE)
		  ENDIF
	    51 CONTINUE

	    ENDIF !0000F
    7122 CONTINUE 

      IF(NSYMY.EQ.1)THEN  !101B
        
	  SM1=FSM+FS1(J,1)-FS1(J,2)                                           
	  SP1=FSP+FS1(J,1)+FS1(J,2)                                           
	  SM2=FS2(J,1)-FS2(J,2)                                                  
	  SP2=FS2(J,1)+FS2(J,2)                                                  
	  VSXP1=VSXP+VSX1(J,1)+VSX1(J,2)                                  
	  VSXM1=VSXM+VSX1(J,1)-VSX1(J,2)                                      
	  VSYP1=VSYP+VSY1(J,1)+VSY1(J,2)                                      
	  VSYM1=VSYM+VSY1(J,1)-VSY1(J,2)                                      
	  VSZP1=VSZP+VSZ1(J,1)+VSZ1(J,2)                                      
	  VSZM1=VSZM+VSZ1(J,1)-VSZ1(J,2)                                      
	  VSXP2=VSX2(J,1)+VSX2(J,2)                                              
	  VSXM2=VSX2(J,1)-VSX2(J,2)                                              
	  VSYP2=VSY2(J,1)+VSY2(J,2)                                              
	  VSYM2=VSY2(J,1)-VSY2(J,2)                                              
	  VSZP2=VSZ2(J,1)+VSZ2(J,2)                                              
	  VSZM2=VSZ2(J,1)-VSZ2(J,2)   
      ELSE      !101E
	  SP1=FSP+FS1(J,1)                                                    
	  SM1=SP1                                                             
	  SP2=FS2(J,1)                                                          
	  SM2=SP2                                                             
	  VSXP1=VSXP+VSX1(J,1)                                                
	  VSXM1=VSXP1                                                         
	  VSYP1=VSYP+VSY1(J,1)                                                
	  VSYM1=VSYP1                                                         
	  VSZP1=VSZP+VSZ1(J,1)                                                
	  VSZM1=VSZP1                                                         
	  VSXP2=VSX2(J,1)                                                        
	  VSXM2=VSXP2                                                         
	  VSYP2=VSY2(J,1)                                                        
	  VSYM2=VSYP2                                                        
	  VSZP2=VSZ2(J,1)                                                        
	  VSZM2=VSZP2                                                              
      ENDIF !101F

      RETURN                                                                    
      END SUBROUTINE  
!------------------------------------------------------------------------

      SUBROUTINE CINT_FD(AKK,N)

      INTEGER::I,J,N,NPIN
      REAL::AKK,CQ(101),QQ(101),PI                      
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
!      WRITE(*,3001)N,AKK
 3001 FORMAT( 2X,'INTEGRATION PAR QG',I2 &
     & /2X,'K0 ADIMENSIONNEL PAR RAPPORT LA PLUS GRANDE DIMENSION = '&
     & ,1PE16.6)
      GOTO 807
  806 CONTINUE

      N=51
      IF(AKK.GT.40.)THEN
!      WRITE(*,822)
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
 !     WRITE(*,3002)N,AKK
 3002 FORMAT(2X,'SIMPSON ',I3,' POINTS'&
     & /2X,'K0 ADIMENSIONNEL PAR RAPPORT A LA PLUS GRANDE DIMENSION = '&
     &,1PE16.6)
  807 CONTINUE
      RETURN
      END SUBROUTINE
!--------------------------------------------------------------------------
                                                                  
      SUBROUTINE LISC(AK0,AM0,NM)
                                     
      INTEGER::NM,I,J,NJ,NPP
      REAL:: AK0,AM0                                                              
      REAL:: POL(31),A,B,H                               
      REAL:: S(4*(31-1),31+1),XT(4*(31-1)+1),YT(4*(31-1)+1)            
      REAL:: SC(31),VR(31),VC(31)
      INTEGER::ISTIR,NMAX,NK,ISOR,NPI,NMO
      REAL:: PRECI,ERMAX,ERMOY,XX,YY,TT,DIF,RT,NEXR 
      COMPLEX:: COM(31)

      S=0.
      SC=0.
      AR=0.
      NEXR=31                                           
      PRECI=1.E-02                                                              
      ISTIR=0                                                                   
      NMAX=4*(NEXR-1)                                                           
      NK=4                                                                      
      A=-0.1                                                                    
      B=20.                                                                     
   62 CONTINUE                                                                  
      NM=NK                                                                     
      NJ=4*NM                                                                                                                                   
      NPP=NJ+1                                                                   
      H=(B-A)/NJ                                                                
      DO 10 I=1,NPP                                                              
      XT(I)=A+(I-1)*H                                                            
      YT(I)=FF(XT(I),AK0,AM0)                                                 
   10 CONTINUE                                                                  
      ISOR=0                                                                         
      CALL EXPORS(XT,YT,NJ,NM,AMBDA,NMAX,S,SC,VR,VC,COM,POL,AR) 
    
      NPI=2                                                                     
      NMO=NPI*NPP-NPI+1                                                          
      ERMAX=0.                                                                  
      ERMOY=0.                                                                  
      DO 20 I=1,NMO                                                             
      XX=(I-1)*B/(NMO-1)                                                        
      YY=FF(XX,AK0,AM0)                                                         
      TT=0.                                                                      
      DO 30 J=1,NM                                                              
      RT=AMBDA(J)*XX                                                            
      IF(RT.GT.-20.)THEN
      TT=TT+AR(J)*EXP(RT)                                                         
      ENDIF
   30 CONTINUE                                                                  
      DIF=YY-TT                                                                  
      ERMOY=ERMOY+DIF                                                           
      ERMAX=AMAX1(ERMAX,ABS(DIF))                                               
      IF(ABS(DIF).GT.PRECI)ISOR=1                                               
   20 CONTINUE                                                                  
      ERMOY=ERMOY/NMO                                                        
 !     WRITE(*,1111)NM,ERMAX,ERMOY                                              
 1111 FORMAT(5X,I2,'EXPONENTIELLES  ECART MAXI = ',E10.3,'ECART MOYEN = ',E10.3/)                                                              
      IF(ISTIR.EQ.1)GOTO 61                                                     
      IF(ISOR)63,61,63                                                          
   63 CONTINUE                                                                  
      NK=NK+2                                                                   
      IF(NK-(NEXR-1))62,62,65                                                   
   65 CONTINUE !WRITE(*,6500)PRECI,NM                                                    
 6500 FORMAT(/5X,'PRECISION = ',E10.3,'  NON ATTEINTE AVEC ',I2,'  EXPONENTIELLES')                                                               
! C      STOP                                                                      
   61 DO 60 J=1,NM                                                              
!      WRITE(*,1100)AR(J),AMBDA(J)                                              
 1100 FORMAT(5X,E16.7,'EXP(',E16.7,')')                                         
      IF(AMBDA(J).GT.0.)STOP                                                    
   60 CONTINUE                                                                  
      RETURN
                                                                    
      END SUBROUTINE            
!-------------------------------------------------------------------------------!                                                                        

      SUBROUTINE EXPORS(XT,YT,NJ,NM,VCOM,NMAX,S,SC,VR,VC,COM,POL,AR)           
      INTEGER::NJ,NM,NMAX
      REAL:: VCOM(31),POL(31),AR(31),SC(31),VR(31),VC(31)
      REAL:: S(4*(31-1),31+1),XT(4*(31-1)+1),YT(4*(31-1)+1)           
      COMPLEX:: COM(31) 
      INTEGER::I,J,K,NPP,JJ,II,IJ,MN,NEXP
      INTEGER::IS,IER
      REAL::H,EPS                                                       
                                                     
      NPP=NJ+1                                                                   
      H=(XT(NPP)-XT(1))/NJ                                                         
      K=NPP-NM                                                                  
      DO 2 I=1,K                                                                
      DO 1 J=1,NM                                                               
      JJ=NM-J+I                                                     
 1    S(I,J)=YT(JJ)                                                              
      II=NM+I                                                         
 2    S(I,NM+1)=-YT(II)                                              
      EPS=1.E-20
                                                                
      CALL HOUSRS(S,NMAX,K,NM,1,EPS)          
      DO 5 I=1,NM                                                               
      IJ=NM-I+1                                                                 
 5    SC(IJ)=S(I,NM+1)                                                          
      MN=NM+1                                                                   
      SC(MN)=1.                                                                 
      CALL SPRBM(SC,MN,VR,VC,POL,IS,IER)                              
      DO 6 I=1,NM                                                               
      COM(I)=CMPLX(VR(I),VC(I))                                             
      COM(I)=CLOG(COM(I))/H                                                     
      VR(I)=REAL(COM(I))                                                  
      VC(I)=AIMAG(COM(I))                                  
    6 CONTINUE                                                                  
      I=1                                                                       
      J=0                                                                       
  100 IF(VC(I))110,111,110                                                      
  111 J=J+1                                                                     
      VCOM(J)=VR(I)                                                             
      I=I+1                                                                     
      GO TO 101                                                                 
  110 IF(ABS(VR(I)-VR(I+1))-1.E-5)120,120,121                                   
  120 J=J+1                                                                     
      VCOM(J)=VR(I)                                                             
      I=I+2                                                                     
      GO TO 101                                                                 
  121 J=J+1                                                                     
      VCOM(J)=VR(I)                                                             
      I=I+1                                                                     
  101 IF(I-NM)100,100,102                                                       
  102 NEXP=J                                                                    
      J=0                                                                       
      DO 300 I=1,NEXP                                                           
      J=J+1                                                                     
      IF(VCOM(I).GE.0.)GOTO 301                                                 
      IF(VCOM(I)+20.)301,301,302                                                
  301 J=J-1                                                                     
      GO TO 300                                                                 
  302 VCOM(J)=VCOM(I)                                                           
  300 CONTINUE                                                                  
      NEXP=J                                                                    
      NM=NEXP                                                                 
      CALL MCAS(NEXP,VCOM,XT,YT,NPP,AR,S,NMAX)                                 
      RETURN                                                                    
      END SUBROUTINE
!----------------------------------------------------------------------------

      SUBROUTINE MCAS(NEXP,TEXP,XT,YT,NPP,AR,A,NMAX)                            
      INTEGER:: NEXP,NPP,NMAX
      REAL::XT(4*(31-1)+1),YT(4*(31-1)+1),A(4*(31-1),31+1),AR(31),TEXP(31)
      INTEGER::I,J,L,M,N
      REAL::S,TT,TTT,EPS
                        
      EPS=1.E-20                                                                
      DO 1 I=1,NEXP                                                             
      DO 1 J=1,NEXP                                                             
      S=0                                                                       
      DO 3 L=1,NPP                                                               
      TT=(TEXP(I)+TEXP(J))*XT(L)                                                 
      IF(TT+30)3,4,4                                                            
    4 S=S+EXP(TT)                                                               
    3 CONTINUE                                                                  
      A(I,J)=S                                                                  
    1 CONTINUE                                                                  
      DO 5 I=1,NEXP                                                             
      S=0                                                                       
      DO 6 L=1,NPP                                                               
      TTT=TEXP(I)*XT(L)                                                          
      IF(TTT+30)6,7,7                                                           
    7 S=S+EXP(TTT)*YT(L)                                                         
    6 CONTINUE                                                                  
      A(I,NEXP+1)=S                                                             
    5 CONTINUE                                                                  
      N=NEXP                                                                    
      M=N+1                                                                  
      CALL HOUSRS(A,NMAX,N,N,1,EPS)                                       
      DO 10 I=1,NEXP                                                            
   10 AR(I)=A(I,NEXP+1)                                                         
      RETURN
                                                                    
      END SUBROUTINE
!----------------------------------------------------------------------

      SUBROUTINE SPRBM(C,IC,RR,RC,POL,IR,IER)                                   
      INTEGER::IC,IR,IER
      REAL:: C(31),RR(31),RC(31),POL(31)                                   
      INTEGER::I,J,L,N,LIM,IST
      REAL::A,B,H
      REAL::EPS,Q1,Q2,Q(4) 

      EPS=1.E-6                                                                 
      LIM=100                                                                   
      IR=IC+1                                                                   
    1 IR=IR-1                                                                   
      IF(IR-1)42,42,2                                                           
    2 IF(C(IR))3,1,3                                                            
    3 IER=0                                                                     
      J=IR                                                                      
      L=0                                                                       
      A=C(IR)                                                                   
      DO 8 I=1,IR                                                               
      IF(L)4,4,7                                                                
    4 IF(C(I))6,5,6                                                             
    5 RR(I)=0.                                                                  
      RC(I)=0.                                                                  
      POL(J)=0.                                                                 
      J=J-1                                                                     
      GO TO 8                                                                   
    6 L=1                                                                       
      IST=I                                                                     
      J=0                                                                       
    7 J=J+1                                                                     
      C(I)=C(I)/A                                                               
      POL(J)=C(I)                                                               
      IF(ABS(POL(J))-1.E27)8,42,42                                              
    8 CONTINUE                                                                  
      Q1=0.                                                                     
      Q2=0.                                                                     
    9 IF(J-2)33,10,14                                                           
   10 A=POL(1)                                                                  
      RR(IST)=-A                                                                
      RC(IST)=0.                                                                
      IR=IR-1                                                                   
      Q2=0.                                                                     
      IF(IR-1)13,13,11                                                          
   11 DO 12 I=2,IR                                                              
      Q1=Q2                                                                     
      Q2=POL(I+1)                                                               
   12 POL(I)=A*Q2+Q1                                                            
   13 POL(IR+1)=A+Q2                                                            
      GO TO 34                                                                  
   14 DO 22 L=1,10                                                              
      N=1                                                                       
   15 Q(1)=Q1                                                                   
      Q(2)=Q2                                                                   
      CALL SPQFB(POL,J,Q,LIM,I)                                                 
      IF(I)16,24,23                                                             
   16 IF(Q1)18,17,18                                                            
   17 IF(Q2)18,21,18                                                            
   18 GOTO(19,20,19,21),N                                                       
   19 Q1=-Q1                                                                    
      N=N+1                                                                     
      GO TO 15                                                                  
   20 Q2=-Q2                                                                    
      N=N+1                                                                     
      GO TO 15                                                                  
   21 Q1=1.+Q1                                                                  
   22 Q2=1.-Q2                                                                  
      IER=3                                                                     
      IR=IR-J                                                                   
      GOTO 45                                                                   
   23 IER=1                                                                     
   24 Q1=Q(1)                                                                   
      Q2=Q(2)                                                                   
      B=0.                                                                      
      A=0.                                                                      
      I=J                                                                       
   25 H=-Q1*B-Q2*A+POL(I)                                                       
      POL(I)=B                                                                  
      B=A                                                                       
      A=H                                                                       
      I=I-1                                                                     
      IF(I-2)26,26,25                                                           
   26 POL(2)=B                                                                  
      POL(1)=A                                                                  
      L=IR-1                                                                    
      IF(J-L)27,27,29                                                           
   27 DO 28 I=J,L                                                               
   28 POL(I-1)=POL(I-1)+POL(I)*Q2+POL(I+1)*Q1                                   
   29 POL(L)=POL(L)+POL(L+1)*Q2+Q1                                              
      POL(IR)=POL(IR)+Q2                                                        
      H=-.5*Q2                                                                  
      A=H*H-Q1                                                                  
      B=SQRT(ABS(A))                                                            
      IF(A)30,30,31                                                             
   30 RR(IST)=H                                                                 
      RC(IST)=B                                                                 
      IST=IST+1                                                                 
      RR(IST)=H                                                                 
      RC(IST)=-B                                                                
      GO TO 32                                                                  
   31 B=H+SIGN(B,H)                                                             
      RR(IST)=Q1/B                                                              
      RC(IST)=0.                                                                
      IST=IST+1                                                                 
      RR(IST)=B                                                                 
      RC(IST)=0.                                                                
   32 IST=IST+1                                                                 
      J=J-2                                                                     
      GO TO 9                                                                   
   33 IR=IR-1                                                                   
   34 A=0.                                                                      
      DO 38 I=1,IR                                                              
      Q1=C(I)                                                                   
      Q2=POL(I+1)                                                               
      POL(I)=Q2                                                                 
      IF(Q1)35,36,35                                                            
   35 Q2=(Q1-Q2)/Q1                                                             
   36 Q2=ABS(Q2)                                                                
      IF(Q2-A)38,38,37                                                          
   37 A=Q2                                                                      
   38 CONTINUE                                                                  
      I=IR+1                                                                    
      POL(I)=1.                                                                 
      RR(I)=A                                                                   
      RC(I)=0.                                                                  
      IF(IER)39,39,41                                                           
   39 IF(A-EPS)41,41,40                                                         
   40 IER=-1                                                                    
   41 GOTO 45                                                                   
   42 IER=2                                                                     
      IR=0                                                                      
   45 IF(IER-2)46,47,46                                                         
   47 WRITE(*,48)IER                                                           
   48 FORMAT(/5X,'IER = ',I3,'  ERREUR DANS SPRBM'/)                            
      STOP                                                                      
   46 RETURN  
                                                                  
      END SUBROUTINE
!----------------------------------------------------------------

      SUBROUTINE SPQFB(C,IC,Q,LIM,IER)                                          
      INTEGER::IC,LIM,IER,I,J,L,LL
      REAL:: C(31),Q(4)    
      REAL:: H,HH,A,A1,AA,B,BB,B1,C1,CA,CB,CC,CD,DQ1,DQ2,EPS,EPS1                                                     
      REAL:: Q1,Q2,QQ1,QQ2,QQQ1,QQQ2  
      
!--------- Value non initialized in previous versions ?!
      H=0.
      HH=0.
!---------      
      IER=0                                                                     
      J=IC+1                                                                    
    1 J=J-1                                                                     
      IF(J-1) 40,40,2                                                           
    2 IF(C(J)) 3,1,3                                                            
    3 A=C(J)                                                                    
      IF(A-1.) 4,6,4                                                            
    4 DO 5 I=1,J                                                                
      C(I)=C(I)/A                                                               
      IF(ABS(C(I))-1.E27)5,40,40                                                
    5 CONTINUE                                                                  
    6 IF(J-3) 41,38,7                                                           
    7 EPS=1.E-14                                                                
      EPS1=1.E-6                                                                
      L=0                                                                       
      LL=0                                                                      
      Q1=Q(1)                                                                   
      Q2=Q(2)                                                                   
      QQ1=0.                                                                    
      QQ2=0.                                                                    
      AA=C(1)                                                                   
      BB=C(2)                                                                   
      CB=ABS(AA)                                                                
      CA=ABS(BB)                                                                
      IF(CB-CA) 8,9,10                                                          
    8 CC=CB+CB                                                                  
      CB=CB/CA                                                                  
      CA=1.                                                                     
      GO TO 11                                                                  
    9 CC=CA+CA                                                                  
      CA=1.                                                                     
      CB=1.                                                                     
      GO TO 11                                                                  
   10 CC=CA+CA                                                                  
      CA=CA/CB                                                                  
      CB=1.                                                                     
   11 CD=CC*.1                                                                  
   12 A=0.                                                                      
      B=A                                                                       
      A1=A                                                                      
      B1=A                                                                      
      I=J                                                                       
      QQQ1=Q1                                                                   
      QQQ2=Q2                                                                   
      DQ1=HH                                                                    
      DQ2=H                                                                     
   13 H=-Q1*B-Q2*A+C(I)                                                         
      IF(ABS(H)-1.E27)14,42,42                                                  
   14 B=A                                                                       
      A=H                                                                       
      I=I-1                                                                     
      IF(I-1) 18,15,16                                                          
   15 H=0.                                                                      
   16 H=-Q1*B1-Q2*A1+H                                                          
      IF(ABS(H)-1.E27)17,42,42                                                  
   17 C1=B1                                                                     
      B1=A1                                                                     
      A1=H                                                                      
      GO TO 13                                                                  
   18 H=CA*ABS(A)+CB*ABS(B)                                                     
      IF(LL) 19,19,39                                                           
   19 L=L+1                                                                     
      IF(ABS(A)-EPS*ABS(C(1))) 20,20,21                                         
   20 IF(ABS(B)-EPS*ABS(C(2))) 39,39,21                                         
   21 IF(H-CC) 22,22,23                                                         
   22 AA=A                                                                      
      BB=B                                                                      
      CC=H                                                                      
      QQ1=Q1                                                                    
      QQ2=Q2                                                                    
   23 IF(L-LIM) 28,28,24                                                        
   24 IF(H-CD) 43,43,25                                                         
   25 IF(Q(1)) 27,26,27                                                         
   26 IF(Q(2)) 27,42,27                                                         
   27 Q(1)=0.                                                                   
      Q(2)=0.                                                                   
      GO TO 7                                                                   
   28 HH=AMAX1(ABS(A1),ABS(B1),ABS(C1))                                         
      IF(HH) 42,42,29                                                           
   29 A1=A1/HH                                                                  
      B1=B1/HH                                                                  
      C1=C1/HH                                                                  
      H=A1*C1-B1*B1                                                             
      IF(H) 30,42,30                                                            
   30 A=A/HH                                                                    
      B=B/HH                                                                    
      HH=(B*A1-A*B1)/H                                                          
      H=(A*C1-B*B1)/H                                                           
      Q1=Q1+HH                                                                  
      Q2=Q2+H                                                                   
      IF(ABS(HH)-EPS*ABS(Q1)) 31,31,33                                          
   31 IF(ABS(H)-EPS*ABS(Q2)) 32,32,33                                           
   32 LL=1                                                                      
      GO TO 12                                                                  
   33 IF(L-1)12,12,34                                                           
   34 IF(ABS(HH)-EPS1*ABS(Q1)) 35,35,12                                         
   35 IF(ABS(H)-EPS1*ABS(Q2)) 36,36,12                                          
   36 IF(ABS(QQQ1*HH)-ABS(Q1*DQ1)) 37,44,44                                     
   37 IF(ABS(QQQ2*H)-ABS(Q2*DQ2)) 12,44,44                                      
   38 Q(1)=C(1)                                                                 
      Q(2)=C(2)                                                                 
      Q(3)=0.                                                                   
      Q(4)=0.                                                                   
      GOTO 45                                                                   
   39 Q(1)=Q1                                                                   
      Q(2)=Q2                                                                   
      Q(3)=A                                                                    
      Q(4)=B                                                                    
      GOTO 45                                                                   
   40 IER=-1                                                                    
      GOTO 45                                                                   
   41 IER=-2                                                                    
      GOTO 45                                                                   
   42 IER=-3                                                                    
      GO TO 44                                                                  
   43 IER=1                                                                     
   44 Q(1)=QQ1                                                                  
      Q(2)=QQ2                                                                  
      Q(3)=AA                                                                   
      Q(4)=BB                                                                   
   45 RETURN                                                                    
      END SUBROUTINE
!---------------------------------------------------------------------

      SUBROUTINE HOUSRS(A,NMAX,NL,NCC,NS,EPS)
      INTEGER::NMAX,NL,NCC,NS
      REAL:: A(NMAX,31+1),EPS
      INTEGER::I,J,K,L,M,NCJ,NTC,KP1
      REAL::E,E0,AR,BA,ETA
      INTEGER :: I1

      NTC=NCC+NS                                                                 
      IF(NCC.GT.NL)THEN                                                          
	WRITE(*,3010)                                                            
  3010 FORMAT(' NBRE DE COLONNES > NBRES DE LIGNES')                             
	STOP                                                             
      ENDIF
      DO 13 K=1,NCC                                                              
	E=0                                                                       
	DO 1101 I=K,NL                                                            
	  E=E+A(I,K)**2                                                            
   1101 CONTINUE                                                                  
	E0=SQRT(E)             
	IF(E0.LT.EPS)THEN
	  WRITE(*,201)EPS
      201 FORMAT(1X,'NORME INFERIEURE A ',1PE16.6/)
	  STOP                                                                    
	ENDIF
	IF(A(K,K).EQ.0)THEN                                                      
	  AR=-E0                                                                    
	ELSE                                                                   
	  AR=-SIGN(E0,A(K,K))
	ENDIF
	ETA=AR*(AR-A(K,K))                                                        
	KP1=K+1                                                                   
	DO 10 J=KP1,NTC                                                           
	  BA=(A(K,K)-AR)*A(K,J)                                                     
	  DO 9 I=KP1,NL                                                            
	    BA=BA+A(I,K)*A(I,J)   
	9 CONTINUE                                                   
	  A(K,J)=A(K,J)+BA/AR
	  DO 11 I=KP1,NL                                                            
	    A(I,J)=A(I,J)-A(I,K)*BA/ETA                                               
       11 CONTINUE
     10 CONTINUE                                                                  
	A(K,K)=AR                                                                 
	DO 12 I=KP1,NL                                                            
  12   A(I,K)=0                                                                  
 13   CONTINUE                                                                  
      DO 1006 J=1,NS                                                            
	NCJ=NCC+J                                                                  
	A(NCC,NCJ)=A(NCC,NCJ)/A(NCC,NCC)                                              
	DO 1005 L=2,NCC                                                            
	  I1=NCC+1-L                                                                 
	  M=I1+1                                                                    
	  DO 1004 I=M,NCC                                                            
      1004 A(I1,NCJ)=A(I1,NCJ)-A(I1,I)*A(I,NCJ)                                      
    1005 A(I1,NCJ)=A(I1,NCJ)/A(I1,I1)                                              
 1006 CONTINUE                                                                  
      RETURN

      END SUBROUTINE
!-----------------------------------------------------------

      SUBROUTINE VSD(XFT,YFT,ZFT,JJ,XNF,YNF,ZNF,A,TD,XGJ,YGJ,ZGJ, &
                   & XGI,YGI,ZIG,FS,VX,VY,VZ)                                              
      INTEGER::JJ
      REAL::XGJ,YGJ,ZGJ,XGI,YGI,ZIG,FS,VX,VY,VZ
      REAL:: A,TD,XNF,YNF,ZNF,XFT(5),YFT(5),ZFT(5)
      REAL:: QPI                    
      INTEGER::MJJ,L
      REAL:: RR(5),DRX(5),DRY(5),DRZ(5)
      REAL:: AIJS,ALDEN,ANL,ANLX,ANLY,ANLZ,ANT,ANTY,ANTX,ANTZ,ARG,ASRO
      REAL:: AT,ATX,ATY,ATZ,DAT,DDK,DEN,DENL,DENT,DK,DLG
      REAL::DNL,DNT,DNTX,DNTY,DNTZ,DR,DS,EPS,GY,GYX,GYZ,GYY,GZ,PJ,QJ 
      REAL::QMP,RJ,RO,SGN,VXS,VZS,VYS


      QPI=4.*4.*ATAN(1.)
      MJJ=(-1)**(JJ+1)
      QMP=MJJ/(2*QPI)
      EPS=1.E-4
      RO=SQRT((XGI-XGJ)**2+(YGI-YGJ*MJJ)**2+(ZIG-ZGJ)**2)
      GZ=(XGI-XGJ)*XNF+(YGI-YGJ*MJJ)*YNF*MJJ+(ZIG-ZGJ)*ZNF
      IF(RO.GT.7*TD)THEN
      FS=-A/(RO*2*QPI)
      ASRO=FS/RO**2                                                       
      VX=-(XGI-XGJ)*ASRO                                              
      VY=-(YGI-YGJ*MJJ)*ASRO                                              
      VZ=-(ZIG-ZGJ)*ASRO                                              
      ELSE
	DO 212 L=1,4
      RR(L)=SQRT((XGI-XFT(L))**2+(YGI-YFT(L)*MJJ)**2+(ZIG-ZFT(L))**2)
      DRX(L)=(XGI-XFT(L))/RR(L)                                                
      DRY(L)=(YGI-YFT(L)*MJJ)/RR(L)                                                
      DRZ(L)=(ZIG-ZFT(L))/RR(L)                                                
  212 CONTINUE
      RR(5)=RR(1)
      DRX(5)=DRX(1)                                                           
      DRY(5)=DRY(1)                                                           
      DRZ(5)=DRZ(1)                                                           
      AIJS=0.
      VXS=0.
      VYS=0.
      VZS=0.
      DO 29 L=1,4
      DK=SQRT((XFT(L+1)-XFT(L))**2+(YFT(L+1)-YFT(L))**2+(ZFT(L+1)-ZFT(L))**2)
      IF(DK.GE.1.E-3*TD)THEN                                             
      PJ=(XFT(L+1)-XFT(L))/DK                                                     
      QJ=(YFT(L+1)-YFT(L))*MJJ/DK                                                     
      RJ=(ZFT(L+1)-ZFT(L))/DK                                                     
      GYX=YNF*MJJ*RJ-ZNF*QJ                                                     
      GYY=ZNF*PJ-XNF*RJ                                                     
      GYZ=XNF*QJ-YNF*MJJ*PJ                                                     
      GY=(XGI-XFT(L))*GYX+(YGI-YFT(L)*MJJ)*GYY+(ZIG-ZFT(L))*GYZ
      SGN=SIGN(1.,GZ)                                                           
      DDK=2.*DK                                                                 
      ANT=GY*DDK                                                                
      DNT=(RR(L+1)+RR(L))**2-DK*DK+2.*ABS(GZ)*(RR(L+1)+RR(L))
      ANL=RR(L+1)+RR(L)+DK                                                      
      DNL=RR(L+1)+RR(L)-DK                                                      
      DEN=ANL/DNL                                                               
      ALDEN=ALOG(DEN)                                                           
      IF(ABS(GZ).GT.1.E-4*TD)THEN
      ARG=ANT/DNT
      AT=ATAN(ARG)
      ELSE
      AT=0.                                                                     
      ENDIF                                                                  
      AIJS=AIJS+GY*ALDEN-2.*ABS(GZ)*AT
      DAT=2.*AT*SGN                                                             
      ANTX=GYX*DDK                                                              
      ANTY=GYY*DDK                                                              
      ANTZ=GYZ*DDK                                                              
      ANLX=DRX(L+1)+DRX(L)                                                      
      ANLY=DRY(L+1)+DRY(L)                                                      
      ANLZ=DRZ(L+1)+DRZ(L)                                                      
      DR=2.*(RR(L+1)+RR(L)+ABS(GZ))                                             
      DS=2.*(RR(L+1)+RR(L))*SGN                                                 
      DNTX=DR*ANLX+XNF*DS                                                     
      DNTY=DR*ANLY+YNF*MJJ*DS                                                     
      DNTZ=DR*ANLZ+ZNF*DS                                                     
      DENL=ANL*DNL                                                              
      DENT=ANT*ANT+DNT*DNT                                                      
      ATX=(ANTX*DNT-DNTX*ANT)/DENT                                              
      ATY=(ANTY*DNT-DNTY*ANT)/DENT                                              
      ATZ=(ANTZ*DNT-DNTZ*ANT)/DENT                                              
      DLG=(DNL-ANL)/DENL                                                       
      VXS=VXS+GYX*ALDEN+GY*ANLX*DLG-2.*ABS(GZ)*ATX-DAT*XNF
      VYS=VYS+GYY*ALDEN+GY*ANLY*DLG-2.*ABS(GZ)*ATY-DAT*YNF*MJJ             
      VZS=VZS+GYZ*ALDEN+GY*ANLZ*DLG-2.*ABS(GZ)*ATZ-DAT*ZNF   
      ENDIF        
   29 CONTINUE
      FS=-AIJS*QMP
      VX=-VXS*QMP
      VY=-VYS*QMP
      VZ=-VZS*QMP
      ENDIF
      RETURN
      END SUBROUTINE
!---------------------------------------------------------------------

      FUNCTION FF(XTT,AK,AM)
      REAL::XTT,AK,AM
      REAL::COEF,TOL,H,A,B,C,D,E,F,FF                                                                     
      COEF=(AM+AK)**2/(AM**2-AK**2+AK)                                        
      H=XTT-AM                                                                    
      TOL=AMAX1(0.1,0.1*AM)
      IF(ABS(H).GT.TOL)THEN                                                  
      FF=(XTT+AK)*EXP(XTT)/(XTT*SINH(XTT)-AK*COSH(XTT))-COEF/(XTT-AM)-2.                                                             
      ELSE                                                                    
      A=AM-TOL                                                                  
      B=AM                                                                      
      C=AM+TOL                                                                  
      D=(A+AK)*EXP(A)/(A*SINH(A)-AK*COSH(A))-COEF/(A-AM)-2                                                                   
      E=COEF/(AM+AK)*(AM+AK+1)-(COEF/(AM+AK))**2*AM-2                           
      F=(C+AK)*EXP(C)/(C*SINH(C)-AK*COSH(C))-COEF/(C-AM)-2                                                                   
      FF=(XTT-B)*(XTT-C)*D/((A-B)*(A-C))+(XTT-C)*(XTT-A)*E/((B-C)*(B-A))+&               
     & (XTT-A)*(XTT-B)*F/((C-A)*(C-B))                                           
      ENDIF
      RETURN
                                                                    
      END FUNCTION
!-----------------------------------------------------------------------------

      REAL FUNCTION X0(AK)
      REAL::AK
      INTEGER::ITOUR,IITER    
      REAL::XS,XI,XM,PAS,EPS,VAL1,VAL2
                                                 
      EPS=5.E-6                                                                 
      ITOUR=0                                                                   
      XI=0.                                                                     
      XS=XI                                                                     
      PAS=AMAX1(AK,SQRT(AK))                                                    
   30 XS=XS+PAS                                                                 
      ITOUR=ITOUR+1                                                             
      IF(ITOUR.LE.1000)THEN
        VAL1= AK-XS*TANH(XS)
        VAL2= AK-XI*TANH(XI)                                                  
        IF(VAL1*VAL2.GT.0) GOTO 30                                  
      ENDIF
      IITER=0                                                                   
   10 CONTINUE
      XM=(XI+XS)*0.5                                                            
      IITER=IITER+1                                                             
      IF(IITER.GT.1000.OR.ITOUR.GT.1000)THEN
	WRITE(*,110)ITOUR,IITER                                                  
    110 FORMAT(2X,'ERREUR DANS LA RECHERCHE DE LA RACINE', &                       
      & /2X,'APPROXIMATION =',I5,'   DICHOTOMIE = ',I5)                           
	STOP                                                                      
      ELSE
	IF(ABS((XS-XI)/XM).GT.EPS)THEN

          VAL1= AK-XM*TANH(XM)
          VAL2= AK-XI*TANH(XI)
        if(abs(val1*val2) .le. 1.E-12)val1=0.
          IF(VAL1*VAL2.LT.0)THEN                                                
	    XS=XM                                                                     
	  ELSE                                                                      
	    XI=XM                                                                     
	  ENDIF                                                     
	  GOTO 10
          
	ELSE                                 
	  X0=XM
	ENDIF
      ENDIF

      RETURN                                                                    
      END FUNCTION
!-----------------------------------------------------------------------------!

      REAL FUNCTION PL2(U1,U2,U3,XU)
	REAL::U1,U2,U3,XU
	PL2=((XU-U1)*(XU-U2))/((U3-U1)*(U3-U2))
	RETURN
      END FUNCTION
END MODULE