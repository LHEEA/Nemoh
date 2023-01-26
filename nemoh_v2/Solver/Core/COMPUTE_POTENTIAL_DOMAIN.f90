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
SUBROUTINE COMPUTE_POTENTIAL_DOMAIN(ID,PHI,XC,YC,ZC,AM0,AMH,NEXP)
!
  USE MIDENTIFICATION
  USE COM_VAR
  USE ELEMENTARY_FNS
  USE COMPUTE_GREEN_FREESURFACE
!  
  IMPLICIT NONE
!
  TYPE(TID) :: ID
  INTEGER::I,J,NEXP
  REAL::XC,YC,ZC,AM0,AMH
  REAL:: S1B,S2B,S1S,S2S,C,D,POT1,POT2
  REAL:: W,PI,DPI
  COMPLEX :: PHI
!
  PI=4.*ATAN(1.)                                                          
  DPI=2.*PI 
  W=DPI/T  
  POT1=0.0
  POT2=0.0
  DO J=1,IMX
        SP1=0.
        SP2=0.
        SM1=0.
        SM2=0.
      IF((Depth .EQ. 0.) .OR. (AMH .GE. 20))  THEN         
	        CALL VVV(1,J,XC,YC,ZC) ! finite part
	        CALL VNV(J,XC,YC,ZC)  !infinite part
      ELSE
	        CALL VVV(2,J,XC,YC,ZC)
 	        CALL VNVF(AM0,AMH,NEXP,J,XC,YC,ZC)
      ENDIF
      S1B=REAL(ZIGB(J)) 
      S1S=REAL(ZIGS(J))
      S2B=AIMAG(ZIGB(J))
      S2S=AIMAG(ZIGS(J))
      IF(NSYMY.EQ.0)THEN
	      C=S1B*SP1
	      D=S2B*SP2
	      POT1=POT1+C-D
	      C=S1B*SP2
	      D=S2B*SP1
	      POT2=POT2+C+D
      ELSE
	      C=S1B*(SP1+SM1)+S1S*(SP1-SM1)
	      D=S2B*(SP2+SM2)+S2S*(SP2-SM2)
	      POT1=POT1+(C-D)*0.5
	      C=S1B*(SP2+SM2)+S1S*(SP2-SM2)
	      D=S2B*(SP1+SM1)+S2S*(SP1-SM1)
	      POT2=POT2+(C+D)*0.5
      ENDIF
  END DO 
  PHI=CMPLX(POT1,POT2)
!      PRE1=-POT2*W/G
!      PRE2=POT1*W/G
!      POMOD=SQRT(PRE1**2+PRE2**2)
!
END SUBROUTINE
