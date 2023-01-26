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
MODULE ELEMENTARY_FNS

  IMPLICIT NONE
 
  CONTAINS
!-----------------------------------------------------------------------

      COMPLEX FUNCTION GG(Z,CEX)                                              
      COMPLEX:: Z,CEX,Y
      REAL::T

      IF(REAL(Z)+16.)2,2,1       ! case 1 pp 368                                          
    2 Y=1./Z                                                                    
      GG=Y*(1.+Y*(-1.+Y*(2.+Y*(-6.+Y*(24.+Y*(-120.))))))                      
      RETURN                                                                    
    1 T=AIMAG(Z)                                                                
      IF(ABS(T)-10.)5,5,6          !case 3 pp 368                                             
    6 GG=0.711093/(Z+0.415775)+0.278518/(Z+2.29428)+0.010389/(Z+6.2900)
      IF(T)33,44,44                                                             
    5 CONTINUE
      IF(REAL(Z)+0.5)7,7,8  !case 2 pp 368
    8 GG=-(CLOG(Z)*(.1E+01+Z*(0.23721365E+00+Z*(0.206543E-01+Z*(0.763297E-03+&
         & Z*0.97087007E-05))))+0.5772156649E+00*(0.99999207E+00+&
         & Z*(-0.149545886E+01+Z*(0.41806426E-01+Z*(-0.3000591E-01+&
         & Z*(0.19387339E-02+Z*(-0.51801555E-03)))))))/(0.1E+01+&
         & Z*(-0.76273617E+00+Z*(0.28388363E+00+Z*(-0.66786033E-01+Z*(0.12982719E-01&
         & +Z*(-0.8700861E-03+Z*0.2989204E-03))))))                                                         
      IF(T)33,44,44   
   33 GG=GG-(0.,3.14159265)*CEX                                                 
      RETURN                                                                    
   44 GG=GG+(0.,3.14159265)*CEX                                                 
      RETURN
    7 CONTINUE

      IF(T)3,4,4   
!   Case 4 pp. 369                                                             
    3 GG=((((((( (1.000000,1.3935496E-06)*Z+ (15.82958,-20.14222))  &
     &*Z+ (-70.52863,-227.9511))*Z+ (-985.4221,-226.6272))*Z        &
     &+ (-1202.318,1580.907))*Z+ (953.2441,1342.447))*Z             &
     &+ (417.3716,-196.6665))*Z+ (-9.881266,-24.24952))/            &
     &(((((((( (1.000000,0.0000000E+00)*Z+ (16.83184,-20.14481))*Z  &
     &+ (-55.66969,-248.1167))*Z+ (-1068.640,-434.4707))*Z          &
     &+ (-2082.250,1522.471))*Z+ (383.3455,2730.378))*Z             &
     &+ (1216.791,351.7189))*Z+ (115.3926,-161.2647))*Z             &
     &+ (-3.777369,-4.510900))-(0.,3.14159265)*CEX  
      RETURN
!   Case 4 pp. 369
    4 GG=((((((( (1.000000,-1.3935496E-06)*Z+ (15.82958,20.14222))  &
     &*Z+ (-70.52863,227.9511))*Z+ (-985.4221,226.6272))*Z          &
     &+ (-1202.318,-1580.907))*Z+ (953.2441,-1342.447))*Z           &
     &+ (417.3716,196.6665))*Z+ (-9.881266,24.24952))/              &
     &(((((((( (1.000000,0.0000000E+00)*Z+ (16.83184,20.14481))*Z   &
     &+ (-55.66969,248.1167))*Z+ (-1068.640,434.4707))*Z            &
     &+ (-2082.250,-1522.471))*Z+ (383.3455,-2730.378))*Z           &
     &+ (1216.791,-351.7189))*Z+ (115.3926,161.2647))*Z             &
     &+ (-3.777369,4.510900))+(0.,3.14159265)*CEX
      RETURN
      END FUNCTION  
!-------------------------------------------------------------------------------!

     REAL FUNCTION CIH(AK,Z,H)
     REAL:: AK,Z,H

	IF((AK*H.LE.20).AND.(AK*H.GT.0.))THEN
	CIH=COSH(AK*(Z+H))/COSH(AK*H)
	ELSE
	CIH=EXP(AK*Z)
	ENDIF
      RETURN
      END FUNCTION
!-------------------------------------------------------------------------------!

      REAL FUNCTION SIH(AK,Z,H)
      REAL:: AK,Z,H

	IF((AK*H.LE.20).AND.(AK*H.GT.0.))THEN
	SIH=SINH(AK*(Z+H))/COSH(AK*H)
	ELSE
	SIH=EXP(AK*Z)
	ENDIF
      RETURN
      END FUNCTION
!-------------------------------------------------------------------------------!
END MODULE