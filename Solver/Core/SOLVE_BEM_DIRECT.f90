MODULE SOLVE_BEM_DIRECT

  USE COM_VAR
  USE COMPUTE_GREEN_INFD
  USE COMPUTE_GREEN_FD
  USE ELEMENTARY_FNS
  USE M_SOLVER
  IMPLICIT NONE

  CONTAINS
!--------------------------------------------------------------------------!
  SUBROUTINE SOLVE_POTENTIAL_DIRECT(AM0,NVEL)
!In this subroutine the linear system Ax=b is constructed

    INTEGER :: BX,ISYM,NJJ
    INTEGER ::I,J
    REAL :: G,W,tdepth,BETA
    REAL :: VSXP(NFA,NFA),VSYP(NFA,NFA),VSZP(NFA,NFA)
    REAL :: VSXM(NFA,NFA),VSYM(NFA,NFA),VSZM(NFA,NFA)
    REAL :: FSP(NFA,NFA),FSM(NFA,NFA)
    REAL:: CB,SB,CP,CR,COEFB,COEFS,PCOS,PSIN
    REAL:: AKAD,AM0,SD1B,SD1S,SD2B,SD2S,PI,DPI,ZERO
    COMPLEX:: ZOL(IMX,2)
!   Body condition
    COMPLEX,DIMENSION(*) :: NVEL    

      NJJ=NSYMY+1
      PI=4.*ATAN(1.)
      DPI=2.*PI
      G=9.81
      W=DPI/T
      ZERO=0.0
!--------------Initilizations---------------
      VSXP=0.
      VSXM=0.
      VSYP=0.
      VSYM=0.
      VSZP=0.
      VSZM=0.
      ZOL=CMPLX(0.,0.)
      ZIJ=CMPLX(0.,0.)


! --- Modification AB 170412 --- 
      AM0=W*W/G
      IF ((Depth .EQ. 0.).OR.(AM0*Depth.GE.20)) then
!     IF (Depth .LE. 0.) THEN
        tdepth=1.0E+20
	    CALL VAVINFD(1,FSP,FSM,VSXP,VSXM,VSYP,VSYM,VSZP,VSZM) 
	    CALL VNSINFD(1,FSP,FSM,VSXP,VSXM,VSYP,VSYM,VSZP,VSZM) 
      ELSE
        AM0=0. ! computed through bisection method in NNSFD
        tdepth=depth
	    CALL VAVFD(2,FSP,FSM,VSXP,VSXM,VSYP,VSYM,VSZP,VSZM) 
	    CALL VNSFD(AM0,FSP,FSM,VSXP,VSXM,VSYP,VSYM,VSZP,VSZM) 
      END IF
! --- End of modification AB 170412 ---

!Construction of the influence matrix 
    DO ISYM=1,NJJ
        BX=(-1)**(ISYM+1)
	    DO I=1,IMX
	        IF(ISYM.EQ.1)THEN
	            DO J=1,IMX
	                PCOS=VSXP1(I,J)*XN(I)+VSYP1(I,J)*YN(I)+VSZP1(I,J)*ZN(I)
	                PSIN=VSXP2(I,J)*XN(I)+VSYP2(I,J)*YN(I)+VSZP2(I,J)*ZN(I)
	                ZIJ(I,J)=CMPLX(PCOS,PSIN)
	            END DO
	        ELSE
	            DO J=1,IMX
	                PCOS=VSXM1(I,J)*XN(I)+VSYM1(I,J)*YN(I)+VSZM1(I,J)*ZN(I)
	                PSIN=VSXM2(I,J)*XN(I)+VSYM2(I,J)*YN(I)+VSZM2(I,J)*ZN(I)
	                ZIJ(I,J)=CMPLX(PCOS,PSIN)           
                END DO
	        ENDIF
        END DO 

! --- Modification AB 180412 --- 
	    DO I=1,IMX	    
	        IF (NSYMY.EQ.1) THEN
	            ZIJ(I,IMX+1)=(NVEL(I)+BX*NVEL(I+NFA))*0.5
	        ELSE
	            ZIJ(I,IMX+1)=NVEL(I)
	        END IF
	    ENDDO
! --- End of modification AB 180412 ---
 
!------------------------------------------------!
        CALL GAUSSZ(ZIJ,NFA,IMX,IMX+1)
!------------------------------------------------!

	    DO I=1,IMX
	        ZOL(I,(ISYM-1)+1)=ZIJ(I,IMX+1)
        END DO
!           stop
    END DO

	ZIGB=(0.,0.)
	ZIGS=(0.,0.)
	DO 2255 I=1,IMX
	ZIGB(I)=(ZOL(I,1)+ZOL(I,2))
	ZIGS(I)=(ZOL(I,1)-ZOL(I,2))
  2255 CONTINUE      

!computation of potential phi=S*sigma
      ZPB=(0.,0.)
      ZPS=(0.,0.)
      DO 2240 I=1,IMX
	  DO 2261 J=1,IMX
	      ZPB(I)=ZPB(I)+0.5*(ZIGB(J)*CMPLX(SP1(I,J)+SM1(I,J),SP2(I,J)+ &
            & SM2(I,J))+ZIGS(J)*CMPLX(SP1(I,J)-SM1(I,J),SP2(I,J)-SM2(I,J)))
	      ZPS(I)=ZPS(I)+0.5*(ZIGS(J)*CMPLX(SP1(I,J)+SM1(I,J),SP2(I,J)+ &
            & SM2(I,J))+ZIGB(J)*CMPLX(SP1(I,J)-SM1(I,J),SP2(I,J)-SM2(I,J)))
    2261 CONTINUE
 2240 CONTINUE

     do i=1,IMX
!     print*,'ZPB,ZPS',I,ZPB(I),ZPS(I)
     enddo
!      read*
  END SUBROUTINE


END MODULE