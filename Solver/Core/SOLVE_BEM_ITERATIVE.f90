MODULE SOLVE_BEM_ITERATIVE

  USE COM_VAR
  USE ELEMENTARY_FNS
  USE COMPUTE_GREEN_INFD
  USE COMPUTE_GREEN_FD
  USE M_SOLVER
  IMPLICIT NONE

  CONTAINS
!--------------------------------------------------------------------------!
  SUBROUTINE SOLVE_POTENTIAL_ITERATIVE(AM0,NVEL,LWORK,WORK)

    LOGICAL:: RHS
    INTEGER:: LWORK
    COMPLEX:: WORK(LWORK)
    INTEGER :: BX,ISYM,NJJ
    INTEGER ::I,J
    REAL :: G,W,tdepth,BETA
    REAL :: VSXP(NFA,NFA),VSYP(NFA,NFA),VSZP(NFA,NFA)
    REAL :: VSXM(NFA,NFA),VSYM(NFA,NFA),VSZM(NFA,NFA)
    REAL :: FSP(NFA,NFA),FSM(NFA,NFA)
    REAL:: CB,SB,CP,CR,COEFB,COEFS
    REAL:: AKAD,AM0,SD1B,SD1S,SD2B,SD2S,PI,DPI,ZERO
    COMPLEX:: ZOL(IMX,2),CZERO
!   Body condition    
    COMPLEX,DIMENSION(*) :: NVEL    
! Variables for GMRES
    INTEGER::revcom,colx, coly, colz, nbscal
    INTEGER:: irc(5), icntl(8), info(3)
    REAL:: cntl(5), rinfo(2)
      NJJ=NSYMY+1
      PI=4.*ATAN(1.)
      DPI=2.*PI
      G=9.81
      W=DPI/T
      ZERO=0.0
      CZERO=(0.0,0.0)
!--------------Initilizations---------------
      VSXP=0.
      VSXM=0.
      VSYP=0.
      VSYM=0.
      VSZP=0.
      VSZM=0.
      ZOL=CMPLX(0.,0.)
      ZIJ=CMPLX(0.,0.)
      WORK=CMPLX(0.,0.)

! --- Modification AB 170412 --- 
      AM0=W*W/G
      IF ((Depth .EQ. 0.).OR.(AM0*Depth.GE.20)) then
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

!Construction of the RHS vector
      DO 2000 ISYM=1,NJJ
        BX=(-1)**(ISYM+1)
        
! --- Modification AB 180412 --- 
	  DO I=1,IMX	    
	        IF (NSYMY.EQ.1) THEN
	            WORK(I+IMX)=(NVEL(I)+BX*NVEL(I+NFA))*0.5
	        ELSE
	            WORK(I+IMX)=NVEL(I)
	        END IF
	  ENDDO
! --- End of modification AB 180412 ---

! Initialize the control parameters for GMRES
! For details on GMRES parameter see at the end of this subroutine
       cntl(1) = TOL_GMRES !Tolerance
       cntl(2) = 0.0e-0
       cntl(3) = 0.0e-0
       cntl(4) = 0.0e-0
       cntl(5) = 0.0e-0
       icntl(1) = 6
       icntl(2) = 6
       icntl(3) = 20       !Save the convergence history in file fort.20
       icntl(4) = 0        !Preconditioning
       icntl(5) = 3        !ICGS orthogonalization
       icntl(6) = 0
       icntl(7) = MAXIT !Maximum number of iterations
       icntl(8) = 1

      10 call drive_cgmres(IMX,IMX,IRES,LWORK,WORK,IRC,ICNTL,CNTL,INFO,RINFO)
      revcom = irc(1)
      colx = irc(2)
      coly = irc(3)
      colz = irc(4)
      nbscal = irc(5)
!             print*,'check revcom',revcom
!       read*
      if (revcom.eq. 1) then
  !     * perform the matrix vector product work(colz) <-- A * work(colx)
        call matvec(ISYM,work(colx),work(colz))
	goto 10    
      else if (revcom.eq. 2) then
  !     * perform the left preconditioning work(colz) <-- M_1^{-1} * work(colx)
	call ccopy(IMX,work(colx),1,work(colz),1)
	goto 10
      else if (revcom.eq. 3) then
  !     * perform the right preconditioning work(colz) <-- M_2^{-1} * work(colx)
	call ccopy(IMX,work(colx),1,work(colz),1)
	goto 10      
      else if (revcom.eq. 4) then
	do i=0,nbscal-1
	  work(colz+i) = cdotc(IMX,work(colx+i*IMX),1,work(coly),1)
	enddo
	goto 10
      endif
      
      
! --- Modification AB 180412 ---   
     
      if (info(1).eq.0) then
            OPEN(13,FILE='gmres.log') 
            write(13,*) ' Normal exit '
            write(13,*) ' Convergence after ', info(2),' iterations'
            write(13,*) ' Backward error - preconditioned system ', rinfo(1)
            write(13,*) ' Backward error - unpreconditioned system', rinfo(2)
            write(13,*) ' Optimal size for workspace ', info(3)
            CLOSE(13)
      else if (info(1).eq.-1) then
            write(*,*) ' Bad value of n'
            STOP
      else if (info(1).eq.-2) then
            write(13,*) ' Bad value of m'
            STOP
      else if (info(1).eq.-3) then
            write(13,*) ' Too small workspace. '
            write(13,*) ' Minimal value should be ', info(2)
            STOP
      else if (info(1).eq.-4) then
            write(13,*) ' No convergence after ', icntl(7),' iterations'
            STOP
      else if (info(1).eq.-5) then
            write(13,*) ' Type of preconditioner not specified'
            STOP
      endif
! --- End modification AB 180412 --- 

	DO 2002 I=1,IMX
	  ZOL(I,(ISYM-1)+1)=work(I) 
!           print*,'sol',I,work(I) 
  2002 CONTINUE
!      stop
 2000 CONTINUE

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

  END SUBROUTINE
!-------------------------------------------------------------------------

  SUBROUTINE MATVEC(ISYM,XX,AXX)
  INTEGER:: ISYM
  COMPLEX:: XX(*),AXX(*)
  INTEGER:: I,J,BX
  REAL:: PCOS,PSIN

  do j=1,IMX
    AXX(j)=0.
  enddo

  BX=(-1)**(ISYM+1)
  DO 1 I=1,IMX
    IF(ISYM.EQ.1)THEN
      DO 31 J=1,IMX
	PCOS=VSXP1(I,J)*XN(I)+VSYP1(I,J)*YN(I)+VSZP1(I,J)*ZN(I)
	PSIN=VSXP2(I,J)*XN(I)+VSYP2(I,J)*YN(I)+VSZP2(I,J)*ZN(I)
	AXX(I)=AXX(I)+CMPLX(PCOS,PSIN)*XX(J)
      31 CONTINUE
    ELSE
      DO 32 J=1,IMX
	PCOS=VSXM1(I,J)*XN(I)+VSYM1(I,J)*YN(I)+VSZM1(I,J)*ZN(I)
	PSIN=VSXM2(I,J)*XN(I)+VSYM2(I,J)*YN(I)+VSZM2(I,J)*ZN(I)
	AXX(I)=AXX(I)+CMPLX(PCOS,PSIN)*XX(J)           
      32 CONTINUE
    ENDIF
  1 CONTINUE
  END SUBROUTINE
!----------------------------------------------------------------------

END MODULE

!*  some control parametrs definition for GMRES
!* icntl    (input) INTEGER array. length 6
!*            icntl(1) : stdout for error messages
!*            icntl(2) : stdout for warnings
!*            icntl(3) : stdout for convergence history
!*            icntl(4) : 0 - no preconditioning
!*                       1 - left preconditioning
!*                       2 - right preconditioning
!*                       3 - double side preconditioning
!*                       4 - error, default set in Init
!*            icntl(5) : 0 - modified Gram-Schmidt
!*                       1 - iterative modified Gram-Schmidt
!*                       2 - classical Gram-Schmidt
!*                       3 - iterative classical Gram-Schmidt
!*            icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
!*                       1 - user supplied initial guess
!*            icntl(7) : maximum number of iterations
!*            icntl(8) : 1 - default compute the true residual at each restart
!*                       0 - use recurence formaula at restart
!*
!* cntl     (input) real array, length 5
!*            cntl(1) : tolerance for convergence
!*            cntl(2) : scaling factor for normwise perturbation on A
!*            cntl(3) : scaling factor for normwise perturbation on b
!*            cntl(4) : scaling factor for normwise perturbation on the
!*                      preconditioned matrix
!*            cntl(5) : scaling factor for normwise perturbation on 
!*                      preconditioned right hand side