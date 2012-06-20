MODULE SOLVE_BEM_FD_ITERATIVE

  USE COM_VAR
  USE ELEMENTARY_FNS
  USE COMPUTE_GREEN_FD
  USE M_SOLVER
  IMPLICIT NONE

  CONTAINS
!--------------------------------------------------------------------------!
  SUBROUTINE SOLVE_POTENTIAL_FD_ITERATIVE(LWORK,WORK,NVEL,AMH,NEXP)

    COMPLEX,DIMENSION(*) :: NVEL
    INTEGER:: LWORK
    COMPLEX:: WORK(LWORK)
    INTEGER:: BX,ISYM,NJJ,N
    INTEGER:: I,J
    INTEGER:: NP1,I1,JJ,NEXP
    REAL:: W,BJ,GM,DIJ,AKK 
!     REAL:: CB,SB,CP,CR,COEFB,COEFS,BETA,AKAD,SD1B,SD1S,SD2B,SD2S,
    REAL:: AM0,AKH,AMH,PI,DPI,ZERO
    COMPLEX:: ZOL(IMX,2)
    
! Variables for GMRES
    INTEGER::revcom,colx, coly, colz, nbscal
    INTEGER:: irc(5), icntl(8), info(3)
    REAL:: cntl(5), rinfo(2)

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
 7106 WRITE(*,7102)                                                            
 7102 FORMAT(/5X,'PROFONDEUR QUASI-INFINIE'/5X, &                                
     & 'LE PROGRAMME EN PROFONDEUR INFINIE SERAIT PLUS ADAPTE')                  
      GOTO 7104                                                                 
 7101 IF(AKH-0.1)7103,7103,7104                                                 
 7103 WRITE(*,7105)                                                            
 7105 FORMAT(/5X,'PROFONDEUR TROP FAIBLE POUR LA LONGUEUR D''ONDE')             
 7104 CONTINUE 

!--------------Initilizations---------------
      CQ=0.0
      QQ=0.0
      AMBDA=0.0
      AR=0.0
      ZOL=CMPLX(0.,0.)
      ZIJ=CMPLX(0.,0.)
      WORK=CMPLX(0.,0.)
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
!-------------------------------------------------
!Construction of the RHS vector
DO 2000 ISYM=1,NJJ
      BX=(-1)**(ISYM+1)
	  DO I=1,IMX
	    IF (NSYMY.EQ.1) THEN
	      WORK(I+IMX)=(NVEL(I)+BX*NVEL(I+NFA))*0.5
	    ELSE
	      WORK(I+IMX)=NVEL(I)
	    END IF
	  ENDDO

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
        call matvec(ISYM,work(colx),work(colz),AM0,AMH,NEXP)
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
      
      if (info(1).eq.0) then
      write(*,*) ' Normal exit '
      write(*,*) ' Convergence after ', info(2),' iterations'
      write(*,*) ' Backward error - preconditioned system ', rinfo(1)
      write(*,*) ' Backward error - unpreconditioned system', rinfo(2)
      write(*,*) ' Optimal size for workspace ', info(3)
      else if (info(1).eq.-1) then
      write(*,*) ' Bad value of n'
      else if (info(1).eq.-2) then
      write(*,*) ' Bad value of m'
      else if (info(1).eq.-3) then
      write(*,*) ' Too small workspace. '
      write(*,*) ' Minimal value should be ', info(2)
      else if (info(1).eq.-4) then
      write(*,*) ' No convergence after ', icntl(7),' iterations'
      else if (info(1).eq.-5) then
      write(*,*) ' Type of preconditioner not specified'
      endif

	DO 2002 I=1,IMX
	  ZOL(I,(ISYM-1)+1)=work(I) 
!           print*,'sol',I,work(I) 
  2002 CONTINUE
!      stop
 2000 CONTINUE

	ZIGB=(0.,0.)
	ZIGS=(0.,0.)
	DO 2255 I=1,IMX
	  IF(NSYMY .EQ. 0)THEN
	    ZIGB(I)=ZOL(I,1)    ! Factor 2 is removed in comparison with previous version of Aquaplus because
	                        ! Normal velocity is halved in Aquaplus (because of symmetry)
	    ZIGS(I)=0.0
	  ELSE
	  ZIGB(I)=(ZOL(I,1)+ZOL(I,2))
	  ZIGS(I)=(ZOL(I,1)-ZOL(I,2))
          ENDIF
  2255 CONTINUE       

!computation of potential phi=S*sigma
      ZPB=(0.,0.)
      ZPS=(0.,0.)
      DO 2240 I=1,IMX
	  DO 2261 J=1,IMX
	      call VAVFD(2,XG(I),YG(I),ZG(I),I,J)

              call VNSFD(AM0,AMH,NEXP,I,J,XG(I),YG(I),ZG(I))  

	      ZPB(I)=ZPB(I)+0.5*(ZIGB(J)*CMPLX(SP1+SM1,SP2+SM2)&
                    +ZIGS(J)*CMPLX(SP1-SM1,SP2-SM2))
	      ZPS(I)=ZPS(I)+0.5*(ZIGS(J)*CMPLX(SP1+SM1,SP2+SM2)&
                    +ZIGB(J)*CMPLX(SP1-SM1,SP2-SM2))
    2261 CONTINUE
 2240 CONTINUE

   END SUBROUTINE
!-------------------------------------------------------------------------

  SUBROUTINE MATVEC(ISYM,XX,AXX,AM0,AMH,NEXP)
  INTEGER:: ISYM,I,NEXP
  COMPLEX:: XX(*),AXX(*)
  INTEGER:: ISP,IFP,BX
  REAL:: AM0,AMH,PCOS,PSIN

  do I=1,IMX
    AXX(I)=0.
  enddo

  BX=(-1)**(ISYM+1)
  DO 1 ISP=1,IMX
    IF(ISYM.EQ.1)THEN
      DO 31 IFP=1,IMX
	call VAVFD(2,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP)

	call VNSFD(AM0,AMH,NEXP,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP))  

	PCOS=VSXP1*XN(ISP)+VSYP1*YN(ISP)+VSZP1*ZN(ISP)
	PSIN=VSXP2*XN(ISP)+VSYP2*YN(ISP)+VSZP2*ZN(ISP)
	AXX(ISP)=AXX(ISP)+CMPLX(PCOS,PSIN)*XX(IFP)
      31 CONTINUE
    ELSE
      DO 32 IFP=1,IMX
	call VAVFD(2,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP)

	call VNSFD(AM0,AMH,NEXP,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP))  

	PCOS=VSXM1*XN(ISP)+VSYM1*YN(ISP)+VSZM1*ZN(ISP)
	PSIN=VSXM2*XN(ISP)+VSYM2*YN(ISP)+VSZM2*ZN(ISP)
	AXX(ISP)=AXX(ISP)+CMPLX(PCOS,PSIN)*XX(IFP)           
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