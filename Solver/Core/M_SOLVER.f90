MODULE M_SOLVER

  IMPLICIT NONE
 
  CONTAINS
!---------------------------------------------------------------------------!
      SUBROUTINE GAUSSZ(A,NMAX,N,M)
        INTEGER:: N,M,NMAX,I,J,K,L,IL
	COMPLEX A(NMAX,1),C,P
        REAL:: EPS
        
	EPS=1.E-20
	DO 10 J=1,N-1
	K=J
	DO 20 I=J+1,N
	IF(ABS(A(K,J))-ABS(A(I,J)))30,20,20
    30 K=I
    20 CONTINUE
	IF(K-J)50,40,50
    50 DO 60 L=J,M
	C=A(J,L)
	A(J,L)=A(K,L)
    60 A(K,L)=C
    40 IF(ABS(A(J,J))-EPS)120,120,70
    70 DO 80 K=J+1,M
	P=A(J,K)/A(J,J)
	DO 80 I=J+1,N
    80 A(I,K)=A(I,K)-A(I,J)*P
    10 CONTINUE
	IF(ABS(A(N,N))-EPS)120,120,90
    90 DO 100 IL=N+1,M
	DO 100 J=N,1,-1
	A(J,IL)=A(J,IL)/A(J,J)
	DO 100 I=1,J-1
    100 A(I,IL)=A(I,IL)-A(I,J)*A(J,IL)
	RETURN
    120 WRITE(*,500)EPS
    500 FORMAT(5X,'PIVOT INFERIEUR A ',1P,E16.6)
	STOP
      END SUBROUTINE
!-------------------------------------------------------------------------

     subroutine drive_cgmres(n,nloc,m,lwork,work, &
     &                       irc,icntl,cntl,info,rinfo)
! *
! *  Purpose
! *  =======
! *    drive_cgmres is the driver routine for solving the linear system 
! *  Ax = b using the *  Generalized Minimal Residual iterative method 
! *  with preconditioning.
! *  This solver is implemented with a reverse communication scheme: control
! *  is returned to the user for computing the (preconditioned) 
! *  matrix-vector product.
! *  See the User's Guide for an example of use.
! *
! *
! * Written : June 1996
! * Authors : Luc Giraud, Serge Gratton, V. Fraysse
! *             Parallel Algorithms - CERFACS
! *
! * Updated : April 1997
! * Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
! *             Parallel Algorithms - CERFACS
! *
! * Updated : June 1998
! * Authors : Valerie Fraysse, Luc Giraud, Serge Gratton
! *             Parallel Algorithms - CERFACS
! * Purpose : Make clear that the warning and error messages come from the
! *           cgmres modules.
! *
! * Updated : December 2002 - L. Giraud, J.Langou
! * Purpose : Add the capability to avoid explicit residual calculation at restart
! *
! * Updated : March 2005 - L. Giraud
! * Purpose : small bug in the computation of the restart if too large vs the
! *           size of the workspace
! *
! *  Arguments
! *  =========
! *
! *  n      (input) INTEGER.
! *          On entry, the dimension of the problem.
! *          Unchanged on exit.
! *
! *  nloc   (input) INTEGER.
! *          On entry, the dimension of the local problem.
! *          In a parallel distributed envirionment, this corresponds
! *          to the size of the subset of entries of the right hand side
! *          and solution allocated to the calling process.
! *          Unchanged on exit.
! *
! *
! *  m      (input) INTEGER
! *          Restart parameter, <= N. This parameter controls the amount
! *          of memory required for matrix H (see WORK and H).
! *          Unchanged on exit.
! *
! *  lwork  (input) INTEGER
! *          size of the workspace
! *          if (icntl(5) = 0 or 1 )
! *            lwork >= m*m + m*(n+5) + 5*n+2, if icntl(8) = 1
! *            lwork >= m*m + m*(n+5) + 6*n+2, if icntl(8) = 0
! *          if (icntl(5) = 2 or 3 )
! *            lwork >= m*m + m*(n+5) + 5*n+m+1, if icntl(8) = 1
! *            lwork >= m*m + m*(n+5) + 6*n+m+1, if icntl(8) = 0
! *
! *  work   (workspace) real/complex array, length lwork
! *          work contains the required vector and matrices stored in the 
! *          following order :
! *            x  (n,1)       : computed solution.
! *            b  (n,1)       : right hand side.
! *            r0 (n,1)       : vector workspace.
! *            w  (n,1)       : vector workspace.
! *            V  (n,m)       : Krylov basis.
! *            H  (m+1,m+1)   : Hessenberg matrix (full storage).
! *            yCurrent (m,1) : solution of the current LS
! *            xCurrent (n,1) : current iterate
! *            rotSin (m,1)   : Sine of the Givens rotation
! *            rotCos (m,1)   : Cosine of the Givens rotation
! *
! *  irc     (input/output) INTEGER array. length 5
! *            irc(1) : REVCOM   used for reverse communication
! *                             (type of external operation)
! *            irc(2) : COLX     used for reverse communication
! *            irc(3) : COLY     used for reverse communication
! *            irc(4) : COLZ     used for reverse communication
! *            irc(5) : NBSCAL   used for reverse communication
! *
! *  icntl   (input) INTEGER array. length 7
! *            icntl(1) : stdout for error messages
! *            icntl(2) : stdout for warnings
! *            icntl(3) : stdout for convergence history
! *            icntl(4) : 0 - no preconditioning
! *                       1 - left preconditioning
! *                       2 - right preconditioning
! *                       3 - double side preconditioning
! *                       4 - error, default set in Init
! *            icntl(5) : 0 - modified Gram-Schmidt
! *                       1 - iterative modified Gram-Schmidt
! *                       2 - classical Gram-Schmidt
! *                       3 - iterative classical Gram-Schmidt
! *            icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
! *                       1 - user supplied initial guess
! *            icntl(7) : maximum number of iterations
! *            icntl(8) : 0 - use recurence formula at restart
! *                       1 - default compute the true residual at each restart
! *
! *  cntl    (input) real array, length 5
! *            cntl(1) : tolerance for convergence
! *            cntl(2) : scaling factor for normwise perturbation on A
! *            cntl(3) : scaling factor for normwise perturbation on b
! *            cntl(4) : scaling factor for normwise perturbation on the
! *                      preconditioned matrix
! *            cntl(5) : scaling factor for normwise perturbation on 
! *                      preconditioned right hand side
! *
! *  info    (output) INTEGER array, length 3
! *            info(1) :  0 - normal exit
! *                      -1 - n < 1
! *                      -2 - m < 1
! *                      -3 - lwork too small
! *                      -4 - convergence not achieved after icntl(7) iterations
! *                      -5 - precondition type not set by user
! *            info(2) : if info(1)=0 - number of iteration to converge
! *                      if info(1)=-3 - minimum workspace size necessary
! *            info(3) : optimal size for the workspace
! *
! *  rinfo   (output) real array, length 2
! *            if info(1)=0 
! *              rinfo(1) : backward error for the preconditioned system
! *              rinfo(2) : backward error for the unpreconditioned system
! *
! * Input variables
! * ---------------
       integer n, nloc, lwork, icntl(*)
       real   cntl(*)
       real   sA, sb, sPA, sPb
! * Output variables
! * ----------------
       integer  info(*)
       real    rinfo(*)
! * Input/Output variables
! * ----------------------
       integer  m, irc(*)
       complex work(*)
! * Local variables
! * ---------------
       integer xptr, bptr, wptr, r0ptr, Vptr, Hptr, dotptr
       integer yCurrent,rotSin, rotCos, xCurrent
       integer sizeWrk, newRestart
       integer iwarn, ierr, ihist, compRsd
       real    rn,rx,rc
       real DZRO
       parameter (DZRO = 0.0e0)

       integer icheck
       save icheck
       DATA icheck /0/

       intrinsic ifix, float
! *
! *       Executable statements :
! *
       ierr    = icntl(1)
       iwarn   = icntl(2)
       ihist   = icntl(3)
       compRsd = icntl(8)

       if (ierr.lt.0) ierr = 6

       if (compRsd.eq.1) then
          sizeWrk  = m*m + m*(nloc+5) + 5*nloc+ 1
       else
          sizeWrk  = m*m + m*(nloc+5) + 6*nloc+ 1
       endif

       if (icheck.eq.0) then
!* Check the value of the arguments
         if ((n.lt.1).or.(nloc.lt.1)) then
            write(ierr,*)
            write(ierr,*)' ERROR GMRES : '
            write(ierr,*)'     N < 1 '
            write(ierr,*)
            info(1) = -1
            irc(1)  = 0
            return
         endif
         if (m.lt.1) then
            write(ierr,*)
            write(ierr,*)' ERROR GMRES :'
            write(ierr,*)'     M < 1 '
            write(ierr,*)
            info(1) = -2
            irc(1)  = 0
            return
         endif
         if ((icntl(4).ne.0).and.(icntl(4).ne.1).and. &
     &     (icntl(4).ne.2).and.(icntl(4).ne.3)) then
            write(ierr,*)
            write(ierr,*)' ERROR GMRES : '
            write(ierr,*)'     Undefined preconditioner '
            write(ierr,*)
            info(1) = -5
            irc(1)  = 0
            return
         endif
!*
         if ((icntl(5).lt.0).or.(icntl(5).gt.3)) then
           icntl(5) = 0
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING  GMRES : '
             write(iwarn,*) '       Undefined orthogonalisation '  
             write(iwarn,*) '       Default MGS '
             write(iwarn,*)
           endif
         endif
!*
         if ((icntl(5).eq.2).or.(icntl(5).eq.3)) then
!* the workspace should be large enough to store the m dot-products
            sizeWrk  = sizeWrk  + m
         else
            sizeWrk  = sizeWrk  + 1
         endif
!*
         if (iwarn.ne.0) then
           write(iwarn,*)
           write(iwarn,*) ' WARNING GMRES : '
           write(iwarn,*) '       For M = ',m,' optimal value '     
           write(iwarn,*) '       for LWORK =  ', sizeWrk
           write(iwarn,*)
         endif
!*
         if ((icntl(6).ne.0).and.(icntl(6).ne.1)) then
           icntl(6) = 0
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING GMRES : '
             write(iwarn,*) '       Undefined intial guess '
             write(iwarn,*) '       Default x0 = 0 '
             write(iwarn,*)
           endif
         endif
         if (icntl(7).le.0) then
           icntl(7) = n
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING GMRES :'
             write(iwarn,*) '       Negative max number of iterations'  
             write(iwarn,*) '       Default N '
             write(iwarn,*)
           endif
         endif
         if ((icntl(8).ne.0).and.(icntl(8).ne.1)) then
           icntl(8) = 1
           write(iwarn,*)
           write(iwarn,*) ' WARNING GMRES :'
           write(iwarn,*) '       Undefined strategy for the residual'
           write(iwarn,*) '       at restart'
           write(iwarn,*) '       Default 1 '
           write(iwarn,*)
         endif
!* Check if the restart parameter is correct and if the size of the
!*  workspace is big enough for the restart.
!* If not try to fix correctly the parameters
!*
         if ((m .gt. n).or.(lwork.lt.sizeWrk)) then
           if (m .gt. n) then
             m = n
             if (iwarn.ne.0) then
               write(iwarn,*)
               write(iwarn,*) ' WARNING GMRES : '
               write(iwarn,*) '       Parameter M bigger than N'  
               write(iwarn,*) '       New value for M ',m
               write(iwarn,*)
             endif
             if (compRsd.eq.1) then
               sizeWrk = m*m + m*(nloc+5) + 5*nloc+1
             else
               sizeWrk = m*m + m*(nloc+5) + 6*nloc+1
             endif
             if ((icntl(5).eq.2).or.(icntl(5).eq.3)) then
!* the workspace should be large enough to store the m dot-products
                sizeWrk  = sizeWrk  + m
             else
                sizeWrk  = sizeWrk  + 1
             endif
           endif
           if ((lwork.lt.sizeWrk).and.(n.eq.nloc)) then
!* Compute the maximum size of the m according to the memory space
             rn         = float(n)
             rx         = rn + 5.0
             rc         = 5.0*rn + 1 - float(lwork)
!*
!* Update the linear part of the second order equation to be solved
             if ((icntl(5).eq.2).or.(icntl(5).eq.3)) then
               rx = rx + 1
             endif
!* Update the constant part of the second order equation to be solved
!*             
             if (icntl(8).eq.0) then
               rc = rc + rn
             endif
!* Solution of the the second order equation
             newRestart = ifix((-rx+sqrt(rx**2-4.0*rc))/2.0)
             if (newRestart.gt.0) then
               m = newRestart
               if (iwarn.ne.0) then
                 write(iwarn,*)
                 write(iwarn,*)' WARNING GMRES : '
                 write(iwarn,*)'       Workspace too small for M'  
                 write(iwarn,*)'       New value for M ',m
                 write(iwarn,*)
               endif
             else
               write(ierr,*)
               write(ierr,*)' ERROR GMRES : '
               write(ierr,*)'     Not enough space for the problem'
               write(ierr,*)'     the space does not permit any m'
               write(ierr,*)
               info(1) = -3
               irc(1)  = 0
               return
             endif
           endif
           if ((lwork.lt.sizeWrk).and.(n.ne.nloc)) then
              write(ierr,*)
              write(ierr,*)' ERROR GMRES : '
              write(ierr,*)'     Not enough space for the problem'
              write(ierr,*)
              info(1) = -3
              irc(1)  = 0
              return
           endif
         endif
!*
         info(3) = sizeWrk
         icheck = 1
!*
!* save the parameters the the history file
!*
         if (ihist.ne.0) then
           write(ihist,'(10x,A39)') 'CONVERGENCE HISTORY FOR GMRES'
           write(ihist,*)
           write(ihist,'(A30,I2)') 'Errors are displayed in unit: ',ierr
           if (iwarn.eq.0) then
             write(ihist,'(A27)') 'Warnings are not displayed:'
           else
             write(ihist,'(A32,I2)') 'Warnings are displayed in unit: ', &
     &                               iwarn
           endif 
           write(ihist,'(A13,I7)') 'Matrix size: ',n
           write(ihist,'(A19,I7)') 'Local matrix size: ',nloc
           write(ihist,'(A9,I7)') 'Restart: ',m
           if (icntl(4).eq.0) then
             write(ihist,'(A18)') 'No preconditioning'
           elseif (icntl(4).eq.1) then
             write(ihist,'(A20)') 'Left preconditioning'
           elseif (icntl(4).eq.2) then
             write(ihist,'(A21)') 'Right preconditioning'
           elseif (icntl(4).eq.3) then
             write(ihist,'(A30)') 'Left and right preconditioning'
           endif
           if (icntl(5).eq.0) then
             write(ihist,'(A21)') 'Modified Gram-Schmidt'
           elseif (icntl(5).eq.1) then
             write(ihist,'(A31)') 'Iterative modified Gram-Schmidt'
           elseif (icntl(5).eq.2) then
             write(ihist,'(A22)') 'Classical Gram-Schmidt'
           else
             write(ihist,'(A32)') 'Iterative classical Gram-Schmidt'
           endif
           if (icntl(6).eq.0) then
             write(ihist,'(A29)') 'Default initial guess x_0 = 0'
           else
             write(ihist,'(A27)') 'User supplied initial guess'
           endif
           if (icntl(8).eq.1) then
             write(ihist,'(A33)') 'True residual computed at restart'
           else
             write(ihist,'(A30)') 'Recurrence residual at restart'
           endif
           write(ihist,'(A30,I5)') 'Maximum number of iterations: ',&
     &                              icntl(7)
           write(ihist,'(A27,E11.4)') 'Tolerance for convergence: ', &
     &                                cntl(1) 
!* 
           write(ihist,'(A53)') &
     &       'Backward error on the unpreconditioned system Ax = b:'
           sA       = cntl(2)
           sb       = cntl(3)
           if ((sA.eq.DZRO).and.(sb.eq.DZRO)) then
             write(ihist,'(A39)') &
     &       '    the residual is normalised by ||b||'
           else
             write(ihist,1) sA,sb
1     format('    the residual is normalised by         ',E11.4, &
     &       ' * ||x|| + ',E11.4)
           endif
           sPA      = cntl(4)
           sPb      = cntl(5)
             write(ihist,2)
2     format('Backward error on the preconditioned system', &
     &        ' (P1)A(P2)y = (P1)b:')
           if ((sPA.eq.DZRO).and.(sPb.eq.DZRO)) then
             write(ihist,3)
3     format('    the preconditioned residual is normalised ', &
     &       'by ||(P1)b||')
           else
             write(ihist,4) sPA, sPb
4     format('    the preconditioned residual is normalised by ', E11.4, &
     &        ' * ||(P2)y|| + ',E11.4)
           endif
!*
           write(ihist,5) info(3)
5     format('Optimal size for the local workspace:',I7)
           write(ihist,*) 
           write(ihist,6)
6     format('Convergence history: b.e. on the preconditioned system')
           write(ihist,7)
7     format(' Iteration   Arnoldi b.e.    True b.e.')
         endif
!*
       endif
!* setup some pointers on the workspace
       xptr     = 1
       bptr     = xptr + nloc
       r0ptr    = bptr + nloc
       wptr     = r0ptr + nloc
       Vptr     = wptr + nloc
       if (compRsd.eq.1) then
          Hptr     = Vptr + m*nloc
       else
          Hptr     = Vptr + (m+1)*nloc
       endif
       dotptr   = Hptr + (m+1)*(m+1)
       if ((icntl(5).eq.2).or.(icntl(5).eq.3)) then
          yCurrent = dotptr + m
       else
         yCurrent = dotptr + 1
       endif
       xCurrent = yCurrent + m
       rotSin   = xCurrent + nloc
       rotCos   = rotSin + m
!*
       call cgmres(nloc,m,work(bptr),work(xptr),                      &
     &            work(Hptr),work(wptr),work(r0ptr),work(Vptr),       &
     &            work(dotptr),work(yCurrent),work(xCurrent),         &
     &            work(rotSin),work(rotCos),irc,icntl,cntl,info,rinfo)
!*
       if (irc(1).eq.0) then
         icheck = 0
       endif
!*
       return

       END SUBROUTINE
!----------------------------------------------------------------------------

        subroutine cgmres(n,m,b,x,H,w,r0,V,dot,yCurrent,xCurrent,rotSin,&
     &                   rotCos,irc,icntl,cntl,info,rinfo)
!*
!*
!*  Purpose
!*  =======
!*  cgmres solves the linear system Ax = b using the
!*  Generalized Minimal Residual iterative method
!*
!* When preconditioning is used we solve :
!*     M_1^{-1} A M_2^{-1} y = M_1^{-1} b
!*     x = M_2^{-1} y
!*
!*   Convergence test based on the normwise backward error for
!*  the preconditioned system
!*
!* Written : June 1996
!* Authors : Luc Giraud, Serge Gratton, V. Fraysse
!*             Parallel Algorithms - CERFACS
!*
!* Updated : April 1997
!* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
!*             Parallel Algorithms - CERFACS
!*
!* Updated : March 1998
!* Purpose : Pb with F90 on DEC ws
!*           cure : remove "ZDSCAL" when used to initialize vectors to zero
!*
!* Updated : May 1998
!* Purpose : r0(1) <-- r0'r0 : pb when used with DGEMV for the dot product
!*           cure : w(1) <--  r0'r0
!*
!* Updated : June 1998
!* Purpose : Make clear that the warning and error messages come from the
!*           cgmres modules.
!*
!* Updated : February 2001 - L. Giraud
!* Purpose : In complex version, initializations to zero performed  in complex
!*           arithmetic to avoid implicit conversion by the compiler.
!*
!* Updated : July 2001 - L. Giraud, J. Langou
!* Purpose : Avoid to compute the approximate solution at each step of
!*           the Krylov space construction when spA is zero.
!*
!* Updated : November 2002 - S. Gratton
!* Purpose : Use Givens rotations conform to the classical definition.
!*           No impact one the convergence history.
!*
!* Updated : November 2002 - L. Giraud
!* Purpose : Properly handle the situation when the convergence is obtained
!*           exactly at the "IterMax" iteration
!*
!* Updated : December 2002 - L. Giraud, J.Langou
!* Purpose : Add the capability to avoid explicit residual calculation at restart
!*
!* Updated : January  2003 - L. Giraud, S. Gratton
!* Purpose : Use Givens rotations from BLAS.
!*
!* Updated : March    2003 - L. Giraud
!* Purpose : Set back retlbl to zero, if initial guess is solution
!*           or right-hand side is zero
!*
!* Updated : September 2003 - L. Giraud
!* Purpose : Include room in the workspace to store the results of the dot products
!*           Fix the bugs that appeared when M > Nloc 
!*
!*  Arguments
!*  =========
!*
!*  n       (input) INTEGER.
!*           On entry, the dimension of the problem.
!*           Unchanged on exit.
!*
!*  m        (input) INTEGER
!*           Restart parameter, <= N. This parameter controls the amount
!*           of memory required for matrix H (see WORK and H).
!*           Unchanged on exit.
!*
!*  b        (input) real/complex
!*           Right hand side of the linear system.
!*
!*  x        (output) real/complex
!*           Computed solution of the linear system.
!*
!*  H        (workspace)  real/complex
!*           Hessenberg matrix built within dgmres
!*
!*  w        (workspace)  real/complex
!*           Vector used as temporary storage
!*
!*  r0       (workspace)  real/complex
!*           Vector used as temporary storage
!*
!*  V        (workspace)  real/complex
!*           Basis computed by the Arnoldi's procedure.
!*  
!*  dot      (workspace) real/complex
!*           Store the results of the dot product calculation
!*
!*  yCurrent (workspace) real/complex
!*           solution of the current LS
!*
!*  xCurrent (workspace) real/complex
!*           current iterate
!*
!*  rotSin   (workspace) real/complex
!*           Sine of the Givens rotation
!*
!*  rotCos   (workspace) real
!*           Cosine of the Givens rotation
!*
!*  irc      (input/output) INTEGER array. length 3
!*             irc(1) : REVCOM   used for reverse communication
!*                              (type of external operation)
!*             irc(2) : COLX     used for reverse communication
!*             irc(3) : COLY     used for reverse communication
!*             irc(4) : COLZ     used for reverse communication
!*             irc(5) : NBSCAL   used for reverse communication
!*
!*  icntl    (input) INTEGER array. length 7
!*             icntl(1) : stdout for error messages
!*             icntl(2) : stdout for warnings
!*             icntl(3) : stdout for convergence history
!*             icntl(4) : 0 - no preconditioning
!*                        1 - left preconditioning
!*                        2 - right preconditioning
!*                        3 - double side preconditioning
!*                        4 - error, default set in Init
!*             icntl(5) : 0 - modified Gram-Schmidt
!*                        1 - iterative modified Gram-Schmidt
!*                        2 - classical Gram-Schmidt
!*                        3 - iterative classical Gram-Schmidt
!*             icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
!*                        1 - user supplied initial guess
!*             icntl(7) : maximum number of iterations
!*             icntl(8) : 1 - default compute the true residual at each restart
!*                        0 - use recurence formula at restart
!*
!*  cntl     (input) real array, length 5
!*             cntl(1) : tolerance for convergence
!*             cntl(2) : scaling factor for normwise perturbation on A
!*             cntl(3) : scaling factor for normwise perturbation on b
!*             cntl(4) : scaling factor for normwise perturbation on the
!*                       preconditioned matrix
!*             cntl(5) : scaling factor for normwise perturbation on
!*                       preconditioned right hand side
!*
!*  info     (output) INTEGER array, length 2
!*             info(1) :  0 - normal exit
!*                       -1 - n < 1
!*                       -2 - m < 1
!*                       -3 - lwork too small
!*                       -4 - convergence not achieved after icntl(7) iterations
!*                       -5 - precondition type not set by user
!*             info(2) : if info(1)=0 - number of iteration to converge
!*                       if info(1)=-3 - minimum workspace size necessary
!*             info(3) : optimal size for the workspace
!*
!* rinfo     (output) real array, length 2
!*             if info(1)=0 
!*               rinfo(1) : backward error for the preconditioned system
!*               rinfo(2) : backward error for the unpreconditioned system
!*
!* Input variables
!* ---------------
        integer  n, m, icntl(*)
        complex b(*)
        real    cntl(*)
!*
!* Output variables
!* ----------------
       integer  info(*)
       real    rinfo(*)
!*
!* Input/Output variables
!* ----------------------
       integer  irc(*)
       complex x(*), H(m+1,*), w(*), r0(*) 
       complex V(n,*),dot(*),yCurrent(*)  
       complex xCurrent(*), rotSin(*)
       complex rotCos(*)
!*
!* Local variables
!* ---------------
       integer  j, jH, iterOut, nOrtho, iterMax, initGuess, iOrthog
       integer  xptr, bptr, wptr, r0ptr, Vptr, Hptr, yptr, xcuptr
       integer  dotptr
       integer  typePrec, leftPrec, rightPrec, dblePrec, noPrec
       integer  iwarn, ihist
       integer  compRsd
       real    beta, bn, sA, sb, sPA, sPb, bea, be, temp
       real    dloo, dnormw, dnormx, dnormres, trueNormRes
       complex dVi, aux
       complex auxHjj, auxHjp1j
!*
       parameter (noPrec = 0, leftPrec = 1)
       parameter (rightPrec = 2, dblePrec = 3)  
!*
       complex ZERO, ONE
       parameter (ZERO = (0.0e0, 0.0e0), ONE = (1.0e0, 0.0e0))
       real DZRO,DONE
       parameter (DZRO = 0.0e0, DONE = 1.0e0)
!*
!*
!* External functions
!* ------------------
!         real    scnrm2
!        external scnrm2
!*
!* Saved variables
!* ---------------
       save iterOut, jH, beta, bn, dnormres, retlbl, j
       save sA, sb, sPA, sPb, dnormx, trueNormRes, bea, be
       save dloo, nOrtho, compRsd
!*
!* Intrinsic function
!* ------------------
       intrinsic abs, sqrt, real, conjg
!*
!*
!* Reverse communication variables
!* -------------------------------
       integer matvec, precondLeft, precondRight, prosca
       parameter(matvec=1, precondLeft=2, precondRight=3, prosca=4)
       integer retlbl
       DATA retlbl /0/
!*
!*       Executable statements
!*
!* setup some pointers on the workspace
       xptr     = 1
       bptr     = xptr + n
       r0ptr    = bptr + n
       wptr     = r0ptr + n
       Vptr     = wptr + n
       if (icntl(8).eq.1) then
         Hptr     = Vptr + m*n
       else
         Hptr     = Vptr + (m+1)*n
       endif
       dotptr   = Hptr + (m+1)*(m+1)
       if ((icntl(5).eq.2).or.(icntl(5).eq.3)) then
         yptr = dotptr + m
       else
         yptr = dotptr + 1
       endif
       xcuptr   = yptr + m
!*
       iwarn      = icntl(2)
       ihist      = icntl(3)
       typePrec   = icntl(4)
       iOrthog    = icntl(5)
       initGuess  = icntl(6)
       iterMax    = icntl(7)
!*
       if (retlbl.eq.0) then
         compRsd    = icntl(8)
       endif
!*
      if (retlbl.ne.0) then
        if (retlbl.eq.5) then
          goto 5
        else if (retlbl.eq.6) then
          goto 6
        else if (retlbl.eq.8) then
          goto 8
        else if (retlbl.eq.11) then
          goto 11
        else if (retlbl.eq.16) then
          goto 16
        else if (retlbl.eq.18) then
          goto 18
        else if (retlbl.eq.21) then
          goto 21
        else if (retlbl.eq.26) then
          goto 26
        else if (retlbl.eq.31) then
          goto 31
        else if (retlbl.eq.32) then
          goto 32
        else if (retlbl.eq.33) then
          goto 33
        else if (retlbl.eq.34) then
          goto 34 
        else if (retlbl.eq.36) then
          goto 36
        else if (retlbl.eq.37) then
          goto 37
        else if (retlbl.eq.38) then
          goto 38
        else if (retlbl.eq.41) then
          goto 41
        else if (retlbl.eq.43) then
          goto 43
        else if (retlbl.eq.46) then
          goto 46
        else if (retlbl.eq.48) then
          goto 48
        else if (retlbl.eq.51) then
          goto 51
        else if (retlbl.eq.52) then
          goto 52
        else if (retlbl.eq.61) then
          goto 61
        else if (retlbl.eq.66) then
          goto 66
        else if (retlbl.eq.68) then
          goto 68
        endif
      endif
!*
!*
!* intialization of various variables
!*
      iterOut  = 0
      beta     = DZRO
!*
      if (initGuess.eq.0) then
        do j=1,n
          x(j) = ZERO
        enddo
      endif
!*
!*        bn = scnrm2(n,b,1)
!*
      irc(1) = prosca
      irc(2) = bptr
      irc(3) = bptr
      irc(4) = dotptr
      irc(5) = 1
      retlbl = 5
      return
 5    continue
      bn = sqrt(real(dot(1)))
!*
      if (bn.eq.DZRO) then
        do j=1,n
          x(j) = ZERO
        enddo  
        if (iwarn.ne.0) then
          write(iwarn,*)
          write(iwarn,*) ' WARNING GMRES : '
          write(iwarn,*) '       Null right hand side'
          write(iwarn,*) '       solution set to zero'
          write(iwarn,*)
        endif
        jH = 0
        bea = DZRO
        be  = DZRO
        write(ihist,'(I5,11x,E11.4,$)') jH,bea
        write(ihist,'(7x,E11.4)') be
        info(1)  = 0
        info(2)  = 0
        rinfo(1) = DZRO
        rinfo(2) = DZRO
        irc(1)   = 0
        retlbl = 0
        return
      endif
!*
!* Compute the scaling factor for the backward error on the 
!*  unpreconditioned sytem
!*
      sA       = cntl(2)
      sb       = cntl(3)
      if ((sA.eq.DZRO).and.(sb.eq.DZRO)) then
        sb = bn
      endif
!* Compute the scaling factor for the backward error on the
!*  preconditioned sytem
!*
       sPA      = cntl(4)
       sPb      = cntl(5)
       if ((sPA.eq.DZRO).and.(sPb.eq.DZRO)) then
         if ((typePrec.eq.noPrec).or.(typePrec.eq.rightPrec)) then
           sPb = bn
         else
           irc(1) = precondLeft
           irc(2) = bptr
           irc(4) = r0ptr
           retlbl = 6
           return
         endif
       endif
 6     continue
       if ((sPA.eq.DZRO).and.(sPb.eq.DZRO)) then
         if ((typePrec.eq.dblePrec).or.(typePrec.eq.leftPrec)) then
!*
!*           sPb = scnrm2(n,r0,1)
!*
           irc(1) = prosca
           irc(2) = r0ptr
           irc(3) = r0ptr
           irc(4) = dotptr
           irc(5) = 1
           retlbl = 8
           return
         endif
       endif
 8     continue
       if ((sPA.eq.DZRO).and.(sPb.eq.DZRO)) then
         if ((typePrec.eq.dblePrec).or.(typePrec.eq.leftPrec)) then
           sPb = sqrt(real(dot(1)))
!* 
         endif
       endif
!*
!*
!* Compute the first residual
!*           Y = AX : r0 <-- A x
!*
!* The residual is computed only if the initial guess is not zero
!*
       if (initGuess.ne.0) then
         irc(1) = matvec
         irc(2) = xptr
         irc(4) = r0ptr
         retlbl = 11
         return
       endif
 11    continue
       if (initGuess.ne.0) then
         do j=1,n
           r0(j) = b(j)-r0(j)
         enddo
       else
         call ccopy(n,b,1,r0,1)
       endif 
!*
!* Compute the preconditioned residual if necessary
!*      M_1Y = X : w <-- M_1^{-1} r0
!*
       if ((typePrec.eq.noPrec).or.(typePrec.eq.rightPrec)) then
         call ccopy(n,r0,1,w,1)
       else
         irc(1) = precondLeft
         irc(2) = r0ptr
         irc(4) = wptr
         retlbl = 16
         return
       endif
 16    continue
!*
!*
!*       beta = scnrm2(n,w,1)
!*
!*
       irc(1) = prosca
       irc(2) = wptr
       irc(3) = wptr
       irc(4) = dotptr
       irc(5) = 1
       retlbl = 18
       return
 18    continue
       beta = sqrt(real(dot(1)))
!*
       if (beta .eq. DZRO) then
!*  The residual is exactly zero : x is the exact solution
         info(1) = 0
         info(2) = 0
         rinfo(1) = DZRO
         rinfo(2) = DZRO
         irc(1)   = 0
         retlbl = 0
         jH = 0
         bea = DZRO
         be  = DZRO
         write(ihist,'(I5,11x,E11.4,$)') jH,bea
         write(ihist,'(7x,E11.4)') be
         if (iwarn.ne.0) then
          write(iwarn,*)
          write(iwarn,*) ' WARNING GMRES : '
          write(iwarn,*) '       Intial residual is zero'
          write(iwarn,*) '       initial guess is solution'
          write(iwarn,*)
         endif
         return
       endif
!*
       aux = ONE/beta
       do j=1,n
         V(j,1) = ZERO
       enddo
       call caxpy(n,aux,w,1,V(1,1),1)
!*
!*       Most outer loop : cgmres iteration
!*
!*       REPEAT
 7     continue
!*
!*
       H(1,m+1)=beta
       do j=1,m
         H(j+1,m+1) = ZERO
       enddo
!*
!*        Construction of the hessenberg matrix WORK and of the orthogonal
!*        basis V such that AV=VH 
!*
       jH = 1
 10    continue
!* Remark : this  do loop has been written with a while do
!*          because the
!*               " do jH=1,restart "
!*         fails with the reverse communication.
!*      do  jH=1,restart
!*
!*
!* Compute the preconditioned residual if necessary
!*
       if ((typePrec.eq.rightPrec).or.(typePrec.eq.dblePrec)) then  
!*
!*           Y = M_2^{-1}X : w <-- M_2^{-1} V(1,jH)
!*
         irc(1) = precondRight
         irc(2) = vptr + (jH-1)*n
         irc(4) = wptr
         retlbl = 21
         return
       else
         call ccopy(n,V(1,jH),1,w,1)
       endif
 21    continue
!*
!*           Y = AX : r0 <-- A w
!*
       irc(1) = matvec
       irc(2) = wptr
       irc(4) = r0ptr
       retlbl = 26
       return
 26    continue
!*
!*      MY = X : w <-- M_1^{-1} r0
!*
       if ((typePrec.eq.noPrec).or.(typePrec.eq.rightPrec)) then
         call ccopy(n,r0,1,w,1)
       else
         irc(1) = precondLeft
         irc(2) = r0ptr
         irc(4) = wptr
         retlbl = 31
         return
       endif
 31    continue
!*
!* Orthogonalization using either MGS or IMGS
!*  
!* initialize the Hessenberg matrix to zero in order to be able to use
!*     IMGS as orthogonalization procedure.
       do j=1,jH
         H(j,jH) = ZERO
       enddo
       nOrtho = 0
 19    continue
       nOrtho = nOrtho +1
       dloo   = DZRO
!*
       if ((iOrthog.eq.0).or.(iOrthog.eq.1)) then
!* MGS
!*
!*           do j=1,jH
!*
         j = 1
!*           REPEAT
       endif
 23    continue
       if ((iOrthog.eq.0).or.(iOrthog.eq.1)) then
!*
!*             dVi     = cdotc(n,V(1,j),1,w,1)
!*
         irc(1) = prosca
         irc(2) = vptr + (j-1)*n
         irc(3) = wptr
         irc(4) = dotptr
         irc(5) = 1
         retlbl = 32
         return
       endif
 32    continue
       if ((iOrthog.eq.0).or.(iOrthog.eq.1)) then
         dVi     = dot(1)
         H(j,jH) = H(j,jH) + dVi
         dloo    = dloo + abs(dVi)**2
         aux = -ONE*dVi
         call caxpy(n,aux,V(1,j),1,w,1)
         j = j + 1
         if (j.le.jH) goto 23
!*          enddo_j
       else
!* CGS
!* Gathered dot product calculation
!*
!*           call cgemv('C',n,jH,ONE,V(1,1),n,w,1,ZERO,r0,1)
!*
         irc(1) = prosca
         irc(2) = vptr
         irc(3) = wptr
         irc(4) = dotptr
         irc(5) = jH
         retlbl = 34
         return
       endif
 34    continue
       if ((iOrthog.eq.2).or.(iOrthog.eq.3)) then
!*
         call caxpy(jH,ONE,dot,1,H(1,jH),1)
         call cgemv('N',n,jH,-ONE,V(1,1),n,dot,1,ONE,w,1)
         dloo = scnrm2(jH,dot,1)**2
       endif
!*
!*         dnormw = scnrm2(n,w,1)
!*
       irc(1) = prosca
       irc(2) = wptr
       irc(3) = wptr
       irc(4) = dotptr
       irc(5) = 1
       retlbl = 33
       return
 33    continue
       dnormw = sqrt(real(dot(1)))
!*
       if ((iOrthog.eq.1).or.(iOrthog.eq.3)) then
!* IMGS / CGS orthogonalisation
         dloo = sqrt(dloo)
!* check the orthogonalization quality
         if (((2.0*dnormw).le.dloo).and.(nOrtho.lt.3)) then
           goto 19
         endif
       endif
!*
       H(jH+1,jH) = dnormw
       if ((jH.lt.m).or.(icntl(8).eq.0)) then
         aux = ONE/dnormw
         do j=1,n
           V(j,jH+1) = ZERO
         enddo
         call caxpy(n,aux,w,1,V(1,jH+1),1)
       endif
!* Apply previous Givens rotations to the new column of H
       do j=1,jH-1
         call crot(1, H(j,jH), 1, H(j+1,jH), 1,real(rotCos(j)),  &
     &             rotSin(j))
       enddo
       auxHjj = H(jH,jH)
       auxHjp1j= H(jH+1,jH)
       call crotg(auxHjj, auxHjp1j,temp,rotSin(jH))
       rotCos(jH)= dcmplx(temp, DZRO)
!* Apply current rotation to the rhs of the least squares problem
       call crot(1, H(jH,m+1), 1, H(jH+1,m+1), 1, real(rotCos(jH)),&
     &            rotSin(jH))
!*
!* zabs(H(jH+1,m+1)) is the residual computed using the least squares
!*          solver
!* Complete the QR factorisation of the Hessenberg matrix by apply the current
!* rotation to the last entry of the collumn
       call crot(1, H(jH,jH), 1, H(jH+1,jH), 1, real(rotCos(jH)),&
     &           rotSin(jH))
       H(jH+1,jH) = ZERO
!*
!* Get the Least square residual
!*
       dnormres = abs(H(jH+1,m+1))
       if (sPa.ne.DZRO) then
!*
!* Compute the solution of the current linear least squares problem
!*
         call ccopy(jH,H(1,m+1),1,yCurrent,1)
         call ctrsv('U','N','N',jH,H,m+1,yCurrent,1)
!*
!* Compute the value of the new iterate 
!*
         call cgemv('N',n,jH,ONE,v,n,                &
     &            yCurrent,1,ZERO,xCurrent,1)
!*
         if ((typePrec.eq.rightPrec).or.(typePrec.eq.dblePrec)) then  
!*
!*         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent
!*
           irc(1) = precondRight
           irc(2) = xcuptr
           irc(4) = r0ptr
           retlbl = 36
           return
         else
           call ccopy(n,xCurrent,1,r0,1)
         endif
       endif
 36    continue
!*
!*
       if (sPa.ne.DZRO) then
!* Update the current solution
         call ccopy(n,x,1,xCurrent,1)
         call caxpy(n,ONE,r0,1,xCurrent,1)
!*
!*         dnormx = scnrm2(n,xCurrent,1)
!*
         irc(1) = prosca
         irc(2) = xcuptr
         irc(3) = xcuptr
         irc(4) = dotptr
         irc(5) = 1
         retlbl = 38
         return
       else
         dnormx    = DONE
       endif
 38    continue
       if (sPa.ne.DZRO) then
         dnormx = sqrt(real(dot(1)))
       endif
!*
       bea = dnormres/(sPA*dnormx+sPb)
!*
!* Check the convergence based on the Arnoldi Backward error for the
!* preconditioned system
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then  
!* 
!* The Arnoldi Backward error indicates that cgmres might have converge
!* enforce the calculation of the true residual at next restart
         compRsd = 1
!*
!*  If the update of X has not yet been performed
         if (sPA.eq.DZRO) then
!*
!* Compute the solution of the current linear least squares problem
!*
           call ccopy(jH,H(1,m+1),1,yCurrent,1)
           call ctrsv('U','N','N',jH,H,m+1,yCurrent,1)
!*
!* Compute the value of the new iterate 
!*
           call cgemv('N',n,jH,ONE,v,n,     &
     &            yCurrent,1,ZERO,xCurrent,1)
!*
           if ((typePrec.eq.rightPrec).or.(typePrec.eq.dblePrec)) then
!* 
!*         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent
!*
             irc(1) = precondRight
             irc(2) = xcuptr
             irc(4) = r0ptr 
             retlbl = 37
             return
           else
             call ccopy(n,xCurrent,1,r0,1)
           endif
         endif
       endif
 37    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         if (sPA.eq.DZRO) then
!* Update the current solution
            call ccopy(n,x,1,xCurrent,1)
            call caxpy(n,ONE,r0,1,xCurrent,1)
         endif
!*
         call ccopy(n,xCurrent,1,r0,1)
!* Compute the true residual, the Arnoldi one may be unaccurate
!*
!*           Y = AX : w  <-- A r0
!*
         irc(1) = matvec
         irc(2) = r0ptr
         irc(4) = wptr
         retlbl = 41
         return
       endif
 41    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
!*
         do j=1,n
           w(j) = b(j) - w(j)
         enddo
!* Compute the norm of the unpreconditioned residual
!*
!*        trueNormRes = scnrm2(n,w,1)
!*
         irc(1) = prosca
         irc(2) = wptr
         irc(3) = wptr
         irc(4) = dotptr
         irc(5) = 1
         retlbl = 43
         return
       endif
 43    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         trueNormRes = sqrt(real(dot(1)))
!*
         if ((typePrec.eq.leftPrec).or.(typePrec.eq.dblePrec)) then  
!*
!*      MY = X : r0 <-- M_1^{-1} w 
!*
           irc(1) = precondLeft
           irc(2) = wptr
           irc(4) = r0ptr
           retlbl = 46
           return
         else 
           call ccopy(n,w,1,r0,1)
         endif
       endif
 46    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
!*
!*        dnormres = scnrm2(n,r0,1) 
!*
         irc(1) = prosca
         irc(2) = r0ptr
         irc(3) = r0ptr
         irc(4) = dotptr
         irc(5) = 1
         retlbl = 48
         return
       endif
 48    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         dnormres = sqrt(real(dot(1)))
!*
         be = dnormres/(sPA*dnormx+sPb)
!* Save the backward error on a file if convergence history requested
         if (ihist.ne.0) then
           write(ihist,1000)iterOut*m+jH,bea,be
1000  format(I5,11x,E11.4,7x,E11.4)
         endif
!*
       endif
!*
!*
!* Check again the convergence
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then   
         if ((be.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then   
!* The convergence has been achieved, we restore the solution in x
!* and compute the two backward errors.
           call ccopy(n,xCurrent,1,x,1)
!*
           if (sA.ne.DZRO) then
!*
!*            dnormx = scnrm2(n,x,1)
!*
             irc(1) = prosca
             irc(2) = xptr
             irc(3) = xptr
             irc(4) = dotptr
             irc(5) = 1
             retlbl = 51
             return
           endif
         endif
       endif
 51    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         if ((be.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
           if (sA.ne.DZRO) then
             dnormx = sqrt(real(dot(1)))
!*
           else
             dnormx = DONE
           endif
!* Return the backward errors
           rinfo(1) = be
           rinfo(2) = trueNormRes/(sA*dnormx+sb)
           if (be.le.cntl(1)) then
             info(1) = 0
             if (ihist.ne.0) then
               write(ihist,*)
               write(ihist,'(A20)') 'Convergence achieved'
             endif
           else if (be.gt.cntl(1)) then
             if (iwarn.ne.0) then
               write(iwarn,*)
               write(iwarn,*) ' WARNING GMRES : '
               write(iwarn,*) '       No convergence after '
               write(iwarn,*) iterOut*m+jH,' iterations '
               write(iwarn,*)
             endif
             if (ihist.ne.0) then
               write(ihist,*)
               write(ihist,*) ' WARNING GMRES :'
               write(ihist,*) '       No convergence after '
               write(ihist,*) iterOut*m+jH,' iterations '
               write(ihist,*)
             endif
             info(1) = -4
           endif
           if (ihist.ne.0) then
             write(ihist,1010) rinfo(1)
             write(ihist,1011) rinfo(2)
1010  format('B.E. on the preconditioned system:   ',E11.4)
1011  format('B.E. on the unpreconditioned system: ',E11.4)
           endif
           info(2) = iterOut*m+jH
           if (ihist.ne.0) then
             write(ihist,'(A10,I2)') 'info(1) = ',info(1)
             write(ihist,'(A32,I5)')  &
     &                'Number of iterations (info(2)): ',info(2)  
           endif
           irc(1)  = 0
           retlbl  = 0
           return
         endif
       else
!* Save the backward error on a file if convergence history requested
         if (ihist.ne.0) then
           write(ihist,1001)iterOut*m+jH,bea
1001  format(I5,11x,E11.4,7x,'--')
         endif
!*
       endif  
!*
       jH = jH + 1
       if (jH.le.m) then
         goto 10
       endif
!*
       iterOut = iterOut + 1
!*
!* we have completed the Krylov space construction, we restart if
!* we have not yet exceeded the maximum number of iterations allowed.
!*
       if ((sPa.eq.DZRO).and.(bea.gt.cntl(1))) then
!*
!* Compute the solution of the current linear least squares problem
!*
         jH = jH - 1
         call ccopy(jH,H(1,m+1),1,yCurrent,1)
         call ctrsv('U','N','N',jH,H,m+1,yCurrent,1)
!*
!* Compute the value of the new iterate
!*
         call cgemv('N',n,jH,ONE,v,n,  &
     &            yCurrent,1,ZERO,xCurrent,1)
!*
         if ((typePrec.eq.rightPrec).or.(typePrec.eq.dblePrec)) then
!*
!*         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent
!*
           irc(1) = precondRight
           irc(2) = xcuptr
           irc(4) = r0ptr
           retlbl = 52
           return
         else
           call ccopy(n,xCurrent,1,r0,1)
         endif
       endif
 52    continue
       if ((sPa.eq.DZRO).and.(bea.gt.cntl(1))) then
!* Update the current solution
         call ccopy(n,x,1,xCurrent,1)
         call caxpy(n,ONE,r0,1,xCurrent,1)
       endif
!*
       call ccopy(n,xCurrent,1,x,1)
!*
       if (compRsd.eq.1) then
!*
!* Compute the true residual
!*
         call ccopy(n,x,1,w,1)
         irc(1) = matvec
         irc(2) = wptr
         irc(4) = r0ptr
         retlbl = 61
         return
       endif
 61    continue
       if (compRsd.eq.1) then
         do j=1,n
           r0(j) = b(j) - r0(j)
         enddo
!*
!* Precondition the new residual if necessary
!*
         if ((typePrec.eq.leftPrec).or.(typePrec.eq.dblePrec)) then
!*
!*      MY = X : w <-- M_1^{-1} r0
!*
           irc(1) = precondLeft
           irc(2) = r0ptr
           irc(4) = wptr
           retlbl = 66
           return
         else
           call ccopy(n,r0,1,w,1)
         endif
       endif
 66    continue
!*
!*           beta = scnrm2(n,w,1)
!*
       if (compRsd.eq.1) then
         irc(1) = prosca
         irc(2) = wptr
         irc(3) = wptr
         irc(4) = dotptr
         irc(5) = 1
         retlbl = 68
         return
       endif
 68    continue
       if (compRsd.eq.1) then
         beta = sqrt(real(dot(1)))
!*
       else
!* Use recurrence to approximate the residual at restart
         beta = abs(H(m+1,m+1))
!* Apply the Givens rotation is the reverse order
         do j=m,1,-1
           H(j,m+1)   = ZERO
           call crot(1, H(j,m+1), 1, H(j+1,m+1), 1,  &
     &               real(rotCos(j)), -rotSin(j))
         enddo
!*
!* On applique les vecteurs V
!*
         call cgemv('N',n,m+1,ONE,v,n,H(1,m+1),1,ZERO,w,1)
!*
       endif
       do j=1,n
         V(j,1) = ZERO
       enddo
       aux = ONE/beta
       call caxpy(n,aux,w,1,V(1,1),1)
!*
       goto 7
!*
       END SUBROUTINE
!---------------------------------------------------------------

      SUBROUTINE CGEMV (TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y,&
     &   INCY)
! C***BEGIN PROLOGUE  CGEMV
! C***PURPOSE  Multiply a complex vector by a complex general matrix.
! C***LIBRARY   SLATEC (BLAS)
! C***CATEGORY  D1B4
! C***TYPE      COMPLEX (SGEMV-S, DGEMV-D, CGEMV-C)
! C***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
! C***AUTHOR  Dongarra, J. J., (ANL)
! C           Du Croz, J., (NAG)
! C           Hammarling, S., (NAG)
! C           Hanson, R. J., (SNLA)
! C***DESCRIPTION
! C
! C  CGEMV  performs one of the matrix-vector operations
! C
! C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
! C
! C     y := alpha*conjg( A' )*x + beta*y,
! C
! C  where alpha and beta are scalars, x and y are vectors and A is an
! C  m by n matrix.
! C
! C  Parameters
! C  ==========
! C
! C  TRANS  - CHARACTER*1.
! C           On entry, TRANS specifies the operation to be performed as
! C           follows:
! C
! C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
! C
! C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
! C
! C              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
! C
! C           Unchanged on exit.
! C
! C  M      - INTEGER.
! C           On entry, M specifies the number of rows of the matrix A.
! C           M must be at least zero.
! C           Unchanged on exit.
! C
! C  N      - INTEGER.
! C           On entry, N specifies the number of columns of the matrix A.
! C           N must be at least zero.
! C           Unchanged on exit.
! C
! C  ALPHA  - COMPLEX         .
! C           On entry, ALPHA specifies the scalar alpha.
! C           Unchanged on exit.
! C
! C  A      - COMPLEX          array of DIMENSION ( LDA, n ).
! C           Before entry, the leading m by n part of the array A must
! C           contain the matrix of coefficients.
! C           Unchanged on exit.
! C
! C  LDA    - INTEGER.
! C           On entry, LDA specifies the first dimension of A as declared
! C           in the calling (sub) program. LDA must be at least
! C           max( 1, m ).
! C           Unchanged on exit.
! C
! C  X      - COMPLEX          array of DIMENSION at least
! C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
! C           and at least
! C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
! C           Before entry, the incremented array X must contain the
! C           vector x.
! C           Unchanged on exit.
! C
! C  INCX   - INTEGER.
! C           On entry, INCX specifies the increment for the elements of
! C           X. INCX must not be zero.
! C           Unchanged on exit.
! C
! C  BETA   - COMPLEX         .
! C           On entry, BETA specifies the scalar beta. When BETA is
! C           supplied as zero then Y need not be set on input.
! C           Unchanged on exit.
! C
! C  Y      - COMPLEX          array of DIMENSION at least
! C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
! C           and at least
! C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
! C           Before entry with BETA non-zero, the incremented array Y
! C           must contain the vector y. On exit, Y is overwritten by the
! C           updated vector y.
! C
! C  INCY   - INTEGER.
! C           On entry, INCY specifies the increment for the elements of
! C           Y. INCY must not be zero.
! C           Unchanged on exit.
! C
! C***REFERENCES  Dongarra, J. J., Du Croz, J., Hammarling, S., and
! C                 Hanson, R. J.  An extended set of Fortran basic linear
! C                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1,
! C                 pp. 1-17, March 1988.
! C***ROUTINES CALLED  LSAME, XERBLA
! C***REVISION HISTORY  (YYMMDD)
! C   861022  DATE WRITTEN
! C   910605  Modified to meet SLATEC prologue standards.  Only comment
! C           lines were modified.  (BKS)
! C***END PROLOGUE  CGEMV
! C     .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
! C     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
! C     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
      LOGICAL            NOCONJ
! C     .. External Functions ..
!       LOGICAL            LSAME
!       EXTERNAL           LSAME
! C     .. External Subroutines ..
!       EXTERNAL           XERBLA
! C     .. Intrinsic Functions ..
       INTRINSIC          CONJG, MAX
! C***FIRST EXECUTABLE STATEMENT  CGEMV
! C
! C     Test the input parameters.
! C

      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.        &
     &         .NOT.LSAME( TRANS, 'T' ).AND.        &
     &         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGEMV ', INFO )
         RETURN
      END IF
! C
! C     Quick return if possible.
! C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.               &
     &    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )&
     &   RETURN
! C
      NOCONJ = LSAME( TRANS, 'T' )
! C
! C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
! C     up the start points in  X  and  Y.
! C
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
! C
! C     Start the operations. In this version the elements of A are
! C     accessed sequentially with one pass through A.
! C
! C     First form  y := beta*y.
! C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO ) RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
! C
! C        Form  y := alpha*A*x + y.
! C
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
! C
! C        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
! C
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 110, J = 1, N
               TEMP = ZERO
               IF( NOCONJ )THEN
                  DO 90, I = 1, M
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
               ELSE
                  DO 100, I = 1, M
                     TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  100             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  110       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               IF( NOCONJ )THEN
                  DO 120, I = 1, M
                     TEMP = TEMP + A( I, J )*X( IX )
                     IX   = IX   + INCX
  120             CONTINUE
               ELSE
                  DO 130, I = 1, M
                     TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
                     IX   = IX   + INCX
  130             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  140       CONTINUE
         END IF
      END IF
! C
      RETURN
! C
! C     End of CGEMV .
! C
      END SUBROUTINE
!-----------------------------------------------------------------------
      COMPLEX FUNCTION CDOTC(N,CX,INCX,CY,INCY)
! *     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
! *     ..
! *     .. Array Arguments ..
      COMPLEX CX(*),CY(*)
! *     ..
! *
! *  Purpose
! *  =======
! *
! *     forms the dot product of two vectors, conjugating the first
! *     vector.
! *
! *  Further Details
! *  ===============
! *
! *     jack dongarra, linpack,  3/11/78.
! *     modified 12/3/93, array(1) declarations changed to array(*)
! *
! *  =====================================================================
! *
! *     .. Local Scalars ..
      COMPLEX CTEMP
      INTEGER I,IX,IY
! *     ..
! *     .. Intrinsic Functions ..
!       INTRINSIC CONJG
! *   

      CTEMP = (0.0,0.0)
      CDOTC = (0.0,0.0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
! *
! *        code for both increments equal to 1
! *
         DO I = 1,N
            CTEMP = CTEMP + CONJG(CX(I))*CY(I)
         END DO
      ELSE

! *
! *        code for unequal increments or equal increments
! *          not equal to 1
! *
         IX = 1
         IY = 1

         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1

         DO I = 1,N
            CTEMP = CTEMP + CONJG(CX(IX))*CY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      CDOTC = CTEMP
      RETURN

      END FUNCTION
!-----------------------------------------------------------------
      SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY)
! *     .. Scalar Arguments ..
      COMPLEX CA
      INTEGER INCX,INCY,N
! *     ..
! *     .. Array Arguments ..
      COMPLEX CX(*),CY(*)
! *     ..
! *
! *  Purpose
! *  =======
! *
! *     CAXPY constant times a vector plus a vector.
! *
! *  Further Details
! *  ===============
! *
! *     jack dongarra, linpack, 3/11/78.
! *     modified 12/3/93, array(1) declarations changed to array(*)
! *
! *  =====================================================================
! *
! *     .. Local Scalars ..
      INTEGER I,IX,IY
! *     ..
! *     .. External Functions ..
!       REAL SCABS1
!       EXTERNAL SCABS1
! *     ..
      IF (N.LE.0) RETURN
      IF (SCABS1(CA).EQ.0.0E+0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
! *
! *        code for both increments equal to 1
! *
         DO I = 1,N
            CY(I) = CY(I) + CA*CX(I)
         END DO
      ELSE
! *
! *        code for unequal increments or equal increments
! *          not equal to 1
! *
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CY(IY) = CY(IY) + CA*CX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN

      END SUBROUTINE
!-----------------------------------------------------------------------

      SUBROUTINE CCOPY(N,CX,INCX,CY,INCY)
! *     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
! *     ..
! *     .. Array Arguments ..
      COMPLEX CX(*),CY(*)
! *     ..
! *
! *  Purpose
! *  =======
! *
! *     CCOPY copies a vector x to a vector y.
! *
! *  Further Details
! *  ===============
! *
! *     jack dongarra, linpack, 3/11/78.
! *     modified 12/3/93, array(1) declarations changed to array(*)
! *
! *  =====================================================================
! *
! *     .. Local Scalars ..
      INTEGER I,IX,IY
! *     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
! *
! *        code for both increments equal to 1
! *
         DO I = 1,N
            CY(I) = CX(I)
         END DO
      ELSE
! *
! *        code for unequal increments or equal increments
! *          not equal to 1
! *
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CY(IY) = CX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN

      END SUBROUTINE
!---------------------------------------------------------------------
      LOGICAL FUNCTION LSAME(CA,CB)
! *
! *  -- LAPACK auxiliary routine (version 3.1) --
! *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
! *     November 2006
! *
! *     .. Scalar Arguments ..
      CHARACTER CA,CB
! *     ..
! *
! *  Purpose
! *  =======
! *
! *  LSAME returns .TRUE. if CA is the same letter as CB regardless of
! *  case.
! *
! *  Arguments
! *  =========
! *
! *  CA      (input) CHARACTER*1
! *
! *  CB      (input) CHARACTER*1
! *          CA and CB specify the single characters to be compared.
! *
! * =====================================================================
! *
! *     .. Intrinsic Functions ..
      INTRINSIC ICHAR
! *     ..
! *     .. Local Scalars ..
      INTEGER INTA,INTB,ZCODE
! *     ..
! *
! *     Test if the characters are equal
! *
      LSAME = CA .EQ. CB
      IF (LSAME) RETURN
! *
! *     Now test for equivalence if both characters are alphabetic.
! *
      ZCODE = ICHAR('Z')
! *
! *     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
! *     machines, on which ICHAR returns a value with bit 8 set.
! *     ICHAR('A') on Prime machines returns 193 which is the same as
! *     ICHAR('A') on an EBCDIC machine.
! *
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
! *
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
! *
! *        ASCII is assumed - ZCODE is the ASCII code of either lower or
! *        upper case 'Z'.
! *
          IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
          IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
! *
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
! *
! *        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
! *        upper case 'Z'.
! *
          IF (INTA.GE.129 .AND. INTA.LE.137 .OR.  &
     &        INTA.GE.145 .AND. INTA.LE.153 .OR.  &
     &        INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
          IF (INTB.GE.129 .AND. INTB.LE.137 .OR.  &
     &        INTB.GE.145 .AND. INTB.LE.153 .OR.  &
     &        INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
! *
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
! *
! *        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
! *        plus 128 of either lower or upper case 'Z'.
! *
          IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
          IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
      LSAME = INTA .EQ. INTB
! *
! *     RETURN
! *
! *     End of LSAME
! *
      END FUNCTION
!-----------------------------------------------------------------
      SUBROUTINE XERBLA( SRNAME, INFO )
! *
! *  -- LAPACK auxiliary routine (preliminary version) --
! *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
! *     November 2006
! *
! *     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
! *     ..
! *
! *  Purpose
! *  =======
! *
! *  XERBLA  is an error handler for the LAPACK routines.
! *  It is called by an LAPACK routine if an input parameter has an
! *  invalid value.  A message is printed and execution stops.
! *
! *  Installers may consider modifying the STOP statement in order to
! *  call system-specific exception-handling facilities.
! *
! *  Arguments
! *  =========
! *
! *  SRNAME  (input) CHARACTER*(*)
! *          The name of the routine which called XERBLA.
! *
! *  INFO    (input) INTEGER
! *          The position of the invalid parameter in the parameter list
! *          of the calling routine.
! *
! * =====================================================================
! *
! *     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
! *     ..
! *     .. Executable Statements ..
! *
      WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
! *
      STOP
! *
 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ', &
     &      'an illegal value' )
! *
! *     End of XERBLA
! *
      END SUBROUTINE
!-------------------------------------------------------------

      REAL FUNCTION SCABS1(Z)
! *     .. Scalar Arguments ..
      COMPLEX Z
! *     ..
! *
! *  Purpose
! *  =======
! *
! *  SCABS1 computes absolute value of a complex number
! *
! *  =====================================================================
! *
! *     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,REAL
! *     ..
      SCABS1 = ABS(REAL(Z)) + ABS(AIMAG(Z))
      RETURN

      END FUNCTION
!---------------------------------------------------------------

      SUBROUTINE CROTG(CA,CB,C,S)
! *     .. Scalar Arguments ..
      COMPLEX CA,CB,S
      REAL C
! *     ..
! *
! *  Purpose
! *  =======
! *
! *  CROTG determines a complex Givens rotation.
! *
! *  =====================================================================
! *
! *     .. Local Scalars ..
      COMPLEX ALPHA
      REAL NORM,SCALE
! *     ..
! *     .. Intrinsic Functions ..
      INTRINSIC CABS,CONJG,SQRT
! *     ..
      IF (CABS(CA).EQ.0.) THEN
         C = 0.
         S = (1.,0.)
         CA = CB
      ELSE
         SCALE = CABS(CA) + CABS(CB)
         NORM = SCALE*SQRT((CABS(CA/SCALE))**2+ (CABS(CB/SCALE))**2)
         ALPHA = CA/CABS(CA)
         C = CABS(CA)/NORM
         S = ALPHA*CONJG(CB)/NORM
         CA = ALPHA*NORM
      END IF
      RETURN

      END SUBROUTINE
!----------------------------------------------------

      SUBROUTINE CTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
! *     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
! *     ..
! *     .. Array Arguments ..
      COMPLEX A(LDA,*),X(*)
! *     ..
! *
! *  Purpose
! *  =======
! *
! *  CTRSV  solves one of the systems of equations
! *
! *     A*x = b,   or   A**T*x = b,   or   A**H*x = b,
! *
! *  where b and x are n element vectors and A is an n by n unit, or
! *  non-unit, upper or lower triangular matrix.
! *
! *  No test for singularity or near-singularity is included in this
! *  routine. Such tests must be performed before calling this routine.
! *
! *  Arguments
! *  ==========
! *
! *  UPLO   - CHARACTER*1.
! *           On entry, UPLO specifies whether the matrix is an upper or
! *           lower triangular matrix as follows:
! *
! *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
! *
! *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
! *
! *           Unchanged on exit.
! *
! *  TRANS  - CHARACTER*1.
! *           On entry, TRANS specifies the equations to be solved as
! *           follows:
! *
! *              TRANS = 'N' or 'n'   A*x = b.
! *
! *              TRANS = 'T' or 't'   A**T*x = b.
! *
! *              TRANS = 'C' or 'c'   A**H*x = b.
! *
! *           Unchanged on exit.
! *
! *  DIAG   - CHARACTER*1.
! *           On entry, DIAG specifies whether or not A is unit
! *           triangular as follows:
! *
! *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
! *
! *              DIAG = 'N' or 'n'   A is not assumed to be unit
! *                                  triangular.
! *
! *           Unchanged on exit.
! *
! *  N      - INTEGER.
! *           On entry, N specifies the order of the matrix A.
! *           N must be at least zero.
! *           Unchanged on exit.
! *
! *  A      - COMPLEX          array of DIMENSION ( LDA, n ).
! *           Before entry with  UPLO = 'U' or 'u', the leading n by n
! *           upper triangular part of the array A must contain the upper
! *           triangular matrix and the strictly lower triangular part of
! *           A is not referenced.
! *           Before entry with UPLO = 'L' or 'l', the leading n by n
! *           lower triangular part of the array A must contain the lower
! *           triangular matrix and the strictly upper triangular part of
! *           A is not referenced.
! *           Note that when  DIAG = 'U' or 'u', the diagonal elements of
! *           A are not referenced either, but are assumed to be unity.
! *           Unchanged on exit.
! *
! *  LDA    - INTEGER.
! *           On entry, LDA specifies the first dimension of A as declared
! *           in the calling (sub) program. LDA must be at least
! *           max( 1, n ).
! *           Unchanged on exit.
! *
! *  X      - COMPLEX          array of dimension at least
! *           ( 1 + ( n - 1 )*abs( INCX ) ).
! *           Before entry, the incremented array X must contain the n
! *           element right-hand side vector b. On exit, X is overwritten
! *           with the solution vector x.
! *
! *  INCX   - INTEGER.
! *           On entry, INCX specifies the increment for the elements of
! *           X. INCX must not be zero.
! *           Unchanged on exit.
! *
! *  Further Details
! *  ===============
! *
! *  Level 2 Blas routine.
! *
! *  -- Written on 22-October-1986.
! *     Jack Dongarra, Argonne National Lab.
! *     Jeremy Du Croz, Nag Central Office.
! *     Sven Hammarling, Nag Central Office.
! *     Richard Hanson, Sandia National Labs.
! *
! *  =====================================================================
! *
! *     .. Parameters ..
      COMPLEX ZERO
      PARAMETER (ZERO= (0.0E+0,0.0E+0))
! *     ..
! *     .. Local Scalars ..
      COMPLEX TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
! *     ..
! *     .. External Functions ..
!       LOGICAL LSAME
!       EXTERNAL LSAME
! *     ..
! *     .. External Subroutines ..
!        EXTERNAL XERBLA
! *     ..
! ! *     .. Intrinsic Functions ..
       INTRINSIC CONJG,MAX
! *     ..
! *
! *     Test the input parameters.
! *
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
     &         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('CTRSV ',INFO)
          RETURN
      END IF
! *
! *     Quick return if possible.
! *
      IF (N.EQ.0) RETURN
! *
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
! *
! *     Set up the start point in X if the increment is not unity. This
! *     will be  ( N - 1 )*INCX  too small for descending loops.
! *
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
! *
! *     Start the operations. In this version the elements of A are
! *     accessed sequentially with one pass through A.
! *
      IF (LSAME(TRANS,'N')) THEN
! *
! *        Form  x := inv( A )*x.
! *
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/A(J,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,1,-1
                              X(I) = X(I) - TEMP*A(I,J)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 40 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                          TEMP = X(JX)
                          IX = JX
                          DO 30 I = J - 1,1,-1
                              IX = IX - INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/A(J,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,N
                              X(I) = X(I) - TEMP*A(I,J)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                          TEMP = X(JX)
                          IX = JX
                          DO 70 I = J + 1,N
                              IX = IX + INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
! *
! *        Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
! *
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 110 J = 1,N
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          DO 90 I = 1,J - 1
                              TEMP = TEMP - A(I,J)*X(I)
   90                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      ELSE
                          DO 100 I = 1,J - 1
                              TEMP = TEMP - CONJG(A(I,J))*X(I)
  100                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(J,J))
                      END IF
                      X(J) = TEMP
  110             CONTINUE
              ELSE
                  JX = KX
                  DO 140 J = 1,N
                      IX = KX
                      TEMP = X(JX)
                      IF (NOCONJ) THEN
                          DO 120 I = 1,J - 1
                              TEMP = TEMP - A(I,J)*X(IX)
                              IX = IX + INCX
  120                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      ELSE
                          DO 130 I = 1,J - 1
                              TEMP = TEMP - CONJG(A(I,J))*X(IX)
                              IX = IX + INCX
  130                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(J,J))
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
  140             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 170 J = N,1,-1
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          DO 150 I = N,J + 1,-1
                              TEMP = TEMP - A(I,J)*X(I)
  150                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      ELSE
                          DO 160 I = N,J + 1,-1
                              TEMP = TEMP - CONJG(A(I,J))*X(I)
  160                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(J,J))
                      END IF
                      X(J) = TEMP
  170             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 200 J = N,1,-1
                      IX = KX
                      TEMP = X(JX)
                      IF (NOCONJ) THEN
                          DO 180 I = N,J + 1,-1
                              TEMP = TEMP - A(I,J)*X(IX)
                              IX = IX - INCX
  180                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      ELSE
                          DO 190 I = N,J + 1,-1
                              TEMP = TEMP - CONJG(A(I,J))*X(IX)
                              IX = IX - INCX
  190                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(J,J))
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
  200             CONTINUE
              END IF
          END IF
      END IF
! *
      RETURN
! *
! *     End of CTRSV .
! *
      END SUBROUTINE
!----------------------------------------------------------

      REAL FUNCTION SCNRM2(N,X,INCX)
! *     .. Scalar Arguments ..
      INTEGER INCX,N
! *     ..
! *     .. Array Arguments ..
      COMPLEX X(*)
! *     ..
! *
! *  Purpose
! *  =======
! *
! *  SCNRM2 returns the euclidean norm of a vector via the function
! *  name, so that
! *
! *     SCNRM2 := sqrt( x**H*x )
! *
! *  Further Details
! *  ===============
! *
! *  -- This version written on 25-October-1982.
! *     Modified on 14-October-1993 to inline the call to CLASSQ.
! *     Sven Hammarling, Nag Ltd.
! *
! *  =====================================================================
! *
! *     .. Parameters ..
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
! *     ..
! *     .. Local Scalars ..
      REAL NORM,SCALE,SSQ,TEMP
      INTEGER IX
! *     ..
! *     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,REAL,SQRT
! *     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE
          SCALE = ZERO
          SSQ = ONE
! *        The following loop is equivalent to this call to the LAPACK
! *        auxiliary routine:
! *        CALL CLASSQ( N, X, INCX, SCALE, SSQ )
! *
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (REAL(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(REAL(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
              IF (AIMAG(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(AIMAG(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
! *
      SCNRM2 = NORM
      RETURN
! *
! *     End of SCNRM2.
! *
      END FUNCTION
!-----------------------------------------------------------

      SUBROUTINE CROT( N, CX, INCX, CY, INCY, C, S )
! 	*
! 	*  -- LAPACK auxiliary routine (version 3.0) --
! 	*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
! 	*     Courant Institute, Argonne National Lab, and Rice University
! 	*     October 31, 1992
! 	*
! 	*     .. Scalar Arguments ..
	      INTEGER            INCX, INCY, N
	      REAL               C
	      COMPLEX            S
! 	*     ..
! 	*     .. Array Arguments ..
	      COMPLEX            CX( * ), CY( * )
! 	*     ..
! 	*
! 	*  Purpose
! 	*  =======
! 	*
! 	*  CROT   applies a plane rotation, where the cos (C) is real and the
! 	*  sin (S) is complex, and the vectors CX and CY are complex.
! 	*
! 	*  Arguments
! 	*  =========
! 	*
! 	*  N       (input) INTEGER
! 	*          The number of elements in the vectors CX and CY.
! 	*
! 	*  CX      (input/output) COMPLEX array, dimension (N)
! 	*          On input, the vector X.
! 	*          On output, CX is overwritten with C*X + S*Y.
! 	*
! 	*  INCX    (input) INTEGER
! 	*          The increment between successive values of CY.  INCX <> 0.
! 	*
! 	*  CY      (input/output) COMPLEX array, dimension (N)
! 	*          On input, the vector Y.
! 	*          On output, CY is overwritten with -CONJG(S)*X + C*Y.
! 	*
! 	*  INCY    (input) INTEGER
! 	*          The increment between successive values of CY.  INCX <> 0.
! 	*
! 	*  C       (input) REAL
! 	*  S       (input) COMPLEX
! 	*          C and S define a rotation
! 	*             [  C          S  ]
! 	*             [ -conjg(S)   C  ]
! 	*          where C*C + S*CONJG(S) = 1.0.
! 	*
! 	* =====================================================================
! 	*
! 	*     .. Local Scalars ..
	      INTEGER            I, IX, IY
	      COMPLEX            STEMP
! 	*     ..
! 	*     .. Intrinsic Functions ..
	      INTRINSIC          CONJG
! 	*     ..
! 	*     .. Executable Statements ..
! 	*
	      IF( N.LE.0 ) &
	     &   RETURN
	      IF( INCX.EQ.1 .AND. INCY.EQ.1 ) &
	     &   GO TO 20
! 	*
! 	*     Code for unequal increments or equal increments not equal to 1
! 	*
	      IX = 1
	      IY = 1
	      IF( INCX.LT.0 ) &
	     &   IX = ( -N+1 )*INCX + 1
	      IF( INCY.LT.0 ) &
	     &   IY = ( -N+1 )*INCY + 1
	      DO 10 I = 1, N
	         STEMP = C*CX( IX ) + S*CY( IY )
	         CY( IY ) = C*CY( IY ) - CONJG( S )*CX( IX )
	         CX( IX ) = STEMP
	         IX = IX + INCX
	         IY = IY + INCY
	   10 CONTINUE
	      RETURN
! 	*
! 	*     Code for both increments equal to 1
! 	*
	   20 CONTINUE
	      DO 30 I = 1, N
	         STEMP = C*CX( I ) + S*CY( I )
	         CY( I ) = C*CY( I ) - CONJG( S )*CX( I )
	         CX( I ) = STEMP
	   30 CONTINUE
	      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------

END MODULE 