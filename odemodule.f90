!#############################################################################
MODULE odemodule
    IMPLICIT NONE
    INTEGER, PARAMETER :: k=SELECTED_REAL_KIND(14)
    PRIVATE :: rk4,rkqcck,rkck
    ABSTRACT INTERFACE
       SUBROUTINE derivation (x,y,dy)
          IMPORT k
          REAL(KIND=k), INTENT(IN), DIMENSION(:)        :: y
          REAL(KIND=k), INTENT(IN)                      :: x
          REAL(KIND=k), INTENT(OUT), DIMENSION(SIZE(y)) :: dy
       END SUBROUTINE derivation
    END INTERFACE
CONTAINS
    SUBROUTINE rk4(x,y,dy,h,yout,dydx)
        PROCEDURE(derivation) :: dydx
        REAL(KIND=k), INTENT(IN), DIMENSION(:) :: y,dy
        REAL(KIND=k), INTENT(OUT), DIMENSION(:) :: yout
        REAL(KIND=k), INTENT(IN) :: h,x
        REAL(KIND=k) :: hh,h6
        REAL(KIND=k) :: xh
        REAL(KIND=k), DIMENSION (SIZE(y)) :: yt,dyt,dym
        hh=0.5_k*h
        h6=h/6._k
        xh=x+hh
        yt=y+hh*dy
        CALL dydx(xh,yt,dyt) ! to produce k2
        yt=y+hh*dyt
        CALL dydx(xh,yt,dym) ! to produce k3
        yt=y+h*dym
        dym=dym+dyt
        CALL dydx(x+h,yt,dyt) ! to produce k4
        yout=y+h6*(dy+dyt+2._k*dym)
    END SUBROUTINE rk4
    !#########################################################################
    SUBROUTINE rkck(x,y,dy,h,yout,yerr,dydx)
        !Cash-Karp Runge-Kutta step
        PROCEDURE(derivation) :: dydx
        REAL(KIND=k), INTENT(IN) :: h,x
        REAL(KIND=k), INTENT(IN), DIMENSION(:) :: dy,y
        REAL(KIND=k), INTENT(OUT), DIMENSION(:) :: yerr,yout
        REAL(KIND=k), DIMENSION (SIZE(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
        REAL(KIND=k), PARAMETER ::  a2=.2_k , a3=.3_k , a4=.6_k , a5=1._k , a6=.875_k , &
            b21=.2_k , b31=3._k/40._k , b32=9._k/40._k , b41=.3_k , b42=-.9_k , &
            b43=1.2_k , b51=-11._k/54._k , b52=2.5_k , b53=-70._k/27._k , &
            b54=35._k/27._k , b61=1631._k/55296._k , b62=175._k/512._k , &
            b63=575._k/13824._k , b64=44275._k/110592._k , b65=253._k/4096._k , &
            c1=37._k/378._k , c3=250._k/621._k , c4=125._k/594._k , c6=512._k/1771._k , &
            dc1=c1-2825._k/27648._k , dc3=c3-18575._k/48384._k , &
            dc4=c4-13525._k/55296._k , dc5=-277._k/14336._k , dc6=c6-.25_k
        ytemp=y+b21*h*dy
        CALL dydx(x+a2*h,ytemp,ak2)
        ytemp=y+h*(b31*dy+b32*ak2)
        CALL dydx(x+a3*h,ytemp,ak3)
        ytemp=y+h*(b41*dy+b42*ak2+b43*ak3)
        CALL dydx(x+a4*h,ytemp,ak4)
        ytemp=y+h*(b51*dy+b52*ak2+b53*ak3+b54*ak4)
        CALL dydx(x+a5*h,ytemp,ak5)
        ytemp=y+h*(b61*dy+b62*ak2+b63*ak3+b64*ak4+b65*ak5)
        CALL dydx(x+a6*h,ytemp,ak6)
        ! accumulate increments with proper weights.
        yout=y+h*(c1*dy+c3*ak3+c4*ak4+c6*ak6)
        ! estimate error as difference between fourth and fifth order methods.
        yerr=h*(dc1*dy+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)
    END SUBROUTINE rkck
    !##########################################################################
    SUBROUTINE rkqcck(y,dy,x,htry,eps,hdid,hnext,dydx)
        ! A subroutine to determine optimum h by Cash-Karp method
        PROCEDURE(derivation) :: dydx
        REAL(KIND=k), INTENT (IN), DIMENSION(:) :: dy
        REAL(KIND=k), INTENT(INOUT), DIMENSION(:) :: y
        REAL(KIND=k), INTENT(IN) :: htry,eps
        REAL(KIND=k), INTENT(INOUT) :: x
        REAL(KIND=k), INTENT(OUT) :: hdid,hnext
        REAL(KIND=k),PARAMETER :: pgrow=-0.2_k , pshrink=-0.25_k , epss=1.e-30_k , &
             safety=0.9_k , errconv=EXP((1._k/pgrow)*LOG(5._k/safety))
        REAL(KIND=k), DIMENSION (SIZE(y)) :: ytemp,yerr,yscal
        REAL(KIND=k) :: h,errmax,xnew,htemp
        h=htry
        DO
            CALL rkck(x,y,dy,h,ytemp,yerr,dydx)
            errmax=0._k
            yscal=ABS(y)+ABS(h*dy)+epss
            errmax=MAXVAL(ABS(yerr/yscal))/eps
            IF (errmax>1._k) THEN
                htemp=safety*h*(errmax**pshrink)
                h=SIGN(MAX(ABS(htemp),0.1_k*ABS(h)),h)
                xnew=x+h
                IF (xnew.eq.x) THEN
                    PRINT*,'Step got too small to be significant at x=', x
                    !h=0.001_k
                    !exit
                ENDIF
            ELSE
                EXIT
            ENDIF
        ENDDO
        hdid=h
        x=x+hdid
        IF (errmax>errconv) THEN
            hnext=safety*h*(errmax**pgrow)
        ELSE
            hnext=5._k*h
        ENDIF
        y=ytemp
    END SUBROUTINE rkqcck
    !##########################################################################
    SUBROUTINE odesolve(x0,xend,y0,dydx,x,y,ep,hm)
        REAL(KIND=k), INTENT(IN), DIMENSION(:) :: y0
        REAL(KIND=k), INTENT(IN) :: x0,xend
        REAL(KIND=k), INTENT(IN), OPTIONAL :: ep , hm
        REAL(KIND=k), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: y
        REAL(KIND=k), INTENT(OUT), DIMENSION(:), ALLOCATABLE :: x
        PROCEDURE (derivation) :: dydx
        INTEGER :: n,i,istat
        LOGICAL :: send
        REAL(KIND=k) :: htry,hdid,hnext,xstep,eps,hmax
        REAL(KIND=k), DIMENSION (SIZE(y0)) :: ystep,dystep
        TYPE :: xlist
            REAL(KIND=k) :: xval
            TYPE(xlist), POINTER :: xvalp
        END TYPE xlist
        TYPE :: ylist
            REAL(KIND=k), DIMENSION(:), ALLOCATABLE :: yval
            TYPE(ylist), POINTER :: yvalp
        END TYPE ylist
        TYPE (xlist), POINTER :: xtemp,xhead,xtail
        TYPE (ylist), POINTER :: ytemp,yhead,ytail

        IF (PRESENT(ep)) THEN
            eps=ep
        ELSE
            eps=1.e-6_k
        ENDIF

        IF (PRESENT(hm)) THEN
            hmax=hm
        ELSE
            hmax=0.1_k
        ENDIF
        
        send=.FALSE.
        htry=SIGN(1.e-4_k,xend-x0)
        n=SIZE(y0)
        
        ! Allocating the head and tail of xlist with only the initial conditions (x0).
        ALLOCATE (xhead,STAT=istat);IF (istat.ne.0) PRINT*,'Allocation of xhead failed.'
        xtail=>xhead
        NULLIFY(xhead%xvalp)
        xhead%xval=x0
        
        ! Allocating the head and tail of ylist with only initial conditions (y0).
        ALLOCATE (yhead,STAT=istat);IF (istat.ne.0) PRINT*,'Allocation of yhead failed.'
        ytail=>yhead
        NULLIFY(yhead%yvalp)
        ALLOCATE(yhead%yval(n),STAT=istat);
        IF (istat.NE.0) PRINT*,'Allocation of yhead%yval failed.'
        yhead%yval=y0

        ystep=y0
        xstep=x0
        CALL dydx (xstep,ystep,dystep)
        
        i=1
        DO WHILE (.NOT. send)
            i=i+1
            CALL rkqcck(ystep,dystep,xstep,htry,eps,hdid,hnext,dydx)
            IF (ABS(xstep) .GT. ABS(xend)) THEN
                CALL rk4(xtail%xval,ytail%yval,dystep,xend-xtail%xval,ystep,dydx)
                xstep=xend
                send=.TRUE.
             ENDIF
            ! Adding the new x to the tail of xlist.
            ALLOCATE(xtail%xvalp,STAT=istat);
            IF (istat.NE.0) PRINT*,'Allocation of new xtail failed.'

            xtail=>xtail%xvalp
            xtail%xval=xstep
            NULLIFY(xtail%xvalp)

            ! Adding the new y to the end of ylist.
            ALLOCATE(ytail%yvalp,STAT=istat);
            IF (istat.NE.0) PRINT*,'Allocation of new ytail failed.'

            ytail=>ytail%yvalp
            ALLOCATE(ytail%yval(n),STAT=istat);
            IF (istat.NE.0) PRINT*,'Allocation of new ytail%yval failed.'
            
            ytail%yval=ystep
            NULLIFY(ytail%yvalp)

            ! Producing dy for the next round.
            CALL dydx(xstep,ystep,dystep)
            htry=MIN(hnext,hmax)
        ENDDO

        ALLOCATE (y(n,i),x(i),STAT=istat);
        IF (istat.NE.0) PRINT*,'Allocation of x & y failed.'

        xtemp=>xhead
        ytemp=>yhead
        i=1
        ! Filling output x & y with the xlist & ylist
        DO WHILE (ASSOCIATED(xtemp))
            x(i)=xtemp%xval
            y(:,i)=ytemp%yval
            xtemp=>xtemp%xvalp
            ytemp=>ytemp%yvalp
            i=i+1
        ENDDO
    END SUBROUTINE odesolve
    !#######################################################################
END MODULE odemodule
