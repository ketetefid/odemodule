PROGRAM testode
   USE odemodule
   IMPLICIT NONE
   REAL(KIND=k), ALLOCATABLE, DIMENSION (:,:) :: yf
   REAL(KIND=k), ALLOCATABLE, DIMENSION (:) :: xf
   REAL(KIND=k), PARAMETER :: pi=ACOS(-1._k)
   INTEGER :: i
   REAL(kind=k) :: x0 , xend , y0(1) 

   x0 = 0._k      ! the starting point
   xend = 2*pi    ! the final point
   y0 = [ 0._k ]  ! the vector of initial conditions

   CALL odesolve ( x0 , xend , y0 , dydx , xf , yf )

   ! check the result:  dy/dx=cos(x) => y=sin(x)
   DO i=1,SIZE(xf)
      WRITE(*,*) ABS(SIN(xf(i))-yf(1,i))
   ENDDO
   
CONTAINS
   SUBROUTINE dydx (x,y,dy)
      REAL(k) , INTENT(in) :: x
      REAL(k) , INTENT(in) , DIMENSION(:) :: y
      REAL(k) , INTENT(out) , DIMENSION(SIZE(y)) :: dy
      
      dy(1)=COS(x)

   END SUBROUTINE dydx
END PROGRAM testode
      
