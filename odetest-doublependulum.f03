PROGRAM testode
   USE odemodule
   IMPLICIT NONE
   REAL(KIND=k), ALLOCATABLE, DIMENSION (:,:) :: yf
   REAL(KIND=k), ALLOCATABLE, DIMENSION (:) :: xf
   REAL(KIND=k), PARAMETER :: pi=ACOS(-1._k)
   INTEGER :: i
   REAL(kind=k) :: x0 , xend , y0(4)
   REAL(kind=k) :: th1 , th2
   

   x0 = 0._k
   xend = 20._k

   ! y(1)=θ1
   ! y(2)=d(θ1)/dt
   ! y(3)=θ2
   ! y(4)=d(θ2)/dt

   th1 = 2*pi/3
   th2 = pi/2
   
   y0 = [ th1 , 0._k , th2 , 0._k ]

   CALL odesolve ( x0 , xend , y0 , dydx , xf , yf , 1.e-9_k , 0.01_k )
   PRINT*, 'Number of nodes:',SIZE(xf)

   OPEN (1, file='results')
   DO i=1,SIZE(xf)
      WRITE(1,*) SIN(yf(1,i))+SIN(yf(3,i)),-COS(yf(1,i))-cos(yf(3,i))
   ENDDO
   
   !CALL system("gnuplot -p -e "//'"plot '//"'results'"// ' with lines'//'"')
   CALL system("gnuplot doublep.gp")
   
CONTAINS
   SUBROUTINE dydx (x,y,dy)
      REAL(k),INTENT(in) :: x
      REAL(k),INTENT(in), DIMENSION(:) :: y
      REAL(k),INTENT(out), DIMENSION(SIZE(y)) :: dy
      REAL(k) :: c1 , c2 , c3 , c4 , c5
      ! y(1)=θ1
      ! y(2)=d(θ1)/dt
      ! y(3)=θ2
      ! y(4)=d(θ2)/dt

      c1 = 2.16_k
      c2 = 5._k/18._k
      c3 = 1._k
      c4 = 172.5_k
      c5 = 588._k/18._k
      
      dy(1)=y(2)
      dy(2)=(2*c2*c4*SIN(y(1))+(c3*y(2))**2*SIN(y(1)-y(3))*COS(y(1)-y(3))&
           +2*c2*c3*y(4)**2*SIN(y(1)-y(3))-c3*c5*COS(y(1)-y(3))*SIN(y(3)))&
           /((c3*COS(y(1)-y(3)))**2-4*c1*c2)

      dy(3)=y(4)
      dy(4)=(2*c1*c5*SIN(y(3))-(c3*y(4))**2*SIN(y(1)-y(3))*COS(y(1)-y(3))&
           -2*c1*c3*y(2)**2*SIN(y(1)-y(3))-c3*c4*COS(y(1)-y(3))*SIN(y(1)))&
           /((c3*COS(y(1)-y(3)))**2-4*c1*c2)

   END SUBROUTINE dydx
END PROGRAM testode
      
