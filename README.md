# Odemodule

Odemodule is a Fortran 2003 module for solving ODE's. It utilizes the fourth-order Runge-Kutta with Cash-Karp implementation as the adaptive step size.

The module uses whole array operations and linked-lists for efficient storage of the intermediate results. It needs a Fortran 2003 compliant compiler such as `gfortran`.

# How to use

In your program define two allocatable arrays like this:
```fortran
   USE odemodule
   REAL(KIND=k), ALLOCATABLE, DIMENSION (:)   :: x_result
   REAL(KIND=k), ALLOCATABLE, DIMENSION (:,:) :: y_result
```
Do not allocate these in your main program. We will feed them into the solver and they will be allocated there.

The derivatives of the system states should be calculated in a subroutine with an interface as follows:
```fortran
   SUBROUTINE dydx (x,y,dy)
      REAL(k),INTENT(in) :: x
      REAL(k),INTENT(in), DIMENSION(:) :: y
      REAL(k),INTENT(out), DIMENSION(SIZE(y)) :: dy
```
Finally, you can get the solution by:

```fortran
CALL odesolve ( x0 , xend , y0 , dydx , x_result , y_result )
```
Where `x0` and `xend` are the start and end points (the independent variable) and `y0` would be the vector of initial conditions.
`odesolve` will fill the allocatable arrays i.e., `x_result` and `y_result` with all the points.
If you would like, you can define the desired accuracy (default is 1.e-6) and the maximum stepsize (default is 0.1) and set them in the solver as well: 
```yaml
CALL odesolve ( x0 , xend , y0 , dydx , x_result , y_result , 1.e-9_k , 0.01_k)
```

# Examples
## Solving a single input, single output ODE
![equation](https://latex.codecogs.com/gif.latex?\dpi{130}\frac{dy}{dx}=cos(x))

Even with one input and one output, the interface of the derivative subroutine must hold. So, we will create `dy` as an array with one element. In the main program, we will define the initial condition and will call `odesolve`.
The whole program would look like this:
 ```fortran
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
```
## The double pendulum 

[The double pendulum](https://www.astro.umd.edu/~adhabal/V1/Reports/Order_and_Chaos.pdf) will have 2 degrees of freedom and 4 elements as the derivatives. Based on the initial conditions, this system can be completely chaotic and impossible to predict.

The derivation subroutine would be:
```fortran
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
```
and the main program where we have set `1.e-9` for the accuracy and a finer maximum step of `0.01`:
```yaml
   USE odemodule
   IMPLICIT NONE
   REAL(KIND=k), ALLOCATABLE, DIMENSION (:,:) :: yf
   REAL(KIND=k), ALLOCATABLE, DIMENSION (:) :: xf
   REAL(KIND=k), PARAMETER :: pi=ACOS(-1._k)
   REAL(kind=k) :: x0 , xend , y0(4)

   x0 = 0._k
   xend = 20._k

   ! y(1)=θ1
   ! y(2)=d(θ1)/dt
   ! y(3)=θ2
   ! y(4)=d(θ2)/dt
   y0 = [ 2*pi/3 , 0._k , pi/2 , 0._k ]

   CALL odesolve ( x0 , xend , y0 , dydx , xf , yf , 1.e-9_k , 0.01_k )
```
The trace of the pendulum end:
![The chaotic pendulum](/double_pendulum.gif)

The full programs have been supplied in the repository.
# Suggestions and bug submit

If you have any suggestions or have found a bug, please let me know.
