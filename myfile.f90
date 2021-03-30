program addNumbers

! This simple program adds two numbers
   implicit none

! Type declarations
   real :: result
   real :: A, B, time, total_time, twist, rho, F_a_0, F_a_1, g_0, g_1, r, delta_t, delta_r, Radius
   real(kind = 8) :: EI_00, EI_01, EI_10, EI_11
   integer, parameter :: N = 100
   integer, parameter :: T = 10000000
   integer :: i,j


 ! Dynamic velocity arrays:
   real, dimension(N+2) :: V_0_new
   real, dimension(N+2) :: V_1_new
   real, dimension(N+2) :: V_0_old
   real, dimension(N+2) :: V_1_old
   
 ! Moments arrays:
   real, dimension(N+2) :: M_0_new
   real, dimension(N+2) :: M_1_new
   real, dimension(N+2) :: M_0_old
   real, dimension(N+2) :: M_1_old
   

 ! Dynamic displacement arrays:
   real, dimension(N+2) :: q_0_new
   real, dimension(N+2) :: q_1_new
   real, dimension(N+2) :: q_0_old
   real, dimension(N+2) :: q_1_old

 ! Radial position in the blade array:
   real, dimension(N+2) :: Radius_position

 ! Time position in the blade array:
   real, dimension(T+1) :: Time_position

 ! Root bending moments array:
   real, dimension(T+1) :: Root_M_0
   real, dimension(T+1) :: Root_M_1
   
 ! Tip velocity of the displacements array:
   real, dimension(T+1) :: Tip_V_0
   real, dimension(T+1) :: Tip_V_1
  
 ! Tip displacements array:
   real, dimension(T+1) :: Tip_q_0
   real, dimension(T+1) :: Tip_q_1

! Executable statements
 ! Get all arrays with initial value 0:
   V_0_new = 0.0
   V_1_new = 0.0
   V_0_old = 0.0
   V_1_old = 0.0
   M_0_new = 0.0
   M_1_new = 0.0
   M_0_old = 0.0
   M_1_old = 0.0
   q_0_new = 0.0
   q_1_new = 0.0
   q_0_old = 0.0
   q_1_old = 0.0
   Radius_position = 0.0
   Time_position = 0.0
   Root_M_0 = 0.0
   Root_M_1 = 0.0
   Tip_V_0 = 0.0
   Tip_V_1 = 0.0
   Tip_q_0 = 0.0
   Tip_q_1 = 0.0

 ! Main program
   total_time = 1.275
   time = 0
   A = -100.0 ! Amplitud aerodynamic force direction 1
   B = 0.0 ! Amplitud aerodynamic force direction 0
   
   ! Properties DEFINITION:
   EI_00 = (10**9)*1.5
   EI_11 = (10**9)*1.5
   EI_01 = 0
   EI_10 = 0
   twist = 0
   rho = 500.0
   F_a_0 = B
   F_a_1 = A
   g_0 = 0.0
   g_1 = 0.0
   Radius = 62.5
   delta_r = Radius/(N-1)

   ! Program loops iteration:

   do j = 1, T
      
      delta_t = total_time/T
      r = 0



      do i = 2, N+1
         ! Calculation of the velocity new values from the old step time:
         V_0_new(i) = V_0_old(i) + delta_t*((-1/rho)*((M_0_old(i+1)-(2*M_0_old(i))+M_0_old(i-1))/delta_r**2)+(F_a_0/rho)+ g_0)
         V_1_new(i) = V_1_old(i) + delta_t*((-1/rho)*((M_1_old(i+1)-(2*M_1_old(i))+M_1_old(i-1))/delta_r**2)+(F_a_1/rho)+ g_1)
         
         ! Boundary conditions at the ROOT position: 
         V_0_new(1) = 0
         V_1_new(1) = 0
         V_0_new(2) = 0
         V_1_new(2) = 0

         ! Calculation of the bending moment from the new velocity values:
         M_0_new(i) = M_0_old(i) + delta_t*((EI_00*((V_0_new(i+1)-(2*V_0_new(i))+V_0_new(i-1))/delta_r**2)))! +EI_01*((V_1_new[i+1]-(2*V_1_new[i])+V_1_new[i-1])/delta_r**2)))
         M_1_new(i) = M_1_old(i) + delta_t*((EI_11*((V_1_new(i+1)-(2*V_1_new(i))+V_1_new(i-1))/delta_r**2)))! +EI_10*((V_0_new[i+1]-(2*V_0_new[i])+V_0_new[i-1])/delta_r**2)))

         ! Boundary conditions at the TIPS position:
         M_0_new(N+2) = 0
         M_1_new(N+2) = 0
         M_0_new(N+1) = 0
         M_1_new(N+1) = 0

         ! Calculation of the deformation value from the velocity calculated:
         q_0_new(i) = q_0_old(i) + delta_t*((V_0_old(i)))
         q_1_new(i) = q_1_old(i) + delta_t*((V_1_old(i)))

         ! Position control:
         Radius_position(i) = r
         r = r + delta_r
      
      end do

      ! Upload data to new time-step:
       V_0_old = V_0_new 
       V_1_old = V_1_new 
       M_0_old = M_0_new 
       M_1_old = M_1_new 
       q_0_old = q_0_new
       q_1_old = q_1_new

      ! Store data to plot: 
       Root_M_0(j) = M_0_new(2)
       Root_M_1(j) = M_1_new(2)
       Tip_V_0(j) = V_0_new(N+1)
       Tip_V_1(j) = V_0_new(N+1)
       Tip_q_0(j) = q_0_new(N+1)
       Tip_q_1(j) = q_1_new(N+1)

      ! Time control:
       Time_position(j) = time
       time = time + delta_t

   end do

   print *, huge(EI_00)
   print *, 'The  ', M_1_new

OPEN(1, FILE='Tip_q_0.csv', FORM='formatted')
Do i = 1, T
   write(1, '(10F5.3)') Tip_q_0(i)
End do

OPEN(1, FILE='Tip_q_1.csv', FORM='formatted')
Do i = 1, T
   write(1, '(10F5.3)') Tip_q_1(i)
End do

OPEN(1, FILE='Root_M_0.csv', FORM='formatted')
Do i = 1, T
   write(1, '(10F5.3)') Root_M_0(i)
End do

OPEN(1, FILE='Root_M_1.csv', FORM='formatted')
Do i = 1, T
   write(1, '(10F5.3)') Root_M_1(i)
End do

OPEN(1, FILE='Radius_position.csv', FORM='formatted')
Do i = 1, N
   write(1, '(10F5.3)') Radius_position(i)
End do

OPEN(1, FILE='Time_position.csv', FORM='formatted')
Do i = 1, T
   write(1, '(10F5.3)') Time_position(i)
End do

OPEN(1, FILE='q_0_new.csv', FORM='formatted')
Do i = 1, N
   write(1, '(10F5.3)') q_0_new(i)
End do

OPEN(1, FILE='q_1_new.csv', FORM='formatted')
Do i = 1, N
   write(1, '(10F5.3)') q_1_new(i)
End do

OPEN(1, FILE='M_0_new.csv', FORM='formatted')
Do i = 1, N
   write(1, '(10F5.3)') M_0_new(i)
End do

OPEN(1, FILE='M_1_new.csv', FORM='formatted')
Do i = 1, N
   write(1, '(10F5.3)') M_1_new(i)
End do



end program addNumbers