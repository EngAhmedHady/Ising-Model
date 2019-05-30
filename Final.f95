! A fortran95 program for G95
program main
implicit none
   integer, dimension (:,:), allocatable :: Lattice
   integer  :: l = 1 ! ------------------------ Lattice Linear Size
   real     :: N ! ---------------------------- Lattice Size
   integer  :: i, j, o,k, Counter ! ---- Loops Prameters
   integer  :: IntRand(2) !-------------------- Random selected element
   real     :: Temp(2) ! ---------------------- Temperature
   real     :: dT !---------------------------- Temperature step
   real     :: T !----------------------------- Loop Temperature
   integer  :: IC ! --------------------------- Intial Condition
   real     :: MeanEn, EnSpin   = 0  ! -------- Energy Retioes 
   real     :: MeanMag, MagSpin =  0 ! -------- Magnetization Retioes 
   Integer  :: Mag, En ! ---------------------- Total Magnetization & Energy 
   real     :: MeanM2, MeanE2, MB
   integer, dimension(:), allocatable :: UN !-- UpdateNeighbors
   integer  :: Itt,Itt2, Ittdy ! -------------- Main Metropolis Itterations
   integer  :: TStep ! ------------------------ Number of Temperature Steps
   real     :: start_time, stop_time
   logical  :: file_exists, TempCheck, SoluSign, Stability
   integer  :: NowDate(3), Nowtime(3)
   character(len = 120):: S ! -------- Ploting String
   character(len = 5) :: FileName,NIC
   character(len = 150) :: PlotText
   
   
do while (l > 0)

   print *, "================ MONTE CARLO SIMULATION OF THE 2D ISING MODEL ================="
   call DateAndTime(NowDate,Nowtime)
    ! ========================= Check Existance And Open Data Files ===========================|

     INQUIRE(FILE="TimeDependant.dat", EXIST=file_exists) ! ------------ Transiant state file data
     if (file_exists) then
         open(1, file = "TimeDependant.dat" , status = 'old') 
     else 
         open(1, file = "TimeDependant.dat" , status = 'new')
     end if
     
     INQUIRE(FILE="TempDependant.dat", EXIST=file_exists)! ------------- Temperature state file data
     if (file_exists) then
         open(2, file = "TempDependant.dat" , status = 'old') 
     else 
         open(2, file = "TempDependant.dat" , status = 'new')  
     end if
     
     INQUIRE(FILE="Lattice.dat", EXIST=file_exists)! ------------------- Lattice Matrix file
     if (file_exists) then
         open(3, file = "Lattice.dat" , status = 'old') 
     else 
         open(3, file = "Lattice.dat" , status = 'new') 
     end if
     close(3) 
   ! == Step 0 ========================= Lattice Builde ======================================|
   !                                  Initiate The Lattic Dimention                           |
   !==========================================================================================|
     
     Lattice = Initiate(l , IC, Mag, En)
   
   ! ============================== Advanced Initiate Options ================================|
     Itt = l**3
     Itt2 = l**4/Itt
     TStep = 60
     T = Temp(1)
     dT = (Temp(2)-Temp(1))/TStep
     
     print 25, " -------------------| Starting Metroplis Equelibrium |--------------------------"
     25 format(/,a)
     call DateAndTime(NowDate,Nowtime)
    
    ! ================================= The Equilibruim Loop =================================|
     k = 0
     Ittdy = Itt
     Stability = .True.
     TempCheck = .FALSE.

     call cpu_time(start_time)
     
     do o =1,Itt2
         call Metropolis (l,Lattice,T,TempCheck,MeanEn, MeanMag,MeanE2,MeanM2,Counter)
         MagSpin = real(Mag) / N
         EnSpin  = real(En)  / N
        
         write(1,*) o*Ittdy, MagSpin, EnSpin
         if (abs(MagSpin) > 0.999 .and. abs(MagSpin) < 1.001 ) then
             IF(MOD(o,5).EQ.0) WRITE(6,42,ADVANCE='NO') "Loading: ", (real(o)/(Itt2))*100, ' %',"Stable  "//CHAR(13)
             42  FORMAT(A, f 6.2, A,10x,A)
             If (Stability) then
             Itt = o * Ittdy
             end if
             Stability = .False.
         else
             IF(MOD(o,5).EQ.0) WRITE(6,42,ADVANCE='NO') "Loading: ", (real(o)/(Itt2))*100, ' %',"Unstable"//CHAR(13)
         End if
     end do
     close(1) 
     call system('gnuplot -p TransientResults.gnu')

     if (MagSpin > 0) then
         SoluSign = .TRUE.
     else
         SoluSign = .FALSE.
     end if
     
    ! =========================================================================================|
     call cpu_time(stop_time)
     print *, "------------------------------------------------------------------------------"
     print 12, "Selected Energy/spin : ", EnSpin, "Selected Mag/spin: ", MagSpin
     12 format(1x,a,2x,f6.3,7x,a,2x,f6.3)
     print 13, "CPU Time (Ts)        : ", stop_time - start_time, "Seconds"
     13 format(1x,a,e12.3,3x,a)
     print 14, "Stability Loops      : ", Itt 
     14 format(1x,a,2x,i10)
     print 16, "Temperature step (ΔT): ", dT
     16 format(1x,a,2x,f6.3)
     if ((T < 1.5) .or. (T > 2.5)) Then 
         if (abs(MagSpin) > 0.999 .and. abs(MagSpin) < 1.001 ) then
             print *, "Simulation Status    :    Stable"
         Else
             print *, "Simulation Status    :    UnStable"
         End if
     End if
     print *, "------------------------------------------------------------------------------"
     call DateAndTime(NowDate,Nowtime)
   ! ==========================================================================================|
     print *, "==============================================================================="
     print *, "Temperature (T) | Mean Mag.<m> | Mean Energy<e> | S.Mean Mag<m²> |S.Mean E <e²>"
     print *, "----------------|--------------|----------------|----------------|-------------"
   ! ====================================== Temperature Turn ==================================|
     call cpu_time(start_time)
     do o = 1, (TStep + 1)
     
         Ittdy = Itt
         if (T > 2.2 .and. T < 2.4) then 
             Ittdy = 5*Itt
         else
             Ittdy = int(N)
         End If 
         
         call Metropolis (l,Lattice,T,TempCheck,MeanEn, MeanMag,MeanE2,MeanM2,Counter)
         
         Ittdy = l**4
         TempCheck = .TRUE.  
             call Metropolis (l,Lattice,T,TempCheck,MeanEn, MeanMag,MeanE2,MeanM2,Counter)

         call TempFunctions(T, En, Mag, MeanM2, MeanMag, MeanE2, MeanEn)
         
         if (mod(o,(TStep+1)/10)==0) then
             print 300, T, MeanMag/real(Counter), MeanEn/real(Counter), MeanM2/real(Counter), MeanE2/real(Counter)
             300 format(6x,f6.4,5x,'|', f10.3, 4x, '|',2x, f10.3, 4x, '|', f10.3,6x, '|',f10.3) 
         end if 
         T = T + dT
         TempCheck = .FALSE.
        !============================ Progress Indecator ====================================
         IF(MOD(o,5).EQ.0) WRITE(6,41,ADVANCE='NO') "Loading: ", (real(o)/(TStep + 1))*100, ' %'//CHAR(13)
         41  FORMAT(A, f 6.2, A)
        !====================================================================================
       
   end do
   close(1) 
   close(2)
   call system('gnuplot -p TemperatureResults.gnu')
   call cpu_time(stop_time)
   print *, "------------------------------------------------------------------------------"
   print *, "Total Ts of teperature loops: ", (stop_time-start_time)/60, "Min"
   print *, "------------------------------------------------------------------------------"
end do 



! =============================================================================|
!                              Functions In Use                                |
! =============================================================================|
contains
function Initiate(l , IC, Mag, En)
     integer , dimension (:,:), allocatable :: Initiate
     integer  :: l! ------------------------ Lattice Linear Size
     integer  :: IC ! ---------------------- Intial Condition
     logical  :: RIC ! --------------------- Random Intial Condition
     integer,intent (out):: Mag !----------- Intial Total Magnetization
     integer,intent (out):: En ! ----------- Intial Total Energy
     real     :: Rnd !---------------------- Random Number
     print *, "Enter the Linear size of the Lattice (L): "
     read  *, l ! -------------------------- Lattic size input
      
     if (l /= 0) then
         N = real(l)**2
         print *, "Enter Intial condition" 
         print *, "1- ⇈     2- ⇊      3- ⇅"
         read  *, IC
         print *, "Enter Temperature Range  (T): "
         read  *, Temp ! ------------------ Lattic size input
     End If
     
     write (NIC,'(I3)') IC
     Mag = 0 ! ---------------------------- Initiate Magnetization For New Lattice
     En = 0 ! ----------------------------- Initiate Energy For New Lattice

     if (IC == 1) then
         IC  = 1
         RIC = .false.
     else if (IC == 2)then
         IC  = -1
         RIC = .false.
     else 
         RIC = .true.
     end if
     allocate (Initiate(l,l))
      
     ! ######## Initiate New clear Lattice #######
      do j=1,l 
         do i=1,l 
             Initiate(i,j) = 0 
         end do 
      end do 
     ! ###########################################
      
     !$OMP DO
     
     do i=1,l
         do j = 1, l
             If (RIC) then
                  call random_number(Rnd)
                  if (Rnd > 0.5) then
                     IC = 1
                  else
                     IC = -1
                  End If
              end if
             Initiate(i, j) = IC
         end do
     end do
     
     do i = 1,l
         do j = 1,l
             Mag  = Mag + Initiate(i,j)
             Un = Neighbors(l , i , j)
             En  = En - Initiate(i,j) * (Initiate(UN(2), j)+Initiate(i, UN(4)))
         end do
     end do
     
     !$OMP END DO

     MagSpin = Mag/N
     EnSpin  = En /N
     
     print 406, "Intial Lattic Magnetization (M):", Mag
     print 405, "Intial Magnetization/spin (m)  :", MagSpin
     print 406, "Intial Lattic Energy (H)       :", En
     print 405, "Intial Energy/spin (e)         :", EnSpin
     print 405, "Intial Temperature (T)         :", Temp(1)
     405 format(a,5x,f10.3)
     406 format(a,5x,i10)
 end function Initiate
! ==========================================================================

   subroutine Metropolis(l, W, T,Check, eM, mM, E2, M2, Counter) 
     integer  :: l, m, u
     integer, intent(out) :: Counter
     real, intent(out)    :: eM, mM, E2, M2
     integer, intent (inout), dimension (:,:) :: W ! ---Lattice Input and Output
     real,    intent (in)     :: T ! ------------------ Temperature
     integer  :: dEnergy = 0! ------------------------- NewEnergy, EnergyDifference = 0
     real     :: r(2) ! ------------------------------- Random Number
     real     :: PiPr !-------------------------------- Pi Probability
     real     :: Rnd !--------------------------------- Random Number
     logical  :: Check !------------------------------- If Metropolis works under temperature loop
     character (len = 2) :: SaveData ! ---------------- Check to save data or not
     eM = 0.0; mM = 0.0;
     M2 = 0.0; E2 = 0.0;
     counter = 0;
     If ((((T > 1.9999).and.(T < 2.0001)).or. ((T > 2.499).and.(T < 2.5001))) .and. (Check .eqv. .False.)) then
         open(1, file = "TimeDependant.dat" , status = 'old') 
     end if
     !$OMP DO Parallel
     do u = 1, Ittdy
         call random_number(r)
         IntRand = int(r * l) + 1 !------------------- Selecting Random Element in side l
         UN = Neighbors(l , IntRand(1), IntRand(2))
         
         ! =|1|= Local Energy Calculation Before Change ============
         dEnergy = (W(UN(1),IntRand(2)) + W(UN(2),IntRand(2)) + W(IntRand(1),UN(3)) + W(IntRand(1),UN(4)))
         dEnergy =   2 * dEnergy* W(IntRand(1),IntRand(2))
         PiPr = exp(-(1/T) * dEnergy)
         call random_number(Rnd)
         ! ######## Energy Change Acceptance ########
         If ((dEnergy <= 0) .or. (PiPr >= Rnd)) Then
             W(IntRand(1),IntRand(2)) = -1 * W(IntRand(1),IntRand(2))
             Mag = Mag + 2 * W(IntRand(1),IntRand(2))
             En = En + dEnergy
         end if
         
         if (mod(u,Ittdy/5) == 0 .and. k < 90) then
             open(3, file = "Lattice.dat" , status = 'old') 
             do m = lbound(W,1), ubound(W,1)
                 write(3,*) (W(m,j), j = lbound(W,2), ubound(W,2))
             end do
             close(3) 
             write (FileName,'(I3)') k
             s ='gnuplot -e "splot ''Lattice.dat'' matrix; set palette negative; set pm3d map; set term png size 1000,1000;'
             PlotText = trim(s)//'set key off;set output'''//trim(NIC)//trim(FileName)//'.png''; replot;"'
             call system(PlotText)
             k = k + 1
         end if
         If ((((T > 1.999).and.(T < 2.001)).or. ((T > 2.499).and.(T < 2.501))) .and. (Check .eqv. .False.)) then
             If (mod(u,int(Ittdy/N)) == 0) then
                 MagSpin = real(Mag) / N
                 EnSpin  = real(En)  / N
                 write(1,*) u, MagSpin, EnSpin
                 if (mod(u,int(Ittdy/45)) == 0) then
                     open(3, file = "Lattice.dat" , status = 'old') 
                     do m = lbound(W,1), ubound(W,1)
                         write(3,*) (W(m,j), j = lbound(W,2), ubound(W,2))
                     end do
                     close(3) 
                     write (FileName,'(I3)') k
                     s ='gnuplot -e "splot ''Lattice.dat'' matrix; set palette negative; set pm3d map; set term png size 1000,1000;'
                     PlotText = trim(s)//'set key off;set output'''//trim(FileName)//'.png''; replot;"'
                     call system(PlotText)
                     k = k + 1
                 end if
             End if
         End If
         
         if (Check .and. (mod(u,(l*l))==0)) then
                eM = eM + real(En)/N
                mM = mM + real(Mag)/N
                E2 = E2 + (real(En)/N)**2
                M2 = M2 + (real(Mag)/N)**2
                Counter = Counter + 1
         end if
     end do
     If ((((T > 1.9999).and.(T < 2.0001)) .or. ((T > 2.499).and.(T < 2.501))) .and. (Check .eqv. .False.)) then
         close (1)
         print 8, "Temperature is :", T
         8 format(1x,a,2x,f6.3)
         print *, "Do you want save Temperature data [Y/N]"
         read *, SaveData
         If (SaveData == 'y') then
              call system('gnuplot -p TransientResults.gnu')
         end if
     End If
     !$End OMP Parallel
   end subroutine Metropolis
   
! ====================================================================================================

   subroutine TempFunctions(T, En, Mag, M2Mean, MeanMag, E2Mean, MeanEn)
     real, intent(in)   :: T
     Integer, intent(in):: En, Mag
     real :: X, Cv ! --------------------------- Magnetic Susceptibility & Spacific Heat
     real :: E2Mean, M2Mean,MeanMag,MeanEn ! --- Energy and Magnetization/Spin

     If ((T < 2.269) .and. SoluSign) then
         MB = (1 - (1/ (sinh(2/T))**4))**(0.125)
     else if ((T < 2.269) .and. (SoluSign .eqv. .False.)) then
         MB = -(1 - (1/ (sinh(2/T))**4))**(0.125)
     else
         MB = 0
     end If
     MagSpin = real(Mag) / N;  EnSpin  = real(En)  / N;
     
     X  = (N/T)*((M2Mean/real(Counter))-(MeanMag/real(Counter))**2)
     Cv = ((real(l)/T)**2)*((E2Mean/real(Counter))-(MeanEn/real(Counter))**2)
     write(2,*) T, MeanMag/real(Counter), MeanEn/real(Counter), X, Cv, MB
     
   end subroutine TempFunctions
   
   
! ====================================================================================================
   
   function Neighbors(l , x , y) 
      integer :: l , x , y
      integer, dimension(4) :: Neighbors
      if (x == 1 ) then
         Neighbors(1) = l
         Neighbors(2) = x + 1 
      else if (x == l) then
         Neighbors(1) = x - 1
         Neighbors(2) = 1
      else 
         Neighbors(1) = x - 1
         Neighbors(2) = x + 1 
      end if
   
      if (y == 1 ) then
         Neighbors(3)= l
         Neighbors(4) = y + 1 
      else if (y == l) then
         Neighbors(3) = y - 1 
         Neighbors(4) = 1
      else
         Neighbors(3) = y - 1
         Neighbors(4) = y + 1 
      end if
     
   end function Neighbors


! ==================================================================================================== 
   subroutine DateAndTime(D, time)
     integer*4, intent(out) :: D(3), time(3)
     
     call idate(D)   
     call itime(time)
     
     print 1000, D(2), D(1), D(3), time
     1000 format (20x,'Date: ', i2.2, '/', i2.2, '/', i4.4,5x,'time: ',i2.2, ':', i2.2, ':', i2.2 )
     
   end subroutine DateAndTime

end program main
