 subroutine sample( step, time, aver, switch )

   use global
   implicit none

   type(averages), intent(in) :: aver


   integer, intent(in) :: step, switch

   real(8) :: totEnergy, potEnergy, kinEnergy, temp
   real(8), intent(in) :: time


   open(unit=7,file='output.dat',status='old',position='append')

   if( switch == 0 ) then

     write(unit=7,fmt="(//'       step',2x,'       time',2x,'  totEnergy',2x,'  potEnergy',2x,&
                      '  kinEnergy',2x,'temperature')") 
     write(unit=7,fmt="('===================================================================&
                         ===================================='),/") 

     close(unit=7)

     return
  
   end if

   totEnergy = aver%totEnergy(1)/aver%counter
   potEnergy = aver%potEnergy(1)/aver%counter
   kinEnergy = aver%kinEnergy(1)/aver%counter
   temp = aver%temp(1)/aver%counter

   write(unit=7,fmt="(i11,2x,f11.5,4(2x,e11.4))") step, time, totEnergy, potEnergy,&
                     kinEnergy, temp

  close(unit=7)

 end subroutine sample
