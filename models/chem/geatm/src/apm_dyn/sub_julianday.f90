
subroutine julian(y,m,d,julianday)
implicit none
integer i
integer y, m, d, dm(12), dmtem(12,2)
integer id400,id100,id4
integer sum,julianday

data dmtem/31,28,31,30,31,30,31,31,30,31,30,31,&
         & 31,29,31,30,31,30,31,31,30,31,30,31/


! check whether year y is a leap year
! if the answer is yes, we should use dm(1:12,2)
id400=mod(y,400)
id100=mod(y,100)
id4=mod(y,4)
if((id400.eq.0).or.((id4.eq.0).and.(id100.ne.0))) then
 do i=1,12
  dm(i)=dmtem(i,2)
 enddo
else
do i=1,12
  dm(i)=dmtem(i,1)
 enddo
endif

sum=0
do i=1,m-1
 sum=sum+dm(i)
enddo
julianday=sum+d

end
