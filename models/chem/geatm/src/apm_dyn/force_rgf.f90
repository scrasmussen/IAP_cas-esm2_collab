
subroutine force_rgfac(label,ii,jj,kk,rgf1,rgf2,def,rgf,iostate)
implicit none
character :: label*200
integer   :: ii,jj,kk
real*8    :: rgf
real      :: rgf1,rgf2,def
integer   :: iostate

iostate=0

if( .not.(rgf.ge.rgf1.and.rgf.le.rgf2) ) then
!  write(6,*) trim(label)
!  write(6,*) ii,jj,kk,rgf
  iostate=1
  rgf=def
endif

end subroutine force_rgfac

