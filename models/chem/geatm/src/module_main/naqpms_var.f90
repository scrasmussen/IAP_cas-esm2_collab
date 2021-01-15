module naqpms_varlist

!====================

logical,parameter :: laerv2=.true.

logical,parameter :: laerv1=.false.

!integer,parameter :: nmoxdt=5,nm12=12

integer,parameter :: naerbin=2,naersp=14,nbinbdy=naerbin+1

integer,parameter :: idx_ppm=1,idx_bc=2,idx_oc=3,idx_so4=4,idx_nh4=5,idx_no3=6

integer,parameter :: idx_soa01=7,idx_soa02=8,idx_soa03=9,idx_soa04=10,idx_soa05=11,idx_soa06=12

real,parameter :: aer_rhop(1:naersp)=(/1.0E06,2.0E06,1.0E06,1.5E06,1.5E06 &
                                      ,1.5E06,1.0E06,1.0E06,1.0E06,1.0E06 &
                                      ,1.0E06,1.0E06,1.0E06,1.0E06  /)

real,parameter :: diambdy(1:nbinbdy)=(/0.1,2.5,10.0/) ! um

real :: diamlgm(naerbin)

integer :: ist_aerom,ied_aerom,ntr_aerom

integer :: iedgas


!====================

! flag vars

logical :: lglbrun

integer :: idifvert,ichemgas,idry,iglobal,imodis,iopt_gaschem

! chenxsh@20181206
real    :: degres
integer :: nxxppt,nxxppt2,nxxpptb

logical :: lrd_lai

character*3 :: ddep_flag


character*10 :: flagadv

character*10 :: agtflag

integer,parameter :: itzon=0

! shun@20161227
logical :: lgaschemsmp
integer :: nhfq_updtmet(5)
integer :: nhfq_output(5)
real    :: dtstep_syn(5)
integer :: ndt_syn(5)
integer :: idt_syn(5)
character :: caveoutclab*20
!integer :: igassmp

!=====================
! ip var for tracers
integer,allocatable,dimension(:) :: ip2mem
integer,allocatable,dimension(:,:) :: ip3mem
integer,allocatable,dimension(:,:,:) :: ip4mem
integer,allocatable,dimension(:,:,:,:) :: ip5mem
integer,allocatable,dimension(:,:,:,:) :: ip5memc,ip5memcs

integer,allocatable,dimension(:,:) :: ip2memGas

integer,allocatable,dimension(:,:,:,:) :: ipSMmem

integer,allocatable,dimension(:,:,:) :: ip3memaer

integer,allocatable,dimension(:,:,:) :: ip_emit2d


integer,allocatable,dimension(:,:,:,:) :: ip4mem_aer

integer,allocatable,dimension(:,:,:) :: ip2memaer

integer,allocatable,dimension(:,:) :: ip4aveout


!integer,allocatable,dimension(:,:,:,:) :: ip4mem_oxdt

!integer, allocatable, dimension(:,:,:) :: ib4mem,jb4mem
!integer,allocatable,dimension(:,:,:,:) :: ib5mem,jb5mem

!=================================================================
! naqpms tracers

real,allocatable,dimension(:) :: avgvar

! global conditions
real,allocatable,dimension(:) :: globalno2,globalo3,globalco
real,allocatable,dimension(:) :: topo3



real,allocatable,dimension(:) :: gas
real,allocatable,dimension(:) :: gasOLD
real,allocatable,dimension(:) :: aer,aer_src

real,allocatable,dimension(:) :: aerom

!real,allocatable,dimension(:) :: mmean_oxdt

! hygroscopic growth factor
real,allocatable,dimension(:) :: hgfac
real,allocatable,dimension(:) :: awc3d

! FOR HETEROGENEOUS CHEMISTRY
real,allocatable,dimension(:) :: DUSTCOMP ! DUST COMPOSITIONS 
! 1: CACO3, 2 MAGCO3; 3. K2CO3, 4:NA2CO3, 5: SIO2,
! 6: AL2O3,  7: SO4, 8: NO3, 9: Fe(II),10: Total
! Fe(III) 11: Coated Fe(III)
real,allocatable,dimension(:) :: SEACOMP
! 1: CL, 2: NA 3 SS-SO4, 4: MG, 5: Ca, 6: K 7: SO4 8: NO3


!=========================
! emit
real,allocatable,dimension(:) :: EmtaGas,EmtpGas,EmttGas,EmtbGas,Emt5Gas 

real,allocatable,dimension(:) :: emit2d

real,allocatable,dimension(:,:) :: emt2d_zfrc

character*2,allocatable,dimension(:) :: cfmode

!==========
! dust
real,allocatable,dimension(:) :: DUSTEMISS,DUSTDRY,DUSTWET,DUSTGRAV
!  DUST EMISS,DRY and WET and GRAVITY  DEPOSITION : kg/m2/hour

real,allocatable,dimension(:) :: EMITFACT,TOTALDUST

! FOR HETEROGENEOUS CHEMISTRY
! kg/m2/hour
real,allocatable,dimension(:) :: DUSTDRYSO4,DUSTDRYNO3,DUSTDRYFeII,DUSTDRYFeIII 
real,allocatable,dimension(:) :: DUSTWETSO4,DUSTWETNO3,DUSTWETFeII,DUSTWETFeIII 
real,allocatable,dimension(:) :: DUSTGRAVSO4,DUSTGRAVNO3,DUSTGRAVFeII,DUSTGRAVFeIII 
real, allocatable, dimension(:) :: DryVeldust



!===========
! sea salt
real,allocatable,dimension(:) :: SEAEMISS

!===========
real,allocatable,dimension(:) :: SourceMark
real,allocatable,dimension(:) :: MapSource
real,allocatable,dimension(:) :: tmpMarkCon

!===========
! drydep
real, allocatable, dimension(:) :: DryVelGas
!real, allocatable, dimension(:) :: dvel_z03

real, allocatable, dimension(:) :: dryvelaer

real, allocatable, dimension(:) :: aerom_wdepfld,aerom_wdepfld2

!===========
! wetdep
real,allocatable,dimension(:) :: WETDEP,WETDEP2
!2-D array of wet deposited mass (mol/ha,g/ha)

real,allocatable,dimension(:) :: gscav_so2,ascav_so4


!===================================================================
! model diags

real,allocatable,dimension(:) :: CPH

real,allocatable,dimension(:) :: jo1d,jno2
real,allocatable,dimension(:) :: EXT,EXTASO4,EXTANO3,EXTANH4,EXTBC,EXTOC
real,allocatable,dimension(:) :: VISIB,UVB,UVBS,UVA
real,allocatable,dimension(:) :: VIS,SSA,AOD,CLDOPD
real,allocatable,dimension(:) :: DUSO2,DUO3,DUNO2
real,allocatable,dimension(:) :: ANA,ASO4,ANH4,ANO3,ACL,OPE

real,allocatable,dimension(:) :: DUSTEXT 
real,allocatable,dimension(:) :: DUSTAOD,PBLAOD


!============================================
! process analysis

!integer,dimension(102) :: PrintTermGas  ! process by lijie

integer, allocatable, dimension(:) :: IGGPOS,IGOpos

real,allocatable,dimension(:) :: GasTermBal

integer, allocatable, dimension(:,:,:,:) :: ipGasTermBal




contains

subroutine allo_naqpms_var(nx,ny,nzz,nest,sx,ex,sy,ey,mem_per_block &
                          ,igas,IAER,ISIZE,ndustcom,nseacom &
                          ,ifsmt,ifsm,idmSet,ismMax,igMark &
                          ,nctg_emt &
                          ,lprocess,iprocess,PrintTermGas &
                          ,mem2d,mem3d,mem4d,mem5d,mem2dgas,mem3daer,mem_emt2d )
 implicit none

 integer :: ne,is,k,ig,ia,iduc,idm,ism,ictg
 integer :: ii,mem5d,mem4d,mem3d,mem2d,mem2dgas,mem3daer

 integer :: mem4d_aerom

 integer :: memtmp

!shun@20161228
! integer :: mem4d_oxdt

 integer :: mem5dc,mem5dterm

 integer :: mem5dsmark

 integer :: nest
 integer :: nx(5),ny(5),nzz
 integer :: sx(5), ex(5), sy(5), ey(5)

 integer :: mem_per_block(5)

 integer :: igas

 integer :: IAER,ISIZE

 integer :: ndustcom,nseacom

 integer :: nctg_emt

 integer :: ifsmt
 integer :: ifsm(5)
 integer :: idmSet,ismMax
 integer :: igMark(idmSet)

 logical :: lprocess

 integer :: PrintTermGas(102)
 integer :: ip
 integer :: iPrintTermGas
 integer :: iprocess


 integer :: mem_emt2d

 integer :: mem2daer

 integer :: ibin

 integer :: imon

!===================================================================

 do ibin=1,naerbin
   diamlgm(ibin)=sqrt(diambdy(ibin)*diambdy(ibin+1)) ! particle diameter : um
 enddo

 ! number of memory for 2d variables

  allocate(ip2mem(nest))
  ii=1
  do ne=1,nest
   ip2mem(ne)=ii
   ii=ii+mem_per_block(ne)
  enddo


 mem2d=0          
 do ne=1,nest
   mem2d=mem2d+mem_per_block(ne)
 enddo

 allocate( topo3(mem2d) &
          ,CLDOPD(mem2d) &
          ,AOD(mem2d) &
          ,DUSTAOD(mem2d) &
          ,PBLAOD(mem2d) &
          ,DUSO2(mem2d) &
          ,DUO3(mem2d) &
          ,DUNO2(mem2d) )

 topo3=0.0



 allocate( SEAEMISS(mem2d))


 allocate( EMITFACT(mem2d) &
          ,TOTALDUST(mem2d) &
          ,DUSTEMISS(mem2d) )

 allocate( DUSTDRY(mem2d) &
          ,DUSTWET(mem2d) &
          ,DUSTGRAV(mem2d))

 allocate( DUSTDRYSO4(mem2d) &
          ,DUSTDRYNO3(mem2d) &
          ,DUSTDRYFeII(mem2d) &
          ,DUSTDRYFeIII(mem2d))

 allocate( DUSTWETSO4(mem2d) &
          ,DUSTWETNO3(mem2d) &
          ,DUSTWETFeII(mem2d) &
          ,DUSTWETFeIII(mem2d) )

 allocate( DUSTGRAVSO4(mem2d) &
          ,DUSTGRAVNO3(mem2d) &
          ,DUSTGRAVFeII(mem2d) &
          ,DUSTGRAVFeIII(mem2d) )

  ! For Source
  if(ifsmt>0)then ! 
     allocate(MapSource(mem2d))
     allocate(tmpMarkCon(mem2d))
  endif


  allocate(ip2memGas(igas,nest))
  ii=1
  do ne=1,nest
  do ig=1,igas
   ip2memGas(ig,ne)=ii
   ii=ii+mem_per_block(ne)
  enddo
  enddo


  mem2dgas=0   !number of memory for 2d variables such emit,dep
  do ne=1,nest
  do ig=1,igas
    mem2dgas=mem2dgas+mem_per_block(ne)
  enddo
  enddo

  allocate( EmtaGas(mem2dgas) &
           ,EmtpGas(mem2dgas) &
           ,EmttGas(mem2dgas) &
           ,EmtbGas(mem2dgas) )

  allocate( Emt5Gas(mem2dgas) )

  allocate( DryVelGas(mem2dgas) )


  allocate(ip2memaer(naerbin,naersp,nest))
  ii=1
  do ne=1,nest
  do ia=1,naersp
  do is=1,naerbin
   ip2memaer(is,ia,ne)=ii
   ii=ii+mem_per_block(ne)
  enddo
  enddo
  enddo

  mem2daer=0
  do ne=1,nest
  do ia=1,naersp
  do is=1,naerbin
    mem2daer=mem2daer+mem_per_block(ne)
  enddo
  enddo
  enddo
  allocate(dryvelaer(mem2daer))
  allocate(aerom_wdepfld(mem2daer));  aerom_wdepfld=0.0
  allocate(aerom_wdepfld2(mem2daer)); aerom_wdepfld2=0.0

  


!  allocate( dvel_z03(mem2dgas) )

  allocate( WETDEP(mem2dgas) &
           ,WETDEP2(mem2dgas) )

  WETDEP=0.0
  WETDEP2=0.0


  allocate(ip_emit2d(igas,nctg_emt,nest))
  ii=1
  do ne=1,nest
  do ictg=1,nctg_emt
  do ig=1,igas
   ip_emit2d(ig,ictg,ne)=ii
   ii=ii+mem_per_block(ne)
  enddo
  enddo
  enddo

  mem_emt2d=0
  do ne=1,nest
  do ictg=1,nctg_emt
  do ig=1,igas
    mem_emt2d=mem_emt2d+mem_per_block(ne)
  enddo
  enddo
  enddo
  allocate(emit2d(mem_emt2d)); emit2d=0.0



!----- FOR DUST AND SEA SALT -----
 mem3daer =0
 DO NE = 1, NEST
 DO IA = 1, IAER
 DO IS = 1, ISIZE
   mem3daer = mem3daer + mem_per_block(ne) 
 ENDDO
 ENDDO
 ENDDO

 allocate(ip3memaer(ISIZE, IAER, NEST))
 II = 1
 DO NE = 1, NEST
 DO IA = 1, IAER
 DO IS  =1,ISIZE
  ip3memaer(IS,IA,NE) = II
  II = II + mem_per_block(ne) 
 ENDDO
 ENDDO
 ENDDO

 allocate(DryVeldust(mem3daer))


 allocate(ip3mem(nzz,nest))
 ii=1
 do ne=1,nest
 do k=1,nzz
   ip3mem(k,ne)=ii
   ii=ii+mem_per_block(ne)
 enddo
 enddo

 mem3d=0
 do ne=1,nest
 do k=1,nzz
    mem3d=mem3d+mem_per_block(ne)
 enddo
 enddo

 allocate( hgfac(mem3d) )
 if(laerv2) then
   allocate( awc3d(mem3d) )
 endif

 allocate( globalno2(mem3d) &
          ,globalo3(mem3d) &
          ,globalco(mem3d) )

 allocate( EXT(mem3d) &
          ,UVA(mem3d) &
          ,UVB(mem3d) &
          ,UVBS(mem3d) &
          ,VIS(mem3d) &
          ,VISIB(mem3d) )

 allocate( DUSTEXT(mem3d) )

 allocate( EXTASO4(mem3d) &
          ,EXTANO3(mem3d) &
          ,EXTANH4(mem3d) &
          ,EXTBC(mem3d) &
          ,EXTOC(mem3d) )

 allocate( jo1d(mem3d) &
          ,jno2(mem3d) &
          ,SSA(mem3d) )

 allocate( ANA(mem3d) &
          ,ASO4(mem3d) &
          ,ANO3(mem3d) &
          ,ANH4(mem3d) &
          ,ACL(mem3d) )

 allocate( CPH(mem3d) )  ! AQUEOUS CHEMISTRY 
 allocate( OPE(mem3d) )  ! OPE

 allocate( gscav_so2(mem3d) )
 allocate( ascav_so4(mem3d) )



  allocate(ip4mem(nzz,igas,nest))
  ii=1
  do ne=1,nest
  do ig=1,igas
  do k=1,nzz
   ip4mem(k,ig,ne)=ii
   ii=ii+mem_per_block(ne)
  enddo
  enddo
  enddo


  mem4d=0        ! number of memory for 4d variables
  do ne=1,nest
  do ig=1,igas
  do k=1,nzz
    mem4d=mem4d+mem_per_block(ne)
  enddo
  enddo
  enddo

  allocate(gas(mem4d)); gas=0.0

  if(lprocess)  then
    allocate(gasOLD(mem4d))
    gasOLD=0.0
  endif


  mem4d_aerom=0
  do ne=1,nest
  do ig=1,naersp
  do is=1,naerbin
  do k=1,nzz
    mem4d_aerom=mem4d_aerom+mem_per_block(ne)
  enddo
  enddo
  enddo
  enddo
  allocate(aerom(mem4d_aerom))
  aerom=0.0


  allocate(ip4mem_aer(nzz,naerbin,naersp,nest))
  ii=1
  do ne=1,nest
  do ia=1,naersp
  do is=1,naerbin
  do k=1,nzz
    ip4mem_aer(k,is,ia,ne)=ii
    ii=ii+mem_per_block(ne)
  enddo
  enddo
  enddo
  enddo


  allocate(ip5mem(nzz,isize,iaer,nest))
  ii=1
  do ne=1,nest
  do ia=1,iaer
  do is=1,isize
  do k=1,nzz
    ip5mem(k,is,ia,ne)=ii
    ii=ii+mem_per_block(ne)
  enddo
  enddo
  enddo
  enddo

  mem5d=0        ! number of memory for 5D variables
  do ne=1,nest
  do ia=1,iaer
  do is=1,isize
  do k=1,nzz
    mem5d=mem5d+mem_per_block(ne)
  enddo
  enddo
  enddo
  enddo

  allocate(aer(mem5d),aer_src(mem5d))
  aer=0.0
  aer_src=0.0


  allocate(ip5memc(nzz,isize,ndustcom,nest))
  ii = 1
  do ne = 1,nest
  do iduc = 1, ndustcom
  do is = 1, isize
  do k = 1, nzz
    ip5memc (k, is, iduc, ne) = ii
    ii = ii + mem_per_block(ne)
  enddo
  enddo
  enddo
  enddo

  mem5dc=0
  do ne = 1, nest
  do iduc = 1, ndustcom
  do is = 1, isize
  do k = 1, nzz
    mem5dc = mem5dc + mem_per_block(ne)
  enddo
  enddo
  enddo
  enddo

  allocate (DUSTCOMP(mem5dc))
  DUSTCOMP=0.0


  allocate(ip5memcs(nzz,isize,nseacom,nest))
  ii = 1
  do ne = 1, nest
  do iduc = 1, nseacom
  do is = 1, isize
  do k = 1, nzz
    ip5memcs(k,is,iduc, ne) = ii
    ii = ii + mem_per_block(ne)
  enddo
  enddo
  enddo
  enddo

  mem5dc = 0
  do ne = 1, nest
  do iduc = 1, nseacom
  do is = 1, isize
  do k = 1, nzz
     mem5dc = mem5dc + mem_per_block(ne)
  enddo
  enddo
  enddo
  enddo

  allocate ( SEACOMP(mem5dc) )
  SEACOMP=0.0




  allocate(ip4aveout(2000,nest))
  ii = 1
  do ne = 1, nest
  do k = 1, 2000
    ip4aveout(k,ne) = ii
    ii = ii + mem_per_block(ne)
  enddo
  enddo

  memtmp = 0
  do ne = 1, nest
  do k = 1, 2000 ! need to set
     memtmp = memtmp + mem_per_block(ne)
  enddo
  enddo
  allocate(avgvar(memtmp))




 !=============================================
 ! For Source Mark
 if(ifsmt>0)then  ! For Source Mark
   mem5dsmark=0        ! number of memory for 5d variables
   do ne=1,nest
     if(ifsm(ne)==1)then
      do idm=1,idmSet
      do ism=1,ismMax
      do k=1,nzz
         mem5dsmark=mem5dsmark+mem_per_block(ne)
      enddo
      enddo
      enddo
    endif
   enddo
   allocate(SourceMark(mem5dsmark))
   SourceMark=0.0

   allocate(ipSMmem(nzz,ismMax,idmSet,nest))
   ii=1
   do ne=1,nest
    if(ifsm(ne)==1)then
      do idm=1,idmSet
      do ism=1,ismMax
      do k=1,nzz
         ipSMmem(k,ism,idm,ne)=ii
         ii=ii+mem_per_block(ne)
      enddo
      enddo
      enddo
    endif
   enddo
 endif
 !====================================

if(lprocess) then

 iPrintTermGas = 0
 do ig=1,igas
    iPrintTermGas = iPrintTermGas + PrintTermGas(ig)
 enddo

 allocate(IGGPOS( iPrintTermGas ))
 allocate(IGOPos( igas ))

 iPrintTermGas = 0
 IGOPos=0
 do ig=1,igas
  if( PrintTermGas(ig) == 1 )then
    iPrintTermGas = iPrintTermGas + 1
    IGGPOS(iPrintTermGas) = ig
    IGOPos(ig) = iPrintTermGas
  endif
 enddo

 mem5dterm=0        ! number of memory for Track -GAS -Balance
 do ne=1,nest
 do ig=1,iPrintTermGas
 do ip=1,iprocess ! this term is very large
 do k=1,nzz
   mem5dterm=mem5dterm+mem_per_block(ne)
 enddo
 enddo
 enddo
 enddo
 allocate(GasTermBal(mem5dterm))
 GasTermBal=0.0


 allocate(ipGasTermBal(nzz,iprocess,iPrintTermGas,nest))
 ii=1
 do ne=1,nest
 do ig=1,iPrintTermGas
 do ip=1,iprocess
 do k=1,nzz
    ipGasTermBal(k,ip,ig,ne)=ii
    ii=ii+mem_per_block(ne)
 enddo
 enddo
 enddo
 enddo

endif

end subroutine allo_naqpms_var



end module naqpms_varlist
