#ifdef CRTM
subroutine da_crtm_sensor_descriptor(nsensor,sensor_descriptor)

 integer,        intent(in)  :: nsensor
 character(len=80), intent(inout) :: sensor_descriptor(nsensor)

 integer :: i,platform_id,satellite_id,sensor_id
 character (len=80) :: crtm_sat, crtm_sensor

   if (trace_use) call da_trace_entry("da_crtm_sensor_descriptor")  

  do i=1,nsensor
     
     platform_id  = rtminit_platform(i)
     satellite_id = rtminit_satid(i)
     sensor_id    = rtminit_sensor(i)

     if (trim(crtm_platform_name(platform_id)) == 'eos') then
        if (satellite_id == 2) crtm_sat='aqua'
        if (satellite_id == 1) crtm_sat='terra'
     else if (trim(crtm_platform_name(platform_id)) == 'metop') then
        if (satellite_id == 1) crtm_sat='metop-b'
        if (satellite_id == 2) crtm_sat='metop-a'
        if (satellite_id == 3) crtm_sat='metop-c'
     else
        write(crtm_sat, '(a,i2.2)')  &
             trim( crtm_platform_name(platform_id) ),satellite_id
     end if

     if ( trim(crtm_sensor_name(sensor_id)) == 'airs' ) then
#ifdef CRTM_1_1
        crtm_sensor='airs281SUBSET'  ! CRTM_1_1
#else
        crtm_sensor='airs281'        ! CRTM_1_2
#endif
     elseif ( trim(crtm_sensor_name(sensor_id)) == 'hirs' ) then
        if (satellite_id <= 14) crtm_sensor='hirs2'
        if (satellite_id >= 15 .and. satellite_id <= 17) crtm_sensor='hirs3'
        if (satellite_id == 18) crtm_sensor='hirs4'
     elseif ( trim(crtm_sensor_name(sensor_id)) == 'avhrr' ) then
        if (satellite_id <= 14) crtm_sensor='avhrr2'
        if (satellite_id >= 15 .and. satellite_id <= 17) crtm_sensor='avhrr3'
        if (satellite_id == 18) crtm_sensor='avhrr4'
     else
        crtm_sensor=crtm_sensor_name(sensor_id)
     end if

     sensor_descriptor(i)=trim(crtm_sensor)//'_'//trim(crtm_sat)

  end do

   if (trace_use) call da_trace_exit("da_crtm_sensor_descriptor")  

end subroutine da_crtm_sensor_descriptor
#endif
