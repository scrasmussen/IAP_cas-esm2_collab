subroutine da_analysis_stats  (grid, stats_unit)
   
   !------------------------------------------------------------------------
   ! Purpose: Calculate min, max, mean and RMS of input 1d field.
   !------------------------------------------------------------------------

   implicit none

   type (domain), intent (in) :: grid
   integer,       intent (in) :: stats_unit ! Output unit for stats.
   
   integer :: i, j, k
   integer :: ij_g, ijk_g   ! ij, ijk for global domain. 
   integer :: kdim  ! k range 

   real    :: um, vm, tm, pm, qm , qcwm, qrnm ! On local domain.
   real    :: rij_g, rijk_g      ! On global domain.

   type (maxmin_field_type) :: max_u(kts:kte), max_v(kts:kte), &
                               max_t(kts:kte), max_p(kts:kte), &
                               max_q(kts:kte), &
                               max_qcw(kts:kte),max_qrn(kts:kte),&
                               min_u(kts:kte), min_v(kts:kte), &
                               min_t(kts:kte), min_p(kts:kte), &
                               min_q(kts:kte),&
                               min_qcw(kts:kte),min_qrn(kts:kte)
 
   real                        uv(kts:kte), vv(kts:kte), &
                               tv(kts:kte), pv(kts:kte), &
                               qv(kts:kte),&
                               qcwv(kts:kte),qrnv(kts:kte)

   call da_trace_entry("da_analysis_stats")

   kdim = kte-kts+1

   ij_g = mix * mjy
   ijk_g = ij_g * mkz
   
   rij_g  = 1.0/real(ij_g)
   rijk_g = 1.0/real(ijk_g)
  
   if (rootproc) then
      write(unit=stats_unit, fmt='(/a/)') &
         ' Minimum of gridded analysis increments'
   if (crtm_cloud .or. use_radar_rf) then
      write(unit=stats_unit, fmt='(8a/)') &
         ' Lvl         ', &
         'u     i    j          ', &
         'v     i    j          ', &
         't     i    j          ', &
         'p     i    j          ', &
         'q     i    j          ',&
         'qcw     i    j          ',&
         'qrn     i    j          '
   else
      write(unit=stats_unit, fmt='(6a/)') &
         ' Lvl         ', &
         'u     i    j          ', &
         'v     i    j          ', &
         't     i    j          ', &
         'p     i    j          ', &
         'q     i    j          '
   end if
   end if

   call da_maxmin_in_field(grid%xa%u(its:ite,jts:jte,kts:kte), max_u, min_u)
   call da_proc_maxmin_combine(kdim, max_u, min_u)
   call da_maxmin_in_field(grid%xa%v(its:ite,jts:jte,kts:kte), max_v, min_v)
   call da_proc_maxmin_combine(kdim, max_v, min_v)
   call da_maxmin_in_field(grid%xa%t(its:ite,jts:jte,kts:kte), max_t, min_t)
   call da_proc_maxmin_combine(kdim, max_t, min_t)
   call da_maxmin_in_field(grid%xa%p(its:ite,jts:jte,kts:kte), max_p, min_p)
   call da_proc_maxmin_combine(kdim, max_p, min_p)
   call da_maxmin_in_field(grid%xa%q(its:ite,jts:jte,kts:kte), max_q, min_q)
   call da_proc_maxmin_combine(kdim, max_q, min_q)

   if (crtm_cloud .or. use_radar_rf) then
      call da_maxmin_in_field(grid%xa%qcw(its:ite,jts:jte,kts:kte), max_qcw, min_qcw)
      call da_proc_maxmin_combine(kdim, max_qcw, min_qcw)
      call da_maxmin_in_field(grid%xa%qrn(its:ite,jts:jte,kts:kte), max_qrn, min_qrn)
      call da_proc_maxmin_combine(kdim, max_qrn, min_qrn)
   end if

        

   um = 999999.0
   vm = 999999.0
   tm = 999999.0
   pm = 999999.0
   qm = 999999.0

   if (crtm_cloud .or. use_radar_rf) then
      qcwm = 999999.0
      qrnm = 999999.0
   end if

   do k = kts, kte   
      if (rootproc) then
         if (crtm_cloud .or. use_radar_rf) then
            write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),3(e12.4,2i5))') k, &
              min_u(k), min_v(k), min_t(k), min_p(k), min_q(k),min_qcw(k),min_qrn(k)
         else
             if ( abs(min_q(k)%value) < 1.e-30 ) min_q(k)%value = 0.0
             write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),e12.4,2i5)') k, &
              min_u(k), min_v(k), min_t(k), min_p(k), min_q(k)
         end if
      end if

      um=minval(min_u(:)%value)
      vm=minval(min_v(:)%value)
      tm=minval(min_t(:)%value)
      pm=minval(min_p(:)%value)
      qm=minval(min_q(:)%value)

      if (crtm_cloud .or. use_radar_rf) then
         qcwm=minval(min_qcw(:)%value)
         qrnm=minval(min_qrn(:)%value)
      end if

   end do

   if (rootproc) then

     if (crtm_cloud .or. use_radar_rf) then
        write(unit=stats_unit, fmt='(a,4(f12.4,10x),3(e12.4))') ' ALL', &
          um, vm, tm, pm, qm, qcwm, qrnm
     else
        write(unit=stats_unit, fmt='(a,4(f12.4,10x),e12.4)') ' ALL', &
          um, vm, tm, pm, qm
     end if
   
  
      write(unit=stats_unit, fmt='(/a/)') &
         ' Maximum of gridded analysis increments'

     if (crtm_cloud .or. use_radar_rf) then
        write(unit=stats_unit, fmt='(8a/)') &
           ' Lvl         ', &
           'u     i    j          ', &
           'v     i    j          ', &
           't     i    j          ', &
           'p     i    j          ', &
           'q     i    j          ',&
           'qcw   i    j        ',&
           'qrn   i    j        '
     else
        write(unit=stats_unit, fmt='(6a/)') &
           ' Lvl         ', &
           'u     i    j          ', &
           'v     i    j          ', &
           't     i    j          ', &
           'p     i    j          ', &
           'q     i    j          '
     end if
   end if

   do k = kts, kte
      if (rootproc) then
        if (crtm_cloud .or. use_radar_rf) then
           write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),3(e12.4,2i5))') k, &
             max_u(k), max_v(k), max_t(k), max_p(k), max_q(k),max_qcw(k),max_qrn(k) 
        else
           if ( abs(max_q(k)%value) < 1.e-30 ) max_q(k)%value = 0.0
           write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),e12.4,2i5)') k, &
             max_u(k), max_v(k), max_t(k), max_p(k), max_q(k)
        end if
      end if

      um=maxval(max_u(:)%value)
      vm=maxval(max_v(:)%value)
      tm=maxval(max_t(:)%value)
      pm=maxval(max_p(:)%value)
      qm=maxval(max_q(:)%value)
      if (crtm_cloud .or. use_radar_rf) then
         qcwm=maxval(max_qcw(:)%value)
         qrnm=maxval(max_qrn(:)%value)
      end if

   end do

   if (rootproc) then
     if (crtm_cloud .or. use_radar_rf) then
        write(unit=stats_unit, fmt='(a,4(f12.4,10x),e12.4)') ' ALL', &
           um, vm, tm, pm, qm, qcwm, qrnm
     else
        write(unit=stats_unit, fmt='(a,4(f12.4,10x),e12.4)') ' ALL', &
           um, vm, tm, pm, qm
     end if
   
      write(unit=stats_unit, fmt='(/a/)') ' Mean of gridded analysis increments'
      

     if (crtm_cloud .or. use_radar_rf) then
        write(unit=stats_unit, fmt='(a/)') &
          ' Lvl        u           v           t           p           q           qcw           qrn'
     else
        write(unit=stats_unit, fmt='(a/)') &
          ' Lvl        u           v           t           p           q'
     end if
   end if

   um = 0.0
   vm = 0.0
   tm = 0.0
   pm = 0.0
   qm = 0.0


   do k = kts, kte
      uv(k) = sum(grid%xa%u(its:ite,jts:jte,k))
      vv(k) = sum(grid%xa%v(its:ite,jts:jte,k))
      tv(k) = sum(grid%xa%t(its:ite,jts:jte,k))
      pv(k) = sum(grid%xa%p(its:ite,jts:jte,k))
      qv(k) = sum(grid%xa%q(its:ite,jts:jte,k))
   end do
   if (crtm_cloud .or. use_radar_rf) then
      qcwm = 0.0
      qrnm = 0.0
       do k = kts, kte
          qcwv(k) = sum(grid%xa%qcw(its:ite,jts:jte,k))
          qrnv(k) = sum(grid%xa%qrn(its:ite,jts:jte,k))
       end do
      call da_proc_sum_real  (qcwv)
      call da_proc_sum_real  (qrnv)
   end if

   call da_proc_sum_real  (uv)
   call da_proc_sum_real  (vv)
   call da_proc_sum_real  (tv)
   call da_proc_sum_real  (pv)
   call da_proc_sum_real  (qv)

   if (rootproc) then
      do k = kts, kte
         write(unit=stats_unit, fmt='(i4,4f12.4,4e12.4)') k, &
            uv(k)*rij_g, vv(k)*rij_g, tv(k)*rij_g, &
            pv(k)*rij_g, qv(k)*rij_g

         um=um+uv(k)
         vm=vm+vv(k)
         tm=tm+tv(k)
         pm=pm+pv(k)
         qm=qm+qv(k)
      end do
   end if

   if (rootproc) then
      write(unit=stats_unit, fmt='(a,4f12.4,4e12.4)') ' ALL', &
         um*rijk_g, vm*rijk_g, tm*rijk_g, pm*rijk_g, qm*rijk_g

      write(unit=stats_unit, fmt='(/a/)') ' RMSE of gridded analysis increments'

      write(unit=stats_unit, fmt='(a/)') &
         ' Lvl        u           v           t           p           q'
   end if

   um = 0.0
   vm = 0.0
   tm = 0.0
   pm = 0.0
   qm = 0.0
   uv = 0.0
   vv = 0.0
   tv = 0.0
   pv = 0.0
   qv = 0.0

   do k = kts, kte
      do j=jts,jte
         do i=its,ite
            uv(k) = uv(k) + grid%xa%u(i,j,k) * grid%xa%u(i,j,k)
            vv(k) = vv(k) + grid%xa%v(i,j,k) * grid%xa%v(i,j,k)
            tv(k) = tv(k) + grid%xa%t(i,j,k) * grid%xa%t(i,j,k)
            pv(k) = pv(k) + grid%xa%p(i,j,k) * grid%xa%p(i,j,k)
            qv(k) = qv(k) + grid%xa%q(i,j,k) * grid%xa%q(i,j,k)
         end do
      end do
   end do

   call da_proc_sum_real  (uv)
   call da_proc_sum_real  (vv)
   call da_proc_sum_real  (tv)
   call da_proc_sum_real  (pv)
   call da_proc_sum_real  (qv)

   if (rootproc) then
      do k = kts, kte
         write(unit=stats_unit, fmt='(i4,4f12.4,4e12.4)') k, &
            sqrt(uv(k)*rij_g), &
            sqrt(vv(k)*rij_g), &
            sqrt(tv(k)*rij_g), &
            sqrt(pv(k)*rij_g), &
            sqrt(qv(k)*rij_g)

         um=um+uv(k)
         vm=vm+vv(k)
         tm=tm+tv(k)
         pm=pm+pv(k)
         qm=qm+qv(k)
      end do
   end if

   if (rootproc) then
      write(unit=stats_unit, fmt='(a,4f12.4,4e12.4)') ' ALL', &
         sqrt(um*rijk_g), sqrt(vm*rijk_g), sqrt(tm*rijk_g), &
         sqrt(pm*rijk_g), sqrt(qm*rijk_g)
   end if

   call da_trace_exit("da_analysis_stats")

contains

#include "da_maxmin_in_field.inc"
   
end subroutine da_analysis_stats


