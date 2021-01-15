

MODULE module_gtiming

   INTEGER, PARAMETER, PRIVATE :: cnmax = 30
   INTEGER, PRIVATE, DIMENSION(cnmax) :: count_int1 , count_rate_int1 , count_max_int1
   INTEGER, PRIVATE, DIMENSION(cnmax) :: count_int2 , count_rate_int2 , count_max_int2
   INTEGER, PRIVATE :: cn = 0
   REAL, PRIVATE :: elapsed_seconds , elapsed_seconds_total = 0
   REAL, PRIVATE :: cpu_1 , cpu_2 , cpu_seconds , cpu_seconds_total = 0

CONTAINS

   SUBROUTINE init_module_gtiming
      cn = 0
   END SUBROUTINE init_module_gtiming


   SUBROUTINE start_gtiming ( string ) 

      IMPLICIT NONE
      CHARACTER *(*) :: string

      cn = cn + 1
      IF ( cn .gt. cnmax ) THEN
        print*,'module_timing: clock nesting error (too many nests)'
        RETURN
      ENDIF
      CALL SYSTEM_CLOCK ( count_int1(cn) , count_rate_int1(cn) , count_max_int1(cn) )


   END SUBROUTINE start_gtiming


   SUBROUTINE end_gtiming ( string, idomain, myid )

      IMPLICIT NONE
      integer :: idomain, myid
      CHARACTER *(*) :: string

      IF ( cn .lt. 1 ) THEN
        if(myid.eq.0) print*,'module_timing: clock nesting error, cn<1'
      ELSE IF ( cn .gt. cnmax ) THEN
        if(myid.eq.0) print*,'module_timing: clock nesting error, cn>cnmax'
      ENDIF

      CALL SYSTEM_CLOCK ( count_int2(cn) , count_rate_int2(cn) , count_max_int2(cn) )


      IF ( count_int2(cn) < count_int1(cn) ) THEN
         count_int2(cn) = count_int2(cn) + count_max_int2(cn)
      ENDIF

      count_int2(cn) = count_int2(cn) - count_int1(cn)
      elapsed_seconds = REAL(count_int2(cn)) / REAL(count_rate_int2(cn))
      elapsed_seconds_total = elapsed_seconds_total + elapsed_seconds

      if(myid.eq.0) then
          WRITE(6,'(A,i2.2,2x,F10.5,A,A)') &
            ' d',idomain,elapsed_seconds,' elapsed seconds % ',TRIM(string)

      endif

      cn = cn - 1

   END SUBROUTINE end_gtiming

END MODULE module_gtiming
