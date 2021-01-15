


! in k loop

             ! For gas speceis
             DO IG=1,iedgas
               ips1 = ips1 + 1
               I04 = IP4MEM(K,IG,NE)
               CALL PUTVALUE(MYID,gas(i04),SX(NE),EX(NE),SY(NE),EY(NE),&
               IPOLARMRK(3,IPS),IPP,ATESTR(K,ips1))
             ENDDO


              ! for aerosols
              DO IA=1,IAER
              DO IS=1,ISIZE
              !IPS1 = igas + (IA-1)*ISIZE + IS
              ips1 = ips1 + 1
              I05= IP5MEM(K,IS,IA,NE)
              CALL PUTVALUE(MYID,AER(I05),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(3,IPS),IPP,ATESTR(K,IPS1))
              ENDDO
              ENDDO

           !!!!!!!!!!!!!!!!!
          ! for Source Mark
          if(ifsmt>0 )then   !checking
            do idm=1,idmSet
            do ism=1,ismMax
             !ips1 = igas + iaer*isize +2+(idm-1)*ismMax + ism
             ips1 = ips1 + 1
             i0=ipSMmem(k,ism,idm,ne)
            call putvalue(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
            enddo
            enddo
           endif
           !!!!!!!!!!!!!!!!!!!!
           ! for sea salt and dust compositions
!          if(ifseacom.eq.1) then
          do iduc = 1, nseacom
          do is = 1, isize
            ips1 = ips1 + 1
            i05c  = ip5memcs(k,is,iduc,ne)
     call putvalue(myid,seacomp(i05c),sx(ne),ex(ne),sy(ne),ey(ne), &
              IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
          enddo
          enddo
!          endif

!          if(ifdustcom.eq.1) then
          do iduc = 1, ndustcom
          do is = 1, isize
            ips1 = ips1 + 1
            i05c  = ip5memc(k,is,iduc,ne)
     call putvalue(myid,dustcomp(i05c),sx(ne),ex(ne),sy(ne),ey(ne), &
              IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
          enddo
          enddo
!          endif

        if(laerv2) then
           do iaersp  = 1, naersp
           do iaerbin = 1, naerbin
             ips1 = ips1 + 1
             i04aer  = ip4mem_aer(k,iaerbin,iaersp,ne) ! juanxiong he
     call putvalue(myid,aerom(i04aer),sx(ne),ex(ne),sy(ne),ey(ne), &
              IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
           enddo
           enddo
        endif



if(1==2) then
              ! for uv
              !IPS1 = igas + IAER*ISIZE + 1
              ips1 = ips1 + 1
              I03 = IP3MEM(K,NE)
             CALL PUTVALUE(MYID,U(I03),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(3,IPS),IPP,ATESTR(K,IPS1))
              !IPS1 = igas + IAER*ISIZE + 2
              ips1 = ips1 + 1
              I03 = IP3MEM(K,NE)
             CALL PUTVALUE(MYID,V(I03),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(3,IPS),IPP,ATESTR(K,IPS1))
endif


