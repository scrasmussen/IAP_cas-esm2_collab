


! in k loop

          DO IG=1,iedgas
            ITSP=ITSP+1
            I04 = IP4MEM(K,IG,NE)
            CALL GETVALUE(MYID,gas(I04),SX(NE),EX(NE),SY(NE),EY(NE),&
                       IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP))
          ENDDO

          !For aerosols
          DO IA=1,IAER
          DO IS=1,ISIZE
           I05= IP5MEM(K,IS,IA,NE)
           ITSP=ITSP+1
          CALL GETVALUE(MYID,AER(I05),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP))
          ENDDO
          ENDDO

        !!!!!!!!!!!!!!!!!!!!
        ! for Source Mark
        if(ifsmt>0)then    ! checking 
         do idm=1,idmSet
         do ism=1,ismMax
            i0=ipSMmem(k,ism,idm,ne)
            itsp=itsp+1
         call getvalue(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),atestS(k,itsp))
         enddo
         enddo
        endif
        !!!!!!!!!!!!!!!!!!!!
        ! for dust and sea salt composition
!        if(ifseacom.eq.1) then
        do iduc = 1, nseacom
        do is =1 , isize
          i05c=ip5memcs(k,is,iduc,ne)
          itsp = itsp + 1
     call getvalue(myid,seacomp(i05c),sx(ne),ex(ne),sy(ne),ey(ne), &
             IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),atestS(k,itsp))
        enddo
        enddo
!        endif

!        if(ifdustcom.eq.1) then
        do iduc = 1, ndustcom
        do is =1 , isize
          i05c=ip5memc(k,is,iduc,ne)
          itsp = itsp + 1
     call getvalue(myid, dustcomp(i05c) ,sx(ne),ex(ne),sy(ne),ey(ne), &
             IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),atestS(k,itsp))
        enddo
        enddo
!        endif


        do iaersp = 1, naersp
        do iaerbin =1, naerbin
          i04aer = ip4mem_aer(k,iaerbin,iaersp,ne)
          itsp = itsp + 1
     call getvalue(myid, aerom(i04aer) ,sx(ne),ex(ne),sy(ne),ey(ne), &
             IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),atestS(k,itsp))
        enddo
        enddo

if(1==2) then
           ! FOR U,V 
           I03= IP3MEM(K,NE)
           ITSP=ITSP+1
           CALL GETVALUE(MYID,U(I03),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP))
           ITSP=ITSP+1
           CALL GETVALUE(MYID,V(I03),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP))

endif
