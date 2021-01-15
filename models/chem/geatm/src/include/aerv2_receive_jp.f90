
           do iaersp = 1, naersp
           do iaerbin = 1, naerbin
             ips1 =  ips1 + 1
             i04aer= ip4mem_aer(k,iaerbin,iaersp,ne+1)
     call putvalue(myid,aerom(i04aer),sx(ne+1),ex(ne+1),sy(ne+1),ey(ne+1), &
              i,jp,atestR(k,ips1))
           enddo
           enddo

