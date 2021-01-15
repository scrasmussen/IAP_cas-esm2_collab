
        do iaersp = 1, naersp
        do iaerbin =1, naerbin
          i04aer = ip4mem_aer(k,iaerbin,iaersp,ne)
          itsp = itsp + 1
     call getvalue(myid, aerom(i04aer) ,sx(ne),ex(ne),sy(ne),ey(ne), &
             i,j,atestS(k,itsp))
        enddo
        enddo

