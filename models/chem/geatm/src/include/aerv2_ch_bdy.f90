
   do iaersp  = 1, naersp
   do iaerbin = 1, naerbin
     i04aer = ip4mem_aer (k,iaerbin, iaersp, ne)
     call exchng2 ( myid, aerom(i04aer), sx(ne), ex(ne), sy(ne),ey(ne),  &
                comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
                nbrtop(ne), nbrbottom(ne))
   enddo
   enddo

