        subroutine getland_usgs2zhang03(landuse0,ilanduse0)

          real :: landuse0 ! the 25 catagories in MM5
          integer :: ilanduse0 ! the 11 catagories in RADM
!--------------------------------------Table----------------------
!C        Change usgs 24 landuse categories to zhang2003 26 categories
!C        -----------------------------------------------
!C        MM5        CATEGORIES                 CMAQ/RADM
!C        -----------------------------------------------
!C         1      URBAN LAND                          1
!C         2      DRYLAND CROPLAND & PASTURE          2
!C         3      IRRIGATED CROPLAND & PASTURE        2
!C         4      MIXED DRYLAND & IRRIGATED CROPLAND  2
!C         5      CROPLAND/GRASSLAND MOSAIC           2
!C         6      CROPLAND/WOODLAND MOSAIC           10
!C         7      GRASSLAND                           3
!C         8      SHRUBLAND                           3
!C         9      MIXED SHRUBLAND/GRASSLAND           3
!C        10      SAVANNAH                           10
!C        11      DECIDUOUS BROADLEAF FOREST          4
!C        12      DECIDUOUS NEEDLELEAF FOREST         4
!C        13      EVERGREEN BROADLEAF FOREST          5
!C        14      EVERGREEN NEEDLELEAF FOREST         5
!C        15      MIXED FOREST                        6
!C        16      WATER                               7
!C        17      HERBACEOUS WETLAND                  9
!C        18      WOODED WETLAND                      5
!C        19      BARREN OR SPARSELY VEGETATED        8
!C        20      HERBACEOUS TUNDRA                  11
!C        21      WOODED TUNDRA                       5
!C        22      MIXED TUNDRA                       10
!C        23      BARE GROUND TUNDRA                  8
!C        24      SNOW OR ICE                        11
!-----------------------------------------------------------
         select case(int(landuse0))
             case(0)
                 ilanduse0=0 ! no data
             case(1)
                 ilanduse0=21
             case(2)
                 ilanduse0=15 ! 
             case(3)
                 ilanduse0=20
             case(4)
                 ilanduse0=15
             case(5)
                 ilanduse0=13
             case(6)
                 ilanduse0=26
             case(7)
                 ilanduse0=14
             case(8)
                 ilanduse0=10
             case(9)
                 ilanduse0=10
             case(10)
                 ilanduse0=14
             case(11)
                 ilanduse0=7
             case(12)
                 ilanduse0=6
             case(13)
                 ilanduse0=5
             case(14)
                 ilanduse0=4
             case(15)
                 ilanduse0=25
             case(16)
                 ilanduse0=1
             case(17)
                 ilanduse0=23
             case(18)
                 ilanduse0=5
             case(19)
                 ilanduse0=24
             case(20)
                 ilanduse0=22
             case(21)
                 ilanduse0=4
             case(22)
                 ilanduse0=22
             case(23)
                 ilanduse0=22
             case(24)
                 ilanduse0=2
          end select
        return

        end



