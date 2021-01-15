c----CAMx v5.40 111010
c  
c     CHMDAT.COM contains all chemistry variables 
c                            
c     Copyright 1996 - 2010
c     ENVIRON International Corporation
c            
c     Modifications:  
c       4/4/00     Added aerosol deposition variables to /aerochm/
c       1/9/02     Aerosol size cut points and density now a function of
c                  species
c       8/20/02    Added minimum CWC to define presence of clouds
c      12/12/02    Expanded species list for Mechanism 4
c       1/10/03    Added array for deposition output species names
c      10/14/04    Modifications for mechanism 10
c      11/30/06    Created by taking Static variables out of chmstry.com
c      12/29/06 -bkoo-     Expanded species list for the updated SOA scheme
c      01/08/07 -bkoo-     Expanded species list for Mechanism 6 (CB05)
c                          Now ipigrxn is a scalar pointer to NO self-reaction
c      07/04/07 -bkoo-     Added pointer to hydrolysis of N2O5
c      04/28/08 -gyarwood- Added EBI solver
c 
c-----------------------------------------------------------------------
c     Parameters for some of the switches:
c
c     CDCMC  -- string for requesting the CMC chemistry solver
c     CDEBI  -- string for requesting the EBI chemistry solver
c     CDIEH  -- string for requesting the IEH chemistry solver
c     CDLSOD -- string for requesting the LSODE chemistry solver
c     IDCMC  -- code for using the CMC chemistry solver
c     IDEBI  -- code for using the EBI chemistry solver
c     IDIEH  -- code for using the IEH chemistry solver
c     IDLSOD -- code for using the LSODE chemistry solver
c
c-----------------------------------------------------------------------
c
      character*10 CDCMC
      character*10 CDEBI
      character*10 CDIEH
      character*10 CDLSOD
      integer      IDCMC
      integer      IDEBI
      integer      IDIEH
      integer      IDLSOD
c
      parameter( CDCMC  = "CMC       " )
      parameter( CDEBI  = "EBI       " )
      parameter( CDIEH  = "IEH       " )
      parameter( CDLSOD = "LSODE     " )
      parameter( IDCMC  = 1 )
      parameter( IDEBI  = 2 )
      parameter( IDIEH  = 3 )
      parameter( IDLSOD = 4 )
c 
c-----------------------------------------------------------------------
c    Variables for the number of species in input files:
c
c    ngas   --  number of gas species being modeled
c    naero  --  number of aersol species being modeled
c    nspec  --  total number of modeled species
c    nrad   --  number of radical species being modeled
c    nreact --  number of chemical reactions
c    nspfst --  number of "fast" species -- handled by the fast solver
c    iessrad--  number of radicals in steady state for IEH solver
c    idmech --  the code which determines which chemical mechanism is used
c    idsolv --  the code which determines which chemstry solver to use
c    navspc --  number of species to write to output average file
c    ndepspc--  number of species to write to deposition file
c    nicspc --  number of species in the initial conditions file
c    nbcspc --  number of species in the boundary conditions file
c    nptspc --  number of species in the point source emissions file
c
c-----------------------------------------------------------------------
c
       integer   ngas
       integer   naero
       integer   nspec
       integer   nrad
       integer   nreact
       integer   nspfst
       integer   iessrad
       integer   idmech
       integer   idsolv
       integer   navspc
       integer   ndepspc
       integer   nicspc
       integer   nbcspc
       integer   nptspc
c
      common /chm1/ ngas, naero, nspec, nrad, nreact, nspfst, iessrad,
     &                         idmech, idsolv, navspc, ndepspc, nicspc,
     &                                                    nbcspc, nptspc
c
c-----------------------------------------------------------------------
c
c    aeropt --  keyword for aerosol scheme ('NONE', 'INERT', 'CF' or 'CMU')
c
c-----------------------------------------------------------------------
c
      character*10 aeropt
c
      common /chm2/ aeropt
c
c-----------------------------------------------------------------------
c     Variables for photolysis data:
c
c     nphot1   -- number of primary photolysis reactions
c     nphot2   -- number of secondary (scaled) photolysis reactions
c
c-----------------------------------------------------------------------
c
      integer nphot1
      integer nphot2
c
      common /photmap/ nphot1, nphot2
c 
c-----------------------------------------------------------------------
c     Variables for vertical profiles of HCl and CL2 for Hg chemistry
c
c     NHTCL   --  number of heights for concentrations 
c     htcl    --  heights (km) for the top of each concentration range
c     hclprof --  profile of HCl (ppm) with height
c     cl2day  --  profile of Cl2 (ppm) with height during day over ocean
c     cl2nite --  profile of Cl2 (ppm) with height at night over ocean
c
c-----------------------------------------------------------------------
c
      integer NHTCL
c
      parameter ( NHTCL = 6 )
c
      real htcl(NHTCL)
      real hclprof(NHTCL)
      real cl2day(NHTCL)
      real cl2nite(NHTCL)
c
      common /clprofs/ htcl, hclprof, cl2day, cl2nite
c 
c-----------------------------------------------------------------------
c     Parameters for mapping input species to internal model order:
c
c      NSPNAM  -- number of recognizable state species
c
      integer   NSPNAM
c
      parameter ( NSPNAM  =   372 )
c 
c-----------------------------------------------------------------------
c     Variables for mapping input species to internal model order:
c
c     kmap     -- mapping of species on chemistry parameters file to
c                 internal order
c
c-----------------------------------------------------------------------
c
      integer kmap(NSPNAM)
c
      common /kname/ kmap
c
      integer       ko1d,      ko,     kcl
      integer       kclo,     koh,    kho2
      integer      kc2o3,    kxo2,   kxo2n
      integer      kcxo3,   kmeo2,    kto2
      integer       kror,    kcro,   kro2r
      integer      kr2o2,   kro2n,   kcco3
      integer      krco3,   kmco3,   kbzco
      integer      kcxo2,   khco3,   ktbuo
      integer       kbzo,   kbzno,      ki
      integer        kno,    kno2,     ko3
      integer       kpan,   kpanx,   kcres
      integer      kpan2,   kmpan,   kpbzn
      integer      knphe,   krno3,   kdcb2
      integer      kdcb3,   khno4,   kbzo2
      integer      kcao2,   kcrn2,   kcrno
      integer      kepx2,   kiso2,   kopan
      integer      kopo3,    kro2,   kxlo2
      integer      kxo2h,   kaacd,   kacet
      integer      kald2,   kaldx,   kalk1
      integer      kalk2,   kalk3,   kalk4
      integer      kalk5,   karo1,   karo2
      integer      kbacl,   kbald,   kbcl1
      integer      kbcl2,   kbenz,   kbuta
      integer      kcat1,   kccho,   kccrs
      integer       kcg1,    kcg2,    kcg3
      integer       kcg4,    kcg5,    kcg6
      integer       kcg7,   kch3i,    kcl2
      integer        kco,   kco2h,   kco3h
      integer      kcooh,   kcprm,   kcron
      integer      kcrpx,   kdcb1,   kepox
      integer       keth,   ketha,   kethe
      integer      kethy,   ketoh,   kfacd
      integer      kfcrs,   kfmcl,   kform
      integer      kfprm,    kgly,   kglyd
      integer      kh2o2,   khc2h,   khcho
      integer       khcl,    khg0,    khg2
      integer       khgp,     khi,   khio3
      integer      khno3,   kho2h,   khocl
      integer       khoi,   khono,     ki2
      integer      kibac,   kibuo,   kicl1
      integer      kicl2,    kino,   kino2
      integer      kino3,   kintr,     kio
      integer      kiole,   kisop,    kisp
      integer      kispd,   kispx,   kixoy
      integer       kket,   kmbut,    kmek
      integer      kmeoh,   kmepx,   kmeth
      integer      kmgly,   kmtbe,    kmvk
      integer      kn2o5,   kn3cl,     kna
      integer      knbui,    knh3,    kno3
      integer      kntcl,    kntr,    koio
      integer       kole,   kole1,   kole2
      integer      kopen,   kpacd,    kpar
      integer       kpcl,    kpec,   kph2o
      integer      kphen,    kpna,   kpnh4
      integer      kpno3,    kpoa,   kprod
      integer      kprpa,   kpso4,   krc2h
      integer      krc3h,   krcho,   krooh
      integer       kso2,   ksoa1,   ksoa2
      integer      ksoa3,   ksoa4,   ksoa5
      integer      ksoa6,   ksoa7,   ksopa
      integer      ksopb,    ksqt,     kss
      integer      ksulf,   kterp,    ktol
      integer      ktola,    ktrp,     kxn
      integer      kxopn,    kxyl,   kxyla
      integer    ksoa1_1, ksoa1_2, ksoa1_3
      integer    ksoa1_4, ksoa1_5, ksoa1_6
      integer    ksoa1_7, ksoa1_8, ksoa1_9
      integer   ksoa1_10, ksoa2_1, ksoa2_2
      integer    ksoa2_3, ksoa2_4, ksoa2_5
      integer    ksoa2_6, ksoa2_7, ksoa2_8
      integer    ksoa2_9,ksoa2_10, ksoa3_1
      integer    ksoa3_2, ksoa3_3, ksoa3_4
      integer    ksoa3_5, ksoa3_6, ksoa3_7
      integer    ksoa3_8, ksoa3_9,ksoa3_10
      integer    ksoa4_1, ksoa4_2, ksoa4_3
      integer    ksoa4_4, ksoa4_5, ksoa4_6
      integer    ksoa4_7, ksoa4_8, ksoa4_9
      integer   ksoa4_10, ksoa5_1, ksoa5_2
      integer    ksoa5_3, ksoa5_4, ksoa5_5
      integer    ksoa5_6, ksoa5_7, ksoa5_8
      integer    ksoa5_9,ksoa5_10, ksoa6_1
      integer    ksoa6_2, ksoa6_3, ksoa6_4
      integer    ksoa6_5, ksoa6_6, ksoa6_7
      integer    ksoa6_8, ksoa6_9,ksoa6_10
      integer    ksoa7_1, ksoa7_2, ksoa7_3
      integer    ksoa7_4, ksoa7_5, ksoa7_6
      integer    ksoa7_7, ksoa7_8, ksoa7_9
      integer   ksoa7_10, ksopa_1, ksopa_2
      integer    ksopa_3, ksopa_4, ksopa_5
      integer    ksopa_6, ksopa_7, ksopa_8
      integer    ksopa_9,ksopa_10, ksopb_1
      integer    ksopb_2, ksopb_3, ksopb_4
      integer    ksopb_5, ksopb_6, ksopb_7
      integer    ksopb_8, ksopb_9,ksopb_10
      integer     kpoa_1,  kpoa_2,  kpoa_3
      integer     kpoa_4,  kpoa_5,  kpoa_6
      integer     kpoa_7,  kpoa_8,  kpoa_9
      integer    kpoa_10,  kpec_1,  kpec_2
      integer     kpec_3,  kpec_4,  kpec_5
      integer     kpec_6,  kpec_7,  kpec_8
      integer     kpec_9, kpec_10, kcrst_1
      integer    kcrst_2, kcrst_3, kcrst_4
      integer    kcrst_5, kcrst_6, kcrst_7
      integer    kcrst_8, kcrst_9,kcrst_10
      integer    kph2o_1, kph2o_2, kph2o_3
      integer    kph2o_4, kph2o_5, kph2o_6
      integer    kph2o_7, kph2o_8, kph2o_9
      integer   kph2o_10,  kpcl_1,  kpcl_2
      integer     kpcl_3,  kpcl_4,  kpcl_5
      integer     kpcl_6,  kpcl_7,  kpcl_8
      integer     kpcl_9, kpcl_10,   kna_1
      integer      kna_2,   kna_3,   kna_4
      integer      kna_5,   kna_6,   kna_7
      integer      kna_8,   kna_9,  kna_10
      integer    kpnh4_1, kpnh4_2, kpnh4_3
      integer    kpnh4_4, kpnh4_5, kpnh4_6
      integer    kpnh4_7, kpnh4_8, kpnh4_9
      integer   kpnh4_10, kpno3_1, kpno3_2
      integer    kpno3_3, kpno3_4, kpno3_5
      integer    kpno3_6, kpno3_7, kpno3_8
      integer    kpno3_9,kpno3_10, kpso4_1
      integer    kpso4_2, kpso4_3, kpso4_4
      integer    kpso4_5, kpso4_6, kpso4_7
      integer    kpso4_8, kpso4_9,kpso4_10
c
      equivalence
     &   (kmap(  1),    ko1d),(kmap(  2),      ko),(kmap(  3),     kcl),
     &   (kmap(  4),    kclo),(kmap(  5),     koh),(kmap(  6),    kho2),
     &   (kmap(  7),   kc2o3),(kmap(  8),    kxo2),(kmap(  9),   kxo2n),
     &   (kmap( 10),   kcxo3),(kmap( 11),   kmeo2),(kmap( 12),    kto2),
     &   (kmap( 13),    kror),(kmap( 14),    kcro),(kmap( 15),   kro2r),
     &   (kmap( 16),   kr2o2),(kmap( 17),   kro2n),(kmap( 18),   kcco3),
     &   (kmap( 19),   krco3),(kmap( 20),   kmco3),(kmap( 21),   kbzco),
     &   (kmap( 22),   kcxo2),(kmap( 23),   khco3),(kmap( 24),   ktbuo),
     &   (kmap( 25),    kbzo),(kmap( 26),   kbzno),(kmap( 27),      ki),
     &   (kmap( 28),     kno),(kmap( 29),    kno2),(kmap( 30),     ko3),
     &   (kmap( 31),    kpan),(kmap( 32),   kpanx),(kmap( 33),   kcres),
     &   (kmap( 34),   kpan2),(kmap( 35),   kmpan),(kmap( 36),   kpbzn),
     &   (kmap( 37),   knphe),(kmap( 38),   krno3),(kmap( 39),   kdcb2),
     &   (kmap( 40),   kdcb3),(kmap( 41),   khno4),(kmap( 42),   kbzo2),
     &   (kmap( 43),   kcao2),(kmap( 44),   kcrn2),(kmap( 45),   kcrno),
     &   (kmap( 46),   kepx2),(kmap( 47),   kiso2),(kmap( 48),   kopan),
     &   (kmap( 49),   kopo3),(kmap( 50),    kro2),(kmap( 51),   kxlo2),
     &   (kmap( 52),   kxo2h),(kmap( 53),   kaacd),(kmap( 54),   kacet),
     &   (kmap( 55),   kald2),(kmap( 56),   kaldx),(kmap( 57),   kalk1),
     &   (kmap( 58),   kalk2),(kmap( 59),   kalk3),(kmap( 60),   kalk4),
     &   (kmap( 61),   kalk5),(kmap( 62),   karo1),(kmap( 63),   karo2),
     &   (kmap( 64),   kbacl),(kmap( 65),   kbald),(kmap( 66),   kbcl1),
     &   (kmap( 67),   kbcl2),(kmap( 68),   kbenz),(kmap( 69),   kbuta),
     &   (kmap( 70),   kcat1),(kmap( 71),   kccho),(kmap( 72),   kccrs),
     &   (kmap( 73),    kcg1),(kmap( 74),    kcg2),(kmap( 75),    kcg3),
     &   (kmap( 76),    kcg4),(kmap( 77),    kcg5),(kmap( 78),    kcg6),
     &   (kmap( 79),    kcg7),(kmap( 80),   kch3i),(kmap( 81),    kcl2),
     &   (kmap( 82),     kco),(kmap( 83),   kco2h),(kmap( 84),   kco3h),
     &   (kmap( 85),   kcooh),(kmap( 86),   kcprm),(kmap( 87),   kcron),
     &   (kmap( 88),   kcrpx),(kmap( 89),   kdcb1),(kmap( 90),   kepox),
     &   (kmap( 91),    keth),(kmap( 92),   ketha),(kmap( 93),   kethe),
     &   (kmap( 94),   kethy),(kmap( 95),   ketoh),(kmap( 96),   kfacd),
     &   (kmap( 97),   kfcrs),(kmap( 98),   kfmcl),(kmap( 99),   kform),
     &   (kmap(100),   kfprm),(kmap(101),    kgly),(kmap(102),   kglyd),
     &   (kmap(103),   kh2o2),(kmap(104),   khc2h),(kmap(105),   khcho),
     &   (kmap(106),    khcl),(kmap(107),    khg0),(kmap(108),    khg2),
     &   (kmap(109),    khgp),(kmap(110),     khi),(kmap(111),   khio3),
     &   (kmap(112),   khno3),(kmap(113),   kho2h),(kmap(114),   khocl),
     &   (kmap(115),    khoi),(kmap(116),   khono),(kmap(117),     ki2),
     &   (kmap(118),   kibac),(kmap(119),   kibuo),(kmap(120),   kicl1),
     &   (kmap(121),   kicl2),(kmap(122),    kino),(kmap(123),   kino2),
     &   (kmap(124),   kino3),(kmap(125),   kintr),(kmap(126),     kio),
     &   (kmap(127),   kiole),(kmap(128),   kisop),(kmap(129),    kisp),
     &   (kmap(130),   kispd),(kmap(131),   kispx),(kmap(132),   kixoy),
     &   (kmap(133),    kket),(kmap(134),   kmbut),(kmap(135),    kmek),
     &   (kmap(136),   kmeoh),(kmap(137),   kmepx),(kmap(138),   kmeth),
     &   (kmap(139),   kmgly),(kmap(140),   kmtbe),(kmap(141),    kmvk),
     &   (kmap(142),   kn2o5),(kmap(143),   kn3cl),(kmap(144),     kna),
     &   (kmap(145),   knbui),(kmap(146),    knh3),(kmap(147),    kno3),
     &   (kmap(148),   kntcl),(kmap(149),    kntr),(kmap(150),    koio),
     &   (kmap(151),    kole),(kmap(152),   kole1),(kmap(153),   kole2),
     &   (kmap(154),   kopen),(kmap(155),   kpacd),(kmap(156),    kpar),
     &   (kmap(157),    kpcl),(kmap(158),    kpec),(kmap(159),   kph2o),
     &   (kmap(160),   kphen),(kmap(161),    kpna),(kmap(162),   kpnh4),
     &   (kmap(163),   kpno3),(kmap(164),    kpoa),(kmap(165),   kprod),
     &   (kmap(166),   kprpa),(kmap(167),   kpso4),(kmap(168),   krc2h),
     &   (kmap(169),   krc3h),(kmap(170),   krcho),(kmap(171),   krooh),
     &   (kmap(172),    kso2),(kmap(173),   ksoa1),(kmap(174),   ksoa2),
     &   (kmap(175),   ksoa3),(kmap(176),   ksoa4),(kmap(177),   ksoa5),
     &   (kmap(178),   ksoa6),(kmap(179),   ksoa7),(kmap(180),   ksopa),
     &   (kmap(181),   ksopb),(kmap(182),    ksqt),(kmap(183),     kss),
     &   (kmap(184),   ksulf),(kmap(185),   kterp),(kmap(186),    ktol),
     &   (kmap(187),   ktola),(kmap(188),    ktrp),(kmap(189),     kxn),
     &   (kmap(190),   kxopn),(kmap(191),    kxyl),(kmap(192),   kxyla),
     &   (kmap(193), ksoa1_1),(kmap(194), ksoa1_2),(kmap(195), ksoa1_3),
     &   (kmap(196), ksoa1_4),(kmap(197), ksoa1_5),(kmap(198), ksoa1_6),
     &   (kmap(199), ksoa1_7),(kmap(200), ksoa1_8),(kmap(201), ksoa1_9),
     &   (kmap(202),ksoa1_10),(kmap(203), ksoa2_1),(kmap(204), ksoa2_2),
     &   (kmap(205), ksoa2_3),(kmap(206), ksoa2_4),(kmap(207), ksoa2_5),
     &   (kmap(208), ksoa2_6),(kmap(209), ksoa2_7),(kmap(210), ksoa2_8),
     &   (kmap(211), ksoa2_9),(kmap(212),ksoa2_10),(kmap(213), ksoa3_1),
     &   (kmap(214), ksoa3_2),(kmap(215), ksoa3_3),(kmap(216), ksoa3_4),
     &   (kmap(217), ksoa3_5),(kmap(218), ksoa3_6),(kmap(219), ksoa3_7),
     &   (kmap(220), ksoa3_8),(kmap(221), ksoa3_9),(kmap(222),ksoa3_10),
     &   (kmap(223), ksoa4_1),(kmap(224), ksoa4_2),(kmap(225), ksoa4_3),
     &   (kmap(226), ksoa4_4),(kmap(227), ksoa4_5),(kmap(228), ksoa4_6),
     &   (kmap(229), ksoa4_7),(kmap(230), ksoa4_8),(kmap(231), ksoa4_9),
     &   (kmap(232),ksoa4_10),(kmap(233), ksoa5_1),(kmap(234), ksoa5_2),
     &   (kmap(235), ksoa5_3),(kmap(236), ksoa5_4),(kmap(237), ksoa5_5),
     &   (kmap(238), ksoa5_6),(kmap(239), ksoa5_7),(kmap(240), ksoa5_8),
     &   (kmap(241), ksoa5_9),(kmap(242),ksoa5_10),(kmap(243), ksoa6_1),
     &   (kmap(244), ksoa6_2),(kmap(245), ksoa6_3),(kmap(246), ksoa6_4),
     &   (kmap(247), ksoa6_5),(kmap(248), ksoa6_6),(kmap(249), ksoa6_7),
     &   (kmap(250), ksoa6_8),(kmap(251), ksoa6_9),(kmap(252),ksoa6_10),
     &   (kmap(253), ksoa7_1),(kmap(254), ksoa7_2),(kmap(255), ksoa7_3),
     &   (kmap(256), ksoa7_4),(kmap(257), ksoa7_5),(kmap(258), ksoa7_6),
     &   (kmap(259), ksoa7_7),(kmap(260), ksoa7_8),(kmap(261), ksoa7_9),
     &   (kmap(262),ksoa7_10),(kmap(263), ksopa_1),(kmap(264), ksopa_2),
     &   (kmap(265), ksopa_3),(kmap(266), ksopa_4),(kmap(267), ksopa_5),
     &   (kmap(268), ksopa_6),(kmap(269), ksopa_7),(kmap(270), ksopa_8)
      equivalence
     &   (kmap(271), ksopa_9),(kmap(272),ksopa_10),(kmap(273), ksopb_1),
     &   (kmap(274), ksopb_2),(kmap(275), ksopb_3),(kmap(276), ksopb_4),
     &   (kmap(277), ksopb_5),(kmap(278), ksopb_6),(kmap(279), ksopb_7),
     &   (kmap(280), ksopb_8),(kmap(281), ksopb_9),(kmap(282),ksopb_10),
     &   (kmap(283),  kpoa_1),(kmap(284),  kpoa_2),(kmap(285),  kpoa_3),
     &   (kmap(286),  kpoa_4),(kmap(287),  kpoa_5),(kmap(288),  kpoa_6),
     &   (kmap(289),  kpoa_7),(kmap(290),  kpoa_8),(kmap(291),  kpoa_9),
     &   (kmap(292), kpoa_10),(kmap(293),  kpec_1),(kmap(294),  kpec_2),
     &   (kmap(295),  kpec_3),(kmap(296),  kpec_4),(kmap(297),  kpec_5),
     &   (kmap(298),  kpec_6),(kmap(299),  kpec_7),(kmap(300),  kpec_8),
     &   (kmap(301),  kpec_9),(kmap(302), kpec_10),(kmap(303), kcrst_1),
     &   (kmap(304), kcrst_2),(kmap(305), kcrst_3),(kmap(306), kcrst_4),
     &   (kmap(307), kcrst_5),(kmap(308), kcrst_6),(kmap(309), kcrst_7),
     &   (kmap(310), kcrst_8),(kmap(311), kcrst_9),(kmap(312),kcrst_10),
     &   (kmap(313), kph2o_1),(kmap(314), kph2o_2),(kmap(315), kph2o_3),
     &   (kmap(316), kph2o_4),(kmap(317), kph2o_5),(kmap(318), kph2o_6),
     &   (kmap(319), kph2o_7),(kmap(320), kph2o_8),(kmap(321), kph2o_9),
     &   (kmap(322),kph2o_10),(kmap(323),  kpcl_1),(kmap(324),  kpcl_2),
     &   (kmap(325),  kpcl_3),(kmap(326),  kpcl_4),(kmap(327),  kpcl_5),
     &   (kmap(328),  kpcl_6),(kmap(329),  kpcl_7),(kmap(330),  kpcl_8),
     &   (kmap(331),  kpcl_9),(kmap(332), kpcl_10),(kmap(333),   kna_1),
     &   (kmap(334),   kna_2),(kmap(335),   kna_3),(kmap(336),   kna_4),
     &   (kmap(337),   kna_5),(kmap(338),   kna_6),(kmap(339),   kna_7),
     &   (kmap(340),   kna_8),(kmap(341),   kna_9),(kmap(342),  kna_10),
     &   (kmap(343), kpnh4_1),(kmap(344), kpnh4_2),(kmap(345), kpnh4_3),
     &   (kmap(346), kpnh4_4),(kmap(347), kpnh4_5),(kmap(348), kpnh4_6),
     &   (kmap(349), kpnh4_7),(kmap(350), kpnh4_8),(kmap(351), kpnh4_9),
     &   (kmap(352),kpnh4_10),(kmap(353), kpno3_1),(kmap(354), kpno3_2),
     &   (kmap(355), kpno3_3),(kmap(356), kpno3_4),(kmap(357), kpno3_5),
     &   (kmap(358), kpno3_6),(kmap(359), kpno3_7),(kmap(360), kpno3_8),
     &   (kmap(361), kpno3_9),(kmap(362),kpno3_10),(kmap(363), kpso4_1),
     &   (kmap(364), kpso4_2),(kmap(365), kpso4_3),(kmap(366), kpso4_4),
     &   (kmap(367), kpso4_5),(kmap(368), kpso4_6),(kmap(369), kpso4_7),
     &   (kmap(370), kpso4_8),(kmap(371), kpso4_9),(kmap(372),kpso4_10)
c
c-----------------------------------------------------------------------
c     Variables to define parameters for each chemical species:
c
c     henso20  -- Henry's Law constant at STP for SO2 (molar/atm)
c     tfactso2 -- Temperature dependence of SO2 Henry's Law constant (1/K)
c     nbin     -- Number of aerosol size bins
c     cwmin    -- Minimum cloud water threshold (g/m3)
c     tamin    -- Cloud water freezing threshold (K)
c
c-----------------------------------------------------------------------
c
      real henso20
      real tfactso2
      real cwmin
      real tamin
c
      common /depchm/ henso20, tfactso2, cwmin, tamin
c
      integer nbin
c
      common /aerochm/ nbin
c
c-----------------------------------------------------------------------
c     Pointer to hydrolysis of N2O5
c                    N2O5 + H2O -> 2 HNO3
c
      integer ihydrxn
      common /hydrxn/ ihydrxn
c
c-----------------------------------------------------------------------
c     Pointer used to lookup pig chemistry rate constant
c
c     ipigrxn  -- pointer to NO self-reaction
c                    NO + NO -> 2 NO2
c
      integer ipigrxn
c
      common /pigrxn/ ipigrxn
c
c----------------------------------------------------------------------
c    Variables for controlling calls to aerosol routines
c
c     dtaero     -- user input time interval (min) between calls to aerosol 
c                   routines
c     dt_aero    -- adjusted (for I/O frequency) time interval (min) between
c                   calls to aerosol routines
c
      real dtaero
      real dt_aero
c
      common /aero_t/ dtaero, dt_aero
