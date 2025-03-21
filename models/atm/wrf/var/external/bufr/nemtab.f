      SUBROUTINE NEMTAB(LUN,NEMO,IDN,TAB,IRET)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    NEMTAB
C   PRGMMR: WOOLLEN          ORG: NP20       DATE: 1994-01-06
C
C ABSTRACT: THIS SUBROUTINE SEARCHES FOR MNEMONIC NEMO WITHIN THE
C   INTERNAL TABLE B AND D ARRAYS HOLDING THE DICTIONARY TABLE (ARRAYS
C   IN COMMON BLOCK /TABABD/) AND, IF FOUND, RETURNS INFORMATION ABOUT
C   THAT MNEMONIC FROM WITHIN THESE ARRAYS.  OTHERWISE, IT CHECKS
C   WHETHER NEMO IS A TABLE C OPERATOR DESCRIPTOR AND, IF SO, DIRECTLY
C   COMPUTES AND RETURNS SIMILAR INFORMATION ABOUT THAT DESCRIPTOR.
C   THIS SUBROUTINE MAY BE USEFUL TO APPLICATION PROGRAMS WHICH WANT
C   TO CHECK WHETHER A PARTICULAR MNEMONIC IS IN THE DICTIONARY.  IN
C   THIS CASE, BUFR ARCHIVE LIBRARY SUBROUTINE OPENBF MUST FIRST BE
C   CALLED TO STORE THE DICTIONARY TABLE INTERNALLY, AND BUFR ARCHIVE
C   LIBRARY SUBROUTINE STATUS MUST BE CALLED TO CONNECT THE LOGICAL
C   UNIT NUMBER FOR THE BUFR FILE OPENED IN OPENBF TO LUN.
C
C PROGRAM HISTORY LOG:
C 1994-01-06  J. WOOLLEN -- ORIGINAL AUTHOR
C 1995-06-28  J. WOOLLEN -- INCREASED THE SIZE OF INTERNAL BUFR TABLE
C                           ARRAYS IN ORDER TO HANDLE BIGGER FILES
C 1999-11-18  J. WOOLLEN -- THE NUMBER OF BUFR FILES WHICH CAN BE
C                           OPENED AT ONE TIME INCREASED FROM 10 TO 32
C                           (NECESSARY IN ORDER TO PROCESS MULTIPLE
C                           BUFR FILES UNDER THE MPI)
C 2000-09-19  J. WOOLLEN -- ADDED CAPABILITY TO ENCODE AND DECODE DATA
C                           USING THE OPERATOR DESCRIPTORS (BUFR TABLE
C                           C) FOR CHANGING WIDTH AND CHANGING SCALE
C 2003-11-04  J. ATOR    -- ADDED DOCUMENTATION
C 2003-11-04  S. BENDER  -- ADDED REMARKS/BUFRLIB ROUTINE
C                           INTERDEPENDENCIES
C 2003-11-04  D. KEYSER  -- UNIFIED/PORTABLE FOR WRF; ADDED HISTORY
C                           DOCUMENTATION
C 2005-11-29  J. ATOR    -- ADDED SUPPORT FOR 207 AND 208 OPERATORS
C
C USAGE:    CALL NEMTAB (LUN, NEMO, IDN, TAB, IRET)
C   INPUT ARGUMENT LIST:
C     LUN      - INTEGER: I/O STREAM INDEX INTO INTERNAL MEMORY ARRAYS
C     NEMO     - CHARACTER*(*): MNEMONIC TO SEARCH FOR
C
C   OUTPUT ARGUMENT LIST:
C     IDN      - INTEGER: BIT-WISE REPRESENTATION OF FXY VALUE
C                CORRESPONDING TO NEMO (IF NEMO WAS FOUND)
C     TAB      - CHARACTER*1: INTERNAL TABLE ARRAY IN WHICH NEMO WAS
C                FOUND:
C                     'B' = Table B array
C                     'C' = Table C array
C                     'D' = Table D array
C     IRET     - INTEGER: POSITIONAL INDEX OF NEMO WITHIN TAB
C                       0 = NEMO was not found within any of the Table
C                           B, C, or D arrays
C
C REMARKS:
C    THIS ROUTINE CALLS:        IFXY
C    THIS ROUTINE IS CALLED BY: CHEKSTAB CMSGINI  ELEMDX   MSGINI
C                               SEQSDX   TABSUB   UFBDMP   UFBQCD
C                               UFDUMP   UPFTBV
C                               Also called by application programs
C                               (see ABSTRACT).
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  PORTABLE TO ALL PLATFORMS
C
C$$$

      INCLUDE 'bufrlib.prm'

      COMMON /TABABD/ NTBA(0:NFILES),NTBB(0:NFILES),NTBD(0:NFILES),
     .                MTAB(MAXTBA,NFILES),IDNA(MAXTBA,NFILES,2),
     .                IDNB(MAXTBB,NFILES),IDND(MAXTBD,NFILES),
     .                TABA(MAXTBA,NFILES),TABB(MAXTBB,NFILES),
     .                TABD(MAXTBD,NFILES)

      CHARACTER*(*) NEMO
      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*8   NEMT
      CHARACTER*1   TAB
      LOGICAL       FOLVAL

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      FOLVAL = NEMO(1:1).EQ.'.'
      IRET = 0
      TAB = ' '

C  LOOK FOR NEMO IN TABLE B
C  ------------------------

      DO 1 I=1,NTBB(LUN)
      NEMT = TABB(I,LUN)(7:14)
      IF(NEMT.EQ.NEMO) THEN
         IDN  = IDNB(I,LUN)
         TAB  = 'B'
         IRET = I
         GOTO 100
      ELSEIF(FOLVAL.AND.NEMT(1:1).EQ.'.') THEN
         DO J=2,LEN(NEMT)
         IF(NEMT(J:J).NE.'.' .AND. NEMT(J:J).NE.NEMO(J:J)) GOTO 1
         ENDDO
         IDN  = IDNB(I,LUN)
         TAB  = 'B'
         IRET = I
         GOTO 100
      ENDIF
1     ENDDO

C  DON'T LOOK IN TABLE D FOR FOLLOWING VALUE-MNEMONICS
C  ---------------------------------------------------

      IF(FOLVAL) GOTO 100

C  LOOK IN TABLE D IF WE GOT THIS FAR
C  ----------------------------------

      DO I=1,NTBD(LUN)
      NEMT = TABD(I,LUN)(7:14)
      IF(NEMT.EQ.NEMO) THEN
         IDN  = IDND(I,LUN)
         TAB  = 'D'
         IRET = I
         GOTO 100
      ENDIF
      ENDDO

C  IF STILL NOTHING, CHECK HERE FOR TABLE C OPERATOR DESCRIPTORS
C  -------------------------------------------------------------

      IF(NEMO(1:3).EQ.'201' .OR.
     .   NEMO(1:3).EQ.'202' .OR.
     .   NEMO(1:3).EQ.'206' .OR.
     .   NEMO(1:3).EQ.'207' .OR.
     .   NEMO(1:3).EQ.'208' ) THEN
         READ(NEMO,'(1X,I2)') IRET
         IDN = IFXY(NEMO)
         TAB = 'C'
         GOTO 100
      ENDIF

C  EXIT
C  ----

100   RETURN
      END
