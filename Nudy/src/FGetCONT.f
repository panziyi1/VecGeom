      SUBROUTINE FGetCONT(LINE, C, NL, MTF)
      CHARACTER*80 LINE
      DOUBLE PRECISION C(2)
      INTEGER NL(4), MTF(4)
      READ(LINE,10) (C(L1),L1=1,2),(NL(L2),L2=1,4),(MTF(L),L=1,4)
 10   FORMAT(2E11.0,4I11,I4,I2,I3,I5) 
      END
