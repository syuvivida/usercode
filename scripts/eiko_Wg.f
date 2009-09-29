      PROGRAM CHK018
C...Test program to generate ttbar events at Tevatron using PYTHIA
C...internal ttbar production subprocesses.
C...Ref: PYTHIA Tutorial, Fermilab, Dec 2004.

C ------ PREAMBLE: COMMON BLOCK DECLARATIONS ETC -----------------
C...All real arithmetic done in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

C ----------------- PYTHIA SETUP ---------------------------------
C...Number of events to generate
      NEV=10000
C...Select type of events to be generated: ttbar (using PYGIVE)
C...And use the new world average mt.
C... Genearting Pythia Wgamma

      CALL PYGIVE('MSEL=0')
      CALL PYGIVE('MSUB(20)=1')
C      CALL PYGIVE('CKIN(3)=10D0')   
      CALL PYGIVE('MDME(190,1)=0')
      CALL PYGIVE('MDME(191,1)=0')
      CALL PYGIVE('MDME(192,1)=0')
      CALL PYGIVE('MDME(194,1)=0')
      CALL PYGIVE('MDME(195,1)=0')
      CALL PYGIVE('MDME(196,1)=0')
      CALL PYGIVE('MDME(198,1)=0')
      CALL PYGIVE('MDME(199,1)=0')
      CALL PYGIVE('MDME(200,1)=0')
      CALL PYGIVE('MDME(206,1)=1')
      CALL PYGIVE('MDME(207,1)=0')
      CALL PYGIVE('MDME(208,1)=0')


C...Initialize PYTHIA for Tevatron ppbar collisions
C      ECM=10000D0
C      CALL PYINIT('CMS','p','p',ECM)
      ECM= 1960D0
      CALL PYINIT('CMS','p','pbar',ECM)
C...Initialize user stuff, e.g. book histograms etc.
      CALL MYSTUF(0,NEV)
      
      OPEN (unit = 1, file = "lgamma.data")
 
C ----------------- EVENT LOOP -----------------------------------
      DO 1000 IEV=1,NEV
C...Generate event
        CALL PYEVNT
C...Print out the event record of the first event
      IF (IEV.EQ.1) THEN 
          CALL PYLIST(2)
      ENDIF
C...Do event-by-event user stuff, e.g. fill histograms.
        CALL MYSTUF(1,IEV)
 1000 CONTINUE

      CLOSE(1)

C ----------------- FINALIZATION ---------------------------------
C...Print some info on cross sections and errors/warnings
      CALL PYSTAT(1)
C...Finalize my user stuff, e.g. close histogram file.
      CALL MYSTUF(2,NEV)
      
      END




C ================================================================



      SUBROUTINE MYSTUF(MODE,IEV)
C...User routine, this does not come with PYTHIA. 
C...To be written by user for whatever analysis is desired.
C...I've decided I want to compute the average pT of the top here.
C...MODE = 0: Initialization
C...     = 1: Event-by-event analysis
C...     = 2: Finalization
C...IEV: Lets us know which event we're analyzing.

C ------ PREAMBLE: COMMON BLOCK DECLARATIONS ETC -----------------
C...All real arithmetic done in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C...PYTHIA stuff (the event record etc)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...LOCAL variables
      DOUBLE PRECISION PTSUM,EVWT, GPX, GPY, GPZ, GE, LPX, LPY, LPZ, LE
      DOUBLE PRECISION NPX, NPY, NPZ, NE
      DOUBLE PRECISION G2PX,G2PY, G2PZ, G2E
      SAVE PTSUM,EVWT

C---------------------------------------------------------------------
C MODE = 0: INITIALIZATION.
C Compute event weight and zero cumulative pT counter.
      IF (MODE.EQ.0) THEN
        EVWT=1D0/DBLE(IEV)
        PTSUM=0D0
 
C---------------------------------------------------------------------
C MODE = 1: Event-by-event analysis
C Loop over the event record...
 
      ELSEIF(MODE.EQ.1) THEN
        I=1
        HASPHOTON=0
        GPX  = 0D0
        GPY  = 0D0
        GPZ  = 0D0
        GE   = 0D0
        LPX  = 0D0
        LPY  = 0D0
        LPZ  = 0D0
        LE   = 0D0
        NPX  = 0D0
        NPY  = 0D0
        NPZ  = 0D0
        NE   = 0D0
        G2PX = 0D0
        G2PY = 0D0
        G2PZ = 0D0
        G2E  = 0D0
        DO 10 I=1, N
C Look for photons and only the first photon         
        
        MOM = K(I,3)


        IF (I.LT.20.AND.K(I,1).EQ.1.AND.K(I,2).EQ.22.AND.HASPHOTON.EQ.0.
     &  AND.K(MOM,1).EQ.21.AND.K(MOM,2).EQ.22) 
     &  THEN
          GPX=PYP(I,1)
          GPY=PYP(I,2)
          GPZ=PYP(I,3)
          GE =PYP(I,4)
          HASPHOTON=1
C Look for photons and only the second photon       
        ELSEIF(I.LT.20.AND.K(I,1).EQ.1.AND.K(I,2).EQ.22.AND.HASPHOTON.
     &   EQ.1.AND.K(MOM,1).EQ.21.AND.(ABS(K(MOM,2)).EQ.11.OR.
     &   ABS(K(MOM,2)).EQ.13)) THEN
          G2PX=PYP(I,1)
          G2PY=PYP(I,2)
          G2PZ=PYP(I,3)
          G2E =PYP(I,4)
          HASPHOTON=2

          
C look for muons
        ELSEIF (I.LT.20.AND.K(I,1).EQ.1.AND.
     &           ((ABS(K(I,2)).EQ.12).OR.(ABS(K(I,2)).EQ.14))
           
     &           ) THEN
           NPX=PYP(I,1)
           NPY=PYP(I,2)
           NPZ=PYP(I,3)
           NE =PYP(I,4)

        ELSEIF (I.LT.20.AND.K(I,1).EQ.1.AND.
     &           ((ABS(K(I,2)).EQ.11).OR.(ABS(K(I,2)).EQ.13))           
     &           ) THEN
           LPX=PYP(I,1)
           LPY=PYP(I,2)
           LPZ=PYP(I,3)
           LE =PYP(I,4)
C           LPT=PYP(I,10)
C           LETA=PYP(I,19)
C        WRITE(1,'("electron: "2x,F10.3,3x,F10.3,3x, F10.3,3x, F10.3)')
        ENDIF  
C        PRINT *, HASPHOTON

 10       CONTINUE 
        WRITE(1,'(F12.5,1x,F12.5,1x, F12.5,1x, F12.5,1x,
     &  F12.5,1x,F12.5,1x,F12.5,1x,F12.5,1x,
     &  F12.5,1x,F12.5,1x,F12.5,1x,F12.5,1x,
     &  F12.5,1x,F12.5,1x,F12.5,1x,F12.5)')
     &  GPX, GPY, GPZ, GE, 
     &  LPX, LPY, LPZ, LE, 
     &  NPX, NPY, NPZ, NE,
     &  G2PX, G2PY, G2PZ, G2E
        
C---------------------------------------------------------------------
C MODE = 2: Finalization
C Compute and print average pT to output.
      ELSEIF(MODE.EQ.2) THEN
        AVGPT=PTSUM*EVWT
        UNCRT=SQRT(EVWT)*AVGPT
 
      ENDIF

      RETURN 
      END
