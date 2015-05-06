

C--------------------------------------------------------
C  COMPUTE PARTIAL CORRELATION COEFFICIENT AND
C       SIGNIFICANCE FOR CENSORED DATA
C
C               AUGUST 1995
C 
C  THE CODE IS BASED ON THE METHODOLOGY PRESENTED IN 
C     'A test for partial correlation with censored 
C                astronomical data'
C                      BY
C            M.G.Akritas and J.Siebert 
C  Monthly Notices of the Royal Astronomical Society
C               278, 919-924, 1996
C-------------------------------------------------------
      program partial_tau

      common /data/ dat(500,3),idat(500,3)
      common ntot
      common /kx/ k1,k2,k3

C-------------------------------------------------------------
C  INPUT DATA FILE CALLED 'DATA'.
C  CURRENTLY THE DATA FORMAT IS FIXED TO 3(f10.4,1x,i2,1x).
C  1ST, 3RD AND 5TH COLUMN ARE INDEPENDENT, DEPENDENT AND TEST
C  VARIABLE, RESPECTIVELY. 2ND, 4TH AND 6TH COLUMN DENOTE
C  CENSORING WITH 1 = DETECTION, 0= UPPER LIMIT
C  EXAMPLE:
C  '   26.9800  1    44.4340  0    -1.0714  1 '
C------------------------------------------------------------- 

      open(10,file='data',status='old')

C-------------------------------------------------
C READ IN DATA:
C     DAT(I,K)  = MEASUREMENT I OF VARIABLE K
C     IDAT(I,K) = CENSORING INDICATOR FOR DATA POINT (I,K)
C                 DETECTION   --> IDAT(I,K)=1
C                 UPPER LIMIT --> IDAT(I,K)=0
C--------------------------------------------------
 
      i=1
1     read(10,*,end=99) dat(i,1),idat(i,1),dat(i,2),
     #     idat(i,2),dat(i,3),idat(i,3)
110   format(3(f10.4,1x,i2,1x))
      dat(i,1)=-dat(i,1)     ! CHANGE TO RIGHT CENSORING
      dat(i,2)=-dat(i,2)     ! CHANGE TO RIGHT CENSORING
      dat(i,3)=-dat(i,3)     ! CHANGE TO RIGHT CENSORING
      i=i+1
      goto 1
99    ntot=i-1

      k1=1       ! INDEPENDENT VARIABLE = 1.COL OF DAT
      k2=2       ! DEPENDENT VARIABLE   = 2.COL OF DAT
      k3=3       ! THIRD VARIABLE       = 3.COL OF DAT

      call tau123(res)     ! COMPUTE PARTIAL KENDALLS TAU

      write(6,*) 'Tau(1,2):',tau(k1,k2)
      write(6,*) 'Tau(1,3):',tau(k1,k3)
      write(6,*) 'Tau(2,3):',tau(k2,k3)
      write(6,*) '--> Partial Kendalls tau:', res
      write(6,*) '  '
      write(6,*) 'Calculating variance...this takes some time....'
      write(6,*) '  '

      call sigma(sig)      ! COMPUTE VARIANCE

      write(6,*) 'Square root of variance (sigma):',sig
      write(6,*) '  '
      write(6,*) 'Z-score: ',abs(res/sig)
      if(abs(res/sig).gt.1.96) then
        write(6,*) 'Zero partial correlation rejected at level 0.05'
      else
       write(6,*) 'Null hypothesis cannot be rejected!'
       write(6,*) '(--> No correlation present, if influence of
     #third variable is excluded)'
      endif

      stop
      end

C------------------------------------------------------
C-------- SUBROUTINES AND FUNCTIONS -------------------
C------------------------------------------------------

C------------ TAU123 ---------------------------------
C-------- PARTIAL KENDALLS TAU -----------------------

      subroutine tau123(res)
      common /kx/ k1,k2,k3

      res= (tau(k1,k2)-tau(k1,k3)*tau(k2,k3))/
     # sqrt((1.-tau(k1,k3)**2)*(1.-tau(k2,k3)**2))
      end

C------------ SIGMA -----------------------------------
C-------- VARIANCE OF STATISTIC -----------------------

      subroutine sigma(sigres)
      common ntot
      common /kx/ k1,k2,k3

      sig2=an( )/(ntot*(1.-tau(k1,k3)**2)*(1.-tau(k2,k3)**2))
      sigres=sqrt(sig2)
      end

C------------ AN ---------------------------------------
C------- COMPUTES VALUE FOR A_N -------------------------

      function an( )
      double precision aasum(500)
      common ntot
      common /data/ dat(500,3),idat(500,3)
      common /kx/ k1,k2,k3
      c1=16./(float(ntot)-1.)
      c2=6./((float(ntot)-1.)*(float(ntot)-2.)*(float(ntot)-3.))
      asum=0.0
      ave = 0.0
      do 5 i=1,ntot
        aasum(i)=0.0
 5    continue
      do 10 i1=1,ntot     ! OUTER SUMMATION (I1)
      write(6,*) i1
        do 11 j1=1,ntot-2         ! INNER SUMMATION WITH
          if(j1.eq.i1) goto 11    ! J1<I2<J2 AND ALL .NE. I1
          do 12 j2=j1+2,ntot      !
            if(j2.eq.i1) goto 12  !
            do 13 i2=j1+1,j2-1    !
            if(i2.eq.i1) goto 13  !
            cj1=- idat(j1,k1)
            if(dat(i1,k1).lt.dat(j1,k1)) cj1=idat(i1,k1)
            cj2=- idat(j1,k2)
            if(dat(i1,k2).lt.dat(j1,k2)) cj2=idat(i1,k2)
            cj3=- idat(j1,k3)
            if(dat(i1,k3).lt.dat(j1,k3)) cj3=idat(i1,k3)
            cj4=- idat(j2,k2)
            if(dat(i2,k2).lt.dat(j2,k2)) cj4=idat(i2,k2)
            cj5=- idat(j2,k3)
            if(dat(i2,k3).lt.dat(j2,k3)) cj5=idat(i2,k3)
            cj6=- idat(i2,k2)
            if(dat(j2,k2).lt.dat(i2,k2)) cj6=idat(j2,k2)
            cj7=- idat(i2,k3)
            if(dat(j2,k3).lt.dat(i2,k3)) cj7=idat(j2,k3)
            gtsum=cj1*(2.0*cj2 - cj3*(cj4*cj5+cj6*cj7) )

            cj1=- idat(j2,k1)
            if(dat(i1,k1).lt.dat(j2,k1)) cj1=idat(i1,k1)
            cj2=- idat(j2,k2)
            if(dat(i1,k2).lt.dat(j2,k2)) cj2=idat(i1,k2)
            cj3=- idat(j2,k3)
            if(dat(i1,k3).lt.dat(j2,k3)) cj3=idat(i1,k3)
            cj4=- idat(j1,k2)
            if(dat(i2,k2).lt.dat(j1,k2)) cj4=idat(i2,k2)
            cj5=- idat(j1,k3)
            if(dat(i2,k3).lt.dat(j1,k3)) cj5=idat(i2,k3)
            cj6=- idat(i2,k2)
            if(dat(j1,k2).lt.dat(i2,k2)) cj6=idat(j1,k2)
            cj7=- idat(i2,k3)
            if(dat(j1,k3).lt.dat(i2,k3)) cj7=idat(j1,k3)
            gtsum=gtsum+cj1*(2.0*cj2 - cj3*(cj4*cj5+cj6*cj7) )

            cj1=- idat(i2,k1)
            if(dat(i1,k1).lt.dat(i2,k1)) cj1=idat(i1,k1)
            cj2=- idat(i2,k2)
            if(dat(i1,k2).lt.dat(i2,k2)) cj2=idat(i1,k2)
            cj3=- idat(i2,k3)
            if(dat(i1,k3).lt.dat(i2,k3)) cj3=idat(i1,k3)
            cj4=- idat(j1,k2)
            if(dat(j2,k2).lt.dat(j1,k2)) cj4=idat(j2,k2)
            cj5=- idat(j1,k3)
            if(dat(j2,k3).lt.dat(j1,k3)) cj5=idat(j2,k3)
            cj6=- idat(j2,k2)
            if(dat(j1,k2).lt.dat(j2,k2)) cj6=idat(j1,k2)
            cj7=- idat(j2,k3)
            if(dat(j1,k3).lt.dat(j2,k3)) cj7=idat(j1,k3)
            gtsum=gtsum+cj1*(2.0*cj2 - cj3*(cj4*cj5+cj6*cj7) )

            cj1=- idat(i1,k1)
            if(dat(j1,k1).lt.dat(i1,k1)) cj1=idat(j1,k1)
            cj2=- idat(i1,k2)
            if(dat(j1,k2).lt.dat(i1,k2)) cj2=idat(j1,k2)
            cj3=- idat(i1,k3)
            if(dat(j1,k3).lt.dat(i1,k3)) cj3=idat(j1,k3)
            cj4=- idat(j2,k2)
            if(dat(i2,k2).lt.dat(j2,k2)) cj4=idat(i2,k2)
            cj5=- idat(j2,k3)
            if(dat(i2,k3).lt.dat(j2,k3)) cj5=idat(i2,k3)
            cj6=- idat(i2,k2)
            if(dat(j2,k2).lt.dat(i2,k2)) cj6=idat(j2,k2)
            cj7=- idat(i2,k3)
            if(dat(j2,k3).lt.dat(i2,k3)) cj7=idat(j2,k3)
            gtsum=gtsum+cj1*(2.0*cj2 - cj3*(cj4*cj5+cj6*cj7) )

            cj1=- idat(i2,k1)
            if(dat(j1,k1).lt.dat(i2,k1)) cj1=idat(j1,k1)
            cj2=- idat(i2,k2)
            if(dat(j1,k2).lt.dat(i2,k2)) cj2=idat(j1,k2)
            cj3=- idat(i2,k3)
            if(dat(j1,k3).lt.dat(i2,k3)) cj3=idat(j1,k3)
            cj4=- idat(j2,k2)
            if(dat(i1,k2).lt.dat(j2,k2)) cj4=idat(i1,k2)
            cj5=- idat(j2,k3)
            if(dat(i1,k3).lt.dat(j2,k3)) cj5=idat(i1,k3)
            cj6=- idat(i1,k2)
            if(dat(j2,k2).lt.dat(i1,k2)) cj6=idat(j2,k2)
            cj7=- idat(i1,k3)
            if(dat(j2,k3).lt.dat(i1,k3)) cj7=idat(j2,k3)
            gtsum=gtsum+cj1*(2.0*cj2 - cj3*(cj4*cj5+cj6*cj7) )

            cj1=- idat(j2,k1)
            if(dat(j1,k1).lt.dat(j2,k1)) cj1=idat(j1,k1)
            cj2=- idat(j2,k2)
            if(dat(j1,k2).lt.dat(j2,k2)) cj2=idat(j1,k2)
            cj3=- idat(j2,k3)
            if(dat(j1,k3).lt.dat(j2,k3)) cj3=idat(j1,k3)
            cj4=- idat(i2,k2)
            if(dat(i1,k2).lt.dat(i2,k2)) cj4=idat(i1,k2)
            cj5=- idat(i2,k3)
            if(dat(i1,k3).lt.dat(i2,k3)) cj5=idat(i1,k3)
            cj6=- idat(i1,k2)
            if(dat(i2,k2).lt.dat(i1,k2)) cj6=idat(i2,k2)
            cj7=- idat(i1,k3)
            if(dat(i2,k3).lt.dat(i1,k3)) cj7=idat(i2,k3)
            gtsum=gtsum+cj1*(2.0*cj2 - cj3*(cj4*cj5+cj6*cj7) )

            cj1=- idat(i1,k1)
            if(dat(i2,k1).lt.dat(i1,k1)) cj1=idat(i2,k1)
            cj2=- idat(i1,k2)
            if(dat(i2,k2).lt.dat(i1,k2)) cj2=idat(i2,k2)
            cj3=- idat(i1,k3)
            if(dat(i2,k3).lt.dat(i1,k3)) cj3=idat(i2,k3)
            cj4=- idat(j2,k2)
            if(dat(j1,k2).lt.dat(j2,k2)) cj4=idat(j1,k2)
            cj5=- idat(j2,k3)
            if(dat(j1,k3).lt.dat(j2,k3)) cj5=idat(j1,k3)
            cj6=- idat(j1,k2)
            if(dat(j2,k2).lt.dat(j1,k2)) cj6=idat(j2,k2)
            cj7=- idat(j1,k3)
            if(dat(j2,k3).lt.dat(j1,k3)) cj7=idat(j2,k3)
            gtsum=gtsum+cj1*(2.0*cj2 - cj3*(cj4*cj5+cj6*cj7) )

            cj1=- idat(j1,k1)
            if(dat(i2,k1).lt.dat(j1,k1)) cj1=idat(i2,k1)
            cj2=- idat(j1,k2)
            if(dat(i2,k2).lt.dat(j1,k2)) cj2=idat(i2,k2)
            cj3=- idat(j1,k3)
            if(dat(i2,k3).lt.dat(j1,k3)) cj3=idat(i2,k3)
            cj4=- idat(j2,k2)
            if(dat(i1,k2).lt.dat(j2,k2)) cj4=idat(i1,k2)
            cj5=- idat(j2,k3)
            if(dat(i1,k3).lt.dat(j2,k3)) cj5=idat(i1,k3)
            cj6=- idat(i1,k2)
            if(dat(j2,k2).lt.dat(i1,k2)) cj6=idat(j2,k2)
            cj7=- idat(i1,k3)
            if(dat(j2,k3).lt.dat(i1,k3)) cj7=idat(j2,k3)
            gtsum=gtsum+cj1*(2.0*cj2 - cj3*(cj4*cj5+cj6*cj7) )

            cj1=- idat(j2,k1)
            if(dat(i2,k1).lt.dat(j2,k1)) cj1=idat(i2,k1)
            cj2=- idat(j2,k2)
            if(dat(i2,k2).lt.dat(j2,k2)) cj2=idat(i2,k2)
            cj3=- idat(j2,k3)
            if(dat(i2,k3).lt.dat(j2,k3)) cj3=idat(i2,k3)
            cj4=- idat(j1,k2)
            if(dat(i1,k2).lt.dat(j1,k2)) cj4=idat(i1,k2)
            cj5=- idat(j1,k3)
            if(dat(i1,k3).lt.dat(j1,k3)) cj5=idat(i1,k3)
            cj6=- idat(i1,k2)
            if(dat(j1,k2).lt.dat(i1,k2)) cj6=idat(j1,k2)
            cj7=- idat(i1,k3)
            if(dat(j1,k3).lt.dat(i1,k3)) cj7=idat(j1,k3)
            gtsum=gtsum+cj1*(2.0*cj2 - cj3*(cj4*cj5+cj6*cj7) )

            cj1=- idat(i1,k1)
            if(dat(j2,k1).lt.dat(i1,k1)) cj1=idat(j2,k1)
            cj2=- idat(i1,k2)
            if(dat(j2,k2).lt.dat(i1,k2)) cj2=idat(j2,k2)
            cj3=- idat(i1,k3)
            if(dat(j2,k3).lt.dat(i1,k3)) cj3=idat(j2,k3)
            cj4=- idat(i2,k2)
            if(dat(j1,k2).lt.dat(i2,k2)) cj4=idat(j1,k2)
            cj5=- idat(i2,k3)
            if(dat(j1,k3).lt.dat(i2,k3)) cj5=idat(j1,k3)
            cj6=- idat(j1,k2)
            if(dat(i2,k2).lt.dat(j1,k2)) cj6=idat(i2,k2)
            cj7=- idat(j1,k3)
            if(dat(i2,k3).lt.dat(j1,k3)) cj7=idat(i2,k3)
            gtsum=gtsum+cj1*(2.0*cj2 - cj3*(cj4*cj5+cj6*cj7) )

            cj1=- idat(j1,k1)
            if(dat(j2,k1).lt.dat(j1,k1)) cj1=idat(j2,k1)
            cj2=- idat(j1,k2)
            if(dat(j2,k2).lt.dat(j1,k2)) cj2=idat(j2,k2)
            cj3=- idat(j1,k3)
            if(dat(j2,k3).lt.dat(j1,k3)) cj3=idat(j2,k3)
            cj4=- idat(i1,k2)
            if(dat(i2,k2).lt.dat(i1,k2)) cj4=idat(i2,k2)
            cj5=- idat(i1,k3)
            if(dat(i2,k3).lt.dat(i1,k3)) cj5=idat(i2,k3)
            cj6=- idat(i2,k2)
            if(dat(i1,k2).lt.dat(i2,k2)) cj6=idat(i1,k2)
            cj7=- idat(i2,k3)
            if(dat(i1,k3).lt.dat(i2,k3)) cj7=idat(i1,k3)
            gtsum=gtsum+cj1*(2.0*cj2 - cj3*(cj4*cj5+cj6*cj7) )

            cj1=- idat(i2,k1)
            if(dat(j2,k1).lt.dat(i2,k1)) cj1=idat(j2,k1)
            cj2=- idat(i2,k2)
            if(dat(j2,k2).lt.dat(i2,k2)) cj2=idat(j2,k2)
            cj3=- idat(i2,k3)
            if(dat(j2,k3).lt.dat(i2,k3)) cj3=idat(j2,k3)
            cj4=- idat(j1,k2)
            if(dat(i1,k2).lt.dat(j1,k2)) cj4=idat(i1,k2)
            cj5=- idat(j1,k3)
            if(dat(i1,k3).lt.dat(j1,k3)) cj5=idat(i1,k3)
            cj6=- idat(i1,k2)
            if(dat(j1,k2).lt.dat(i1,k2)) cj6=idat(j1,k2)
            cj7=- idat(i1,k3)
            if(dat(j1,k3).lt.dat(i1,k3)) cj7=idat(j1,k3)
            gtsum=gtsum+cj1*(2.0*cj2 - cj3*(cj4*cj5+cj6*cj7) )

       aasum(i1)=aasum(i1)+1./24.*gtsum !ADD SUMMATION OVER PERMUTATIONS
13          continue              !
12        continue                !
11      continue                  !
      ave = ave + c2*aasum(i1)
10    continue
      ave=ave/float(ntot)
      do 20 i=1,ntot
      asum=asum+(c2*aasum(i)-ave)**2
 20   continue
      an=asum*c1
      return
      end

C------------- TAU -------------------------------------------
C------- COMPUTES KENDALLS TAU -------------------------------

      function tau(k,l)
      common ntot
      ac=2./(float(ntot)*(float(ntot)-1))
      sum=0.0
      do 11 j=1,ntot
        do 12 i=1,ntot
          if (i.ge.j) goto 11
          sum=sum+h(k,l,i,j)
12      continue
11    continue
      tau=sum*ac
      return
      end

C-------------- H --------------------------------------------
C------- COMPUTE VALUE FOR H (SEE FORMULA) --------------------

      function h(k,l,i,j)
      common /data/ dat(500,3),idat(500,3)
      cj1=-idat(j,k)
      if(dat(i,k).lt.dat(j,k)) cj1=idat(i,k)
      cj2=-idat(j,l)
      if(dat(i,l).lt.dat(j,l)) cj2=idat(i,l)
      h=cj1*cj2
      return
      end
