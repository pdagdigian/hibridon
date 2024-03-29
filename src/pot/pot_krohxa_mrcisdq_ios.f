*system:  routine for OH(A/X)-Kr system
*         with MRCISD+Q/AVTZDK J. Klos's PESs
*
*  this is a modification of pot_krohxa_mrcisdq_basgpi1.f
*  pot routine
*
*  in this version, a 3-state calculation (sigma A', pi A', pi A")
*  will be set up.  this is appropriate to an IOS treatment of the
*  dynamics, as described by mha and g.c.corey, jcp 84, 100 (1986)
*
*  divided V1 PES by sqrt(2).  see mha and corey, jcp 84, 100 (1986)
*
*  routine set up for one vibrational level in each electronic state
*
*  written by P.j.dagdigian
*  current revision date;  17-apr-2012
*   
      subroutine driver
      implicit double precision (a-h,o-z)
      character *48 potnam
      logical csflag, ljunk, ihomo, lljunk
      include "common/parbas"
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(4)
      common /covpot/ numvib,ivibpi(5)
      common /cosysr/ isrcod, junkr, e1, e2, thta
      potnam='OH(A/X)-Kr Quenching Klos MRCISDQ PESs'
      print *, potnam
*  consider only v=0 and 1 vib levels of 2pi state
      numvib = 1
      ivibpi(1)=0
      nterm=1
      write (6,89) numvib, (ivibpi(i),i=1,numvpi)
89    format(' Number of 2pi vibrational levels =',i3)
      print *
1     print *, ' theta (angle), r (bohr) '
      read (5, *, end=93) thta, r
      if (r .lt. 0.d0) go to 93
      call pot(vv0,r)
      write (6, 100) vvl 
* 219474.6
100   format(' vA''',3(pe16.8),
     :    '    vA" ',1pe16.8)
      goto 1
93    continue
99    end
* --------------------------------------------------------------------------
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine syusr (irpot, readpt, iread)
*  subroutine to read in system dependent parameters for sigma-sigma
*   + atom scattering (collisional transfer between two sigma states)
*  we consider here one vibrational level in each sigma state
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 16-mar-2012 by p.dagdigian
*  -----------------------------------------------------------------------
*  variables in common /cosysr/
*    isrcod:  number of real parameters
*    e1:       energy of 1st state in cm-1
*    e2:       energy of 2nd state in cm-1
*    thta:     angle between R and r, in degrees
*  variable in common /cosysi/
*    nscode:   total number of system dependent parameters
*              nscode = isicod + isrcod + 3
*    isicod:   number of integer system dependent parameters
*    nterm:    number of different diabatic potentials, set to 1 in this model
*  variable in common /cosys/
*    scod:    character*8 array of dimension nscode, which contains names
*             of all system dependent parameters.  note that the ordering
*             of the variable names in scod must correspond to the ordering
*             of the variable names in cosysi followed by the ordering of
*             variable names in cosysr followed by lammin, lammax, and mproj
*
*  subroutines called: loapot(iunit,filnam)
*  -----------------------------------------------------------------------
*
*  variables in common/cobspt/ must be set in loapot!!
*
      implicit double precision (a-h,o-z)
      include "common/parsys"
      logical readpt, existf
      character*8 scod, rcod
      character*8 char
      character*(*) fname
      character*1 dot
      character*60 filnam, line, potfil, filnm1
      include "common/parbas"
      common /cosys/ scod(maxpar)
      common /cosyr/ rcod(maxpar)
      common /cosysi/ nscode, isicod, ispar(maxpar)
      common /cosysr/ isrcod, junkr, e1, e2, thta
      logical airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar
      common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo, lpar(18)
      common /coskip/ nskip,iskip
      common /conlam/ nlam
      save potfil
      include "common/comdot"
      irpot = 1
*     number and names of system dependent parameters
*  first all the system dependent integer variables
*  in the same order as in the common block /cosysi/
*  variable names must be in uppercase of maximum length 6 characters
      scod(1) = 'NTERM'
      nterm = 1
      ispar(1) = nterm
      isicod = 1
*  then all the system dependent real variables
*  in the same order as in the common block /cosysr/
      scod(isicod+1) = 'E1'
      scod(isicod+2) = 'E2'
      scod(isicod+3) = 'THETA'
      isrcod = 3
*  set default values
      potfil = ' '
      if (iread .eq. 0) then
        niout=2
        indout(1)=1
        indout(2)=2
*  parameters for pes's
        mproj(1) = 0
        lammin(1) = 1
        lammax(1) = 4
      end if
      if(iread.eq.0) return
*  read e1 and e2
      read (8, *, err=888) e1, e2
*  read angle between R and r
      read (8, *, err=888) thta
      nscode = isicod + isrcod + 3
*  read name of file containing potential parameters
      if(.not.readpt.or.iread.eq.0) then
        call loapot(1,' ')
        return
      end if
      read (8, 285, end=888) line
      potfil=line
285   format (a)
      goto 1000
* here if read error occurs
888   write(6,900)
900   format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
*
      entry ptrusr (fname,readpt)
      line = fname
      readpt = .true.
1000  if (readpt) then
        l=1
        call parse(line,l,filnam,lc)
        if(lc.eq.0) then
          write(6,1020)
1020      format(' FILENAME MISSING FOR POTENTIAL INPUT')
        end if
        j=index(filnam(1:lc),dot)
        if(j.eq.0) then
          call gennam(potfil,filnam,0,'BIN',lc)
          filnam = potfil
        end if
        potfil=filnam
        filnm1 = 'potdata/'//filnam
        inquire(file=filnm1,exist=existf)
        if(.not.existf) then
          write(6,1025) filnam(1:lc)
1025      format(' FILE ',(a),' NOT FOUND')
          return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
      end if
* check for consistency
      close (8)
      irpot=1
      return
*
      entry savusr (readpt)
*  save parameters
      write (8, 201) e1, e2,
     :   ' e1, e2'
      write (8, 203) thta,
     :   ' thta'
201   format(2g14.6,t55,a)
203   format(g14.6,t55,a)
      write (8, 300) potfil
300   format(a)
      return
      end
* --------------------------------------------------------------------------
      subroutine bausr (j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential
*  for collisional coupling of 2sigma and 2pi state at a particular
*  angle thta between R and r, for ios calculation
*  author:  paul dagdigian
*  current revision date:  18-apr-2012
*  variables in call list:
*    j:        on return contains rotational quantum numbers for each
*              channel (dummy)
*    l:        on return contains orbital angular momentum for each
*              channel (zero)
*    is:       on return contains the state label for each channel
*              equals 1 for first 1Sigma+ state
*              equals 2 for second 1Sigma+state
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level (dummy)
*    ehold:    on return contains energy in hartrees of each level
*    ishold:   on return contains symmetry index of each rotational level
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    sc1,sc2:  scratch vectors of length at least nmax
*    sc3,sc4:  scratch vectors of length at least nmax
*              these scratch vectors are not used here
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              note that this is ignored in molecule-surface collisions!!!
*    jtot:     total angular momentum
*    flaghf:   if .true., then system with half-integer spin
*              if .false., then system with integer spin (this is the case
*              here)
*    flagsu:   if .true., then molecule-surface collisons
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*              NOTE:  csflag not used in this basis routine
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true. , then homonuclear molecule
*              if .false., then heteronuclear molecule
*              if the molecule is homonuclear (ihomo = .true.), the
*              rotational levels included go from jmin to jmax in steps
*              of 2 and only even lambda terms in the anisotropy are
*              included. Not used here.
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              only those channels are included for which
*                  (-1)**(parity+l-jtot)=jlpar
*              where parity is defined for each electronic, rotationnal
*              vibrationnal state, Here both sigma 2 surfaces have the
*              same symmetry and no rotation is taken into account.
*              We just define jlpar=1 for all channels involved in
*              calculation.
*              in cs calculation jlpar is set equal to 1 in calling
*              program
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*    note!!!   we mst have here flaghf = .false.
*  variables in common /cosysr/
*    isrcod:  number of real parameters
*    e1:        energy of 1st state in cm-1
*    e2:        energy of 2nd state in cm-1
*    thta:      angle between R and r in degrees
*  variable in common /cosysi/
*    nscode:   total number of system dependent parameters
*              nscode = isicod + isrcod +3
*  variables in common block /cobspt/
*    lammin:   array containing minimum value of lambda for each term
*    lammax:   array containing maximum value of lambda for each term
*    mproj:    array containing the order of the reduced rotation matrix
*              elements for each term.  lammin can not be less than mproj.
*              for homonuclear molecules, the allowed values of lambda for
*              each term range from lammin to lammax in steps of 2
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coint/
*    eint:      array containing channel energies (in hartree)
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units
*               (mass of electron = 1)
*  variable in common block /conlam/
*    nlam:      the number of angular coupling terms actually used
*    nlammx:    the maximum number of angular coupling terms
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               only nonzero elements are stored
*  variable in common block /coiv2/
*   iv2:        matrix address of v2 matrix for each non-zero element
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ihomo, flaghf, csflag, clist, flagsu, bastst
      include "common/parbas"
      include "common/parbasl"
      include "common/parsys"
      common /coipar/ iiipar(9), iprint
      common /cosysi/ nscode, isicod, nterm
      common /cosysr/ isrcod, junkr, e1, e2, thta
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/  nlam, nlammx,lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      dimension j(1), l(1), jhold(1), ehold(1),
     :          sc1(1), sc2(1), sc3(1),
     :          sc4(1), ishold(1), is(1)
*   econv is conversion factor from cm-1 to hartrees
*   xmconv is conversion factor from amu to atomic units
      tzero = 1.d-9
*  check for consistency in the value of flaghf
      if (flaghf) then
        write (6, 1)
        write (9, 1)
1       format (' *** FLAGHF MUST BE FALSE; ABORT ***' )
        if (bastst) then
          return
        else
          call exit
        end if
      end if
*  do not set up basis routine for gas-surface scattering
      if (flagsu) then
        write (6, 8)
        write (9, 8)
8      format
     :   ('  *** CSFLAG = .FALSE. FOR SURFACE CALCULATION; ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
*      
      if (clist) then
        write (6,20) rmu * xmconv, ered * econv, jtot, jlpar
        write (9,20) rmu * xmconv, ered * econv, jtot, jlpar
20      format(/,' **  IOS 2SIGMA-2PI TRANSITIONS  ** RMU=', f9.4,
     :           '   E=',f10.3,'   JTOT=', i5, 2x,' JLPAR=',i2)
        if (.not. flagsu) write (9,30) rcut
30      format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :          f6.2)
*  print out spectroscopic parameters for the states
        write(6,31) ' State    Energy'
        write(9,31) ' State    Energy'
31      format(/a,8x,a)
        write(6,35) e1, e2
        write(9,35) e1, e2
35      format (1x,'1    ',f12.6/
     :    1x,'2    ',f12.6)
      end if
*
*  list the 3 states
*  2sigma state
      n=1
      j(n) = 0
      cent(n) = jtot * (jtot + 1)
      is(n) = 1
      l(n) = jtot
      eint(n) = e1 / econv
*  2pi(A') state
      n = n + 1
      j(n) = 0
      cent(n) = jtot * (jtot + 1)
      is(n) = 2
      l(n) = jtot
      eint(n) = e2 / econv
*  2pi(A") state
      n = n + 1
      j(n) = 1
      cent(n) = jtot * (jtot + 1)
      is(n) = 2
      l(n) = jtot
      eint(n) = e2 / econv
*       
130   if (n .gt. nmax) then
        write (9, 140) n, nmax
        write (6, 140) n, nmax
140     format(/' *** NCHANNELS=', i3,' .GT. MAX DIMENSION OF',
     :       i3,'; ABORT')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
*  return if no channels
      if (n .eq. 0) then
        write (6, 1140)
        write (9, 1140)
1140    format(/' *** NCHANNELS = 0; ABORT')
        return
      end if
*
*  we will not do rcut test here
*
*  return if no channels
      if (n .eq. 0) then
        write (6, 1140)
        write (9, 1140)
        return
      end if
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
      ntop = 3
      nlevel = 0
*  form list of all energetically distinct rotational levels included
*  calculations and their energies
      do 200 i=1,n
        nlevel = nlevel + 1
        ehold(nlevel) = eint(i)
        jhold(nlevel) = j(i)
        ishold(nlevel) = is(i)
200   continue
*  now sort this list to put closed levels at end
*  also determine number of levels which are open
      nlevop = 0
      if (nlevel .gt. 1) then
        do 80  i = 1, nlevel - 1
          if (ehold(i) .le. ered) then
             nlevop = nlevop + 1
          else
            do 75 ii = i + 1, nlevel
              if (ehold(ii) .le. ered) then
                nlevop = nlevop + 1
                ikeep = jhold(i)
                jhold(i) = jhold(ii)
                jhold(ii) = ikeep
                ikeep = ishold(i)
                ishold(i) = ishold(ii)
                ishold(ii) = ikeep
                ekeep = ehold(i)
                ehold(i) = ehold(ii)
                ehold(ii) = ekeep
                go to 80
              end if
75          continue
          end if
80      continue
      else
* here for only one level
        if (ehold(1) .le. ered) then
          nlevop=1
        else
          nlevop=0
        endif
      endif
      if (nlevop .le. 0) then
        write (9, 85)
        write (6, 85)
85      format('*** NO OPEN LEVELS IN BAUSR; ABORT')
        if (bastst) return
        call exit
      endif
      if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
      ntop = max(n, nlevop)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1
*  now list channels if requested
      if (clist) then
        write (6, 255)
        write (9, 255)
255     format(/'   N   State   J    L    EINT(CM-1)')
        do 265  i = 1, n
          write (6, 260) i, is(i), j(i), l(i), eint(i) * econv
          write (9, 260) i, is(i), j(i), l(i), eint(i) * econv
260       format (i4, i6, 1x, 2i5, f13.3)
265     continue
        write (6, 256)
256     format(/' OPEN LEVELS:'//
     1          '   N   State   J      EINT(CM-1)')
        do 266  i = 1, nlevop
          write (6, 261) i, ishold(i), jhold(i),  ehold(i) * econv
261       format (i4, i6, i6, f15.3)
266     continue
      end if
*
*  now calculate coupling matrix elements
      if (bastst .and. iprint.gt.1) then
        write (6, 280)
        write (9, 280)
280     format (/' ILAM   LAMBDA   MU     I   ICOL  IROW',
     :           '  IV2      VEE')
      end if
*
* i counts v2 elements
* inum counts v2 elements for given lambda
* ilam counts numver of v2 matrices
* ij is address of given v2 element in present v2 matrix
      nlam = 0
      i = 0
      inum = 0
      nlami = lammax(1) - lammin(1) + 1
      do 600 it=1,nlami
        nlam = nlam + 1
        lb = lammin(1) + (it - 1)
        mu = mproj(1)
        inum = 0
        do 450 icol = 1, n
        do 450 irow = icol, n
          ij = ntop * (icol - 1) + irow
          vee = zero
          goto (455, 460, 465, 470), it
*  vsig
455       if (it.eq.1 .and. irow.eq.1 .and.
     :        icol.eq.1) then
            vee = 1.d0
          end if
          go to 440
*  vpi(A')
460       if (it.eq.2 .and. irow.eq.2 .and.
     :        icol.eq.2) then
            vee = 1.d0
          end if
          go to 440
*  v1
465       if (it.eq.3 .and. ((irow.eq.1 .and.
     :        icol.eq.2) .or. (irow.eq.2 .and.
     :        icol.eq.1)) ) then
            vee = 1.d0
          end if
          go to 440
*  vpi(A")
470       if (it.eq.4 .and. irow.eq.3 .and.
     :        icol.eq.3) then
            vee = 1.d0
          end if
440       if (abs(vee).lt.tzero) goto 450
          i = i + 1
          inum = inum + 1
          if (i.gt.nv2max) goto 450
          v2(i) = vee
          iv2(i) = ij
          if (.not.bastst .or. iprint.le.1) goto 450
          write (6,495) nlam, lb, mu, i, icol, irow, iv2(i), vee
          write (9,495) nlam, lb, mu, i, icol, irow, iv2(i), vee
495       format (i4, 3i7, 2i6, i6, g17.8)
450     continue
        if (i.le.nv2max) lamnum(nlam) = inum
        if (bastst) then
            write (6, 360) it, nlam, lb, mu, lamnum(nlam)
            write (9, 360) it, nlam, lb, mu, lamnum(nlam)
360         format ('ITERM=',i3,' ILAM=', i3,' LAMBDA=',i3,
     :        ' MU=',i3,' LAMNUM(ILAM)=',i8)
          end if
          lamsum = lamsum + lamnum(il)
500     continue
600   continue
      if (nlam.gt.nlammx) then
        write(6,610) nlam,nlammx
        write(9,610) nlam,nlammx
610     format (' *** NLAM = ',i6,' .GT. NLAMMX=',i6,'; ABORT ***')
        call exit
      end if
      if (i.gt.nv2max) then
        write (6, 620) i, nv2max
        write (9, 620) i, nv2max
620     format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6,
     :           ' .GT. NV2MAX=',i6,'; ABORT ***')
        call exit
      end if
      if (clist) then
        write (6, 630) i
        write (9, 630) i
630     format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :           i9)
      end if
      return
1000  write (6, 1010) n, nmax
      write (9, 1010) n, nmax
1010  format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF',
     :         i4,' ABORT ***')
      call exit
      return
      end
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      common /covpot/ numvib,ivibpi(5)
      common /cosysr/ isrcod, junkr, thta
       potnam='OH(A/X)-Kr Quenching Klos MRCISDQ PESs'
      ibasty=99
      numvib=1
*  consider only v=0 and 1 vib level of X state
      ivibpi(1)=0
      ivibpi(2)=1
      nterm=1
      lammin(1)=1
      lammax(1)=4
      mproj(1)=0
      return
      end
* --------------------------------------------------------------------------
      subroutine pot (vv0, r)
* ----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients 
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*    vv0      dummy term here
*  variable in common block /covvl/
*    vvl:     vector to store r-dependence of each term
*             in potential expansion
*    vvl(1)   value of vsigma for the given angle
*    vvl(2)   value of vpi(A')
*    vvl(3)   value of v1
*    vvl(4)   value of vpi(A")
*
* authors:  paul dagdigian and jacek klos
* current revision:  16-apr-2012 by p. dagdigian
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      dimension FCF(2) ! stores Franck-Condon factors
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(4)
      common /covpot/ numvib,ivibpi(5)
      common /cosysr/ isrcod, junkr, e1, e2, thta

* hyperbolic tangent scaling factor
      data alph /1.2d0/
      data rmax /13d0/
      aa=zero
* Here are Franck-Condon factors (as squares of the <Av=0|Xv=v'>
* overlap integrals) between transitions from vibrational
* state v=0 of A state
* to vibrational states v' of the X state:
* FCF_0_0 = 0.917664288481494
* FCF_0_1 = 0.0800257653502253
* FCF_0_2 =0.00227074305148058
* FCF_0_3 =3.83580435449925e-05
*Obtained from MRCISD+Q/AVQZ diatomic X and A potentials.
*J. Klos 2011 November
       FCF(1)=0.917664288481494D0
       FCF(2)=0.0800257653502253D0
*
      vsig = VCI_H22(r,thta)
      vpiap = VCI_H11(r,thta)
      vpiapp = VCI_XAbis(r,thta)
*   divide V1 by sqrt(2)
      v1 = VCI_H12(r,thta) / sqrt(2.d0)
      if (r .gt. rmax) then
        damp=-half*(tanh(3.d0*(r-rmax))-one)
        v1=v1*damp
      endif
      vvl(1) = vsig
      vvl(2) = vpiap
      vvl(3) = v1 * FCF(1)
      vvl(4) = vpiapp
* convert to hartree
      econv=1./219474.6
      call dscal(4,econv,vvl,1)
      vv0 = 0.d0
      return
      end

C******************************************************

      FUNCTION PLGNDR(L,M,x)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)pause
     *'bad arguments in plgndr'
      pmm=1.
      if(m.gt.0) then
        somx2=dsqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
         do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END
C******************************************************
C    BELOW Kr-OH(X-A) MRCISD+Q/AVTZDK PES ROUTINES
C******************************************************

      FUNCTION VCI_H11(R,theta)
C*********************************
C System: Kr-OH(X2Pi - A2Sigma) MIXING
C Method: MRCISD+Q
C Basis: AVTZ-DK
C with no bond functions
C Size Consistency Corrected
C PES: H11 CI DIABAT
C Variables:
C R in Bohr
C Theta n Degrees
C Output in cm-1
C Theta=0 for Kr---HO geom.
C********************************
C Authors:
C Dr Jacek Klos
C Department of Chemistry and Biochemistry
C University of Maryland
C College Park, MD 20742
C email:jklos@umd.edu
C**********************************
      implicit double precision(a-h, o-z)
      dimension V0(11)
      dimension T(11)
      pi=dacos(-1.d0)
      conv=1.D0 ! IN CM-1
      V0(1)=H11CI_VL0_0(r)*conv    
      V0(2)=H11CI_VL0_1(r)*conv
      V0(3)=H11CI_VL0_2(r)*conv
      V0(4)=H11CI_VL0_3(r)*conv
      V0(5)=H11CI_VL0_4(r)*conv
      V0(6)=H11CI_VL0_5(r)*conv
      V0(7)=H11CI_VL0_6(r)*conv
      V0(8)=H11CI_VL0_7(r)*conv
      V0(9)=H11CI_VL0_8(r)*conv
      V0(10)=H11CI_VL0_9(r)*conv
      V0(11)=H11CI_VL0_10(r)*conv
       do j=1,11
       T(j)=PLGNDR((j-1),0,DCOS(theta*pi/180.D0))
       enddo
       s=0.D0
       do i=1,11
       s=s+V0(i)*T(i)
       enddo
       VCI_H11=s
       return
       end

C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H11 DIABAT WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H11CI_VL0_0(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.754022545523574257D+10,
     *   -0.188885854017529640D+11,
     *    0.304973445473534317D+11,
     *   -0.235179915963814392D+11,
     *   -0.137423787556856680D+10,
     *    0.231768964888736038D+11,
     *   -0.206715019149821739D+11,
     *    0.873723032954555893D+10,
     *   -0.113282242877030106D+11,
     *    0.123893124458545895D+11,
     *   -0.962926769739731598D+10,
     *    0.497369998465775108D+10,
     *   -0.147157345009616566D+10,
     *    0.621664517864005327D+09,
     *   -0.127765052327549648D+10,
     *    0.263975303284722328D+09,
     *   -0.987744340910186768D+09,
     *   -0.180437241938340664D+09,
     *   -0.660646635749313354D+09,
     *   -0.117160994504904866D+09,
     *   -0.126687371394032300D+09,
     *   -0.256609981959662437D+09,
     *    0.221323955851361275D+09,
     *   -0.640724794580335617D+08,
     *    0.105095260322260857D+09,
     *    0.342697512039466143D+09,
     *    0.141603144376274347D+08,
     *    0.395963099439535618D+09,
     *    0.107174495115200520D+09,
     *    0.427634794234447598D+09,
     *    0.157776395526118279D+08,
     *    0.197311467338355064D+09,
     *    0.188876198678314209D+09,
     *   -0.248757443464279175D+08,
     *    0.126676712649784088D+09,
     *   -0.196583943733751297D+09,
     *    0.198043617345201492D+09,
     *   -0.200652484682610512D+09,
     *    0.159089512843903542D+09,
     *   -0.815070216723299026D+08,
     *    0.388953236640494347D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H11CI_VL0_0=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H11 DIABAT WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H11CI_VL0_10(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.451182067686770248D+10,
     *   -0.109648030854760399D+11,
     *    0.921162244282719040D+10,
     *   -0.412180102918363476D+10,
     *    0.204403672116384354D+11,
     *   -0.517274375167659683D+11,
     *    0.437261909249506760D+11,
     *   -0.122740233026348629D+11,
     *    0.868610618109011650D+09,
     *   -0.665688225320239544D+10,
     *    0.426728080571465302D+11,
     *   -0.685035176155561371D+11,
     *    0.359280709915393906D+11,
     *   -0.465374081436768150D+10,
     *    0.237660444102289581D+10,
     *   -0.128795780794557142D+10,
     *    0.533025176829717159D+09,
     *   -0.129146449620189369D+08,
     *   -0.169976350610788703D+09,
     *    0.201398572551484674D+09,
     *   -0.193797439807785153D+09,
     *    0.133287330372264504D+09,
     *   -0.308282254593070745D+08,
     *   -0.397678056801165342D+08,
     *   -0.843127037745648026D+08,
     *    0.420922042865592957D+09,
     *   -0.508288241198408008D+09,
     *    0.155571066384999752D+09,
     *    0.204824476882104039D+09,
     *    0.569843037957491279D+08,
     *   -0.754176072962012172D+09,
     *    0.692349981286295414D+09,
     *   -0.215715378544527054D+09,
     *    0.155826671810019970D+09,
     *   -0.737824646725330353D+08,
     *   -0.782784677470560074D+08,
     *    0.987405133275794983D+07,
     *    0.280064688170995712D+08,
     *    0.192358175810528398D+09,
     *   -0.332068277006895840D+09,
     *    0.170728747444962233D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H11CI_VL0_10=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H11 DIABAT WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H11CI_VL0_1(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.830388126960833931D+10,
     *   -0.729687168057940960D+10,
     *   -0.424782528216062117D+10,
     *   -0.281748669787919712D+10,
     *    0.730639575179757385D+11,
     *   -0.131805420495141739D+12,
     *    0.766388544478160858D+11,
     *   -0.123428199885251770D+11,
     *   -0.579758356266753197D+09,
     *   -0.379859246768474579D+09,
     *    0.314986540039419632D+11,
     *   -0.555753354902990875D+11,
     *    0.283607460511757813D+11,
     *   -0.383970175676017141D+10,
     *    0.283653155040002060D+10,
     *   -0.231538818815261745D+10,
     *    0.539290627190757036D+09,
     *   -0.156770476658982217D+09,
     *   -0.101799891025373971D+10,
     *    0.600828648291476011D+09,
     *   -0.317549775457067847D+09,
     *   -0.616394994305056691D+09,
     *    0.632501688545850873D+09,
     *   -0.651348605134928703D+09,
     *    0.451647129016574860D+09,
     *   -0.664319310415291786D+08,
     *    0.322338052501032352D+09,
     *   -0.227864656490840912D+08,
     *    0.131457999389681578D+09,
     *    0.264092194600946426D+09,
     *   -0.125748999939052343D+09,
     *    0.244312075558607578D+09,
     *    0.254337375979538918D+09,
     *   -0.122888989310960323D+08,
     *   -0.584534033302323818D+08,
     *   -0.327037777965903759D+09,
     *    0.744612559537177801D+09,
     *   -0.267760071339043856D+09,
     *    0.142003532436867714D+09,
     *    0.181188215378812790D+09,
     *   -0.460401108892367363D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H11CI_VL0_1=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H11 DIABAT WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H11CI_VL0_2(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.410723197014157677D+10,
     *    0.506689944720027828D+10,
     *   -0.154440423309340801D+11,
     *    0.126939494319864254D+11,
     *    0.335514547282957726D+11,
     *   -0.821361540854758453D+11,
     *    0.498088524164014359D+11,
     *   -0.105559793721869926D+11,
     *    0.362981176883468437D+10,
     *   -0.472645597316711617D+10,
     *    0.384012503683604279D+11,
     *   -0.583301527661400375D+11,
     *    0.263402032034487419D+11,
     *   -0.230395123565426540D+10,
     *    0.110066076602727008D+10,
     *   -0.144595151514403129D+10,
     *   -0.539098457904319763D+08,
     *   -0.323882813806885779D+09,
     *   -0.753713412393162131D+09,
     *    0.962182565826382935D+08,
     *   -0.558394725976745129D+09,
     *   -0.121416889170882940D+09,
     *   -0.540433904816019535D+08,
     *   -0.131504786818957329D+09,
     *    0.241187483892197609D+09,
     *    0.249764473024056911D+09,
     *    0.206639494773453236D+09,
     *    0.413757425026777506D+09,
     *    0.867946096652550697D+08,
     *    0.353138927211316228D+09,
     *   -0.120443162527495325D+09,
     *    0.315829632068882108D+09,
     *    0.265771125274375319D+09,
     *    0.315551104204180717D+09,
     *   -0.475071663331336021D+09,
     *    0.541717347237930298D+08,
     *    0.327515289140102386D+09,
     *   -0.222645768151259422D+09,
     *   -0.533186532045941353D+08,
     *    0.391068611823366165D+09,
     *   -0.208879066976512432D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H11CI_VL0_2=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H11 DIABAT WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H11CI_VL0_3(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.992851907650230408D+10,
     *   -0.123021873137497444D+11,
     *    0.528643781225123024D+10,
     *    0.204572061610103202D+10,
     *   -0.746832379888492432D+11,
     *    0.123360523679246964D+12,
     *   -0.599290450794751434D+11,
     *    0.978422179450627136D+10,
     *   -0.619488264249820614D+10,
     *    0.218860642894198799D+10,
     *    0.143144310799655609D+11,
     *   -0.383166609031813507D+11,
     *    0.250124761505249252D+11,
     *    0.841737059358996391D+09,
     *   -0.433116040464500999D+10,
     *    0.438576998983644867D+10,
     *   -0.171946923477388716D+10,
     *   -0.173010801201131105D+09,
     *    0.758537372477397680D+09,
     *   -0.903778177667283177D+09,
     *   -0.127693845466809273D+09,
     *    0.448666102151245594D+09,
     *   -0.100774780368533111D+10,
     *    0.560160777115999460D+09,
     *   -0.288476928764286995D+09,
     *    0.147368386186193794D+09,
     *    0.139614741773515701D+09,
     *    0.407759379931960106D+09,
     *    0.174228457067628860D+09,
     *   -0.277216981869235039D+08,
     *    0.579011793293105841D+09,
     *   -0.626220574198467374D+09,
     *    0.266252836004965544D+09,
     *   -0.174148819981927872D+08,
     *   -0.183495646226751328D+09,
     *    0.169598350686431885D+09,
     *   -0.159193173787382126D+09,
     *    0.445064189102616310D+09,
     *   -0.651015481688286781D+09,
     *    0.866148128024418354D+09,
     *   -0.461161080520584106D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H11CI_VL0_3=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H11 DIABAT WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H11CI_VL0_4(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.120341543222196884D+11,
     *   -0.128863529975174656D+11,
     *   -0.830959422539284611D+10,
     *    0.153167128937894478D+11,
     *   -0.245488854455680885D+11,
     *    0.109437039128152390D+11,
     *    0.139321314882062950D+11,
     *   -0.129880526733412552D+11,
     *    0.146433238262807484D+11,
     *   -0.185917011277836304D+11,
     *    0.473463994773544769D+11,
     *   -0.728080719212940826D+11,
     *    0.400270938101350708D+11,
     *   -0.535804305404485989D+10,
     *    0.138113171182600760D+10,
     *   -0.575079868253699064D+09,
     *    0.770347913640305758D+09,
     *   -0.808335078946387887D+09,
     *    0.554016089538181186D+09,
     *   -0.259348996479093075D+09,
     *   -0.350770436526856899D+09,
     *    0.534729178031826019D+09,
     *   -0.673994915547287941D+09,
     *    0.325325582670249462D+09,
     *   -0.168418134695476532D+09,
     *   -0.192545722035191536D+09,
     *    0.519147976541693687D+09,
     *   -0.273096854725134850D+09,
     *    0.837146929944561481D+09,
     *   -0.456275113768770218D+09,
     *    0.391741134843969822D+09,
     *   -0.339250064237517357D+09,
     *    0.539062873531436920D+08,
     *    0.292963763780741692D+08,
     *   -0.979049318684167862D+08,
     *   -0.138314825983955145D+09,
     *   -0.538012179377819598D+08,
     *    0.523355764679781914D+09,
     *   -0.539172825590674400D+09,
     *    0.368290924615418434D+09,
     *   -0.941859932006244659D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H11CI_VL0_4=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H11 DIABAT WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H11CI_VL0_5(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.125498875200508480D+11,
     *   -0.221503403816442986D+11,
     *    0.132539674282388535D+11,
     *   -0.165793403579035034D+11,
     *    0.157343111845110626D+12,
     *   -0.281168266528799011D+12,
     *    0.159512174716666046D+12,
     *   -0.315481090133780441D+11,
     *    0.106055458607758827D+11,
     *    0.193187596078476572D+10,
     *   -0.528254646651716232D+11,
     *    0.106110758341830627D+12,
     *   -0.638119272805368500D+11,
     *    0.688615982860062790D+10,
     *    0.291880709682986259D+09,
     *   -0.130332529343302608D+10,
     *    0.982397829352440834D+09,
     *   -0.304427693545155048D+09,
     *    0.125949776113070965D+09,
     *   -0.619641495195833445D+08,
     *   -0.807575978833029270D+08,
     *    0.239442522578252792D+09,
     *   -0.208530747546630561D+09,
     *    0.975347506528282166D+08,
     *   -0.163613499888999939D+09,
     *   -0.975819671306536198D+08,
     *    0.321140929496298015D+09,
     *   -0.326216914014689803D+09,
     *    0.608283250523697853D+09,
     *   -0.431383007504911423D+08,
     *   -0.584550073340663910D+08,
     *    0.625494956916294098D+08,
     *   -0.421317146245628357D+09,
     *    0.353143190766916037D+09,
     *   -0.166516648737144470D+09,
     *    0.140089913414128780D+09,
     *   -0.403137252403277278D+09,
     *    0.299042502016615868D+09,
     *    0.263776727560466766D+09,
     *   -0.512032699311965942D+09,
     *    0.279634935369014740D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H11CI_VL0_5=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H11 DIABAT WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H11CI_VL0_6(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.158429580428204422D+11,
     *   -0.392475834234195862D+11,
     *    0.422151220569631348D+11,
     *   -0.303552084302411041D+11,
     *    0.468877584058168640D+11,
     *   -0.834394184982988434D+11,
     *    0.594234901310046921D+11,
     *   -0.122456231719285278D+11,
     *   -0.177518515673693156D+10,
     *    0.156617360966087608D+11,
     *   -0.776533705533361206D+11,
     *    0.116982040399426025D+12,
     *   -0.572407676292626877D+11,
     *    0.607333842402847958D+10,
     *   -0.300860377378458786D+10,
     *    0.219078298267440701D+10,
     *   -0.529262671896023095D+09,
     *   -0.346481606663372517D+08,
     *    0.456309691841410160D+09,
     *   -0.595418867991861463D+09,
     *    0.565133197921074986D+09,
     *   -0.319167990095046759D+09,
     *    0.144135739119715214D+09,
     *    0.546489081183366776D+08,
     *   -0.744571156180148125D+08,
     *   -0.222036245032580853D+09,
     *    0.158716503765999556D+09,
     *   -0.335215427914468050D+08,
     *   -0.131001884711058140D+09,
     *    0.242683891304682732D+09,
     *    0.326813428381870270D+09,
     *   -0.256688731965733528D+09,
     *   -0.952716355257225037D+08,
     *   -0.101038120841587067D+09,
     *    0.461775469861051559D+09,
     *   -0.638569115879884720D+09,
     *    0.478106648700817108D+09,
     *   -0.239690063689750671D+09,
     *    0.114291720626981735D+09,
     *   -0.333090136301479340D+08,
     *   -0.131279910900278091D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H11CI_VL0_6=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H11 DIABAT WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H11CI_VL0_7(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.235705377045790215D+11,
     *   -0.707763062242668762D+11,
     *    0.906643620608644562D+11,
     *   -0.790252435957259979D+11,
     *    0.227206004164486053D+12,
     *   -0.372937872785397339D+12,
     *    0.213753420908294495D+12,
     *   -0.424440290543384705D+11,
     *    0.813920189415389252D+10,
     *    0.540402342093138313D+10,
     *   -0.229477035685399857D+11,
     *    0.360244331071565170D+11,
     *   -0.177739444476916733D+11,
     *    0.458625385328266144D+09,
     *    0.688043577526339531D+09,
     *   -0.558088772686892033D+09,
     *    0.407551855232027054D+09,
     *    0.303533904754509449D+09,
     *   -0.438404843449235916D+09,
     *    0.327886302706208944D+09,
     *   -0.913340912721986771D+08,
     *   -0.333664098970842361D+08,
     *    0.136450649229245186D+09,
     *   -0.832710079521808624D+08,
     *    0.363866914219799042D+08,
     *    0.893775473571090698D+08,
     *   -0.324690129095209599D+09,
     *    0.360308138888592005D+09,
     *   -0.348225215918651104D+09,
     *    0.184494372394289017D+09,
     *    0.219756066507015228D+08,
     *    0.122633635142797470D+09,
     *   -0.137241285174808502D+09,
     *   -0.945872674624633789D+08,
     *    0.148515535638841629D+09,
     *    0.128920027976055145D+09,
     *   -0.365487387449590683D+09,
     *    0.315589443536647797D+09,
     *   -0.243981173941619873D+09,
     *    0.177866917339542389D+09,
     *   -0.408803739876136780D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H11CI_VL0_7=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H11 DIABAT WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H11CI_VL0_8(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.564135212594949818D+10,
     *   -0.870175514122092056D+10,
     *    0.518302019313124657D+10,
     *   -0.657418361603513145D+10,
     *    0.325669964508629684D+11,
     *   -0.898362795632160187D+11,
     *    0.865583933756932678D+11,
     *   -0.349099306140149918D+11,
     *    0.162285338789425850D+11,
     *   -0.493001090936865139D+10,
     *   -0.383693729488838959D+11,
     *    0.762035127312288361D+11,
     *   -0.434583064592741852D+11,
     *    0.511576356421382809D+10,
     *   -0.118885112676238918D+10,
     *    0.360485073675737262D+09,
     *   -0.237766431331490159D+09,
     *    0.608460177465176940D+09,
     *   -0.342785895009413362D+09,
     *   -0.386683113200291395D+08,
     *    0.250793126312446356D+09,
     *   -0.272186699472331166D+09,
     *    0.228266375937395602D+09,
     *   -0.191331292893972754D+09,
     *    0.250080128822263896D+09,
     *   -0.284591702572794914D+09,
     *    0.302676732244371355D+09,
     *   -0.189082345528565526D+09,
     *   -0.210901296271149516D+08,
     *   -0.336359348894661307D+09,
     *    0.825686994071604252D+09,
     *   -0.520005164816785336D+09,
     *    0.154538871033468485D+09,
     *   -0.123763414209801197D+09,
     *   -0.922189677623782158D+08,
     *    0.283369809211835384D+09,
     *   -0.116038943630863667D+09,
     *    0.255746593313138485D+08,
     *   -0.284857344918354511D+09,
     *    0.458643099147464037D+09,
     *   -0.237642054573683023D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H11CI_VL0_8=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H11 DIABAT WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H11CI_VL0_9(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *   -0.695710226488211751D+08,
     *    0.106095930436345444D+11,
     *   -0.235440872827960205D+11,
     *    0.130683798313802719D+11,
     *    0.115062156176361938D+12,
     *   -0.250283761214301331D+12,
     *    0.160351006294754486D+12,
     *   -0.379433215279756012D+11,
     *    0.226735107592784958D+11,
     *   -0.194620819987529449D+11,
     *    0.345071721569521179D+11,
     *   -0.505415384452243729D+11,
     *    0.289313905817205772D+11,
     *   -0.417759420931861687D+10,
     *    0.675636523884788036D+09,
     *    0.359106481407902241D+08,
     *   -0.273161614073288441D+08,
     *    0.146823977753425598D+09,
     *    0.117505589793963447D+09,
     *   -0.347096754128953338D+09,
     *    0.392307900113504171D+09,
     *   -0.331157337774027586D+09,
     *    0.228375172078185081D+09,
     *   -0.127710055394392967D+09,
     *    0.765544435905094147D+08,
     *    0.416059905197749138D+08,
     *   -0.715463654455238581D+08,
     *    0.342687834537581503D+08,
     *   -0.126254295446875095D+09,
     *    0.243052854200771928D+09,
     *   -0.324338575275112629D+09,
     *    0.189800997133445144D+09,
     *    0.333861041360360384D+08,
     *   -0.775859587310514450D+08,
     *    0.203473883804778099D+09,
     *   -0.551976660357192993D+09,
     *    0.799864970193072319D+09,
     *   -0.535338215010530472D+09,
     *    0.916798744582834244D+08,
     *    0.152509423100502491D+09,
     *   -0.137718269556090593D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H11CI_VL0_9=SUMA
       RETURN
       END

      FUNCTION VCI_H22(R,theta)
C*********************************
C System: Kr-OH(X2Pi - A2Sigma) MIXING
C Method: MRCISD+Q
C Basis: AVTZ-DK
C with no bond functions
C Size Consistency Corrected
C PES: H22 CI DIABAT
C ASYMPTOTICS IS E=0 for R=>INFTY
C TO GET ASYMPTOTICS TOWARDS OH(A) STATE
C ADD 32592.4 cm-1
C Variables:
C R in Bohr
C Theta n Degrees
C Output in cm-1
C Theta=0 for Kr---HO geom.
C********************************
C Authors:
C Dr Jacek Klos
C Department of Chemistry and Biochemistry
C University of Maryland
C College Park, MD 20742
C email:jklos@umd.edu
C**********************************
      implicit double precision(a-h, o-z)
      dimension V0(11)
      dimension T(11)
      pi=dacos(-1.d0)
      conv=1.D0 ! IN CM-1
      V0(1)=H22CI_VL0_0(r)*conv    
      V0(2)=H22CI_VL0_1(r)*conv
      V0(3)=H22CI_VL0_2(r)*conv
      V0(4)=H22CI_VL0_3(r)*conv
      V0(5)=H22CI_VL0_4(r)*conv
      V0(6)=H22CI_VL0_5(r)*conv
      V0(7)=H22CI_VL0_6(r)*conv
      V0(8)=H22CI_VL0_7(r)*conv
      V0(9)=H22CI_VL0_8(r)*conv
      V0(10)=H22CI_VL0_9(r)*conv
      V0(11)=H22CI_VL0_10(r)*conv
       do j=1,11
       T(j)=PLGNDR((j-1),0,DCOS(theta*pi/180.D0))
       enddo
       s=0.D0
       do i=1,11
       s=s+V0(i)*T(i)
       enddo
       VCI_H22=s
       return
       end
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H22 DIABAT WITH
C RIGID OH(A) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H22CI_VL0_0(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.331257857517450523D+10,
     *    0.101269218665445995D+10,
     *   -0.774039525187470245D+10,
     *    0.579853429975492191D+10,
     *    0.160603959743558693D+11,
     *   -0.335425173313512306D+11,
     *    0.221876780530350838D+11,
     *   -0.933596530737809944D+10,
     *    0.109739978507358189D+11,
     *   -0.128686257308045197D+11,
     *    0.926889848394235992D+10,
     *   -0.494754885139596272D+10,
     *    0.889321949314870834D+09,
     *   -0.130710401965448976D+10,
     *    0.602884442645624638D+09,
     *   -0.991620220321126461D+09,
     *    0.172998191789619565D+09,
     *   -0.614346722910949469D+09,
     *   -0.531405761974209845D+08,
     *   -0.369387381963611901D+09,
     *   -0.234673783705638528D+09,
     *    0.979342736551548243D+08,
     *   -0.253719434933477163D+09,
     *    0.173917352487463713D+09,
     *   -0.119304812169454575D+09,
     *    0.134676644045323670D+09,
     *    0.100313349983981133D+09,
     *    0.242224814890018582D+09,
     *    0.221688770573425770D+09,
     *    0.517684244518345952D+09,
     *   -0.419106346985851169D+09,
     *    0.659372483184395313D+09,
     *   -0.332821109977515340D+08,
     *    0.427178763582198620D+08,
     *    0.576147185015592575D+08,
     *   -0.197105399693387270D+09,
     *    0.852698674388120174D+08,
     *   -0.156123538805770636D+09,
     *    0.632643188327219248D+09,
     *   -0.266298012930672169D+10,
     *    0.372970228077689123D+10/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H22CI_VL0_0=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H22 DIABAT WITH
C RIGID OH(A) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H22CI_VL0_10(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *   -0.508021105796177959D+10,
     *    0.180759218787667427D+11,
     *   -0.218703044268266258D+11,
     *    0.140617889891763554D+11,
     *   -0.463479929372833710D+11,
     *    0.888135592320015106D+11,
     *   -0.607735467736958923D+11,
     *    0.145110284681325989D+11,
     *   -0.139721952130599403D+10,
     *    0.687952304439463043D+10,
     *   -0.426489964393608551D+11,
     *    0.684168558187577286D+11,
     *   -0.358542978372004089D+11,
     *    0.465539849143638420D+10,
     *   -0.239795573677840996D+10,
     *    0.127031328495244241D+10,
     *   -0.488671767963995516D+09,
     *    0.156621317605828404D+09,
     *    0.105437228670125395D+09,
     *   -0.828047366822471917D+08,
     *    0.273147664980888367D+05,
     *   -0.958044161089720726D+08,
     *    0.800744140039055347D+08,
     *    0.113015239686178446D+09,
     *   -0.162632582294077873D+08,
     *   -0.394350898896825790D+09,
     *    0.510247599942516685D+09,
     *    0.414476773334602416D+08,
     *   -0.731954697039530277D+09,
     *    0.642977286095053911D+09,
     *   -0.785513515611159801D+08,
     *   -0.165490276123747826D+08,
     *   -0.235776045520183563D+09,
     *    0.453295658245712280D+09,
     *   -0.487402755625935555D+09,
     *    0.102274736467597008D+09,
     *    0.493763395650296688D+09,
     *   -0.628066274795193672D+09,
     *    0.522369576102302074D+09,
     *   -0.478183110835051537D+09,
     *    0.207404096611095428D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H22CI_VL0_10=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H22 DIABAT WITH
C RIGID OH(A) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H22CI_VL0_1(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.297939180427505493D+10,
     *    0.705363651238331604D+10,
     *   -0.191633892280233421D+11,
     *    0.219623756996385956D+11,
     *   -0.859324250921464386D+11,
     *    0.138107926145799561D+12,
     *   -0.758935897185085144D+11,
     *    0.127104043980402164D+11,
     *    0.100866084956903744D+10,
     *    0.682381015780406952D+09,
     *   -0.313542213154001274D+11,
     *    0.556194746964961090D+11,
     *   -0.285681959818380737D+11,
     *    0.279650387006481171D+10,
     *   -0.436874696463905621D+10,
     *    0.114947995049565172D+10,
     *   -0.123627525765800524D+10,
     *   -0.149942460171507090D+09,
     *    0.101123490848018396D+10,
     *   -0.400043399003867149D+09,
     *    0.450705562197506905D+09,
     *    0.718037732136734009D+09,
     *   -0.546480432904430389D+09,
     *    0.671300402503577471D+09,
     *    0.467272537356419563D+08,
     *   -0.340277624435527325D+08,
     *    0.494541733336481810D+09,
     *   -0.104305047786185265D+09,
     *    0.525777787933301449D+09,
     *   -0.805839406779500723D+09,
     *    0.119934573249440980D+10,
     *   -0.764806320635955572D+09,
     *    0.477945284368544579D+09,
     *   -0.595275758457131863D+09,
     *    0.402324321655319154D+09,
     *   -0.511367879858423471D+09,
     *    0.356836003404825211D+09,
     *    0.203581148054740906D+09,
     *   -0.102659221418133736D+09,
     *    0.372970683054441452D+09,
     *   -0.544796535278690338D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H22CI_VL0_1=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H22 DIABAT WITH
C RIGID OH(A) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H22CI_VL0_2(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.141618659134991493D+11,
     *   -0.200125652820648804D+11,
     *    0.949165603351620102D+10,
     *   -0.472217659784441757D+10,
     *   -0.511630874531631622D+11,
     *    0.102430419769087463D+12,
     *   -0.585008036456834106D+11,
     *    0.108039994393630714D+11,
     *   -0.477667102654557419D+10,
     *    0.401046480287586975D+10,
     *   -0.390904676683183365D+11,
     *    0.584314732260088501D+11,
     *   -0.271508005279040909D+11,
     *    0.192506374420333481D+10,
     *   -0.730057928690143824D+09,
     *    0.206352329486082602D+10,
     *    0.567520696868871450D+09,
     *    0.633675045931585670D+09,
     *    0.107976529583816195D+10,
     *    0.142549038772753656D+09,
     *    0.669233477085536003D+09,
     *    0.153323507083946705D+09,
     *    0.792490285603260994D+07,
     *    0.342832463500521541D+09,
     *    0.379201347861531973D+08,
     *    0.144188334262307763D+09,
     *    0.168776698264168620D+09,
     *    0.115724085216129541D+09,
     *   -0.243537940738795280D+09,
     *    0.262332248972253799D+09,
     *   -0.269825131743741035D+09,
     *   -0.430229734743182182D+09,
     *   -0.856774094231843948D+08,
     *    0.506431862946033478D+07,
     *   -0.179485159350826263D+09,
     *   -0.167334903096955299D+09,
     *   -0.297922388442811012D+09,
     *    0.364080944401950836D+09,
     *   -0.321421922298557758D+09,
     *    0.393513390238523901D+09,
     *   -0.294006256666805744D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H22CI_VL0_2=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H22 DIABAT WITH
C RIGID OH(A) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H22CI_VL0_3(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.125399950630822334D+11,
     *   -0.282750637581199608D+11,
     *    0.321028112832101898D+11,
     *   -0.327457383938403702D+11,
     *    0.112728322679523193D+12,
     *   -0.164012987174408264D+12,
     *    0.741692403929079742D+11,
     *   -0.117103241193093681D+11,
     *    0.693951753369677258D+10,
     *   -0.213972251626498842D+10,
     *   -0.140079119838218708D+11,
     *    0.382258094343669434D+11,
     *   -0.247418779647072601D+11,
     *   -0.964918549357278824D+09,
     *    0.376473225104815340D+10,
     *   -0.490777831996900463D+10,
     *    0.137179001610179377D+10,
     *   -0.429550889528766274D+08,
     *   -0.705696452732154846D+09,
     *    0.106481670265207505D+10,
     *    0.507125005200473309D+09,
     *   -0.681957433312105179D+09,
     *    0.126367865165457869D+10,
     *   -0.486258971188908696D+09,
     *    0.183846760863800049D+09,
     *    0.735447901918906927D+09,
     *   -0.553423098464929104D+09,
     *    0.769121579070619583D+09,
     *   -0.756905534531563520D+08,
     *    0.487178153822612762D+09,
     *   -0.103529220938797474D+10,
     *    0.903529871832629204D+09,
     *   -0.573163139597532272D+09,
     *    0.108745037300357342D+09,
     *    0.282775327229192257D+08,
     *    0.605400033642954826D+08,
     *   -0.101331029580907536D+10,
     *    0.118330780108559990D+10,
     *   -0.841379889555204391D+09,
     *    0.774075751450897217D+09,
     *   -0.396229998455173492D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H22CI_VL0_3=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H22 DIABAT WITH
C RIGID OH(A) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H22CI_VL0_4(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.904257655095584106D+10,
     *   -0.148000563520735207D+11,
     *    0.156796679058230953D+11,
     *   -0.174155932344628830D+11,
     *    0.114435252896499825D+11,
     *    0.142873994615531635D+11,
     *   -0.267199319497228317D+11,
     *    0.150054889171965446D+11,
     *   -0.150044832501446075D+11,
     *    0.188047971468357887D+11,
     *   -0.474936961449331512D+11,
     *    0.729715418723686218D+11,
     *   -0.403677647010647354D+11,
     *    0.483140933649440765D+10,
     *   -0.182468987169125652D+10,
     *    0.456151167308328748D+09,
     *   -0.811731513209398627D+09,
     *    0.734857058907907128D+09,
     *   -0.535155008072624683D+09,
     *    0.463477865385331392D+09,
     *    0.822166516935206652D+09,
     *   -0.661185576503141284D+09,
     *    0.970241353205680370D+09,
     *   -0.678436029019505382D+09,
     *    0.796281453994968534D+09,
     *   -0.420984203306871176D+09,
     *    0.500847319539754391D+09,
     *    0.188570518194947660D+09,
     *   -0.416988227208676457D+09,
     *    0.867889190817394733D+09,
     *    0.301955892800684929D+09,
     *   -0.117799162961494541D+10,
     *    0.685152685603351116D+09,
     *   -0.608825895866374493D+09,
     *    0.470259115645734310D+09,
     *   -0.406842380828254223D+09,
     *   -0.590035501488684177D+09,
     *    0.104746458659519386D+10,
     *   -0.663377182609004974D+09,
     *    0.313300214849209785D+09,
     *   -0.777402479579076767D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H22CI_VL0_4=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H22 DIABAT WITH
C RIGID OH(A) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H22CI_VL0_5(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.630476659952422047D+10,
     *    0.122888744652420521D+10,
     *   -0.163489747512033081D+11,
     *    0.201390294848125534D+11,
     *   -0.140238881256418457D+12,
     *    0.248729697621020660D+12,
     *   -0.142250469961313782D+12,
     *    0.293751955400679169D+11,
     *   -0.984275863178364563D+10,
     *   -0.208815736570602798D+10,
     *    0.530084855111134415D+11,
     *   -0.106249190451144989D+12,
     *    0.639701664923305969D+11,
     *   -0.677029769559156418D+10,
     *   -0.216800290553425312D+09,
     *    0.140539679597218513D+10,
     *   -0.962138062209702492D+09,
     *    0.198175688987132579D+09,
     *   -0.305274936654112816D+09,
     *    0.264775068052698761D+09,
     *    0.182692353638484001D+09,
     *    0.177254287540200233D+09,
     *   -0.104374069175613403D+09,
     *    0.326989475377894521D+09,
     *   -0.135975602129821777D+09,
     *   -0.897397573859150410D+08,
     *    0.154191926988011360D+09,
     *    0.212807394888031006D+09,
     *   -0.431633580050520897D+08,
     *    0.442224179408937097D+08,
     *    0.436737486252130747D+09,
     *   -0.125935454803923607D+09,
     *   -0.815033408600414753D+09,
     *    0.770316099803897500D+09,
     *   -0.220230536492660046D+09,
     *    0.798598831374373436D+08,
     *   -0.742257961119610548D+09,
     *    0.603510484845246077D+09,
     *    0.188445317443236351D+09,
     *   -0.493584037299928665D+09,
     *    0.258809324605190277D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H22CI_VL0_5=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H22 DIABAT WITH
C RIGID OH(A) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H22CI_VL0_6(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *   -0.760247402395385933D+10,
     *    0.432155050262567749D+11,
     *   -0.691925385270851135D+11,
     *    0.543563177516326981D+11,
     *   -0.101524452021851471D+12,
     *    0.155725525938389343D+12,
     *   -0.914761190715946960D+11,
     *    0.163208301294383392D+11,
     *    0.684654244743764877D+09,
     *   -0.150364948521558533D+11,
     *    0.778347641108581390D+11,
     *   -0.117087137576079041D+12,
     *    0.576848630497441406D+11,
     *   -0.573423449013644695D+10,
     *    0.315823444244393730D+10,
     *   -0.203014682272218513D+10,
     *    0.663965985445343494D+09,
     *    0.280813237949562073D+07,
     *   -0.636542733463556767D+09,
     *    0.579133906826222181D+09,
     *   -0.363331483043745041D+09,
     *    0.327974933616406441D+09,
     *    0.936080772866020203D+08,
     *   -0.143116374987736702D+09,
     *    0.137166061981201172D+06,
     *    0.594919251567821503D+09,
     *   -0.869735884937392235D+09,
     *    0.282753250794470787D+09,
     *    0.856435829147543460D+08,
     *    0.658662026033361435D+09,
     *   -0.143020684351741409D+10,
     *    0.129306472517292595D+10,
     *   -0.388679714408124924D+09,
     *   -0.254954067138753891D+09,
     *    0.566236278803352356D+09,
     *   -0.716232167310655594D+09,
     *    0.770890832057651520D+09,
     *   -0.637247130014369965D+09,
     *    0.344920278181163788D+09,
     *   -0.122641211597633362D+09,
     *   -0.375641061753082275D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H22CI_VL0_6=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H22 DIABAT WITH
C RIGID OH(A) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H22CI_VL0_7(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.138366879515616918D+10,
     *    0.116946304897078819D+11,
     *   -0.318405540963202896D+11,
     *    0.355610927098372879D+11,
     *   -0.145082723893493164D+12,
     *    0.266905282314706970D+12,
     *   -0.168087509219721619D+12,
     *    0.366732488811457672D+11,
     *   -0.631243089389179420D+10,
     *   -0.583615482725802612D+10,
     *    0.233666330227527046D+11,
     *   -0.362430267433386307D+11,
     *    0.180479290338031311D+11,
     *   -0.242495324493152618D+09,
     *   -0.548406225854614258D+09,
     *    0.699365134940030217D+09,
     *   -0.136884900215273857D+09,
     *   -0.196877025348369598D+09,
     *    0.391716674245104790D+09,
     *   -0.491463004230932891D+09,
     *    0.139563407044875622D+09,
     *    0.881439323175950050D+08,
     *    0.464208202999744415D+08,
     *    0.146305306109333038D+08,
     *    0.128569558383989334D+07,
     *   -0.317201769934878349D+08,
     *    0.382762088182897568D+08,
     *   -0.542259606711673737D+07,
     *   -0.110862062478926182D+09,
     *   -0.166441271228904009D+09,
     *    0.955761277195309401D+09,
     *   -0.861847675338696480D+09,
     *    0.248021630060190201D+09,
     *   -0.170828341306333542D+08,
     *   -0.257709516584201813D+09,
     *    0.167796480722738266D+09,
     *    0.708509414111170292D+09,
     *   -0.992075011941078186D+09,
     *    0.437527096979688644D+09,
     *   -0.133096256363800049D+09,
     *    0.320088498442859650D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H22CI_VL0_7=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H22 DIABAT WITH
C RIGID OH(A) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H22CI_VL0_8(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.148365170431118317D+11,
     *   -0.385853655083746185D+11,
     *    0.373155452978020630D+11,
     *   -0.157633244097904243D+11,
     *   -0.494255224850299072D+11,
     *    0.140093870876024658D+12,
     *   -0.117183014025952637D+12,
     *    0.387895514072212982D+11,
     *   -0.172874191463605843D+11,
     *    0.531620443171011257D+10,
     *    0.383527594731639709D+11,
     *   -0.762450617358952789D+11,
     *    0.435731680044964294D+11,
     *   -0.501019532401720142D+10,
     *    0.125878817281413579D+10,
     *   -0.291606269695208311D+09,
     *    0.484822068646267772D+09,
     *   -0.475159153109042466D+09,
     *    0.437258718947526991D+09,
     *   -0.294279204009179473D+09,
     *   -0.147027761838858604D+09,
     *    0.460324410096245289D+09,
     *   -0.382532023951336384D+09,
     *    0.346872345834131241D+09,
     *   -0.329078745084552526D+09,
     *    0.147530532591603756D+09,
     *    0.352237249335796475D+09,
     *   -0.133869956800350380D+10,
     *    0.265614911629975605D+10,
     *   -0.333845065299302292D+10,
     *    0.214562625468936062D+10,
     *   -0.213365519376773834D+09,
     *   -0.492323151255650520D+09,
     *    0.506470776753686905D+09,
     *   -0.750512759098478317D+09,
     *    0.438939081040187836D+09,
     *    0.946578898973796844D+09,
     *   -0.129152497195089149D+10,
     *    0.388073909977273941D+09,
     *    0.145255813989465714D+09,
     *   -0.156656852642167091D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H22CI_VL0_8=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X-A) H22 DIABAT WITH
C RIGID OH(A) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION H22CI_VL0_9(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .   2.750D0,
     .   3.000D0,
     .   3.250D0,
     .   3.500D0,
     .   3.700D0,
     .   3.750D0,
     .   3.800D0,
     .   3.900D0,
     .   4.000D0,
     .   4.100D0,
     .   4.200D0,
     .   4.250D0,
     .   4.300D0,
     .   4.500D0,
     .   4.750D0,
     .   5.000D0,
     .   5.250D0,
     .   5.500D0,
     .   5.750D0,
     .   6.000D0,
     .   6.250D0,
     .   6.500D0,
     .   6.750D0,
     .   7.000D0,
     .   7.250D0,
     .   7.500D0,
     .   7.750D0,
     .   8.000D0,
     .   8.250D0,
     .   8.500D0,
     .   8.750D0,
     .   9.000D0,
     .   9.500D0,
     .  10.000D0,
     .  10.500D0,
     .  11.000D0,
     .  11.500D0,
     .  12.000D0,
     .  13.000D0,
     .  14.000D0,
     .  15.000D0/
      DATA VCOEF/
     *    0.116601421993027258D+10,
     *    0.366358963381005645D+09,
     *   -0.180588166024535990D+10,
     *    0.290856147203425503D+10,
     *   -0.964039278364802704D+11,
     *    0.204118906715059021D+12,
     *   -0.133327067705672699D+12,
     *    0.345565290124969406D+11,
     *   -0.216966899464147491D+11,
     *    0.191407035605465775D+11,
     *   -0.342804973939449692D+11,
     *    0.503768648560288544D+11,
     *   -0.288343566477055664D+11,
     *    0.423680106257506704D+10,
     *   -0.675277586933554530D+09,
     *   -0.232933645426276922D+08,
     *    0.158345645615461588D+09,
     *    0.799898840430558324D+08,
     *   -0.148067933729498506D+09,
     *    0.271166894923068047D+09,
     *   -0.453907215313519955D+09,
     *    0.375028960447992206D+09,
     *   -0.626919785495548248D+08,
     *   -0.212001377055881023D+09,
     *    0.320375858205129981D+09,
     *   -0.165112304690139771D+09,
     *   -0.801115021161909103D+08,
     *    0.363112057685965300D+09,
     *   -0.795121203225088954D+09,
     *    0.150337274056097078D+10,
     *   -0.202694162437012315D+10,
     *    0.118267341322977924D+10,
     *   -0.392394846222305298D+07,
     *   -0.878075787564081252D+08,
     *   -0.177672456320217609D+09,
     *   -0.393229664459240913D+09,
     *    0.173027824256668091D+10,
     *   -0.170092687603166437D+10,
     *    0.663520494917922378D+09,
     *   -0.110269121655395031D+09,
     *   -0.646616443126935959D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       H22CI_VL0_9=SUMA
       RETURN
       END
      FUNCTION VCI_XAbis(R,theta)
C*********************************
C System: Kr-OH(X2Pi)
C Method: MRCISD+Q
C Basis: AVTZ-DK
C with no bond functions
C Size Consistency Corrected
C PES: 1A"
C Variables:
C R in Bohr
C Theta n Degrees
C Output in cm-1
C Theta=0 for Kr---HO geom.
C********************************
C Authors:
C Dr Jacek Klos
C Department of Chemistry and Biochemistry
C University of Maryland
C College Park, MD 20742
C email:jklos@umd.edu
C**********************************
      implicit double precision(a-h, o-z)
      dimension V0(11)
      dimension T(11)
      pi=dacos(-1.d0)
      conv=1.D0 ! IN CM-1
      V0(1)=VL0_0(r)*conv    
      V0(2)=VL0_1(r)*conv
      V0(3)=VL0_2(r)*conv
      V0(4)=VL0_3(r)*conv
      V0(5)=VL0_4(r)*conv
      V0(6)=VL0_5(r)*conv
      V0(7)=VL0_6(r)*conv
      V0(8)=VL0_7(r)*conv
      V0(9)=VL0_8(r)*conv
      V0(10)=VL0_9(r)*conv
      V0(11)=VL0_10(r)*conv
       do j=1,11
       T(j)=PLGNDR((j-1),0,DCOS(theta*pi/180.D0))
       enddo
       s=0.D0
       do i=1,11
       s=s+V0(i)*T(i)
       enddo
       VCI_XAbis=s
       return
       end

C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X) 1A'' PES WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION VL0_0(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .     2.500D0,
     .     2.750D0,
     .     3.000D0,
     .     3.250D0,
     .     3.500D0,
     .     3.700D0,
     .     3.750D0,
     .     3.800D0,
     .     3.900D0,
     .     4.000D0,
     .     4.100D0,
     .     4.200D0,
     .     4.250D0,
     .     4.300D0,
     .     4.500D0,
     .     4.750D0,
     .     5.000D0,
     .     5.250D0,
     .     5.500D0,
     .     5.750D0,
     .     6.000D0,
     .     6.250D0,
     .     6.500D0,
     .     6.750D0,
     .     7.000D0,
     .     7.250D0,
     .     7.500D0,
     .     7.750D0,
     .     8.000D0,
     .     8.250D0,
     .     8.500D0,
     .     8.750D0,
     .     9.000D0,
     .     9.500D0,
     .    10.000D0,
     .    10.500D0,
     .    11.000D0,
     .    12.000D0,
     .    13.000D0,
     .    14.000D0,
     .    15.000D0/
      DATA VCOEF/
     *    0.444247880320973873D+10,
     *   -0.560340406058358765D+10,
     *    0.387283470428195095D+10,
     *    0.764889080519841909D+09,
     *    0.111268963736257315D+09,
     *   -0.163557736496400690D+10,
     *    0.426049298336184168D+10,
     *   -0.353119123472344637D+10,
     *    0.748000552164716005D+09,
     *   -0.421945555365309417D+09,
     *   -0.116026424818402290D+09,
     *   -0.250227135494054079D+09,
     *    0.929722054089851081D+08,
     *   -0.405593321142051578D+09,
     *   -0.616668169763036013D+09,
     *   -0.745074661577533484D+09,
     *   -0.763320925223823309D+09,
     *   -0.760742623446693301D+09,
     *   -0.704267995581511259D+09,
     *   -0.600707862776609778D+09,
     *   -0.439556564962812066D+09,
     *   -0.286026622326476157D+09,
     *   -0.146452363726243734D+09,
     *   -0.409539871993731260D+08,
     *    0.574862229128122330D+06,
     *    0.906989542317146063D+08,
     *    0.138199419082749248D+09,
     *    0.221944410127321601D+09,
     *    0.308571557519117117D+09,
     *    0.361712318698771596D+09,
     *    0.342374057644880533D+09,
     *    0.178933602453395605D+09,
     *    0.399875773426320553D+09,
     *    0.208912956530405879D+09,
     *    0.128367162238637686D+09,
     *   -0.807873800811969042D+08,
     *    0.355929259386969805D+08,
     *   -0.852882693945827484D+08,
     *    0.135196295233482599D+09,
     *    0.148595452770279378D+09,
     *    0.272803559566966236D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_0=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X) 1A'' PES WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION VL0_10(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .     2.500D0,
     .     2.750D0,
     .     3.000D0,
     .     3.250D0,
     .     3.500D0,
     .     3.700D0,
     .     3.750D0,
     .     3.800D0,
     .     3.900D0,
     .     4.000D0,
     .     4.100D0,
     .     4.200D0,
     .     4.250D0,
     .     4.300D0,
     .     4.500D0,
     .     4.750D0,
     .     5.000D0,
     .     5.250D0,
     .     5.500D0,
     .     5.750D0,
     .     6.000D0,
     .     6.250D0,
     .     6.500D0,
     .     6.750D0,
     .     7.000D0,
     .     7.250D0,
     .     7.500D0,
     .     7.750D0,
     .     8.000D0,
     .     8.250D0,
     .     8.500D0,
     .     8.750D0,
     .     9.000D0,
     .     9.500D0,
     .    10.000D0,
     .    10.500D0,
     .    11.000D0,
     .    12.000D0,
     .    13.000D0,
     .    14.000D0,
     .    15.000D0/
      DATA VCOEF/
     *    0.389617816256470776D+10,
     *   -0.692707135882719040D+10,
     *    0.269001560400565052D+10,
     *    0.122913782296267557D+10,
     *   -0.223127947925571346D+10,
     *   -0.203103428343693018D+09,
     *    0.596830536349816990D+10,
     *   -0.582412953799997711D+10,
     *    0.159022271586371970D+10,
     *   -0.434892598288800716D+09,
     *    0.171067240302111149D+09,
     *   -0.110668960609132484D+09,
     *    0.147412500071243227D+09,
     *   -0.376148970531623960D+08,
     *    0.404192297770532370D+08,
     *   -0.856723685244109482D+07,
     *   -0.245788896678175926D+08,
     *    0.305296722807315439D+08,
     *    0.422324398002459705D+08,
     *    0.219163176086535454D+08,
     *   -0.381823625194964111D+08,
     *    0.297897621706753969D+06,
     *   -0.121173991655623317D+08,
     *    0.462495095227683783D+08,
     *   -0.252665469852033854D+08,
     *   -0.935217418690514565D+07,
     *    0.412762914817135334D+08,
     *   -0.672736670853142738D+08,
     *    0.155356306066228151D+09,
     *   -0.153088134428171158D+09,
     *   -0.101201638553743124D+09,
     *    0.159478778743889213D+09,
     *    0.103931114639715910D+09,
     *   -0.407526284208165050D+09,
     *    0.680773662733382583D+09,
     *   -0.752267422086614609D+09,
     *    0.481763736009705305D+09,
     *   -0.266801627328372836D+09,
     *    0.249074460059555352D+09,
     *   -0.114956213131740808D+09,
     *   -0.185531027235329151D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_10=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X) 1A'' PES WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION VL0_1(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .     2.500D0,
     .     2.750D0,
     .     3.000D0,
     .     3.250D0,
     .     3.500D0,
     .     3.700D0,
     .     3.750D0,
     .     3.800D0,
     .     3.900D0,
     .     4.000D0,
     .     4.100D0,
     .     4.200D0,
     .     4.250D0,
     .     4.300D0,
     .     4.500D0,
     .     4.750D0,
     .     5.000D0,
     .     5.250D0,
     .     5.500D0,
     .     5.750D0,
     .     6.000D0,
     .     6.250D0,
     .     6.500D0,
     .     6.750D0,
     .     7.000D0,
     .     7.250D0,
     .     7.500D0,
     .     7.750D0,
     .     8.000D0,
     .     8.250D0,
     .     8.500D0,
     .     8.750D0,
     .     9.000D0,
     .     9.500D0,
     .    10.000D0,
     .    10.500D0,
     .    11.000D0,
     .    12.000D0,
     .    13.000D0,
     .    14.000D0,
     .    15.000D0/
      DATA VCOEF/
     *    0.479807236989410019D+10,
     *   -0.463919785129395390D+10,
     *    0.188726283929634023D+10,
     *   -0.226095092276942682D+10,
     *    0.495854666895620108D+09,
     *   -0.446780207808389187D+10,
     *    0.104916671150613766D+11,
     *   -0.799942157592443371D+10,
     *    0.212061829331425261D+10,
     *   -0.533131493599532843D+09,
     *    0.250833772067033172D+09,
     *   -0.358332364736093283D+08,
     *    0.465423932509451509D+08,
     *    0.628954669374675751D+08,
     *    0.588887404662029147D+08,
     *   -0.336945788608616590D+07,
     *   -0.666754348890333176D+08,
     *   -0.158830183738280952D+09,
     *   -0.225413392116856813D+09,
     *   -0.215051705389508396D+09,
     *   -0.175674571248754591D+09,
     *   -0.175596962732050449D+09,
     *   -0.193357813122147202D+09,
     *   -0.160260690570041299D+09,
     *   -0.500501830586522222D+08,
     *    0.139968958517110348D+07,
     *    0.166846487988615751D+09,
     *    0.901734825486818552D+08,
     *    0.112774194899466276D+09,
     *    0.232340616867304087D+09,
     *    0.795104849857158661D+08,
     *   -0.722470399758901596D+08,
     *    0.274872858868443847D+09,
     *    0.149449662945100546D+09,
     *    0.405262778044469357D+08,
     *   -0.297738570481367469D+09,
     *    0.281541938466142774D+09,
     *    0.191622829021427333D+09,
     *    0.354722362572936416D+08,
     *    0.340402713939232588D+09,
     *   -0.567549556781368017D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_1=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X) 1A'' PES WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION VL0_2(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .     2.500D0,
     .     2.750D0,
     .     3.000D0,
     .     3.250D0,
     .     3.500D0,
     .     3.700D0,
     .     3.750D0,
     .     3.800D0,
     .     3.900D0,
     .     4.000D0,
     .     4.100D0,
     .     4.200D0,
     .     4.250D0,
     .     4.300D0,
     .     4.500D0,
     .     4.750D0,
     .     5.000D0,
     .     5.250D0,
     .     5.500D0,
     .     5.750D0,
     .     6.000D0,
     .     6.250D0,
     .     6.500D0,
     .     6.750D0,
     .     7.000D0,
     .     7.250D0,
     .     7.500D0,
     .     7.750D0,
     .     8.000D0,
     .     8.250D0,
     .     8.500D0,
     .     8.750D0,
     .     9.000D0,
     .     9.500D0,
     .    10.000D0,
     .    10.500D0,
     .    11.000D0,
     .    12.000D0,
     .    13.000D0,
     .    14.000D0,
     .    15.000D0/
      DATA VCOEF/
     *    0.815945354062853622D+10,
     *   -0.844254544423920918D+10,
     *    0.218627941845488501D+10,
     *   -0.287139636112785339D+10,
     *    0.141260274153377485D+10,
     *   -0.425567041350677347D+10,
     *    0.916024798592223549D+10,
     *   -0.679735883666688156D+10,
     *    0.178493001043066120D+10,
     *   -0.465974368930305839D+09,
     *    0.173915332490313053D+09,
     *   -0.489296681093040407D+08,
     *    0.276769587545170784D+08,
     *   -0.282000504059414864D+08,
     *   -0.582756208425753713D+08,
     *   -0.619414592238248587D+08,
     *   -0.749616015486735106D+08,
     *   -0.146936169979472399D+09,
     *   -0.221990577350172043D+09,
     *   -0.196477393373019338D+09,
     *   -0.173562168198243499D+09,
     *   -0.164041025108932614D+09,
     *   -0.199350294118442774D+09,
     *   -0.931654397839841843D+08,
     *    0.276586109935564995D+08,
     *    0.193423212332655191D+09,
     *    0.193893921641383410D+09,
     *    0.267662450211387157D+09,
     *    0.304131740158772707D+09,
     *    0.595959430226509571D+08,
     *    0.222751144980378389D+09,
     *   -0.121748205906846046D+09,
     *    0.259399646350972652D+08,
     *    0.358243528502864122D+09,
     *   -0.405540539399476051D+08,
     *   -0.221853355434982777D+09,
     *    0.840694511087305546D+08,
     *    0.164991508578667641D+09,
     *   -0.395141699508907318D+09,
     *    0.577334423244143724D+09,
     *   -0.335371764707724631D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_2=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X) 1A'' PES WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION VL0_3(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .     2.500D0,
     .     2.750D0,
     .     3.000D0,
     .     3.250D0,
     .     3.500D0,
     .     3.700D0,
     .     3.750D0,
     .     3.800D0,
     .     3.900D0,
     .     4.000D0,
     .     4.100D0,
     .     4.200D0,
     .     4.250D0,
     .     4.300D0,
     .     4.500D0,
     .     4.750D0,
     .     5.000D0,
     .     5.250D0,
     .     5.500D0,
     .     5.750D0,
     .     6.000D0,
     .     6.250D0,
     .     6.500D0,
     .     6.750D0,
     .     7.000D0,
     .     7.250D0,
     .     7.500D0,
     .     7.750D0,
     .     8.000D0,
     .     8.250D0,
     .     8.500D0,
     .     8.750D0,
     .     9.000D0,
     .     9.500D0,
     .    10.000D0,
     .    10.500D0,
     .    11.000D0,
     .    12.000D0,
     .    13.000D0,
     .    14.000D0,
     .    15.000D0/
      DATA VCOEF/
     *    0.930564831881555367D+10,
     *   -0.114889138860875282D+11,
     *    0.377102227263774633D+10,
     *   -0.174040440651176643D+10,
     *   -0.352528258430121064D+09,
     *    0.149067195807482862D+10,
     *   -0.183815258735292983D+10,
     *    0.672149261705048203D+09,
     *   -0.163740885582811326D+09,
     *    0.493277706259444058D+08,
     *   -0.380734215891093016D+07,
     *    0.316401046513360739D+08,
     *   -0.119311906500871181D+08,
     *    0.404256246233400106D+08,
     *    0.488042433141206503D+08,
     *    0.661830130534069538D+08,
     *    0.573487910557239652D+08,
     *   -0.453174337495635748D+08,
     *   -0.164864480368007064D+09,
     *   -0.189428339851999760D+09,
     *   -0.119004368725137115D+09,
     *   -0.152607858517853856D+09,
     *   -0.156557709176206112D+09,
     *   -0.150423407980467558D+09,
     *   -0.221703207903122902D+08,
     *    0.222613157539944649D+08,
     *    0.118186476258069515D+09,
     *    0.157551100395653248D+09,
     *    0.401296348939333916D+09,
     *    0.498039294265785217D+08,
     *    0.251977823431519985D+09,
     *    0.247241493493816853D+09,
     *   -0.337501660288780451D+09,
     *    0.105340030851531029D+09,
     *    0.110962398429538965D+09,
     *   -0.118954416793255091D+09,
     *   -0.233726552770478725D+08,
     *    0.233118393319381118D+09,
     *   -0.304544461717576146D+09,
     *    0.365044082027730405D+09,
     *   -0.198512854833940148D+09/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_3=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X) 1A'' PES WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION VL0_4(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .     2.500D0,
     .     2.750D0,
     .     3.000D0,
     .     3.250D0,
     .     3.500D0,
     .     3.700D0,
     .     3.750D0,
     .     3.800D0,
     .     3.900D0,
     .     4.000D0,
     .     4.100D0,
     .     4.200D0,
     .     4.250D0,
     .     4.300D0,
     .     4.500D0,
     .     4.750D0,
     .     5.000D0,
     .     5.250D0,
     .     5.500D0,
     .     5.750D0,
     .     6.000D0,
     .     6.250D0,
     .     6.500D0,
     .     6.750D0,
     .     7.000D0,
     .     7.250D0,
     .     7.500D0,
     .     7.750D0,
     .     8.000D0,
     .     8.250D0,
     .     8.500D0,
     .     8.750D0,
     .     9.000D0,
     .     9.500D0,
     .    10.000D0,
     .    10.500D0,
     .    11.000D0,
     .    12.000D0,
     .    13.000D0,
     .    14.000D0,
     .    15.000D0/
      DATA VCOEF/
     *    0.929310837641872978D+10,
     *   -0.120226499897965584D+11,
     *    0.452286214353290272D+10,
     *   -0.199208097081481743D+10,
     *   -0.990678897782604694D+09,
     *    0.609390615861050606D+10,
     *   -0.134599332160586605D+11,
     *    0.963287310516804504D+10,
     *   -0.247992782563880253D+10,
     *    0.719534315337842226D+09,
     *   -0.168568470846518993D+09,
     *    0.892652895244128704D+08,
     *    0.772290758456313610D+07,
     *    0.371927179610063434D+08,
     *    0.757322782591098547D+08,
     *    0.136925104649899602D+09,
     *    0.186888335993679345D+09,
     *    0.110454992075304210D+09,
     *   -0.363840873287358880D+08,
     *   -0.837293787898122966D+08,
     *   -0.723821300468546152D+07,
     *   -0.161907169309473038D+08,
     *   -0.414561973988742828D+08,
     *   -0.130782482324384212D+09,
     *    0.160412444515299797D+07,
     *   -0.558759675305151939D+08,
     *    0.120103547533417940D+09,
     *    0.402228182308630943D+08,
     *    0.341848674664684534D+09,
     *    0.625908012154107094D+08,
     *    0.131706629573152065D+09,
     *    0.248027981471798658D+09,
     *   -0.344171744238330364D+09,
     *   -0.189563197491960049D+09,
     *    0.586439266307759523D+09,
     *   -0.669500195844280720D+09,
     *    0.139178059261439085D+09,
     *    0.377936745914343357D+09,
     *   -0.488195952532669544D+09,
     *    0.285011284286151528D+09,
     *   -0.454174977701806426D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_4=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X) 1A'' PES WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION VL0_5(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .     2.500D0,
     .     2.750D0,
     .     3.000D0,
     .     3.250D0,
     .     3.500D0,
     .     3.700D0,
     .     3.750D0,
     .     3.800D0,
     .     3.900D0,
     .     4.000D0,
     .     4.100D0,
     .     4.200D0,
     .     4.250D0,
     .     4.300D0,
     .     4.500D0,
     .     4.750D0,
     .     5.000D0,
     .     5.250D0,
     .     5.500D0,
     .     5.750D0,
     .     6.000D0,
     .     6.250D0,
     .     6.500D0,
     .     6.750D0,
     .     7.000D0,
     .     7.250D0,
     .     7.500D0,
     .     7.750D0,
     .     8.000D0,
     .     8.250D0,
     .     8.500D0,
     .     8.750D0,
     .     9.000D0,
     .     9.500D0,
     .    10.000D0,
     .    10.500D0,
     .    11.000D0,
     .    12.000D0,
     .    13.000D0,
     .    14.000D0,
     .    15.000D0/
      DATA VCOEF/
     *    0.904004852429895401D+10,
     *   -0.129878174443920021D+11,
     *    0.626574989020529366D+10,
     *   -0.389446867202347422D+10,
     *    0.181239748839883137D+10,
     *    0.483260857921675587D+10,
     *   -0.176098058487145157D+11,
     *    0.147818825265091152D+11,
     *   -0.389207622163572550D+10,
     *    0.112788072581564927D+10,
     *   -0.323303941176759720D+09,
     *    0.169768781715528607D+09,
     *   -0.752149621972426176D+08,
     *    0.377264112935075462D+08,
     *    0.460365913057785630D+08,
     *    0.114563667718437076D+09,
     *    0.207387527408397198D+09,
     *    0.160296131187769115D+09,
     *    0.185171967605620623D+08,
     *   -0.648578419227396846D+08,
     *    0.190121752754724026D+08,
     *    0.640268741663485765D+08,
     *    0.132015170125974417D+08,
     *   -0.262889039201059341D+08,
     *   -0.689186135770659447D+08,
     *    0.360155089229154587D+08,
     *   -0.124447656435595512D+09,
     *    0.133105976891431808D+09,
     *   -0.661278290199134350D+08,
     *    0.429734345670295238D+09,
     *   -0.966652353629703522D+08,
     *    0.292187115689948559D+09,
     *   -0.306496094332467079D+09,
     *   -0.189611673154765606D+09,
     *    0.186939580683349848D+09,
     *    0.376435412608709335D+08,
     *   -0.216490173479629278D+09,
     *    0.165834586991403818D+09,
     *   -0.316236505400052071D+08,
     *   -0.434429941211109161D+08,
     *    0.312583187300014496D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_5=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X) 1A'' PES WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION VL0_6(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .     2.500D0,
     .     2.750D0,
     .     3.000D0,
     .     3.250D0,
     .     3.500D0,
     .     3.700D0,
     .     3.750D0,
     .     3.800D0,
     .     3.900D0,
     .     4.000D0,
     .     4.100D0,
     .     4.200D0,
     .     4.250D0,
     .     4.300D0,
     .     4.500D0,
     .     4.750D0,
     .     5.000D0,
     .     5.250D0,
     .     5.500D0,
     .     5.750D0,
     .     6.000D0,
     .     6.250D0,
     .     6.500D0,
     .     6.750D0,
     .     7.000D0,
     .     7.250D0,
     .     7.500D0,
     .     7.750D0,
     .     8.000D0,
     .     8.250D0,
     .     8.500D0,
     .     8.750D0,
     .     9.000D0,
     .     9.500D0,
     .    10.000D0,
     .    10.500D0,
     .    11.000D0,
     .    12.000D0,
     .    13.000D0,
     .    14.000D0,
     .    15.000D0/
      DATA VCOEF/
     *    0.100759706728617020D+11,
     *   -0.160691616760629120D+11,
     *    0.782762612250686359D+10,
     *   -0.317852159413108349D+10,
     *    0.120618413542768025D+10,
     *    0.405198504740569258D+10,
     *   -0.136766678784151917D+11,
     *    0.114382158579803257D+11,
     *   -0.300626671341010237D+10,
     *    0.869218244948612928D+09,
     *   -0.239876617197723269D+09,
     *    0.141055063243318677D+09,
     *   -0.597043410202339888D+08,
     *    0.510982283457837105D+08,
     *    0.300447837950269580D+08,
     *    0.563325506348527074D+08,
     *    0.161084918740831077D+09,
     *    0.166161090814083427D+09,
     *    0.649990978194611371D+08,
     *   -0.525734270981119275D+08,
     *   -0.961139563135194778D+07,
     *    0.868178070121011138D+08,
     *   -0.521407431301712990D+06,
     *    0.946027535465812683D+08,
     *   -0.807226182871661186D+08,
     *    0.290848893214812279D+08,
     *   -0.627901847283496857D+08,
     *   -0.111039479500036240D+08,
     *   -0.222480107742202759D+09,
     *    0.379478024184921265D+09,
     *   -0.224299103278502941D+09,
     *    0.385584442919106007D+09,
     *   -0.896961690174956322D+08,
     *   -0.241048468192775726D+09,
     *    0.640238772032732964D+08,
     *    0.260834640291483402D+09,
     *   -0.282038811679893494D+09,
     *    0.327295814822568893D+08,
     *    0.368467775306375027D+08,
     *    0.678157379622745514D+08,
     *   -0.763065674926862717D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_6=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X) 1A'' PES WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION VL0_7(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .     2.500D0,
     .     2.750D0,
     .     3.000D0,
     .     3.250D0,
     .     3.500D0,
     .     3.700D0,
     .     3.750D0,
     .     3.800D0,
     .     3.900D0,
     .     4.000D0,
     .     4.100D0,
     .     4.200D0,
     .     4.250D0,
     .     4.300D0,
     .     4.500D0,
     .     4.750D0,
     .     5.000D0,
     .     5.250D0,
     .     5.500D0,
     .     5.750D0,
     .     6.000D0,
     .     6.250D0,
     .     6.500D0,
     .     6.750D0,
     .     7.000D0,
     .     7.250D0,
     .     7.500D0,
     .     7.750D0,
     .     8.000D0,
     .     8.250D0,
     .     8.500D0,
     .     8.750D0,
     .     9.000D0,
     .     9.500D0,
     .    10.000D0,
     .    10.500D0,
     .    11.000D0,
     .    12.000D0,
     .    13.000D0,
     .    14.000D0,
     .    15.000D0/
      DATA VCOEF/
     *    0.114500711247574978D+11,
     *   -0.197456295175286179D+11,
     *    0.102524484671034431D+11,
     *   -0.335813528957645416D+10,
     *    0.904583576154235840D+09,
     *   -0.214306781872536063D+09,
     *   -0.915403097080819488D+09,
     *    0.126758062523236871D+10,
     *   -0.321830419513134599D+09,
     *    0.144052968915025592D+09,
     *    0.270082561156666800D+08,
     *    0.735252448733273745D+08,
     *   -0.587072384511522651D+08,
     *    0.926510908190626502D+08,
     *    0.282529886791144609D+08,
     *    0.713991204979151487D+07,
     *    0.684263440580661595D+08,
     *    0.161439647380557179D+09,
     *    0.104136119696494371D+09,
     *   -0.166586929740441889D+08,
     *   -0.314278545559059978D+08,
     *    0.265942514878582954D+07,
     *    0.846562504611854553D+08,
     *   -0.620014889443945885D+07,
     *    0.954772365936779976D+07,
     *    0.252086053651237488D+08,
     *   -0.288761609910874367D+08,
     *   -0.514448738522992134D+08,
     *    0.662346102052121162D+08,
     *   -0.210917424398010254D+09,
     *    0.168142506747240543D+09,
     *    0.691137438617730141D+08,
     *    0.117876992959890366D+08,
     *   -0.821320313152875900D+08,
     *   -0.470315410271868706D+08,
     *    0.116788805798274994D+09,
     *   -0.252599590231947899D+08,
     *   -0.465616111265759468D+08,
     *    0.187342182750673294D+08,
     *    0.842550384613752365D+07,
     *    0.124426372140359879D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_7=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X) 1A'' PES WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION VL0_8(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .     2.500D0,
     .     2.750D0,
     .     3.000D0,
     .     3.250D0,
     .     3.500D0,
     .     3.700D0,
     .     3.750D0,
     .     3.800D0,
     .     3.900D0,
     .     4.000D0,
     .     4.100D0,
     .     4.200D0,
     .     4.250D0,
     .     4.300D0,
     .     4.500D0,
     .     4.750D0,
     .     5.000D0,
     .     5.250D0,
     .     5.500D0,
     .     5.750D0,
     .     6.000D0,
     .     6.250D0,
     .     6.500D0,
     .     6.750D0,
     .     7.000D0,
     .     7.250D0,
     .     7.500D0,
     .     7.750D0,
     .     8.000D0,
     .     8.250D0,
     .     8.500D0,
     .     8.750D0,
     .     9.000D0,
     .     9.500D0,
     .    10.000D0,
     .    10.500D0,
     .    11.000D0,
     .    12.000D0,
     .    13.000D0,
     .    14.000D0,
     .    15.000D0/
      DATA VCOEF/
     *    0.978447770752510452D+10,
     *   -0.182453761473609886D+11,
     *    0.121699425379775467D+11,
     *   -0.775953160069848633D+10,
     *    0.642140059812753677D+10,
     *   -0.104628958867475338D+11,
     *    0.126818375045029049D+11,
     *   -0.615147854767738152D+10,
     *    0.143547357792194223D+10,
     *   -0.353080228477867484D+09,
     *    0.156970756524426162D+09,
     *    0.713104080287766457D+08,
     *   -0.127575588596078277D+09,
     *    0.140410049866054773D+09,
     *    0.306113241797862053D+08,
     *   -0.119428610262348652D+08,
     *    0.102435302392476797D+08,
     *    0.101461409092207730D+09,
     *    0.101499999955560863D+09,
     *    0.620050592034339905D+07,
     *   -0.451440323106168509D+08,
     *   -0.220690703605403900D+08,
     *    0.102971905184883326D+09,
     *   -0.706303168862035275D+08,
     *    0.918474008209607601D+08,
     *   -0.553167942398841381D+08,
     *    0.962237547287940979D+06,
     *    0.371593543939793110D+08,
     *   -0.102130840270210266D+09,
     *   -0.962859213309240341D+07,
     *    0.219022483850303650D+09,
     *   -0.269167811936429501D+09,
     *   -0.345924615064330101D+08,
     *    0.601958861791147232D+09,
     *   -0.973328636085651398D+09,
     *    0.889451853970012665D+09,
     *   -0.483444472312618732D+09,
     *    0.272531090733235359D+09,
     *   -0.250393536983942986D+09,
     *    0.730207298939580917D+08,
     *    0.346101344531655312D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_8=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN Legendre POLYNOMIALS
C ANGULAR EXPANSION OF Kr-OH(X) 1A'' PES WITH
C RIGID OH(X) 
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO:MRCISD+Q/AVTZ-DK SC CORRECTED
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C 2012.02.29
C
C------SAMPLE PROGRAM CALLING PES FUNCTION-------
      FUNCTION VL0_9(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=41)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .     2.500D0,
     .     2.750D0,
     .     3.000D0,
     .     3.250D0,
     .     3.500D0,
     .     3.700D0,
     .     3.750D0,
     .     3.800D0,
     .     3.900D0,
     .     4.000D0,
     .     4.100D0,
     .     4.200D0,
     .     4.250D0,
     .     4.300D0,
     .     4.500D0,
     .     4.750D0,
     .     5.000D0,
     .     5.250D0,
     .     5.500D0,
     .     5.750D0,
     .     6.000D0,
     .     6.250D0,
     .     6.500D0,
     .     6.750D0,
     .     7.000D0,
     .     7.250D0,
     .     7.500D0,
     .     7.750D0,
     .     8.000D0,
     .     8.250D0,
     .     8.500D0,
     .     8.750D0,
     .     9.000D0,
     .     9.500D0,
     .    10.000D0,
     .    10.500D0,
     .    11.000D0,
     .    12.000D0,
     .    13.000D0,
     .    14.000D0,
     .    15.000D0/
      DATA VCOEF/
     *    0.754650577438139057D+10,
     *   -0.136328403744863052D+11,
     *    0.741382861798679638D+10,
     *   -0.296287623683514357D+10,
     *    0.225947659146189547D+10,
     *   -0.758455384539548111D+10,
     *    0.144783174801535511D+11,
     *   -0.977670210602620125D+10,
     *    0.248615610516424561D+10,
     *   -0.672210097211129665D+09,
     *    0.253505365142645121D+09,
     *   -0.592249594944287539D+08,
     *    0.479623096698665470D+08,
     *    0.479724691834328175D+08,
     *    0.510773633609924316D+08,
     *   -0.153173013747112751D+08,
     *   -0.226329315020355359D+08,
     *    0.591442927572842166D+08,
     *    0.917405867161158621D+08,
     *    0.146506748416869640D+08,
     *   -0.564460726706182659D+08,
     *    0.227051235812735558D+07,
     *   -0.183300596404278278D+07,
     *    0.609489500037976503D+08,
     *   -0.103423387479720116D+08,
     *   -0.144113086155002117D+08,
     *   -0.148106635797567368D+08,
     *    0.928027703063044548D+08,
     *   -0.167734455492631912D+09,
     *    0.181399512670068741D+09,
     *   -0.163020242442821980D+09,
     *   -0.693547759598622322D+08,
     *    0.179787767039282322D+09,
     *   -0.386918293533477783D+08,
     *    0.504320152460837364D+08,
     *   -0.127968815694409609D+09,
     *    0.682814708938987255D+08,
     *    0.333545053853590488D+08,
     *   -0.170995617069973946D+08,
     *   -0.604162517952039242D+08,
     *    0.535945490994033813D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_9=SUMA
       RETURN
       END

      FUNCTION VCI_H12(R,theta)
C*********************************
C System: Kr-OH(X2Pi-A2Sigma) Mixing 
C Diabatic H12CI Potential That 
C Couples 1A'(X-state) and 2A'(A-state) 
C Method: MRCISD+Q H12CI 
C Basis:aug-cc-pvtz-DK
C Variables:
C R in Bohr
C Theta n Degrees
C Output in cm-1
C Theta=0 for Kr---HO geom.
C********************************
C Authors:
C Dr Jacek Klos
C Department of Chemistry and Biochemistry
C University of Maryland
C College Park, MD 20742
C email:jklos@umd.edu
C 2012.02.28
C**********************************
C Needs link with LAPACK
C*********************************
      implicit double precision(a-h, o-z)
      dimension V0(8)
      dimension T(8)
      pi=dacos(-1.d0)
      call Vlpes(R,V0)
       do j=1,8
       T(j)=PLGNDR(j,1,DCOS(theta*pi/180.D0))
       enddo
       s=0.D0
       do i=1,8
       s=s+V0(i)*T(i)
       enddo
       RDAMP=7.D0
       damp=-0.5D0*(tanh(2.d0*(R-RDAMP))-1.D0)
       IF (R.GT.3.25D0) THEN
        VCI_H12=s*damp
       ELSE
       call Vlpes(3.25D0,V0)
       do j=1,8
       T(j)=PLGNDR(j,1,DCOS(theta*pi/180.D0))
       enddo
       s=0.D0
       do i=1,8
       s=s+V0(i)*T(i)
       enddo
       dampsr=-0.5D0*(tanh(10.D0*(R-3.25D0))-1.D0)
       VCI_H12=s*R/3D0*dampsr+(1.D0-dampsr)*s
       ENDIF
       return
       end


      subroutine Vlpes(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(8,8),ipvt(8),work(100)
      dimension V0(8),theta(8)
      LWORK=100
      pi=dacos(-1.d0)
      theta(1)=20.D0
      theta(2)=40.D0
      theta(3)=60.D0
      theta(4)=90.D0
      theta(5)=100.D0
      theta(6)=120.D0
      theta(7)=140.D0
      theta(8)=160.D0
      do i=1,8
       do j=1,8
       T(i,j)=PLGNDR(j,1,DCOS(theta(i)*pi/180.D0)) 
       enddo
      enddo
      do k=1,8
      V0(k)=VPESN(R,k)
      enddo
      call dgesv(8,1,T,8,ipvt,V0,8,info)
c      call dgels('N',11,9,1,T,11,V0,11,WORK,LWORK,info)
      return
      end


      FUNCTION VPESN(R,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(15),XX1(15),XX2(15),XX3(15),
     * XX4(15),XX5(15),XX6(15),XX7(15),
     * XX8(15)
C theta=20
      DATA XX1/
     *   0.278425374601325215D+01,
     *  -0.964989991984151629D+01,
     *   0.111321932171445782D+07,
     *  -0.250378401839806407D+06,
     *  -0.503634740456341897D+05,
     *   0.100726611142919592D+05,
     *   0.487573094216573918D+04,
     *   0.862267003447068987D+02,
     *  -0.218464971074001795D+03,
     *  -0.127967574552116803D+01,
     *   0.501359077698944766D+01,
     *   0.482719905126656704D+04,
     *  -0.297569680437203944D+09,
     *   0.143670599314466834D+10,
     *   0.964306617031438637D+10/
C theta=40
      DATA XX2/
     *   0.293165071127555343D+01,
     *  -0.937968524879164711D+01,
     *   0.105647923412486655D+07,
     *  -0.255007550874013046D+06,
     *  -0.493321806383625590D+05,
     *   0.102159722402185489D+05,
     *   0.468842539111288806D+04,
     *   0.153497247629266163D+02,
     *  -0.221243284823940854D+03,
     *   0.606929928028193366D+01,
     *   0.317692776893101492D+01,
     *   0.482719905126656704D+04,
     *  -0.686179539641419833D+05,
     *  -0.835543424035223842D+09,
     *   0.705578493742819691D+10/
C theta=60
      DATA XX3/       
     *   0.264970445541508060D+01,
     *  -0.895017124697448629D+01,
     *   0.110793066416331963D+07,
     *  -0.231575719713749451D+06,
     *  -0.445841153366188883D+05,
     *   0.913944972816614063D+04,
     *   0.428275573631500538D+04,
     *   0.161363723464199438D+03,
     *  -0.153300352633726902D+03,
     *  -0.164320283786271837D+02,
     *   0.623272777152817348D+01,
     *   0.482719905126656704D+04,
     *  -0.553044570367415309D+09,
     *   0.489750513187515354D+10,
     *   0.786854716780362427D+08/
C theta=90
      DATA XX4/       
     *   0.290453423918802134D+01,
     *  -0.109183084245323059D+02,
     *   0.124354159126578667D+07,
     *  -0.275411741325174982D+06,
     *  -0.585087469163042115D+05,
     *   0.111696798548448296D+05,
     *   0.560602484378459212D+04,
     *   0.118154153798641417D+03,
     *  -0.259563013301619492D+03,
     *  -0.302409873222872938D+00,
     *   0.600089201432861330D+01,
     *   0.482719905126656704D+04,
     *  -0.964120732970950305D+08,
     *  -0.453483592952518272D+10,
     *   0.456181407003679581D+11/
C theta=100
      DATA XX5/       
     *   0.282910330769598417D+01,
     *  -0.105961502533886716D+02,
     *   0.123645399602043629D+07,
     *  -0.266771546856021858D+06,
     *  -0.569350719581452067D+05,
     *   0.106114238639096493D+05,
     *   0.544415110365345390D+04,
     *   0.144942750197369350D+03,
     *  -0.244165062610752784D+03,
     *  -0.309364779934446910D+01,
     *   0.608604061799635776D+01,
     *   0.482719905126656704D+04,
     *  -0.449916998545581460D+09,
     *   0.218425574126848578D+07,
     *   0.332319031756438675D+11/
C theta=120
      DATA XX6/
     *   0.257732055609478916D+01,
     *  -0.896335141000426106D+01,
     *   0.115378366535044950D+07,
     *  -0.186255544414435455D+06,
     *  -0.416928521548271456D+05,
     *   0.762319696659298188D+04,
     *   0.485238075589710297D+04,
     *  -0.183904603831572082D+02,
     *  -0.131987360532601201D+02,
     *  -0.438636634171842559D+02,
     *   0.878034891353057034D+01,
     *   0.482719905126656704D+04,
     *  -0.132907584815306568D+10,
     *   0.130316385724559822D+11,
     *  -0.100171011411218586D+11/
C theta=140
      DATA XX7/
     *   0.494381911330640023D+01,
     *  -0.169693692997838959D+02,
     *   0.989438613150500809D+06,
     *  -0.468584221245118475D+06,
     *  -0.483536843874391925D+05,
     *   0.244607480687114694D+05,
     *   0.733675421452483079D+04,
     *  -0.507238640034652519D+03,
     *  -0.506558713702059663D+03,
     *   0.386755294333027317D-03,
     *   0.109819691655251095D+02,
     *   0.482719905126656704D+04,
     *  -0.244026707530737631D+08,
     *   0.479560423110483229D+09,
     *  -0.243103135204228401D+10/
C theta=160
      DATA XX8/
     *   0.636804559982606921D+01,
     *  -0.235703859663975308D+02,
     *   0.957622999410666758D+06,
     *  -0.445523946353751933D+06,
     *  -0.479872879253503634D+05,
     *   0.229318021966122469D+05,
     *   0.706143341432095713D+04,
     *  -0.437348301390015422D+03,
     *  -0.467649834832280646D+03,
     *   0.621768085165382889D+01,
     *   0.849752679277828982D+01,
     *   0.482719905126656704D+04,
     *   0.898188024302440695D+07,
     *   0.987784443895764649D+07,
     *  -0.707388799860129356D+09/
      IF(N.EQ.1) THEN
      DO I=1,15
      XX(I)=XX1(I)
      ENDDO
      ENDIF
      IF(N.EQ.2) THEN
      DO I=1,15
      XX(I)=XX2(I)
      ENDDO
      ENDIF
      IF(N.EQ.3) THEN
      DO I=1,15
      XX(I)=XX3(I)
      ENDDO
      ENDIF
      IF(N.EQ.4) THEN
      DO I=1,15
      XX(I)=XX4(I)
      ENDDO
      ENDIF
      IF(N.EQ.5) THEN
      DO I=1,15
      XX(I)=XX5(I)
      ENDDO
      ENDIF
      IF(N.EQ.6) THEN
      DO I=1,15
      XX(I)=XX6(I)
      ENDDO
      ENDIF
      IF(N.EQ.7) THEN
      DO I=1,15
      XX(I)=XX7(I)
      ENDDO
      ENDIF
      IF(N.EQ.8) THEN
      DO I=1,15
      XX(I)=XX8(I)
      ENDDO
      ENDIF


      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**6)+XX(14)/(R**7)+XX(15)/(R**8)
       VPESN = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END


C------RKHS ROUTINES

C------------------------RKHS---------------------------------
C Fortran  Deck with distance-like and angular-like 
C kernels for RKHS interpolation method
C****************************************************
C* By J. KLOS  12/05/2007   University of Maryland  *
C* REVISON: 1.0                                     * 
C* jasztol@gmail.com                                *
C****************************************************
C-------------------------------------------------------------
C FUNCTIONS:
C
C RKHS_DK(X,Y,N,M): Reproducing kernel for  distance-like (DK) 
C                   variables, see Ho, Rabitz Eq(17) of
C                   J. Chem. Phys. 104, 2584 (1996) 
C                   Range = [0, inf]
C
C  INPUT: X,Y-DOUBLE PRECISION
C          N:INTEGER:ORDER OF SMOOTHNESS OF THE FUNCTION
C          M>=0:INTEGER:RECIPROCAL POWER OF THE WEIGHT W(X)=X**(-M)
C  OUTPUT: DOUBLE PRECISION
C -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -            
C RKHS_AK(X,Y,N): Reproducing kernel for  angular-like (AK) 
C                   variables, see Ho, Rabitz Eq(23) of
C                   J. Chem. Phys. 104, 2584 (1996) 
C                    Range = [0, 1]
C    
C  INPUT: X,Y-DOUBLE PRECISION
C          N:INTEGER:ORDER OF SMOOTHNESS OF THE FUNCTION
C  OUPUT: DOUBLE PRECISION
C  HINT: SCALE ANGLES TO [0,1] INTEGRAL FOR ANGLES, FOR EXAMPLE: 
C          X=(1.d0-cosd(ANGLE))/2.d0
C -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C RKHS_OP(X,Y,N,M):Reproducing kernel for angular like 
C                    variables with a use of orthogonal 
C                    polynomials Plm(L,M,X) where X is 
C                    K=SUM_L=M^N PLM(X)PLM(X')
C                    variable transformed from angles as
C                    x=cos(theta)
C-------------------------------------------------------------
       DOUBLE PRECISION FUNCTION RKHS_DKN2(X,Y,M)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IF (X .LT. Y) THEN
           XB=X
           XU=Y
       ELSE
           XB=Y
           XU=X
       ENDIF

C       M=ABS(M)! FORCE M>=0
       NSQR=4
       XUPOWM=XU**(-(M+1))

C      TERM WITH BETA FUNCTION
C       CALL BETA(DFLOAT(M+1),DFLOAT(N),BETATERM)
C      N=2 case direct expression
       BETATERM=1.D0/DFLOAT((M+1)*(M+2))

C      GET TERM WITH 2F1 HYPERGEOMETRIC  GAUSS FUNCTION
C       CALL HYGFX(DFLOAT(-N+1),DFLOAT(M+1),DFLOAT(N+M+1),XB/XU,HF)
C      DIRECT EXPRESSION FOR 2F1 WITH N=2 case
       HF=1.D0-(XB/XU)*(DFLOAT(M+1)/DFLOAT(M+3))
       RKHS_DKN2=NSQR*XUPOWM*BETATERM*HF
       RETURN
       END

       DOUBLE PRECISION FUNCTION RKHS_DK(X,Y,N,M)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IF (X .LT. Y) THEN
           XB=X
           XU=Y
       ELSE
           XB=Y
           XU=X
       ENDIF

C       M=ABS(M)! FORCE M>=0
       NSQR=N*N
       XUPOWM=XU**(-(M+1))

C      TERM WITH BETA FUNCTION
       CALL BETA(DFLOAT(M+1),DFLOAT(N),BETATERM)

C      GET TERM WITH 2F1 HYPERGEOMETRIC  GAUSS FUNCTION
       CALL HYGFX(DFLOAT(-N+1),DFLOAT(M+1),DFLOAT(N+M+1),XB/XU,HF)
       RKHS_DK=NSQR*XUPOWM*BETATERM*HF
       RETURN
       END

C-------------------------------------------------------------
       DOUBLE PRECISION FUNCTION RKHS_EXPR(X,Y,XE,N)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)
       
       SUMA=ZERO
       Z1=(X-XE)/XE
       Z2=(Y-XE)/XE
       DO I=0,N
       SUMA=SUMA+(Z1**I)*(Z2**I)
       ENDDO

       RKHS_EXPR=SUMA
       RETURN
       END


C-------------------------------------------------------------
       DOUBLE PRECISION FUNCTION RKHS_AK(X,Y,N)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)
       
       IF ((X.EQ.ZERO).OR.(Y.EQ.ZERO)) THEN
          RKHS_AK=ONE
          RETURN
       ENDIF
       IF (X .LT. Y) THEN
           XB=X
           XU=Y
       ELSE
           XB=Y
           XU=X
       ENDIF

       SUMA=ZERO
C      GET POLYNOMIAL TERM
       DO I=0,N-1
       SUMA=SUMA+(XU**I)*(XB**I)
       END DO

C      GET TERM WITH 2F1 HYPERGEOMETRIC  GAUSS FUNCTION
       CALL HYGFX(DFLOAT(1),DFLOAT(-N+1),DFLOAT(N+1),XB/XU,HF)
       
       RKHS_AK=SUMA+N*(XB**N)*(XU**(N-1))*HF
       RETURN
       END

       DOUBLE PRECISION FUNCTION RKHS_AKN2(X,Y)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)
       
       IF ((X.EQ.ZERO).OR.(Y.EQ.ZERO)) THEN
          RKHS_AKN2=ONE
          RETURN
       ENDIF 
       IF (X .LT. Y) THEN
           XB=X
           XU=Y
       ELSE
           XB=Y
           XU=X
       ENDIF

C      GET POLYNOMIAL TERM
       SUMA=1.D0+XU*XB

C      GET TERM WITH 2F1 HYPERGEOMETRIC  GAUSS FUNCTION
C       CALL HYGFX(DFLOAT(1),DFLOAT(-N+1),DFLOAT(N+1),XB/XU,HF)
       HF=1.D0-(XB/XU)/3.D0
       
       RKHS_AKN2=SUMA+2.D0*(XB**2)*XU*HF
       RETURN
       END


       DOUBLE PRECISION FUNCTION RKHS_OP(X,Y,N,M)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)

       SUMA=ZERO
C      GET ORTHOGONAL POLYNOMIAL TERM
       DO I=ABS(M),N
       SUMA=SUMA+PLGNDR(I,M,X)*PLGNDR(I,M,Y)
       END DO
       RKHS_OP=SUMA
       RETURN
       END

       DOUBLE PRECISION FUNCTION RKHS_OP2L(X,Y,N,M)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)
 
       SUMA=ZERO
C      GET ORTHOGONAL POLYNOMIAL TERM
       DO I=ABS(M),N,2
       SUMA=SUMA+PLGNDR(I,M,X)*PLGNDR(I,M,Y)
       END DO
       RKHS_OP2L=SUMA
       RETURN
       END


            



        SUBROUTINE BETA(P,Q,BT)
C
C       ==========================================
C       Purpose: Compute the beta function B(p,q)
C       Input :  p  --- Parameter  ( p > 0 )
C                q  --- Parameter  ( q > 0 )
C       Output:  BT --- B(p,q)
C       Routine called: GAMMA for computing �(x)
C       ==========================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        CALL GAMMA(P,GP)
        CALL GAMMA(Q,GQ)
        PPQ=P+Q
        CALL GAMMA(PPQ,GPQ)
        BT=GP*GQ/GPQ
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function (x)
C       Input :  x  --- Argument of ( x )
C                       ( x is not equal to 0,-1,-2,...)
C       Output:  GA --- (x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END



        SUBROUTINE HYGFX(A,B,C,X,HF)
C
C       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C       Routines called:
C            (1) GAMMA for computing gamma function
C            (2) PSI_AAA for computing psi function
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        L0=C.EQ.INT(C).AND.C.LT.0.0
        L1=1.0D0-X.LT.1.0D-15.AND.C-A-B.LE.0.0
        L2=A.EQ.INT(A).AND.A.LT.0.0
        L3=B.EQ.INT(B).AND.B.LT.0.0
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0
        IF (L0.OR.L1) THEN
           WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        EPS=1.0D-15
        IF (X.GT.0.95) EPS=1.0D-8
        IF (X.EQ.0.0.OR.A.EQ.0.0.OR.B.EQ.0.0) THEN
           HF=1.0D0
           RETURN
        ELSE IF (1.0D0-X.EQ.EPS.AND.C-A-B.GT.0.0) THEN
           CALL GAMMA(C,GC)
           CALL GAMMA(C-A-B,GCAB)
           CALL GAMMA(C-A,GCA)
           CALL GAMMA(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (1.0D0+X.LE.EPS.AND.DABS(C-A+B-1.0).LE.EPS) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMA(C,G1)
           CALL GAMMA(1.0D0+A/2.0-B,G2)
           CALL GAMMA(0.5D0+0.5*A,G3)
           HF=G0*G1/(G2*G3)
           RETURN
        ELSE IF (L2.OR.L3) THEN
           IF (L2) NM=INT(ABS(A))
           IF (L3) NM=INT(ABS(B))
           HF=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
10            HF=HF+R
           RETURN
        ELSE IF (L4.OR.L5) THEN
           IF (L4) NM=INT(ABS(C-A))
           IF (L5) NM=INT(ABS(C-B))
           HF=1.0D0
           R=1.0D0
           DO 15 K=1,NM
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X
15            HF=HF+R
           HF=(1.0D0-X)**(C-A-B)*HF
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0D0) THEN
           X=X/(X-1.0D0)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        IF (X.GE.0.75D0) THEN
           GM=0.0D0
           IF (DABS(C-A-B-INT(C-A-B)).LT.1.0D-15) THEN
              M=INT(C-A-B)
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(A+M,GAM)
              CALL GAMMA(B+M,GBM)
              CALL PSI_AAA(A,PA)
              CALL PSI_AAA(B,PB)
              IF (M.NE.0) GM=1.0D0
              DO 30 J=1,ABS(M)-1
30               GM=GM*J
              RM=1.0D0
              DO 35 J=1,ABS(M)
35               RM=RM*J
              F0=1.0D0
              R0=1.0D0
              R1=1.0D0
              SP0=0.D0
              SP=0.0D0
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=-GC*(X-1.0D0)**M/(GA*GB*RM)
                 DO 40 K=1,M-1
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X)
40                  F0=F0+R0
                 DO 45 K=1,M
45                  SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K
                 F1=PA+PB+SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 55 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 50 J=1,M
50                     SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/
     &                    (B+J+K-1.0)
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 60
55                  HW=F1
60               HF=F0*C0+F1*C1
              ELSE IF (M.LT.0) THEN
                 M=-M
                 C0=GM*GC/(GA*GB*(1.0D0-X)**M)
                 C1=-(-1)**M*GC/(GAM*GBM*RM)
                 DO 65 K=1,M-1
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X)
65                  F0=F0+R0
                 DO 70 K=1,M
70                  SP0=SP0+1.0D0/K
                 F1=PA+PB-SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 80 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 75 J=1,M
75                     SM=SM+1.0D0/(J+K)
                    RP=PA+PB+2.0D0*EL+SP-SM+DLOG(1.0D0-X)
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 85
80                  HW=F1
85               HF=F0*C0+F1*C1
              ENDIF
           ELSE
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(C-A,GCA)
              CALL GAMMA(C-B,GCB)
              CALL GAMMA(C-A-B,GCAB)
              CALL GAMMA(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)
              HF=0.0D0
              R0=C0
              R1=C1
              DO 90 K=1,250
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X)
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))
     &              *(1.0-X)
                 HF=HF+R0+R1
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 95
90               HW=HF
95            HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0D0
           IF (C.GT.A.AND.C.LT.2.0D0*A.AND.
     &         C.GT.B.AND.C.LT.2.0D0*B) THEN
              A0=(1.0D0-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0D0
           R=1.0D0
           DO 100 K=1,250
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
              HF=HF+R
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 105
100           HW=HF
105        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0D0) THEN
           X=X1
           C0=1.0D0/(1.0D0-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.120) WRITE(*,115)
115     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END



        SUBROUTINE PSI_AAA(X,PS)
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
