* ----------------------------------------------------------------------
      subroutine basgpi (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, iom, iohold, nvhold, nlv,
     :                  rcut, jtot, flaghf, flagsu, csflag, clist,
     :                  bastst, ihomo, nu, numin, jlpar, n, nmax, ntop)
* ----------------------------------------------------------------------
*  subroutine to determine angular coupling potential for collision of a
*  molecule in either a 2pi electronic state in an intermediate coupling
*  basis or in a 2sigma electronic state in a hund's case (a) coupling
*  with a structureless atom or with an uncorrugated surface
*  authors:  original version by millard alexander and didier lemoine
*            rewritten and extended to vibration and sigma-pi coupling
*            by h.-j. werner
*  current revision date:  14-nov-1994 by hjw
*  slightly modified    :  20-mar-1992 by ab
*  rcut bug fixed:         31-mar-1992 by mha
*  no rcut for boundc   :  13-may-1997 by mha
*  increased format for NCHANNELS= :  7-apr-2003 by mha

*  --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains rotational quantum numbers for each
*              channel
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains symmetry index of each channel
*  note that we have adopted the following convention for the symmetry
*  index "is" so that on return the doublet pi molecular levels can be
*  uniquely identified by the two labels "j" and "is":
*           for omega = 1/2   is = (100+ivp)*eps
*           for omega = 3/2   is = (200+ivp)*eps
*           for sigma         is = (300+ivs)*eps
*        where  eps = +1 or -1 is the "true" case (a) symmetry index
*        end ivs, ivp are the vibrational quantum numbers
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level
*    ehold:    on return contains energy in hartrees of each rotational
*              level
*    ishold:   on return contains symmetry index of each energetically
*              distinct level, similar to is
*    iohold:   on return, contains nominal omega values for each level
*    nvhold:   on return, contains vibr. quantum numbers for each level
*    nlv:      on return, contains vibr. quantum numbers for each chan.
*    nlvp:     scratch array which holds perturbing vibrational levels
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    iom:      on return contains nominal omega values for each channel
*              for 2sigma states iom is assigned the value of 2
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              note that this is ignored in molecule-surface collisions!
*    jtot:     total angular momentum
*              in cc calculation jtot+1/2 is the total angular momentum
*              in cs calculation jtot is the l-bar quantum number
*    flaghf:   if .true., then system has half-integer spin
*              if .false., then system has integer spin
*    flagsu:   if .true., then molecule-surface collisons
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true., then homonuclear molecule
*              only the s or a levels will be included depending on the
*              value of the parameter isa in common /cosysi/ (see below)
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              only those channels are included for which
*                  (-1)**(parity+l-jtot)=jlpar
*              where parity is +1 for e levels and -1 for f levels
*              with the standard definition that e levels are those for
*              which eps*(-1)**(j-1/2-s) = 1
*              and f levels are those for which eps(-1)**(j-1/2-s) = -1
*              here s=1 for sigma-minus and s=0 otherwise
*              in cs calculation jlpar is set equal to 1 in calling program
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*    note!!!   if flaghf = .true., then the true values of the rotational
*    quantum numbers, the total angular momentum, and the coupled-states
*    projection index are equal to the values stored in j, jtot, and nu
*    plus 1/2
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*              rpar: spectroscopic constants for sigma and pi states,
*              defined in subroutine hisysgpi and used in subroutine hsgpi
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential.  This can be 1, 2, or 4 here and
*              must be consistent with the pot routine.
*    nvmins:   lowest vibrational level for sigma state
*    nvmaxs:   higest vibrational level for sigma state
*    nvminp:   lowest vibrational level for pi state
*    nvmaxp:   higest vibrational level for pi state
*    jmin:     the minimum rotational angular momenta for each 2pi channel
*    jmax:     the maximum rotational angular momenta for each 2pi channel
*              in each spin-orbit manifold with convention
*              omega .le. j .le. jmax+0.5
*    nmin:     the minimum case (b) rotational angular momenta for the
*              2sigma state
*    nmax:     the maximum case (b) rotational angular momenta for the
*              2sigma state
*    jmin,jmax,nmin,nmax are defined separately for each vibartional level
*              Each vibrational Pi level ivp may be perturbed by one sigma
*              vibrational level ivs. This level is given by ivs=ipert(ivp)
*              (see hisysgpi)
*    igupi:    permutation inversion symmetry of 2pi electronic state
*              igupi=1 for gerade states, igupi=-1 for ungerade states
*              for heteronuclear molecules igu should be +1
*    igusg:    permutation inversion symmetry of 2sig electronic state
*              igusg=1 for gerade states, igusg=-1 for ungerade states
*              for heteronuclear molecules igu should be +1
*    nparpi:   number of 2pi symmetry doublets included (npar=2 will ensure
*              both lambda doublets; npar=1 just eps=1 doublets;
*              npar=-1, just eps=-1 doublets
*    nparsg:   number of spin doublets in 2sigma state included (npar=2 will
*              ensure both spin doublets)
*    isymsg:   if isym=+1, then the electronic symmetry of the 2sigma state
*              is sigma-plus
*              if isym=-1, then the electronic symmetry is sigma-minus
*    isa:      s/a symmetry index, if the molecule is homonuclear (ihomo=t)
*              then, if isa=+1 then only the s-levels (both for sigma and
*              pi) are included in the basis, if isa=-1, then only the
*              a-levels are included
*    isg:      if isg=1 and ipi=0 then 2sig + atom scattering
*    ipi:      if ipi=1 and isg=0 then 2pi + atom scattering
*              if isg=1 and ipi=1 then 2pi-2sig + atom scattering
*  variables in common block /cobsp2/
*    nvt:      Number of vibrational blocks for each term. All of these
*              use same lammin, lammax, mproj. These numbers as well
*              as the corresponding ivcol, ivrow (see below) should be
*              set in loapot and must be consistent with the pot routine
*              Otherwise unrecognized chaos!!!
*              Each block corresponds to a pair of vibrational quantum
*              numbers (ivrow,ivcol), which must correspond to lower triangle
*              of potential matrix (i.e., ivrow.ge.ivcol for
*              iterm=1 (Vsig), iterm=2 (Vpi), iterm=4 (V2) and
*              for iterm=3 (V1) ivrow=ivs, ivcol=ivp)
*    ivrow(ivb,iterm): row vibrational state for vibrational block ivb
*                      in term iterm
*    ivcol(ivb,iterm): column vibrational state for vibrational block ivb
*                      in term iterm
*  variables in common block /cosgpi/
*    nvibs:    number of vibrational terms for sigma state in input file
*    ivibs:    vibrational quantum number for each of these
*    nvibp:    number of vibrational terms for pi state in input file
*    ivibp:    vibrational quantum number for each of these
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units
*               (mass of electron = 1)
*  variable in common block /conlam/
*    nlam:      the number of angular coupling terms actually used
*    nlammx:    the maximum number of angular coupling terms
*  variables in common block /cosgpi/
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               only nonzero elements are stored
*  variable in common block /coiv2/
*   iv2:        matrix address of v2 matrix for each non-zero element
*  variable in common block /coconv/
*     econv:    conversion factor from cm-1 to hartrees
*     xmconv:   converson factor from amu to atomic units
*  subroutines called:
*   vlsgpi:    returns angular coupling coefficient for particular
*              choice of channel index
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical csflag, clist, flaghf, flagsu, ihomo, bastst
      include "common/parbas"
      include "common/parbasl"
      common /cosysr/ isrcod, junkr, rpar(10)
      common /coipar/ iiipar(9), iprint
      common /cosysi/ nscode, isicod, ispar(10)
*  these parameters must be the same as in hisysgpi
      common /covib/ nvibs,ivibs(maxvib),nvibp,ivibp(maxvib)
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /cotq1/ vec(3,3,1)
      common /coisc1/ ivec(1)
      common /coisc2/ ivhold(1)
      common /coisc3/ nlvp(1)
      common /coisc4/ nvphol(1)
      common /conlam/ nlam, nlammx,lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      dimension e(3,3), ieps(2), iepp(2), iomc(4), iomr(4), eig(3)
      dimension jhold(1), ishold(1), iohold(1), nvhold(1), ehold(1),
     :          j(1), is(1), iom(1), nlv(1), l(1)
      parameter (nvmax=20)
      dimension ipvs(0:nvmax),ipvp(0:nvmax)
      data iomc/2,1,1,1/
      data iomr/2,1,2,1/
* this determines which eps level is first in channel list (arbitrary)
      data ieps / -1, 1 /
      data iepp / 1, -1 /
      data izero, ione, itwo
     :   / 0,     1,    2 /
* nprsg is number of real parameters per vib state for sigma state
* nprpi is number of real parameters per vib state for pi state
      data nprsg,nprpi/7,14/
* recover system parameters
      nterm=ispar(1)
      isg=ispar(2)
      ipi=ispar(3)
      if(nterm.eq.1) then
        isg=1
        ipi=0
      else if(nterm.eq.2) then
        isg=0
        ipi=1
      end if
      isgpi=ispar(4)
* nprpi is number of real system parameters for pi state
* this is 12 for just pi, 14 if sg/pi coupling, in which case alpha and
* beta are added for mixed levels defined by ipert
* isgpi.lt.0: mixing of pi(3/2) and pi(1/2) suppressed
      if (isgpi .le. 0) then
        nprpi=12
      else if (isgpi .eq. 1) then
        nprpi=14
      end if
      isa=ispar(5)
      if(.not.ihomo) isa=1
      iofi=5
      iofr=0
      if(isg.ne.0) then
* paramters for sigma state
        igusg=ispar(iofi+1)
        if(.not.ihomo) igusg=1
        nparsg=ispar(iofi+2)
        isymsg=ispar(iofi+3)
        nvmins=ispar(iofi+4)
        nvmaxs=ispar(iofi+5)
*  check that requested vib levels have been defined
        if(nvmaxs.gt.nvmax) stop 'nvmax'
        do 12 ivs=nvmins,nvmaxs
        do 10 i=1,nvibs
10      if(ivibs(i).eq.ivs) goto 12
          write(6,11) ivs
11        format(/' NO PARAMETERS DEFINED FOR VIB STATE',i3)
          call exit
12      ipvs(ivs)=i
        isplus = igusg * isymsg * isa
        nmins=iofi+4
        nmaxs=iofi+5
        noffs=-nprsg
        iofi=iofi+5+2*nvibs
        iofr=nprsg*nvibs
* these are offsets:
*   nmin(iv)=ispar(nmins+2*iv)
*   nmax(iv)=ispar(nmaxs+2*iv)
*   bsg(iv)=rpar(noffs+nprsg*iv+1)
*   dsg(iv)=rpar(noffs+nprsg*iv+2)  etc.
      end if
      if(ipi.ne.0) then
* parameters for pi state
        igupi=ispar(iofi+1)
        if(.not.ihomo) igupi=1
        nparpi=ispar(iofi+2)
        nvminp=ispar(iofi+3)
        nvmaxp=ispar(iofi+4)
*  check that requested vib levels have been defined
        if(nvmaxp.gt.nvmax) stop 'nvmax'
        do 22 ivp=nvminp,nvmaxp
        do 20 i=1,nvibp
20      if(ivibp(i).eq.ivp) goto 22
          write(6,11) ivp
          call exit
22      ipvp(ivp)=i
        jminp=iofi+2
        jmaxp=iofi+3
        nvoff=iofi+4
        noffp=iofr-nprpi
        iofi=iofi+4+3*nvibp
        iofr=iofr+nprpi*nvibp
* these are offsets:
*   jmin(iv)=ispar(jminp+3*iv)
*   jmax(iv)=ispar(jmaxp+3*iv)
*   bpi(iv)=rpar(noffp+nprpi*iv+1)
*   dpi(iv)=rpar(noffp+nprpi*iv+2)  etc.
      end if
*
      zero = 0.d0
      tzero=1.d-12
      one = 1.d0
      two = 2.d0
      four = 4.d0
      half = 0.5d0
      quart = 0.25d0
      xjtot = jtot + half
      xnu = nu + half
*  check for consistency in the values of flaghf and csflag
      if (.not.flaghf) then
        write (6, 2)
        write (9, 2)
2       format (' *** FLAGHF = .FALSE. FOR DOUBLET SYSTEM; ABORT ***')
        call exit
      end if
      if (flagsu .and. .not. csflag) then
        write (6, 3)
        write (9, 3)
3       format
     :   ('  *** CSFLAG = .FALSE. FOR SURFACE CALCULATION; ABORT ***')
        call exit
      end if
      istep = 1
      if (ihomo) istep = 2
      do 6  i = 1, nterm
        if (ihomo) then
          if (mod(lammax(i)-lammin(i),2) .ne. 0) then
            write (6, 4) i, lammin(i), lammax(i)
            write (9, 4) i, lammin(i), lammax(i)
4           format (' *** IHOMO=T BUT ODD NO. OF TERMS FOR I=', i2/,
     :              '     LAMMIN=', i2, ' LAMMAX=', i2, '; ABORT ***')
            call exit
          end if
        end if
6     continue
      if (clist) then
        if (flagsu) then
          write (9, 7) rmu * xmconv, ered * econv, jtot, xnu
          write (6, 7) rmu * xmconv, ered * econv, jtot, xnu
7         format (/,' ** 2SIG + 2PI INT. COUPLING UNCORR. SURFACE **',
     :            /,'    RMU=', f9.4,'  E=', f7.2, '  LBAR=', i5,
     :              '  NU=', f5.1)
        else
          if (ipi.eq.0.0) then
            write (9, 25) igusg,isymsg,isa
            write (6, 25) igusg,isymsg,isa
          else if (isg.eq.0) then
            write (9, 26) igupi,isa
            write (6, 26) igupi,isa
          else
            write (9, 27) igusg,isymsg,igupi,isa
            write (6, 27) igusg,isymsg,igupi,isa
          end if
25        format (/,' **  2SIG CASE (A) COUPLING **',
     :              '   g/u=', i2,'  +/- = ', i2,'   s/a = ',i2)
26        format (/,' **  2PI INT. COUPLING **',
     :              '   g/u=', i2,'   s/a = ',i2)
27        format (/,' **  2SIG + 2PI INT. COUPLING **'/
     :              '     2SIG:   g/u=', i2,'  +/- = ', i2/
     :              '     2PI:    g/u=', i2,'  s/a = ',i2)
          if (csflag) then
            write (9, 30) rmu * xmconv, ered * econv, jtot, xnu
            write (6, 30) rmu * xmconv, ered * econv, jtot, xnu
30          format (/,' ** CS **    RMU=', f9.4,'  E=', f8.2,
     :              '  LBAR=', i5, '  NU=', f5.1)
          else
            write (9, 31) rmu * xmconv, ered * econv, xjtot, jlpar
            write (6, 31) rmu * xmconv, ered * econv, xjtot, jlpar
31          format (/,' ** CC **    RMU=', f9.4,'  E=', f8.2,
     :              '  JTOT=', f5.1, ' JLPAR=', i2)
          end if
          if(ipi.ne.0) then
            write(9,32)
            write(6,32)
32          format(/' 2PI-PARAMETERS:'/
     :        '  V  JMIN JMAX',5x,'BROT',8x,'DROT',8x,'HROT',9x,'ASO',
     :           9x,'P',/,20x,'PD',10x,'Q',11x,'QD',11x,'QH',7x,'EVIB',
     :           /,20x,'AD',9x,'GPI')
            do 34 ivp=nvminp,nvmaxp
              i=ipvp(ivp)
              jmin=ispar(jminp+3*i)
              jmax=ispar(jmaxp+3*i)
              write (9, 33) ivp,jmin,jmax,(rpar(noffp+i*nprpi+k),k=1,12)
              write (6, 33) ivp,jmin,jmax,(rpar(noffp+i*nprpi+k),k=1,12)
33            format(1x,i2,2i5,5(1pe12.4),/,13x,4(1pe12.4),0pf12.2,
     :               /,13x,2(1pe12.4))
34          continue
          end if
c          write(6,33)
c          write(9,33)
c          do 342 ivp=nvminp,nvmaxp
c             i=ipvp(ivp)
c             gam=rpar(noffp+i*nprpi+12)
c             if(gam.ne.0) then
c               write(6,341) ivp,gam
c               write(9,341) ivp,gam
c             end if
c341          format(1x,i2,'  SPIN-ROTATION CONSTANT=',9x,d12.4)
c342       continue
          if(isg.ne.0) then
            write(9,35)
            write(6,35)
35          format(/' 2SIG-PARAMETERS:'/
     :         '  V  JMIN JMAX',5x,'BROT',8x,'DROT',8x,'HROT',10x,'GS',
     :                9x,'GSD',/,20x,'GSH',8x,'EVIB')
            do 37 ivs=nvmins,nvmaxs
              i=ipvs(ivs)
              nmns=ispar(nmins+2*i)
              nmxs=ispar(nmaxs+2*i)
              write (9, 36) ivs,nmns,nmxs,(rpar(noffs+i*nprsg+k),
     :                      k=1,nprsg)
              write (6, 36) ivs,nmns,nmxs,(rpar(noffs+i*nprsg+k),
     :                      k=1,nprsg)
36            format(1x,i2,2i5,5(1pe12.4),/,13x,1pe12.4,0pf12.2)
37          continue
          end if
        end if
        if(isgpi.gt.0) then
          write(9,38)
          write(6,38)
38        format(/' SIGMA-PI COUPLING PARAMETERS:'/
     :            ' VPI VSIG',9x,'ALPHA',8x,'BETA')
          do 39 ivp=nvminp,nvmaxp
            i=ipvp(ivp)
            ivs=ispar(nvoff+3*i)
            if(ivs.lt.nvmins.or.ivs.gt.nvmaxs) goto 39
            write(9,391) ivp,ivs,(rpar(noffp+i*nprpi+k),k=13,14)
            write(6,391) ivp,ivs,(rpar(noffp+i*nprpi+k),k=13,14)
39        continue
391       format(1x,i2,i5,5x,2g12.5)
        end if
        if (.not. flagsu) write (9, 392) rcut
392     format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :         f6.2)
      end if
*
*    every vibrational energy level is referred to the lowest
*    vibrational level
*
      n=0
*  first transformation matrix for unperturbed sigma states
      nvec=1
      do 85 i=1,3
      do 85 k=1,3
85    vec(k,i,1)=0
      vec(3,3,1)=1.0d0
      if (ipi.eq.0) then
        nlevel = 0
        npi = 0
      else
*  first set up list of all case (a) 2pi levels for omega=1/2
*  for homonuclear molecules in gerade electronic states the e levels
*  are s for j=1/2, a for j=3/2, etc. while the f levels are a for
*  j=1/2, s for j=3/2, etc.  this reverses for ungerade states.
*  see table i of alexander and corey, j. chem. phys. 84, 100 (1986)
*  loop over pi-state vibrational levels
        do 80 iip=nvminp,nvmaxp
        ivp=ipvp(iip)
        jmin=ispar(jminp+3*ivp)
        jmax=ispar(jmaxp+3*ivp)
        n1=n
        do 40 ji = jmin, jmax
        nptop=2
        npbot=1
        if (nparpi .eq. 1) npbot=2
        if (nparpi .eq. -1) nptop=1
        do 40 ip = npbot,nptop
          ipar = iepp(ip) * (-1) ** ji
          if (ihomo .and. isa.ne.0
     :              .and. ipar*igupi .ne. isa) goto 40
*  omega = 1/2 Pi levels
          n = n + 1
          if(n.gt.nmax) goto 1000
          j(n) = ji
          iom(n) = 0
          is (n) = iepp(ip)
*  dimension of h-matrix is temporarily stored in ivec
          ivec(n) = 1
          nlv(n)=iip
*  default is unperturbed, indicated by nlvp(n)=0
          nlvp(n)=0
40      continue
        n2=n
*  now add the omega = 3/2 pi levels
        do 50 i = n1+1, n2
          if (j(i) .eq. 0) goto 50
          n = n + 1
          if(n.gt.nmax) goto 1000
*  pointer to this level temporarily stored in ishold
          ishold(i) = n
          j(n) = j(i)
          iom(n) = 1
          is (n) = is(i)
          ivec(i) = 2
          ivec(n) = nvec+i-n1
          nlv(n)=nlv(i)
*  default is unperturbed, indicated by nlvp(n)=0
          nlvp(n)=0
50      continue
        if(isg.ne.0.and.isgpi.gt.0) then
*  perturbing sigma levels
          ivps=ispar(nvoff+ivp*3)
          do 51 iis=nvmins,nvmaxs
          ivs=ipvs(iis)
51        if(ivps.eq.iis) goto 52
          ivs=0
52        if(ivs.gt.0) then
            nmn=ispar(nmins+2*ivs)
            nmx=ispar(nmaxs+2*ivs)
            do 60 nnrot = nmn,nmx
            do 60 ip = 1, nparsg
*  include only the eps=+1 level for n=0
              if (nnrot.eq.0 .and. ieps(ip).eq.-1) go to 60
*  now calculate j for each level
*  actual half integer value of j is ji + 1/2
              if (ieps(ip) .eq. -1) then
                ji = nnrot-1
              else
                ji = nnrot
              end if
              if(ihomo .and. isa.ne.0
     :                 .and. ieps(ip)*(-1)**ji.ne.isplus) goto 60
              n = n + 1
              if(n.gt.nmax) goto 1000
              is(n) = ieps(ip)
              j(n) = ji
*  assign iom=2 to sigma levels to distinguish them from pi levels
              iom(n) = 2
              nlv(n) = ivibs(ivs)
*  default is unperturbed, indicated by nlvp(n)=0
              nlvp(n) = 0
*  now find perturbing Pi level
              do 55 i=n1+1,n2
                if(j(i).eq.ji.and.is(i).eq.is(n)) then
                  ivec(i)=3
                  ivec(n)=nvec+i-n1
*  pointer to perturbing p-level temporarily stored in jhold
                  jhold(i)=n
*  nlvp are the perturbing vibrational levels
                  nlvp(i)=nlv(n)
                  nlvp(ishold(i))=nlv(n)
                  nlvp(n)=nlv(i)
                  goto 60
                end if
55            continue
*  no level found, so calculate sigma energy
              call hsgpi(j(n),is(n),e(1,1),rpar(noffs+ivs*nprsg+1),
     :                   rpar(noffp+ivp*nprpi+1),3,-1)
              eint(n)=e(3,3)
              ivec(n)=1
60          continue
          end if
        end if
*
*  now set up and diagonalize the h-matrix for each j
*  the matrix elements are given by a. j. kotlar, r. w. field,
*  and j. i. steinfeld, j. mol. spectr. 80, 86 (1980)
*
        do 70 i=n1+1,n2
          idim=ivec(i)
          nvec=nvec+1
          ivec(i)=nvec
*  hsgpi returns in e the idim*idim h-matrix
          call hsgpi(j(i),is(i),e(1,1),rpar(noffs+ivs*nprsg+1),
     :               rpar(noffp+ivp*nprpi+1),-idim,isgpi)
          call jacobi(idim,e(1,1),3,vec(1,1,nvec),3,eig)
          eint(i)=eig(1)
          if(j(i).ne.0) then
            eint(ishold(i))=eig(2)
            if(ivec(ishold(i)).ne.nvec) stop 'chaos 1'
          else
            do 64 k=1,3
            vec(k,2,nvec)=0
64          vec(2,k,nvec)=0
          end if
          if(idim.eq.3) then
            eint(jhold(i))=eig(3)
            if(ivec(jhold(i)).ne.nvec) stop 'chaos 2'
          else
            do 65 k=1,3
            vec(k,3,nvec)=0
65          vec(3,k,nvec)=0
          end if
70        continue
80      continue
        npi=n
      end if
      n=npi
* ----------------------------------------------------------------------
* now add remaining sigma levels
*  first set up list of all case (a) levels
*  the levels are
*     s for sigma-g-plus or for sigma-u-minus
*     a for sigma-g-minus or for sigma-g-plus
      if(isg.ne.0) then
        do 110 iis=nvmins,nvmaxs
        ivs=ipvs(iis)
        if(isgpi.gt.0) then
*  skip levels already included
          do 90 iip=nvminp,nvmaxp
          ivp=ipvp(iip)
          if(ispar(nvoff+ivp*3).eq.iis) goto 110
90        continue
          ivp=1
        end if
        nmn=ispar(nmins+2*ivs)
        nmx=ispar(nmaxs+2*ivs)
        do 100 nnrot = nmn,nmx
          do 100 ip = 1, nparsg
*  include only the eps=+1 level for n=0
            if (nnrot.eq.0 .and. ieps(ip).eq.-1) go to 100
*  now calculate j for each level
*  actual half integer value of j is ji + 1/2
            if (ieps(ip) .eq. -1) then
              ji = nnrot -1
            else
              ji = nnrot
            end if
            if(ihomo .and. isa.ne.0
     :               .and. ieps(ip)*(-1)**ji.ne.isplus) goto 100
            n = n + 1
            if(n.gt.nmax) goto 1000
            is(n) = ieps(ip)
            j(n) = ji
*  assign iom=2 to sigma levels to distinguish them from pi levels
            iom(n) = 2
            ivec(n) = 1
            nlv(n) = ivibs(ivs)
            nlvp(n) = 0
*  now assign energies for case (a) level and store in array eint
*  the matrix elements are given by a. j. kotlar, r. w. field,
*  and j. i. steinfeld, j. mol. spectr. 80, 86 (1980)
            call hsgpi(j(n),is(n),e(1,1),rpar(noffs+ivs*nprsg+1),
     :                rpar(noffp+ivp*nprpi+1),3,-1)
            eint(n)=e(3,3)
100       continue
110     continue
      end if
      nlevel = n
*  nlevel now contains the total number of distinct levels
*
*  form list of all energetically distinct rotational levels  (case(a))
*  included in the  channel basis and their energies.
*  sort this list according to omega and v
*
      n=0
      nvm=nvminp
      nvx=nvmaxp
      emin=1.d10
      do 120 io=0,2
      if(io.eq.2) then
        nvm=nvmins
        nvx=nvmaxs
      end if
      do 120 ivv=nvm,nvx
      n1=n+1
      do 115  i = 1, nlevel
        if(iom(i).eq.io.and.nlv(i).eq.ivv) then
          n=n+1
          ehold(n) = eint(i) / econv
          emin=min(emin,ehold(n))
          jhold(n) = j(i)
          ishold(n) = is(i)
          iohold(n) = iom(i)
          ivhold(n) = ivec(i)
          nvhold(n) = nlv(i)
          nvphol(n) = nlvp(i)
        end if
115   continue
120   continue
      if(n.ne.nlevel) stop 'chaos 3'
      do 125 i=1,n
        ehold(i)=ehold(i)-emin
        eint(i)=ehold(i)
        j(i)=jhold(i)
        is(i)=ishold(i)
        iom(i)=iohold(i)
        ivec(i)=ivhold(i)
        nlv(i)=nvhold(i)
        nlvp(i)=nvphol(i)
125   continue

*
*  set up coupled-states channel basis (if desired)
      if (csflag) then
*
*  delete channels with j less than coupled states projection index
*
        n=0
        do 130 i=1,nlevel
          if(j(i).ge.nu) then
            n=n+1
            j(n)=j(i)
            eint(n)=eint(i)
            is(n)=is(i)
            iom(n)=iom(i)
            ivec(n)=ivec(i)
            nlv(n)=nlv(i)
            nlvp(n)=nlvp(i)
            l(n) = jtot
            cent(n) = jtot * (jtot + 1)
          end if
130     continue
*
*  set up close-coupled channel basis (if desired)
*
      else if (.not. csflag) then
*
        n=0
        do 145  i = 1, nlevel
*  now duplicate  as many times as is required for rotational degeneray
          ji = jhold(i)
          lmax = jtot + ji + 1
          lmin = iabs (jtot - ji)
          do 140  li = lmin, lmax
            ix = (-1) ** (ji + li - jtot) * ishold(i)
            if (ix .eq. jlpar) then
              n = n + 1
              if (n .gt. nmax) goto 1000
              eint(n) = ehold(i)
              j(n) = jhold(i)
              is(n) = ishold(i)
              iom(n) = iohold(i)
              ivec(n) = ivhold(i)
              nlv(n) = nvhold(i)
              nlvp(n) = nvphol(i)
              l(n) = li
              cent(n) = li * ( li + 1)
            end if
140       continue
145     continue
      end if
*  n now contains the number of 2pi and 2sigma channels
*  now sort the complete (sigma+pi) channel list to put closed levels
*  at end. Also determine number of levels which are open
      nlevop = 0
      do 160  i = 1, nlevel - 1
        if (ehold(i) .le. ered) then
          nlevop = nlevop + 1
        else
          do 155 ii = i + 1, nlevel
            if (ehold(ii) .le. ered) then
              nlevop = nlevop + 1
              call iswap(jhold(i),jhold(ii))
              call iswap(ishold(i),ishold(ii))
              call iswap(nvhold(i),nvhold(ii))
              call iswap(iohold(i),iohold(ii))
              call rswap(ehold(i),ehold(ii))
              go to 160
            end if
155       continue
        end if
160   continue
      if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
* ----------------------------------------------------------------------
*  now check to see if any of the open channels are closed at r=rcut
*  this is not done for molecule-surface collisions of for rcut < 0
      if (.not.flagsu .and. rcut .gt. 0.d0 .and..not.boundc) then
        emin = 1. e+7
        do 220  i = 1, n
          if (eint(i) .le. ered) then
*  here if channel is open asymptotically
            if ( jtot * (jtot + 1) / (two * rmu * rcut * rcut)
     :          .gt. (ered - eint(i)) ) then
*  here if channel is open asymptotically but closed at r = rcut
              if (eint(i) .lt. emin) emin = eint(i)
*  emin now contains the lowest channel energy for which this condition
*  is met
            end if
          end if
220     continue
*  now eliminate all channels with eint .ge. emin if any of the channels
*  are open asymptotically but closed at r = rcut
        if (emin.lt.ered) then
          nn=n
          n = 0
          do 230 i = 1, nn
            if (eint(i) .lt. emin) then
*  here if this channel is to be included
              n = n + 1
              eint(n) = eint(i)
              j(n) = j(i)
              is(n) = is(i)
              cent(n) = cent(i)
              iom(n) = iom(i)
              ivec(n) = ivec(i)
              nlv(n) = nlv(i)
              nlvp(n) = nlvp(i)
              l(n) = l(i)
            end if
230       continue
*  reset number of channels
        end if
      end if
*  if no channels, return
      if (n .eq. 0) return
      if (nu .eq. numin) then
        ntop = max(n, nlevop)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1
      else
        if (n.gt.ntop) then
          write (6, 232) nu, n, ntop
          write (9, 232) nu, n, ntop
232       format (' *** NCH = ',i3, ' AT NU = ',i2, ' .GT. NTOP = ',i3,
     :            '; ABORT **',/,
     :    '     CHECK RCUT')
          call exit
        endif
      end if
*  now list channels if requested
      if (clist .or. bastst .and. iprint .ge. 1) then
        write (6, 250)
        write (9, 250)
250     format(/1x,'   N   V  VP   J  EPS OMEGA    L  EINT(CM-1)',
     :   '     C-1/2       C-3/2       C-SIG')
        do 280  i = 1, n
          fj = j(i) + half
          xmga = iom(i) + half
          ecm = eint(i) * econv
          ive = ivec(i)
          ico = iom(i) + 1
          write (6, 270) i,nlv(i),nlvp(i),fj,is(i),xmga,l(i), ecm,
     :                               (vec(k,ico,ive),k=1,3)
          write (9, 270) i, nlv(i),nlvp(i),fj, is(i), xmga, l(i), ecm,
     :                               (vec(k,ico,ive),k=1,3)
270       format (1x,3i4, f5.1, i4, f5.1, i6, f10.3, 3f12.6)
280     continue
      end if
*  now calculate coupling matrix elements
      if (bastst .and.iprint.gt.1) then
        write (6, 285)
        write (9, 285)
285     format (/'  IT VR VC  LAMBDA ILAM    I   ICOL  IROW',
     :           '  IV2       VEE')
      end if
*
*  the following only for safety. Checks that pot routine
*  contains all required vib. levels in correct order
*
      do 295 irow=1,n
      ior=iom(irow)
      ivr=nlv(irow)
      do 290 icol=1,irow
      ioc=iom(icol)
      ivc=nlv(icol)
      iterm=0
      if(ior.eq.ioc) then
        if(ior.eq.2) iterm=1
        if(ior.le.1) iterm=2
      else
        if(ior.eq.1) iterm=3
        if(ior.eq.2) iterm=4
      end if
      if(iterm.eq.0) stop 'chaos 4'
      do 290 ivb=1,ntv(iterm)
290   if(ivrow(ivb,iterm).eq.ivr.and.
     :   ivcol(ivb,iterm).eq.ivc) goto 295
      write(6,291) ivr,ivc
      write(9,291) ivr,ivc
291   format(/' VIBRATIONAL STATE NOT DEFINED IN POT',2i5)
      call exit
295   continue
      i=0
*
*  i counts v2 elements
*  ilam counts number of v2 matrices
*  inum counts v2 elements for given ilam
*  ij is address of given v2 element in present v2 matrix
*
      ilam=0
      do 600 it=1,nterm
        iterm=it
*  formally reassign iterm if no sigma state present
        if(nterm.eq.2) then
           if(it.eq.1) iterm=2
           if(it.eq.2) iterm=4
        end if
*  ior, ioc will be 1 for pi and 2 for sigma
        ior=iomr(iterm)
        ioc=iomc(iterm)
*  loop over row vibrational levels
        do 500 iv=1,ntv(it)
*  loop over column vibrational levels
          ivc=ivcol(iv,it)
          ivr=ivrow(iv,it)
          do 450 il=lammin(it),lammax(it),istep
            ilam=ilam+1
            inum = 0
            do 400  icol = 1, n
            if(max(1,iom(icol)).ne.ioc) then
*  skip unless level of channel icol is perturbed by vibr.level ivc
              if(nlvp(icol).ne.ivc) goto 400
            else
*  skip unless vibrational level of channel icol equals ivc
              if(nlv(icol).ne.ivc) goto 400
            end if
            do 350  irow = icol, n
            if(max(1,iom(irow)).ne.ior) then
              if(nlvp(irow).ne.ivr) goto 350
            else
              if(nlv(irow).ne.ivr) goto 350
            end if
            ij = ntop * (icol - 1) +irow
            lrow = l(irow)
            if (csflag) lrow = nu
*  always initialize potential to zero
            vee = 0
            do 300 mm=1,3
            do 300 nn=1,3
300         e(nn,mm)=0
*  ilam is the angular expansion label
*  here for the sigma-sigma potential
*  here only if initial and final states are 2 sigma
*  lb is the actual value of lambda
               if(ipi.ne.0) then
                 if(iterm.eq.2) then
                   call vlsgpi (j(irow), lrow, j(icol), l(icol), jtot,
     :                         izero, izero, il,
     :                         is(irow), is(icol), e(1,1), csflag)
                   call vlsgpi (j(irow), lrow, j(icol), l(icol), jtot,
     :                         ione, ione, il,
     :                         is(irow), is(icol), e(2,2), csflag)
                 else if(iterm.eq.4) then
                   call vlsgpi (j(irow), lrow, j(icol), l(icol), jtot,
     :                         ione, izero, il,
     :                         is(irow), is(icol), e(2,1), csflag)
                   call vlsgpi (j(irow), lrow, j(icol), l(icol), jtot,
     :                         izero, ione, il,
     :                         is(irow), is(icol), e(1,2), csflag)
                 end if
               end if
               if(isg.ne.0) then
                 if(iterm.eq.1) then
                   call vlsgpi (j(irow), lrow, j(icol), l(icol), jtot,
     :                         itwo, itwo, il,
     :                         is(irow), is(icol), e(3,3), csflag)
                 else if(iterm.eq.3.and.ipi.ne.0) then
                   call vlsgpi (j(irow), lrow, j(icol), l(icol), jtot,
     :                         itwo, izero, il,
     :                         is(irow), is(icol), e(3,1), csflag)
                   call vlsgpi (j(irow), lrow, j(icol), l(icol), jtot,
     :                         itwo, ione, il,
     :                         is(irow), is(icol), e(3,2), csflag)
                   call vlsgpi (j(irow), lrow, j(icol), l(icol), jtot,
     :                         izero, itwo, il,
     :                         is(irow), is(icol), e(1,3), csflag)
                   call vlsgpi (j(irow), lrow, j(icol), l(icol), jtot,
     :                         ione, itwo, il,
     :                         is(irow), is(icol), e(2,3), csflag)
                 end if
               end if
               icc=iom(icol)+1
               icr=iom(irow)+1
               itc=ivec(icol)
               itr=ivec(irow)
c               call mxva(e,1,3,vec(1,icc,itc),1,eig,1,3,3)
c JK replaced with dgemv call
              call dgemv('n',1,3,1.d0,e,vec(1,icc,itc),1,0.d0,eig,1)
               vee=ddot(3,eig,1,vec(1,icr,itr),1)
             if (abs(vee) .gt. tzero) then
               i = i + 1
               inum = inum + 1
               if (bastst .and. iprint.gt.1) then
                 write (6, 340) it,ivr,ivc,il,ilam,i,icol,irow,
     :                          ij,vee
                 write (9, 340) it,ivr,ivc,il,ilam,i,icol,irow,
     :                          ij,vee
340              format (1x,3i3,6i6, g17.8)
               end if
               if (i .le. nv2max) then
                 v2(i) = vee
                 iv2(i) = ij
               end if
             end if
350        continue
400        continue
           if(ilam.le.nlammx) lamnum(ilam) = inum
           if (bastst .and.iprint.ge.1) then
             write (6, 420) it,ivr,ivc,il,inum
             write (9, 420) it,ivr,ivc,il,inum
420          format(' ITERM=',i1,'  IVR=',i2,'  IVC=',i2,'  LAMBDA=',i2,
     :              ' NUMBER OF NONZERO MATRIX ELEMENTS',i8)
           end if
           if(inum.ne.0) nlam=ilam
450      continue
500    continue
      continue
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
*  set epsilon equal to +/- 2 for nominally omega=3/2 pi states
*  and equal to +/- 3 for sigma states
      do 640 i = 1, n
        is(i) = ((1 + iom(i))*100 + nlv(i)) * is(i)
640   continue
      do 650 i=1,nlevel
        ishold(i) = ((1 + iohold(i))*100 + nvhold(i)) * ishold(i)
650   continue
      if (bastst .and. iprint .ge. 1) then
        write (6, 660)
660     format(/1x,'   N   V  VP   J   EPS    L  EINT(CM-1)',
     :   '    C-1/2       C-3/2       C-SIG')
        do 680  i = 1, n
          fj = j(i) + half
          ecm = eint(i) * econv
          ive = ivec(i)
          ico = iom(i) + 1
          write (6, 670) i,nlv(i),nlvp(i),fj,is(i),l(i), ecm,
     :                               (vec(k,ico,ive),k=1,3)
670       format (1x,3i4, f5.1, i5, i6, f10.3, 3f12.6)
680     continue
      end if
      return
1000  write (6, 1010) n, nmax
      write (9, 1010) n, nmax
1010  format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF',
     :         i4,' ABORT ***')
      call exit
      end
* --------------------------------------------------------------------
      subroutine vlsgpi (jp, lp, j, l, jtot, iomegp, iomeg, lambda,
     :                   iepsp, ieps, v, csflag)
*  subroutine to calculate v-lambda matrices for close-coupled and
*  coupled-states treatments of collisions of a molecule in a 2pi or 2sigma
*  electronic state
*  the cc matrix elements are given in eq. (29) of m.h. alexander, chem. phys.
*  92, 337 (1985)
*  the cs matrix elements are given in eq. (14) of t. orlikowski and m.h.
*  alexander, j. chem. phys. 79, 6006 (1983)
*  note that for cc collisions of a molecule with a flat surface, the
*  coupling matrix elements [m.h. alexander, j. chem. phys. 80, 3485 (1984)]
*  are identical to the cs matrix elements here
*  author:  millard alexander
*  current revision date:  6-oct-87
*  -----------------------------------------------------------------------
*  variables in call list:
*    jp:       rotational quantum number of left side of matrix element (bra)
*              minus 1/2
*    lp:       orbital angular momentum of left side of matrix element (bra)
*    j:        rotational quantum number of right side of matrix element (ket)
*              minus 1/2
*    l:        orbital angular momentum of right side of matrix element (ket)
*    jtot:     total angular momentum
*    iomegp:   omega quantum number of bra
*    iomeg:    omega quantum number of ket
*              we assume that
*              iomega=0 corresponds to the omega=1/2 2pi state
*              iomega=1 corresponds to the omega=3/2 2pi state
*              iomega=2 corresponds to the 2sigma state
*    lambda:   order of legendre term in expansion of potential
*    iepsp:    symmetry index of bra
*    ieps:     symmetry index of ket
*    v:        on return, contains matrix element
*    csflag:   if .true., then cs calculation; in this case:
*                jtot in call list is cs lbar
*                lp is nu (cs projection index) minus 1/2
*                l is not used
*              if .false., then cc calculation; in this case:
*                jtot is total angular momentum minus 1/2
*    for collisions of a 2pi molecule with a surface, nu is equivalent
*    to m (the projection of j along the surface normal) minus 1/2
*  subroutines called:
*     xf3j, xf6j
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
*      real half, halfm, one, onem, threhf, thrhfm, two, v, x, xj,
*     :     xjp, xjtot, xl, xlamda, xlp, xnorm, xnu, xnum, xomeg,
*     :     xomegm, xomegp, xomgpm, zero
*      real xf3j, xf6j
      integer ieps, iepsp, iomeg, iomegp, ione, iphase, j, jp, jtot,
     :        l, lambda, lp, nu
      logical csflag
      data thrhfm, onem, halfm, zero, half, ione, one, threhf, two
     :   /-1.5d0, -1.d0,-0.5d0,0.0d0,0.5d0,   1, 1.d0, 1.5d0,  2.d0/
      v = zero
      xjp = jp + half
      xj = j + half
      xomegp = iomegp + half
      xomeg = iomeg + half
      if (iomegp .eq. 2) xomegp = half
      if (iomeg .eq. 2) xomeg = half
      xjtot = jtot + half
      xomegm = - xomeg
      xomgpm = - xomegp
      if (csflag) then
        nu = lp
        xnu = nu + half
      end if
      xlp = lp
      xl = l
      xlamda = lambda
      iphase = ieps * iepsp * ((-1) ** (jp + j + lambda + 1))
      if (iphase .eq. 1) return
      if (csflag) then
        iphase = nu - iomeg
        xnorm = (2. * xjp + 1.) * (2. * xj + 1.)
        xnu = nu + half
        xnum = - xnu
        x = xf3j (xjp, xlamda, xj, xnum, zero, xnu)
      else
        iphase = jp + j + ione + jtot - iomeg
        xnorm = (2. * xjp + 1.) * (2. * xj + 1.) * (2. * lp + 1.)
     :        * (2. * l + 1.)
        x = xf3j (xlp, xlamda, xl, zero, zero, zero)
        if  (x .eq. zero) return
        x = x * xf6j (xjp, xlp, xjtot, xl, xj, xlamda)
      end if
      if (x .eq. zero) return
      if (iomegp .eq. iomeg) then
*  here for omega=omega' for 2pi-2pi or 2sigma-2sigma coupling
        x = x * xf3j (xjp, xlamda, xj, xomgpm, zero, xomegp)
      else if ((iomegp .eq. 0 .and. iomeg .eq. 1) .or.
     :         (iomegp .eq. 1 .and. iomeg .eq. 0) ) then
*  here for omega.ne.omega' for 2pi-2pi coupling
*  note that phase is the negative of what was erroneously given in eq. (29)
*  of m.h. alexander, chem. phys. 92, 337 (1985)
        x = - ieps * x * xf3j (xjp, xlamda, xj, xomgpm, two, xomegm)
      else if ((iomegp. eq. 0 .and. iomeg .eq. 2) .or.
     :         (iomegp. eq. 2 .and. iomeg .eq. 0) ) then
*  here for coupling between 2sigma and 2pi (omega=1/2)
*  see eq. 36 of m.h. alexander and g.c. corey, j. chem. phys. 84, 100
*  (1986)
        x = x * ieps * xf3j (xjp, xlamda, xj, halfm, one, halfm)
      else if (iomegp. eq. 1 .and. iomeg .eq. 2) then
*  here for <2sigma/ v/ 2pi (omega)=3/2>
*  see eq. 37 of m.h. alexander and g.c. corey, j. chem. phys. 84, 100
*  (1986)
        x = x * xf3j (xjp, xlamda, xj, thrhfm, one, half)
      else if (iomegp. eq. 2 .and. iomeg .eq. 1) then
*  here for <2pi (omega)=3/2/ v/ 2sigma>
*  see eq. 38 of m.h. alexander and g.c. corey, j. chem. phys. 84, 100
*  (1986)
        x = x * xf3j (xjp, xlamda, xj, halfm, onem, threhf)
      end if
      iphase = (-1) ** iabs(iphase)
      v = iphase * x * sqrt(xnorm)
      return
      end
c
c       subroutine to diagonalize square matrix
c       n is the first dimension of the matrices h and u
c       h is the matrix to be diagonalized, u is the matrix
c       of eigenvectors as columns and e is the vector
c       of eigenvalues.
c
        subroutine jacobi(n,h,nh,u,nu,e)
        implicit double precision (a-h,s-z)
        implicit integer (n,p-r)
        dimension h(nh,n),u(nu,n),e(n)
        logical rot
        th=0.0d0
        do 10 p=1,n
          do 20 q=1,p-1
            u(p,q)=0.0d0
            u(q,p)=0.0d0
            th=th+h(p,q)*h(p,q)
   20     continue
          u(p,p)=1.0d0
          e(p)=h(p,p)
   10   continue
        th=dsqrt(th)/(n*n)
        thmin=th*1.0d-10
   30   rot=.false.
        do 40 p=1,n
          do 50 q=p+1,n
            if (dabs(h(p,q)).gt.th) then
              rot=.true.
              theta=(e(q)-e(p))/h(p,q)*0.5d0
              t=dsign(1.0d0,theta)/
     1        (dabs(theta)+dsqrt(theta*theta+1.0d0))
              c=1.0d0/dsqrt(1.0d0+t*t)
              s=t*c
              do 60 r=1,n
                a=u(r,p)
                u(r,p)=c*a-s*u(r,q)
                u(r,q)=s*a+c*u(r,q)
   60         continue
              do 70 r=1,p-1
                a=h(r,p)
                h(r,p)=c*a-s*h(r,q)
                h(r,q)=s*a+c*h(r,q)
   70         continue
              do 80 r=p+1,q-1
                a=h(p,r)
                h(p,r)=c*a-s*h(r,q)
                h(r,q)=s*a+c*h(r,q)
   80         continue
              do 90 r=q+1,n
                a=h(p,r)
                h(p,r)=c*a-s*h(q,r)
                h(q,r)=s*a+c*h(q,r)
   90         continue
              a=e(p)
              e(p)=c*c*a+s*s*e(q)-2.0d0*s*c*h(p,q)
              e(q)=s*s*a+c*c*e(q)+2.0d0*s*c*h(p,q)
              h(p,q)=0.0d0
            endif
   50     continue
   40   continue
        if (rot) goto 30
        th=0.1d0*th
        if (th.gt.thmin) goto 30
        return
        end
* ----------------------------------------------------------------------
      subroutine iswap(ia,ib)
      implicit integer (a-z)
      ii=ia
      ia=ib
      ib=ii
      return
      end
* ----------------------------------------------------------------------
      subroutine rswap(a,b)
      implicit double precision (a-h,o-z)
      c=a
      a=b
      b=c
      return
      end
* ----------------------------------------------------------------------
      subroutine hsgpi(ji,ise,e,rp1,rp2,iflag,iperpi)
* ----------------------------------------------------------------------
* subroutine to determine full doublet-sigma and/or doublet-pi
* hamiltonian
* written by alessandra degli-esposti
* current  revision:  06-mar-1992 by pjd
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical if1,if2,if3
*
*        pi state constants in input, they depend on v  (array rp2)
*    1-3 brot,drot,hrot
*      4 aso
*    5-6 p,pd
*    7-9 q,qd,qh
*     10 evib
*     11 ad
*     12 gpi
*        sigma state constants in input, they depend on v (array rp1)
*    1-3 brotsg,drotsg,hrotsg
*    4-6 gsg,gsgd,gsgh
*      7 evib
*        sigma!pi mixing parameters (array rp2)
*     13 alpha
*     14 beta
*
*        1 = pi 1/2
*        2 = pi 3/2
*        3 = sigma
*
       dimension e(3,3),rp1(7),rp2(14)
*
         if1=(iflag.eq.1.or.iflag.le.-1)
         if2=(iflag.eq.2.or.iflag.le.-2).and.ji.ne.0
         if3=(iflag.eq.3.or.iflag.le.-3)
         x =  float(ji) + 1.0d0
         x2 = x**2
         x3 = x**3
         x4 = x2 * x2
         do 10 i=1,3
         do 10 j=1,3
10       e(j,i)=0
*
*
*        <pi1/2!pi1/2>  constants
*
         if(if1) then
           tepi = rp2(10)
           brotp1 = rp2(1) *x2
           drotp1 = rp2(2) *(1.d0-x2-x4)
           hrotp1 = rp2(3) *(x**6 + 3.d0*(x2-1.d0)**2 +x2-1.d0)
           asop1  = rp2(4) *(-0.5d0)
           pp1    = rp2(5) *0.5d0 *(1.d0-ise*x)
           pdp1   = rp2(6) *0.5d0 *(1.d0-ise*x) *(x2-.25d0)
           qp1    = rp2(7) *0.5d0 *(1.d0-ise*x)**2
           qdp1   = rp2(8) *0.5d0 *(1.d0-ise*x)**2 *(x2-.25d0)
           qhp1   = rp2(9) *0.5d0 *(1.d0-ise*x)**2 *(x2-.25d0)**2
           ad1    = rp2(11) * (-x2)
           gp1    = rp2(12) * (-0.5d0)
           e(1,1) = asop1 + brotp1 + drotp1 + hrotp1  + pp1 + pdp1
     >           + qp1 + qdp1 +  qhp1 + ad1 + tepi + gp1
         end if
*
*        <sigma!sigma>  constants
*
         if(if3) then
           tesg = rp1(7)
           brotsg = rp1(1) *x  *(x-ise)
           drotsg = rp1(2) * (-x2) *(x-ise)**2
           hrotsg = rp1(3) * x3 *(x-ise)**3
           gsg    = rp1(4) *(-0.5d0) *(1.d0-ise*x)
           gsgd   = rp1(5) *(-0.5d0) *(1.d0-ise*x) *(x2-.25d0)
           gsgh   = rp1(6) *(-0.5d0) *(1.d0-ise*x) *(x2-.25d0)**2
           e(3,3) = brotsg + drotsg + hrotsg + gsg + gsgd + gsgh + tesg
         end if
*
*         <sigma!pi1/2> constants
*
         if(if1.and.if3) then
           alphas1= rp2(13)
           betas1 = rp2(14)*(1.d0-ise*x)
           e(1,3) = alphas1 + betas1
           e(3,1) = e(1,3)
         end if
*
         if(if2) then
*
*        <pi3/2!pi3/2>  constants
*
             tepi = rp2(10)
             brotp2 = rp2(1) *(x2-2.d0)
             drotp2 = rp2(2) *(1.d0-x2 -(x2- 2.d0)**2)
             hrotp2 = rp2(3) *((x2-2.d0)**3 + 3.d0*(x2-1.d0)**2
     >                       -x3+1.d0)
             asop2  = rp2(4) *0.5d0
             qp2    = rp2(7) *0.5d0*(x2-1.d0)
             qdp2   = rp2(8) *0.5d0*(x2-1.d0) *(x2-.25d0)
             qhp2   = rp2(9) *0.5d0*(x2-1.d0) *(x2-.25d0)**2
             ad2    = rp2(11) * (x2-2.d0)
             e(2,2) = asop2+brotp2+drotp2+hrotp2+qp2+qdp2+qhp2
     >                       +ad2+tepi
*
*        <sigma!pi3/2>
*
           if(if3) then
             betas2 = rp2(14) *(-sqrt(x2-1.d0))
             e(2,3) = betas2
             e(3,2) = e(2,3)
           end if
*
*          <pi1/2!pi3/2> constants
*
           if(if1.and.iperpi.ge.0) then
             brot12 = rp2(1) *(-sqrt(x2-1.d0))
             drot12 = rp2(2) *2.d0*sqrt(x2-1.d0)**3
             hrot12 = rp2(3) *(-sqrt(x2-1.d0) *(3.d0*(x2-1.d0)**2 +x2))
             p12  = rp2(5) *(-.25d0)*sqrt(x2-1.d0)
             pd12 = rp2(6) *(-.25d0)*sqrt(x2-1.d0) *(x2-.25d0)
             q12  = rp2(7) *(-.5d0)*(1.d0-ise*x)*sqrt(x2-1.d0)
             qd12 = rp2(8) *(-.5d0)*(1.d0-ise*x)*sqrt(x2-1.d0)
     >                     *(x2-.25d0)
             qh12 = rp2(9) *(-.5d0)*(1.d0-ise*x)*sqrt(x2-1.d0)
     >                     *(x2-.25d0)**2
             gp12 = rp2(12) * (-0.5d0) * sqrt(x2 - 1.0d0)
             e(1,2) = brot12 + drot12 + hrot12 + p12 + pd12 +
     >              q12 + qd12 + qh12 + gp12
             e(2,1) = e(1,2)
           end if
         end if
         return
         end
