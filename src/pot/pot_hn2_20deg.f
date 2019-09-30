*  System:  H + N2(X) RCCSD(T)-F12a
*
*   calculation of potential energy curves by p.dagdigian
*   theta = 10 deg cut
*  
*   written by p. dagdigian
*   current revision date:  15-sep-2016
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(1)
      include "common/parpot"
      econv=219474.6d0
      potnam='H-N2 RCCSD(T)-F12a -- 20 deg'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      if (r.le.0.d0) goto 99
      call pot(vv0,r)
*     potential is returned in atomic units (hartree)
*     convert from atomic units for printout
      write (6, 100) vv0*econv
100   format(9(1pe16.8))
      goto 1
c
99    rr=3.0d0
      dr=0.1d0
      open (unit=12,file='hn2_pot_0deg.dat')
      write(12,109)
109   format(' %R/bohr V00  ...')
      do i=1,250
        call pot(vv0,rr)
        write (12,110) rr,vv0*econv
110     format(f7.2,1pe16.8)
        rr = rr + dr
      enddo
      close(12)
      end
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*(*) filnam
      include "common/parbas"
      include "common/parpot" 
      common /conlam/ nlam, nlammx, lamnum(2)
      common /cosysi/ nscode, isicod, nterm
      potnam='H-N2 RCCSD(T)-F12a -- 20 deg'
      npot=1
      nterm=1
      lammin(1)=2
      lammax(1)=2
      mproj(1)=0
      ipotsy = 2
*
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
*  count number of anisotropic terms
      nlam = 1
      nlammx = nlam
      return
      end
* --------------------------------------------------------------------------
      subroutine pot (vv0, r)


*  subroutine to calculate the r-dependent coefficients in the
*  collision of a homonuclear diatomic with a structureless target
*  in units of hartree for energy and bohr for distance

*  on return:
*  vv0 contains the the potential
*  vvl(1) is set to zero and returned in common block /covvl/

*  variable in common block /conlam/
*    nlammx:    the maximum number of anisotropic terms
*    nlam:      the total number of angular coupling terms
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential
* 
* author:  paul dagdigian
* latest revision date:  15-sep-2016
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      dimension v(9)
      dimension rr(37), vl(37), csplin(37)
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(1)
      data ifirst / 0 /
*
*  37 values or R
      data rr /
     +  2.000,  2.250,  2.500,  2.750,  3.000,
     +  3.250,  3.500,  3.750,  4.000,  4.250,
     +  4.500,  4.750,  5.000,  5.250,  5.500,
     +  5.750,  6.000,  6.250,  6.500,  6.750,
     +  7.000,  7.250,  7.500,  7.750,  8.000,
     +  8.500,  9.000,  9.500, 10.000, 10.500,
     + 11.000, 11.500, 12.000, 13.000, 14.000,
     + 15.000, 20.000 /
* pot:
      data vl /
     +  1.0986321e+05,  4.1598908e+04,  1.4314121e+04,  5.9595979e+03,
     +  6.1416526e+03,  8.7382279e+03,  9.5341797e+03,  7.9701707e+03,
     +  6.0360065e+03,  4.3597956e+03,  3.0546553e+03,  2.0881906e+03,
     +  1.3936912e+03,  9.0667458e+02,  5.7251314e+02,  3.4796642e+02,
     +  2.0045909e+02,  1.0592382e+02,  4.7083175e+01,  1.1859416e+01,
     + -8.0492651e+00, -1.8271581e+01, -2.2599305e+01, -2.3483197e+01,
     + -2.2490508e+01, -1.8304136e+01, -1.3821009e+01, -1.0126731e+01,
     + -7.3598457e+00, -5.3667040e+00, -3.9488869e+00, -2.9380092e+00,
     + -2.2157298e+00, -1.3113719e+00, -8.1757573e-01, -5.3764664e-01,
     + -9.0513507e-02 /
*
* spline fit
      if (ifirst .eq. 0) then
* spline fit
*    evaluate derivative at first point
         der1=(vl(2)-vl(1))/(rr(2)-rr(1))
         call dspline(rr,vl,37,der1,0d0,csplin)
         ifirst = 1
       end if
* r^-6 fit to at R = 15 bohr for long-range part of the potential
       c6sum = -5.2883574764d+06
* switching function for long-range
       switch_lr=0.5*(tanh(0.5*(r - 18.d0)) + 1.d0)
* determine splined coefficients at R=r
       call dsplint(rr,vl,csplin,37,r,vvx)
* kill anisotropic terms at large R
       vvx = (1.d0 - switch_lr)*vvx
* merge with asymptotic form
       vvx = vvx + switch_lr*c6sum/(r**6)
       vvl(1) = 0.d0
* convert to hartree
       econv=1.d0/219474.6d0
       vv0 = vvx * econv
*
       return
       end
*===========================eof===============================