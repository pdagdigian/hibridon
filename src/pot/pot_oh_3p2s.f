*  System:  O(3P) + H(2S)
*
*   calculation of potential energy curves and spin-orbit matrix
*   elements, taken from Yarkony [j. chem. phys. 97, 1838 (1992)] and 
*   Parlant and Yarkony [j. chem. phys. 110, 363 (1999)]
*
*   3P and 1D-3P spin-orbit matrix elements extend to large R.
*   compute 3P spin-orbit splitting using large R values
*  
*   written by p. dagdigian
*   current revision date:  13-feb-2015
*
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      common /coselb/ ibasty
      include "common/parbas"
      include "common/parpot"
      potnam='O(3P)-H(2S)'
      npot=1
      ibasty=23
      lammin(1)=1
      lammax(1)=18
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysi/ junk(5), npot
      common /cosysr/ isrcod, junkr, en(4)
      common /covvl/ vvl(18)
      include "common/parpot"
      potnam='O(3P)-H(2S)'
      econv=219474.6d0
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      if (r.le.0.d0) goto 99
      call pot(vv0,r)
      write (6, 100) (econv*vvl(i), i=1,18)
100   format('  ',18(1pe16.8))
      goto 1
99    rr=1.075d0
      dr=0.2d0
      open (unit=12,file='oh_vlms.txt')
      write(12,109)
109   format(' R/bohr X2Pi  4Pi  2Sig-  4Sig-  A   B   C',
     :  '  D   E   F   G   A(R=inf)   B(R=inf)',
     :  '  C(R=inf)  D(R=inf)  E(R=inf)  F(R=inf)',
     :  '  G(R=inf)')
      do i=1,276
        call pot(vv0,rr)
        write (12,110) rr,(econv*vvl(j),j=1,18)
110     format(f7.2,18(1pe16.8))
        rr = rr + dr
      enddo
      close(12)
      end
* --------------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*    vv0        (not used)
*    vvl(1,4)   contains the X2Pi, 4Pi, 2Sig-, 4Sig- energies
*    vvl(5,11)  contains the spin-orbit matrix elements [A - G]
*    vvl(12,18) contains the spin-orbit matrix elements A - G
*               at R=Rinf 
*  variable in common block /conlam/ used here
*    nlam:      the number of angular coupling terms actually used
*  variable in common block /covvl/
*    vvl:       array to store r-dependence of each
*               angular term in the potential
*  variable in common block /coconv/
*    econv:     conversion factor from cm-1 to hartrees
*  variables in common /cosysr/
*    isrcod:    number of real parameters
*    en1d:      asymptotic energy of the 1D state (cm-1)
* 
* author:  paul dagdigian
* latest revision date:  19-sep-2014
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/parbas"
      common /covvl/ vvl(18)

      dimension rrx(53),vlx(53),rr(54),vl4p(54),vl2s(54),vl4s(54),
     +  rso(31),somat(217),vl(54)
      dimension csx(53),cs4p(54),cs42s(54),cd4s(54),csplso(31,11),vecelx(53),
     +  vecl4p(54),vecl2s(54),vecl4s(54),vecso(31)
      dimension v(18)
*
*  values of R for X2Pi PE curve
      data rrx /
     +  1.05835,  1.19065,  1.32294,  1.45524,  1.58753, 
     +  1.71983,  1.85212,  1.98441,  2.11671,  2.24900, 
     +  2.38130,  2.51359,  2.64589,  2.77818,  2.91047, 
     +  3.04277,  3.17506,  3.30736,  3.43965,  3.57195, 
     +  3.70424,  3.83653,  3.96883,  4.10112,  4.20000, 
     +  4.33229,  4.46459,  4.59688,  4.72918,  4.86147, 
     +  4.99377,  5.12606,  5.25835,  5.39065,  5.52294, 
     +  5.65524,  5.78753,  5.91983,  6.05212,  6.18441, 
     +  6.31671,  6.44900,  6.58130,  6.71359,  6.84589, 
     +  6.97818,  7.11047,  7.24277,  7.37506,  7.50736, 
     +  7.63965,  7.77195,  7.90424 /   
      data vlx /
     +  5.8134742d+04, 1.8325793d+04, -8.1629391d+03, -2.4237418d+04, 
     + -3.2803607d+04, -3.6605972d+04, -3.7442428d+04, -3.6439429d+04, 
     + -3.4350387d+04, -3.1630953d+04, -2.8609891d+04, -2.5493023d+04, 
     + -2.2417476d+04, -1.9478743d+04, -1.6738930d+04, -1.4233948d+04, 
     + -1.1985965d+04, -1.0008076d+04, -8.2950919d+03, -6.8341434d+03, 
     + -5.6111554d+03, -4.6030315d+03, -3.7824933d+03, -3.1222395d+03, 
     + -2.7463846d+03, -2.2800764d+03, -1.9035682d+03, -1.5976310d+03, 
     + -1.3475461d+03, -1.1419595d+03, -9.7204897d+02, -8.3091197d+02, 
     + -7.1311254d+02, -6.1434315d+02, -5.3117037d+02, -4.6084218d+02, 
     + -4.0114093d+02, -3.5027055d+02, -3.0676935d+02, -2.6944226d+02, 
     + -2.3730786d+02, -2.0955680d+02, -1.8551890d+02, -1.6463711d+02, 
     + -1.4644660d+02, -1.3055820d+02, -1.1664490d+02, -1.0443102d+02, 
     + -9.3683439d+01, -8.4204388d+01, -7.5825612d+01, -6.8403549d+01, 
     + -6.1815371d+01 /   
*  values of R for excited-state curves
*  eliminated highly repulsive R = 2.00 values
      data rr /
     +  1.05835,  1.19065,  1.32294,  1.45524,  1.58753, 
     +  1.71983,  1.85212,  1.98441,  2.11671,  2.24900, 
     +  2.38130,  2.51359,  2.64589,  2.77818,  2.91047, 
     +  3.04277,  3.17506,  3.30736,  3.43965,  3.57195, 
     +  3.70424,  3.83653,  3.96883,  4.10112,  4.23342, 
     +  4.36571,  4.49801,  4.63030,  4.76259,  4.89489, 
     +  5.02718,  5.05364,  5.18594,  5.31823,  5.45052, 
     +  5.58282,  5.71511,  5.84741,  5.97970,  6.11200, 
     +  6.24429,  6.37658,  6.50888,  6.64117,  6.77347, 
     +  6.90576,  7.03806,  7.17035,  7.30264,  7.43494, 
     +  7.56723,  7.69953,  7.83182,  7.96412 /  
      data vl4p /
     + 1.5393394e+05, 1.2166936e+05, 9.7800152e+04, 8.0576682e+04, 
     + 6.8249300e+04, 5.9068362e+04, 5.1403585e+04, 4.4324744e+04, 
     + 3.7568858e+04, 3.1685842e+04, 2.6424280e+04, 2.1804497e+04, 
     + 1.7920743e+04, 1.4584723e+04, 1.1584804e+04, 9.7509938e+03, 
     + 7.8335862e+03, 6.3586853e+03, 5.1334088e+03, 4.1197401e+03, 
     + 3.2909826e+03, 2.6205798e+03, 2.0836803e+03, 1.6565706e+03, 
     + 1.3155576e+03, 1.0387715e+03, 8.1499963e+02, 6.3689061e+02, 
     + 4.9709850e+02, 3.8827733e+02, 3.0308115e+02, 2.9368403e+02, 
     + 2.5150354e+02, 2.1622414e+02, 1.8658521e+02, 1.6157933e+02, 
     + 1.4039706e+02, 1.2238463e+02, 1.0701136e+02, 9.3844463e+01, 
     + 8.2529367e+01, 7.2774340e+01, 6.4338346e+01, 5.7021443e+01, 
     + 5.0657149e+01, 4.5106330e+01, 4.0252301e+01, 3.5996874e+01, 
     + 3.2257153e+01, 2.8962938e+01, 2.6054595e+01, 2.3481325e+01, 
     + 2.1199727e+01, 1.9172624e+01  / 
      data vl2s /
     + 1.3284182e+05, 9.4997105e+04, 6.8000849e+04, 4.9826345e+04, 
     + 3.8446890e+04, 3.1835784e+04, 2.7975384e+04, 2.5238874e+04, 
     + 2.2742050e+04, 2.0253046e+04, 1.7756642e+04, 1.5355822e+04, 
     + 1.3157343e+04, 1.1189809e+04, 9.4631439e+03, 7.9638482e+03, 
     + 6.6772355e+03, 5.5780986e+03, 4.6444740e+03, 3.8561140e+03, 
     + 3.1928636e+03, 2.6349423e+03, 2.1671308e+03, 1.7772524e+03, 
     + 1.4531859e+03, 1.1834626e+03, 9.6042972e+02, 7.7781711e+02, 
     + 6.2935640e+02, 5.0877928e+02, 4.0981743e+02, 3.9711092e+02, 
     + 3.4007570e+02, 2.9237194e+02, 2.5229504e+02, 2.1848282e+02, 
     + 1.8984079e+02, 1.6548491e+02, 1.4469762e+02, 1.2689373e+02, 
     + 1.1159379e+02, 9.8403326e+01, 8.6996422e+01, 7.7102721e+01, 
     + 6.8497109e+01, 6.0991454e+01, 5.4427979e+01, 4.8673916e+01, 
     + 4.3617176e+01, 3.9162834e+01, 3.5230259e+01, 3.1750758e+01, 
     + 2.8665648e+01, 2.5924659e+01 / 
      data vl4s /
     + 1.2979041e+05, 9.2298066e+04, 6.5373259e+04, 4.7017079e+04, 
     + 3.5230613e+04, 2.8014950e+04, 2.3384083e+04, 1.9861202e+04, 
     + 1.6808901e+04, 1.4073202e+04, 1.1633224e+04, 9.5131712e+03, 
     + 7.7310305e+03, 6.2551948e+03, 5.0470225e+03, 4.0682097e+03, 
     + 3.2782826e+03, 2.6396873e+03, 2.1234134e+03, 1.7088281e+03, 
     + 1.3757851e+03, 1.1047182e+03, 8.8311744e+02, 7.0318061e+02, 
     + 5.5719112e+02, 4.3794643e+02, 3.4124890e+02, 2.6398965e+02, 
     + 2.0306130e+02, 1.5535648e+02, 1.1776781e+02, 1.1411639e+02, 
     + 9.7726373e+01, 8.4017908e+01, 7.2501149e+01, 6.2784650e+01, 
     + 5.4553887e+01, 4.7554823e+01, 4.1581253e+01, 3.6465010e+01, 
     + 3.2068319e+01, 2.8277822e+01, 2.4999859e+01, 2.2156741e+01, 
     + 1.9683776e+01, 1.7526902e+01, 1.5640779e+01, 1.3987254e+01, 
     + 1.2534116e+01, 1.1254088e+01, 1.0123997e+01, 9.1241050e+00, 
     + 8.2375477e+00, 7.4498793e+00 / 
      data rso /
     +  2.00000,  2.25000,  2.50000,  2.75000,  3.00000, 
     +  3.25000,  3.50000,  3.75000,  4.00000,  4.25000, 
     +  4.50000,  4.75000,  5.00000,  5.25000,  5.50000, 
     +  5.75000,  6.00000,  6.25000,  6.50000,  6.75000, 
     +  7.00000,  7.25000,  7.50000,  7.75000,  8.00000, 
     +  8.25000,  8.50000,  8.75000,  9.00000,  9.25000, 
     +  9.50000 /
      data somat /
     + -2.6178000d+01, -2.5624206d+01, -2.5322333d+01, -2.5194460d+01, 
     + -2.5162667d+01, -2.5162000d+01, -2.5177000d+01, -2.5197333d+01, 
     + -2.5217333d+01, -2.5234771d+01, -2.5249000d+01, -2.5259872d+01, 
     + -2.5268000d+01, -2.5274095d+01, -2.5278503d+01, -2.5281475d+01, 
     + -2.5283265d+01, -2.5284127d+01, -2.5284312d+01, -2.5284074d+01, 
     + -2.5283667d+01, -2.5283298d+01, -2.5282999d+01, -2.5282759d+01, 
     + -2.5282564d+01, -2.5282403d+01, -2.5282264d+01, -2.5282134d+01, 
     + -2.5282000d+01, -2.5281853d+01, -2.5281692d+01, 
     + -1.7460000d+00, -3.5041328d+00, -5.5650000d+00, -8.0248672d+00, 
     + -1.0980000d+01, -1.4480000d+01, -1.8361000d+01, -2.2287000d+01, 
     + -2.5845000d+01, -2.8750040d+01, -3.0962000d+01, -3.2528256d+01, 
     + -3.3624000d+01, -3.4427786d+01, -3.5003801d+01, -3.5387640d+01, 
     + -3.5614896d+01, -3.5721166d+01, -3.5742043d+01, -3.5713123d+01, 
     + -3.5670000d+01, -3.5641725d+01, -3.5631179d+01, -3.5634696d+01, 
     + -3.5648613d+01, -3.5669266d+01, -3.5692990d+01, -3.5716123d+01, 
     + -3.5735000d+01, -3.5746650d+01, -3.5750871d+01, 
     + -5.7509000d+01, -5.6802859d+01, -5.6806000d+01, -5.7294141d+01, 
     + -5.8043000d+01, -5.8849000d+01, -5.9598000d+01, -6.0231000d+01, 
     + -6.0732000d+01, -6.1111844d+01, -6.1391000d+01, -6.1591449d+01, 
     + -6.1736000d+01, -6.1844550d+01, -6.1924531d+01, -6.1980254d+01, 
     + -6.2016035d+01, -6.2036186d+01, -6.2045022d+01, -6.2046855d+01, 
     + -6.2046000d+01, -6.2045998d+01, -6.2047306d+01, -6.2049608d+01, 
     + -6.2052589d+01, -6.2055935d+01, -6.2059328d+01, -6.2062455d+01, 
     + -6.2065000d+01, -6.2066709d+01, -6.2067576d+01, 
     + -3.9973000d+01, -4.0247413d+01, -4.0802000d+01, -4.1482087d+01, 
     + -4.2133000d+01, -4.2639000d+01, -4.3027000d+01, -4.3313000d+01, 
     + -4.3512000d+01, -4.3648057d+01, -4.3738000d+01, -4.3795743d+01, 
     + -4.3833000d+01, -4.3859289d+01, -4.3877551d+01, -4.3889078d+01, 
     + -4.3895165d+01, -4.3897107d+01, -4.3896197d+01, -4.3893730d+01, 
     + -4.3891000d+01, -4.3889057d+01, -4.3887975d+01, -4.3887583d+01, 
     + -4.3887710d+01, -4.3888186d+01, -4.3888840d+01, -4.3889502d+01, 
     + -4.3890000d+01, -4.3890197d+01, -4.3890081d+01, 
     + -2.3533000d+01, -3.0389322d+01, -3.5561000d+01, -3.9798428d+01, 
     + -4.3852000d+01, -4.8294000d+01, -5.2979000d+01, -5.7537000d+01, 
     + -6.1526000d+01, -6.4678410d+01, -6.7005000d+01, -6.8606145d+01, 
     + -6.9697000d+01, -7.0486843d+01, -7.1046669d+01, -7.1412902d+01, 
     + -7.1621965d+01, -7.1710283d+01, -7.1714279d+01, -7.1670377d+01, 
     + -7.1615000d+01, -7.1577802d+01, -7.1561351d+01, -7.1561445d+01, 
     + -7.1573881d+01, -7.1594457d+01, -7.1618970d+01, -7.1643219d+01, 
     + -7.1663000d+01, -7.1674902d+01, -7.1678675d+01, 
     + -6.8146000d+01, -6.8684892d+01, -6.9187000d+01, -6.9449858d+01, 
     + -6.9271000d+01, -6.8442000d+01, -6.6779000d+01, -6.4343000d+01, 
     + -6.1500000d+01, -5.8736533d+01, -5.6374000d+01, -5.4608382d+01, 
     + -5.3327000d+01, -5.2365232d+01, -5.1659336d+01, -5.1170789d+01, 
     + -5.0861069d+01, -5.0691652d+01, -5.0624017d+01, -5.0619640d+01, 
     + -5.0640000d+01, -5.0653406d+01, -5.0655499d+01, -5.0648753d+01, 
     + -5.0635641d+01, -5.0618636d+01, -5.0600212d+01, -5.0582842d+01, 
     + -5.0569000d+01, -5.0560684d+01, -5.0557995d+01, 
     + 4.8837000d+01, 5.8669650d+01, 6.4021000d+01, 6.6035600d+01, 
     + 6.5858000d+01, 6.4428000d+01, 6.1960000d+01, 5.8839000d+01, 
     + 5.5554000d+01, 5.2559885d+01, 5.0074000d+01, 4.8198276d+01, 
     + 4.6817000d+01, 4.5781239d+01, 4.5026825d+01, 4.4510783d+01, 
     + 4.4190138d+01, 4.4021912d+01, 4.3963132d+01, 4.3970819d+01, 
     + 4.4002000d+01, 4.4021440d+01, 4.4024873d+01, 4.4015775d+01, 
     + 4.3997624d+01, 4.3973894d+01, 4.3948063d+01, 4.3923606d+01, 
     + 4.3904000d+01, 4.3892055d+01, 4.3887919d+01 / 
*
*  set vv0 term to zero
      vv0 = 0.d0
*
*  spline fits
      if (ifirst .eq. 0) then
*  spline fit of the coefficients for the PE curves
*  evaluate derivative at first point
        der1=(vlx(2)-vlx(1))/(rrx(2)-rrx(1))
        call dspline(rrx,vlx,53,der1,0d0,csx)
*
        der1=(vl4p(2)-vl4p(1))/(rr(2)-rr(1))
        call dspline(rr,vl4p,54,der1,0d0,cs4p)
*
        der1=(vl2s(2)-vl2s(1))/(rr(2)-rr(1))
        call dspline(rr,vl2s,54,der1,0d0,cs2s)
*
        der1=(vl4s(2)-vl4s(1))/(rr(2)-rr(1))
        call dspline(rr,vl4s,54,der1,0d0,cs4s)
*
*  spline fit of the coefficients for the spin-orbit matrix elements
        indso=1
        do ilam=5,11
          call dcopy(31,somat(indso),1,vecso,1)
*  evaluate derivative at first point
          der1=(vecso(2)-vecso(1))/(rso(2)-rso(1))
          call dspline(rso,vecso,31,der1,0d0,csplso(1,ilam))
          indso = indso + 31
        enddo
        ifirst = 1
      end if
*
*  determine splined coefficients for PE curves at r=R
      call dsplint(rrx,vlx,csx,53,r,vvx)
      v(1) = vvx
*  fix long-range behavior of PE curve as c6/r^6 from last point
      if (r.gt.rrx(53)) then
        c6 = vlx(53) * rrx(53)**6
        v(1) = c6 / r**6
      endif
*
      call dsplint(rr,vl4p,cs4p,54,r,vvx)
      v(2) = vvx
*  fix long-range behavior of PE curve as c6/r^6 from last point
      if (r.gt.rr(54)) then
        c6 = vl4p(54) * rr(54)**6
        v(2) = c6 / r**6
      endif
*
      call dsplint(rr,vl2s,cs2s,54,r,vvx)
      v(3) = vvx
*  fix long-range behavior of PE curve as c6/r^6 from last point
      if (r.gt.rr(54)) then
        c6 = vl2s(54) * rr(54)**6
        v(3) = c6 / r**6
      endif
*
      call dsplint(rr,vl4s,cs4s,54,r,vvx)
      v(4) = vvx
*  fix long-range behavior of PE curve as c6/r^6 from last point
      if (r.gt.rr(54)) then
        c6 = vl4s(54) * rr(54)**6
        v(4) = c6 / r**6
      endif
*
*  determine splined coefficients for spin-orbit matrix elements at r=R
      indso=1
      do ilam=5,11
        call dcopy(31,somat(indso),1,vecso,1)
        call dsplint(rso,vecso,csplso(1,ilam),31,r,vvx)
        v(ilam) = vvx
        call dcopy(31,somat(indso),1,vecso,1)
        indso = indso + 31
*  if r≥8.5 bohr, set equal to SO matrix element at R=9.5
        if (r.ge.9.5d0) v(ilam) = somat(indso-1)
*  get value of spin-orbit parameter at R=inf
        v(ilam+7) = somat(indso-1)
      enddo
      call dcopy(18,v,1,vvl,1)
* convert to hartree
      econv=1.d0/219474.6d0
      call dscal(18,econv,vvl,1)
      return
      end
*===================================eof====================================