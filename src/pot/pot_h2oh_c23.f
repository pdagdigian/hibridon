*  system:  H2O-H, PES
*  written by P. Dagdigian, Feb 2013
*
*  H2O defined to lie on the xz plane, with the origin at center of mass
*  x axis is the C2 symmetry axis of the CH2 molecule
*  z axis is the a inertial axis of the molecule (perpendicular to C2 axis)
*  when theta = 90, phi = 0, He is on C side of molecule
*  phi = 0 has all 4 atoms coplanar
*
*  the PES is fitted with 23 angular terms

      include "common/syusr"
      include "common/ground"
      include "common/bausr"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(22)
      include "common/parpot"
      s4pi = sqrt ( 4.d0 * acos(-1.d0) )
      potnam='H2O-H CCSD(T) PES-23 coeffs'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      call pot(vv0,r)
*     vlm coefficient is returned in atomic units (hartree)
*     convert from atomic units for printout
      econv=219474.6d0
      write (6, 100) vv0*econv*s4pi, (econv*vvl(i), i=1,22)
100   format(' v(lam,0):',5(1pe16.8)/' v(lam,1):',4(1pe16.8)/
     :    ' v(lam,2):',4(1pe16.8)/' v(lam,3):',3(1pe16.8)/
     :    ' v(lam,4):',3(1pe16.8)/' v(lam,5):',2(1pe16.8)/
     :    ' v(lam,6):',2(1pe16.8))
99    goto 1
      end
* ------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
      implicit double precision (a-h,o-z)
      character*(*) filnam
      include "common/parbas"
      include "common/parpot" 
      common /conlam/ nlam, nlammx, lamnum(2)
      common /cosysi/ nscode, isicod, nterm
      potnam='H2O-H CCSD(T) PES-23 coeffs'
      ibasty = 16
*
      nterm = 7
      do i=1,nterm
        mproj(i)=i-1
      enddo
      lammin(1) = 2
      do i=2,nterm
        lammin(i)=i-1
      enddo
      lammax(1) = 8
      lammax(2) = 7
      lammax(3) = 8
      lammax(4) = 7
      lammax(5) = 8
      lammax(6) = 7
      lammax(7) = 8
* 
      ipotsy = 2
*
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
*  calculate total number of anisotropic terms
      nlam = 0
      do i=1,nterm
        lmin = lammin(i)
        lmax = lammax(i)
        lamnum(i)=(lammax(i)-lammin(i))/2+1
        do lb = lmin, lmax, ipotsy
          nlam = nlam + 1
        end do
      end do
      nlammx = nlam
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:          interparticle distance
*  on return:
*    vv0         contains isotropic term (Y00)
*  variable in common block /covvl/
*    vvl:        vector of length 22 to store r-dependence of each term
*                in potential expansion
*    vvl(1-4):   expansion coefficients in Yl0 (l=2:8:2) of v(lam,0)
*    vvl(5-9):   expansion coefficients in [Yl1 - Y(l,-1)] (l=1:7:2) of v(lam,1)
*    vvl(9-12):  expansion coefficients in [Yl2 + Y(l,-2)] (l=2:8:2) of v(lam,2)
*    vvl(13-15): expansion coefficients in [Yl3 - Y(l,-3)] (l=3:7:2) of v(lam,3)
*    vvl(16-18): expansion coefficients in [Yl4 + Y(l,-4)] (l=4:8:2) of v(lam,4)
*    vvl(19-20): expansion coefficients in [Yl5 - Y(l,-5)] (l=5:7:2) of v(lam,5)
*    vvl(21-22): expansion coefficients in [Yl6 + Y(l,-6)] (l=6:8:2) of v(lam,6)
*  variable in common block /coloapot/
*    s4pi:       normalization factor for isotropic potential
*
*  uses linear least squares routines from lapack
*
* author:  paul dagdigian
* latest revision date:  20-feb-2013
* ----------------------------------------------------------------------
* this pes is a fit to 19 values of R and 190 orientations, namely
* R = [3:0.5:10 11 12 13 15 20]
* theta=[0:10:90] and phi=[0:10:180]
*
* in this version, 23 terms are employed in the angular fit
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      dimension v(23)
      dimension csplin(20,23)
      dimension rr(20), vl(460),vec(23)
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(22)
* hyperbolic tangent scaling factor
      data alph /1.2d0/
      data rmax /13d0/
* sqrt(4*pi)
      s4pi = 3.544907701811032d0
* here are the 20 values of R
      data rr/3d0,3.5d0,4d0,4.5d0,5d0,5.5d0,6d0,6.5d0,7d0,7.5d0,
     +  8d0,8.5d0,9d0,9.5d0,10d0,11d0,12d0,13d0,15d0,20d0/
* here are column ordered legendre expansion coefficients (12 total) for the
* potential at each of 20 values of R (460 values)
      data vl/
     + 2.3867375e+04, 1.0479715e+04, 4.3951279e+03, 1.6145198e+03,
     + 4.2166277e+02, -2.8575651e+01, -1.5812633e+02, -1.6545872e+02,
     + -1.3579671e+02, -1.0160035e+02, -7.3038153e+01, -5.1733103e+01,
     + -3.6605805e+01, -2.6046142e+01, -1.8711374e+01, -9.9899304e+00,
     + -5.5754766e+00, -3.2365744e+00, -1.1437447e+00, 1.3908087e-01,
     + 6.6144371e+03, 2.4988441e+03, 9.7053659e+02, 3.4383585e+02,
     + 9.8475820e+01, 1.4119269e+01, -8.5302228e+00, -1.0801998e+01,
     + -7.9940602e+00, -5.0029987e+00, -2.8908646e+00, -1.5925657e+00,
     + -8.4590724e-01, -4.4428663e-01, -2.3518792e-01, -5.7041624e-02,
     + -1.7902631e-03, 1.6136017e-02, -3.3315009e-03, 5.1830847e-03,
     + 3.2407288e+01, 1.4313296e+02, 5.2340311e+01, 8.6488061e+00,
     + -2.3010777e+00, -2.8655893e+00, -1.5900299e+00, -5.9217405e-01,
     + -7.9517415e-02, 1.0899866e-01, 1.4862420e-01, 1.2706809e-01,
     + 9.9008282e-02, 6.9986699e-02, 5.5779115e-02, 2.2020279e-02,
     + 1.1286601e-02, 6.5722204e-03, 1.0923885e-02, 1.3719462e-03,
     + -7.0014807e+02, -9.7641184e+01, -1.3900238e+01, -2.0822736e+00,
     + -2.3003901e-02, 3.2274767e-01, 3.2455958e-01, 2.5735370e-01,
     + 1.7554233e-01, 1.1181826e-01, 6.7247136e-02, 4.2841409e-02,
     + 2.4081053e-02, 1.5463959e-02, 3.5650117e-03, 4.1080196e-03,
     + 1.3703501e-03, -5.6965896e-04, -5.7938896e-03, 8.9070078e-03,
     + -7.5160667e+01, -1.7672336e+01, -4.7086773e+00, -1.1618813e+00,
     + -2.8313182e-01, -5.2857892e-02, -9.3385640e-03, 7.4772400e-03,
     + 6.8283140e-03, 9.9578926e-04, 1.7680098e-03, -2.6008230e-03,
     + 4.0583296e-04, 1.9599226e-03, 5.5880291e-03, 5.7007456e-04,
     + 1.4903037e-03, 4.3607106e-03, 2.1879985e-03, -2.8998749e-03,
     + 4.1318121e+03, 1.6085626e+03, 7.3812846e+02, 3.3415415e+02,
     + 1.3544685e+02, 4.4566125e+01, 7.7529645e+00, -4.5671856e+00,
     + -7.0855004e+00, -6.2685694e+00, -4.7302684e+00, -3.3346720e+00,
     + -2.2727499e+00, -1.5474607e+00, -1.0655658e+00, -5.1278969e-01,
     + -2.5163274e-01, -1.3330513e-01, -4.3146729e-02, 4.6673087e-03,
     + 3.7087877e+03, 1.2649136e+03, 5.3642480e+02, 2.2791952e+02,
     + 8.5555713e+01, 2.4581686e+01, 1.7963106e+00, -4.8564335e+00,
     + -5.5342635e+00, -4.4772670e+00, -3.2206827e+00, -2.2046063e+00,
     + -1.4809664e+00, -9.8631025e-01, -6.5612628e-01, -3.1758995e-01,
     + -1.5953254e-01, -8.0719982e-02, -4.1053355e-02, -1.1474939e-03,
     + 7.1129597e+02, 1.8187546e+02, 5.4795759e+01, 1.6390136e+01,
     + 3.9953844e+00, 3.1057746e-01, -5.4392610e-01, -5.7108911e-01,
     + -4.0913770e-01, -2.6021866e-01, -1.5638533e-01, -9.1809400e-02,
     + -5.3729539e-02, -3.1046848e-02, -2.4565258e-02, -6.2641365e-03,
     + -3.8931057e-03, 2.6677535e-03, 6.0687671e-03, -7.2316023e-03,
     + -3.3140048e+02, -3.1858161e+01, -2.7265012e-01, 1.0575406e+00,
     + 5.2804574e-01, 2.1574238e-01, 9.6777199e-02, 4.9587430e-02,
     + 2.8682795e-02, 1.5640487e-02, 9.6514051e-03, 3.6564026e-03,
     + 3.2989783e-03, 2.2264487e-03, 5.8459514e-03, 5.8244868e-04,
     + 1.7819517e-04, -1.6758482e-03, -3.4946417e-03, 1.6381672e-03,
     + 2.3603213e+03, 9.1956932e+02, 3.8110188e+02, 1.4678554e+02,
     + 4.6523571e+01, 8.0344247e+00, -4.1429274e+00, -6.3384433e+00,
     + -5.4315294e+00, -3.9934178e+00, -2.7746114e+00, -1.8931521e+00,
     + -1.2997050e+00, -8.9954937e-01, -6.2820469e-01, -3.4634533e-01,
     + -2.0456066e-01, -1.1742879e-01, -3.7317996e-02, -1.1167013e-02,
     + 1.7651722e+03, 4.0263585e+02, 1.2968562e+02, 4.8113703e+01,
     + 1.5490981e+01, 2.6171637e+00, -1.5845453e+00, -2.3490168e+00,
     + -2.0004031e+00, -1.4454973e+00, -9.7200552e-01, -6.3941876e-01,
     + -4.1160668e-01, -2.6235547e-01, -1.7427930e-01, -7.1491921e-02,
     + -3.4449119e-02, -1.5357808e-02, -6.4412183e-03, -3.0976554e-03,
     + 8.0260531e+02, 1.4612188e+02, 3.1594921e+01, 7.5773471e+00,
     + 1.4996272e+00, -4.1878101e-02, -3.5034750e-01, -3.0302744e-01,
     + -2.0968750e-01, -1.3379688e-01, -7.9425075e-02, -4.3061825e-02,
     + -2.7065484e-02, -1.5478143e-02, -3.9482310e-03, -3.3400756e-03,
     + -2.0777968e-03, 2.5295564e-03, 2.2151412e-03, -3.9439540e-03,
     + 4.3378345e+00, 8.5114230e+00, 3.6330081e+00, 1.1372589e+00,
     + 3.0485307e-01, 7.7492023e-02, 6.3973638e-03, -2.0211957e-03,
     + -2.8977408e-03, -2.5228111e-03, -1.2018287e-03, -1.0939482e-04,
     + 1.1883665e-03, -3.2855831e-04, -3.5414118e-03, 1.2791217e-04,
     + 6.1805992e-04, 4.9879771e-04, -9.6092562e-04, -5.9609304e-03,
     + 6.6443637e+02, 2.3779548e+02, 1.0250167e+02, 4.3179929e+01,
     + 1.5475158e+01, 3.7319484e+00, -5.3964045e-01, -1.7321751e+00,
     + -1.6636712e+00, -1.2942493e+00, -9.2547276e-01, -6.4382725e-01,
     + -4.3996421e-01, -3.0012452e-01, -2.1540476e-01, -1.0033345e-01,
     + -5.4164876e-02, -3.2118187e-02, -7.5207439e-03, 3.0562332e-04,
     + 9.4470158e+02, 1.7392752e+02, 4.0954386e+01, 1.1767839e+01,
     + 2.9498794e+00, 8.1052843e-02, -6.5321096e-01, -7.0442033e-01,
     + -5.3763355e-01, -3.6451853e-01, -2.3153738e-01, -1.3873534e-01,
     + -8.9307011e-02, -5.5352678e-02, -3.0540208e-02, -1.3906493e-02,
     + -7.1642893e-03, -3.6441367e-03, 1.2685766e-03, -1.7293759e-02,
     + 6.0631609e+02, 9.5564732e+01, 1.7062171e+01, 3.3796136e+00,
     + 5.6738258e-01, -3.2608404e-02, -1.3414037e-01, -1.0426671e-01,
     + -6.8784613e-02, -4.2201496e-02, -2.4772747e-02, -1.6664584e-02,
     + -7.3015300e-03, -3.4292499e-03, -4.8099010e-03, -3.3711061e-04,
     + 1.2168899e-03, 1.6051002e-03, 2.0214925e-03, -7.9619821e-03,
     + 2.3610838e+02, 6.8438367e+01, 2.4872522e+01, 9.5105161e+00,
     + 3.1663071e+00, 7.0070504e-01, -9.9882355e-02, -3.3554074e-01,
     + -3.0862277e-01, -2.2970949e-01, -1.5770199e-01, -1.0727775e-01,
     + -6.9821901e-02, -4.8240370e-02, -2.6881318e-02, -1.2625537e-02,
     + -2.3001480e-03, -5.5615629e-03, -2.6148588e-03, -3.9044855e-03,
     + 4.6632954e+02, 8.0393171e+01, 1.6227606e+01, 3.8056679e+00,
     + 7.6147137e-01, -3.1900138e-02, -1.6246464e-01, -1.6859964e-01,
     + -1.2079007e-01, -7.7966956e-02, -4.9221288e-02, -2.9169423e-02,
     + -1.6078986e-02, -9.2532978e-03, -9.1991606e-03, -2.2266989e-03,
     + 2.5763837e-03, -3.0006375e-03, -3.4084504e-03, 1.5090299e-03,
     + 3.5531882e+02, 5.1553790e+01, 8.4082272e+00, 1.5198844e+00,
     + 2.4990420e-01, -2.9017836e-03, -2.8672767e-02, -2.8431270e-02,
     + -1.7265629e-02, -1.1056356e-02, -6.4107717e-03, -5.8746530e-03,
     + -7.4480816e-04, 1.5321225e-03, 2.4354447e-03, -6.1041435e-04,
     + -9.7367902e-04, -5.5351998e-04, -2.4152212e-03, -3.3283992e-03,
     + 1.0701559e+02, 2.3766096e+01, 7.0452895e+00, 2.3840817e+00,
     + 8.5791402e-01, 2.5640048e-01, -1.8839603e-03, -4.8311379e-02,
     + -4.9662348e-02, -3.7551719e-02, -2.6707778e-02, -1.0703634e-02,
     + -1.0243543e-02, -7.3219718e-03, -2.0609179e-03, -1.5247441e-03,
     + -5.4245599e-03, -2.2223820e-03, -5.2300404e-03, -2.8217197e-03,
     + 2.1403989e+02, 3.4219824e+01, 6.2419511e+00, 1.2670130e+00,
     + 2.5356407e-01, 1.6106472e-02, -4.1992719e-02, -3.4823626e-02,
     + -2.4178959e-02, -1.4910869e-02, -9.0848876e-03, -6.2898248e-03,
     + -1.6305706e-03, -4.7751728e-04, 1.1513493e-03, -4.2390450e-04,
     + -2.2539002e-03, -2.3446758e-03, -2.7792986e-03, -1.8797771e-03,
     + 4.2545141e+01, 8.6083920e+00, 2.2553202e+00, 6.9477115e-01,
     + 1.9025420e-01, 3.4134276e-02, 3.8359364e-04, -5.1967873e-03,
     + -6.8642916e-03, -5.1725933e-03, -4.0068661e-03, -4.2400736e-03,
     + 1.4712544e-04, -2.1426727e-03, -4.9189979e-03, 9.3447261e-04,
     + 4.1814715e-03, -1.6082540e-03, -5.8795160e-03, -2.6006740e-03,
     + 8.7539175e+01, 1.3110208e+01, 2.2870035e+00, 4.2492265e-01,
     + 7.6641761e-02, -6.0607281e-03, -1.1398397e-02, -6.2043831e-03,
     + -4.0193799e-03, -2.2967645e-03, -2.2670020e-03, -1.5252144e-03,
     + -3.0883808e-04, 5.7068790e-04, 3.0427554e-03, 1.3447144e-03,
     + 3.1152965e-04, -9.0832410e-04, -4.0714277e-03, -5.6433540e-03 /
*
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         ind=1
         do ilam=1,23
           call dcopy(23,vl(ind),1,vec,1)
*    evaluate derivative at first point
           der1=(vec(2)-vec(1))/(rr(2)-rr(1))
           call dspline(rr,vec,20,der1,0d0,csplin(1,ilam))
           ind = ind + 20
         enddo
         ifirst = 1
       end if
* r^-6 fit to isotropic part of potential
       c6sum = -1.9457075e+07
* switching function for long-range
       switch_lr=0.5*(tanh(0.5*(r - 25.d0)) + 1.d0)
* determine splined coefficients at R=r
       ind=1
       do ilam=1,23
         call dcopy(20,vl(ind),1,vec,1)
         call dsplint(rr,vec,csplin(1,ilam),20,r,vvx)
* kill anisotropic terms at large R
         vvx = (1.d0 - switch_lr)*vvx
         if (ilam.eq.1) then
* merge with asymptotic form
            vvx = vvx + switch_lr*c6sum/(r**6)
         endif
         v(ilam)=vvx
         call dcopy(20,vl(ind),1,vec,1)
         ind = ind + 20
       enddo
       call dcopy(22,v(2),1,vvl(1),1)
* convert to hartree
       econv=1.d0/219474.6d0
       call dscal(22,econv,vvl,1)
* isotropic term - divide by sqrt(4*pi)
       vv0 = v(1)*econv/s4pi
*
       return
       end
