*  system:  H2O-H, PES
*  written by P. Dagdigian, Feb 2013
*
*  H2O defined to lie on the xz plane, with the origin at center of mass
*  x axis is the C2 symmetry axis of the CH2 molecule
*  z axis is the a inertial axis of the molecule (perpendicular to C2 axis)
*  when theta = 90, phi = 0, He is on C side of molecule
*  phi = 0 has all 4 atoms coplanar
*
*  the PES is fitted with 34 angular terms

      include "common/syusr"
      include "common/ground"
      include "common/bausr"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(33)
      include "common/parpot"
      s4pi = sqrt ( 4.d0 * acos(-1.d0) )
      econv=219474.6d0
      potnam='H2O-H CCSD(T) PES-34 coeffs'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      call pot(vv0,r)
      if (r.le.0.d0) goto 99
*     vlm coefficient is returned in atomic units (hartree)
*     convert from atomic units for printout
      write (6, 100) vv0*econv*s4pi, (econv*vvl(i), i=1,33)
100   format(' v(lam,0):',6(1pe16.8)/' v(lam,1):',5(1pe16.8)/
     :    ' v(lam,2):',5(1pe16.8)/' v(lam,3):',4(1pe16.8)/
     :    ' v(lam,4):',4(1pe16.8)/' v(lam,5):',3(1pe16.8)/
     :    ' v(lam,6):',3(1pe16.8)/' v(lam,7):',2(1pe16.8)/
     :    ' v(lam,8):',2(1pe16.8))
      goto 1
99    rr=3.5d0
      dr=0.5d0
      open (unit=12,file='h2oh_c34_vlms.dat')
      write(12,109)
109   format(' %R/bohr V00  ...')
      do i=1,40
        call pot(vv0,rr)
        write (12,110) rr,vv0*econv*s4pi, (econv*vvl(j),j=1,33)
110     format(f7.2,34(1pe16.8))
        rr = rr + dr
      enddo
      close(12)
      end
* ------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
      implicit double precision (a-h,o-z)
      character*(*) filnam
      include "common/parbas"
      include "common/parpot" 
      common /conlam/ nlam, nlammx, lamnum(2)
      common /cosysi/ nscode, isicod, nterm
      potnam='H2O-H CCSD(T) PES-34 coeffs'
      ibasty = 16
*
      nterm = 9
      do i=1,nterm
        mproj(i)=i-1
      enddo
      lammin(1) = 2
      do i=2,nterm
        lammin(i)=i-1
      enddo
      lammax(1) = 10
      lammax(2) = 9
      lammax(3) = 10
      lammax(4) = 9
      lammax(5) = 10
      lammax(6) = 9
      lammax(7) = 10
      lammax(8) = 9
      lammax(9) = 10
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
*    vvl:        vector of length 33 to store r-dependence of each term
*                in potential expansion
*    vvl(1-5):   expansion coefficients in Yl0 (l=2:10:2) of v(lam,0)
*    vvl(6-10):  expansion coefficients in [Yl1 - Y(l,-1)] (l=1:9:2) of v(lam,1)
*    vvl(11-15): expansion coefficients in [Yl2 + Y(l,-2)] (l=2:10:2) of v(lam,2)
*    vvl(16-19): expansion coefficients in [Yl3 - Y(l,-3)] (l=3:9:2) of v(lam,3)
*    vvl(20-23): expansion coefficients in [Yl4 + Y(l,-4)] (l=4:10:2) of v(lam,4)
*    vvl(24-26): expansion coefficients in [Yl5 - Y(l,-5)] (l=5:9:2) of v(lam,5)
*    vvl(27-29): expansion coefficients in [Yl6 + Y(l,-6)] (l=6:10:2) of v(lam,6)
*    vvl(30-31): expansion coefficients in [Yl7 - Y(l,-7)] (l=5:9:2) of v(lam,7)
*    vvl(32-33): expansion coefficients in [Yl8 + Y(l,-8)] (l=6:10:2) of v(lam,8)
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
* in this version, 34 terms are employed in the angular fit
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      dimension v(34)
      dimension csplin(20,34)
      dimension rr(20), vl(680),vec(34)
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(33)
* hyperbolic tangent scaling factor
      data alph /1.2d0/
      data rmax /13d0/
* sqrt(4*pi)
      s4pi = 3.544907701811032d0
* here are the 20 values of R
      data rr/3d0,3.5d0,4d0,4.5d0,5d0,5.5d0,6d0,6.5d0,7d0,7.5d0,
     +  8d0,8.5d0,9d0,9.5d0,10d0,11d0,12d0,13d0,15d0,20d0/
* here are column ordered legendre expansion coefficients (12 total) for the
* potential at each of 20 values of R (680 values)
      data vl/
     + 2.3870357e+04, 1.0479954e+04, 4.3951288e+03, 1.6145148e+03,
     + 4.2166022e+02, -2.8578999e+01, -1.5812621e+02, -1.6545889e+02,
     + -1.3579678e+02, -1.0160047e+02, -7.3038058e+01, -5.1733076e+01,
     + -3.6605774e+01, -2.6046065e+01, -1.8711015e+01, -9.9899932e+00,
     + -5.5752958e+00, -3.2362363e+00, -1.1431448e+00, 1.3926819e-01,
     + 6.5945741e+03, 2.4966265e+03, 9.7027462e+02, 3.4379864e+02,
     + 9.8474861e+01, 1.4122315e+01, -8.5300887e+00, -1.0801656e+01,
     + -7.9939589e+00, -5.0027685e+00, -2.8909459e+00, -1.5925829e+00,
     + -8.4591172e-01, -4.4430255e-01, -2.3538002e-01, -5.7025197e-02,
     + -2.0529232e-03, 1.5552458e-02, -3.9549155e-03, 5.2946772e-03,
     + 3.2815318e+01, 1.4313960e+02, 5.2329497e+01, 8.6449818e+00,
     + -2.3027417e+00, -2.8665969e+00, -1.5900439e+00, -5.9224643e-01,
     + -7.9491003e-02, 1.0897188e-01, 1.4862921e-01, 1.2688137e-01,
     + 9.8954212e-02, 6.9874235e-02, 5.6185584e-02, 2.1880592e-02,
     + 1.1184641e-02, 6.6277472e-03, 1.0765280e-02, 1.2969071e-03,
     + -7.2211276e+02, -9.9956270e+01, -1.4133640e+01, -2.1137174e+00,
     + -1.8205421e-02, 3.2687589e-01, 3.2306115e-01, 2.5784282e-01,
     + 1.7573417e-01, 1.1212496e-01, 6.7057125e-02, 4.2626516e-02,
     + 2.3933369e-02, 1.5280690e-02, 3.5828045e-03, 3.9623435e-03,
     + 9.9993524e-04, -1.3739419e-03, -7.3574357e-03, 8.1916486e-03,
     + -1.1888105e+02, -2.2304172e+01, -5.1818132e+00, -1.2225677e+00,
     + -2.7460755e-01, -4.5955290e-02, -1.0770685e-02, 8.3378611e-03,
     + 7.1265903e-03, 1.6400504e-03, 1.4025022e-03, -2.8942577e-03,
     + 2.7901035e-04, 1.6354373e-03, 6.1549555e-03, 7.3970663e-05,
     + 5.3047931e-04, 2.6264147e-03, -3.8402808e-05, -3.3738083e-03,
     + 9.0726065e+01, 9.7785138e+00, 1.0500586e+00, 1.4362767e-01,
     + -1.0902452e-02, -1.2520244e-02, 3.0137675e-03, -1.6615099e-03,
     + -6.3400261e-04, -1.1688460e-03, 6.3081264e-04, 7.3046530e-04,
     + 3.1380228e-04, 6.8613384e-04, -1.1225914e-03, 8.9348601e-04,
     + 1.7389320e-03, 3.0930305e-03, 4.7299319e-03, 1.2608656e-03,
     + 4.1492296e+03, 1.6107937e+03, 7.3844131e+02, 3.3420189e+02,
     + 1.3544839e+02, 4.4566678e+01, 7.7523494e+00, -4.5674650e+00,
     + -7.0856025e+00, -6.2684751e+00, -4.7302867e+00, -3.3347077e+00,
     + -2.2728072e+00, -1.5474925e+00, -1.0656438e+00, -5.1275317e-01,
     + -2.5163552e-01, -1.3342274e-01, -4.2522217e-02, 4.4165250e-03,
     + 3.6856244e+03, 1.2618417e+03, 5.3596531e+02, 2.2784386e+02,
     + 8.5545504e+01, 2.4581164e+01, 1.7966076e+00, -4.8560144e+00,
     + -5.5339676e+00, -4.4772209e+00, -3.2205903e+00, -2.2048458e+00,
     + -1.4808433e+00, -9.8633122e-01, -6.5597340e-01, -3.1762592e-01,
     + -1.5957772e-01, -8.0612061e-02, -4.0785575e-02, -5.5381518e-04,
     + 7.3733334e+02, 1.8526100e+02, 5.5291350e+01, 1.6469287e+01,
     + 4.0038441e+00, 3.0772529e-01, -5.4269335e-01, -5.7151348e-01,
     + -4.0914623e-01, -2.6012936e-01, -1.5631953e-01, -9.1922564e-02,
     + -5.3781701e-02, -3.0911229e-02, -2.4291195e-02, -6.2491982e-03,
     + -3.8801269e-03, 2.3382587e-03, 6.3338413e-03, -7.3901053e-03,
     + -2.6745497e+02, -2.3508083e+01, 9.5170275e-01, 1.2545119e+00,
     + 5.5111587e-01, 2.1471438e-01, 9.8447168e-02, 4.8429832e-02,
     + 2.8189142e-02, 1.5506070e-02, 9.4849351e-03, 3.4948617e-03,
     + 3.1384361e-03, 2.3781437e-03, 7.0315726e-03, 6.5829819e-04,
     + -5.9477198e-05, -1.9705594e-03, -3.4318202e-03, 1.2452795e-03,
     + -1.7574744e+02, -2.3255657e+01, -3.5026115e+00, -5.8547197e-01,
     + -7.4358817e-02, 6.0325522e-04, -4.5605949e-03, 3.3565470e-03,
     + 1.2159983e-03, 2.0381438e-04, 2.1688137e-04, 1.1865142e-03,
     + 2.5067490e-04, -8.5604575e-04, -3.4206942e-03, 1.3866986e-04,
     + 5.4104171e-04, 8.0354047e-04, -1.4751721e-03, 1.1733999e-04,
     + 2.3714506e+03, 9.2078884e+02, 3.8123332e+02, 1.4679735e+02,
     + 4.6524992e+01, 8.0315942e+00, -4.1437639e+00, -6.3385923e+00,
     + -5.4315640e+00, -3.9934560e+00, -2.7746284e+00, -1.8933581e+00,
     + -1.2998894e+00, -8.9947708e-01, -6.2797599e-01, -3.4641741e-01,
     + -2.0448616e-01, -1.1739129e-01, -3.6901815e-02, -1.1147050e-02,
     + 1.7518333e+03, 4.0109389e+02, 1.2949275e+02, 4.8087752e+01,
     + 1.5488016e+01, 2.6177690e+00, -1.5837300e+00, -2.3488959e+00,
     + -2.0004232e+00, -1.4454854e+00, -9.7193186e-01, -6.3912487e-01,
     + -4.1140523e-01, -2.6226567e-01, -1.7442935e-01, -7.1403531e-02,
     + -3.4361278e-02, -1.5310903e-02, -6.2917742e-03, -2.8301339e-03,
     + 8.1274072e+02, 1.4719449e+02, 3.1696528e+01, 7.5856049e+00,
     + 1.4986681e+00, -4.3030562e-02, -3.5030746e-01, -3.0314467e-01,
     + -2.0970223e-01, -1.3375940e-01, -7.9470390e-02, -4.3304137e-02,
     + -2.7254938e-02, -1.5396753e-02, -3.9965403e-03, -3.3150781e-03,
     + -2.1723126e-03, 2.3472013e-03, 2.6784601e-03, -3.5086593e-03,
     + 2.9268637e+01, 1.1260233e+01, 3.9330896e+00, 1.1691341e+00,
     + 3.0619488e-01, 7.6634138e-02, 5.0688230e-03, -2.1905800e-03,
     + -2.8958588e-03, -2.5358774e-03, -1.2843649e-03, -5.3286898e-04,
     + 7.4579911e-04, -1.6767264e-04, -4.2887366e-03, 4.2278337e-04,
     + 7.6373293e-04, 4.2005518e-04, -6.8090395e-04, -6.0485019e-03,
     + -1.0334845e+02, -1.1658602e+01, -1.3596275e+00, -1.6868026e-01,
     + -1.5297521e-02, 1.9085558e-03, 3.9859865e-03, 6.5122065e-04,
     + 1.3378927e-04, -8.3119250e-05, 4.2651936e-04, 1.1508797e-03,
     + 1.4183458e-03, -7.2446345e-04, 2.7098220e-03, -9.3159652e-04,
     + -3.1174513e-04, 7.4184324e-04, -2.0898838e-03, -5.3805151e-04,
     + 6.7881586e+02, 2.3996426e+02, 1.0289382e+02, 4.3261219e+01,
     + 1.5489647e+01, 3.7393535e+00, -5.4214971e-01, -1.7326251e+00,
     + -1.6641474e+00, -1.2943173e+00, -9.2558708e-01, -6.4393295e-01,
     + -4.3997609e-01, -2.9998130e-01, -2.1569628e-01, -1.0052807e-01,
     + -5.4107832e-02, -3.1841533e-02, -6.9434653e-03, 2.3794158e-04,
     + 9.2321943e+02, 1.7073732e+02, 4.0384662e+01, 1.1652764e+01,
     + 2.9256554e+00, 7.4987252e-02, -6.5180169e-01, -7.0380688e-01,
     + -5.3706707e-01, -3.6434168e-01, -2.3141400e-01, -1.3873848e-01,
     + -8.9266476e-02, -5.5630285e-02, -3.0416500e-02, -1.3682222e-02,
     + -7.2526950e-03, -3.8539518e-03, 1.4495812e-03, -1.7114122e-02,
     + 5.7429831e+02, 9.0882864e+01, 1.6243871e+01, 3.2167128e+00,
     + 5.3464832e-01, -4.3270034e-02, -1.3216745e-01, -1.0332596e-01,
     + -6.7875801e-02, -4.1815597e-02, -2.4479644e-02, -1.6426570e-02,
     + -7.2854071e-03, -3.7535315e-03, -5.3761954e-03, -3.7320129e-05,
     + 1.3382847e-03, 1.0851986e-03, 2.0993845e-03, -8.1042585e-03,
     + 1.1998832e+02, 1.8443184e+01, 3.4555994e+00, 7.3374303e-01,
     + 1.5744947e-01, 5.1109576e-02, -9.3728478e-03, -4.0738095e-03,
     + -3.6118072e-03, -1.3609739e-03, -8.0523490e-04, -2.2146245e-03,
     + 3.0109837e-04, 2.1705092e-03, 2.2125603e-03, -1.9654785e-03,
     + -2.0243677e-04, 2.3746230e-03, 1.7808657e-03, 2.3257961e-03,
     + 2.5098425e+02, 7.0327468e+01, 2.5163368e+01, 9.5598735e+00,
     + 3.1766519e+00, 6.9997508e-01, -1.0105017e-01, -3.3563636e-01,
     + -3.0868084e-01, -2.2974221e-01, -1.5778551e-01, -1.0721136e-01,
     + -6.9827725e-02, -4.8295432e-02, -2.6544149e-02, -1.2844869e-02,
     + -2.2779111e-03, -5.5421259e-03, -2.5160836e-03, -4.3322721e-03,
     + 4.4982698e+02, 7.8288698e+01, 1.5901865e+01, 3.7476744e+00,
     + 7.5117610e-01, -3.3466185e-02, -1.6165248e-01, -1.6853514e-01,
     + -1.2078091e-01, -7.7984141e-02, -4.9084901e-02, -2.9144774e-02,
     + -1.6045315e-02, -9.0927656e-03, -9.3273266e-03, -2.0526870e-03,
     + 2.7335575e-03, -2.8674620e-03, -3.1245271e-03, 1.9694385e-03,
     + 3.3483093e+02, 4.8958231e+01, 8.0114909e+00, 1.4511474e+00,
     + 2.3734357e-01, -5.3871003e-03, -2.7159026e-02, -2.8416759e-02,
     + -1.7242181e-02, -1.1035733e-02, -6.2832439e-03, -5.9898984e-03,
     + -6.9100009e-04, 1.5797174e-03, 2.9452577e-03, -6.8021589e-04,
     + -1.0880688e-03, -5.5066781e-04, -2.0640439e-03, -2.5156350e-03,
     + 1.1447097e+02, 1.4685005e+01, 2.2951474e+00, 4.1279583e-01,
     + 7.7443035e-02, 1.4146005e-02, -5.5252954e-03, -1.0890799e-04,
     + -3.3356169e-04, 1.2811512e-04, -8.1953839e-04, 1.2254260e-03,
     + 1.4351044e-04, -9.4427609e-05, -2.1044931e-03, -3.6639672e-05,
     + 6.8399162e-05, -7.8639822e-04, -6.5368847e-05, -2.4770032e-03,
     + 7.8351971e+01, 1.9863070e+01, 6.4333653e+00, 2.2775385e+00,
     + 8.3617008e-01, 2.5404561e-01, -2.5265326e-03, -4.7671448e-02,
     + -4.9025088e-02, -3.7128781e-02, -2.6406123e-02, -1.1002844e-02,
     + -1.0107075e-02, -7.3986968e-03, -2.7521017e-03, -1.5442475e-03,
     + -5.2602394e-03, -2.3957439e-03, -3.7385979e-03, -2.3616899e-03,
     + 2.0566666e+02, 3.3090034e+01, 6.0687986e+00, 1.2383281e+00,
     + 2.4772390e-01, 1.6814833e-02, -4.1986414e-02, -3.4688794e-02,
     + -2.4036704e-02, -1.4873275e-02, -9.0452640e-03, -6.6238920e-03,
     + -1.5361628e-03, -4.8452147e-04, 1.5058890e-03, -4.6053696e-04,
     + -2.3611404e-03, -2.2434226e-03, -2.1717362e-03, -1.3402891e-03,
     + 1.8071414e+02, 2.4242215e+01, 3.6744204e+00, 6.1470970e-01,
     + 1.0613833e-01, 1.9628383e-02, -4.2797093e-03, -3.6806873e-03,
     + -4.2033833e-03, -2.0047935e-03, -1.9863224e-03, 3.2867311e-03,
     + -1.2540643e-03, -3.4730196e-04, 3.5182102e-04, 5.0230681e-04,
     + -3.9538253e-04, 8.1912914e-04, -7.6549491e-03, -5.6620913e-03,
     + 2.9223336e+01, 6.9512621e+00, 2.0112322e+00, 6.5308330e-01,
     + 1.8478337e-01, 3.0020730e-02, 1.5513745e-03, -5.2485026e-03,
     + -6.9436047e-03, -5.1577043e-03, -3.9084129e-03, -4.1274331e-03,
     + 3.0044649e-04, -2.0053608e-03, -4.0556923e-03, 6.7705908e-04,
     + 4.1200572e-03, -1.6319371e-03, -4.8710920e-03, -1.6303397e-03,
     + 8.5605839e+01, 1.2870179e+01, 2.2515125e+00, 4.1856552e-01,
     + 7.6171786e-02, -6.1497967e-03, -1.1388175e-02, -6.1622720e-03,
     + -4.0532369e-03, -2.3074166e-03, -2.2297735e-03, -1.3574466e-03,
     + -2.6423165e-04, 7.1657413e-04, 2.7009891e-03, 1.4978903e-03,
     + 4.9306655e-04, -8.4797936e-04, -3.7977928e-03, -5.5366721e-03,
     + 8.3971657e+01, 1.0322748e+01, 1.4914429e+00, 2.3600339e-01,
     + 3.5840371e-02, 9.5247101e-03, -1.0056573e-02, -2.2331711e-04,
     + 4.3103945e-04, -3.9274384e-04, -4.9651105e-04, -1.4282960e-03,
     + -1.1444543e-03, -1.2445884e-03, -1.2481127e-03, 1.2814045e-04,
     + 3.1062222e-04, 8.5971286e-04, -5.8290829e-03, -6.7735005e-03,
     + 1.5333941e+01, 2.7447892e+00, 6.5798213e-01, 1.6028167e-01,
     + 6.7012532e-02, -6.0381982e-03, 1.5298838e-02, -1.0149721e-03,
     + -4.8719226e-05, -1.4742346e-03, -5.2730522e-05, -2.5362559e-03,
     + 6.7240673e-04, 1.5294850e-03, 7.9558297e-03, -6.8939119e-04,
     + -1.3575514e-03, 7.0846540e-04, -3.7545723e-03, 5.0022990e-03,
     + 3.4141077e+01, 4.7365271e+00, 7.6808936e-01, 1.3044238e-01,
     + 3.7875081e-02, -1.7308256e-02, 3.2923217e-03, -3.2655709e-04,
     + -1.0572710e-04, 1.6880895e-04, 2.9101392e-04, 2.5813484e-03,
     + -6.0481861e-04, 3.9228504e-04, -3.6950550e-03, 1.8174963e-04,
     + 1.3470112e-03, -1.5117922e-03, -4.9383277e-03, -3.9287929e-03,
     + 5.6399463e+00, 9.3447871e-01, 1.9261742e-01, 6.8160793e-02,
     + -8.2025029e-06, 3.3224875e-02, 4.4947822e-03, 1.1008833e-03,
     + 1.4720834e-04, 5.5057327e-04, -2.5945420e-04, 1.4289963e-03,
     + 2.8954785e-04, 7.5615216e-04, -8.6700489e-03, 3.0866868e-03,
     + 3.1569093e-04, -1.2882760e-03, -1.3299036e-03, 8.0905234e-04,
     + 1.2852571e+01, 1.6329128e+00, 2.5173804e-01, 5.2879618e-02,
     + 6.9979716e-04, 2.8408487e-03, 1.7863049e-03, -4.0412308e-04,
     + 3.8426673e-04, 2.4677419e-04, -4.2832584e-04, -1.8070501e-03,
     + -3.7516888e-04, -1.6165791e-03, 3.7599178e-03, -1.6959082e-03,
     + -2.3468636e-03, -1.0232667e-03, -2.7698835e-03, -2.8582541e-04/
*
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         ind=1
         do ilam=1,34
           call dcopy(20,vl(ind),1,vec,1)
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
       do ilam=1,34
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
       call dcopy(33,v(2),1,vvl(1),1)
* convert to hartree
       econv=1.d0/219474.6d0
       call dscal(33,econv,vvl,1)
* isotropic term - divide by sqrt(4*pi)
       vv0 = v(1)*econv/s4pi
*
       return
       end
