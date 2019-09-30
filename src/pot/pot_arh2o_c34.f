*  system:  H2O + Ar PES
*  written by P. Dagdigian, May 2013
*
*  program to compute potential vs geometry from:
*  J. Makarewicz, J. Chem. Phys. 129, 184310 (2008)
*  CCSD(T) calculation, basis set included bond functions
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
* sqrt(4*pi)
      sq4pi = 3.544907701811032d0
      potnam='Ar-H2O Makarewicz PES-34 coeffs'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      call pot(vv0,r)
*  vlm coefficient is returned in atomic units (hartree)
*  convert from atomic units for printout
      econv=219474.6d0
      write (6, 100) vv0*econv*sq4pi, (econv*vvl(i), i=1,33)
100   format(' v(lam,0):',6(1pe16.8)/' v(lam,1):',5(1pe16.8)/
     :    ' v(lam,2):',5(1pe16.8)/' v(lam,3):',4(1pe16.8)/
     :    ' v(lam,4):',4(1pe16.8)/' v(lam,5):',3(1pe16.8)/
     :    ' v(lam,6):',3(1pe16.8)/' v(lam,7):',2(1pe16.8)/
     :    ' v(lam,8):',2(1pe16.8))
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
      potnam='Ar-H2O Makarewicz PES-34 coeffs'
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
      dimension v(33)
      dimension csplin(20,34)
      dimension rr(23), vl(680),vec(34)
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(33)
* hyperbolic tangent scaling factor
      data alph /1.2d0/
      data rmax /13d0/
* sqrt(4*pi)
      data sq4pi / 3.544907701811032d0 /
* here are the 20 values of R
      data rr/3d0,3.5d0,4d0,4.5d0,5d0,5.5d0,6d0,6.5d0,7d0,7.5d0,
     +  8d0,8.5d0,9d0,9.5d0,10d0,11d0,12d0,13d0,15d0,20d0/
      data ifirst /0/
* here are column ordered legendre expansion coefficients (12 total) for the
* potential at each of 34 values of R (680 values)
      data vl/
     + 5.7076075e+05, 1.8184581e+05, 6.3246861e+04, 2.2265894e+04,
     + 7.2914885e+03, 1.9022779e+03, 1.2117649e+02, -3.5169989e+02,
     + -3.9333508e+02, -3.2067501e+02, -2.3703018e+02, -1.6931365e+02,
     + -1.1997831e+02, -8.5311667e+01, -6.1162314e+01, -3.2540248e+01,
     + -1.8689702e+01, -1.8689702e+01, -4.4462018e+00, -6.7227567e-01,
* 
     + 6.7928480e+04, 2.6194125e+04, 8.8674681e+03, 2.6838803e+03,
     + 6.8763458e+02, 1.0910328e+02, -2.8170910e+01, -4.2560593e+01,
     + -3.0425919e+01, -1.7633593e+01, -9.1227757e+00, -4.3910286e+00,
     + -2.0522819e+00, -9.9442314e-01, -5.3211102e-01, -1.5576289e-01,
     + 9.3094454e-02, 9.3094454e-02, 4.5021599e-02, 9.3033912e-05,
* 
     + -5.2213030e+03, -1.8519426e+03, -7.7319423e+02, -3.4135413e+02,
     + -1.4825668e+02, -6.1534032e+01, -2.4251529e+01, -9.1479132e+00,
     + -3.3941018e+00, -1.3122722e+00, -5.7067559e-01, -2.8662163e-01,
     + -1.5491733e-01, -7.9683046e-02, -3.2605394e-02, 1.1269886e-02,
     + 1.7748329e-02, 1.7748329e-02, 1.1139089e-05, 4.3445586e-09,
* 
     + -3.6503141e+03, -1.0026389e+03, -3.2112854e+02, -1.0975177e+02,
     + -3.6743825e+01, -1.1182962e+01, -2.7441095e+00, -3.1060373e-01,
     + 1.9827388e-01, 1.9378137e-01, 1.0905897e-01, 4.7322438e-02,
     + 1.6198075e-02, 4.1905746e-03, 1.1698663e-03, 2.0448348e-03,
     + 2.4971093e-03, 2.4971093e-03, 2.9461152e-05, -4.1468983e-08,
* 
     + -3.1526073e+02, -6.9488217e+01, -1.9618065e+01, -6.3955892e+00,
     + -2.1763721e+00, -7.3030352e-01, -2.3594713e-01, -7.2587466e-02,
     + -2.1338494e-02, -6.3336939e-03, -2.3164328e-03, -1.1689488e-03,
     + -7.3748668e-04, -4.8443675e-04, -3.3033422e-04, -1.0374221e-04,
     + -1.6887692e-05, -1.6887692e-05, 7.0538920e-07, -1.0530885e-07,
* 
     + 1.8286185e+02, 3.7348333e+01, 8.7233830e+00, 2.0326037e+00,
     + 4.2077615e-01, 5.7820080e-02, -8.7463124e-03, -1.2706581e-02,
     + -7.3059150e-03, -3.4545804e-03, -1.4535895e-03, -5.4734158e-04,
     + -1.6560988e-04, -2.4864303e-05, 4.4605146e-06, 8.1791055e-06,
     + 3.4718844e-06, 3.4718844e-06, 6.3994595e-07, 1.2229491e-08,
* 
     + -7.6860414e+04, -3.2010615e+03, 4.3958518e+03, 2.7356298e+03,
     + 1.1749152e+03, 4.2191668e+02, 1.2328749e+02, 1.9510179e+01,
     + -1.0296369e+01, -1.4910932e+01, -1.2461040e+01, -9.0433337e+00,
     + -6.3017829e+00, -4.3908851e+00, -3.0838182e+00, -1.4148978e+00,
     + -3.0744572e-01, -3.0744572e-01, 2.2483817e-01, 1.0215044e-02,
* 
     + 7.3013693e+04, 2.4546053e+04, 8.3364201e+03, 2.8767386e+03,
     + 9.8918137e+02, 3.2310576e+02, 8.9957150e+01, 1.3271939e+01,
     + -7.8438773e+00, -1.0705514e+01, -8.6954959e+00, -6.1636142e+00,
     + -4.2070845e+00, -2.8873379e+00, -2.0265903e+00, -1.0550031e+00,
     + -5.4098433e-01, -5.4098433e-01, -4.8180307e-03, -8.7068302e-07,
* 
     + 5.0827712e+03, 1.3517706e+03, 3.8447284e+02, 1.1613131e+02,
     + 3.5082040e+01, 9.2810344e+00, 1.2831480e+00, -7.9070690e-01,
     + -1.0073180e+00, -7.6306640e-01, -4.9410299e-01, -3.0109143e-01,
     + -1.8169171e-01, -1.1243573e-01, -7.2802631e-02, -3.4923959e-02,
     + -1.8198214e-02, -1.8198214e-02, -2.4749389e-04, -6.4341310e-10,
* 
     + -1.0419903e+03, -2.9517434e+02, -9.5464930e+01, -3.1173735e+01,
     + -9.5050697e+00, -2.4901202e+00, -4.4225411e-01, 3.9328952e-02,
     + 9.1051037e-02, 5.7493327e-02, 2.5591341e-02, 7.9247835e-03,
     + 4.4644048e-04, -1.7697965e-03, -1.8640609e-03, -7.6607143e-04,
     + -1.2636112e-04, -1.2636112e-04, -1.2825257e-06, 2.7103725e-08,
* 
     + -4.1042716e+02, -9.2270183e+01, -2.4713597e+01, -6.9210944e+00,
     + -1.8687045e+00, -4.4999843e-01, -8.0544460e-02, 1.8659048e-04,
     + 9.8862018e-03, 6.5106267e-03, 2.9365672e-03, 9.7807447e-04,
     + 1.5798423e-04, -1.0374938e-04, -1.5594860e-04, -7.6496800e-05,
     + -2.5171905e-05, -2.5171905e-05, -7.2208174e-07, 2.3993036e-08,
* 
     + 6.0444332e+04, 1.8091651e+04, 5.2368788e+03, 1.4578763e+03,
     + 3.7253704e+02, 7.0316711e+01, -6.7612122e+00, -2.0411769e+01,
     + -1.7734301e+01, -1.2352905e+01, -7.8592606e+00, -4.8138489e+00,
     + -2.9463618e+00, -1.8611885e+00, -1.2406554e+00, -6.4489755e-01,
     + -3.5723147e-01, -3.5723147e-01, -1.3883035e-02, -1.1358896e-04,
* 
     + 2.1405631e+04, 6.6855283e+03, 2.2446266e+03, 7.9306539e+02,
     + 2.7976013e+02, 9.1701202e+01, 2.4481992e+01, 2.7222445e+00,
     + -2.7583250e+00, -3.1070658e+00, -2.2813784e+00, -1.4785304e+00,
     + -9.3369470e-01, -6.0695581e-01, -4.1671854e-01, -2.2509326e-01,
     + -1.2553460e-01, -1.2553460e-01, -1.0369587e-03, 4.2856071e-08,
* 
     + 3.7258534e+03, 9.9741072e+02, 3.1093703e+02, 1.0464481e+02,
     + 3.4932436e+01, 1.0712532e+01, 2.6863322e+00, 3.3763000e-01,
     + -1.7388850e-01, -1.8375012e-01, -1.0823906e-01, -5.0180421e-02,
     + -1.9692387e-02, -7.1097665e-03, -3.2446319e-03, -2.8822198e-03,
     + -2.7932902e-03, -2.7932902e-03, -4.2435487e-05, 3.3799804e-08,
* 
     + 4.4971888e+01, 1.7690905e+00, -2.7853498e-01, 3.2693646e-01,
     + 4.0119870e-01, 2.6787527e-01, 1.4402903e-01, 6.9232904e-02,
     + 3.0890123e-02, 1.3128125e-02, 5.3532515e-03, 2.0782978e-03,
     + 7.7482218e-04, 2.5588111e-04, 7.5094245e-05, -9.3456115e-06,
     + -1.1609040e-05, -1.1609040e-05, -5.2162776e-07, -3.4720711e-08,
* 
     + -2.1104256e+02, -4.3131743e+01, -1.0222033e+01, -2.4525934e+00,
     + -5.3697999e-01, -8.7633847e-02, 2.1290256e-03, 1.1632191e-02,
     + 7.5277470e-03, 3.6966583e-03, 1.5552540e-03, 5.7996292e-04,
     + 1.7475712e-04, 2.5498276e-05, -1.4037104e-05, -1.0498240e-05,
     + -4.6903682e-06, -4.6903682e-06, 4.7518777e-07, -2.3906848e-08,
* 
     + 9.4629582e+03, 3.2986339e+03, 1.1827878e+03, 4.3277538e+02,
     + 1.5500063e+02, 5.0253383e+01, 1.1994879e+01, -4.9214874e-01,
     + -3.4354228e+00, -3.2677985e+00, -2.3894043e+00, -1.5760103e+00,
     + -9.9613584e-01, -6.2501994e-01, -3.9876688e-01, -1.7857862e-01,
     + -8.7936753e-02, -8.7936753e-02, -3.4838333e-04, 1.0523165e-06,
* 
     + 6.2376490e+03, 1.8364211e+03, 6.1355645e+02, 2.1771009e+02,
     + 7.5758878e+01, 2.3883702e+01, 5.9435340e+00, 5.3039329e-01,
     + -6.4002074e-01, -6.1009739e-01, -3.7613415e-01, -1.9577014e-01,
     + -9.5469327e-02, -4.8939582e-02, -3.0098723e-02, -1.9566474e-02,
     + -1.3985778e-02, -1.3985778e-02, -1.1549297e-04, 8.8290524e-09,
* 
     + 1.9797591e+03, 5.2198757e+02, 1.6180762e+02, 5.2772794e+01,
     + 1.6673522e+01, 4.7661335e+00, 1.0944228e+00, 1.1308539e-01,
     + -7.1500649e-02, -6.2918039e-02, -3.0549149e-02, -9.4176387e-03,
     + -3.2030553e-05, 2.7286469e-03, 2.6900026e-03, 9.9767826e-04,
     + 6.7185081e-05, 6.7185081e-05, -3.1897003e-06, -1.1001547e-07,
* 
     + 2.6953636e+02, 5.8986208e+01, 1.5873134e+01, 4.5971530e+00,
     + 1.3238227e+00, 3.5755640e-01, 8.4308058e-02, 1.4194608e-02,
     + -3.7291525e-04, -1.7910404e-03, -9.1563649e-04, -2.6189732e-04,
     + 4.6245564e-05, 1.2278199e-04, 1.1333607e-04, 5.1709919e-05,
     + 1.4190950e-05, 1.4190950e-05, 4.6952268e-07, -7.5624029e-09,
* 
     + 1.9171861e+03, 5.6708891e+02, 1.8852816e+02, 6.8004097e+01,
     + 2.4521481e+01, 7.9629186e+00, 1.8642019e+00, -1.0160344e-01,
     + -5.3622344e-01, -4.8744591e-01, -3.4107436e-01, -2.1498221e-01,
     + -1.2998700e-01, -7.8403129e-02, -4.8572372e-02, -2.1293969e-02,
     + -1.0545539e-02, -1.0545539e-02, -5.0155495e-05, 2.8032899e-08,
* 
     + 2.0020197e+03, 5.7726590e+02, 1.9248403e+02, 6.6969759e+01,
     + 2.2397357e+01, 6.7085195e+00, 1.5815878e+00, 1.4241358e-01,
     + -1.3640347e-01, -1.1800444e-01, -6.0601422e-02, -2.1953894e-02,
     + -3.7675199e-03, 2.2651559e-03, 2.9252237e-03, 6.4230055e-04,
     + -5.9302830e-04, -5.9302830e-04, -9.7376259e-06, 5.1880742e-09,
* 
     + 9.3066897e+02, 2.3265721e+02, 6.7827586e+01, 2.0388992e+01,
     + 5.8266566e+00, 1.4557067e+00, 2.5476923e-01, -1.5291156e-02,
     + -4.4484217e-02, -2.8255422e-02, -1.2731290e-02, -4.1732538e-03,
     + -5.2796557e-04, 6.2441713e-04, 7.5621816e-04, 3.7053620e-04,
     + 9.9677915e-05, 9.9677915e-05, -1.0621518e-06, -1.8400341e-08,
* 
     + 2.2836095e+02, 4.6084345e+01, 1.1089825e+01, 2.7826573e+00,
     + 6.6698063e-01, 1.3633801e-01, 1.5693469e-02, -5.1951222e-03,
     + -5.2988496e-03, -2.9243120e-03, -1.2704655e-03, -4.6813704e-04,
     + -1.3220315e-04, -1.2666988e-05, 2.1408753e-05, 1.4441294e-05,
     + 4.9711081e-06, 4.9711081e-06, -6.3779798e-07, 6.4216992e-08,
* 
     + 4.3930499e+02, 1.3775811e+02, 5.1227626e+01, 2.0587796e+01,
     + 8.2612870e+00, 3.1708342e+00, 1.1328549e+00, 3.6482447e-01,
     + 9.8796689e-02, 1.7213481e-02, -2.7567630e-03, -5.1432134e-03,
     + -3.8961419e-03, -2.5936421e-03, -1.8238496e-03, -1.1554061e-03,
     + -7.5735946e-04, -7.5735946e-04, -4.2686019e-06, 9.0245665e-08,
* 
     + 6.9135912e+02, 1.9055743e+02, 6.0629833e+01, 1.9778575e+01,
     + 6.1154007e+00, 1.6642843e+00, 3.3489720e-01, 3.4130034e-03,
     + -4.4928174e-02, -3.1343550e-02, -1.4480224e-02, -4.5499358e-03,
     + -2.2801550e-04, 1.1024646e-03, 1.1591699e-03, 5.1039955e-04,
     + 1.0414664e-04, 1.0414664e-04, -9.5136059e-07, 1.8693077e-07,
* 
     + 4.2148602e+02, 9.5959879e+01, 2.5314215e+01, 6.8061412e+00,
     + 1.7028910e+00, 3.4845073e-01, 2.9321602e-02, -2.3751642e-02,
     + -2.0280167e-02, -1.0825241e-02, -4.7367692e-03, -1.7196109e-03,
     + -4.5704524e-04, -1.1861351e-05, 1.0494277e-04, 7.4557482e-05,
     + 2.4852652e-05, 2.4852652e-05, 5.9997600e-07, 7.2245530e-08,
* 
     + 1.2110632e+02, 3.9299739e+01, 1.4989301e+01, 6.0109276e+00,
     + 2.3929042e+00, 9.2951837e-01, 3.5165680e-01, 1.3097404e-01,
     + 4.9081216e-02, 1.9115443e-02, 8.0964077e-03, 3.7994168e-03,
     + 1.9487450e-03, 1.0308599e-03, 5.5065473e-04, 1.1189160e-04,
     + -5.6892230e-06, -5.6892230e-06, 1.5738779e-07, -4.6041260e-08,
* 
     + 2.3649212e+02, 6.0721567e+01, 1.7873495e+01, 5.3252488e+00,
     + 1.4827632e+00, 3.5026551e-01, 4.9686478e-02, -1.2358007e-02,
     + -1.5564258e-02, -9.1141241e-03, -4.0675445e-03, -1.4145291e-03,
     + -2.8171928e-04, 9.9101447e-05, 1.8315339e-04, 1.0053588e-04,
     + 3.0786014e-05, 3.0786014e-05, 5.0476829e-08, -8.5978576e-10,
* 
     + 1.7217534e+02, 3.5223266e+01, 8.3186462e+00, 1.9820878e+00,
     + 4.2740047e-01, 6.5562756e-02, -4.5418516e-03, -1.0768851e-02,
     + -6.8434273e-03, -3.2671997e-03, -1.3918049e-03, -5.3276990e-04,
     + -1.7991552e-04, -3.6469150e-05, -2.7164181e-06, 1.3957303e-05,
     + 4.3260283e-06, 4.3260283e-06, 7.8501866e-07, -8.9565472e-08,
* 
     + 4.5304975e+01, 1.3108425e+01, 4.5112101e+00, 1.6464275e+00,
     + 6.0477651e-01, 2.2066161e-01, 8.0280369e-02, 2.9585676e-02,
     + 1.1422292e-02, 4.6264259e-03, 2.0844159e-03, 1.0203282e-03,
     + 5.4036953e-04, 2.9492888e-04, 1.6082077e-04, 4.1584374e-05,
     + 8.2653066e-06, 8.2653066e-06, -9.8290309e-07, -9.4344843e-08,
* 
     + 8.0550215e+01, 1.8765853e+01, 4.9898262e+00, 1.3331397e+00,
     + 3.2703112e-01, 6.3159727e-02, 2.9466083e-03, -5.9559392e-03,
     + -4.6539055e-03, -2.4210390e-03, -1.0401412e-03, -3.9620250e-04,
     + -1.1642752e-04, -1.6932033e-05, 1.0132047e-05, 1.2505181e-05,
     + 5.4628301e-06, 5.4628301e-06, -1.0880958e-06, 1.2576045e-08,
* 
     + 1.4037954e+01, 3.6976866e+00, 1.1539829e+00, 3.8835403e-01,
     + 1.3043483e-01, 4.4414443e-02, 1.5204512e-02, 5.3874746e-03,
     + 2.0391542e-03, 8.0182719e-04, 3.6075960e-04, 1.7057586e-04,
     + 1.0104257e-04, 4.9297807e-05, 2.5282847e-05, 6.0020874e-06,
     + 3.0952690e-06, 3.0952690e-06, 1.6565398e-06, 2.5392664e-08,
* 
     + 2.6685191e+01, 5.5563330e+00, 1.3279359e+00, 3.1351889e-01,
     + 6.6891294e-02, 9.5503049e-03, -1.1290009e-03, -1.8851960e-03,
     + -1.0289040e-03, -5.7212976e-04, -2.4092348e-04, -8.4547552e-05,
     + -2.5463518e-05, -6.3634533e-06, 4.1170102e-07, -7.0056667e-07,
     + 1.4277967e-06, 1.4277967e-06, -7.1600105e-07, 5.6039543e-08  /
*
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         ind=1
         do ilam=1,34
           call dcopy(34,vl(ind),1,vec,1)
*    evaluate derivative at first point
           der1=(vec(2)-vec(1))/(rr(2)-rr(1))
           call dspline(rr,vec,20,der1,0d0,csplin(1,ilam))
           ind = ind + 20
         enddo
         ifirst = 1
      end if
* r^-6 fit to isotropic part of potential
      c6sum = -6.376139065322550e+07
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
      vv0 = v(1)*econv/sq4pi
*
      return
      end
