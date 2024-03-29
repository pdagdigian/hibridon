cstart unix-ifort
cdec$ fixedformlinesize:132
cend

*  System:  O(1D,3P) + Kr
*  Theory  Level: MRCISD+Q(Davidson) aug-cc-pvqz-DK BSSE and SC corrected
*   Scalar relativistic effects included by using all-electron
*   Douglas-Kroll Hamiltonian and integrals
*   calculation of potential energy curves and spin-orbit matrix
*   elements by j.klos 
*
*   written by p. dagdigian
*   adapted to O-Kr by j. klos
*   current revision date:  6-oct-2015 by pjd
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
      potnam='O(1D,3P)-Kr MRCISD+Q'
      ibasty=22
      lammin(1)=1
      lammax(1)=12
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
      common /covvl/ vvl(19)
      include "common/parpot"
      potnam='O(1D,3P)-Kr MRCISD+Q'
      econv=219474.6d0
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      if (r.le.0.d0) goto 99
      call pot(vv0,r)
      write (6, 100) (econv*vvl(i), i=1,19)
100   format('  ',6(1pe16.8)/'  ',6(1pe16.8)/
     +  '  ',7(1pe16.8))
      goto 1
99    rr=2.6d0
      dr=0.25d0
      open (unit=12,file='okr_vlms.txt')
      write(12,109)
109   format(' R/bohr 1Sig  1Pi  1Del  3Sig  3Pi  Axy   Azy',
     :  '  Bss  Byx  Bxs  Bsy  Bxd  Axy(R=inf)   Azy(R=inf)',
     :  '  Bss(R=inf)  Byx(R=inf)  Bxs(R=inf)  Bsy(R=inf)',
     :  '  Bxd(R=inf)')
      do i=1,201
        call pot(vv0,rr)
        write (12,110) rr,(econv*vvl(j),j=1,19)
110     format(f7.2,19(1pe16.8))
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
*    vvl(1,3)   contains the 1Sigma+, 1Pi, 1Delta energies
*    vvl(4,5)   contains the 3Sigma- and 3Pi PE energies
*    vvl(6,7)   contains Axy and Axz, the spin-orbit matrix elements
*               within the 3P state
*    vvl(8,12)  contains Bss, Byx, Bxs, Bsy, and Bxd, the spin-orbit
*               matrix elements coupling the 1D and 3P states
*    vvl(13,19) values of the Axy, Azy, Bss, Byx, Bxs, Bsy, and Bxd
*               matrix elements at R=inf
**  variable in common block /conlam/ used here
*    nlam:    the number of angular coupling terms actually used
*  variable in common block /covvl/
*    vvl:     array to store r-dependence of each
*             angular term in the potential
*  variable in common block /coconv/
*    econv:     conversion factor from cm-1 to hartrees
*  variables in common /cosysr/
*    isrcod:    number of real parameters
*    en1d:      asymptotic energy of the 1D state (cm-1)
* 
* author:  paul dagdigian
* latest revision date:  6-oct-2015
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/parbas"
      common /covvl/ vvl(19)

      dimension rr(85),vl(425),rso(22),somat(154)
      dimension csplel(85,5),csplso(22,12),vecel(85),vecso(22)
      dimension v(19)
*
*  values or R for PE curves
      data rr /
     +  2.600,  2.700,  2.800,  2.900,  3.000,
     +  3.100,  3.200,  3.300,  3.400,  3.500,
     +  3.600,  3.700,  3.800,  3.900,  4.000,
     +  4.100,  4.200,  4.300,  4.400,  4.500,
     +  4.600,  4.700,  4.800,  5.000,  5.100,
     +  5.200,  5.300,  5.400,  5.500,  5.600,
     +  5.700,  5.800,  5.900,  6.000,  6.100,
     +  6.200,  6.300,  6.400,  6.500,  6.600,
     +  6.700,  6.800,  6.900,  7.000,  7.100,
     +  7.200,  7.300,  7.400,  7.500,  7.600,
     +  7.700,  7.800,  7.900,  8.000,  8.250,
     +  8.500,  8.750,  9.000,  9.250,  9.500,
     +  9.750, 10.000, 10.250, 10.500, 10.750,
     + 11.000, 11.250, 11.500, 11.750, 12.000,
     + 12.500, 13.000, 13.500, 14.000, 15.000,
     + 16.000, 17.000, 18.000, 19.000, 20.000,
     + 22.000, 24.000, 26.000, 28.000, 30.000 /
*  values of the 5 PE curves
      data vl /  
C       1Sigma
     +  4.5922918e+04, 3.5208888e+04, 2.7442494e+04, 2.1895398e+04, 1.8012062e+04,
     +  1.5368450e+04, 1.3648335e+04, 1.2618736e+04, 1.2105361e+04, 1.1973199e+04,
     +  1.2113401e+04, 1.2435477e+04, 1.2863541e+04, 1.3335279e+04, 1.3802306e+04,
     +  1.4230581e+04, 1.4599899e+04, 1.4902062e+04, 1.5138026e+04, 1.5314657e+04,
     +  1.5441779e+04, 1.5529930e+04, 1.5588928e+04, 1.5651227e+04, 1.5666248e+04,
     +  1.5675852e+04, 1.5682563e+04, 1.5688040e+04, 1.5693310e+04, 1.5698949e+04,
     +  1.5705230e+04, 1.5712226e+04, 1.5719889e+04, 1.5728101e+04, 1.5736713e+04,
     +  1.5745568e+04, 1.5754517e+04, 1.5763428e+04, 1.5772188e+04, 1.5780709e+04,
     +  1.5788924e+04, 1.5796786e+04, 1.5804264e+04, 1.5811343e+04, 1.5818019e+04,
     +  1.5824294e+04, 1.5830178e+04, 1.5835685e+04, 1.5840829e+04, 1.5845627e+04,
     +  1.5850096e+04, 1.5854252e+04, 1.5858113e+04, 1.5861693e+04, 1.5869520e+04,
     +  1.5875932e+04, 1.5881153e+04, 1.5885389e+04, 1.5888824e+04, 1.5891613e+04,
     +  1.5893888e+04, 1.5895753e+04, 1.5897290e+04, 1.5898784e+04, 1.5899884e+04,
     +  1.5900806e+04, 1.5901583e+04, 1.5902240e+04, 1.5902798e+04, 1.5903274e+04,
     +  1.5904030e+04, 1.5904591e+04, 1.5905013e+04, 1.5905333e+04, 1.5905772e+04,
     +  1.5906063e+04, 1.5906248e+04, 1.5906368e+04, 1.5906449e+04, 1.5906504e+04,
     +  1.5906570e+04, 1.5906604e+04, 1.5906623e+04, 1.5906635e+04, 1.5906642e+04,
C      1Pi  
     +  8.8745234e+04, 7.6174438e+04, 6.6053732e+04, 5.7731464e+04, 5.0783665e+04,
     +  4.4940411e+04, 4.0022470e+04, 3.5897212e+04, 3.2454172e+04, 2.9594574e+04,
     +  2.7228403e+04, 2.5274585e+04, 2.3661842e+04, 2.2329124e+04, 2.1225381e+04,
     +  2.0308722e+04, 1.9545239e+04, 1.8907695e+04, 1.8374261e+04, 1.7927393e+04,
     +  1.7552912e+04, 1.7239252e+04, 1.6976890e+04, 1.6575676e+04, 1.6424570e+04,
     +  1.6299824e+04, 1.6197368e+04, 1.6113719e+04, 1.6045896e+04, 1.5991348e+04,
     +  1.5947891e+04, 1.5913666e+04, 1.5887086e+04, 1.5866808e+04, 1.5851694e+04,
     +  1.5840785e+04, 1.5833276e+04, 1.5828494e+04, 1.5825878e+04, 1.5824964e+04,
     +  1.5825368e+04, 1.5826777e+04, 1.5828934e+04, 1.5831632e+04, 1.5834702e+04,
     +  1.5838013e+04, 1.5841458e+04, 1.5844955e+04, 1.5848440e+04, 1.5851867e+04,
     +  1.5855198e+04, 1.5858409e+04, 1.5861483e+04, 1.5864407e+04, 1.5871039e+04,
     +  1.5876709e+04, 1.5881491e+04, 1.5885490e+04, 1.5888820e+04, 1.5891586e+04,
     +  1.5893884e+04, 1.5895796e+04, 1.5897388e+04, 1.5898719e+04, 1.5899833e+04,
     +  1.5900770e+04, 1.5901561e+04, 1.5902229e+04, 1.5902797e+04, 1.5903282e+04,
     +  1.5904051e+04, 1.5904622e+04, 1.5905051e+04, 1.5905377e+04, 1.5905823e+04,
     +  1.5906096e+04, 1.5906270e+04, 1.5906383e+04, 1.5906459e+04, 1.5906511e+04,
     +  1.5906573e+04, 1.5906606e+04, 1.5906624e+04, 1.5906636e+04, 1.5906642e+04,
C      1Delta 
     +  1.3875953e+05, 1.2453599e+05, 1.1237558e+05, 1.0168596e+05, 9.2081971e+04,
     +  8.3316161e+04, 7.5252733e+04, 6.7840608e+04, 6.1074345e+04, 5.4959374e+04,
     +  4.9491932e+04, 4.4652402e+04, 4.0406490e+04, 3.6709491e+04, 3.3510939e+04,
     +  3.0758471e+04, 2.8400626e+04, 2.6388669e+04, 2.4677638e+04, 2.3226846e+04,
     +  2.1999999e+04, 2.0965077e+04, 2.0094074e+04, 1.8749777e+04, 1.8237348e+04,
     +  1.7809842e+04, 1.7453989e+04, 1.7158471e+04, 1.6913660e+04, 1.6711385e+04,
     +  1.6544722e+04, 1.6407816e+04, 1.6295725e+04, 1.6204288e+04, 1.6130002e+04,
     +  1.6069929e+04, 1.6021604e+04, 1.5982969e+04, 1.5952301e+04, 1.5928168e+04,
     +  1.5909376e+04, 1.5894936e+04, 1.5884025e+04, 1.5875966e+04, 1.5870199e+04,
     +  1.5866263e+04, 1.5863780e+04, 1.5862440e+04, 1.5861990e+04, 1.5862226e+04,
     +  1.5862980e+04, 1.5864118e+04, 1.5865533e+04, 1.5867138e+04, 1.5871575e+04,
     +  1.5876103e+04, 1.5880358e+04, 1.5884177e+04, 1.5887513e+04, 1.5890376e+04,
     +  1.5892809e+04, 1.5894862e+04, 1.5896590e+04, 1.5898041e+04, 1.5899260e+04,
     +  1.5900286e+04, 1.5901151e+04, 1.5901881e+04, 1.5902501e+04, 1.5903028e+04,
     +  1.5903864e+04, 1.5904482e+04, 1.5904945e+04, 1.5905295e+04, 1.5905772e+04,
     +  1.5906063e+04, 1.5906248e+04, 1.5906368e+04, 1.5906449e+04, 1.5906504e+04,
     +  1.5906570e+04, 1.5906604e+04, 1.5906623e+04, 1.5906635e+04, 1.5906642e+04,
C     3Sigma 
     +  1.3127529e+05, 1.1674323e+05, 1.0410102e+05, 9.2735367e+04, 8.2285669e+04,
     +  7.2585856e+04, 6.3602670e+04, 5.5363046e+04, 4.7898465e+04, 4.1218418e+04,
     +  3.5305176e+04, 3.0118649e+04, 2.5604016e+04, 2.1698765e+04, 1.8338052e+04,
     +  1.5458361e+04, 1.2999770e+04, 1.0907198e+04, 9.1309564e+03, 7.6268351e+03,
     +  6.3559148e+03, 5.2842094e+03, 4.3822295e+03, 2.9891722e+03, 2.4574329e+03,
     +  2.0132545e+03, 1.6429530e+03, 1.3348805e+03, 1.0791403e+03, 8.6733756e+02,
     +  6.9236303e+02, 5.4820550e+02, 4.2978988e+02, 3.3283783e+02, 2.5374820e+02,
     +  1.8949437e+02, 1.3753656e+02, 9.5746911e+01, 6.2345485e+01, 3.5845849e+01,
     +  1.5008775e+01,-1.1966870e+00,-1.3628151e+01,-2.2995620e+01,-2.9886962e+01,
     + -3.4787560e+01,-3.8097268e+01,-4.0144692e+01,-4.1198944e+01,-4.1479795e+01,
     + -4.1166161e+01,-4.0403029e+01,-3.9307356e+01,-3.7973006e+01,-3.4046862e+01,
     + -2.9854176e+01,-2.5817591e+01,-2.2137562e+01,-1.8889235e+01,-1.6078924e+01,
     + -1.3678373e+01,-1.1644291e+01,-9.9291460e+00,-8.4867370e+00,-7.2748660e+00,
     + -6.2563310e+00,-5.3990640e+00,-4.6758450e+00,-4.0638580e+00,-3.5441580e+00,
     + -2.7221000e+00,-2.1159870e+00,-1.6630140e+00,-1.3203670e+00,-8.5536900e-01,
     + -5.7087600e-01,-3.8957900e-01,-2.7066200e-01,-1.9083200e-01,-1.3623700e-01,
     + -7.1199000e-02,-3.7480000e-02,-1.8955000e-02,-6.4160000e-03,-3.3233000e-04,
C       3Pi
     +  8.4108649e+04, 7.1092677e+04, 6.0428124e+04, 5.1480502e+04, 4.3851738e+04,
     +  3.7297655e+04, 3.1658417e+04, 2.6814200e+04, 2.2663693e+04, 1.9116318e+04,
     +  1.6090482e+04, 1.3513588e+04, 1.1322027e+04, 9.4607276e+03, 7.8823208e+03,
     +  6.5461573e+03, 5.4173351e+03, 4.4658362e+03, 3.6657989e+03, 2.9949203e+03,
     +  2.4339647e+03, 1.9663543e+03, 1.5778224e+03, 9.9072304e+02, 7.7266479e+02,
     +  5.9427463e+02, 4.4902990e+02, 3.3139467e+02, 2.3668336e+02, 1.6094115e+02,
     +  1.0083952e+02, 5.3585080e+01, 1.6840541e+01,-1.1344086e+01,-3.2590482e+01,
     + -4.8243151e+01,-5.9412621e+01,-6.7014741e+01,-7.1799737e+01,-7.4381214e+01,
     + -7.5258662e+01,-7.4836848e+01,-7.3442186e+01,-7.1336446e+01,-6.8728258e+01,
     + -6.5782695e+01,-6.2629303e+01,-5.9368717e+01,-5.6078082e+01,-5.2815680e+01,
     + -4.9624725e+01,-4.6536312e+01,-4.3571979e+01,-4.0745745e+01,-3.4327338e+01,
     + -2.8836849e+01,-2.4210029e+01,-2.0345511e+01,-1.7133211e+01,-1.4468424e+01,
     + -1.2257936e+01,-1.0421848e+01,-8.8933190e+00,-7.6172120e+00,-6.5483920e+00,
     + -5.6500850e+00,-4.8923820e+00,-4.2509290e+00,-3.7058830e+00,-3.2410550e+00,
     + -2.5014860e+00,-1.9522980e+00,-1.5391740e+00,-1.2246980e+00,-7.9456500e-01,
     + -5.3013300e-01,-3.6201200e-01,-2.5199500e-01,-1.7818900e-01,-1.2758100e-01,
     + -6.7007000e-02,-3.5425000e-02,-1.8002000e-02,-6.0670000e-03,-2.6599000e-04/
*  values or R for spin-orbit matrix elements
      data rso /  
     +  2.500,  2.600,  2.700,  2.800,  2.900,
     +  3.000,  3.250,  3.500,  3.750,  4.000,
     +  4.250,  4.500,  4.750,  5.000,  5.250,
     +  5.500,  6.000,  6.500,  7.000,  7.500,
     +  8.000,  8.500 /
*  values of Axy and Azy - spin-orbit matrix elements within 3P state
*  then values of the B spin-orbit matrix elements, coupling the 1D and 3P states
       data somat /  
C   Axy
     +   7.4310653e+02, 6.7069772e+02, 5.8499264e+02, 4.9384413e+02,
     +   4.0835129e+02, 3.3649932e+02, 2.1936549e+02, 1.5899654e+02,
     +   1.2577625e+02, 1.0638951e+02, 9.4737589e+01, 8.7655842e+01,
     +   8.3344890e+01, 8.0720606e+01, 7.9126056e+01, 7.8158596e+01,
     +   7.7217622e+01, 7.6873876e+01, 7.6749232e+01, 7.6704719e+01,
     +   7.6689289e+01, 7.6684179e+01,
C Axz 
     +  -2.9501700e+02,-3.0661668e+02,-3.0961370e+02,-3.0126331e+02,
     +  -2.8164291e+02,-2.5383538e+02,-1.7321228e+02,-1.0275692e+02,
     +  -5.0794812e+01,-1.4949288e+01, 9.0159980e+00, 2.4794027e+01,
     +   3.5100136e+01, 4.1810931e+01, 4.6172225e+01, 4.9004151e+01,
     +   5.2035036e+01, 5.3310240e+01, 5.3844529e+01, 5.4066456e+01,
     +   5.4157748e+01, 5.4195184e+01,
C Bss  
     +   6.3524968e+01, 5.9001382e+01, 5.4527924e+01, 5.0411023e+01,
     +   4.7062503e+01, 4.4759402e+01, 4.3486554e+01, 4.6877739e+01,
     +   5.2906769e+01, 5.9804882e+01, 6.6061206e+01, 7.0907372e+01,
     +   7.4346759e+01, 7.6705692e+01, 7.8313134e+01, 7.9408379e+01,
     +   8.0657031e+01, 8.1219565e+01, 8.1461746e+01, 8.1561134e+01,
     +   8.1601113e+01, 8.1617131e+01,
C Byx 
     +   7.3145837e+02, 6.6155009e+02, 5.7857642e+02, 4.8996035e+02,
     +   4.0635842e+02, 3.3557367e+02, 2.1852386e+02, 1.5693340e+02,
     +   1.2250657e+02, 1.0219487e+02, 8.9896350e+01, 8.2388454e+01,
     +   7.7809168e+01, 7.5021818e+01, 7.3331022e+01, 7.2308187e+01,
     +   7.1319160e+01, 7.0961836e+01, 7.0834021e+01, 7.0789102e+01,
     +   7.0773817e+01, 7.0768850e+01,
C Bxs  
     +  -3.5198701e+02,-3.6050550e+02,-3.5875895e+02,-3.4437728e+02,
     +  -3.1856874e+02,-2.8544109e+02,-1.9607861e+02,-1.1918869e+02,
     +  -6.1382029e+01,-2.2022547e+01, 2.0227820e+00, 1.5265326e+01,
     +   2.2023046e+01, 2.5349495e+01, 2.6982317e+01, 2.7803296e+01,
     +   2.8471179e+01, 2.8648335e+01, 2.8782352e+01, 2.8822150e+01,
     +   2.8840177e+01, 2.8804200e+01,
C Bsy  
     +  -3.1592250e+02,-3.2903760e+02,-3.3329724e+02,-3.2569918e+02,
     +  -3.0616116e+02,-2.7781898e+02,-1.9441099e+02,-1.2059040e+02,
     +  -6.5267140e+01,-2.6385142e+01, 2.6130000e-02, 1.7582232e+01,
     +   2.9083386e+01, 3.6552682e+01, 4.1375917e+01, 4.4480400e+01,
     +   4.7755768e+01, 4.9104632e+01, 4.9658630e+01, 4.9884853e+01,       
     +   4.9976544e+01, 5.0013573e+01,
C Bxd 
     +  -3.0048536e+02,-3.1448679e+02,-3.2065163e+02,-3.1598533e+02,
     +  -2.9990396e+02,-2.7470490e+02,-1.9464669e+02,-1.1959461e+02,
     +  -6.2784783e+01,-2.3520276e+01, 2.5131320e+00, 1.9451327e+01,
     +   3.0377200e+01, 3.7405670e+01, 4.1921468e+01, 4.4822711e+01,
     +   4.7886435e+01, 4.9153771e+01, 4.9677103e+01, 4.9891716e+01,
     +   4.9978957e+01, 5.0014327e+01/
*
*  set vv0 term to zero
      vv0 = 0.d0
*
*  spline fits
      if (ifirst .eq. 0) then
*  spline fit of the coefficients for the PE curves
        indel=1
        do ilam=1,5
          call dcopy(85,vl(indel),1,vecel,1)
*  evaluate derivative at first point
          der1=(vecel(2)-vecel(1))/(rr(2)-rr(1))
          call dspline(rr,vecel,85,der1,0d0,csplel(1,ilam))
          indel = indel + 85
        enddo
*  spline fit of the coefficients for the spin-orbit matrix elements
        indso=1
        do ilam=6,12
          call dcopy(22,somat(indso),1,vecso,1)
*  evaluate derivative at first point
          der1=(vecso(2)-vecso(1))/(rso(2)-rso(1))
          call dspline(rso,vecso,22,der1,0d0,csplso(1,ilam))
          indso = indso + 22
        enddo
        ifirst = 1
      end if
*  determine splined coefficients for PE curves at r=R
      indel=1
      do ilam=1,5
        call dcopy(85,vl(indel),1,vecel,1)
        call dsplint(rr,vecel,csplel(1,ilam),85,r,vvx)
        v(ilam) = vvx
        call dcopy(85,vl(indel),1,vecel,1)
        indel = indel + 85
*  shift singlet PE curves to E=0 at large R
        if (ilam .le. 3) then
          v(ilam) = v(ilam) - 15906.642d0
*          if (r .gt. 24.d0) v(ilam) = 0.d0
        endif
*  fix long-range behavior of PE curves to go to zero at large R
        if (r.gt.25.d0) then
          switch_lr = 0.5*(tanh(0.5*(r - 30.d0)) + 1.d0)
          v(ilam) = (1.d0 - switch_lr) * v(ilam)
        endif
      enddo
*  determine splined coefficients for spin-orbit matrix elements at r=R
      indso=1
      do ilam=6,12
        call dcopy(22,somat(indso),1,vecso,1)
        call dsplint(rso,vecso,csplso(1,ilam),22,r,vvx)
        v(ilam) = vvx
        call dcopy(22,somat(indso),1,vecso,1)
        indso = indso + 22
*  if r≥8.5 bohr, set equal to SO matrix element at R=8.5
        if (r.ge.8.5d0) v(ilam) = somat(indso-1)
*  get value of spin-orbit parameter at R=inf
        v(ilam+7) = somat(indso-1)
      enddo
      call dcopy(19,v,1,vvl,1)
* convert to hartree
      econv=1.d0/219474.6d0
      call dscal(19,econv,vvl,1)
      return
      end
*===================================eof====================================
