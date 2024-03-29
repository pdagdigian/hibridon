cstart unix-ifort
cdec$ fixedformlinesize:132
cend

*  Multiplied spin-orbit matrix elements (i.e. the difference
*  from their asymtotic values) by 5 X
*  revision:  23-feb-2018 (p.dagdigian)


*  System:  O(1D,3P) + Kr
*  Theory  Level for 3Pi, 3Sigma, 1Sigma, 1Delta curves: 
*                UCCSD(T)/aug-cc-pv5z-dk BSSE and SC corrected
*  Theory  Level for 1Pi curve: 
*       MRCISD+Q(Davidson) aug-cc-pvqz-DK BSSE and SC corrected
*   Scalar relativistic effects included by using all-electron
*   Douglas-Kroll Hamiltonian and integrals
*  Theory Level for SO curves: 
*       MRCISD+Q(Davidson) aug-cc-pvqz-DK
*   calculation of potential energy curves and spin-orbit matrix
*   elements by j.klos 2013 and 2016 (UCCSD(T)) 
*
*   written by p. dagdigian
*   adapted to O-Kr by j. klos
*   past revision date:  1-march-2016 by jklos
*   current revision date:  28-nov2017 by pdagdigian
*   fixed long-range behavior of potentials
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
      potnam='O(1D,3P)-Kr UCCSD(T)'
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
      potnam='O(1D,3P)-Kr UCCSD(T)'
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
      open (unit=12,file='okr_uccsdt_so5X_vlms.txt')
      open (unit=13,file='okr_uccsdt_so5X_pots_so.txt')
      write(12,109)
109   format('% R/bohr 1Sig  1Pi  1Del  3Sig  3Pi  Axy   Azy',
     :  '  Bss  Byx  Bxs  Bsy  Bxd  Axy(R=inf)   Azy(R=inf)',
     :  '  Bss(R=inf)  Byx(R=inf)  Bxs(R=inf)  Bsy(R=inf)',
     :  '  Bxd(R=inf)')
      write(13,119)
119   format(' R/bohr 1Sig  1Pi  1Del  3Sig  3Pi  Axy   Azy',
     :  '  Bss  Byx  Bxs  Bsy  Bxd')
      do i=1,201
        call pot(vv0,rr)
        write (12,110) rr,(econv*vvl(j),j=1,19)
110     format(f7.2,19(1pe16.8))
        write (13,120) rr,(econv*vvl(j),j=1,12)
120     format(f7.2,12(1pe16.8))
        rr = rr + dr
      enddo
      close(12)
      close(13)
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
     +  2.6000000e+00, 2.7000000e+00, 2.8000000e+00, 2.9000000e+00, 3.0000000e+00,
     +  3.1000000e+00, 3.2000000e+00, 3.3000000e+00, 3.4000000e+00, 3.5000000e+00,
     +  3.6000000e+00, 3.7000000e+00, 3.8000000e+00, 3.9000000e+00, 4.0000000e+00,
     +  4.1000000e+00, 4.2000000e+00, 4.3000000e+00, 4.4000000e+00, 4.5000000e+00,
     +  4.6000000e+00, 4.7000000e+00, 4.8000000e+00, 5.0000000e+00, 5.1000000e+00,
     +  5.2000000e+00, 5.3000000e+00, 5.4000000e+00, 5.5000000e+00, 5.6000000e+00,
     +  5.7000000e+00, 5.8000000e+00, 5.9000000e+00, 6.0000000e+00, 6.1000000e+00,
     +  6.2000000e+00, 6.3000000e+00, 6.4000000e+00, 6.5000000e+00, 6.6000000e+00,
     +  6.7000000e+00, 6.8000000e+00, 6.9000000e+00, 7.0000000e+00, 7.1000000e+00,
     +  7.2000000e+00, 7.3000000e+00, 7.4000000e+00, 7.5000000e+00, 7.6000000e+00,
     +  7.7000000e+00, 7.8000000e+00, 7.9000000e+00, 8.0000000e+00, 8.2500000e+00,
     +  8.5000000e+00, 8.7500000e+00, 9.0000000e+00, 9.2500000e+00, 9.5000000e+00,
     +  9.7500000e+00, 1.0000000e+01, 1.0250000e+01, 1.0500000e+01, 1.0750000e+01,
     +  1.1000000e+01, 1.1250000e+01, 1.1500000e+01, 1.1750000e+01, 1.2000000e+01,
     +  1.2500000e+01, 1.3000000e+01, 1.3500000e+01, 1.4000000e+01, 1.5000000e+01,
     +  1.6000000e+01, 1.7000000e+01, 1.8000000e+01, 1.9000000e+01, 2.0000000e+01,
     +  2.2000000e+01, 2.4000000e+01, 2.6000000e+01, 2.8000000e+01, 3.0000000e+01/     
*  values of the 5 PE curves
      data vl /  
C       1Sigma UCCSD(T)
     +  4.2363210e+04, 3.1715699e+04, 2.4002923e+04, 1.8490815e+04, 1.4629209e+04,
     +  1.2006428e+04, 1.0314225e+04, 9.3207036e+03, 8.8495919e+03, 8.7647431e+03,
     +  8.9590339e+03, 9.3468252e+03, 9.8591159e+03, 1.0440545e+04, 1.1047521e+04,
     +  1.1646914e+04, 1.2214949e+04, 1.2736056e+04, 1.3201605e+04, 1.3608531e+04,
     +  1.3957933e+04, 1.4253737e+04, 1.4501529e+04, 1.4878345e+04, 1.5019652e+04,
     +  1.5136814e+04, 1.5234350e+04, 1.5316024e+04, 1.5384903e+04, 1.5443448e+04,
     +  1.5493616e+04, 1.5536948e+04, 1.5574656e+04, 1.5607694e+04, 1.5636815e+04,
     +  1.5662616e+04, 1.5685578e+04, 1.5706088e+04, 1.5724466e+04, 1.5740975e+04,
     +  1.5755837e+04, 1.5769239e+04, 1.5781343e+04, 1.5792287e+04, 1.5802192e+04,
     +  1.5811166e+04, 1.5819300e+04, 1.5826678e+04, 1.5833375e+04, 1.5839455e+04,
     +  1.5844978e+04, 1.5849998e+04, 1.5854563e+04, 1.5858715e+04, 1.5867542e+04,
     +  1.5874546e+04, 1.5880133e+04, 1.5884616e+04, 1.5888238e+04, 1.5891187e+04,
     +  1.5893602e+04, 1.5895592e+04, 1.5897239e+04, 1.5898607e+04, 1.5899748e+04,
     +  1.5900702e+04, 1.5901502e+04, 1.5902178e+04, 1.5902750e+04, 1.5903236e+04,
     +  1.5904010e+04, 1.5904585e+04, 1.5905017e+04, 1.5905346e+04, 1.5905796e+04,
     +  1.5906072e+04, 1.5906248e+04, 1.5906364e+04, 1.5906442e+04, 1.5906496e+04,
     +  1.5906560e+04, 1.5906594e+04, 1.5906612e+04, 1.5906623e+04, 1.5906642e+04,
C      1Pi  MRCISD+Q
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
C      1Delta  UCCSD(T)
     +  1.2301518e+05, 1.1102502e+05, 1.0048505e+05, 9.1015516e+04, 8.2372428e+04,
     +  7.4421887e+04, 6.7113533e+04, 6.0441997e+04, 5.4412636e+04, 4.9021738e+04,
     +  4.4249962e+04, 4.0063492e+04, 3.6418251e+04, 3.3264479e+04, 3.0550545e+04,
     +  2.8225711e+04, 2.6241921e+04, 2.4554844e+04, 2.2827265e+04, 2.1817954e+04,
     +  2.0901916e+04, 2.0093974e+04, 1.9393560e+04, 1.8283147e+04, 1.7852345e+04,
     +  1.7490478e+04, 1.7187911e+04, 1.6935985e+04, 1.6727056e+04, 1.6554464e+04,
     +  1.6412456e+04, 1.6296095e+04, 1.6201171e+04, 1.6124111e+04, 1.6061888e+04,
     +  1.6011952e+04, 1.5972161e+04, 1.5940717e+04, 1.5916118e+04, 1.5897112e+04,
     +  1.5882657e+04, 1.5871891e+04, 1.5864095e+04, 1.5858681e+04, 1.5855164e+04,
     +  1.5853143e+04, 1.5852295e+04, 1.5852351e+04, 1.5853098e+04, 1.5854360e+04,
     +  1.5855998e+04, 1.5857899e+04, 1.5859973e+04, 1.5862152e+04, 1.5867722e+04,
     +  1.5873070e+04, 1.5877930e+04, 1.5882206e+04, 1.5885893e+04, 1.5889032e+04,
     +  1.5891684e+04, 1.5893914e+04, 1.5895785e+04, 1.5897355e+04, 1.5898672e+04,
     +  1.5899779e+04, 1.5900713e+04, 1.5901502e+04, 1.5902171e+04, 1.5902740e+04,
     +  1.5903642e+04, 1.5904308e+04, 1.5904808e+04, 1.5905185e+04, 1.5905699e+04,
     +  1.5906012e+04, 1.5906210e+04, 1.5906339e+04, 1.5906425e+04, 1.5906484e+04,
     +  1.5906554e+04, 1.5906590e+04, 1.5906610e+04, 1.5906622e+04, 1.5906642e+04,
C     3Sigma UCCSD(T)
     +  1.2815304e+05, 1.1379656e+05, 1.0132958e+05, 9.0149067e+04, 7.9918253e+04,
     +  7.0478712e+04, 6.1778150e+04, 5.3815430e+04, 4.6601509e+04, 4.0136813e+04,
     +  3.4402879e+04, 2.9363139e+04, 2.4967767e+04, 2.1159270e+04, 1.7877256e+04,
     +  1.5061926e+04, 1.2656365e+04, 1.0607881e+04, 8.8686520e+03, 7.3959179e+03,
     +  6.1518684e+03, 5.1033627e+03, 4.2215530e+03, 2.8615931e+03, 2.3434593e+03,
     +  1.9112662e+03, 1.5515330e+03, 1.2527855e+03, 1.0052783e+03, 8.0074984e+02,
     +  6.3220884e+02, 4.9374828e+02, 3.8038529e+02, 2.8792264e+02, 2.1283021e+02,
     +  1.5214379e+02, 1.0337782e+02, 6.4451456e+01, 3.3625549e+01, 9.4488510e+00,
     + -9.2875680e+00,-2.3588933e+01,-3.4289586e+01,-4.2080511e+01,-4.7532678e+01,
     + -5.1116885e+01,-5.3220193e+01,-5.4160063e+01,-5.4195756e+01,-5.3538610e+01,
     + -5.2359956e+01,-5.0798095e+01,-4.8964184e+01,-4.6946957e+01,-4.1527009e+01,
     + -3.6100178e+01,-3.1043115e+01,-2.6519629e+01,-2.2574038e+01,-1.9187173e+01,
     + -1.6309445e+01,-1.3879663e+01,-1.1835467e+01,-1.0118371e+01,-8.6761070e+00,
     + -7.4634330e+00,-6.4417760e+00,-5.5788840e+00,-4.8475520e+00,-4.2265330e+00,
     + -3.2431940e+00,-2.5183150e+00,-1.9766470e+00,-1.5673730e+00,-1.0120140e+00,
     + -6.7416800e-01,-4.6135300e-01,-3.2320100e-01,-2.3113700e-01,-1.6835100e-01,
     + -9.3646000e-02,-5.4952000e-02,-3.3701000e-02,-2.1465000e-02,-1.4099000e-02,
C       3Pi UCCSD(T)
     +  8.1767549e+04, 6.8841639e+04, 5.8292023e+04, 4.9481742e+04, 4.2002624e+04,
     +  3.5599161e+04, 3.0104045e+04, 2.5394767e+04, 2.1370489e+04, 1.7942164e+04,
     +  1.5029378e+04, 1.2559919e+04, 1.0469896e+04, 8.7036213e+03, 7.2130883e+03,
     +  5.9571953e+03, 4.9008640e+03, 4.0141789e+03, 3.2716103e+03, 2.6513474e+03,
     +  2.1347379e+03, 1.7058229e+03, 1.3509498e+03, 8.1835503e+02, 6.2217932e+02,
     +  4.6270015e+02, 3.3379044e+02, 2.3026412e+02, 1.4774279e+02, 8.2539872e+01,
     +  3.1559676e+01,-7.7898320e+00,-3.7672528e+01,-5.9889441e+01,-7.5934839e+01,
     + -8.7043972e+01,-9.4234728e+01,-9.8342576e+01,-1.0005052e+02,-9.9914524e+01,
     + -9.8385037e+01,-9.5825082e+01,-9.2525600e+01,-8.8718105e+01,-8.4585347e+01,
     + -8.0270388e+01,-7.5883857e+01,-7.1510246e+01,-6.7212725e+01,-6.3037690e+01,
     + -5.9017963e+01,-5.5175582e+01,-5.1524333e+01,-4.8071447e+01,-4.0312642e+01,
     + -3.3747668e+01,-2.8256211e+01,-2.3693761e+01,-1.9916578e+01,-1.6793317e+01,
     + -1.4209547e+01,-1.2068479e+01,-1.0289725e+01,-8.8073670e+00,-7.5676430e+00,
     + -6.5269790e+00,-5.6500600e+00,-4.9082000e+00,-4.2780020e+00,-3.7414020e+00,
     + -2.8876810e+00,-2.2542890e+00,-1.7780050e+00,-1.4160890e+00,-9.2130800e-01,
     + -6.1760800e-01,-4.2485800e-01,-2.9894400e-01,-2.1458500e-01,-1.5678900e-01,
     + -8.7658000e-02,-5.1634000e-02,-3.1761000e-02,-2.0276000e-02,-1.3350000e-02/
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
C Axy
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
        endif
*
*  fix long-range behavior of PE curves to go to zero as R^{-6} at large R
        if (r.gt.20.d0) then
          vvv = vl(indel-4)
          if (ilam .le. 3) vvv = (vvv - 15906.642d0)
          c6 = vvv * rr(81)**6
          switch_lr = 0.5*(tanh(0.5*(r - 25.d0)) + 1.d0)
          v(ilam) = (1.d0 - switch_lr) * v(ilam) 
     :      + switch_lr * c6 / (r**6)
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
*
*  Multiply spin-orbit matrix elements (i.e. the difference
*  from their asymtotic values) by 5 X
      do ilam=6,12
        dso = v(ilam) - v(ilam+7)
        v(ilam) = v(ilam+7) + 5.d0 * dso
      enddo
*   
      call dcopy(19,v,1,vvl,1)
* convert to hartree
      econv=1.d0/219474.6d0
      call dscal(19,econv,vvl,1)
      return
      end
*===================================eof====================================
