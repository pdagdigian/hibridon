*  System:  O(1D,3P) + Xe - UCCSD(T)
*
*   calculation of potential energy curves and spin-orbit matrix
*   elements by j.klos 
*
*   in this pot routine, the 3Sigma, 3Pi, 1Sigma, and 1Delta states
*   were obtained in RHF/UCCSD(T) calculations.  the 1Pi state potential
*   energy curve is from the earlier MRCISD+Q calculations
*
*   3P and 1D-3P spin-orbit matrix elements extend to large R.
*   compute 3P spin-orbit splitting using large R values
*  
*   written by p. dagdigian
*   current revision date:  22-oct-2014
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
      potnam='O(1D,3P)-Xe UCCSDT'
      npot=1
      ibasty=22
      lammin(1)=1
      lammax(1)=19
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
      potnam='O(1D,3P)-Xe UCCSDT'
      econv=219474.6d0
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      if (r.le.0.d0) goto 99
      call pot(vv0,r)
      write (6, 100) (econv*vvl(i), i=1,19)
100   format('  ',6(1pe16.8)/'  ',6(1pe16.8)
     +  '  ',7(1pe16.8))
      goto 1
99    rr=2.6d0
      dr=0.25d0
      open (unit=12,file='oxe_uccsdt_vlms.txt')
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
*    vvl(6,7)   contains Axy and Azy, the spin-orbit matrix elements
*               within the 3P state
*    vvl(8,12)  contains Bss, Byx, Bxs, Bsy, and Bxd, the spin-orbit
*               matrix elements coupling the 1D and 3P states
*    vvl(13,19) values of the Axy, Azy, Bss, Byx, Bxs, Bsy, and Bxd
*               matrix elements at R=inf
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
* latest revision date:  22-oct-2014
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
     + 6.3630338e+04, 4.1451440e+04, 2.4913192e+04, 1.2703020e+04,
     + 3.8057737e+03, -2.5622581e+03, -7.0049780e+03, -9.9868054e+03,
     + -1.1864182e+04, -1.2910148e+04, -1.3333565e+04, -1.3294307e+04,
     + -1.2915370e+04, -1.2292412e+04, -1.1501006e+04, -1.0601840e+04,
     + -9.6442341e+03, -8.6683387e+03, -7.7064556e+03, -6.7838228e+03,
     + -5.9191616e+03, -5.1251879e+03, -4.4092136e+03, -3.2180301e+03,
     + -2.7377846e+03, -2.3266904e+03, -1.9777255e+03, -1.6835799e+03,
     + -1.4366560e+03, -1.2296474e+03, -1.0564638e+03, -9.1135933e+02,
     + -7.8952010e+02, -6.8706396e+02, -6.0039667e+02, -5.2691355e+02,
     + -4.6427697e+02, -4.1060203e+02, -3.6445654e+02, -3.2452202e+02,
     + -2.8985604e+02, -2.5961705e+02, -2.3312061e+02, -2.0983944e+02,
     + -1.8928502e+02, -1.7109655e+02, -1.5494751e+02, -1.4056609e+02,
     + -1.2773514e+02, -1.1625172e+02, -1.0596016e+02, -9.6717375e+01,
     + -8.8398869e+01, -8.0898699e+01, -6.5150421e+01, -5.2888120e+01,
     + -4.3218042e+01, -3.5557603e+01, -2.9437390e+01, -2.4529505e+01,
     + -2.0559908e+01, -1.7328473e+01, -1.4669956e+01, -1.2484734e+01,
     + -1.0689588e+01, -9.2013000e+00, -7.9456600e+00, -6.8844900e+00,
     + -5.9886220e+00, -5.2288880e+00, -4.0220190e+00, -3.1349060e+00,
     + -2.4688440e+00, -1.9627110e+00, -1.2670050e+00, -8.5235500e-01,
     + -5.8603000e-01, -4.1139800e-01, -2.9442300e-01, -2.1534400e-01,
     + -1.2022300e-01, -7.0741000e-02, -4.3024000e-02, -2.7735000e-02,
     + -1.8262000e-02,
     + 1.0443048e+05, 8.1394452e+04, 6.3783380e+04, 5.0260751e+04,
     + 3.9802459e+04, 3.1634447e+04, 2.5183372e+04, 2.0039117e+04,
     + 1.5916360e+04, 1.2612667e+04, 9.9759146e+03, 7.8847872e+03,
     + 6.2382349e+03, 4.9501355e+03, 3.9468551e+03, 3.1662141e+03,
     + 2.5568641e+03, 2.0775585e+03, 1.6961268e+03, 1.3881891e+03,
     + 1.1357447e+03, 9.2578618e+02, 7.4904894e+02, 4.7087009e+02,
     + 3.6133476e+02, 2.6779000e+02, 1.8820437e+02, 1.2089897e+02,
     + 6.4428738e+01, 1.7511098e+01, -2.1012539e+01, -5.2198705e+01,
     + -7.7011573e+01, -9.6328356e+01, -1.1094169e+02, -1.2156215e+02,
     + -1.2882163e+02, -1.3327717e+02, -1.3541622e+02, -1.3566237e+02,
     + -1.3438145e+02, -1.3188758e+02, -1.2844905e+02, -1.2429374e+02,
     + -1.1961407e+02, -1.1457135e+02, -1.0929956e+02, -1.0391005e+02,
     + -9.8492160e+01, -9.3118324e+01, -8.7846288e+01, -8.2720083e+01,
     + -7.7773004e+01, -7.3029099e+01, -6.2149682e+01, -5.2725319e+01,
     + -4.4700714e+01, -3.7935472e+01, -3.2254700e+01, -2.7483343e+01,
     + -2.3465852e+01, -2.0073285e+01, -1.7201886e+01, -1.4767939e+01,
     + -1.2702522e+01, -1.0947853e+01, -9.4552180e+00, -8.1838580e+00,
     + -7.0997330e+00, -6.1746490e+00, -4.7110970e+00, -3.6422850e+00,
     + -2.8539610e+00, -2.2611300e+00, -1.4469640e+00, -9.3721000e-01,
     + -6.1407400e-01, -4.0944700e-01, -2.7884500e-01, -1.9392100e-01,
     + -9.7891000e-02, -4.9746000e-02, -2.3471000e-02, -8.3250000e-03,
     + 7.7600000e-04,
     + 1.4164956e+05, 1.1986441e+05, 1.0199441e+05, 8.7448966e+04,
     + 7.5637463e+04, 6.5969296e+04, 5.7889746e+04, 5.1002587e+04,
     + 4.5015400e+04, 3.9678610e+04, 3.4707513e+04, 2.9636762e+04,
     + 2.4903105e+04, 2.0899748e+04, 1.7552123e+04, 1.4785661e+04,
     + 1.2525796e+04, 1.0697960e+04, 9.2275845e+03, 8.0401027e+03,
     + 7.0609468e+03, 6.2155491e+03, 5.4314003e+03, 3.9935293e+03,
     + 3.3811484e+03, 2.8450699e+03, 2.3790989e+03, 1.9770201e+03,
     + 1.6325973e+03, 1.3395905e+03, 1.0917479e+03, 8.8284928e+02,
     + 7.0747193e+02, 5.6099064e+02, 4.3902818e+02, 3.3806074e+02,
     + 2.5479978e+02, 1.8646191e+02, 1.3076883e+02, 8.5587956e+01,
     + 4.9281555e+01, 2.0346339e+01, -2.4744980e+00, -2.0191283e+01,
     + -3.3745761e+01, -4.3848221e+01, -5.1146259e+01, -5.6176429e+01,
     + -5.9364242e+01, -6.1103869e+01, -6.1683421e+01, -6.1362641e+01,
     + -6.0358672e+01, -5.8846056e+01, -5.3715495e+01, -4.7641663e+01,
     + -4.1623850e+01, -3.6032675e+01, -3.1042021e+01, -2.6667502e+01,
     + -2.2878374e+01, -1.9616726e+01, -1.6819983e+01, -1.4434075e+01,
     + -1.2407058e+01, -1.0686989e+01, -9.2251920e+00, -7.9860580e+00,
     + -6.9372480e+00, -6.0464200e+00, -4.6319350e+00, -3.5967970e+00,
     + -2.8271210e+00, -2.2463620e+00, -1.4446730e+00, -9.6363100e-01,
     + -6.5696000e-01, -4.5800700e-01, -3.2489500e-01, -2.3521400e-01,
     + -1.2978200e-01, -7.6004000e-02, -4.6599000e-02, -2.9668000e-02,
     + -1.9505000e-02,
     + 1.5464347e+05, 1.3057827e+05, 1.1178918e+05, 9.7109358e+04,
     + 8.5371918e+04, 7.5558697e+04, 6.7216370e+04, 5.9954016e+04,
     + 5.3495001e+04, 4.7648515e+04, 4.2294558e+04, 3.7369471e+04,
     + 3.2847349e+04, 2.8720835e+04, 2.4986430e+04, 2.1636366e+04,
     + 1.8655905e+04, 1.6023862e+04, 1.3714459e+04, 1.1699374e+04,
     + 9.9494996e+03, 8.4362396e+03, 7.1323814e+03, 5.0538066e+03,
     + 4.2347771e+03, 3.5376078e+03, 2.9452455e+03, 2.4429993e+03,
     + 2.0185412e+03, 1.6602455e+03, 1.3588846e+03, 1.1058870e+03,
     + 8.9397108e+02, 7.1714556e+02, 5.6979249e+02, 4.4756320e+02,
     + 3.4645407e+02, 2.6309929e+02, 1.9477094e+02, 1.3892067e+02,
     + 9.3607680e+01, 5.7055876e+01, 2.7784094e+01, 4.6061160e+00,
     + -1.3582877e+01, -2.7613390e+01, -3.8241674e+01, -4.6093450e+01,
     + -5.1663909e+01, -5.5412612e+01, -5.7679299e+01, -5.8771313e+01,
     + -5.8939801e+01, -5.8379719e+01, -5.4857831e+01, -4.9682972e+01,
     + -4.3970684e+01, -3.8345987e+01, -3.3139031e+01, -2.8468763e+01,
     + -2.4390815e+01, -2.0875700e+01, -1.7871105e+01, -1.5318538e+01,
     + -1.3157961e+01, -1.1329337e+01, -9.7768280e+00, -8.4613960e+00,
     + -7.3482050e+00, -6.4024180e+00, -4.8983600e+00, -3.7938020e+00,
     + -2.9691090e+00, -2.3491640e+00, -1.5108430e+00, -1.0033900e+00,
     + -6.8495000e-01, -4.7886200e-01, -3.4186300e-01, -2.4863500e-01,
     + -1.3799600e-01, -8.0841000e-02, -4.9527000e-02, -3.1502000e-02,
     + -2.0703000e-02,
     + 1.1831505e+05, 9.9613722e+04, 8.0846065e+04, 6.3708831e+04,
     + 4.9898775e+04, 4.0511291e+04, 3.4236322e+04, 2.9207000e+04,
     + 2.4581113e+04, 2.0541111e+04, 1.7243853e+04, 1.4565638e+04,
     + 1.2316097e+04, 1.0384791e+04, 8.7412032e+03, 7.3555209e+03,
     + 6.1868289e+03, 5.1922485e+03, 4.3442236e+03, 3.6216790e+03,
     + 3.0060056e+03, 2.4817428e+03, 2.0361062e+03, 1.3383300e+03,
     + 1.0688620e+03, 8.4268068e+02, 6.5365849e+02, 4.9652995e+02,
     + 3.6679234e+02, 2.6017141e+02, 1.7322654e+02, 1.0293457e+02,
     + 4.6620492e+01, 2.0205390e+00, -3.2819919e+01, -5.9569883e+01,
     + -7.9650238e+01, -9.4265378e+01, -1.0443188e+02, -1.1100427e+02,
     + -1.1469777e+02, -1.1610849e+02, -1.1573112e+02, -1.1397423e+02,
     + -1.1117354e+02, -1.0760369e+02, -1.0348783e+02, -9.9006052e+01,
     + -9.4302591e+01, -8.9491814e+01, -8.4663382e+01, -7.9886684e+01,
     + -7.5214174e+01, -7.0685138e+01, -6.0150011e+01, -5.0886061e+01,
     + -4.2916161e+01, -3.6157103e+01, -3.0477683e+01, -2.5732187e+01,
     + -2.1778590e+01, -1.8487647e+01, -1.5746549e+01, -1.3459467e+01,
     + -1.1546432e+01, -9.9414430e+00, -8.5904170e+00, -7.4491270e+00,
     + -6.4815810e+00, -5.6582010e+00, -4.3497170e+00, -3.3877710e+00,
     + -2.6648800e+00, -2.1174790e+00, -1.3723560e+00, -9.1717900e-01,
     + -6.2937100e-01, -4.4193200e-01, -3.1667500e-01, -2.3104800e-01,
     + -1.2887900e-01, -7.5788000e-02, -4.6564000e-02, -2.9689000e-02,
     + -1.9541000e-02 /
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
     + 2.0028983e+03, 1.8717172e+03, 1.6880283e+03, 1.4705986e+03,
     + 1.2515656e+03, 1.0540202e+03, 6.8657523e+02, 4.5633663e+02,
     + 2.9901971e+02, 2.0995807e+02, 1.5842905e+02, 1.2729100e+02,
     + 1.0810864e+02, 9.6205394e+01, 8.8789232e+01, 8.4156174e+01,
     + 7.9445760e+01, 7.7602781e+01, 7.6882855e+01, 7.6596871e+01,
     + 7.6457772e+01, 7.6425669e+01,
     + -6.9421451e+02, -7.4252306e+02, -7.8459119e+02, -8.2167250e+02,
     + -8.5407790e+02, -8.8045923e+02, -9.0016087e+02, -7.7858487e+02,
     + -4.4368091e+02, -2.4291545e+02, -1.3367381e+02, -6.7442900e+01,
     + -2.5689880e+01, 1.2867770e+00, 1.9025303e+01, 3.0781268e+01,
     + 4.3797579e+01, 4.9557749e+01, 5.2095309e+01, 5.3165425e+01,
     + 5.3682220e+01, 5.3870478e+01,
*
     + 2.4995862e+02, 2.3716130e+02, 2.2107208e+02, 2.0248595e+02,
     + 1.8317523e+02, 1.6491391e+02, 1.2997172e+02, 1.0876538e+02,
     + 8.9580946e+01, 8.3015121e+01, 8.2397485e+01, 8.3276335e+01,
     + 8.3415156e+01, 8.2821179e+01, 8.2130193e+01, 8.1615689e+01,
     + 8.1173338e+01, 8.1136899e+01, 8.1189058e+01, 8.1192303e+01,
     + 8.1149602e+01, 8.1163280e+01,
     + 1.8896122e+03, 1.7511517e+03, 1.5670768e+03, 1.3588484e+03,
     + 1.1561932e+03, 9.7722518e+02, 6.4775454e+02, 4.4066966e+02,
     + 2.9831662e+02, 2.1212484e+02, 1.5924548e+02, 1.2613151e+02,
     + 1.0534871e+02, 9.2326060e+01, 8.4169815e+01, 7.9065829e+01,
     + 7.3886486e+01, 7.1874421e+01, 7.1096256e+01, 7.0778941e+01,
     + 7.0612061e+01, 7.0580349e+01,
     + -8.5387530e+02, -8.8124564e+02, -8.8296336e+02, -8.6115306e+02,
     + -8.2410091e+02, -7.8020960e+02, -6.6516488e+02, -5.4369896e+02,
     + -3.9263364e+02, -2.5960203e+02, -1.5740431e+02, -8.4961459e+01,
     + -3.8242154e+01, -1.0338903e+01, 5.7182330e+00, 1.4925300e+01,
     + 2.3498947e+01, 2.6693668e+01, 2.7982582e+01, 2.8517873e+01,
     + 2.8766647e+01, 2.8911813e+01,
     + -7.6359582e+02, -8.1312324e+02, -8.5695588e+02, -8.9616351e+02,
     + -9.3102330e+02, -9.6028175e+02, -9.8856364e+02, -8.7214444e+02,
     + -5.2233533e+02, -3.0193438e+02, -1.7628052e+02, -9.7565111e+01,
     + -4.6631230e+01, -1.3308540e+01, 8.5210580e+00, 2.2828722e+01,
     + 3.8350340e+01, 4.5003005e+01, 4.7844194e+01, 4.9004504e+01,
     + 4.9536857e+01, 4.9737006e+01,
     + -6.9361133e+02, -7.4280537e+02, -7.8764245e+02, -8.2959833e+02,
     + -8.6882597e+02, -9.0393510e+02, -9.5737478e+02, -9.0181011e+02,
     + -5.7865445e+02, -3.1827478e+02, -1.7512821e+02, -9.1406828e+01,
     + -4.0533266e+01, -8.6877000e+00, 1.1703110e+01, 2.4918723e+01,
     + 3.9184392e+01, 4.5317580e+01, 4.7959444e+01, 4.9073217e+01,
     + 4.9548613e+01, 4.9745891e+01 /
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
