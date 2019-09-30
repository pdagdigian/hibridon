*  System:  H-CO - Re(CO) - Nijmegen PES
*
*   L. Song, A. van der Avoird, and G. C. Groenenboom,
*   JPCA 1117, 7571 (2013)
*
*   lambda_max of angular fit = 18
*  
*   written by p. dagdigian
*   current revision date:  21-jan-2015
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(18)
      include "common/parpot"
      econv=219474.6d0
      potnam='HCO PES-Nijmegen-lmax=18'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      call pot(vv0,r)
      if (r.le.0.d0) goto 99
*     vlm coefficient is returned in atomic units (hartree)
*     convert from atomic units for printout
      write (6, 100) vv0*econv, (econv*vvl(i), i=1,18)
100   format(19(1pe16.8))
      goto 1
99    rr=3.0d0
      dr=0.2d0
      open (unit=12,file='hco_vlms.txt')
      write(12,109)
109   format(' %R/bohr V0  V1  V2  V3  V4  V5  V6  V7  V8',
     +  '  V9  V10  V11  V12  V13  V14  V15  V16  V17  V18')
      do i=1,250
        call pot(vv0,rr)
        write (12,110) rr,vv0*econv, (econv*vvl(j),j=1,18)
110     format(f7.2,19(1pe16.8))
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
      potnam='HCO PES-Nijmegen-lmax=18'
      npot=1
      nterm=1
      lammin(1)=1
      lammax(1)=18
      mproj(1)=0
      ipotsy = 1
*
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
*  count number of anisotropic terms
      nlam = 0
      do il=lammin(1),lammax(1),ipotsy
        nlam = nlam + 1
      enddo
      nlammx = nlam
      return
      end
* --------------------------------------------------------------------------
      subroutine pot (vv0, r)


*  subroutine to calculate the r-dependent coefficients in the
*  collision of a diatomic with a structureless target
*  in units of hartree for energy and bohr for distance

*  on return:
*  vv0 contains the isotropic term (n=0) in the potential
*  the coefficients for each angular term in the coupling potential
*  [ vvl(i) for i = 1, nlam ] are returned in common block /covvl/

*  variable in common block /conlam/
*    nlammx:    the maximum number of anisotropic terms
*    nlam:      the total number of angular coupling terms
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential
* 
* author:  paul dagdigian
* latest revision date:  21-jan-2015
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      dimension v(19)
      dimension csplin(47,19)
      dimension rr(47), vl(893),vec(47)
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(18)
*
*  47 values or R
      data rr /
     +  2.000,  2.250,  2.500,  2.750,  3.000, 
     +  3.250,  3.500,  3.750,  4.000,  4.250, 
     +  4.500,  4.750,  5.000,  5.250,  5.500, 
     +  5.750,  6.000,  6.250,  6.500,  6.750, 
     +  7.000,  7.500,  8.000,  8.500,  9.000, 
     +  9.500, 10.000, 11.000, 12.000, 13.000, 
     + 14.000, 15.000, 16.000, 17.000, 18.000, 
     + 19.000, 20.000, 21.000, 22.000, 23.000, 
     + 24.000, 25.000, 26.000, 27.000, 28.000, 
     + 29.000, 30.000 /
*  19 values of the vlam coefficients
      data vl /  
     + 3.9706079e+04, 2.1692465e+04, 1.3338522e+04, 8.8301095e+03,
     + 6.3727160e+03, 4.7349829e+03, 3.6541442e+03, 2.8961201e+03,
     + 2.2638124e+03, 1.6936795e+03, 1.2110582e+03, 8.3186609e+02,
     + 5.4899488e+02, 3.4648758e+02, 2.0639953e+02, 1.1255636e+02,
     + 5.1768792e+01, 1.3994132e+01, -8.2599619e+00, -2.0297316e+01,
     + -2.5828449e+01, -2.6658263e+01, -2.2194468e+01, -1.6990090e+01,
     + -1.2549664e+01, -9.1543559e+00, -6.6746180e+00, -3.6173559e+00,
     + -2.0544957e+00, -1.2232495e+00, -7.6035081e-01, -4.9034402e-01,
     + -3.2638179e-01, -2.2327457e-01, -1.5642863e-01, -1.1190710e-01,
     + -8.1531851e-02, -6.0360238e-02, -4.5334743e-02, -3.4498636e-02,
     + -2.6568637e-02, -2.0687111e-02, -1.6270899e-02, -1.2917130e-02,
     + -1.0343374e-02, -8.3489375e-03, -6.7894079e-03, 
     + -2.7045095e+04, -6.3707818e+03, 5.2200675e+03, 1.0075601e+04,
     + 1.0442431e+04, 8.1040304e+03, 5.1909896e+03, 2.7435126e+03,
     + 1.1338101e+03, 2.8897951e+02, -8.7140686e+01, -2.2292633e+02,
     + -2.4332160e+02, -2.1495690e+02, -1.7092989e+02, -1.2678940e+02,
     + -8.8894339e+01, -5.8987967e+01, -3.6760648e+01, -2.1014004e+01,
     + -1.0361351e+01, 7.4916135e-01, 4.3101603e+00, 4.6539288e+00,
     + 3.9639196e+00, 3.0739172e+00, 2.2698991e+00, 1.1760625e+00,
     + 6.1598366e-01, 3.3401237e-01, 1.8923984e-01, 1.1216644e-01,
     + 6.9500922e-02, 4.4510651e-02, 2.9079811e-02, 1.9314530e-02,
     + 1.3073547e-02, 9.0478868e-03, 6.3938025e-03, 4.6024940e-03,
     + 3.3680392e-03, 2.5013485e-03, 1.8825874e-03, 1.4341094e-03,
     + 1.1045661e-03, 8.5936837e-04, 6.7483013e-04, 
     + 7.6828350e+04, 2.3954712e+04, 3.8081192e+03, -1.1879020e+03,
     + -1.4983996e+03, -9.3883989e+02, 2.9341464e+01, 9.8506703e+02,
     + 1.4909031e+03, 1.5059373e+03, 1.2863814e+03, 1.0121673e+03,
     + 7.5630670e+02, 5.4483227e+02, 3.8091617e+02, 2.5909926e+02,
     + 1.7144174e+02, 1.0995106e+02, 6.7838252e+01, 3.9732784e+01,
     + 2.1471357e+01, 2.9510447e+00, -3.1953182e+00, -4.3765050e+00,
     + -3.8980766e+00, -3.0419106e+00, -2.2405106e+00, -1.1371103e+00,
     + -5.8493497e-01, -3.1365491e-01, -1.7709667e-01, -1.0468735e-01,
     + -6.4428394e-02, -4.1178953e-02, -2.7303620e-02, -1.8696354e-02,
     + -1.3127572e-02, -9.3837094e-03, -6.8112166e-03, -5.0146525e-03,
     + -3.7408623e-03, -2.8249173e-03, -2.1575986e-03, -1.6654217e-03,
     + -1.2982502e-03, -1.0213861e-03, -8.1051597e-04, 
     + -4.5844892e+04, -2.6276376e+04, -1.4206360e+04, -4.9733613e+03,
     + 1.5086698e+02, 1.5168254e+03, 1.1485859e+03, 3.3190646e+02,
     + -2.5277771e+02, -4.3738094e+02, -4.2900782e+02, -3.6423239e+02,
     + -2.8878071e+02, -2.2024340e+02, -1.6341357e+02, -1.1847744e+02,
     + -8.4232802e+01, -5.8673488e+01, -4.0040727e+01, -2.6774299e+01,
     + -1.7513882e+01, -6.9048586e+00, -2.3012921e+00, -5.0656542e-01,
     + 1.1830803e-01, 2.7587746e-01, 2.6096867e-01, 1.4280830e-01,
     + 6.7101789e-02, 3.0795257e-02, 1.4820059e-02, 7.3937490e-03,
     + 3.7036374e-03, 1.7782136e-03, 7.3093262e-04, 1.7501631e-04,
     + -8.8946002e-05, -1.8426745e-04, -2.0123126e-04, -1.8657759e-04,
     + -1.6138538e-04, -1.3484104e-04, -1.1061663e-04, -8.9875191e-05,
     + -7.2700549e-05, -5.8728266e-05, -4.7469460e-05, 
     + 7.0536713e+04, 3.4865858e+04, 1.7496875e+04, 9.0960114e+03,
     + 5.6624429e+03, 3.5241523e+03, 2.3968667e+03, 1.9140574e+03,
     + 1.5438328e+03, 1.1150497e+03, 7.5687258e+02, 5.0728988e+02,
     + 3.3935262e+02, 2.2775348e+02, 1.5324880e+02, 1.0315733e+02,
     + 6.9323934e+01, 4.6381497e+01, 3.0829814e+01, 2.0324820e+01,
     + 1.3269177e+01, 5.4580807e+00, 2.1041570e+00, 7.4254781e-01,
     + 2.1030872e-01, 2.0325230e-02, -3.2150156e-02, -2.9618913e-02,
     + -1.3950520e-02, -5.6282885e-03, -2.3331774e-03, -1.0263191e-03,
     + -4.6145850e-04, -1.2877011e-04, 1.2751191e-04, 2.7965687e-04,
     + 3.2371777e-04, 2.9280515e-04, 2.3860884e-04, 1.8566528e-04,
     + 1.4131471e-04, 1.0646668e-04, 7.9926775e-05, 6.0011362e-05,
     + 4.5169491e-05, 3.4130693e-05, 2.5912028e-05, 
     + -2.4318679e+04, -7.7854036e+03, -2.9961712e+03, -2.5743255e+03,
     + -1.7384551e+03, -1.4785522e+03, -1.3825395e+03, -1.2940559e+03,
     + -1.0549653e+03, -6.9712752e+02, -4.1607710e+02, -2.4391427e+02,
     + -1.4373288e+02, -8.6458084e+01, -5.3298339e+01, -3.3650683e+01,
     + -2.1568231e+01, -1.4038892e+01, -9.2293027e+00, -6.0679187e+00,
     + -3.9711899e+00, -1.6742449e+00, -6.7513404e-01, -2.5621251e-01,
     + -8.5719478e-02, -2.1296803e-02, -8.5006297e-04, 4.6745055e-03,
     + 3.2970118e-03, 1.5964149e-03, 3.4300639e-04, -3.7167655e-04,
     + -6.1249039e-04, -4.1473476e-04, 3.2774915e-05, 3.7550117e-04,
     + 4.9960947e-04, 4.5720713e-04, 3.6435437e-04, 2.7418145e-04,
     + 2.0092337e-04, 1.4548374e-04, 1.0489606e-04, 7.5640862e-05,
     + 5.4697813e-05, 3.9725254e-05, 2.9005560e-05, 
     + 2.7317960e+04, 1.1358722e+04, 6.1646790e+03, 3.7986854e+03,
     + 2.8478145e+03, 2.1643563e+03, 1.7482084e+03, 1.4785367e+03,
     + 1.1021416e+03, 6.7367877e+02, 3.7275453e+02, 2.0169258e+02,
     + 1.0917228e+02, 5.9994428e+01, 3.3709628e+01, 1.9497687e+01,
     + 1.1394802e+01, 6.8498048e+00, 4.2355093e+00, 2.6383203e+00,
     + 1.6403321e+00, 6.4654650e-01, 2.5252122e-01, 9.4773297e-02,
     + 3.3710626e-02, 1.0634743e-02, 2.4118129e-03, -1.0149970e-03,
     + -1.0532034e-03, -3.1951223e-04, 1.8795630e-04, -4.2881301e-05,
     + -4.9164639e-04, -5.3598360e-04, -1.9679424e-04, 1.4323867e-04,
     + 3.1472476e-04, 3.3228651e-04, 2.8807977e-04, 2.3117749e-04,
     + 1.7888460e-04, 1.3595502e-04, 1.0246356e-04, 7.6995177e-05,
     + 5.7870810e-05, 4.3596572e-05, 3.2957352e-05, 
     + -1.4504828e+04, -8.1246406e+03, -3.8533466e+03, -9.5826440e+02,
     + -3.4195994e+02, -5.5052637e+02, -8.0455954e+02, -9.0987013e+02,
     + -7.4177471e+02, -4.4394279e+02, -2.3209465e+02, -1.1700661e+02,
     + -5.8341292e+01, -2.9267115e+01, -1.4807777e+01, -7.7898424e+00,
     + -4.0541426e+00, -2.1733090e+00, -1.2336156e+00, -7.1560281e-01,
     + -4.1048134e-01, -1.4695438e-01, -5.1553686e-02, -1.6564787e-02,
     + -5.5699016e-03, -1.7983901e-03, -2.3619374e-04, 3.8235571e-04,
     + 1.1994948e-04, 7.6850614e-05, 1.8572337e-04, 2.7893024e-04,
     + 3.0335102e-04, 2.3544137e-04, 1.1060120e-04, 9.9990760e-06,
     + -3.9744132e-05, -4.9786187e-05, -4.4311923e-05, -3.5276133e-05,
     + -2.6734961e-05, -1.9794853e-05, -1.4495126e-05, -1.0572311e-05,
     + -7.7076953e-06, -5.6346653e-06, -4.1343342e-06, 
     + 1.0997692e+04, 2.7895059e+03, 5.8581855e+02, 7.1888118e+02,
     + 7.2066682e+02, 5.5007254e+02, 5.8397361e+02, 6.5570020e+02,
     + 5.3787725e+02, 3.0743284e+02, 1.5208908e+02, 7.2844511e+01,
     + 3.4546857e+01, 1.6495377e+01, 7.8913469e+00, 3.9234819e+00,
     + 1.9400839e+00, 9.7083143e-01, 5.0180486e-01, 2.6359644e-01,
     + 1.3533750e-01, 4.1471461e-02, 1.1688682e-02, 1.9844699e-03,
     + 7.9312605e-05, 3.3946401e-05, 2.5460174e-04, 4.7812449e-04,
     + 3.3946543e-04, 1.1149069e-05, -2.6043679e-04, -3.5454994e-04,
     + -3.4030263e-04, -2.8615319e-04, -2.2018561e-04, -1.6166393e-04,
     + -1.1486761e-04, -7.8435432e-05, -5.2413947e-05, -3.4778189e-05,
     + -2.3089078e-05, -1.5397041e-05, -1.0337151e-05, -6.9988358e-06,
     + -4.7785629e-06, -3.2917448e-06, -2.2878555e-06, 
     + -4.5740889e+03, 2.4495669e+02, 1.4254599e+02, -2.9878465e+02,
     + -1.4476989e+01, -9.9053253e+01, -2.5785501e+02, -3.8919246e+02,
     + -3.4573848e+02, -1.8883736e+02, -8.7525197e+01, -3.9354525e+01,
     + -1.7582041e+01, -7.9297311e+00, -3.6207581e+00, -1.6857366e+00,
     + -7.9016515e-01, -3.8777561e-01, -1.8847388e-01, -8.6451624e-02,
     + -4.2380856e-02, -1.1188884e-02, -1.4641013e-03, 2.9126628e-04,
     + 3.7310140e-04, 4.7634704e-04, 6.1777888e-04, 2.2417599e-04,
     + -3.2287837e-04, -2.0511135e-04, 3.6577380e-05, 3.8536983e-05,
     + -6.6550547e-05, -1.3438694e-04, -1.5501111e-04, -1.4532523e-04,
     + -1.1872862e-04, -8.5710614e-05, -5.8104750e-05, -3.8386942e-05,
     + -2.5125088e-05, -1.6428924e-05, -1.0782593e-05, -7.1182620e-06,
     + -4.7406875e-06, -3.1793083e-06, -2.1524270e-06, 
     + 4.1721205e+03, 2.9268391e+03, 1.5000053e+03, 3.3179147e+02,
     + 3.3015388e+02, 2.1301720e+02, 2.1235077e+02, 2.6264603e+02,
     + 2.3291091e+02, 1.1837977e+02, 5.0814602e+01, 2.1188893e+01,
     + 8.8786834e+00, 3.8071263e+00, 1.7358373e+00, 7.2972543e-01,
     + 3.1193690e-01, 1.5069732e-01, 7.1962546e-02, 3.7076175e-02,
     + 2.7830636e-02, 2.5468723e-03, -6.3344530e-04, 8.4603277e-04,
     + 6.7208179e-04, 1.9541463e-04, -1.1887698e-04, 1.4951587e-04,
     + 5.1068676e-04, 7.4512040e-05, -1.2011826e-04, 4.4884282e-04,
     + 1.0242357e-03, 1.0092809e-03, 6.0544335e-04, 2.5321940e-04,
     + 6.8350154e-05, 4.1848523e-06, -1.2199740e-05, -1.3595673e-05,
     + -1.1054541e-05, -8.0929470e-06, -5.6601647e-06, -3.8754408e-06,
     + -2.6265323e-06, -1.7790250e-06, -1.2063208e-06, 
     + -2.3104104e+03, -1.5818041e+03, -5.3005660e+02, -5.6854141e+01,
     + 3.4177486e+00, -8.6249084e+01, -1.1570337e+02, -1.5748875e+02,
     + -1.4415748e+02, -6.8732907e+01, -2.7385809e+01, -1.0547047e+01,
     + -4.1177530e+00, -1.6788211e+00, -7.5917765e-01, -2.6123912e-01,
     + -1.0625099e-01, -3.2480974e-02, -7.0005451e-03, -1.6838147e-02,
     + -2.7624849e-02, -6.6196756e-04, -2.6238375e-04, -2.9258465e-03,
     + -1.8734274e-03, -1.0420800e-04, 1.1230711e-03, 5.8627437e-04,
     + -5.9322223e-04, 4.8993801e-05, 5.5381773e-04, -2.5370807e-05,
     + -7.2891654e-04, -8.5185970e-04, -5.9061627e-04, -3.2900324e-04,
     + -1.6934369e-04, -8.7507562e-05, -4.6284605e-05, -2.5047232e-05,
     + -1.3849492e-05, -7.8147102e-06, -4.4991590e-06, -2.6395852e-06,
     + -1.5758266e-06, -9.5430228e-07, -5.8771984e-07, 
     + 7.0853835e+02, -1.8366075e+03, -8.4217791e+02, 9.0541292e+01,
     + 8.1540236e+01, 7.6503840e+01, 7.1351488e+01, 7.8829407e+01,
     + 7.1976372e+01, 3.2705363e+01, 1.2355112e+01, 4.4789966e+00,
     + 1.6487651e+00, 6.5847882e-01, 2.4662746e-01, 7.3537709e-02,
     + 4.3729075e-02, -5.2073031e-03, -1.5009795e-02, 6.2147080e-03,
     + 2.0127722e-02, 6.5902430e-04, 1.9293412e-04, 1.9779877e-04,
     + -1.0856008e-03, -1.3083701e-03, -8.7113279e-04, 2.6588653e-04,
     + 4.6436925e-04, -7.4750191e-04, -1.3977741e-03, -7.9911133e-04,
     + 4.4791838e-05, 4.9777281e-04, 6.5219456e-04, 6.2128751e-04,
     + 4.8777834e-04, 3.2147676e-04, 1.9505255e-04, 1.1475554e-04,
     + 6.6887820e-05, 3.9024361e-05, 2.2916753e-05, 1.3583649e-05,
     + 8.1328930e-06, 4.9292808e-06, 3.0227969e-06, 
     + -5.5076107e+01, 5.7512437e+00, -4.6384947e-01, -3.5690678e+00,
     + 5.3767999e-01, 5.2248414e-04, 8.3292872e-04, -1.6192717e-04,
     + 1.3221838e-04, -4.2631992e-05, -2.8621009e-04, 2.0306945e-04,
     + 1.0018538e-05, 6.9260047e-05, -1.4197698e-05, 1.4544927e-05,
     + -1.0007609e-05, -1.0044053e-05, -1.1812599e-06, -1.9529587e-06,
     + -1.6203141e-06, 2.4289379e-06, 2.8998226e-06, -2.6359078e-06,
     + -4.3099360e-08, 3.9908420e-08, 2.9034604e-07, -1.1062500e-08,
     + -1.7779856e-07, -2.7827855e-07, 1.4749176e-07, 1.3278105e-08,
     + 1.4974804e-08, 3.5654881e-09, -1.8870062e-08, -4.5026339e-08,
     + -3.1118043e-09, 3.3449592e-09, -2.7043901e-09, 1.4650868e-09,
     + 4.8693227e-09, -3.0163024e-10, -3.7098245e-10, 1.1169042e-09,
     + 2.9289312e-09, 7.7148573e-11, 2.7134337e-10, 
     + 2.3981660e+01, -1.5340476e+01, -6.5993425e-01, -2.8517891e+00,
     + 5.2794976e-01, 4.9049448e-04, 1.0785882e-03, 4.0612564e-04,
     + 1.2691058e-04, 1.1862927e-04, -1.0359655e-04, -4.5057613e-05,
     + -4.5457109e-05, -4.1925174e-05, 1.1255996e-05, -2.2244613e-05,
     + -3.4257454e-06, 1.5951042e-06, 1.1132122e-05, -3.2460837e-06,
     + -4.9337885e-06, -6.2274649e-07, 2.5325822e-06, 2.1832820e-06,
     + -4.5577521e-06, -1.0024862e-06, -1.0359690e-06, 4.1991125e-07,
     + 3.8173700e-07, 1.6502271e-07, -1.4147477e-07, -2.9923247e-08,
     + -4.2006530e-09, 2.0477383e-09, -2.4041101e-08, 1.3826542e-08,
     + 7.6463667e-09, -2.1994870e-09, -1.3011088e-09, 6.4747428e-10,
     + 1.3483812e-09, 1.1470789e-09, -7.0836002e-09, -1.8756976e-09,
     + -1.2856591e-09, -5.3685268e-11, 4.3699716e-10, 
     + 1.3064922e+01, 6.4056154e-01, -2.6953571e+00, -2.4754882e+00,
     + 6.0635294e-01, 6.1446961e-04, 5.3371588e-04, 3.6757636e-04,
     + 9.9567380e-05, -6.8550615e-05, 1.0125145e-04, 2.7190699e-04,
     + 3.1797649e-05, -3.4345742e-05, 2.6499125e-05, 1.9906209e-05,
     + 1.6854579e-05, -3.7250941e-06, -7.7605941e-06, -1.8418123e-07,
     + 1.7509851e-06, -5.6408107e-06, 1.0890245e-06, 5.4940370e-06,
     + 3.9367952e-06, 1.4882233e-06, 4.2876188e-07, -3.4311826e-08,
     + -8.4278928e-08, -2.0293350e-08, 7.5364759e-08, -3.5317007e-08,
     + -2.9477039e-09, -1.3030313e-08, -6.2358770e-09, 2.1150638e-08,
     + -1.9146238e-09, 2.1663723e-10, 8.1489553e-10, 9.7772909e-10,
     + -9.7764008e-10, -3.9249172e-09, -3.3109588e-10, 1.9754120e-09,
     + 3.7161355e-10, 2.3929834e-10, 2.3136818e-10, 
     + -9.7779098e+00, -7.8677822e-01, 7.5630234e-01, -1.6972033e+00,
     + 5.7202273e-01, 3.4234805e-05, 1.0629101e-04, 3.4405457e-05,
     + -5.4590435e-05, 3.4100641e-04, -2.4272392e-04, -1.4728571e-04,
     + 4.1004274e-05, -7.2645792e-06, -2.2636581e-05, -3.4429246e-06,
     + -2.0766994e-05, -6.7114751e-06, 7.5935902e-06, 1.8172724e-06,
     + 1.6763597e-06, 1.9302567e-06, 6.4475998e-07, -2.6257622e-06,
     + 6.2136672e-07, -1.3998971e-06, -3.2090001e-07, -5.5239841e-07,
     + -3.2171168e-07, 9.4658600e-08, 1.7119082e-08, 2.9254262e-08,
     + 5.4858113e-08, -5.2703399e-08, 8.5127270e-09, 1.1232014e-09,
     + 2.5871826e-09, 5.3785303e-10, 1.7361391e-09, -4.4214363e-09,
     + -8.5271518e-10, -2.2216000e-10, -3.8878121e-10, -2.3881929e-10,
     + 2.6564858e-09, -4.9750043e-10, -1.2919072e-10, 
     + 9.1215222e+00, 1.8741956e+01, 9.9075414e-02, -1.6378814e+00,
     + 9.0078784e-01, -7.0006052e-04, 5.8333849e-04, -8.6917074e-06,
     + 2.3522476e-04, 7.1300855e-05, 3.3750433e-04, -9.8742723e-05,
     + 2.2371267e-05, -1.1122742e-04, 9.8796575e-05, -7.5143405e-06,
     + 2.4895172e-05, 1.8625610e-05, -7.1926544e-06, 2.5464359e-06,
     + -1.3986551e-06, 4.8492777e-06, -5.9243431e-06, 1.3883503e-06,
     + -1.4794466e-06, 2.3732852e-07, 1.0312876e-06, 3.0810255e-07,
     + 1.1355078e-07, -3.4562057e-08, -1.4335744e-07, 6.4484133e-09,
     + -1.0761133e-08, 1.3592861e-09, 3.5277316e-08, 1.4271277e-09,
     + -1.5972010e-09, -8.5268080e-10, 1.0429676e-09, -2.2397648e-10,
     + -4.1472842e-09, 1.2469876e-09, 1.0817390e-09, -6.4445445e-09,
     + 1.2939930e-09, 8.2052034e-10, -2.1959283e-10, 
     + 6.8232642e+00, -1.5893785e+01, -5.2331981e-01, -6.4870901e-01,
     + 8.0602392e-01, -3.0476473e-04, 4.8463534e-04, -2.2600039e-04,
     + -2.8099808e-04, -4.6676858e-04, -1.6662500e-04, -1.0966190e-04,
     + 1.5903127e-04, 1.2294094e-04, -1.2978607e-04, 4.2432127e-05,
     + -8.6985535e-06, -1.2116257e-05, -9.7504383e-07, 1.7538501e-06,
     + 1.7603551e-06, -1.2341054e-06, -2.2806276e-06, 1.1244381e-06,
     + -6.8642954e-07, 3.4376466e-07, -4.0782589e-07, 1.1526134e-07,
     + -2.7663482e-08, -5.4985245e-07, 1.8242179e-07, -1.8973306e-08,
     + -1.4983914e-08, 2.2622658e-08, 1.6644059e-08, 1.7560664e-08,
     + 1.6437289e-09, -1.4781804e-09, -3.1344970e-09, 6.6125183e-09,
     + 9.2607039e-10, 1.9781657e-09, 4.3685992e-09, -2.1530878e-09,
     + 2.5492652e-10, -9.7628903e-10, -2.5405791e-10 /
*
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         ind=1
         do ilam=1,19
           call dcopy(47,vl(ind),1,vec,1)
*    evaluate derivative at first point
           der1=(vec(2)-vec(1))/(rr(2)-rr(1))
           call dspline(rr,vec,47,der1,0d0,csplin(1,ilam))
           ind = ind + 47
         enddo
         ifirst = 1
       end if
* r^-6 fit to at R = 25 bohr for isotropic part of potential
       c6sum = -5.05056421875e+06
* switching function for long-range
       switch_lr=0.5*(tanh(0.5*(r - 27.d0)) + 1.d0)
* determine splined coefficients at R=r
       ind=1
       do ilam=1,19
         call dcopy(47,vl(ind),1,vec,1)
         call dsplint(rr,vec,csplin(1,ilam),47,r,vvx)
         if (ilam.eq.1) then
* merge with asymptotic form for V0
           vvx = (1.d0 - switch_lr)*vvx 
     +       + switch_lr*c6sum/(r**6)
         else
* kill anisotropic terms at large R
           vvx = (1.d0 - switch_lr)*vvx
         endif
         v(ilam)=vvx
         call dcopy(47,vl(ind),1,vec,1)
         ind = ind + 47
       enddo
       call dcopy(18,v(2),1,vvl(1),1)
* convert to hartree
       econv=1.d0/219474.6d0
       call dscal(18,econv,vvl,1)
* isotropic term 
       vv0 = v(1)*econv
*
       return
       end
*===========================eof===============================