* System:  CN(B)+Ar, scaled ab initio MRCI+Q PES's
* Reference:
* M. H. Alexander, Xin Yang, P J. Dagdigian, A Berning,
* and H-J Werner J. Chem. Phys. 112, 781 (2000)

      include "common/syusr"
      include "common/bausr"
      include "common/ground"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(6)
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,10(1pe16.8))
      goto 1
99    do i=1,200
         r=(i-1)*.1+4.5d0
         call pot(vv0,r)
         write (6,101) r,vv0
101      format(f6.2,1pe16.8)
      enddo
99    end
      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='ALEXANDER CN(B)-Ar avqz+bondf scaled MRCI(D)'
      lammin(1)=1
      lammax(1)=5
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients in the
*  Ar-CN(X,B) potential of Alexander
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0        contains isotropic term (d00 term in vsum)
*  variable in common block /covvl/
*    vvl:     vector of length 6 to store r-dependence of each term
*             in potential expansion
*    vvl(1-6) expansion coefficients in dl0 (l=1:5) of vsum

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  11-Aug-2006
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      parameter (kmax=19,kang=10)
*     dimension rr(kmax),vx(kmax*kang),vb(kmax*kang),z(kmax),cx(kmax),
*    :    coefx(kmax,kang),coefb(kmax,kang),xa(kmax)
      dimension rr(19),vx(190),vb(190),z(19),cx(19),
     :    coefb(19,10),xa(19),coefc(19,10),coefd(19,10),
     :    cb(19),cc(19),cd(19)
      dimension d0(100),aa(100)
      dimension kpvt(10),qraux(10),work(10),rsd(10),xsum(10),vang(10)

      common /covvl/ vvl(6)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 10 angles and for l=0:9
* angles are 0 15 30 60 90 120 135 150 165 180
* NOTE:  theta = 0 corresponds to CNAr
      data d0/
     : 1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,9.6592583d-1,
     : 8.6602540d-1,5d-1,0d0,-5d-1,-7.0710678d-1,-8.6602540d-1,
     : -9.6592583d-1,-1d0,1d0,8.9951905d-1,6.25d-1,-1.25d-1,-5d-1,
     : -1.25d-1,2.50d-1,6.25d-1,8.9951905d-1,1d0,1d0,8.0416392d-1,
     : 3.2475953d-1,-4.375d-1,0d0,4.375d-1,1.7677670d-1,-3.2475953d-1,
     : -8.0416392d-1,-1d0,1d0,6.8469544d-1,2.3437500d-2,-2.8906250d-1,
     : 3.75d-1,-2.8906250d-1,-4.0625000d-1,2.3437500d-2,6.8469544d-1,
     : 1d0,1d0,5.4712587d-1,-2.2327217d-1,8.9843750d-2,0d0,
     : -8.9843750d-2,3.7565048d-1,2.2327217d-1,-5.4712587d-1,-1d0,1d0,
     : 3.9830599d-1,-3.7402344d-1,3.2324219d-1,-3.125d-1,3.2324219d-1,
     : -1.4843750d-1,-3.7402344d-1,3.9830599d-1,1d0,1d0,2.4554105d-1,
     : -4.1017805d-1,2.2314453d-1,0d0,-2.2314453d-1,-1.2705825d-1,
     : 4.1017805d-1,-2.4554105d-1,-1d0,1d0,9.6184327d-2,-3.3877563d-1,
     : -7.3638916d-2,2.7343750d-1,-7.3638916d-2,2.9833984d-1,
     : -3.3877563d-1,9.6184327d-2,1d0,1d0,-4.2767847d-2,-1.8957520d-1,
     : -2.6789856d-1,0d0,2.6789856d-1,-2.8553579d-1,1.8957520d-1,
     : 4.2767847d-2,-1d0/

* coefficients for expansion of potential, ordered by angles, each column corresponds
* to theta from 0 (1st column) to 180 (last column)
* R=[4.75 5 5.25, 5.5, 5.75, 6.00, 6.25, 6.5, 6.75, 7.00, 
*    7.25, 7.5, 7.75, 8.00, 8.5, 9.00, 10, 12, 40]
* independent variable is x=exp(-0.2*r)
      
      data xa /
     ; 3.3546263d-4, 9.0717953d-2, 1.3533528d-1, 1.6529889d-1,
     : 1.8268352d-1, 2.0189652d-1, 2.1224797d-1, 2.2313016d-1,
     : 2.3457029d-1, 2.4659696d-1, 2.5924026d-1, 2.7253179d-1,
     : 2.8650480d-1, 3.0119421d-1, 3.1663677d-1, 3.3287108d-1,
     : 3.4993775d-1, 3.6787944d-1, 3.8674102d-1/
* scaled (s=1.15) avtz+bondf mrci pes
      data vx  /
* theta = 0
     : 0d0,-8.9490000d0,-3.0847500d1,-5.9866500d1,-8.0578000d1,
     : -9.9553500d1,-1.0279750d2,-9.5049500d1,-6.6614000d1,
     : -1.3730000d0,1.2668050d2,3.5868550d2,7.5851750d2,1.4226275d3,
     : 2.4915415d3,4.1575090d3,6.6510960d3,1.0170417d4,1.4772608d4,
* theta = 15 
     : 0d0,-8.7307111d0,-3.0028054d1,-5.8431889d1,-7.9068250d1,
     : -9.8994049d1,-1.0370252d2,-9.8763281d1,-7.5365207d1,
     : -1.8884433d1,9.4386690d1,3.0209372d2,6.6276097d2,1.2651663d3,
     : 2.2393950d3,3.7653662d3,6.0653024d3,9.3466736d3,1.3690754d4,
* theta = 30
     : 0d0,-8.1704595d0,-2.7843858d1,-5.4432554d1,-7.4499588d1,
     : -9.6069126d1,-1.0386371d2,-1.0509899d2,-9.3105773d1,
     : -5.6664397d1,2.2560261d1,1.7394471d2,4.4356864d2,9.0213191d2,
     : 1.6550241d3,2.8526494d3,4.6938825d3,7.4005895d3,1.1127761d4,
* theta = 60
     : 0d0,-6.7582583d0,-2.2375086d1,-4.3825025d1,-6.1196744d1,
     : -8.3177551d1,-9.4802458d1,-1.0506507d2,-1.1111640d2,
     : -1.0790470d2,-8.6828090d1,-3.3649231d1,7.4674186d1,2.7474019d2,
     : 6.2367681d2,1.2088303d3,2.1602165d3,3.6632541d3,5.9596207d3,
* theta = 90
     : 0d0,-6.2655618d0,-2.0428734d1,-3.9913972d1,-5.5904009d1,
     : -7.6884330d1,-8.8704276d1,-1.0032248d2,-1.0987537d2,
     : -1.1398469d2,-1.0677309d2,-7.8394425d1,-1.2767171d1,1.1587958d2,
     : 3.4810675d2,7.4666426d2,1.4067300d3,2.4686902d3,4.1295285d3,
* theta = 120
     : 0d0,-7.3094087d0,-2.4592067d1,-4.8447897d1,-6.7651981d1,
     : -9.1737814d1,-1.0434093d2,-1.1537914d2,-1.2184997d2,
     : -1.1861765d2,-9.7228915d1,-4.4220255d1,6.1304170d1,2.5089910d2,
     : 5.7080166d2,1.0863298d3,1.8856654d3,3.0807144d3,4.8043621d3,
* theta = 135
     : 0d0,-8.2157844d0,-2.8294715d1,-5.5907765d1,-7.7466354d1,
     : -1.0272744d2,-1.1437009d2,-1.2221442d2,-1.2172398d2,
     : -1.0545854d2,-6.1645927d1,2.7780285d1,1.8969537d2,4.6284044d2,
     : 9.0110532d2,1.5766819d3,2.5833897d3,4.0444124d3,6.1358359d3,
* theta = 150
     : 0d0,-9.1371031d0,-3.2040834d1,-6.3067418d1,-8.6059616d1,
     : -1.1001044d2,-1.1824706d2,-1.1908429d2,-1.0577337d2,
     : -6.7653463d1,1.1616861d1,1.5636961d2,4.0176173d2,7.9716173d2,
     : 1.4103297d3,2.3338754d3,3.6993816d3,5.7091726d3,8.6944001d3,
* theta = 165
     : 0d0,-9.8237914d0,-3.4854119d1,-6.8036402d1,-9.1237719d1,
     : -1.1198927d2,-1.1554997d2,-1.0786398d2,-7.9983066d1,
     : -1.8067598d1,9.8754976d1,3.0113345d2,6.3304398d2,1.1565704d3,
     : 1.9595241d3,3.1705255d3,4.9905328d3,7.7541384d3,1.2014814d4,
* theta = 180
     : 0d0,-1.0038000d1,-3.5811000d1,-6.9530000d1,-9.2374500d1,
     : -1.1092700d2,-1.1183400d2,-9.9566000d1,-6.4122500d1,9.9380000d0,
     : 1.4571800d2,3.7710850d2,7.5280300d2,1.3421090d3,2.2448530d3,
     : 3.6109825d3,5.6814795d3,8.8651305d3,1.3848195d4/
      data vb / 
* theta = 0
     : 0d0,-8.6895000d0,-3.0392500d1,-6.1504000d1,-8.7762500d1,
     : -1.2346050d2,-1.4460150d2,-1.6695450d2,-1.8871750d2,
     : -2.0673850d2,-2.1574750d2,-2.0721500d2,-1.6670100d2,
     : -6.7232000d1,1.4909800d2,6.1905600d2,1.6729000d3,4.0685820d3,
     : 9.3070940d3,
* theta = 15
     : 0d0,-8.5573492d0,-2.9682575d1,-5.9943618d1,-8.5504299d1,
     : -1.2024280d2,-1.4079779d2,-1.6248372d2,-1.8351985d2,
     : -2.0071280d2,-2.0873925d2,-1.9887595d2,-1.5641290d2,
     : -5.4185649d1,1.6469230d2,6.3256707d2,1.6714009d3,4.0518402d3,
     : 9.3989254d3,
* theta = 30
     : 0d0,-8.0578606d0,-2.7516984d1,-5.5157836d1,-7.8293700d1,
     : -1.0930540d2,-1.2728945d2,-1.4572427d2,-1.6253889d2,
     : -1.7400933d2,-1.7387387d2,-1.5180563d2,-9.0747533d1,3.8655613d1,
     : 2.9092227d2,7.7703935d2,1.7382389d3,3.6943366d3,7.6130126d3,
* theta = 60
     : 0d0,-6.7349453d0,-2.2232082d1,-4.3563409d1,-6.0978014d1,
     : -8.3474325d1,-9.5852282d1,-1.0758459d2,-1.1639392d2,
     : -1.1821206d2,-1.0611616d2,-6.8658578d1,1.2705631d1,1.6767939d2,
     : 4.4365152d2,9.1637881d2,1.7078736d3,3.0158930d3,5.1577175d3,
* theta = 90
     : 0d0,-6.4868194d0,-2.1002603d1,-4.0529606d1,-5.6073696d1,
     : -7.5575202d1,-8.5881481d1,-9.5109042d1,-1.0098317d2,
     : -9.9547434d1,-8.4140397d1,-4.3876899d1,3.8650883d1,1.9085052d2,
     : 4.5522349d2,8.9679906d2,1.6137835d3,2.7524224d3,4.5274042d3,
* theta = 120
     : 0d0,-8.1102836d0,-2.7068348d1,-5.1754276d1,-7.0012952d1,
     : -8.9900264d1,-9.7929910d1,-1.0160473d2,-9.6622361d1,
     : -7.6047789d1,-2.8965539d1,6.1440664d1,2.2068564d2,4.8696959d2,
     : 9.1687683d2,1.5936334d3,2.6392983d3,4.2319956d3,6.6278238d3,
* theta = 135
     : 0d0,-9.4625855d0,-3.2299011d1,-6.0838989d1,-7.9845632d1,
     : -9.5503608d1,-9.6942558d1,-8.8417505d1,-6.2369104d1,
     : -7.0108295d0,9.5676277d1,2.7281663d2,5.6469288d2,1.0306052d3,
     : 1.7574160d3,2.8721755d3,4.5597913d3,7.0834746d3,1.0796382d4,
* theta = 150
     : 0d0,-1.0889666d1,-3.7735031d1,-6.9395606d1,-8.7129926d1,
     : -9.2943529d1,-8.2164857d1,-5.2469275d1,8.8316453d0,1.2123639d2,
     : 3.1430908d2,6.3240933d2,1.1415689d3,1.9397926d3,3.1725488d3,
     : 5.0550761d3,7.8990088d3,1.2128540d4,1.8249432d4,
* theta = 165
     : 0d0,-1.2005329d1,-4.1925428d1,-7.5284899d1,-9.0490263d1,
     : -8.4615006d1,-6.0360856d1,-7.8646794d0,9.0920631d1,2.6368329d2,
     : 5.5254313d2,1.0211675d3,1.7656024d3,2.9311013d3,4.7381293d3,
     : 7.5191525d3,1.1759651d4,1.8116026d4,2.7337666d4,
* theta = 180
     : 0d0,-1.2394500d1,-4.3419500d1,-7.7162000d1,-9.0945500d1,
     : -7.9455000d1,-4.8936500d1,1.4087000d1,1.3012900d2,3.3071350d2,
     : 6.6389800d2,1.2026065d3,2.0575065d3,3.3973325d3,5.4805665d3,
     : 8.7005885d3,1.3638460d4,2.1100614d4,3.2092486d4/
      data rr /
     : 40d0, 12d0, 10d0, 9d0, 8.5d0, 8d0, 7.75d0, 7.5d0, 7.25d0, 7d0,
     : 6.75d0, 6.5d0, 6.25d0, 6d0, 5.75d0, 5.5d0, 5.25d0, 5d0, 4.75d0 /
* spline fit to obtain coefficients

      data ispline /0/
      nspline=kmax 
      lmax=6
      nang=kang
      if (ispline.eq.0) then
        ind=0
        do ith=1,nang
            call dcopy(nspline,vb(ind+1),1,z,1)
            call spline(nspline,xa,z,cb,cc,cd)
            call dcopy(nspline,cb,1,coefb(1,ith),1)
            call dcopy(nspline,cc,1,coefc(1,ith),1)
            call dcopy(nspline,cd,1,coefd(1,ith),1)
            ind=ind+nspline
        enddo
        ispline=1
      endif

* determine potentials at angles
      x=exp(-0.2d0*r)
      ind=0
      do ith=1,nang
         if (r.lt.rr(nspline)) then
            vspline=vb(ind+nspline)
         else
            call dcopy(nspline,vb(ind+1),1,z,1)
            call dcopy(nspline,coefb(1,ith),1,cb,1)
            call dcopy(nspline,coefc(1,ith),1,cc,1)
            call dcopy(nspline,coefd(1,ith),1,cd,1)
            vspline=seval(nspline,x,xa,z,cb,cc,cd)
         endif
         ind=ind+nspline
         vang(ith)=vspline
         if (ith.eq.1) then
            vv0=vspline
         else 
            vvl(ith-1)=vspline
         endif
      enddo
* solve simultaneous equations for solutions
* lmax is maximum order of legendre polynomials included in expansion (l=0:lmax)
      tol=1.e-10
      call dcopy(10*(lmax+1),d0,1,aa,1)
      call dqrank(aa,10,10,lmax+1,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,10,10,lmax+1,kr,vang,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6
      call dscal(7,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(6,xsum(2),1,vvl,1)
      end