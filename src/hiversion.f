      subroutine version(iunit)
      character*100 profile
      character*140 build
      common /bld_config/ build(5)
      integer(4) i
* to output version number of hibridon code
* current revision date:  15-aug-2009
      write (iunit, 10)
10    format
     : (/' ---------------------------------------------------',
     :  '-----------------------',
     :  /,9x,
     :   '  HIBRIDON SCATTERING CODE V 4.4 BUILD DATE 09/25/19 17:47:19 EDT',
     : //'      AUTHORS: M. ALEXANDER, D. MANOLOPOULOS,',
     :   ' H.-J. WERNER, B. FOLLMEG,'/
     :   '               P. DAGDIGIAN',
     :  /' CONTRIBUTORS: D. LEMOINE, P. VOHRALIK,',
     :  ' G. COREY, R. JOHNSON, T. ORLIKOWSKI,',
     :  /'          A. BERNING, A. DEGLI-ESPOSTI,',
     :  ' C. RIST, B. POUILLY, J. KLOS, Q. MA,',/,
     :  '          G. VAN DER SANDEN, M. YANG, F. DE WEERD',
     :  ', S. GREGURICK, F. LIQUE',
     : /' ---------------------------------------------------',
     :  '-----------------------')
      write (iunit,15)
15    format(11x,'BUILD CONFIGURATION:')
      do i=1,5
          call linout(build(i),iunit)
      enddo
      write (iunit,20)
20    format (' ---------------------------------------------------',
     :  '-----------------------')
      i=system("hib_sysconfig")
      open(unit=8,file='sysprofile',access='sequential')
      read (8,25) profile
25    format(a100)
      write (iunit,30) profile
30    format(' CURRENT HARDWARE CONFIGURATION:',/,9x,a100,/,
     : ' ---------------------------------------------------',
     :  '-----------------------')
       close(8)
       i=system("rm -f sysprofile")
      return
      end
      subroutine acknow(iunit,ipos)
* to acknowledge authors of hibridon
* current revision date:  24-feb-2004
      logical ipos
      if (ipos) write (iunit, 10)
10    format
     : (/' --------------------------------------------------------',
     :  '----------------------------------------------------------',
     :  '-------------',
     :  /,' All publications resulting from use of the integrators',
     :  ' included in the Hibridon code must include',
     :   ' the following reference:',
     : //,' D. E. Manolopoulos, J. Chem. Phys. 85, 6425 (1986);',
     :  ' M. H. Alexander and D. E. Manolopoulos, J. Chem. Phys.',
     :  ' 80, 2044 (1987).',
     :  //,' All publications resulting from use of the',
     :  ' Hibridon package must include',
     :  ' the following reference:',
     : //,' HIBRIDON is a package of programs for the time-independent',
     :  ' quantum treatment',
     :  ' of inelastic collisions and photodissociation',
     :  /,' written by',
     :  ' M. H. Alexander, D. E.  Manolopoulos, H.-J. Werner, and',
     : ' B. Follmeg,',
     :  ' with contributions by',
     :  ' P. F. Vohralik, D. Lemoine,',
     :  /,' G. Corey, B. Johnson, T. Orlikowski, A. Berning,',
     :  ' A. Degli-Esposti, C. Rist, P. Dagdigian, B. Pouilly,',
     :  ' G. van der Sanden, M. Yang, F. de Weerd, S. Gregurick,',
     :  ' J. Klos, and F. Lique',
     : /,' --------------------------------------------------------',
     :  '-------------',
     :  '----------------------------------------------------------')
      if (.not.ipos) write (iunit, 20)
20    format
     : (/' --------------------------------------------------------',
     :  '---------------------',
     :  /,' All publications resulting from use of the integrators',
     :  ' included in the',/,' Hibridon code must include',
     :   ' the following reference:',
     : //,' D. E. Manolopoulos, J. Chem. Phys. 85, 6425 (1986);',
     :  ' M. H. Alexander and D. E.',/,
     :    ' Manolopoulos, J. Chem. Phys.',
     :  ' 80, 2044 (1987).',
     :  //,' All publications involving the determination of',
     :     'photodissociation cross sections',/,' must also include',
     :     ' the following reference:',
     :  //,' M. H. Alexander, Comput. Phys. Commun, 75, 87 (1993).',
     :  //,' All publications investigating flux redistribution must',
     :  ' include the following references:',
     :  //, ' M. H. Alexander, J. Chem. Phys. 95, 8931 (1991); D. E.',
     :      ' Manolopoulos and',/,' M. H. Alexander, J. Chem. Phys.',
     :      '  97, 2527 (1992).',
     :  //,' All publications resulting from use of the',
     :  ' Hibridon package must include',
     :  /,' the following reference:',
     : //,' HIBRIDON is a package of programs for the time-independent',
     :  ' quantum treatment',
     :  /,' of inelastic collisions and photodissociation',
     :  ' written by',
     :  ' M. H. Alexander,',/,' D. E.  Manolopoulos, H.-J. Werner, and',
     : ' B. Follmeg,',
     :  ' with contributions by',
     :  /,' P. F. Vohralik, D. Lemoine,',
     :  ' G. Corey, B. Johnson, T. Orlikowski, A. Berning,',
     :  /,' A. Degli-Esposti, C. Rist, P. Dagdigian, B. Pouilly,',
     :  ' G. van der Sanden, M. Yang, F. de Weerd, S. Gregurick,',
     :  ' J. Klos, and F. Lique',
     : /,' --------------------------------------------------------',
     :  '---------------------')
      return
      end
      subroutine linout(l,iunit)
      character*(*) l
      write (iunit,'(a)') l(1:lenstr(l))
      return
      end
      block data config
c variables in data config
      character*140 build
      common /bld_config/ build(5)
      data build(1) /
     :  'SYS=RedHat 7.6.1810 ( 3.10.0-862.9.1.el7.x86_64); MACH=unix-x86 unix-ifort'/
      data build(2) /
     :  'CC=; FC=ifort -O3 -save -i8'/
      data build(3) /
     :  ''/
      data build(4) /
     : ''/
      data build(5) /
     : ''/
cend
      end
