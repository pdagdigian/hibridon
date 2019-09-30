**************************************************************************
*                                                                        *
*                    vector routines library supplement             *
*                                                                        *
**************************************************************************
*                        routines included                               *
*                                                                        *
*  1. matmov   puts a nr x nc matrix a into nr x nc matrix b             *
*  2. maxmgv   finds largest value in a vector                           *
*  3. vsmul    multiplies the elements of a vector by a scalar           *
*  4. vmul     multiplies the elements of two vectors                    *
*  4a.dsum     sum of elements of a vector
*  5. blas     basis linear algebra routines from lapack
*     lsame,ilaenv,xerbla
*  6. blas extensions from lapack
*     dlaev2, dlasyf
*     isamax, saxpy, scopy, sdot, sscal, sswap, srot (grot)
*  7. fzero (vector zero) (no longer part of code)
*  8. blas extensions from ibm essl
*     idamin, idmin
*                                                                        *
**************************************************************************
      subroutine matmov (a, b, nr, nc, na, nb)
*  subroutine to put nr x nc matrix a into nr x nc matrix b
*  author:  millard alexander
*  current revision date: 24-sept-87
c ------------------------------------------------------------------
*  variables in call list
*    a,b:     input matrices, stored as one-dimensional arrays
*    nr:      actual row dimension of matrix a
*    nc:      actual column dimension of matrix a
*    na:      maximum row dimension of matrix a
*    nb:      maximum row dimension of matrix b
c ------------------------------------------------------------------
*  the two matrices are treated as vectors here, with column-by-column
*  ordering assumed
*  the coding uses the blas routine scopy
c ------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer ia, ib, j, na, nb, nr, nc
      dimension a(1), b(1)
      ia = 0
      ib = 0
*  ia and ib point to one position before the top of the jth
*  column of each matrix
      do 20  j = 1, nc
        call dcopy (nr, a(ia+1), 1, b(ib+1), 1)
        ia = ia + na
        ib = ib + nb
20    continue
      return
      end
*  -----------------------------------------------------------------------
      subroutine maxmgv (a, na, c, nc, n)
*  subroutine to scan a vector for its maximum magnitude (absolute value)
*  element
*  current revision date: 24-sept-87
*  -----------------------------------------------------------------------
*  variables in call list:
*  a:   floating point input vector
*  na:  integer element step for a
*  c:   floating point output scalar: on return contains value of
*       maximum magnitude (absolute value) element
*  nc:  integer index of maximum magnitude element
*  n:   integer element count
*  subroutines called:
*  isamax: blas routine to find index of maximum magnitude (absolute value)
*          element
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a(1)
      nc = ( idamax (n, a, na) - 1) * na + 1
      c = abs( a(nc) )
      return
      end
*  -----------------------------------------------------------------------
      subroutine vsmul (a, na, b, c, nc, n)
*  subroutine to multiply the elements of a vector by a scalar
*  current revision date: 24-sept-87
*  -----------------------------------------------------------------------
*  variables in call list:
*  a:   floating point input vector
*  na:  integer element step for a
*  b:   floating point input scalar
*  c:   floating point output vector
*  nc:  integer element step for c
*  n:   integer element count
*  c(m) = a(m) * b for m=1 to n
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer n, na, nc
      dimension a(1), c(1)
*  first copy vector a into vector c
*  then multiply by scalar constant
      call dcopy (n, a, na, c, nc)
      call dscal (n, b, c, nc)
      return
      end
*  -----------------------------------------------------------------------
      subroutine vmul (a, na, b, nb, c, nc, n)
*  subroutine to multiply the elements of two vectors
*  current revision date: 23-sept-87
*  -----------------------------------------------------------------------
*  variables in call list:
*  a:   floating point input vector
*  na:  integer element step for a
*  b:   floating point input vector
*  nb:  integer element step for b
*  c:   floating point output vector
*  nc:  integer element step for c
*  n:   integer element count
*  c(m) = a(m) * b(m) for m=1 to n
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer i, inda, indb, indc, n, na, nb, nc
      dimension a(1), b(1), c(1)
      inda = 1
      indb = 1
      indc = 1
      do 4  i = 1, n
        c(indc) = b(indb) * a(inda)
        inda = inda + na
        indb = indb + nb
        indc = indc + nc
4     continue
      return
      end
c JK commented out vadd subrooutine
c*  -----------------------------------------------------------------------
c      subroutine vadd (ic,a, na, b, nb, n)
c*  subroutine to add or subtract the elements of two vectors
c*  current revision date: 6-dec-1991
c*  -----------------------------------------------------------------------
c*  variables in call list:
c*  ic:  factor
c*  a:   floating point input vector
c*  na:  integer element step for a
c*  b:   floating point input vector
c*  nb:  integer element step for b
c*  n:   integer element count
c*  a(m) = a(m) + ic*b(m) for m=1 to n
c*  -----------------------------------------------------------------------
c      implicit double precision (a-h,o-z)
c      integer i,ic, inda, indb, n, na, nb
c      dimension a(1), b(1)
c      inda = 1
c      indb = 1
c      if (ic .gt. 0) then
c        do 4  i = 1, n
c          a(inda) = a(inda) + b(indb)
c          inda = inda + na
c          indb = indb + nb
c4       continue
c      else
c        do 5  i = 1, n
c          a(inda) = a(inda) - b(indb)
c          inda = inda + na
c          indb = indb + nb
c5       continue
c      endif
c      return
c      end
      double precision function dsum(n,dx,incx)
c
c     returns sum of double precision dx
c     dasum = sum from 0 to n-1 of dx(1+i*incx)
c     adapted from blas dasum by mha  4-apr-1996
c
      double precision dx(1)
      dsum = 0.d0
      if(n.le.0)return
      if(incx.eq.1)goto 20
c
c        code for increments not equal to 1
c
      ns = n*incx
          do 10 i=1,ns,incx
          dsum = dsum + dx(i)
   10     continue
      return
c
c        code for increments equal to 1.
c
c
c        clean-up loop so remaining vector length is a multiple of 6.
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         dsum = dsum + dx(i)
   30 continue
      if( n .lt. 6 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,6
         dsum = dsum + dx(i) + dx(i+1) + dx(i+2)
     1   + dx(i+3) + dx(i+4) + dx(i+5)
   50 continue
      return
      end
      subroutine dset(n,da,dx,incx)
c
c     sets a vector equal to a constant
c     uses unrolled loops for increment equal to one.
c     modified by mha from linpack dscal, written originally by
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      double precision dx(1), da
      integer i,incx,ix,m,mp1,n

c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
        dx(ix) = da
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da
        dx(i + 1) = da
        dx(i + 2) = da
        dx(i + 3) = da
        dx(i + 4) = da
   50 continue
      return
      end
      integer function idamin(n, sx, incx)
*
*     find smallest index of minimum magnitude of double precision s
*     isamax =  first i,  i = 1 to n,  to minimize  abs(sx(1-incx+i*incx)
*
      double precision sx(1), smin, xmag
      idamin = 0
      if (n .le. 0) return
      idamin = 1
      if (n .le. 1) return
      if (incx .eq. 1) go to 20
*
*        code for increments not equal to 1.
*
      smin = abs(sx(1))
      ns = n * incx
      ii = 1
          do 10 i = 1, ns, incx
          xmag = abs(sx(i))
          if (xmag .ge. smin) go to 5
          idamin = ii
          smin = xmag
    5     ii = ii + 1
   10     continue
      return
*
*        code for increments equal to 1.
*
   20 smin = abs(sx(1))
      do 30 i = 2, n
         xmag = abs(sx(i))
         if (xmag .ge. smin) go to 30
         idamin = i
         smin = xmag
   30 continue
      return
      end
      integer function idmin(n, sx, incx)
*
*     find smallest index of minimum element in double precision s
*     isamax =  first i,  i = 1 to n,  to minimize  sx(1-incx+i*incx)
*
      double precision sx(1), smin, xmag
      idmin = 0
      if (n .le. 0) return
      idmin = 1
      if (n .le. 1) return
      if (incx .eq. 1) go to 20
*
*        code for increments not equal to 1.
*
      smin = sx(1)
      ns = n * incx
      ii = 1
          do 10 i = 1, ns, incx
          xmag = sx(i)
          if (xmag .ge. smin) go to 5
          idmin = ii
          smin = xmag
    5     ii = ii + 1
   10     continue
      return
*
*        code for increments equal to 1.
*
   20 smin = abs(sx(1))
      do 30 i = 2, n
         xmag = sx(i)
         if (xmag .ge. smin) go to 30
         idmin = i
         smin = xmag
   30 continue
      return
      end
