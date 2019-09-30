      subroutine syminv (a,lda,n,ierr)
      implicit double precision (a-h,o-z)
c
c     -----------------------------------------------------------------
c     This subroutine uses LAPACK DGETRF and DGETRI
c     to invert a real symmetric indefinite matrix.
c     latest revision:  18-feb-2008 by millard alexander (replace calls
c      to "sy" with "ge"
c         lapack routines
c     -----------------------------------------------------------------
c
      dimension a(lda,n)
      dimension ipiv(4*n),work(256*n)

c
      lwork = 256*n
*  form full matrix
      do j = 2,n
         do i = 1,j-1
            a(i,j) = a(j,i)
         enddo
      enddo


*      write(6,444) lda,n
*444   format ('calling dgetrf:  lda=',i8,'  n=',i8)


      call dgetrf (n,n,a,lda,ipiv,ierr)
      if (ierr .ne. 0) stop 'error in dgetrf'
      call dgetri(n,a,lda,ipiv,work,lwork,ierr)
      if (ierr .ne. 0) stop 'error in dgetri'
      return
      end

