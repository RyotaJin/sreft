      real (kind=dpsize), allocatable, save :: offt(:)
    contains

      subroutine initialize(id, meanx, meany, coun, a, b, c, bm, serial, covt, covy)
        real (kind=dpsize), intent(in) :: meanx(*), meany(*), coun(*), a(*), b(*), c(*), covt, covy(*)
        real (kind=dpsize) :: sum, wsum, ww, ma, mi, offsett, m, sigma
        integer (kind=isize), intent(in) :: id, bm, serial
        integer (kind=isize) :: i

        sum = 0.0d0
        wsum = 0.0d0
        do i = 1, numbm
          if (coun(i)==0 .or. b(i)==0.0d0) cycle

          sigma = omega(netas_+i, netas_+i)
          ma = plotmax - meanx(i)
          mi = plotmin + meanx(i)
          if (bm==1) then
!main biomarker
            if (c(i)==0.0d0) then
              m = (meany(i)-a(i))/(b(i)+covy(i))/covt
            else if ((meany(i)-a(i))/(b(i)/c(i)+covy(i))+1.0<0.0d0) then
              if (c(i)>0.0d0) then
                m = mi
              else
                m = ma
              end if
            else
              m = log((meany(i)-a(i))/(b(i)/c(i)+covy(i))+1.0)/c(i)/covt
            end if
          else
!other biomarker
            if (c(i)==0.0d0) then
              m = (meany(i)-a(i)-covy(i))/b(i)/covt
            else if ((meany(i)-a(i)-covy(i))*c(i)/b(i)+1.0<0.0d0) then
              if (c(i)>0.0d0) then
                m = mi
              else
                m = ma
              end if
            else
              m = log((meany(i)-a(i)-covy(i))*c(i)/b(i)+1.0)/c(i)/covt
            end if
          end if
          if (m<mi) then
            m = mi
          else if (m>ma) then
            m = ma
          end if

          ww = sigma/coun(i)
          ww = ww/(b(i)*covt*exp(c(i)*covt*m))**2
          ww = 1/sqrt(ww)
          m = m - meanx(i)
          sum = sum + ww*m
          wsum = wsum + ww
        end do
        if (wsum>0.0d0) then
          offsett = sum/wsum
        else
          offsett = 0.0d0
        end if
        offt(serial) = offsett
        realid(serial) = id
        return
      end subroutine initialize

      subroutine eval(a, b, c, f, g, time, bm, serial, covt, covy)
        use :: prdims, only: gprd
        real (kind=dpsize), intent(in) :: a(*), b(*), c(*), covt, covy(*)
        real (kind=dpsize), intent(out) :: f, g(gprd, 1)
        real (kind=dpsize), intent(inout) :: time
        integer (kind=isize), intent(in) :: bm, serial

        time = time + offt(serial)
!main biomarker
        if (bm==1) then
          if (c(bm)==0) then
            f = a(bm) + (covy(bm)+b(bm))*covt*time
            g = 0.0d0
            g(bm, 1) = 1.0d0
            g(bm+numbm, 1) = covt*time
          else
            f = a(bm) + (b(bm)/c(bm)-covy(bm))*(exp(c(bm)*covt*time)-1.0d0)
            g = 0.0d0
            g(bm, 1) = 1.0d0
            g(bm+numbm, 1) = (exp(c(bm)*covt*time)-1.0d0)/c(bm)
            g(bm+2*numbm, 1) = (b(bm)/c(bm)-covy(bm))*covt*time*exp(c(bm)*covt*time) - &
              b(bm)/c(bm)**2*(exp(c(bm)*covt*time)-1.0d0)
          end if
        else
!other biomarker
          if (c(bm)==0) then
            f = a(bm) + covy(bm) + b(bm)*covt*time
            g = 0.0d0
            g(bm, 1) = 1.0d0
            g(bm+numbm, 1) = covt*time
          else
            f = a(bm) + covy(bm) + b(bm)/c(bm)*(exp(c(bm)*covt*time)-1.0d0)
            g = 0.0d0
            g(bm, 1) = 1.0d0
            g(bm+numbm, 1) = (exp(c(bm)*covt*time)-1.0d0)/c(bm)
            g(bm+2*numbm, 1) = b(bm)/c(bm)*covt*time*exp(c(bm)*covt*time) - b(bm)/c(bm)**2*(exp(c(bm)*covt*time)-1.0d0)
          end if
        end if
        return
      end subroutine eval

      subroutine makecsv
        integer (kind=isize) :: i

        open (13, file='../../offsetT.csv')
        write (13, '(a2, a1, a7)') 'ID', ',', 'offsetT'
        do i = 1, numsj
          write (13, '(i5, a1, e10.4)') realid(i), ',', offt(i)
        end do
        close (13)
        return
      end subroutine makecsv

    end module sreft

    subroutine pred(icall, newind, theta, datrec, indxs, f, g, h)
      use :: sreft
      use :: prdims, only: gprd, hprd
      use :: nmprd_real, only: eta, eps
      use :: nmprd_int, only: newl2, iquit
      implicit none
      real (kind=dpsize), intent(in) :: theta(*), datrec(*)
      real (kind=dpsize), intent(out) :: f, g(gprd, 1), h(hprd, 1)
      real (kind=dpsize) :: covt, time
      real (kind=dpsize), allocatable :: a(:), b(:), c(:), meanx(:), meany(:), coun(:), covy(:)
      integer (kind=isize), intent(in) :: icall, newind, indxs(*)
      integer (kind=isize) :: id, bm, serial

      if (icall==4) then
        if (newind/=2) then
          call simeta(eta)
          if (iquit==1) return
        end if
        if (newl2==1) then
          call simeps(eps)
          if (iquit==1) return
        end if
      else
        if (newind/=2) then
          call geteta(eta)
          if (iquit==1) return
          eps = 0.0d0
        end if
      end if

      if (icall==0 .and. newind==0) then
        numbm = nepss_
        numsj = nindr
        if (.not. allocated(offt)) then
          allocate (offt(numsj), realid(numsj))
        end if
      end if

      allocate (a(numbm), b(numbm), c(numbm), meanx(numbm), meany(numbm), coun(numbm), covy(numbm))
