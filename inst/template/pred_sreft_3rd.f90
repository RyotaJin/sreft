      if (newind<=1) call initialize(id, meanx, meany, coun, a, b, c, bm, serial, covt, covy)

      call eval(a, b, c, f, g, time, bm, serial, covt, covy)

      f = f + eps(bm)
      h = 0.0d0
      h(bm, 1) = 1.0d0

      if (icall==3) call makecsv

      deallocate (a, b, c, meanx, meany, coun, covy)

      return
    end subroutine pred
