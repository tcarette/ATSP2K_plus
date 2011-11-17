     program pppp
     integer inptr(2000), ibe, n_read
     double precision coeff(2000)


       n_read = 1100
      open(unit=39,file='/tmp/georgio/c.lst.00',status='unknown',
     :     form='unformatted');


            read(39) n, (coeff(j),j = ibe+1,ibe+n_read),
     :                    (inptr(j), j = ibe+1,ibe+n_read)
            print*, n, (coeff(j),j = ibe+1,ibe+n_read),
     :                    (inptr(j), j = ibe+1,ibe+n_read)

      end program pppp
