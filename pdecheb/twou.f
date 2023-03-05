      double precision twou, temp
      integer :: ia(4)
      equivalence(ia(1),i1),(ia(2),i2),(ia(3),i3),(ia(4),i4)

      TWOU = 0.1D0
   40 TEMP = 1.0D0 + TWOU
      IF (1.0D0.EQ.TEMP) THEN
         TWOU = TWOU*2.0D0
      ELSE
         TWOU = TWOU*0.5D0
         GO TO 40
      END IF
      print *, TWOU, epsilon(TWOU)

      ia = [(i,i=1,4)]
      print *, i1, i2 ,i3, i4
      end