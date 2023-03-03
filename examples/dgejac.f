C  jac -- if you have set info(5)=0, you can ignore this parameter
C               by treating it as a dummy argument. otherwise, you must
C               provide a subroutine of the form
C               jac(t,y,yprime,pd,cj,rpar,ipar)
C               to define the matrix of partial derivatives
C               pd=dg/dy+cj*dg/dyprime
C               cj is a scalar which is input to jac.
C               for the given values of t,y,yprime, the
C               subroutine must evaluate the non-zero partial
C               derivatives for each equation and each solution
C               compowent, and store these values in the
C               matrix pd. the elements of pd are set to zero
C               before each call to jac so only non-zero elements
C               need to be defined.
C
C               subroutine jac must not alter t,y,(*),yprime(*),or cj.
C               you must declare the name jac in an
C               external statement in your program that calls
C               ddassl. you must dimension y, yprime and pd
C               in jac.
C
C               the way you must store the elements into the pd matrix
C               depends on the structure of the matrix which you
C               indicated by info(6).
C               *** info(6)=0 -- full (dense) matrix ***
C                   when you evaluate the (non-zero) partial derivative
C                   of equation i with respect to variable j, you must
C               store it in pd according to
C                   pd(i,j) = * dg(i)/dy(j)+cj*dg(i)/dyprime(j)*
C               *** info(6)=1 -- banded jacobian with ml lower and mu
C                   upper diagonal bands (refer to info(6) description o
C                   ml and mu) ***
C                   when you evaluate the (non-zero) partial derivative
C                   of equation i with respect to variable j, you must
C                   store it in pd according to
C                   irow = i - j + ml + mu + 1
C                   pd(irow,j) = *dg(i)/dy(j)+cj*dg(i)/dyprime(j)*
C               rpar and ipar are real and integer parameter arrays whic
C               you can use for communication between your calling
C               program and your jacobian subroutine jac. they are not
C               altered by ddassl. if you do not need rpar or ipar, igno
C               these parameters by treating them as dummy arguments. if
C               you do choose to use them, dimension them in your callin
C               program and in jac as arrays of appropriate length.

      subroutine dgejac(t,y,yprime,pd,cj,rpar,ipar)
c     intent(in)
      double precision t, y(*), yprime(*), cj
      double precision rpar(*)
      integer ipar(*)
c     intent(inout)
      double precision pd(*)
c
c     dummy routine
c
      end
