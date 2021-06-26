MODULE bessel0


 contains

function u(x,k,d)
        use nrtype
        implicit none
        real(dp) x,u,d,k1
        integer k
        
        k1=float(k)
        u=(x-d)**(k1-1.0)
        return
        end function

 function u1(x,k,d)
        use nrtype
        implicit none
        real(dp) x,u1,d,k1
        integer k
        
         k1=float(k)
        u1=(k1-1.0)*((x-d)**(k1-2.0))
        return
        end function
 
 
 
 
   function amat(x,fdata,i,j)
        use nrtype
        implicit none
        
        real(dp) x,aa,in1,in2,in3,mm,nn,pp,mnp,amat
        integer i,j
        real(dp) fdata(6)
        in1=fdata(1)
        in2=fdata(2)
        in3=fdata(3)
        mm=fdata(4)
        nn=fdata(5)
        pp=fdata(6)
        
        mnp=mm+nn+pp
        if((i.eq.1).and.(j.eq.1))then
           aa=-(in1+in2+in3)/x
        else if((i.eq.1).and.(j.eq.2))then
            aa=mm
        else if((i.eq.1).and.(j.eq.3))then    
            aa=nn
         else if((i.eq.1).and.(j.eq.4))then  
            aa=pp
            
        else if((i.eq.2).and.(j.eq.1))then    
            aa=-mm
        else if((i.eq.2).and.(j.eq.2))then    
            aa=(in1-in2-in3-1.0)/x
        else if((i.eq.2).and.(j.eq.5))then   
            aa=nn
        else if((i.eq.2).and.(j.eq.7))then
            aa=pp
            
        else if((i.eq.3).and.(j.eq.1))then    
            aa=-nn
        else if((i.eq.3).and.(j.eq.3))then   
            aa=(in2-in1-in3-1.0)/x
         else if((i.eq.3).and.(j.eq.5))then    
            aa=mm
         else if((i.eq.3).and.(j.eq.6))then 
            aa=pp
            
        else if((i.eq.4).and.(j.eq.1))then   
            aa=-pp
        else if((i.eq.4).and.(j.eq.4))then 
            aa=(in3-in1-in2-1.0)/x
        else if((i.eq.4).and.(j.eq.6))then     
            aa=nn
        else if((i.eq.4).and.(j.eq.7))then 
            aa=mm
            
        else if((i.eq.5).and.(j.eq.2))then   
            aa=-nn
        else if((i.eq.5).and.(j.eq.3))then 
            aa=-mm
        else if((i.eq.5).and.(j.eq.5))then     
            aa=(in1+in2-in3-2.0)/x
        else if((i.eq.5).and.(j.eq.8))then 
            aa=pp
            
        else if((i.eq.6).and.(j.eq.3))then   
            aa=-pp
        else if((i.eq.6).and.(j.eq.4))then 
            aa=-nn
        else if((i.eq.6).and.(j.eq.6))then     
            aa=(in2+in3-in1-2.0)/x
        else if((i.eq.6).and.(j.eq.8))then 
            aa=mm     
            
            
        else if((i.eq.7).and.(j.eq.2))then   
            aa=-pp
        else if((i.eq.7).and.(j.eq.4))then 
            aa=-mm
        else if((i.eq.7).and.(j.eq.7))then     
            aa=(in1+in3-in2-2.0)/x
        else if((i.eq.6).and.(j.eq.8))then 
            aa=nn 
            
            
        else if((i.eq.8).and.(j.eq.5))then   
            aa=-pp
        else if((i.eq.8).and.(j.eq.7))then 
            aa=-nn
        else if((i.eq.8).and.(j.eq.6))then     
            aa=-mm
        else if((i.eq.8).and.(j.eq.8))then 
            aa=(in1+in2+in3-3.0)/x   
        
        else if((i.gt.8).or.(j.gt.8))then 
            print *, "A matrix index out of range"
            
        else
           aa=0.0
        endif   
        amat=aa
        return
        end function
        
        
    function amat22(x,r,i,j)
        use nrtype
        implicit none
        real(dp) x,aa,r,amat22
        integer i,j
        
     
        if((i.eq.1).and.(j.eq.1))then
           aa=0.0
        else if((i.eq.1).and.(j.eq.2))then
            aa=-r
        else if((i.eq.2).and.(j.eq.1))then    
            aa=r
         else if((i.eq.2).and.(j.eq.2))then     
            aa=-1.0/x
        else
            print *, "A matrix index out of range"
        endif   
        amat22=aa
        return
        end function
 
  function ff(x)
        use nrtype
        implicit none
        real(dp) x,ff
        
        ff=1.0/(x**2+1.0)
        return
        end function
        
        
     function ffunc(i,x,ind)
        use nrtype
        implicit none
        real(dp) x,ffunc,ind
        integer i
        
        if(i.eq.1)then
        ffunc=x**ind
        else
        ffunc=0.0
        endif
        return
        end function    
        
        function wfunc22(i,x)
        use nrtype
        implicit none
        real(dp) x,wfunc22
        integer i
        
        if(i.eq.1)then
        wfunc22=j0(x)
        else if(i.eq.2)then
         wfunc22=j1(x)
        else
         print *, "w indx out of range"
        endif
        return
        end function 
        
    function wfunc(i,x,fdata)
        use nrtype
        implicit none
        real(dp) x,wfunc,in1,in2,in3,mm,nn,pp
        integer i
         real(dp) fdata(6)
        in1=fdata(1)
        in2=fdata(2)
        in3=fdata(3)
        mm=fdata(4)
        nn=fdata(5)
        pp=fdata(6)
        
        if(i.eq.1)then
        wfunc=beslj(in1,mm*x)*beslj(in2,nn*x)*beslj(in3,pp*x)
        else if(i.eq.2)then
         wfunc=beslj(in1-1.0,mm*x)*beslj(in2,nn*x)*beslj(in3,pp*x)
        else if(i.eq.3)then
         wfunc=beslj(in1,mm*x)*beslj(in2-1.0,nn*x)*beslj(in3,pp*x)
        else if(i.eq.4)then
         wfunc=beslj(in1,mm*x)*beslj(in2,nn*x)*beslj(in3-1.0,pp*x)
        else if(i.eq.5)then
         wfunc=beslj(in1-1.0,mm*x)*beslj(in2-1.0,nn*x)*beslj(in3,pp*x) 
        else if(i.eq.6)then
         wfunc=beslj(in1,mm*x)*beslj(in2-1.0,nn*x)*beslj(in3-1.0,pp*x)
        else if(i.eq.7)then
         wfunc=beslj(in1-1.0,mm*x)*beslj(in2,nn*x)*beslj(in3-1.0,pp*x) 
        else if(i.eq.8)then
         wfunc=beslj(in1-1.0,mm*x)*beslj(in2-1.0,nn*x)*beslj(in3-1.0,pp*x) 
        else
         print *, "w indx out of range"
        endif
        return
        end function 
        
            function wfuncY(i,x,fdata)
        use nrtype
        implicit none
        real(dp) x,wfuncY,in1,in2,in3,mm,nn,pp
        integer i
         real(dp) fdata(6)
        in1=fdata(1)
        in2=fdata(2)
        in3=fdata(3)
        mm=fdata(4)
        nn=fdata(5)
        pp=fdata(6)
        
        if(i.eq.1)then
        wfuncY=beslY(in1,mm*x)*beslj(in2,nn*x)*beslj(in3,pp*x)
        else if(i.eq.2)then
         wfuncY=beslY(in1-1.0,mm*x)*beslj(in2,nn*x)*beslj(in3,pp*x)
        else if(i.eq.3)then
         wfuncY=beslY(in1,mm*x)*beslj(in2-1.0,nn*x)*beslj(in3,pp*x)
        else if(i.eq.4)then
         wfuncY=beslY(in1,mm*x)*beslj(in2,nn*x)*beslj(in3-1.0,pp*x)
        else if(i.eq.5)then
         wfuncY=beslY(in1-1.0,mm*x)*beslj(in2-1.0,nn*x)*beslj(in3,pp*x) 
        else if(i.eq.6)then
         wfuncY=beslY(in1,mm*x)*beslj(in2-1.0,nn*x)*beslj(in3-1.0,pp*x)
        else if(i.eq.7)then
         wfuncY=beslY(in1-1.0,mm*x)*beslj(in2,nn*x)*beslj(in3-1.0,pp*x) 
        else if(i.eq.8)then
         wfuncY=beslY(in1-1.0,mm*x)*beslj(in2-1.0,nn*x)*beslj(in3-1.0,pp*x) 
        else
         print *, "w indx out of range"
        endif
        return
        end function 
        
        
   function gg(x)
        use nrtype
        implicit none
        real(dp) x,gg
        
        gg=0.0
        return
        end function
 
   function j0(x)
        use nrtype
        implicit none
        real(dp) x,j0
        
        j0=bessel_j0(x)
        return
        end function
    
    function j1(x)
        use nrtype
        implicit none
        real(dp) x,j1
        
        j1=bessel_j1(x)
        return
        end function

     
         function beslJ ( alpha, x )
  
  use nrtype
  implicit none
  integer ( kind = 4 ), parameter :: nb_max = 100000
  
  real(dp) alpha
  real(dp) alpha_frac
  real(dp) b(nb_max)
  real(dp) fx
  real(dp) fx3
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) ncalc
  real(dp) x,nu,beslJ
  
  if(x.le.1e4)then
     nb = int ( alpha ) + 1
   alpha_frac = alpha - real ( int ( alpha ), kind = 8 )
  call rjbesl ( x, alpha_frac, nb, b, ncalc )
     beslJ= b(nb) 
 else
    beslj=beslJ_asymp ( alpha, x )
  endif
  return
  end
  
      function beslY ( alpha, x )
  
  use nrtype
  implicit none
  integer ( kind = 4 ), parameter :: nb_max = 100000
  
  real(dp) alpha
  real(dp) alpha_frac
  real(dp) b(nb_max)
  real(dp) fx
  real(dp) fx3
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) ncalc
  real(dp) x,nu,beslY
  
   if(x.le.1e4)then
     nb = int ( alpha ) + 1
   alpha_frac = alpha - real ( int ( alpha ), kind = 8 )
  call rybesl ( x, alpha_frac, nb, b, ncalc )
     beslY= b(nb) 
  else
    besly=beslY_asymp ( alpha, x )
  endif
  
  return
  end
  
     function beslJ_asymp ( nu, x )
     
     use nrtype
     implicit none
     real(dp) beslJ_asymp, nu, x
     
     beslJ_asymp=(Sqrt(2.0/Pi)*Cos(((1.0 + 2.0*nu)*Pi)/4.0 - x))/Sqrt(x) &
     +((-1.0 + 4.0*nu**2)*Sin(((1.0 + 2.0*nu)*Pi)/4.0 &
     - x))/(4.0*Sqrt(2.0*Pi)*x**1.5)
     
     return
     end function
     
    function beslY_asymp ( nu, x )
         use nrtype
     implicit none
     real(dp) beslY_asymp, nu, x
     
     beslY_asymp=((-1.0 + 4.0*nu**2)*Cos(((1.0 + 2.0*nu)*Pi)/4.0 &
     - x))/(4.0*Sqrt(2.0*Pi)*x**1.5) - (Sqrt(2.0/Pi)*Sin(((1.0 &
     + 2.0*nu)*Pi)/4.0 - x))/Sqrt(x)
     
     return
     end function
     
  

  
  subroutine rjbesl ( x, alpha, nb, b, ncalc )

!*****************************************************************************80
!
!! RJBESL calculates J Bessel function with non-integer orders.
!
!  Discussion:
!
!    This routine calculates Bessel functions J sub(N+ALPHA) (X)
!    for non-negative argument X, and non-negative order N+ALPHA.
!
!    This program is based on a program written by David Sookne
!    that computes values of the Bessel functions J or I of real
!    argument and integer order.  Modifications include the restriction
!    of the computation to the J Bessel function of non-negative real
!    argument, the extension of the computation to arbitrary positive
!    order, and the elimination of most underflow.
!
!    In case of an error,  NCALC .NE. NB, and not all J's are
!    calculated to the desired accuracy.
!
!    NCALC < 0:  An argument is out of range. For example,
!    NBES <= 0, ALPHA < 0 or .GT. 1, or X is too large.
!    In this case, B(1) is set to zero, the remainder of the
!    B-vector is not calculated, and NCALC is set to
!    MIN(NB,0)-1 so that NCALC .NE. NB.
!
!    NB .GT. NCALC .GT. 0: Not all requested function values could
!    be calculated accurately.  This usually occurs because NB is
!    much larger than ABS(X).  In this case, B(N) is calculated
!    to the desired accuracy for N <= NCALC, but precision
!    is lost for NCALC < N <= NB.  If B(N) does not vanish
!    for N .GT. NCALC (because it is too small to be represented),
!    and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
!    significant figures of B(N) can be trusted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Frank Olver, David Sookne,
!    A Note on Backward Recurrence Algorithms,
!    Mathematics of Computation,
!    Volume 26, 1972, pages 941-947.
!
!    David Sookne,
!    Bessel Functions of Real Argument and Integer Order,
!    NBS Journal of Res. B,
!    Volume 77B, 1973, pages 125-132.
!
!  Parameters:
!
!    Input, real(dp) X, the argument for which the
!    J's are to be calculated.
!
!    Input, real(dp) ALPHA, the fractional part of order for which
!    the J's or exponentially scaled J's (J*exp(X)) are to be calculated.
!    0 <= ALPHA < 1.0.
!
!    Input, integer ( kind = 4 ) NB, the number of functions to be calculated.
!    0 < NB.  The first function calculated is of order ALPHA, and the
!    last is of order (NB - 1 + ALPHA).
!
!    Output, real(dp) B(NB).  If RJBESL terminates normally, with
!    NCALC = NB, then B contains the functions J/ALPHA/(X) through
!    J/NB-1+ALPHA/(X), or the corresponding exponentially scaled functions.
!
!    Output, integer ( kind = 4 ) NCALC, error indicator.  If NCALC = NB, 
!    then all the requested values were calculated to the desired accuracy.
!
!  Local Parameters:
!
!    IT, the number of bits in the mantissa of a working precision
!    variable.
!
!    NSIG, the decimal significance desired.  Should be set to
!    INT(LOG10(2)*IT+1).  Setting NSIG lower will result
!    in decreased accuracy while setting NSIG higher will
!    increase CPU time without increasing accuracy.  The
!    truncation error is limited to a relative error of
!    T=.5*10^(-NSIG).
!
!    ENTEN = 10.0^K, where K is the largest integer such that
!    ENTEN is machine-representable in working precision
!
!    ENSIG = 10.0^NSIG
!
!    RTNSIG = 10.0^(-K) for the smallest integer K such that
!    K .GE. NSIG/4
!
!    ENMTEN, the smallest ABS(X) such that X/4 does not underflow
!
!    XLARGE, the upper limit on the magnitude of X.  If ABS(X)=N,
!    then at least N iterations of the backward recursion
!    will be executed.  The value of 10000.0 is used on
!    every machine.
  use nrtype
  implicit none

  integer ( kind = 4 ) nb

  real(dp) alpha
  real(dp) alpem
  real(dp) alp2em
  real(dp) b(nb)
  real(dp) capp
  real(dp) capq
  real(dp) eighth
  real(dp) em
  real(dp) en
  real(dp) enmten
  real(dp) ensig
  real(dp) enten
  real(dp) fact(25)
  real(dp) gnu
  real(dp) halfx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) magx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbmx
  integer ( kind = 4 ) ncalc
  integer ( kind = 4 ) nend
  integer ( kind = 4 ) nstart
  real(dp) one30
  real(dp) p
  real(dp) pi2
  real(dp) plast
  real(dp) pold
  real(dp) psave
  real(dp) psavel
!  real(dp) r8_gamma
  real(dp) rtnsig
  real(dp) s
  real(dp) sum
  real(dp) t
  real(dp) t1
  real(dp) tempa
  real(dp) tempb
  real(dp) tempc
  real(dp) test
  real(dp) three
  real(dp) three5
  real(dp) tover
  real(dp) twofiv
  real(dp) twopi1
  real(dp) twopi2
  real(dp) x
  real(dp) xc
  real(dp) xin
  real(dp) xk
  real(dp) xlarge
  real(dp) xm
  real(dp) vcos
  real(dp) vsin
  real(dp) z
!
!  Mathematical constants
!
!   PI2    - 2 / PI
!   TWOPI1 - first few significant digits of 2 * PI
!   TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,
!            TWOPI1 + TWOPI2 = 2 * PI to extra precision.
!
  data pi2 / 0.636619772367581343075535d0 /
  data twopi1 / 6.28125d0 /
  data twopi2 / 1.935307179586476925286767d-3 /
  data eighth / 0.125d0 /
  data three / 3.0d0 /
  data twofiv /25.0d0/
  data one30 /130.0d0 /
  data three5 / 35.0d0/
!
!  Machine-dependent parameters
!
  data enten /1.0d38 /
  data ensig / 1.0d17 /
  data rtnsig / 1.0d-4/
  data enmten /1.2d-37 /
  data xlarge / 1.0d4/
!
!  Factorial(N)
!
  data fact / &
   1.0d0, &
   1.0d0, &
   2.0d0, &
   6.0d0, &
   24.0d0, &
   1.2d2, &
   7.2d2, &
   5.04d3, &
   4.032d4, &
   3.6288d5,3.6288d6,3.99168d7,4.790016d8,6.2270208d9, &
   8.71782912d10,1.307674368d12,2.0922789888d13,3.55687428096d14, &
   6.402373705728d15,1.21645100408832d17,2.43290200817664d18, &
   5.109094217170944d19,1.12400072777760768d21, &
   2.585201673888497664d22, &
   6.2044840173323943936d23/
!
!  Check for out of range arguments.
!
  magx = int ( x )

  if ( &
    0 < nb .and. &
    0.0D+00 <= x .and. &
    x <= xlarge .and. &
    0.0D+00 <= alpha .and. &
    alpha < 1.0D+00 ) then
!
!  Initialize result array to zero.
!
    ncalc = nb
    b(1:nb) = 0.0D+00
!
!  Branch to use 2-term ascending series for small X and asymptotic
!  form for large X when NB is not too large.
!
    if ( x < rtnsig ) then
!
!  Two-term ascending series for small X.
!
      tempa = 1.0D+00
      alpem = 1.0D+00 + alpha
      halfx = 0.0D+00

      if ( enmten < x ) then
        halfx = 0.5D+00 * x
      end if

      if ( alpha /= 0.0D+00 ) then
        tempa = halfx ** alpha / ( alpha * r8_gamma ( alpha ) )
      end if

      tempb = 0.0D+00

      if ( 1.0D+00 < x + 1.0D+00 ) then
        tempb = -halfx * halfx
      end if

      b(1) = tempa + tempa * tempb / alpem

      if ( x /= 0.0D+00 .and. b(1) == 0.0D+00 ) then
        ncalc = 0
      end if

      if ( nb /= 1 ) then

        if ( x <= 0.0D+00 ) then

          do n = 2, nb
            b(n) = 0.0D+00
          end do
!
!  Calculate higher order functions.
!
        else

          tempc = halfx
          tover = ( enmten + enmten ) / x

          if ( tempb /= 0.0D+00 ) then
            tover = enmten / tempb
          end if

          do n = 2, nb

            tempa = tempa / alpem
            alpem = alpem + 1.0D+00
            tempa = tempa * tempc

            if ( tempa <= tover * alpem ) then
              tempa = 0.0D+00
            end if

            b(n) = tempa + tempa * tempb / alpem

            if ( b(n) == 0.0D+00 .and. n < ncalc ) then
              ncalc = n - 1
            end if

          end do

        end if
      end if
!
!  Asymptotic series for 21 < X.
!
    else if ( twofiv < x .and. nb <= magx + 1 ) then

      xc = sqrt ( pi2 / x )
      xin = ( eighth / x )**2
      m = 11

      if ( three5 <= x ) then
        m = 8
      end if

      if ( one30 <= x ) then
        m = 4
      end if

      xm = 4.0D+00 * real ( m, kind = 8 )
!
!  Argument reduction for SIN and COS routines.
!
      t = aint ( x / ( twopi1 + twopi2 ) + 0.5D+00 )
      z = ( ( x - t * twopi1 ) - t * twopi2 ) &
        - ( alpha + 0.5D+00 ) / pi2
      vsin = sin ( z )
      vcos = cos ( z )
      gnu = alpha + alpha

      do i = 1, 2

        s = ( ( xm - 1.0D+00 ) - gnu ) * ( ( xm - 1.0D+00 ) + gnu ) &
          * xin * 0.5D+00
        t = ( gnu - ( xm - three ) ) * ( gnu + ( xm - three ) )
        capp = s * t / fact(2*m+1)
        t1 = ( gnu - ( xm + 1.0D+00 ) ) * ( gnu + ( xm + 1.0D+00 ) )
        capq = s * t1 / fact(2*m+2)
        xk = xm
        k = m + m
        t1 = t

        do j = 2, m
          xk = xk - 4.0D+00
          s = ( ( xk - 1.0D+00 ) - gnu ) * ( ( xk - 1.0D+00 ) + gnu )
          t = ( gnu - ( xk - three ) ) * ( gnu + ( xk - three ) )
          capp = ( capp + 1.0D+00 / fact(k-1) ) * s * t * xin
          capq = ( capq + 1.0D+00 / fact(k) ) * s * t1 * xin
          k = k - 2
          t1 = t
        end do

        capp = capp + 1.0D+00
        capq = ( capq + 1.0D+00 ) * ( gnu * gnu - 1.0D+00 ) * ( eighth / x )
        b(i) = xc * ( capp * vcos - capq * vsin )

        if ( nb == 1 ) then
          return
        end if

        t = vsin
        vsin = -vcos
        vcos = t
        gnu = gnu + 2.0D+00

      end do
!
!  If 2 < NB, compute J(X,ORDER+I)  I = 2, NB-1.
!
      if ( 2 < nb ) then
        gnu = alpha + alpha + 2.0D+00
        do j = 3, nb
          b(j) = gnu * b(j-1) / x - b(j-2)
          gnu = gnu + 2.0D+00
        end do
      end if
!
!  Use recurrence to generate results.  First initialize the
!  calculation of P's.
!
    else

      nbmx = nb - magx
      n = magx + 1
      en = real ( n + n, kind = 8 ) + ( alpha + alpha )
      plast = 1.0D+00
      p = en / x
!
!  Calculate general significance test.
!
      test = ensig + ensig
!
!  Calculate P's until N = NB-1.  Check for possible overflow.
!
      if ( 3 <= nbmx ) then

        tover = enten / ensig
        nstart = magx + 2
        nend = nb - 1
        en = real ( nstart + nstart, kind = 8 ) - 2.0D+00 + ( alpha + alpha )

        do k = nstart, nend

          n = k
          en = en + 2.0D+00
          pold = plast
          plast = p
          p = en * plast / x - pold
!
!  To avoid overflow, divide P's by TOVER.  Calculate P's until
!  1 < ABS(P).
!
          if ( tover < p ) then

            tover = enten
            p = p / tover
            plast = plast / tover
            psave = p
            psavel = plast
            nstart = n + 1

            do

              n = n + 1
              en = en + 2.0D+00
              pold = plast
              plast = p
              p = en * plast / x - pold

              if ( 1.0D+00 < p ) then
                exit
              end if

            end do

            tempb = en / x
!
!  Calculate backward test and find NCALC, the highest N such that
!  the test is passed.
!
            test = pold * plast &
              * ( 0.5D+00 - 0.5D+00 / ( tempb * tempb ) )
            test = test / ensig
            p = plast * tover
            n = n - 1
            en = en - 2.0D+00
            nend = min ( nb, n )

            do l = nstart, nend
              pold = psavel
              psavel = psave
              psave = en * psavel / x - pold
              if ( test < psave * psavel ) then
                ncalc = l - 1
                go to 190
              end if
            end do

            ncalc = nend
            go to 190

          end if

        end do

        n = nend
        en = real ( n + n, kind = 8 ) + ( alpha + alpha )
!
!  Calculate special significance test for 2 < NBMX.
!
        test = max ( test, sqrt ( plast * ensig ) * sqrt ( p + p ) )

      end if
!
!  Calculate P's until significance test passes.
!
      do

        n = n + 1
        en = en + 2.0D+00
        pold = plast
        plast = p
        p = en * plast / x - pold

        if ( test <= p ) then
          exit
        end if

      end do
!
!  Initialize the backward recursion and the normalization sum.
!
  190     continue

      n = n + 1
      en = en + 2.0D+00
      tempb = 0.0D+00
      tempa = 1.0D+00 / p
      m = 2 * n - 4 * ( n / 2 )
      sum = 0.0D+00
      em = real ( n / 2, kind = 8 )
      alpem = ( em - 1.0D+00 ) + alpha
      alp2em = ( em + em ) + alpha

      if ( m /= 0 ) then
        sum = tempa * alpem * alp2em / em
      end if

      nend = n - nb
!
!  Recur backward via difference equation, calculating (but not
!  storing) B(N), until N = NB.
!
      if ( 0 < nend ) then

        do l = 1, nend

          n = n - 1
          en = en - 2.0D+00
          tempc = tempb
          tempb = tempa
          tempa = ( en * tempb ) / x - tempc
          m = 2 - m

          if ( m /= 0 ) then
            em = em - 1.0D+00
            alp2em = ( em + em ) + alpha
            if ( n == 1 ) then
              exit
            end if
            alpem = ( em - 1.0D+00 ) + alpha
            if ( alpem == 0.0D+00 ) then
              alpem = 1.0D+00
            end if
            sum = ( sum + tempa * alp2em ) * alpem / em
          end if

        end do

      end if
!
!  Store B(NB).
!
! 210     continue

      b(n) = tempa

      if ( 0 <= nend ) then

        if ( nb <= 1 ) then

          alp2em = alpha
          if ( alpha + 1.0D+00 == 1.0D+00 ) then
            alp2em = 1.0D+00
          end if
          sum = sum + b(1) * alp2em
          go to 250

        else
!
!  Calculate and store B(NB-1).
!
          n = n - 1
          en = en - 2.0D+00
          b(n) = ( en * tempa ) / x - tempb

          if ( n == 1 ) then
            go to 240
          end if

          m = 2 - m

          if ( m /= 0 ) then
            em = em - 1.0D+00
            alp2em = ( em + em ) + alpha
            alpem = ( em - 1.0D+00 ) + alpha
            if ( alpem == 0.0D+00 ) then
              alpem = 1.0D+00
            end if
            sum = ( sum + b(n) * alp2em ) * alpem / em
          end if

        end if

      end if

      nend = n - 2
!
!  Calculate via difference equation and store B(N), until N = 2.
!
      if ( nend /= 0 ) then

        do l = 1, nend
          n = n - 1
          en = en - 2.0D+00
          b(n) = ( en * b(n+1) ) / x - b(n+2)
          m = 2 - m
          if ( m /= 0 ) then
            em = em - 1.0D+00
            alp2em = ( em + em ) + alpha
            alpem = ( em - 1.0D+00 ) + alpha
            if ( alpem == 0.0D+00 ) then
              alpem = 1.0D+00
            end if
            sum = ( sum + b(n) * alp2em ) * alpem / em
          end if
        end do

      end if
!
!  Calculate B(1).
!
      b(1) = 2.0D+00 * ( alpha + 1.0D+00 ) * b(2) / x - b(3)

  240     continue

      em = em - 1.0D+00
      alp2em = ( em + em ) + alpha

      if ( alp2em == 0.0D+00 ) then
        alp2em = 1.0D+00
      end if

      sum = sum + b(1) * alp2em
!
!  Normalize.  Divide all B(N) by sum.
!
  250     continue

      if ( alpha + 1.0D+00 /= 1.0D+00 ) then
        sum = sum * r8_gamma ( alpha ) * ( x * 0.5D+00 ) ** ( -alpha )
      end if

      tempa = enmten
      if ( 1.0D+00 < sum ) then
        tempa = tempa * sum
      end if

      do n = 1, nb
        if ( abs ( b(n) ) < tempa ) then
          b(n) = 0.0D+00
        end if
        b(n) = b(n) / sum
      end do

    end if
!
!  Error return: X, NB, or ALPHA is out of range.
!
  else
    b(1) = 0.0D+00
    ncalc = min ( nb, 0 ) - 1
  end if

  return
end subroutine


subroutine rybesl ( x, alpha, nb, by, ncalc )

!*****************************************************************************80
!
!! RYBESL calculates Y Bessel function with non-integer orders.
!
!  Discussion:
!
!    This routine calculates Bessel functions Y SUB(N+ALPHA) (X)
!    for non-negative argument X, and non-negative order N+ALPHA.
!
!    This program draws heavily on Temme's Algol program for Y(a,x)
!    and Y(a+1,x) and on Campbell's programs for Y_nu(x).  Temme's
!    scheme is used for x < THRESH, and Campbell's scheme is used
!    in the asymptotic region.  Segments of code from both sources
!    have been translated into Fortran77, merged, and heavily modified.
!    Modifications include parameterization of machine dependencies,
!    use of a new approximation for ln(gamma(x)), and built-in
!    protection against over/underflow.
!
!    In case of an error, NCALC .NE. NB, and not all Y's are
!    calculated to the desired accuracy.
!
!    NCALC < -1:  An argument is out of range. For example,
!    NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
!    XMAX.  In this case, BY(1) = 0.0, the remainder of the
!    BY-vector is not calculated, and NCALC is set to
!    MIN0(NB,0)-2  so that NCALC .NE. NB.
!
!    NCALC = -1:  Y(ALPHA,X) .GE. XINF.  The requested function
!    values are set to 0.0.
!
!    1 < NCALC < NB: Not all requested function values could
!    be calculated accurately.  BY(I) contains correct function
!    values for I <= NCALC, and and the remaining NB-NCALC
!    array elements contain 0.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2013
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    JB Campbell,
!    Bessel functions J_nu(x) and Y_nu(x) of real order and real argument,
!    Computational Physics Communications,
!    Volume 18, 1979, pages 133-142.
!
!    NM Temme,
!    On the numerical evaluation of the ordinary Bessel function
!    of the second kind,
!    Journal of Computational Physics,
!    Volume 21, 1976, pages 343-350.
!
!  Parameters:
!
!    Input, real(dp) X, the argument.  0 <= X.
!
!    Input, real(dp) ALPHA, the fractional part of the order
!    for which the Y's are to be calculated.  0 <= ALPHA < 1.0.
!
!    Input, integer ( kind = 4 ) NB, the number of functions to be calculated, 
!    NB .GT. 0.  The first function calculated is of order ALPHA, and the
!    last is of order (NB - 1 + ALPHA).
!
!    Output, real(dp) BY(NB).  If the routine terminates normally
!    (NCALC=NB), the vector BY contains the functions Y(ALPHA,X) through
!    Y(NB-1+ALPHA,X),  If (0 < NCALC < NB), BY(I) contains correct
!    function values for I <= NCALC, and contains the ratios
!    Y(ALPHA+I-1,X)/Y(ALPHA+I-2,X) for the rest of the array.
!
!    Output, integer ( kind = 4 ) NCALC, error flag.  Before using the vector 
!    BY, the user should check that NCALC=NB, i.e., all orders have been 
!    calculated to the desired accuracy.
!
  use nrtype
  implicit none

  integer ( kind = 4 ) nb

  real(dp) alfa
  real(dp) alpha
  real(dp) aye
  real(dp) b
  real(dp) by(nb)
  real(dp) c
  real(dp) ch(21)
  real(dp) cosmu
  real(dp) d
  real(dp) del
  real(dp) den
  real(dp) ddiv
  real(dp) div
  real(dp) dmu
  real(dp) d1
  real(dp) d2
  real(dp) e
  real(dp) eight
  real(dp) en
  real(dp) enu
  real(dp) en1
  real(dp) even
  real(dp) ex
  real(dp) f
  real(dp) fivpi
  real(dp) g
  real(dp) gamma
  real(dp) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  integer ( kind = 4 ) ncalc
  real(dp) odd
  real(dp) onbpi
  real(dp) one5
  real(dp) p
  real(dp) pa
  real(dp) pa1
  !real(dp), parameter :: pi = 3.1415926535897932384626434D+00
  real(dp) piby2
  real(dp) pim5
  real(dp) q
  real(dp) qa
  real(dp) qa1
  real(dp) q0
  real(dp) r
  real(dp) s
  real(dp) sinmu
  real(dp) sq2bpi
  real(dp) ten9
  real(dp) term
  real(dp) three
  real(dp) thresh
  real(dp) twobyx
  real(dp) x
  real(dp) xinf
  real(dp) xlarge
  real(dp) xmin
  real(dp) xna
  real(dp) x2
  real(dp) ya
  real(dp) ya1
!
!  Mathematical constants
!    FIVPI = 5*PI
!    PIM5 = 5*PI - 15
!    ONBPI = 1/PI
!    PIBY2 = PI/2
!    SQ2BPI = SQUARE ROOT OF 2/PI
!
  data three / 3.0d0 /
  data eight /8.0d0 /
  data one5 / 15.0d0 /
  data ten9 / 1.9d1/
  data fivpi /1.5707963267948966192d1 /
  data piby2 / 1.5707963267948966192d0/
  data sq2bpi / 7.9788456080286535588d-1/
  data pim5 /7.0796326794896619231d-1/
  data onbpi / 3.1830988618379067154d-1/
!
!  Machine-dependent constants
!
  data del / 1.0d-8 /
  data xmin / 4.46d-308 /
  data xinf / 1.79d308 /
  data thresh / 16.0d0 /
  data xlarge / 1.0d8 /
!
!  Coefficients for Chebyshev polynomial expansion of
!  1/gamma(1-x), abs(x) <= .5
!
  data ch/-0.67735241822398840964d-23,-0.61455180116049879894d-22, &
           0.29017595056104745456d-20, 0.13639417919073099464d-18, &
           0.23826220476859635824d-17,-0.90642907957550702534d-17, &
          -0.14943667065169001769d-14,-0.33919078305362211264d-13, &
          -0.17023776642512729175d-12, 0.91609750938768647911d-11, &
           0.24230957900482704055d-09, 0.17451364971382984243d-08, &
          -0.33126119768180852711d-07,-0.86592079961391259661d-06, &
          -0.49717367041957398581d-05, 0.76309597585908126618d-04, &
           0.12719271366545622927d-02, 0.17063050710955562222d-02, &
          -0.76852840844786673690d-01,-0.28387654227602353814d+00, &
           0.92187029365045265648d+00/

  ex = x
  enu = alpha

  if ( 0 < nb .and. &
    xmin <= x .and. &
    ex < xlarge .and. &
    0.0D+00 <= enu .and. &
    enu < 1.0D+00 )  then

    xna = aint ( enu + 0.5D+00 )
    na = int ( xna )

    if ( na == 1 ) then
      enu = enu - xna
    end if

    if ( enu == - 0.5D+00 ) then

      p = sq2bpi / sqrt ( ex )
      ya = p * sin ( ex )
      ya1 = -p * cos ( ex )
!
!  Use Temme's scheme for small X.
!
    else if ( ex < three ) then

      b = ex * 0.5D+00
      d = - log ( b )
      f = enu * d
      e = b**( -enu )

      if ( abs ( enu ) < del ) then
        c = onbpi
      else
        c = enu / sin ( enu * pi )
      end if
!
!  Computation of sinh(f)/f.
!
      if ( abs ( f ) < 1.0D+00 ) then
        x2 = f * f
        en = ten9
        s = 1.0D+00
        do i = 1, 9
          s = s * x2 / en / ( en - 1.0D+00 ) + 1.0D+00
          en = en - 2.0D+00
        end do
      else
        s = ( e - 1.0D+00 / e ) * 0.5D+00 / f
      end if
!
!  Computation of 1/gamma(1-a) using Chebyshev polynomials.
!
      x2 = enu * enu * eight
      aye = ch(1)
      even = 0.0D+00
      alfa = ch(2)
      odd = 0.0D+00

      do i = 3, 19, 2
        even = -( aye + aye + even )
        aye = -even * x2 - aye + ch(i)
        odd = -( alfa + alfa + odd )
        alfa = -odd * x2 - alfa + ch(i+1)
      end do

      even = ( even * 0.5D+00 + aye ) * x2 - aye + ch(21)
      odd = ( odd + alfa ) * 2.0D+00
      gamma = odd * enu + even
!
!  End of computation of 1/gamma(1-a).
!
      g = e * gamma
      e = ( e + 1.0D+00 / e ) * 0.5D+00
      f = 2.0D+00 * c * ( odd * e + even * s * d )
      e = enu * enu
      p = g * c
      q = onbpi / g
      c = enu * piby2

      if ( abs ( c ) < del ) then
        r = 1.0D+00
      else
        r = sin ( c ) / c
      end if

      r = pi * c * r * r
      c = 1.0D+00
      d = - b * b
      h = 0.0D+00
      ya = f + r * q
      ya1 = p
      en = 0.0D+00

      do

        en = en + 1.0D+00

        if ( abs ( g / ( 1.0D+00 + abs ( ya ) ) ) &
          + abs ( h / ( 1.0D+00 + abs ( ya1 ) ) ) <= epsilon ( g ) ) then
          exit
        end if

        f = ( f * en + p + q ) / ( en * en - e )
        c = c * d / en
        p = p / ( en - enu )
        q = q / ( en + enu )
        g = c * ( f + r * q )
        h = c * p - en * g
        ya = ya + g
        ya1 = ya1 + h

      end do

      ya = -ya
      ya1 = -ya1 / b

    else if ( ex < thresh ) then
!
!  Use Temme's scheme for moderate X.
!
      c = ( 0.5D+00 - enu ) * ( 0.5D+00 + enu )
      b = ex + ex
      e = ( ex * onbpi * cos ( enu * pi ) / epsilon ( e ) )
      e = e * e
      p = 1.0D+00
      q = -ex
      r = 1.0D+00 + ex * ex
      s = r
      en = 2.0D+00

      do while ( r * en * en < e )
        en1 = en + 1.0D+00
        d = ( en - 1.0D+00 + c / en ) / s
        p = ( en + en - p * d ) / en1
        q = ( -b + q * d ) / en1
        s = p * p + q * q
        r = r * s
        en = en1
      end do

      f = p / s
      p = f
      g = -q / s
      q = g

      do

        en = en - 1.0D+00

        if ( en <= 0.0D+00 ) then
          exit
        end if

        r = en1 * ( 2.0D+00 - p ) - 2.0D+00
        s = b + en1 * q
        d = ( en - 1.0D+00 + c / en ) / ( r * r + s * s )
        p = d * r
        q = d * s
        e = f + 1.0D+00
        f = p * e - g * q
        g = q * e + p * g
        en1 = en

      end do

      f = 1.0D+00 + f
      d = f * f + g * g
      pa = f / d
      qa = -g / d
      d = enu + 0.5D+00 - p
      q = q + ex
      pa1 = ( pa * q - qa * d ) / ex
      qa1 = ( qa * q + pa * d ) / ex
      b = ex - piby2 * ( enu + 0.5D+00 )
      c = cos ( b )
      s = sin ( b )
      d = sq2bpi / sqrt ( ex )
      ya = d * ( pa * s + qa * c )
      ya1 = d * ( qa1 * s - pa1 * c )
!
!  Use Campbell's asymptotic scheme.
!
    else

      na = 0
      d1 = aint ( ex / fivpi )
      i = int ( d1 )
      dmu = (( ex - one5 * d1 ) - d1 * pim5 ) &
        - ( alpha + 0.5D+00 ) * piby2

      if ( i - 2 * ( i / 2 ) == 0 ) then
        cosmu = cos ( dmu )
        sinmu = sin ( dmu )
      else
        cosmu = -cos ( dmu )
        sinmu = -sin ( dmu )
      end if

      ddiv = eight * ex
      dmu = alpha
      den = sqrt ( ex )

      do k = 1, 2

        p = cosmu
        cosmu = sinmu
        sinmu = -p
        d1 = ( 2.0D+00 * dmu - 1.0D+00 ) * ( 2.0D+00 * dmu + 1.0D+00 )
        d2 = 0.0D+00
        div = ddiv
        p = 0.0D+00
        q = 0.0D+00
        q0 = d1 / div
        term = q0

        do i = 2, 20

          d2 = d2 + eight
          d1 = d1 - d2
          div = div + ddiv
          term = -term * d1 / div
          p = p + term
          d2 = d2 + eight
          d1 = d1 - d2
          div = div + ddiv
          term = term * d1 / div
          q = q + term

          if ( abs ( term ) <= epsilon ( term ) ) then
            exit
          end if

        end do

        p = p + 1.0D+00
        q = q + q0

        if ( k == 1 ) then
          ya = sq2bpi * ( p * cosmu - q * sinmu ) / den
        else
          ya1 = sq2bpi * ( p * cosmu - q * sinmu ) / den
        end if

        dmu = dmu + 1.0D+00

      end do

    end if

    if ( na == 1 ) then
      h = 2.0D+00 * ( enu + 1.0D+00 ) / ex
      if ( 1.0D+00 < h ) then
        if ( xinf / h < abs ( ya1 ) ) then
          h = 0.0D+00
          ya = 0.0D+00
        end if
      end if
      h = h * ya1 - ya
      ya = ya1
      ya1 = h
    end if
!
!  Now have first one or two Y's.
!
    by(1) = ya
    by(2) = ya1

    if ( ya1 == 0.0D+00 ) then

      ncalc = 1

    else

      aye = 1.0D+00 + alpha
      twobyx = 2.0D+00 / ex
      ncalc = 2

      do i = 3, nb

        if ( twobyx < 1.0D+00 ) then
          if ( xinf / aye <= abs ( by(i-1) ) * twobyx ) then
            exit
          end if
        else
          if ( xinf / aye / twobyx <= abs ( by(i-1) ) ) then
            exit
          end if
        end if

        by(i) = twobyx * aye * by(i-1) - by(i-2)
        aye = aye + 1.0D+00
        ncalc = ncalc + 1

      end do

    end if

    do i = ncalc + 1, nb
      by(i) = 0.0D+00
    end do

  else

    by(1) = 0.0D+00
    ncalc = min ( nb, 0 ) - 1

  end if

  return
end subroutine

function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates the gamma function.
!
!  Discussion:
!
!    This function was originally named DGAMMA.
!
!    However, a number of compilers include a library function of this name.
!    To avoid conflicts, this function was renamed R8_GAMMA.
!
!    This routine calculates the GAMMA function for a real argument X.
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the GAMMA
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real(dp) X, the argument of the function.
!
!    Output, real(dp) R8_GAMMA, the value of the function.
!
  use nrtype
  implicit none
!
!  Coefficients for minimax approximation over (12, INF).
!
  real(dp), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real(dp) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real(dp), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
 ! real(dp), parameter :: pi = 3.1415926535897932384626434D+00
  real(dp), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real(dp) r8_gamma
  real(dp) res
  real(dp), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real(dp) sum
  real(dp) x
  real(dp), parameter :: xbig = 171.624D+00
  real(dp) xden
  real(dp), parameter :: xinf = 1.79D+308
  real(dp), parameter :: xminin = 2.23D-308
  real(dp) xnum
  real(dp) y
  real(dp) y1
  real(dp) ysq
  real(dp) z

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < epsilon ( y ) ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = 1.0D+00 / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + 1.0D+00
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5D+00 ) * log ( y )
      res = exp ( sum )

    else

      res = huge ( res )
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end function


    SUBROUTINE bessjy_s(x,xnu,rj,ry,rjp,ryp)
	use nrtype; USE nrutil, ONLY : assert,nrerror
	USE nr, ONLY : beschb
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x,xnu
	REAL(DP), INTENT(OUT) :: rj,ry,rjp,ryp
	INTEGER(I4B), PARAMETER :: MAXIT=10000
	REAL(DP), PARAMETER :: XMIN=2.0_dp,EPS=1.0e-10_dp,FPMIN=1.0e-30_dp
	INTEGER(I4B) :: i,isign,l,nl
	REAL(DP) :: a,b,c,d,del,del1,e,f,fact,fact2,fact3,ff,gam,gam1,gam2,&
		gammi,gampl,h,p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,&
		ry1,rymu,rymup,rytemp,sum,sum1,w,x2,xi,xi2,xmu,xmu2
	COMPLEX(DPC) :: aa,bb,cc,dd,dl,pq
	call assert(x > 0.0, xnu >= 0.0, 'bessjy args')
	nl=merge(int(xnu+0.5_dp), max(0,int(xnu-x+1.5_dp)), x < XMIN)
	xmu=xnu-nl
	xmu2=xmu*xmu
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	w=xi2/PI_D
	isign=1
	h=xnu*xi
	if (h < FPMIN) h=FPMIN
	b=xi2*xnu
	d=0.0
	c=h
	do i=1,MAXIT
		b=b+xi2
		d=b-d
		if (abs(d) < FPMIN) d=FPMIN
		c=b-1.0_dp/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_dp/d
		del=c*d
		h=del*h
		if (d < 0.0) isign=-isign
		if (abs(del-1.0_dp) < EPS) exit
	end do
	if (i > MAXIT) call nrerror('x too large in bessjy; try asymptotic expansion')
	rjl=isign*FPMIN
	rjpl=h*rjl
	rjl1=rjl
	rjp1=rjpl
	fact=xnu*xi
	do l=nl,1,-1
		rjtemp=fact*rjl+rjpl
		fact=fact-xi
		rjpl=fact*rjtemp-rjl
		rjl=rjtemp
	end do
	if (rjl == 0.0) rjl=EPS
	f=rjpl/rjl
	if (x < XMIN) then
		x2=0.5_dp*x
		pimu=PI_D*xmu
		if (abs(pimu) < EPS) then
			fact=1.0
		else
			fact=pimu/sin(pimu)
		end if
		d=-log(x2)
		e=xmu*d
		if (abs(e) < EPS) then
			fact2=1.0
		else
			fact2=sinh(e)/e
		end if
		call beschb(xmu,gam1,gam2,gampl,gammi)
		ff=2.0_dp/PI_D*fact*(gam1*cosh(e)+gam2*fact2*d)
		e=exp(e)
		p=e/(gampl*PI_D)
		q=1.0_dp/(e*PI_D*gammi)
		pimu2=0.5_dp*pimu
		if (abs(pimu2) < EPS) then
			fact3=1.0
		else
			fact3=sin(pimu2)/pimu2
		end if
		r=PI_D*pimu2*fact3*fact3
		c=1.0
		d=-x2*x2
		sum=ff+r*q
		sum1=p
		do i=1,MAXIT
			ff=(i*ff+p+q)/(i*i-xmu2)
			c=c*d/i
			p=p/(i-xmu)
			q=q/(i+xmu)
			del=c*(ff+r*q)
			sum=sum+del
			del1=c*p-i*del
			sum1=sum1+del1
			if (abs(del) < (1.0_dp+abs(sum))*EPS) exit
		end do
		if (i > MAXIT) call nrerror('bessy series failed to converge')
		rymu=-sum
		ry1=-sum1*xi2
		rymup=xmu*xi*rymu-ry1
		rjmu=w/(rymup-f*rymu)
	else
		a=0.25_dp-xmu2
		pq=cmplx(-0.5_dp*xi,1.0_dp,kind=dpc)
		aa=cmplx(0.0_dp,xi*a,kind=dpc)
		bb=cmplx(2.0_dp*x,2.0_dp,kind=dpc)
		cc=bb+aa/pq
		dd=1.0_dp/bb
		pq=cc*dd*pq
		do i=2,MAXIT
			a=a+2*(i-1)
			bb=bb+cmplx(0.0_dp,2.0_dp,kind=dpc)
			dd=a*dd+bb
			if (absc(dd) < FPMIN) dd=FPMIN
			cc=bb+a/cc
			if (absc(cc) < FPMIN) cc=FPMIN
			dd=1.0_dp/dd
			dl=cc*dd
			pq=pq*dl
			if (absc(dl-1.0_dp) < EPS) exit
		end do
		if (i > MAXIT) call nrerror('cf2 failed in bessjy')
		p=real(pq)
		q=aimag(pq)
		gam=(p-f)/q
		rjmu=sqrt(w/((p-f)*gam+q))
		rjmu=sign(rjmu,rjl)
		rymu=rjmu*gam
		rymup=rymu*(p+q/gam)
		ry1=xmu*xi*rymu-rymup
	end if
	fact=rjmu/rjl
	rj=rjl1*fact
	rjp=rjp1*fact
	do i=1,nl
		rytemp=(xmu+i)*xi2*ry1-rymu
		rymu=ry1
		ry1=rytemp
	end do
	ry=rymu
	ryp=xnu*xi*rymu-ry1
	CONTAINS
!BL
	FUNCTION absc(z)
	use nrtype
	IMPLICIT NONE
	COMPLEX(DPC), INTENT(IN) :: z
	REAL(DP) :: absc
	absc=abs(real(z))+abs(aimag(z))
	END FUNCTION absc
	END SUBROUTINE bessjy_s

	SUBROUTINE bessjy_v(x,xnu,rj,ry,rjp,ryp)
	use nrtype; USE nrutil, ONLY : assert,nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: xnu
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
	INTEGER(I4B), PARAMETER :: MAXIT=1000000
	REAL(DP), PARAMETER :: XMIN=2.0_dp,EPS=1.0e-16_dp,FPMIN=1.0e-30_dp
	INTEGER(I4B) :: i,n
	INTEGER(I4B), DIMENSION(size(x)) :: isign
	REAL(DP), DIMENSION(size(x)) :: b,c,d,del,&
		h,xi,xi2,&
		rj_lt,ry_lt,rjp_lt,ryp_lt,rj_ge,ry_ge,rjp_ge,ryp_ge
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	call assert(all(x > 0.0), xnu >= 0.0, 'bessjy args')
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	isign=1
	h=xnu*xi
	where (h < FPMIN) h=FPMIN
	b=xi2*xnu
	d=0.0
	c=h
	converged=.false.
	do i=1,MAXIT
		where (.not. converged)
			b=b+xi2
			d=b-d
			d=merge(FPMIN,d, abs(d) < FPMIN )
			c=b-1.0_dp/c
			c=merge(FPMIN,c, abs(c) < FPMIN )
			d=1.0_dp/d
			del=c*d
			h=del*h
			isign=merge(-isign,isign, d < 0.0 )
			converged=(abs(del-1.0_dp) < EPS)
		end where
		if (all(converged)) exit
	end do
	if (i > MAXIT) call nrerror('bessjy: x too large; try asymptotic expansion')
	converged=(x < XMIN)
	n=count(converged)
	call bessjy_xltxmin(pack(x,converged),xnu,pack(isign*FPMIN,converged),pack(h,converged),&
		rj_lt(1:n),ry_lt(1:n),rjp_lt(1:n),ryp_lt(1:n))
	n=size(x)-n
	call bessjy_xgexmin(pack(x,.not. converged),xnu,pack(isign,.not. converged),pack(h,.not. converged),&
		rj_ge(1:n),ry_ge(1:n),rjp_ge(1:n),ryp_ge(1:n))
	rj=unpacked(converged,rj_lt,rj_ge)
	rjp=unpacked(converged,rjp_lt,rjp_ge)
	ry=unpacked(converged,ry_lt,ry_ge)
	ryp=unpacked(converged,ryp_lt,ryp_ge)
	CONTAINS
!BL
	SUBROUTINE bessjy_xltxmin(x,xnu,rjli,hi,rj,ry,rjp,ryp)
	use nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : beschb
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), INTENT(IN) :: xnu
	REAL(DP), DIMENSION(:), INTENT(IN) :: hi,rjli
	REAL(DP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
	INTEGER(I4B) :: i,nl
	REAL(DP) :: xmu,xmu2,gam1,gam2,gammi,gampl
	REAL(DP), DIMENSION(size(x)) :: &
		c,d,del,del1,e,f,fact,fact2,fact3,ff,h,&
		p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjpl,rjp1,rjtemp,&
		ry1,rymu,rymup,rytemp,sum,sum1,w,x2,xi,xi2
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	nl=int(xnu+0.5_dp)
	xmu=xnu-nl
	xmu2=xmu*xmu
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	w=xi2/PI_D
	h=hi
	rjl=rjli
	rjpl=h*rjl
	rjl1=rjl
	rjp1=rjpl
	fact=xnu*xi
	do i=int(xnu+0.5_dp),1,-1
		rjtemp=fact*rjl+rjpl
		fact=fact-xi
		rjpl=fact*rjtemp-rjl
		rjl=rjtemp
	end do
	where (rjl == 0.0) rjl=EPS
	f=rjpl/rjl
	x2=0.5_dp*x
	pimu=PI_D*xmu
	where (abs(pimu) < EPS)
		fact=1.0
	elsewhere
		fact=pimu/sin(pimu)
	end where
	d=-log(x2)
	e=xmu*d
	where (abs(e) < EPS)
		fact2=1.0
	elsewhere
		fact2=sinh(e)/e
	end where
	call beschb(xmu,gam1,gam2,gampl,gammi)
	ff=2.0_dp/PI_D*fact*(gam1*cosh(e)+gam2*fact2*d)
	e=exp(e)
	p=e/(gampl*PI_D)
	q=1.0_dp/(e*PI_D*gammi)
	pimu2=0.5_dp*pimu
	where (abs(pimu2) < EPS)
		fact3=1.0
	elsewhere
		fact3=sin(pimu2)/pimu2
	end where
	r=PI_D*pimu2*fact3*fact3
	c=1.0
	d=-x2*x2
	sum=ff+r*q
	sum1=p
	converged=.false.
	do i=1,MAXIT
		where (.not. converged)
			ff=(i*ff+p+q)/(i*i-xmu2)
			c=c*d/i
			p=p/(i-xmu)
			q=q/(i+xmu)
			del=c*(ff+r*q)
			sum=sum+del
			del1=c*p-i*del
			sum1=sum1+del1
			converged=(abs(del) < (1.0_dp+abs(sum))*EPS)
		end where
		if (all(converged)) exit
	end do
	if (i > MAXIT) call nrerror('bessy series failed to converge')
	rymu=-sum
	ry1=-sum1*xi2
	rymup=xmu*xi*rymu-ry1
	rjmu=w/(rymup-f*rymu)
	do i=1,nl
		rytemp=(xmu+i)*xi2*ry1-rymu
		rymu=ry1
		ry1=rytemp
	end do
	fact=rjmu/rjl
	rj=rjl1*fact
	rjp=rjp1*fact
	ry=rymu
	ryp=xnu*xi*rymu-ry1
	END SUBROUTINE bessjy_xltxmin
!BL
	SUBROUTINE bessjy_xgexmin(x,xnu,isign,h,rj,ry,rjp,ryp)
	use nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), INTENT(IN) :: xnu
	INTEGER(I4B),DIMENSION(:), INTENT(IN) :: isign
	REAL(DP), DIMENSION(:), INTENT(IN) :: h
	REAL(DP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
	INTEGER(I4B) :: i,nlmax
	INTEGER(I4B),DIMENSION(size(x)) :: nl
	REAL(DP), DIMENSION(size(x)) :: &
		a,f,fact,gam,p,q,rjl,rjl1,rjmu,rjpl,rjp1,rjtemp,ry1,rymu,rymup,rytemp,w,xi,xi2,xmu,xmu2
	COMPLEX(DPC), DIMENSION(size(x)) :: aa,bb,cc,dd,dl,pq
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	nl=max(0,int(xnu-x+1.5_dp))
	nlmax=maxval(nl)
	xmu=xnu-nl
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	w=xi2/PI_D
	xmu2=xmu*xmu
	rjl=isign*FPMIN
	rjpl=h*rjl
	rjl1=rjl
	rjp1=rjpl
	fact=xnu*xi
	do i=nlmax,1,-1
		converged=(i > nl)
		if (all(converged)) exit
		where (.not. converged)
			rjtemp=fact*rjl+rjpl
			fact=fact-xi
			rjpl=fact*rjtemp-rjl
			rjl=rjtemp
		end where
	end do
	where (rjl == 0.0) rjl=EPS
	f=rjpl/rjl
	a=0.25_dp-xmu2
	pq=cmplx(-0.5_dp*xi,1.0_dp,kind=dpc)
	aa=cmplx(0.0_dp,xi*a,kind=dpc)
	bb=cmplx(2.0_dp*x,2.0_dp,kind=dpc)
	cc=bb+aa/pq
	dd=1.0_dp/bb
	pq=cc*dd*pq
	converged=.false.
	do i=2,MAXIT
		where (.not. converged)
			a=a+2*(i-1)
			bb=bb+cmplx(0.0_dp,2.0_dp,kind=dpc)
			dd=a*dd+bb
			dd=merge(cmplx(FPMIN,kind=dpc),dd, absc(dd) < FPMIN )
			cc=bb+a/cc
			cc=merge(cmplx(FPMIN,kind=dpc),cc, absc(cc) < FPMIN )
			dd=1.0_dp/dd
			dl=cc*dd
			pq=pq*dl
			converged=(absc(dl-1.0_dp) < EPS)
		end where
		if (all(converged)) exit
	end do
	if (i > MAXIT) call nrerror('bessjy: cf2 section failed')
	p=real(pq,dp)
	q=aimag(pq)
	gam=(p-f)/q
	rjmu=sqrt(w/((p-f)*gam+q))
	rjmu=sign(rjmu,rjl)
	rymu=rjmu*gam
	rymup=rymu*(p+q/gam)
	ry1=xmu*xi*rymu-rymup
	do i=1,nlmax
		converged=(i > nl)
		if (all(converged)) exit
		where (.not. converged)
			rytemp=(xmu+i)*xi2*ry1-rymu
			rymu=ry1
			ry1=rytemp
		end where
	end do
	fact=rjmu/rjl
	rj=rjl1*fact
	rjp=rjp1*fact
	ry=rymu
	ryp=xnu*xi*rymu-ry1
	END SUBROUTINE bessjy_xgexmin
!BL
	FUNCTION absc(z)
	use nrtype
	IMPLICIT NONE
	COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: z
	REAL(DP), DIMENSION(size(z)) :: absc
	absc=abs(real(z))+abs(aimag(z))
	END FUNCTION absc
!BL
	FUNCTION unpacked(mask,vtrue,vfalse)
	use nrtype
	IMPLICIT NONE
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
	REAL(DP), DIMENSION(:), INTENT(IN) :: vtrue, vfalse
	REAL(DP), DIMENSION(size(mask)) :: unpacked
	unpacked=unpack(vtrue,converged,0.0_dp)
	unpacked=unpack(vfalse,.not. converged,unpacked)
	END FUNCTION unpacked
	END SUBROUTINE bessjy_v


end module bessel0

	SUBROUTINE beschb_s(x,gam1,gam2,gampl,gammi)
	USE nrtype
	USE nr, ONLY : chebev
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	REAL(DP), INTENT(OUT) :: gam1,gam2,gampl,gammi
	INTEGER(I4B), PARAMETER :: NUSE1=5,NUSE2=5
	REAL(SP) :: xx
	REAL(SP), DIMENSION(7) :: c1=(/-1.142022680371168_sp,&
		6.5165112670737e-3_sp,3.087090173086e-4_sp,-3.4706269649e-6_sp,&
		6.9437664e-9_sp,3.67795e-11_sp,-1.356e-13_sp/)
	REAL(SP), DIMENSION(8) :: c2=(/1.843740587300905_sp,&
		-7.68528408447867e-2_sp,1.2719271366546e-3_sp,&
		-4.9717367042e-6_sp, -3.31261198e-8_sp,2.423096e-10_sp,&
		-1.702e-13_sp,-1.49e-15_sp/)
	xx=8.0_dp*x*x-1.0_dp
	gam1=chebev(-1.0_sp,1.0_sp,c1(1:NUSE1),xx)
	gam2=chebev(-1.0_sp,1.0_sp,c2(1:NUSE2),xx)
	gampl=gam2-x*gam1
	gammi=gam2+x*gam1
	END SUBROUTINE 


	SUBROUTINE beschb_v(x,gam1,gam2,gampl,gammi)
	USE nrtype
	USE nr, ONLY : chebev
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(:), INTENT(OUT) :: gam1,gam2,gampl,gammi
	INTEGER(I4B), PARAMETER :: NUSE1=5,NUSE2=5
	REAL(SP), DIMENSION(size(x)) :: xx
	REAL(SP), DIMENSION(7) :: c1=(/-1.142022680371168_sp,&
		6.5165112670737e-3_sp,3.087090173086e-4_sp,-3.4706269649e-6_sp,&
		6.9437664e-9_sp,3.67795e-11_sp,-1.356e-13_sp/)
	REAL(SP), DIMENSION(8) :: c2=(/1.843740587300905_sp,&
		-7.68528408447867e-2_sp,1.2719271366546e-3_sp,&
		-4.9717367042e-6_sp, -3.31261198e-8_sp,2.423096e-10_sp,&
		-1.702e-13_sp,-1.49e-15_sp/)
	xx=8.0_dp*x*x-1.0_dp
	gam1=chebev(-1.0_sp,1.0_sp,c1(1:NUSE1),xx)
	gam2=chebev(-1.0_sp,1.0_sp,c2(1:NUSE2),xx)
	gampl=gam2-x*gam1
	gammi=gam2+x*gam1
	END SUBROUTINE 
	
    FUNCTION chebev_s(a,b,c,x)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,b,x
	REAL(DP), DIMENSION(:), INTENT(IN) :: c
	REAL(SP) :: chebev_s
	INTEGER(I4B) :: j,m
	REAL(DP) :: d,dd,sv,y,y2
	if ((x-a)*(x-b) > 0.0) call nrerror('x not in range in chebev_s')
	m=size(c)
	d=0.0
	dd=0.0
	y=(2.0_sp*x-a-b)/(b-a)
	y2=2.0_sp*y
	do j=m,2,-1
		sv=d
		d=y2*d-dd+c(j)
		dd=sv
	end do
	chebev_s=y*d-dd+0.5_sp*c(1)
	END FUNCTION 


	FUNCTION chebev_v(a,b,c,x)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,b
	REAL(DP), DIMENSION(:), INTENT(IN) :: c,x
	REAL(DP), DIMENSION(size(x)) :: chebev_v
	INTEGER(I4B) :: j,m
	REAL(DP), DIMENSION(size(x)) :: d,dd,sv,y,y2
	if (any((x-a)*(x-b) > 0.0)) call nrerror('x not in range in chebev_v')
	m=size(c)
	d=0.0
	dd=0.0
	y=(2.0_sp*x-a-b)/(b-a)
	y2=2.0_sp*y
	do j=m,2,-1
		sv=d
		d=y2*d-dd+c(j)
		dd=sv
	end do
	chebev_v=y*d-dd+0.5_sp*c(1)
	END FUNCTION 

