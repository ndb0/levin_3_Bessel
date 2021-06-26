program main

use nrtype
use bessel0
use solve
implicit none

real(dp) a,b,r,d,res1,res2,res,x,aterm,u1term,xarg


integer, parameter :: nn7=80


integer jin,n,n2,i,l,j

character(200) filename,func

real(dp) point(nn7),rhs(nn7),mat(nn7,nn7),rhs1(nn7)

integer indx(nn7),ni,nj,im,jm,norder,iv
real(dp) ddd

real(dp) fdata(6),in1,in2,in3,mm,nn,pp,ind,beta

real(dp) rj,ry,rjp,ryp
	
norder =8

beta=.4
mm=1.0
nn= 8.661690959*1e-5!(8.661e-5)*Sqrt(0.3333)
pp=10.0!1.0*Sqrt(0.3333)
in1=beta
in2=beta+2.0
in3=beta+2.0
ind=1.0-beta
	
	fdata(1)=in1
	fdata(2)=in2
	fdata(3)=in3
	fdata(4)=mm
	fdata(5)=nn
	fdata(6)=pp
	

a=1e-2
b=100.0


!call bessjy_s(b,in1,rj,ry,rjp,ryp)

print *, "r",rj,ry
print *, beslj(0.2_dp,3e4_dp),besly(0.2_dp,3e8_dp)

!##################################
!  go to 552
  
 !#################################### 
n=nn7/norder

d=(a+b)/2.0+1e-15

print *, "d",d




 do i=1,n
 point(i)=a+(i-1)*(b-a)/(n-1)
! print *,"point(",i,")=",point(i)
 enddo
 
 
do iv=1,norder   
do i=1,n
 x=point(i)
 rhs(i+(iv-1)*n)=ffunc(iv,x,ind)
 !print *,"rhs(",i,")=",rhs(i)
 enddo    
enddo   
   
!print *, rhs
 
 

  do ni=1,norder
  do nj=1,norder
  do i=1,n
  do j=1,n
  
  x=point(i)
 aterm=amat(x,fdata,nj,ni)
  !aterm=amat22(x,r,nj,ni)
  if(ni.eq.nj)then
  u1term=u1(x,j,d)
  else
  u1term=0.0
  endif
  
  im=i+(ni-1)*n
  jm=j+(nj-1)*n
  mat(im,jm)=u(x,j,d)*aterm+u1term
 
 end do 
 end do
 end do
 end do

  

 	! do i=1,n2
    !    Write(*,"(10F10.4)")(mat(i,l),l=1,n2)

     !end do
 

 
 ! call svdcmp_dp(mat,ww,vV)
 call ludcmp(mat,indx,ddd)
 call lubksb(mat,indx,rhs)
  !call solveqn(mat,rhs,rhs1)
  
!call DGESV(m, 1, real(mat,8), m, pivot,rhs1 , m, ok)
 !call dgels('N', n2,n2, 1, real(mat,8),n2, rhs,n2, pivot, 200, ok)
 !   call gaussj(mat,b11)
    
       
    
 !   do i=1, nn7
!	   write(*,*) rhs(i)
!	end do
	
	! do i=1,n2
    !    Write(*,"(10F10.4)")(mat(i,l),l=1,n2)

     !end do
	print *,"######################"
	

    res1=0.0
    
    do i=1,n
    do iv=1,norder
        xarg=b
        !print *, iv
        res1=res1+rhs(i+(iv-1)*n)*u(b,i,d)*wfunc(iv,xarg,fdata)
        print *, beslj(in1,mm*xarg),beslj(in2,nn*xarg),beslj(in3,pp*xarg)
         print *, xarg,wfunc(iv,xarg,fdata)
     end do
   !  print *, res1!,rhs(i),u(b,i,d),j0(r*b)
     end do
     print *,"######################"
    res2=0.0
    
    do i=1,n
    do iv=1,norder
        xarg=a
        !print *, iv
        res2=res2+rhs(i+(iv-1)*n)*u(a,i,d)*wfunc(iv,xarg,fdata)
       ! print *, res2,wfunc(iv,xarg,fdata)
   !     print *, beslj(in1,a),beslj(in2,nn*a),beslj(in3,pp*a)
    !    print *, xarg,wfunc(iv,xarg,fdata)
     end do
    ! print *, res2!,rhs(i),u(b,i,d),j0(r*b)
     end do

     
     res=res1-res2
    print *, "final result=",res
   
    !################
    res1=0.0
    
    do i=1,n
    do iv=1,norder
        xarg=b
        !print *, iv
        res1=res1+rhs(i+(iv-1)*n)*u(b,i,d)*wfuncY(iv,xarg,fdata)
      !  print *, res1,beslj(in1,b),beslj(in2,nn*b),beslj(in3,pp*b)
     end do
   !  print *, res1!,rhs(i),u(b,i,d),j0(r*b)
     end do
    ! print *,"######################"
    res2=0.0
    
    do i=1,n
    do iv=1,norder
        xarg=a
        !print *, iv
        res2=res2+rhs(i+(iv-1)*n)*u(a,i,d)*wfuncY(iv,xarg,fdata)
       ! print *, res2,wfunc(iv,xarg,fdata)
     end do
    ! print *, res2!,rhs(i),u(b,i,d),j0(r*b)
     end do

     
     res=res1-res2
   552 print *, "final result=",res
    !################
	
end program main





 
 
