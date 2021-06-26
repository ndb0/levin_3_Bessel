module solve
use nrtype
use bessel0
USE nrutil
implicit none

contains
	
		
		
SUBROUTINE ludcmp(a,indx,d)
	use nrtype
	USE bessel0; USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
	IMPLICIT NONE
	REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
	REAL(dp), INTENT(OUT) :: d
	REAL(dp), DIMENSION(size(a,1)) :: vv
	REAL(dp), PARAMETER :: TINY=1.0e-20_sp
	INTEGER(I4B) :: j,n,imax
	n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
	d=1.0
	vv=maxval(abs(a),dim=2)
	if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
	vv=1.0_sp/vv
	do j=1,n
		imax=(j-1)+imaxloc_r(vv(j:n)*abs(a(j:n,j)))
		if (j /= imax) then
			call swap_rv(a(imax,:),a(j,:))
			d=-d
			vv(imax)=vv(j)
		end if
		indx(j)=imax
		if (a(j,j) == 0.0) a(j,j)=TINY
		a(j+1:n,j)=a(j+1:n,j)/a(j,j)
		a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
	end do
	END SUBROUTINE ludcmp
	
		SUBROUTINE lubksb(a,indx,b)
		use nrtype
	USE bessel0; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(dp), DIMENSION(:,:), INTENT(IN) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: i,n,ii,ll
	REAL(dp) :: summ
	n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
	ii=0
	do i=1,n
		ll=indx(i)
		summ=b(ll)
		b(ll)=b(i)
		if (ii /= 0) then
			summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
		else if (summ /= 0.0) then
			ii=i
		end if
		b(i)=summ
	end do
	do i=n,1,-1
		b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
	end do
	END SUBROUTINE lubksb
end module solve 
