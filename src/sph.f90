! calculate the value at any given point using spherical harmonic 
! expansion coefficient 
subroutine spherical_harmonic_expansion(phi,xlm,cilm,psi)
  use const
  implicit none
  real*8, intent(out) :: psi
  real*8, intent(in) :: phi
  real*8, intent(in) :: cilm(:,:,:)
  real*8, intent(in) :: xlm(:)
  integer :: l, m, idx
  integer :: lmax
  !real*8, dimension (2, lmax+1, lmax+1) :: cilm

  write(*,*) phi
  lmax = size(cilm(1,:,1))
  write(*,*) "lmax is",lmax
  write(*,*) size(cilm(1,:,1)),cilm(:,1,1),xlm(1)

  psi = 0.0d0
  do l = 0, lmax
      do m = 0, l
        idx = (l * (l + 1))/ 2 + m + 1
        write(*,*) l,m, xlm(idx)
        psi = psi +  xlm(idx) * (cilm(1,l+1,m+1)*dcos(m*phi)+ cilm(2,l+1,m+1)*dsin(m*phi))
      enddo
  enddo

  end subroutine
