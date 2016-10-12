!   ------------------------
!   09/23/2015 chukren
!   Fix the bug in surface_laplacian_theta()
 
  module const

  character(len=80), parameter :: sphcoef = "MarsTopo2600.shape.sh"
  real*8,  parameter :: pi    = 3.1415926535898d0
  real*8,  parameter :: twopi = 2.0 * pi
  real*8,  parameter :: d2r   = pi/180.d0
  real*8,  parameter :: r2d   = 180.d0/pi
  real*8,  parameter :: R_MARS = 3389500.d0 
  real*8,  parameter :: ell = 1./169.8d0
  !===

  !=== 1 minite resolution ===
  !real*8,  parameter :: step = 1.0d0/60.0   ! 1 minute
  !integer, parameter :: nlat = 10800        ! 180 * 60
  !integer, parameter :: nlon = 21600        ! 360 * 60
  !=== 2 minite resolution ===
  !real*8,  parameter :: step = 1.0d0/30.0   ! 2 minute
  !integer, parameter :: nlat = 5400         ! 180 
  !integer, parameter :: nlon = 10800        ! 360 
  !=== 6 minite resolution ===
  real*8,  parameter :: step = 1.0d0/6.0     ! 6 minute
  integer, parameter :: nlat = 1080          ! 180 * 60
  integer, parameter :: nlon = 2160          ! 360 * 60
  !=== half degree resolution ===
  !real*8,  parameter :: step = 1.0d0/2.0    ! 0.5 degree
  !integer, parameter :: nlat = 360          ! 180 
  !integer, parameter :: nlon = nlat*2       ! 360 
  !=== 1 degree resolution ===
  !real*8,  parameter :: step = 1.0d0        ! 1 degree
  !integer, parameter :: nlat = 180          ! 180 
  !integer, parameter :: nlon = 360          ! 360 


  end module

  program readtopo
  use const
  use SHTOOLS, only : SHRead, PlmBar, MakeGridDH
  implicit none      

  !=== spherical degree ===
!  integer, parameter :: lmax  = 2600        ! degree of sph harmonics
!  integer, parameter :: paramax = 3383901   ! number of alm, blm
  integer, parameter :: lmax  = 539       ! degree of sph harmonics
!  integer, parameter :: paramax =     ! number of alm, blm
!  integer, parameter :: lmax  = 1000        ! degree of sph harmonics
!  integer, parameter :: paramax = 501501    ! number of alm, blm
!  integer, parameter :: lmax  = 20          ! degree of sph harmonics 20 
!  integer, parameter :: paramax = 231       ! number of alm, blm 231
!  integer, parameter :: lmax  = 40          ! degree of sph harmonics 20 
!  integer, parameter :: paramax = 861       ! number of alm, blm 231
!  integer, parameter :: lmax  = 120         ! degree of sph harmonics 20 
!  integer, parameter :: paramax = 7381      ! number of alm, blm 231
!  integer, parameter :: lmax  = 160          ! degree of sph harmonics 20 
!  integer, parameter :: paramax = 13041      ! number of alm, blm 231
!  integer, parameter :: lmax  = 200         ! degree of sph harmonics 20 
!  integer, parameter :: paramax = 20301     ! number of alm, blm 231

  !===
  character (len = 70) :: outfile, file_map, ell_map

  real*8 :: P20,areoid_factor
  real*8 :: R_areoid, R_ell 
  real*8 :: elev,elev_km 
  real*8 :: alat, alon, costheta    
  real*8 :: lat, colat, lon  ! station and source 
  real*8 :: inclon, inclat

  real :: start,finish ! cpu time consumed 
  integer :: lmax_use ! dummy para 
  integer :: i, j, ielev
  logical :: isxlm
  integer, parameter :: gridmax = nlat*nlon    ! 
  integer :: ngrid, norm, sampling, csphase, lmax_calc 
  integer :: l,m,idx
  real*8 :: phi

  !real*8, dimension (2, lmax+1, lmax+1) :: topo_cilm
  !real*8, dimension (paramax) :: xlm 
  !real*8, dimension (nlat,nlon) :: topo

  real*8, allocatable :: topo_cilm(:,:,:)
  real*8, allocatable :: topo(:,:)

  call cpu_time(start)
  
  allocate(topo_cilm(2, lmax+1, lmax+1))
  allocate(topo(nlat,nlon))


  call SHRead(sphcoef, topo_cilm, lmax_use)
  write(*,*)lmax,lmax_use,topo_cilm(:,lmax+1,lmax+1)

  if(lmax_use .ne. lmax) write(*,*) 'check the lmax!' 

  !---------- creat output file ------------       
  outfile = 'topo.10d540.Mars.dat'
  open(100,file=outfile)

  file_map = 'topo.10d540.Mars.xyz'
  open(101,file=file_map)

  ell_map = 'areoid.Mars.xyz'
  open(102,file=ell_map)

  inclon = 360.d0/nlon
  inclat = 180.d0/nlat

  write(*,*) 'check size of cilm before calling MakeGridDH'
  write(*,*) size(topo_cilm(:,1,1)), size(topo_cilm(1,:,1)),size(topo_cilm(1,1,:))
  call MakeGridDH(topo, ngrid, topo_cilm, lmax, norm=1,sampling=2, csphase=1, lmax_calc=lmax)
  write(*,*) size(topo(:,1)), size(topo(1,:))
  if(gridmax.ne.ngrid) write(*,*) 'check the grid size!'

  do i = 1, nlat

    colat = i*inclat*d2r
    lat = 90.0 - i*inclat
   
    costheta = dcos(colat)
    P20 = 0.5d0*(3.0d0*costheta*costheta-1.0d0) 
    areoid_factor = 1.0d0 - (2.0d0/3.0d0)*ell*P20
    
    do j = 1, nlon
      lon = (j-1)*inclon
      elev = topo(i,j)
      ! geographic to geocentric     
      R_areoid = R_MARS*areoid_factor
   
      ielev = IDNINT(elev-R_areoid)
      write(100,1991) ielev

      elev_km = (elev-R_areoid) / 1000.d0   
      write(101,1990) lon,lat,elev_km

      R_ell = (R_areoid - R_MARS) / 1000.d0
      write(102,1990) lon,lat,R_ell

    enddo
  enddo

1990 format(F9.2,F9.2,F8.2)
1991 format(I10)

  deallocate(topo_cilm)
  deallocate(topo)

!  stop

  !!---------- loop over stations ----------  
!  do i = 1, nlat
!     colat = i*step
!     lat   = 90 - colat
!     alat = colat * d2r
!     costheta = dcos(alat)
     
!     isxlm = .true.
!     do j = 1, nlon
!       lon  = j*step
!       alon =  lon * d2r

!       if(alon > twopi)      alon = alon - twopi
!       if(alon < (-1)*twopi) alon = alon + twopi
       
!       if(isxlm) then 
!         ! call PlmBar(xlm, lmax, costheta)
!          isxlm = .false.
!       end if

!       phi = alon

!       elev = 0.0d0
!       do l = 0, lmax
!         do m = 0, l
!             idx = (l * (l + 1))/ 2 + m + 1
!             elev = elev +  xlm(idx) * (cilm(1,l+1,m+1)*dcos(m*phi)+ cilm(2,l+1,m+1)*dsin(m*phi))
!         enddo
!       enddo      

!     end do
!  enddo 
!==== relic code ===


  close(100)
  close(101)
  write(*,*)'output ', outfile
  write(*,*)'map_file ', file_map

  call cpu_time(finish)
  print '("time = ",f10.3," minutes.")',(finish-start)/60.0

  end program 

