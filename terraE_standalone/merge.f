program annual_ij
real*4 ij_data(72,46,242)
character*80 TITLE
OBS(:,:)=1.E-30
open(10,file='ANN1951-1955.ijE8piM20',STATUS='UNKNOWN',&
  form='unformatted')
read(10) TITLE,ij_data
close(10)

!do i=1,72
! do j=1,46
!   if(MASK(i,j).gt.-99.)then
!     OBS(i,j)=O(i,j)
!   else
!     OBS(i,j)=L(i,j)
!   endif
! enddo
!enddo

end program annual_ij
