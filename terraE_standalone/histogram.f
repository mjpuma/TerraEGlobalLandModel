      SUBROUTINE Histogram(X,frac,N,sum_frac)
!------------------------------------------------------------
! This subroutine takes a real array X(N) and forms a
! plain text 20-bin histogram of the distribution of values.
! Each bin entry output represents SCALE input entries.
!------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL*4, INTENT(IN) :: X(N),frac(N)
      INTEGER, PARAMETER :: BINS=20
      REAL*4 :: min_x, max_x, BINWIDTH,sum_frac
      real*4,dimension(N):: pdf
      real*4,dimension(BINS)::H_weight
      integer,dimension(BINS)::H
      integer :: BIN, I, J

      H=0
      H_weight=0
      min_x=MINVAL(X)
      max_x=MAXVAL(X)

      BINWIDTH=(max_x-min_x)/REAL(BINS)
      PRINT *, "Number of values = ", N
      PRINT *, "Minimum value = ", min_x
      PRINT *, "Maximum value = ", max_x

      DO I=1,N
         BIN=INT(1+(X(I)-min_x)/BINWIDTH)
         IF ( BIN < 1    ) BIN=1    ! Check for underflows
         IF ( BIN > BINS ) BIN=BINS ! Check for overflows
         H(BIN)=H(BIN)+1
         H_weight(BIN) = H_weight(BIN)+frac(I)
      END DO

      DO I=1,BINS
        pdf(I) = H_weight(I)/sum_frac/BINWIDTH
      END DO

      END SUBROUTINE Histogram
