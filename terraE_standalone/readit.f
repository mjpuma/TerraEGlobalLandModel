c ifort readit.f -o readit -I/usr/local/other/netcdf/3.6.1/include
c    -L/usr/local/other/netcdf/3.6.1/lib -lnetcdf
      program readit
      implicit none
      include 'netcdf.inc'
      character*80 :: infile,varname
      integer :: fid,varid,status
      integer, parameter :: x=360,y=180,land=15238,ntime=248
      real*4, dimension(x,y) :: arr2d
      integer, dimension(land) :: xy
      real*4, dimension(land,ntime) :: arrpack
c open file
      infile='Tair_cru198207.nc'
      status = nf_open(trim(infile),nf_nowrite,fid)
c read in x y indices
      varname='land'
      status = nf_inq_varid(fid,trim(varname),varid)
      status = nf_get_var_int(fid,varid,xy)
c read in data over land points only
      varname='Tair'
      status = nf_inq_varid(fid,trim(varname),varid)
      status = nf_get_var_real(fid,varid,arrpack)
      status = nf_close(fid)
c spread data to full x y grid

      end program readit
