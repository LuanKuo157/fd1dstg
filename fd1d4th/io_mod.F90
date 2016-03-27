Module io_mod
PUBLIC read_conf,read_media, stf_ricker,write_netcdf


contains
subroutine get_conf_name(fnm_input)
character (*),intent(out) :: fnm_input

#ifdef GETARG
integer numarg,i,iargc
   numarg = iargc( )
   do i=1,numarg
      call getarg(i,fnm_input)
      !print *, "the ",i," input arg is ", trim(fnm_input)
      print *, "the conf file name is "//trim(fnm_input)
   enddo
   if ( numarg/=1 ) then
      print *, "numarg should be equal to 1, but now it is ", numarg
      stop 1
   endif
#else
   fnm_input="fd1d.conf"
#endif
end subroutine get_conf_name

subroutine read_conf(fnm_conf,nx,nt,h,dt,t0,fc,sid,recvid)
integer:: nx,nt,sid,recvid
real:: t0,fc,h,dt
character *80 fnm_conf,str1,str2
open(10,file=fnm_conf,status='old')
read(10,*),str1, str2, nx
read(10,*),str1, str2, nt
read(10,*),str1, str2, h
read(10,*),str1, str2, dt
read(10,*)
read(10,*),str1, str2, t0
read(10,*),str1, str2, fc
read(10,*),str1, str2, sid
read(10,*)
read(10,*),str1, str2, recvid
close(10)
end subroutine


subroutine read_media( fnm_media, p, c, nx )
integer:: nx,i
real::p(nx),c(nx)
character (*) :: fnm_media
open(20,file=fnm_media,status='old')
do i=1,nx
  read(20,*), p(i), c(i)
enddo
close(20)
end subroutine

subroutine stf_ricker(tau,fc,t0,f)
    real:: tau,t0,fc
    real:: u,f0,f,PI=3.141592653
    if (tau<=0.0) then 
        f=0.0 
    else
        f0=sqrt(PI)/2.0
        u=(tau-t0)*2.0*PI*fc
        f=(u**2.0/4.0-0.5)*exp(-u**2.0/4.0)*f0
    endif
end subroutine

subroutine write_netcdf(outfile,xpos,dat,nt)
use netcdf
implicit none
real, dimension(nt), intent(in) :: xpos
real, dimension(nt), intent(in) :: dat
integer(kind=4):: ncid, x_dimid, x_varid, varid
integer(kind=4),intent(in)::  nt
character(*), intent(in) ::outfile
!Create a netCDF file
print*,"1"
call check( nf90_create( outfile, NF90_CLOBBER, ncid ) )
!Define the dimensions
print*,"2"
call check( nf90_def_dim( ncid, "t",nt, x_dimid ))

!Define the coordinate variables
call check( nf90_def_var( ncid, "t",x_dimid, x_varid ) )

!Difine variable
call check( nf90_def_var( ncid, "Vel", NF90_FLOAT,x_dimid, varid ) )
call check( nf90_enddef(ncid) )

!write data
call check(nf90_put_var(ncid, x_varid, xpos))
call check(nf90_put_var(ncid, varid, dat))
call check(nf90_close(ncid))
end subroutine write_netcdf

subroutine check(status)
use netcdf
implicit none 
integer,intent(in):: status
if (status .ne. nf90_noerr) then 
   write(*,*) trim(adjustl(nf90_strerror(status)))
endif 
end subroutine check

end module


















