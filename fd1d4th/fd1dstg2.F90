!---------------------------------------------------------------------
!  This program is used for synthetize seismogram in 1D complex media
!        Author: ZHOU Li        Email:zhouli35@mail.ustc.edu.cn 
!        Date:            2014-06-08 Sun 15:46:26 CST
!---------------------------------------------------------------------

program  fd1dstg
implicit none
integer:: nx,nt,recvid,sid  !网格点数，时间步数
!接收点位置,震源位置
integer:: i,j,k
real::h,dt              !网格间距，时间步长
real:: t0,fc,f        !recker 震源参数
real(8)::c0=9.0/8,c1=-1./24             !差分系数
real,dimension(:),allocatable:: p,c,T,V,Vg,Tg,tser,vrec
character *80 fnm_conf, fnm_media, outfile
fnm_conf='fd1d.conf'
fnm_media='media.dat'
write(*,*) '控制文件和介质参数文件',fnm_conf,fnm_media
call read_conf(fnm_conf,nx,nt,h,dt,t0,fc,sid,recvid)

allocate(p(nx),c(nx),T(nx+2),V(0:nx+1),Vg(nx),Tg(nx),tser(nt),vrec(nt) )
print*,nx,nt,h,dt
print*,t0,fc,sid,recvid
call read_media(fnm_media,p,c,nx)
T=0.; V=0.;  !Initial condition
open(30,file='seismo.dat')
open(40,file='source.dat')
print*,'c0=',c0,'c1=',c1
do i=0, nt
  Tg(1) = (c0*T(2)+c1*(T(3)-T(2)) )/h
  v(0)=-v(1)
  T(nx+1)=-p(nx)*c(nx)*v(nx)
  v(nx+1)=-T(nx+1)/(p(nx)*c(nx))
  T(nx+2)=-p(nx)*c(nx)*v(nx+1)
  do  j= 2, nx
       vg(j) = ( c0* ( v(j)-v(j-1) ) + c1*( v(j+1)-v(j-2) ) )/h 
       Tg(j) = ( c0* ( T(j+1)-T(j) ) + c1*( T(j+2)-T(j-1) ) )/h
  !     write(40,*),j,vg(j),tg(j)
       T(j)=T(j)+p(j)*c(j)**2*vg(j)*dt
       V(j)=v(j)+1./p(j)*Tg(j)*dt
  enddo
      T(1)=0;V(1)=V(1)+1./p(1)*Tg(1)*dt
      call stf_ricker(i*dt,fc,t0,f)
      v(sid)=v(sid)+f/p(sid)*dt
  write(30,*),i*dt,v(recvid)
  write(40,*),i*dt,f
      vrec(i)=v(recvid)
 enddo
 close(30)
 close(40)
 outfile='seismo.nc'
 do i=1,nt
    tser(i)=i*dt
 end do
 call write_netcdf(outfile,tser,vrec,nt)

end

subroutine read_conf(fnm_conf,nx,nt,h,dt,t0,fc,sid,recvid)
integer:: nx,nt,sid,recvid
real:: to,fc,h,dt
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
call check( nf90_def_dim( ncid, "time",nt, x_dimid ))

!Define the coordinate variables
call check( nf90_def_var( ncid, "time",NF90_REAL,x_dimid, x_varid ) )

!Difine variable
call check( nf90_def_var( ncid, "Velocity field", NF90_FLOAT,x_dimid, varid ) )
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



