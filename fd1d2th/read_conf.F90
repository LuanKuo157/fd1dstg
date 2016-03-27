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
