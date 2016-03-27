subroutine read_media( fnm_media, p, c, nx )
integer:: nx,i
real::p(nx),c(nx)
open(20,file=fnm_media,status='old')
do i=1,nx
  read(20,*), p(i), c(i)
  enddo
  close(20)
end subroutine
