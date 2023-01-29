module iomod
 
 implicit none
 contains

subroutine countlines_f(filename,chan,nlines)

 character(*),intent(in) :: filename
 integer,intent(in)      :: chan
 integer,intent(inout)   :: nlines
 integer                 :: io

 open(chan,file=filename)
 do
    read(chan,*,iostat=io)
    if (io/=0) exit
    nlines = nlines + 1
  enddo
  close(chan)
end subroutine

subroutine countlines_u(filename,chan,nlines)

 character(*),intent(in) :: filename
 integer,intent(in)      :: chan
 integer,intent(inout)   :: nlines
 integer                 :: io

 open(chan,file=filename,form='unformatted')
 do
    read(chan,*,iostat=io)
    if (io/=0) exit
    nlines = nlines + 1
  enddo
  close(chan)
end subroutine
end module
