program ml_ci
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  integer :: i

  read_wf = .True.
  TOUCH read_wf

  open(unit=10, file='n_det', FORM='FORMATTED')
  write(10,*) n_det
  close(10)

  open(unit=10, file='coefs', FORM='FORMATTED')
  do i=1,n_det
    write(10,*) psi_coef(i,1)
  end do
  close(10)

  open(unit=10, file='dets', FORM='FORMATTED')
  do i=1,n_det
    write(10,'(100(I20,X))') psi_det(:,:,i)
  end do
  close(10)
end
