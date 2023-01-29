program ml_ci
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  integer :: i

  read_wf = .False.
  TOUCH read_wf

  open(unit=10, file='n_det', FORM='FORMATTED')
  read(10,*) n_det
  close(10)
  TOUCH n_det
  print *, 'Ndet= ', n_det

  open(unit=10, file='coefs', FORM='FORMATTED')
  do i=1,n_det
    read(10,*) psi_coef(i,1)
!    print *, psi_coef(i,1)
  end do
  close(10)

  open(unit=10, file='dets', FORM='FORMATTED')
  do i=1,n_det
    read(10,*) psi_det(:,:,i)
  end do
  close(10)

  TOUCH psi_det psi_coef
  call save_wavefunction

end
