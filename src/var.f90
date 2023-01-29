module var

        use tools
        use omp_lib
        use iomod

        implicit none

        integer :: algo,N,Ne,M_hiddens=0,frozen=2,ir_t=0, &
Nsample=25,it_up=10,el=0
        double precision :: fgen,tempgen=1,tol
        character(10) :: symg ! symmetry point group
        logical :: sym

        integer,dimension(100) :: so_ir_in
        integer, dimension(:), allocatable :: so_ir
        integer, dimension(8,8) :: D2h
        integer, dimension(4,4) :: C2v
        integer, dimension(:,:), allocatable :: M_sym

        contains 

subroutine init_var()

        implicit none

        integer :: i,j,k,l,cnt=0

        D2h(1,:)=[1,2,3,4,5,6,7,8]
        D2h(2,:)=[2,1,4,3,6,5,8,7]
        D2h(3,:)=[3,4,1,2,7,8,5,6]
        D2h(4,:)=[4,3,2,1,8,7,6,5]
        D2h(5,:)=[5,6,7,8,1,2,3,4]
        D2h(6,:)=[6,5,8,7,2,1,4,3]
        D2h(7,:)=[7,8,5,6,3,4,1,2]
        D2h(8,:)=[8,7,6,5,4,3,2,1]

        C2v(1,:)=[1,2,3,4]
        C2v(2,:)=[2,1,4,3]
        C2v(3,:)=[3,4,1,2]
        C2v(4,:)=[4,3,2,1]

        call read_namelist('input.nml',algo,N,Ne,M_hiddens,tempgen,Nsample,fgen,it_up,tol, &
el,frozen,sym,symg,so_ir_in,ir_t)


!       print *, "Input parameters:"
!       print *, "M_hiddens = ",M_hiddens
!       print *, "tempgen = ", tempgen
!       print *, "Nsample = ", Nsample
!       print *, "fgen = ", fgen 
!       print *, "it_up = ", it_up
!       print *, "tol = ", tol
!       print *, "el = ", el
!       print *, "frozen = ", frozen
!       print *, "sym = ", sym
!       print *, "symg = ", symg
!       print *, "so_ir_in = ", so_ir_in
!       print *, "ir_t = ", ir_t
       
        allocate(so_ir(N))
        so_ir=so_ir_in(:N)
        
        if (symg == 'D2h') then
                allocate(M_sym(size(D2h,1),size(D2h,2)))
                M_sym=D2h
        endif
        if (symg == 'C2v') then
                allocate(M_sym(size(C2v,1),size(C2v,2)))
                M_sym=C2v
        endif

  end subroutine         
end module
