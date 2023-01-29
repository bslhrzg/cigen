subroutine read_namelist(file_path,algo,N,Ne,M_hiddens,tempgen,Nsample,fgen,&
it_up,tol,el,frozen,sym,symg,so_ir_in,ir_t)

        use, intrinsic :: iso_fortran_env, only: stderr => error_unit

        !! Reads Namelist from given file.
        character(len=*),                        intent(in)  :: file_path
        integer,                                 intent(out) :: algo,N,Ne,M_hiddens,Nsample,it_up,el,frozen,ir_t

        integer, dimension(100),                 intent(out) :: so_ir_in
        double precision,                        intent(out) :: tempgen,tol,fgen
        logical,                                 intent(out) :: sym
        character(10),                           intent(out) :: symg
        integer                                              :: fu, rc

        ! Namelist definition.
        namelist / params / algo,N,Ne,M_hiddens,tempgen,Nsample,fgen,it_up,tol,el,frozen,sym,symg,so_ir_in,&
ir_t,itcig,maxiter,diagtol

        ! Check whether file exists.
        inquire (file=file_path, iostat=rc)

        if (rc /= 0) then
            write (stderr, '(3a)') 'Error: input file "', trim(file_path), '"does not exist.'
            return
        end if

        ! Open and read Namelist file.
        open (action='read', file=file_path, iostat=rc, newunit=fu)
        read (nml=params, iostat=rc, unit=fu)

        if (rc /= 0) then
            write (stderr, '(a)') 'Error: invalid Namelist format.'
        end if

        close (fu)
end subroutine read_namelist
