program CI


      ! MAIN PROGRAM

      use var
      use tools
      use omp_lib
      use rbm
      use iomod

      implicit none

      integer :: i,j,k,cnt,n_dets,ndn,Ngen,Nsamp
      integer,dimension(:),allocatable ::dets,Dnew,dets_n
      integer,dimension(:,:),allocatable :: dets_bin
      double precision :: eps=0.01,lreg=0.
      double precision, dimension(:),allocatable :: eigenvectors
      double precision, dimension(:),allocatable :: a,b
      double precision, dimension(:,:),allocatable :: W

      ! initialize some variables
      call init_var()

      ! RBM parameters allocation
      allocate(a(N))
      allocate(b(M_hiddens))
      allocate(W(N,M_hiddens))
 
      ! Read current wavefunction
      n_dets=0
      call countlines_f('dets',1,n_dets)      
      allocate(dets(n_dets))
      open(unit=13, file='dets')
      do i=1,n_dets
              read(13,*) dets(i)
      end do
      close(13)

      ! Read training set of determinants
      Nsamp=0
      call countlines_f('dets_sample',1,Nsamp)        
      allocate(dets_bin(Nsamp,N))
      open(unit=13, file='dets_sample')
      do i=1,Nsamp
              read(13,*) (dets_bin(i,j),j=1,N)
      end do
      close(13)


      
      ndn=0

      ! This is not used at the moment **********************************************
      ! I used to do the pruning and generation of training set 
      ! here but now using python to convert QP output and prune
      ! and create trainig set at the same time
      !allocate(eigenvectors(n_dets))
      
      !open(unit=90,file='coefs')
      !read(90,*) eigenvectors
      !close(90)
      
     
      !sort coefs by probability 
      !toto=rargsort(eigenvectors(:,1)*eigenvectors(:,1))
      !toto = reverse(toto)
      
      !eigenvectors(:,1)=eigenvectors(toto,1)
      !dets(:,:)=dets(toto,:)


      !prunning
      !print *, "prunning with tol = ",tol
      !c_p=count_prun(eigenvectors(:,1)*eigenvectors(:,1),tol)

      !allocate(cc_pruned(c_p))
      !allocate(Dp_pruned(c_p))
      
      !call pruning(eigenvectors(:,1)*eigenvectors(:,1),dets,c_p,cc_pruned,Dp_pruned)
      ! *****************************************************************************

      print *, '***************************************************************************'
      print *,''

      print *,"   ____ ___                 "
      print *,"  / ___|_ _|__ _  ___ _ __  "
      print *," | |    | |/ _` |/ _ \ '_ \ "
      print *," | |___ | | (_| |  __/ | | |"
      print *,"  \____|___\__, |\___|_| |_|"
      print *,"           |___/            "

      ! RBM Initialization
      call rand_init(a,b,W)
       
      !RBM update
      print *, ''
      print *, ' ----------------------------------------------------------------------------'
      print *,  "RBM Training:"
      print *,  "Training set size             = ", Nsamp
      print *,  "Number of hidden cells        = ", M_hiddens
      print *,  "Batch size                    = 10"
      print *,  "Learning rate                 = ", eps
      print *,  "Number of learning iterations = ", it_up
      do k=1,it_up
        write(*,*) 'learning iteration : ', k ,  'is complete'
        call update_rbm(a,b,W,dets_bin,10,Nsamp,eps,lreg)
      enddo
      
      print *, ''
      print *, ' ----------------------------------------------------------------------------'
      print *, 'Generative procedure with temperature                     = ',tempgen
      if (sym) then
        print *, 'Spatial symmetry enabled in generation with point group : ',symg 
        print *, 'Irreducible representation of the spin-orbitals         : '
        print *, so_ir 
      endif

      print *, 'fgen                                                      = ',fgen
      print *, 'Number of frozen spin-orbitals in the generation          = ',frozen

      Ngen=int(fgen*n_dets)
      allocate(Dnew(Ngen-1))
      Dnew=0

      if (algo .eq. 1) then
        print *, 'Algo : Direct generation of determinants'
        call generate(Dnew,W,a,b,dets,Ngen,Ngen*1000,cnt)
      else if (algo .eq. 2) then 
        print *, 'Algo : Single/Double substitution generation of determinants'
        call generate_SD(Dnew,W,a,b,dets,Ngen,Ngen*1000,cnt)
      else 
        print *, 'Please enter a correct determinant generation algorithm number:'
        print *, '1: Direct generation'
        print *, '2: Single/Double generations'
      endif


      if (cnt .le. Ngen) Dnew=remove_zeros(Dnew)
      ndn=size(Dnew)
      print *, 'Number of determinants generated : ', ndn

      allocate(dets_n(n_dets+ndn))
      dets_n(:n_dets)=dets
      dets_n(n_dets+1:)=Dnew
      
      dets=dets_n
      deallocate(dets_n)
      n_dets=n_dets+ndn    
      
      
       open(unit=100,file='dets')
         do i=1,n_dets
                write(100,*) dets(i)
        end do
       close(100)

      print *, '***************************************************************************' 

end program CI
