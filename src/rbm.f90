module rbm
        use omp_lib
        use tools
        use var
        implicit none
contains

function normal(mu,s) result(r)

        ! sample number from normal distribution 
        ! with mean mu and std s
        
        implicit none
        
        double precision,intent(in) :: mu,s
        double precision :: r,u,v
        double precision :: pi=3.14159265359
        
        call random_number(u)
        call random_number(v)
        
        r=sqrt(-2*log(u))*cos(2*pi*v)*s+mu
        
        return
end function
        

subroutine rand_init(a,b,W)

        ! initialize the RBM parameters from a normal distribution
        ! N(mu=0,s=0.05)

        implicit none
       
        double precision, intent(inout), dimension(:,:) :: W
        double precision, intent(inout), dimension(:) :: a,b
        integer :: i,j
        double precision :: mu=0.,s=0.05
        
        
        do i=1,size(a)
          a(i)=normal(mu,s)
        end do
        do j=1,size(b)
          b(j)=normal(mu,s)
        end do
        do i=1,size(a)
         do j=1,size(b)
           W(i,j)=normal(mu,s)
         end do
        end do
         
end subroutine

function sigmoid(x) result(y)

        implicit none

        double precision, intent(in) :: x
        double precision :: y

        y=1./(1.+exp(-x))

        return
end function sigmoid

function outer(x,y) result(O)

        ! for double float vectors x,y
        ! return the outer product of vector x and y

        implicit none

        double precision,intent(in), dimension(:) :: x,y
        double precision,dimension(size(x),size(y)) :: O
        integer :: i,j

        O=0
        do j=1,size(y)
                do i=1,size(x)
                        O(i,j)=x(i)*y(j)
                end do
        end do 
        return
end function outer

function outer_int(x,y) result(O)

        ! for integer vectors x,y
        ! return the outer product of vector x and y

        implicit none

        integer,intent(in), dimension(:) :: x,y
        integer,dimension(size(x),size(y)) :: O
        integer :: i,j

        O=0
        do j=1,size(y)
                do i=1,size(x)
                        O(i,j)=x(i)*y(j)
                end do
        end do 
        return
end function outer_int

function SampleHidden(v,W,b,temp) result(h)

        ! sample the hidden layer of the RBM

        implicit none

        double precision, intent(in), dimension(:,:) :: W
        double precision, intent(in), dimension(:) :: b
        integer, intent(in), dimension(:) :: v
        integer, dimension(size(b)) :: h
        real,dimension(size(b)) :: htmp
        double precision, dimension(size(b)) :: mh
        double precision, intent(in) :: temp
        integer :: n,i

        n=size(b)

        mh=0
        mh=matmul(transpose(W),v)
        do i=1,n
                mh(i)=sigmoid((mh(i)+b(i))/temp)
        end do
        h=0
        call random_number(htmp)
        do i=1,n
                if (htmp(i) < mh(i)) then
                        h(i)=1
                endif
        end do

        return
end function SampleHidden

function SampleHidden_p(W,b,temp) result(h)

        ! sample the hidden layer of the RBM
        ! from a uniformly distributed visible vector
        ! v with real values in [0,1]

        implicit none

        double precision, intent(in), dimension(:,:) :: W
        double precision, intent(in), dimension(:) :: b
        double precision, dimension(size(W,1)) :: v
        integer, dimension(size(b)) :: h
        real,dimension(size(b)) :: htmp
        double precision, dimension(size(b)) :: mh
        double precision, intent(in) :: temp
        integer :: n,i

        n=size(b)
        call random_number(v)
        mh=0
        mh=matmul(transpose(W),v)
        do i=1,n
                mh(i)=sigmoid((mh(i)+b(i))/temp)
        end do
        h=0
        call random_number(htmp)
        do i=1,n
                if (htmp(i) < mh(i)) then
                        h(i)=1
                endif
        end do

        return
end function SampleHidden_p

function SampleVisible(h,W,a,temp) result(v)

        ! sample the visible layer of the RBM
        ! use tower sampling to fulfill constraint
        ! on spin up/down electrons number

        implicit none

        double precision, intent(in), dimension(:,:) :: W
        double precision, intent(in), dimension(:) :: a
        integer,dimension(:), intent(in) :: h
        integer, dimension(size(a)) :: v
        double precision, dimension(size(a)) :: mv,mv_
        double precision, intent(in) :: temp
        double precision :: p,c
        integer :: i,k

        mv=0
        mv=matmul(W,h)
        do i=(1+frozen),N
                mv(i)=sigmoid((mv(i)+a(i))/temp)
        end do

        v=0
        v(:frozen)=1
        mv(:frozen)=0
        mv_=mv

        !tower sampling
        do k=1,(Ne-frozen)/2
          c=0 !cumulative
          call random_number(p)
          p=p*sum(mv(1::2)) !draw number inside cumulative
          do i=frozen+1,N,2
            c=c+mv(i)
            if (p <= c) then
              v(i)=1
              mv(i)=0 ! set to zero for next draws
              exit
            endif         
          enddo
        enddo

        do k=1,(Ne-frozen)/2
          c=0 !cumulative
          call random_number(p)
          p=p*sum(mv(2::2)) !draw number inside cumulative
          do i=frozen+2,N,2
            c=c+mv(i)
            if (p <= c) then
              v(i)=1
              mv(i)=0 ! set to zero for next draws
              exit
            endif
          enddo
        enddo

       return
end function SampleVisible

subroutine Sampling(v,h,W,a,b,temp)

        ! one step Gibbs sampling
        
        implicit none

        double precision,intent(in),dimension(:,:) :: W
        double precision,intent(in),dimension(:) :: a,b
        integer,intent(out),dimension(size(a)) :: v
        integer,intent(out),dimension(size(b)) :: h
        double precision, intent(in) :: temp


        h=SampleHidden_p(W,b,temp)
        v=SampleVisible(h,W,a,temp)

end subroutine Sampling

subroutine Sample_SD(dets,v_,h,W,a,b,temp)

        implicit none

        double precision,intent(in),dimension(:,:) :: W
        double precision,intent(in),dimension(:)   :: a,b
        integer,intent(in),dimension(:)            :: dets
        integer,dimension(size(a))                 :: v
        integer,dimension(size(a))                 :: vtmp
        integer,intent(inout),dimension(size(a))   :: v_
        double precision, dimension(size(a))       :: mv
        integer,dimension(size(b))                 :: h
        integer,dimension(Ne)                      :: o !occupied so
        integer,dimension((N-Ne))                  :: u !unocc so
        double precision, dimension(Ne,(N-Ne))     :: jump
        integer                                    :: dec,div,i,j,k,io,iu,sd,id
        double precision, intent(in)               :: temp ! temperature
        real                                       :: p,c,cc


        call random_number(p)
        id=nint(p*(int(size(dets))-1)+1) !+2 <- +1 for index and +1 for HF
        dec=dets(id)
        v=0
        div=dec
        do i=1,N
        v(N+1-i)=mod(div,2)
        div=div/2
        enddo

        !!! -----------------------------
        !!! sampling of the proba of so
        !!! ----------------------------
        v_=v ! out determinant
        vtmp=0 ! temp copy of base det
 
        h=SampleHidden_p(W,b,temp)
        vtmp=SampleVisible(h,W,a,temp)
        h=SampleHidden(vtmp,W,b,temp)

        mv=0
        mv=matmul(W,h)

        do i=(1+frozen),N
          mv(i)=sigmoid((mv(i)+a(i))/temp)
        end do
        mv(:frozen)=1
      
        !------------------------------
        !!! single or double 
        sd=2
        call random_number(p)
        if (p<0.5) sd=1
       
        !------------------------------
        ! occupied & unnoccupied so
        !-----------------------------
        io=1
        iu=1
        do i=1,N
         if (v(i)==0) then 
           u(iu)=i
           iu=iu+1
         else
           o(io)=i
           io=io+1
         endif
        enddo

        !----------------------------
        ! jump probabilities
        ! ---------------------------
        cc=0 ! cumsum
        jump=0
        do i=1,Ne
           do j=1,(N-Ne)
             if (mod(o(i)+u(j),2)==0) then
               jump(i,j)=(1-mv(o(i)))*(mv(u(j)))
               cc=cc+jump(i,j)
             endif
           enddo
         enddo

        !--------------------------
        ! sample jump 
        !--------------------------
        do k=1,sd
         c=0
         call random_number(p)
         p=p*sum(jump(:,:))
         do i=1,Ne
           do j=1,N-Ne
             c=c+jump(i,j) !cumulative(i,j)
             if (p < c) then
               v_(u(j))=1
               v_(o(i))=0
               jump(i,:)=0 !cant jump from i anymore
               jump(:,j)=0 !cant jump to j anymore
               goto 111 
             endif
           enddo
         enddo
        111 continue
        enddo

end subroutine Sample_SD

subroutine generate(Dn,W,a,b,dets,Ngen,Ntry,cnt)

        ! Generation of new determinants : direct algorithm
        ! Dn : empty array to hold new determinants
        ! W,a,b : RBM parameters
        ! dets : determinants already in the pruned wavefunction
        ! Ngen : number of new determinants we want to generate
        ! Ntry : max number of trials
        ! cnt : actual number of determinant generated

        implicit none

        double precision,intent(in),dimension(:,:) :: W
        double precision,intent(in),dimension(:) :: a,b
        double precision :: tempgegen
        integer,intent(in) :: Ngen,Ntry
        integer,intent(inout) :: cnt
        integer,intent(in),dimension(:) :: dets
        integer,dimension(size(dets)) :: dets_sorted
        integer,dimension(:),allocatable,intent(inout) :: Dn
        integer,dimension(Ntry/1000+mod(Ntry/1000,2)) :: Dn_tmp !Dn_tmp even
        integer,dimension(size(Dn_tmp)*2) :: Dn_
        integer :: i,k,ir,nu
        integer :: dec,dec_bar
        integer, dimension(size(a)) :: v_,v_bar
        integer, dimension(size(b)) :: h
        integer, dimension(N) :: d
        logical :: sdcond,isin

        Dn=0
        dets_sorted=dets
       
        ! Here I sort the determinants (their decimal rep.) in a
        ! new array dets_sorted
        ! when dets are generated it is then much quicker to search if
        ! they are already inside the wavefunction
        call quicksort(dets_sorted,1,size(dets_sorted))
       

        tempgegen=tempgen
        cnt=1

        ! do by blocks (1000 times Ntry/1000 tries)
        ! For each block we use Dn_tmp to hold new determinants
        ! The parallelisation is done inside those blocks
        do k=1,1000
        
                Dn_tmp=0 
                !$omp parallel private(v_,v_bar,sdcond,isin,dec,dec_bar,d,h,ir) 
                v_=0
                h=0
                !$omp do
                do i=1,Ntry/1000,2 ! increment by 2 : we generate spin symmetric dets
                  
                  ! generate a det in v_ 
                  call Sampling(v_,h,W,a,b,tempgegen)

                  if (sym) then 
                    ir=s2ir(v_,so_ir,M_sym) ! compute irred. rep. of det

                    !if ir <= ir_t : ideally should be if ir == ir_t but
                    !MO from QP and ir from Orca doesn't always agree
                    if  ((ir .le. ir_t)) then 
                      goto 900 ! go to next step
                    else 
                      goto 1000 ! try again
                    endif
                  else 
                      goto 900 ! go to next step
                  endif

900               dec=bin2dec(v_) ! transform generated det to decimal representation

                  ! check if dec belong to the set of single/doubles substitutions of the det 
                  ! inside the wavefunction (much much muuuuch faster using decimal rep)
                  call is_SD_of_L(dec,dets,sdcond) 
                   if (sdcond) then
                     ! check if dec belong to the current wavefunction 
                     ! again much faster using sorted decimal rep
                     isin=binary_search(dets_sorted, dec)

                     ! it passes all tests, we add it to Dn_tmp
                     if (.not. isin) then
                       Dn_tmp(i)=dec

                       ! For open-shell systems, in a spin up/down representation
                       ! det=(a,b), the symmetric determinant (b,a) 
                       ! possesses the same coeff in the wavefunction
                       ! so I directly add this det in D_new
                       ! (Actually it doesn't seem change a lot the efficiency)
                       v_bar(::2)=v_(2::2)
                       v_bar(2::2)=v_(::2)
                       dec_bar=bin2dec(v_bar)
                       if (dec_bar .ne. dec .and. i .ne. size(Dn_tmp)) then
                         Dn_tmp(i+1)=dec_bar
                       endif
                     endif
                   endif
1000            continue 
                enddo
                !$omp enddo
                !$omp end parallel
                
                !to see where we at for big systems
                !open(unit=69,file='progression')
                !write(69,*) "current iteration k = ",k,"/1000, current ngen produced =", cnt
                !close(69)
                
                ! the temperature in the sampling is raised by 1 every 50 block of Ntry/1000
                ! The idea is that in the begining, you need only ~ 1 block of Ntry/1000 to
                ! generate all your Ngen dets
                ! Near convergence it becomes hard to find the remaining important determinants
                ! So this helps the RBM to sample more out of equilibrium
                ! This is also much more computationally efficient  
                tempgegen=tempgegen+(k/50)

                ! Sort the generated determinants
                call quicksort(Dn_tmp,1,size(Dn_tmp))
                i=1
                do while (Dn_tmp(i)==0)
                  i=i+1
                enddo

                ! The parallelization can produce doublons so we check for them here
                ! Dn_ holds both Dn and Dn_tmp
                ! I_unista(Dn_,nu) return Dn_ with nu unique elements
                ! in first position
      
                Dn_(:size(Dn_tmp)-i)=Dn_tmp(i+1:)
                Dn_( size(Dn_tmp)+1-i : size(Dn_tmp)-i+cnt )=Dn(:cnt)
                call I_unista(Dn_(:size(Dn_tmp)+cnt),nu)
                cnt=nu
                if (cnt >= Ngen) then 
                  Dn=Dn_(:Ngen)
                  goto 1001
                else
                  Dn(:nu)=Dn_(:nu)
                endif

        enddo
1001    continue
!        print *, 'total number of dets passed :    ',cnt
!        print *, 'Ntry max :                       ',Ntry
end subroutine

subroutine generate_SD(Dn,W,a,b,dets,Ngen,Ntry,cnt)

        ! Generation of new determinants : SD algorithm
        ! Dn : empty array to hold new determinants
        ! W,a,b : RBM parameters
        ! dets : determinants already in the pruned wavefunction
        ! Ngen : number of new determinants we want to generate
        ! Ntry : max number of trials
        ! cnt : actual number of determinant generated


        implicit none

        double precision,intent(in),dimension(:,:) :: W
        double precision,intent(in),dimension(:) :: a,b
        double precision :: tempgegen
        integer,intent(in) :: Ngen,Ntry
        integer,intent(inout) :: cnt
        integer,intent(in),dimension(:) :: dets
        integer,dimension(size(dets)) :: dets_sorted
        integer,dimension(:),allocatable,intent(inout) :: Dn
        integer,dimension(Ntry/1000+mod(Ntry/1000,2)) :: Dn_tmp !Dn_tmp even
        integer,dimension(size(Dn_tmp)*2) :: Dn_
        integer :: i,k,ir,nu
        integer :: dec,dec_bar
        integer, dimension(size(a)) :: v_,v_bar
        integer, dimension(size(b)) :: h
        integer, dimension(N) :: d
        logical :: sdcond,isin

        ! dets are sorted in dets_sorted for fast comparison
        dets_sorted=dets
        call quicksort(dets_sorted,1,size(dets_sorted))

        tempgegen=tempgen
        cnt=1

        ! do by blocks (1000 times Ntry/1000 tries)
        ! For each block we use Dn_tmp to hold new determinants
        ! The parallelisation is done inside those blocks
        do k=1,1000
        
                Dn_tmp=0
                !$omp parallel private(v_,v_bar,sdcond,isin,dec,dec_bar,d,h,ir) 
                v_=0
                h=0
                !$omp do
                do i=1,Ntry/1000,2
                  
                  ! generate a single/double susbtitution of a det inside dets
                  call Sample_SD(dets,v_,h,W,a,b,tempgegen)

                  if (sym) then 
                    ir=s2ir(v_,so_ir,M_sym)
                    if  ((ir .le. ir_t)) then
                      goto 900 ! go to next step
                    else 
                      goto 1000 ! try again
                    endif
                  else 
                      goto 900 ! go to next step
                  endif

900               dec=bin2dec(v_)
                  isin=binary_search(dets_sorted, dec)
                  if (.not. isin) then
                       Dn_tmp(i)=dec
                       v_bar(::2)=v_(2::2)
                       v_bar(2::2)=v_(::2)
                       dec_bar=bin2dec(v_bar)
                       if (dec_bar .ne. dec .and. i .ne. size(Dn_tmp)) then
                         Dn_tmp(i+1)=dec_bar
                       endif
                  endif
1000            continue 
                enddo
                !$omp enddo
                !$omp end parallel
!                open(unit=69,file='progression')
!                write(69,*) "current iteration k = ",k,"/1000, current ngen produced =", cnt
!                close(69)
                
                tempgegen=(k+50)/50 

                call quicksort(Dn_tmp,1,size(Dn_tmp))
                i=1
                do while (Dn_tmp(i)==0)
                  i=i+1
                enddo

                !*****************************************************!
                ! Dn_ holds both Dn and Dn_tmp
                ! I_unista(Dn_,nu) return Dn_ with nu unique elements
                ! in first position

                Dn_(:size(Dn_tmp)-i)=Dn_tmp(i+1:)
                Dn_( size(Dn_tmp)+1-i : size(Dn_tmp)-i+cnt )=Dn(:cnt)
                call I_unista(Dn_(:size(Dn_tmp)+cnt),nu)
                cnt=nu
                if (cnt >= Ngen) then
                  Dn=Dn_(:Ngen)
                  goto 1001
                else
                  Dn(:nu)=Dn_(:nu)
                endif


        enddo
1001    continue
!        print *, 'total number of dets passed :    ',cnt
!        print *, 'Ntry max :                       ',Ntry
end subroutine
 
        
subroutine update_rbm(a,b,W,V,bs,Ns,eps,l)

        ! one training step of the RBM

        implicit none

        double precision, intent(inout), dimension(:,:) :: W
        integer, intent(in), dimension(:,:) :: V
        double precision, intent(inout), dimension(:) :: a,b
        integer,intent(in) :: Ns,bs
        double precision, intent(in) :: l,eps
        double precision, dimension(size(a)) :: da
        double precision, dimension(size(b)) :: db
        double precision, dimension(size(W,1),size(W,2)) :: dw
        integer :: i,k,q
        integer, dimension(size(a)) :: vm,vd
        integer, dimension(size(b)) :: hm,h
        double precision :: temp=1.,error

        error=0       
 
        vm=0
        hm=0
        h=0
        do i=1,Ns/bs
          dA=0
          dB=0
          dW=0
          do k=1,bs
            vd=V(i,:)
            h=SampleHidden(vd,W,b,temp)
            hm=SampleHidden_p(W,b,temp)
            vm=SampleVisible(hm,W,a,temp)

            !20 steps Gibss sampling
            do q=1,20
              hm=SampleHidden(vm,W,b,temp)
              vm=SampleVisible(hm,W,a,temp)
            enddo     

            ! update       
            dA=dA + ( real( vd-vm )/real(bs) )
            dB=dB + ( real(h-hm)/real(bs) )
            dW=dW + ( real(outer_int(vd,h) - outer_int(vm,hm))/real(bs) )
          enddo

          ! to use for regularization
          !dA=dA*(1+l*eps) 
          !dB=dB*(1+l*eps)
          !dW=dW*(1+l*eps)
          A=A+dA*eps
          B=B+dB*eps
          W=W+dW*eps
        enddo
        !error=0
        !error=norm2(dA)/size(a)
        !print *, "(<V_d> - <V_rbm>)/ #param : ",error!/(size(a)+size(b)+size(a)*size(b))
end subroutine
       
function s2ir(s,soir,M) result(ir)

        ! return the irreducible representation of a determinant
        ! s : determinant in spin-orbital occupation form ([1,1,0,..])
        ! soir : irreducible representation of the spin orbitals
        ! M : multiplication table for the group

        implicit none
        integer,dimension(N),intent(in) :: s,soir
        integer,dimension(:,:),intent(in) :: M
        integer,dimension(N) :: irs
        integer :: ir,cnt,i

        cnt=1 ! from so to irreducible representation
        do i=1,size(s)
         if (s(i)==1) then
          irs(cnt)=soir(i)
          cnt=cnt+1
         endif
        enddo
     
        ir=M(irs(1),irs(2)) ! Multiplication table
        do i=3,sum(s)
          ir=M(ir,irs(i))
        enddo

        return
end function

end module rbm
