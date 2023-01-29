module tools
        implicit none
contains

function count_prun(cc,tol) result(c)

        ! pruning : count determinants to keep
        ! not used rn, handled by python scripts

        implicit none
        double precision,intent(in),dimension(:) :: cc
        double precision,intent(in) :: tol
        integer :: c,i
        
        c=0
        do i=1,size(cc)
          if (cc(i) > tol) then 
            c=c+1
          endif
        enddo
        
        return
end function

subroutine pruning(cc,Dp,c,cc_pruned,Dp_pruned)

        ! not used rn, handled by python scripts

        implicit none
        
        integer,intent(in) :: c
        double precision, intent(in), dimension(:) :: cc
        integer, intent(in), dimension(:) :: Dp
        double precision, intent(out), dimension(:) :: cc_pruned
        integer, intent(out), dimension(:) :: Dp_pruned
        
        cc_pruned=cc(:c)
        Dp_pruned=Dp(:c)
           
end subroutine

function remove_zeros(D_) result(D)

        ! before Determinants generation
        ! a bigger array is defined and initialized to zero
        ! to avoid loosing time in memory reallocation
        ! this count the actual number of generated determinants
        ! and return the corresponding array

        implicit none

        integer, intent(in), dimension(:) :: D_
        integer, dimension(:),allocatable :: D
        integer :: i,N,cnt

        N=size(D_)
        cnt=0
        do i=1,N
         if (D_(i) > 0) then
          cnt=cnt+1
         else
          exit
         endif
        enddo
        allocate( D( cnt ))
        D=D_(:cnt)
        return
end function
        
logical function equal( array1, array2 )
        integer, dimension(:), intent(in) :: array1, array2
        integer :: i

        equal =size(array1) == size(array2)
        if ( equal ) then
        do i = 1,size(array1)
        equal = array1(i) == array2(i)
        if ( .not. equal )exit
        enddo
        endif
end function equal

subroutine  is_SD_of_L(a,L,isit)

        ! this checks if determinant a belongs
        ! to the set of single/double substituted 
        ! determinants in list L
        ! use bitwise operators to compare decimal rep of determinants

        implicit none 

        integer,intent(in) :: a
        integer, dimension(:), intent(in) :: L
        integer :: i,dif
        logical,intent(inout) :: isit

        isit=.false.
        do i=1,size(L)
          dif=ieor(a,L(i)) !xor
          if (popcnt(dif) .le. 4) then !number of ones in dif
            isit=.true.
            exit
          endif
        enddo
end subroutine

function binary_search(A, T) result(res)

    ! is number T inside list A ?
    ! binary search algorithm 
    ! to check if generated determinant is already 
    ! in the wavefunction
    
    integer,dimension(:) :: A 
    integer :: T,L,R,m
    logical :: res
    
    res=.false.
    L=1
    R = size(A)
    do while (L .le. R)
        m = ((L + R) / 2)
        if (A(m) .lt. T) then
            L = m + 1
        else if (A(m) > T) then
            R = m-1
        else
            res=.true.
            exit
        endif
    enddo
    return
end function

function rargsort(a) result(b)
        ! Returns the indices that would sort an array.
        ! from smallest to greatest 
        ! Arguments
        ! ---------
        !
        double precision, intent(in):: a(:)   ! array of numbers
        integer :: b(size(a))         ! indices into the array 'a' that sort it
        !
        ! Example
        ! -------
        !
        ! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]

        integer :: N                           ! number of numbers/vectors
        integer :: i,imin                      ! indices: i, i of smallest
        integer :: temp1                       ! temporary
        real*8 :: temp2
        real*8 :: a2(size(a))
        a2 = a
        N=size(a)
        do i = 1, N
            b(i) = i
        end do
        do i = 1, N-1
            ! find ith smallest in 'a'
            imin = minloc(a2(i:),1) + i - 1
            ! swap to position i in 'a' and 'b', if not already there
            if (imin /= i) then
                temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
                temp1 = b(i); b(i) = b(imin); b(imin) = temp1
            end if
        end do
        return
end function

function rargsort_int(a) result(b)

        ! not mine : Based on code written by John E. Pask, LLNL.
        ! see : https://github.com/certik/fortran-utils/blob/master/src/sorting.f90

        ! Returns the indices that would sort an array.
        ! from smallest to greatest 
        ! Arguments
        ! ---------
        !
        integer, intent(in):: a(:)   ! array of numbers
        integer :: b(size(a))         ! indices into the array 'a' that sort it
        !
        ! Example
        ! -------
        !
        ! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]

        integer :: N                           ! number of numbers/vectors
        integer :: i,imin                      ! indices: i, i of smallest
        integer :: temp1                       ! temporary
        integer :: temp2
        integer :: a2(size(a))
        a2 = a
        N=size(a)
        do i = 1, N
            b(i) = i
        end do
        do i = 1, N-1
            ! find ith smallest in 'a'
            imin = minloc(a2(i:),1) + i - 1
            ! swap to position i in 'a' and 'b', if not already there
            if (imin /= i) then
                temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
                temp1 = b(i); b(i) = b(imin); b(imin) = temp1
            end if
        end do
        return
end function


function reverse(a) result (b)

        ! reverse the order of an array

        integer,intent(in) :: a(:)
        integer :: b(size(a))
        integer :: N,i
        N=size(a)
        do i=1,N
          b(i)=a(N+1-i)
        end do
       
        return
end function


function s2d(s,Ne) result(d)

        ! not used anymore
        ! transform from spinorbital occupation 
        ! to occupied spinorbital number
        ! eg. [1,0,1,0,1,1] --> [0,2,4,5]

        implicit none

        integer,intent(in),dimension(:) :: s
        integer,intent(in) :: Ne
        integer,allocatable,dimension(:) :: d
        integer :: i,cnt
       
        allocate(d(Ne))
        d=0
        cnt=1
        do i=1,size(s)
          if (s(i)==1) then 
           d(cnt)=(i-1)
           cnt=cnt+1
          endif
        enddo
        return
        
end function s2d

function d2s(d,N) result(s)

        ! not used anymore
        ! transform from occupied spinorbital number
        ! to spinorbital occupation
        ! eg. [0,2,4,5] --> [1,0,1,0,1,1]

        implicit none

        integer,intent(in),dimension(:) :: d
        integer,intent(in) :: N
        integer,dimension(:),allocatable :: s
        integer :: i

        allocate(s(N))
        s=0
        
        do i=1,size(d)
                s(d(i)+1)=1
        end do
        
        return
end function d2s

function sample_n_from_p(p,n,el) result(ind)

        ! not used anymore
        ! handled by python scripts
        ! to create the training sample for the RBM

        ! el allows to choose either a log probability (el > 0)
        ! or (el < 0)  the actual probability without the first biggest |el| elements
        ! eg. for el = -1 we exclude the Hartree-Fock determinant

        implicit none

        double precision,intent(in) :: p(:)
        integer,intent(in) :: el
        double precision,dimension(size(p)) :: cumsum,q
        integer :: n,i,j,shift
        integer, dimension(n) :: ind
        double precision, dimension(n) :: d

        if (el > 0) then
          do i=1,size(p)
                q(i)=exp(log(p(i))*(1.0/dble(el)))
          enddo
        endif
        if (el < 0) then 
          shift=(-el)+1
        q=p
        print *, "one two one two, size(q), size(p) :", size(q), size(p)
        q(:shift)=0
        endif
        if (el .eq. 0) q=p
        q=q/sum(q)

        do i=1,size(q)
          cumsum(i)=sum(q(1:i))/sum(q)
        end do
        call random_number(d)
        
        do j=1,n
          do i=1,size(p)
           if (d(j) < cumsum(i)) then
             ind(j)=i
             exit
           endif
         enddo
       enddo
       return
end function

function bin2dec(bin) result(dec)

  ! transform a binary array [1,1,0,..,1]
  ! to the decimal representation of the number 110..1

  implicit none
  integer,dimension(:) :: bin
  integer :: dec,i,j

  dec=0
  do i=1,size(bin)
    j=size(bin)-i
    dec=dec+bin(i)*(2**j)
  enddo
  return
end function bin2dec

recursive subroutine quicksort(a, first, last)
 
  ! not mine : Author: t-nisie
  ! see : https://gist.github.com/t-nissie/479f0f16966925fa29ea

  ! sort an array a

  implicit none
  integer  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort

subroutine unique(la,lb)

 ! given la, lb, changes lb to the
 ! intersection set of la and lb
 ! until size(lb)

 implicit none
 integer( kind = 8 ) :: kx,nu,la(:),lb(:)
 integer( kind = 8 ),allocatable :: indices(:),list(:)
 logical :: mask(size(la)+size(lb))

 allocate(list(size(la)+size(lb)))
 list(:size(la))=la
 list(size(la)+1:)=lb

 mask(1)=.true.
 do kx=size(list),2,-1
   mask(kx)= .not.(any(list(:kx-1)==list(kx)))
 end do
 indices=pack([( kx,kx=1,size(list) ) ],mask)

 nu=size(indices)
 list(:nu)=list(indices)
 list(nu+1:)=0

 lb=list

end subroutine

! The following routines uniinv, nearless and I_unista are not mine:
! written by Michel Olagnon and part of:
! "ORDERPACK 2.0 -- Unconditional, Unique, and Partial Ranking, Sorting, and Permutation Downloadable Fortran 90 source code"
! "http://www.fortran-2000.com/rank/

Subroutine uniinv (XDONT, IGOEST)
! __________________________________________________________
!   UNIINV = Merge-sort inverse ranking of an array, with removal of
!   duplicate entries.
!   The routine is similar to pure merge-sort ranking, but on
!   the last pass, it sets indices in IGOEST to the rank
!   of the value in the ordered set with duplicates removed.
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer, Dimension (:), Intent (Out) :: IGOEST
! __________________________________________________________
      Integer :: XTST, XDONA, XDONB
!
! __________________________________________________________
      Integer, Dimension (SIZE(IGOEST)) :: JWRKT, IRNGT
      Integer :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IGOEST))
!
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IGOEST (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) < XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo (NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 4) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (2*LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!  One steps in the C subset, that we create in the final rank array
!
!  Make a copy of the rank array for the iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
            XDONA = XDONT (JWRKT(IINDA))
            XDONB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XDONA > XDONB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XDONB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XDONA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
!   Last merge of A and B into C, with removal of duplicates.
!
      IINDA = 1
      IINDB = LMTNA + 1
      NUNI = 0
!
!  One steps in the C subset, that we create in the final rank array
!
      JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
      If (IINDB <= NVAL) Then
        XTST = NEARLESS (Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
      Else
        XTST = NEARLESS (XDONT(JWRKT(1)))
      Endif
      Do IWRK = 1, NVAL
!
!  We still have unprocessed values in both A and B
!
         If (IINDA <= LMTNA) Then
            If (IINDB <= NVAL) Then
               If (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) Then
                  IRNG = IRNGT (IINDB)
                  IINDB = IINDB + 1
               Else
                  IRNG = JWRKT (IINDA)
                  IINDA = IINDA + 1
               End If
            Else
!
!  Only A still with unprocessed values
!
               IRNG = JWRKT (IINDA)
               IINDA = IINDA + 1
            End If
         Else
!
!  Only B still with unprocessed values
!
            IRNG = IRNGT (IWRK)
         End If
         If (XDONT(IRNG) > XTST) Then
            XTST = XDONT (IRNG)
            NUNI = NUNI + 1
         End If
         IGOEST (IRNG) = NUNI
!
      End Do
!
      Return
!
End Subroutine uniinv

Function nearless (XVAL) result (I_nl)
!  Nearest value less than given value
! __________________________________________________________
      Integer, Intent (In) :: XVAL
      Integer :: I_nl
! __________________________________________________________
      I_nl = XVAL - 1
      return
!
End Function nearless

Subroutine I_unista (XDONT, NUNI)
!   UNISTA = (Stable unique) Removes duplicates from an array,
!            leaving unique entries in the order of their first
!            appearance in the initial set.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (InOut)  :: XDONT
      Integer, Intent (Out) :: NUNI
! __________________________________________________________
!
      Integer, Dimension (Size(XDONT)) :: IWRKT
      Logical, Dimension (Size(XDONT)) :: IFMPTYT
      Integer :: ICRS
! __________________________________________________________
      Call UNIINV (XDONT, IWRKT)
      IFMPTYT = .True.
      NUNI = 0
      Do ICRS = 1, Size(XDONT)
         If (IFMPTYT(IWRKT(ICRS))) Then
            IFMPTYT(IWRKT(ICRS)) = .False.
            NUNI = NUNI + 1
            XDONT (NUNI) = XDONT (ICRS)
         End If
      End Do
      Return
!
End Subroutine I_unista


end module tools
