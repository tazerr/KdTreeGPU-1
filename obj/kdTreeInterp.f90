module call_kd_tree
    use iso_c_binding
    implicit none
    interface

        subroutine custom_funct(coordinates, numPoints, numQuerys, query, results, numResults, rootIdx) bind(C, name="custom_funct")
            use iso_c_binding
            real(c_float), dimension(*) :: coordinates
            integer(c_int), value :: numPoints
            integer(c_int), value :: numQuerys
            real(c_float), dimension(*) :: query
            real(c_float), dimension(*) :: results
			integer(c_int), value :: numResults
			integer(c_int), intent(inout) :: rootIdx
        end subroutine

		subroutine search_funct(coordinates, numDimensions, numQuerys, query, results, numResults, rootIdx) bind(C, name="search_funct")
            use iso_c_binding
            real(c_float), dimension(*) :: coordinates
            integer(c_int), value :: numDimensions
            integer(c_int), value :: numQuerys
            real(c_float), dimension(*) :: query
			real(c_float), dimension(*) :: results
			integer(c_int), value :: numResults
			integer(c_int), value:: rootIdx
        end subroutine

    end interface
end module call_kd_tree


program main

use call_kd_tree
implicit none

integer(c_int) :: numPoints, numResults, rootIdx    
integer :: numDimensions = 3
integer :: sizeOfarray 
real(c_float), allocatable :: coordinates(:), results(:)
real(c_float), allocatable :: query(:)
integer(c_int) :: numQuerys
integer :: nblocks, i, j, k, index
integer, allocatable :: NI(:), NJ(:), NK(:)
real, allocatable :: xgrid(:,:,:,:), ygrid(:,:,:,:), zgrid(:,:,:,:)
character(100) :: filename="grid.xyz"
real :: Lx, Ly, Lz
Real,ALLOCATABLE :: Igrid(:,:,:,:), Jgrid(:,:,:,:), Kgrid(:,:,:,:)
Real,ALLOCATABLE :: Xgridnew(:,:,:,:), Ygridnew(:,:,:,:), Zgridnew(:,:,:,:)
Real,ALLOCATABLE :: Igridnew(:,:,:,:), Jgridnew(:,:,:,:), Kgridnew(:,:,:,:)
integer nbl, NImax, NJmax, NKmax
integer nblocksnew,NImaxnew, NJmaxnew, NKmaxnew
Integer,ALLOCATABLE :: NInew(:), NJnew(:), NKnew(:)


    !**********************************Background block*****************************************

nblocks = 1


Lx = 12.d0
Ly = 12.d0
Lz = 3.d0

ALLOCATE(NI(nblocks))
ALLOCATE(NJ(nblocks))
ALLOCATE(NK(nblocks)) 
 
nbl =1
NI(nbl) = 256
NJ(nbl) = 256
NK(nbl) = 1

NImax = maxval(NI)
NJmax = maxval(NJ)
NKmax = maxval(NK)	
  
ALLOCATE(Xgrid(NImax,NJmax,NKmax,nblocks))
ALLOCATE(Ygrid(NImax,NJmax,NKmax,nblocks))
ALLOCATE(Zgrid(NImax,NJmax,NKmax,nblocks))
   
ALLOCATE(Igrid(NImax,NJmax,NKmax,nblocks))
ALLOCATE(Jgrid(NImax,NJmax,NKmax,nblocks))
ALLOCATE(Kgrid(NImax,NJmax,NKmax,nblocks))
   
DO nbl = 1,nblocks
    DO k = 1,NK(nbl)
        DO j = 1,NJ(nbl)
            DO i = 1,NI(nbl)

        Xgrid(i,j,k,nbl) = -6.d0 + (Lx/(NI(nbl)-1.d0)) * ((i-1.d0) + sin(6.d0*3.14d0*(j-1)*(1.d0/(NJ(nbl)-1))))
        Ygrid(i,j,k,nbl) = -6.d0 + (Ly/(NJ(nbl)-1.d0)) * ((j-1.d0) + 2.d0 * sin(6.d0*3.14d0*(i-1)*(1.d0/(NI(nbl)-1))))
        Zgrid(i,j,k,nbl) = Lz*(k-1.d0)
		
		Igrid(i,j,k,nbl) = 1.d0 * i
		Jgrid(i,j,k,nbl) = 1.d0 * j
		Kgrid(i,j,k,nbl) = 1.d0 * k
		
            END DO
        END DO
    END DO
END DO

! Print the values for verification
  	!do nbl = 1, nblocks
    !	write(*,*) "Block:", nbl
    !	do k = 1, NK(nbl)
    !	  do j = 1, NJ(nbl)
    !	    do i = 1, NI(nbl)
    !	      write(*,*) "xgrid(", i, j, k, nbl, ") =", xgrid(i, j, k, nbl)
    !	      write(*,*) "ygrid(", i, j, k, nbl, ") =", ygrid(i, j, k, nbl)
    !	      write(*,*) "zgrid(", i, j, k, nbl, ") =", zgrid(i, j, k, nbl)
    !	    end do
    !	  end do
    !	end do
  	!end do
!***************************************Second grid******************************************
nblocksnew = 1


Lx = 6.d0
Ly = 6.d0
Lz = 3.d0

ALLOCATE(NInew(nblocksnew))
ALLOCATE(NJnew(nblocksnew))
ALLOCATE(NKnew(nblocksnew))
 
nbl =1
NInew(nbl) = 5
NJnew(nbl) = 5
NKnew(nbl) = 1
   	
NImaxnew = maxval(NInew)
NJmaxnew = maxval(NJnew)
NKmaxnew = maxval(NKnew)	
  
ALLOCATE(Xgridnew(NImaxnew,NJmaxnew,NKmaxnew,nblocksnew))
ALLOCATE(Ygridnew(NImaxnew,NJmaxnew,NKmaxnew,nblocksnew))
ALLOCATE(Zgridnew(NImaxnew,NJmaxnew,NKmaxnew,nblocksnew))
   
ALLOCATE(Igridnew(NImaxnew,NJmaxnew,NKmaxnew,nblocksnew))
ALLOCATE(Jgridnew(NImaxnew,NJmaxnew,NKmaxnew,nblocksnew))
ALLOCATE(Kgridnew(NImaxnew,NJmaxnew,NKmaxnew,nblocksnew))

DO nbl = 1,nblocksnew
    DO k = 1,NKnew(nbl)
		DO j = 1,NJnew(nbl)
            DO i = 1,NInew(nbl)

        Xgridnew(i,j,k,nbl) = -3.d0 + Lx*(i-1.d0)/(NInew(nbl)-1.d0)
        Ygridnew(i,j,k,nbl) = -3.d0 + Ly*(j-1.d0)/(NJnew(nbl)-1.d0)
        Zgridnew(i,j,k,nbl) = Lz*(k-1.d0)

            END DO
		END DO
    END DO
END DO

!***************************************Prepartion to call creation of KD tree******************************************

numPoints = size(xgrid)
sizeOfarray = numPoints*numDimensions
numQuerys = size(Xgridnew)
numResults = 2
allocate(results(numResults*numDimensions*numQuerys))

allocate(coordinates(sizeOfarray))
do nbl = 1, nblocks
    do k = 1, NK(nbl)
        do j = 1, NJ(nbl)
            do i = 1, NI(nbl)
            
            ! Get the current index
            !index = (nbl-1)*NI(nbl)*NJ(nbl)*NK(nbl)*numDimensions + (k-1)*NI(nbl)*NJ(nbl)*numDimensions + (j-1)*NI(nbl)*numDimensions +  (i-1)*numDimensions + 1

            index = (nbl-1)*NI(nbl)*NJ(nbl)*NK(nbl)*numDimensions + &
                (k-1)*NI(nbl)*NJ(nbl)*numDimensions + &
                (j-1)*NI(nbl)*numDimensions + &
                (i-1)*numDimensions + 1

            
            !store in appropriate indices
            coordinates(index) = Xgrid(i,j,k,nbl)*100000
            coordinates(index+1) = Ygrid(i,j,k,nbl)*100000
            coordinates(index+2) = Zgrid(i,j,k,nbl)*100000
            
            end do
        end do
    end do
end do
!***************************************call creation of KD tree******************************************

call custom_funct(coordinates, numPoints, numQuerys, query, results, numResults, rootIdx)

!***************************************Prepartion to call searching of KD tree******************************************

allocate(query(numQuerys*numDimensions))

do nbl = 1, nblocksnew
    do k = 1, NKnew(nbl)
        do j = 1, NJnew(nbl)
            do i = 1, NInew(nbl)
            
            ! Get the current index
            index = (nbl-1)*NInew(nbl)*NJnew(nbl)*NKnew(nbl)*numDimensions + &
                (k-1)*NInew(nbl)*NJnew(nbl)*numDimensions + &
                (j-1)*NInew(nbl)*numDimensions + &
                (i-1)*numDimensions + 1

            
            !store in appropriate indices
            query(index) = Xgridnew(i,j,k,nbl)
            query(index+1) = Ygridnew(i,j,k,nbl)
            query(index+2) = Zgridnew(i,j,k,nbl)
            
            end do
        end do
    end do
end do

!***************************************call searching of KD tree******************************************

call search_funct(coordinates, numDimensions, numQuerys, query, results, numResults, rootIdx) 	

!***************************************Print Result******************************************
!i=1
!j=1

!write(*,*) "Closest points points found are as follows : "
!do while (i<=size(results))
!	write(*, *) "X (", j, ") : ", results(i)
!	write(*, *) "Y (", j, ") : ", results(i+1)
!	write(*, *) "Z (", j, ") : ", results(i+2)
!	j=j+1
!	i=i+3
!end do

end program main

!*************************************Local Index subroutine***************************************
SUBROUTINE get_loc_index(global_index,blk,iloc,jloc,kloc)
use arrs
    integer global_index, rest_index
    integer blk, iloc, jloc, kloc, nbl
    integer sblk,eblk
    sblk=0
    !$acc data present(startblkl, endblkl, NI, NJ, NK)

    DO nbl = 1,nblocks
        eblk=NI(nbl)*NJ(nbl)*NK(nbl)
        IF ((global_index>sblk).and.(global_index<=eblk)) THEN 
            blk = nbl
            EXIT
        END IF
        sblk=eblk
    END DO

	rest_index = global_index - sblk
	kloc = int(rest_index/(NJ(blk)*NI(blk)))
	rest_index = rest_index - kloc*NI(blk)*NJ(blk)
	jloc = int(rest_index/(NI(blk)))
	rest_index = rest_index - jloc*NI(blk)
	iloc = int(rest_index)
	
	iloc = iloc+1
	jloc = jloc+1
	kloc = kloc+1
	
END	