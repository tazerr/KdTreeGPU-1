MODULE call_kd_tree
    USE iso_c_binding
    IMPLICIT NONE

    !*************************************C INTERFACE***************************************
    INTERFACE
        SUBROUTINE custom_funct(coordinates, numPoints, numQuerys, query, results, numResults, rootIdx) bind(C, name="custom_funct")
            USE iso_c_binding
            REAL(c_float), dimension(*) :: coordinates
            INTEGER(c_int), value :: numPoints
            INTEGER(c_int), value :: numQuerys
            REAL(c_float), dimension(*) :: query
            REAL(c_float), dimension(*) :: results
			INTEGER(c_int), value :: numResults
			INTEGER(c_int), intent(inout) :: rootIdx
        END SUBROUTINE

		SUBROUTINE search_funct(coordinates, numDimensions, numQuerys, query, results, numResults, rootIdx, mltip) &
         bind(C, name="search_funct")
            USE iso_c_binding
            REAL(c_float), dimension(*) :: coordinates
            INTEGER(c_int), value :: numDimensions
            INTEGER(c_int), value :: numQuerys
            REAL(c_float), dimension(*) :: query
			REAL(c_float), dimension(*) :: results
			INTEGER(c_int), value :: numResults
			INTEGER(c_int), value:: rootIdx
            REAL(c_float), value :: mltip
        END SUBROUTINE
    END INTERFACE
END MODULE call_kd_tree

MODULE arrs
    INTEGER,ALLOCATABLE :: NI(:), NJ(:), NK(:), startblk(:), endblk(:)
    INTEGER nblocks
END MODULE arrs


PROGRAM main

USE call_kd_tree
USE arrs
IMPLICIT NONE

!***********************************declare variables**********************************

INTEGER(c_int) :: numPoints, numResults, rootIdx    
INTEGER :: numDimensions = 3
INTEGER :: sizeOfarray 
REAL(c_float):: mltip=100000.0
REAL(c_float), ALLOCATABLE :: coordinates(:), results(:)
REAL(c_float), ALLOCATABLE :: query(:)
INTEGER(c_int) :: numQuerys
INTEGER :: i, j, k, index
REAL, ALLOCATABLE :: xgrid(:,:,:,:), ygrid(:,:,:,:), zgrid(:,:,:,:)
CHARACTER(100) :: filename="grid.xyz"
REAL :: Lx, Ly, Lz
REAL,ALLOCATABLE :: Igrid(:,:,:,:), Jgrid(:,:,:,:), Kgrid(:,:,:,:)
REAL,ALLOCATABLE :: Xgridnew(:,:,:,:), Ygridnew(:,:,:,:), Zgridnew(:,:,:,:)
REAL,ALLOCATABLE :: Igridnew(:,:,:,:), Jgridnew(:,:,:,:), Kgridnew(:,:,:,:)
INTEGER nbl, NImax, NJmax, NKmax
INTEGER nblocksnew,NImaxnew, NJmaxnew, NKmaxnew
INTEGER,ALLOCATABLE :: NInew(:), NJnew(:), NKnew(:)


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
ALLOCATE(results(numResults*numDimensions*numQuerys))

ALLOCATE(coordinates(sizeOfarray))
DO nbl = 1, nblocks
    DO k = 1, NK(nbl)
        DO j = 1, NJ(nbl)
            DO i = 1, NI(nbl)
            
            ! Get the current index
            !index = (nbl-1)*NI(nbl)*NJ(nbl)*NK(nbl)*numDimensions + (k-1)*NI(nbl)*NJ(nbl)*numDimensions + (j-1)*NI(nbl)*numDimensions +  (i-1)*numDimensions + 1

            index = (nbl-1)*NI(nbl)*NJ(nbl)*NK(nbl)*numDimensions + &
                (k-1)*NI(nbl)*NJ(nbl)*numDimensions + &
                (j-1)*NI(nbl)*numDimensions + &
                (i-1)*numDimensions + 1

            
            !store in appropriate indices
            coordinates(index) = Xgrid(i,j,k,nbl)*mltip
            coordinates(index+1) = Ygrid(i,j,k,nbl)*mltip
            coordinates(index+2) = Zgrid(i,j,k,nbl)*mltip
            
            END DO
        END DO
    END DO
END DO
!***************************************call creation of KD tree******************************************

CALL custom_funct(coordinates, numPoints, numQuerys, query, results, numResults, rootIdx)

!***************************************Prepartion to call searching of KD tree******************************************

ALLOCATE(query(numQuerys*numDimensions))

DO nbl = 1, nblocksnew
    DO k = 1, NKnew(nbl)
        DO j = 1, NJnew(nbl)
            DO i = 1, NInew(nbl)
            
            ! Get the current index
            index = (nbl-1)*NInew(nbl)*NJnew(nbl)*NKnew(nbl)*numDimensions + &
                (k-1)*NInew(nbl)*NJnew(nbl)*numDimensions + &
                (j-1)*NInew(nbl)*numDimensions + &
                (i-1)*numDimensions + 1

            
            !store in appropriate indices
            query(index) = Xgridnew(i,j,k,nbl)
            query(index+1) = Ygridnew(i,j,k,nbl)
            query(index+2) = Zgridnew(i,j,k,nbl)
            
            END DO
        END DO
    END DO
END DO

!***************************************call searching of KD tree******************************************

CALL search_funct(coordinates, numDimensions, numQuerys, query, results, numResults, rootIdx, mltip) 	

!***************************************Print Result******************************************

i=1
j=1
WRITE(*,*) "Closest points points found are as follows : "
DO WHILE (i<=size(results))
	WRITE(*, *) "X (", j, ") : ", results(i)
	WRITE(*, *) "Y (", j, ") : ", results(i+1)
	WRITE(*, *) "Z (", j, ") : ", results(i+2)
	j=j+1
	i=i+3
END DO

END PROGRAM main

!*************************************Local Index subroutine***************************************
SUBROUTINE get_loc_index(global_index,blk,iloc,jloc,kloc)
USE arrs
    INTEGER global_index, rest_index
    INTEGER blk, iloc, jloc, kloc, nbl
    INTEGER sblk,eblk
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
        
END SUBROUTINE