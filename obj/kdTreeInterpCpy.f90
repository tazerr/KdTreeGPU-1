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

		SUBROUTINE search_funct(coordinates, numPoints, numDimensions, numQuerys, query, results, numResults, rootIdx &
            , mltip, gindices, dists, u, v, uinter, vinter) bind(C, name="search_funct")
            !@@$acc routine seq
            USE iso_c_binding
            REAL(c_float), dimension(*) :: coordinates
            INTEGER(c_int), value :: numPoints
            INTEGER(c_int), value :: numDimensions
            INTEGER(c_int), value :: numQuerys
            REAL(c_float), dimension(*) :: query
			REAL(c_float), dimension(*) :: results
			INTEGER(c_int), value :: numResults
			INTEGER(c_int), value:: rootIdx
            REAL(c_float), value :: mltip
            INTEGER(c_int), dimension(*) :: gindices
            REAL(c_double), dimension(*) :: dists
            REAL(c_double), dimension(*) :: u
            REAL(c_double), dimension(*) :: v
            REAL(c_double), dimension(*) :: uinter
            REAL(c_double), dimension(*) :: vinter
        END SUBROUTINE

        SUBROUTINE temp_funct(tarry, tx, ty, tz) bind(C, name="temp_funct")
            USE iso_c_binding
            REAL(c_double), dimension(*) :: tarry
            INTEGER(c_int), value :: tx
            INTEGER(c_int), value :: ty
            INTEGER(c_int), value :: tz
        END SUBROUTINE
    END INTERFACE
END MODULE call_kd_tree

MODULE arrs
    INTEGER,ALLOCATABLE :: NI(:), NJ(:), NK(:), startblk(:), endblk(:)
    INTEGER nblocks
    !$acc declare create(NI, NJ, NK, nblocks)
END MODULE arrs


PROGRAM main

USE call_kd_tree
USE arrs
IMPLICIT NONE

!***********************************declare variables**********************************

INTEGER(c_int) :: numPoints, numResults, rootIdx    
INTEGER :: numDimensions 
INTEGER :: sizeOfarray 
REAL(c_float):: mltip=100000000.0
REAL(c_float), ALLOCATABLE :: coordinates(:), results(:)
INTEGER, ALLOCATABLE :: gindices(:)
REAL(c_double), ALLOCATABLE :: dists(:)
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
INTEGER var,loc_blk,loc_i,loc_j,loc_k
INTEGER interptsx, interptsy, indx, indy, interplace
INTEGER indexing, flag, convertindex
REAL, ALLOCATABLE :: xstencil(:,:), ystencil(:,:), zstencil(:,:), istencil(:), jstencil(:)
REAL, ALLOCATABLE :: Lag_weights(:,:)
REAL iinter,jinter
REAL dw
real time_s, time_einter, time_esearch, time_sinter
REAL(c_double), ALLOCATABLE :: tarry(:,:,:)
INTEGER(c_int) :: tx,ty,tz
REAL(c_double), ALLOCATABLE :: u(:,:,:,:), v(:,:,:,:), uinter(:,:,:,:), vinter(:,:,:,:)
REAL, ALLOCATABLE :: error_u(:,:,:,:), error_v(:,:,:,:)


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


!***************************************Second grid******************************************
nblocksnew = 1


Lx = 6.d0
Ly = 6.d0
Lz = 3.d0

ALLOCATE(NInew(nblocksnew))
ALLOCATE(NJnew(nblocksnew))
ALLOCATE(NKnew(nblocksnew))
 
nbl =1
NInew(nbl) = 128
NJnew(nbl) = 128
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

numDimensions=3
numPoints = size(xgrid)
sizeOfarray = numPoints*numDimensions
numQuerys = size(Xgridnew)
numResults = 8
ALLOCATE(results(numResults*numDimensions*numQuerys))

ALLOCATE(u(NImax,NJmax,NKmax,nblocks))
ALLOCATE(v(NImax,NJmax,NKmax,nblocks))

ALLOCATE(uinter(NImaxnew,NJmaxnew,NKmaxnew,nblocksnew))
ALLOCATE(vinter(NImaxnew,NJmaxnew,NKmaxnew,nblocksnew))

ALLOCATE(error_u(NImaxnew,NJmaxnew,NKmaxnew,nblocksnew))
ALLOCATE(error_v(NImaxnew,NJmaxnew,NKmaxnew,nblocksnew))

ALLOCATE(coordinates(sizeOfarray))
DO nbl = 1, nblocks
    DO k = 1, NK(nbl)
        DO j = 1, NJ(nbl)
            DO i = 1, NI(nbl)
            
            ! Calculate what the global index will be
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

!***************************************Prepartion to call searching of KD tree******************************************

ALLOCATE(query(numQuerys*numDimensions))
ALLOCATE(gindices(numQuerys*numResults))
ALLOCATE(dists(numQuerys*numResults))

DO nbl = 1, nblocksnew
    DO k = 1, NKnew(nbl)
        DO j = 1, NJnew(nbl)
            DO i = 1, NInew(nbl)
            
            ! Calculate what the global index will be
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

!***************************************call creating and searching of KD tree******************************************
call cpu_time(time_s)
print*, time_s

CALL custom_funct(coordinates, numPoints, numQuerys, query, results, numResults, rootIdx)

CALL search_funct(coordinates, numPoints, numDimensions, numQuerys, query, results, numResults, rootIdx, mltip &
    , gindices, dists, u, v, uinter, vinter) 	

!pause
call cpu_time(time_esearch)
print*, time_esearch
print*, "TOTAL TIME TAKEN FOR SEARCH + GETTING RESULTS BACK: ", time_esearch-time_s

!***************************************Write out results******************************************
open(1,file='grid.xyz',form='unformatted')
	  
write(1) nblocks
write(1) ( NI(nbl),NJ(nbl), NK(nbl), nbl=1, nblocks) 

DO nbl=1,nblocks

write(1) (((Xgrid(i,j,k,nbl), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))  &
	        ,(((Ygrid(i,j,k,nbl), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))  &
	        ,(((Zgrid(i,j,k,nbl), i=1,NI(nbl)),j=1,NJ(nbl)),k=1,NK(nbl))
ENDDO!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(1,file='grid2.xyz',form='unformatted')
	  
write(1) nblocksnew
write(1) ( NInew(nbl),NJnew(nbl), NKnew(nbl), nbl=1, nblocks) 
DO nbl=1,nblocksnew

write(1) (((Xgridnew(i,j,k,nbl), i=1,NInew(nbl)),j=1,NJnew(nbl)),k=1,NKnew(nbl))  &
	        ,(((Ygridnew(i,j,k,nbl), i=1,NInew(nbl)),j=1,NJnew(nbl)),k=1,NKnew(nbl))  &
	        ,(((Zgridnew(i,j,k,nbl), i=1,NInew(nbl)),j=1,NJnew(nbl)),k=1,NKnew(nbl))
ENDDO!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open (7, form = 'unformatted', file = 'flow_main.xyz')
		write(7) nblocks
		write(7) (NI(nbl), NJ(nbl), NK(nbl),2, nbl = 1, nblocks)
		Do nbl = 1,nblocks
		write(7) ((( u(i,j,k,nbl), i = 1,NI(nbl)), j = 1,NJ(nbl)),k=1,NK(nbl)) 		&
			, 		((( v(i,j,k,nbl), i = 1,NI(nbl)), j = 1,NJ(nbl)),k=1,NK(nbl)) 
		END DO
		close(7)
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO nbl = 1, nblocksnew
    DO k = 1, NKnew(nbl)
        DO j = 1, NJnew(nbl)
            DO i = 1, NInew(nbl)
			
				error_u(i,j,k,nbl) = abs(uinter(i,j,k,nbl) - (Xgridnew(i,j,k,nbl)**2 + Ygridnew(i,j,k,nbl)**2))
				error_v(i,j,k,nbl) = abs(vinter(i,j,k,nbl) - (sin(Xgridnew(i,j,k,nbl)) * cos(Ygridnew(i,j,k,nbl)**2)))
            
            END DO
        END DO
    END DO
END DO



open (7, form = 'unformatted', file = 'flow_query.xyz')
		write(7) nblocksnew
		write(7) (NInew(nbl), NJnew(nbl), NKnew(nbl),4, nbl = 1, nblocksnew)
		Do nbl = 1,nblocksnew
		write(7) ((( uinter(i,j,k,nbl), i = 1,NInew(nbl)), j = 1,NJnew(nbl)),k=1,NKnew(nbl)) 		&
			, 		((( vinter(i,j,k,nbl), i = 1,NInew(nbl)), j = 1,NJnew(nbl)),k=1,NKnew(nbl))  		&
			, 		((( error_u(i,j,k,nbl), i = 1,NInew(nbl)), j = 1,NJnew(nbl)),k=1,NKnew(nbl))  		&
			, 		((( error_v(i,j,k,nbl), i = 1,NInew(nbl)), j = 1,NJnew(nbl)),k=1,NKnew(nbl)) 
		END DO
		close(7)

!*******************************************************END*****************************************
END PROGRAM main


