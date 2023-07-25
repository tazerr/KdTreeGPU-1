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

		SUBROUTINE search_funct(coordinates, numDimensions, numQuerys, query, results, numResults, rootIdx &
            , mltip, gindices, dists) bind(C, name="search_funct")

            !@@$acc routine seq
            USE iso_c_binding
            REAL(c_float), dimension(*) :: coordinates
            INTEGER(c_int), value :: numDimensions
            INTEGER(c_int), value :: numQuerys
            REAL(c_float), dimension(*) :: query
			REAL(c_float), dimension(*) :: results
			INTEGER(c_int), value :: numResults
			INTEGER(c_int), value:: rootIdx
            REAL(c_float), value :: mltip
            INTEGER(c_int), dimension(*) :: gindices
            REAL(c_double), dimension(*) :: dists
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
INTEGER :: numDimensions 
INTEGER :: sizeOfarray 
REAL(c_float):: mltip=100000.0
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

numDimensions=3
numPoints = size(xgrid)
sizeOfarray = numPoints*numDimensions
numQuerys = size(Xgridnew)
numResults = 1
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
ALLOCATE(gindices(numQuerys*numResults))
ALLOCATE(dists(numQuerys*numResults))

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

CALL search_funct(coordinates, numDimensions, numQuerys, query, results, numResults, rootIdx, mltip &
    , gindices, dists) 	

!***************************************Print Result******************************************

!i=1
!j=1
!WRITE(*,*) "Closest points points found are as follows : "
!DO WHILE (i<=size(results))
!	WRITE(*, *) "X (", j, ") : ", results(i)
!	WRITE(*, *) "Y (", j, ") : ", results(i+1)
!	WRITE(*, *) "Z (", j, ") : ", results(i+2)
!	j=j+1
!	i=i+3
!END DO

!**************************************KD Tree Search******************************************
interptsx = 5
interptsy = 3
interplace = 2

!ALLOCATE(query(numQuerys*numDimensions))
!ALLOCATE(gindices(numQuerys*numResults))
!ALLOCATE(dists(numQuerys*numResults))

ALLOCATE (xstencil(interptsx,interptsy))
ALLOCATE (ystencil(interptsx,interptsy))
ALLOCATE (zstencil(interptsx,interptsy))
ALLOCATE (istencil(interptsx))
ALLOCATE (jstencil(interptsy))

!$acc data copyin (Xgrid, Ygrid, Zgrid, Xgridnew, Ygridnew, Zgridnew)
!$acc data copyin (Igrid, Jgrid, Igridnew, Jgridnew)
!$acc data copyin (coordinates, results, query, gindices, dists)
!$acc data copyin (istencil, jstencil, xstencil, ystencil, zstencil)

	DO nbl = 1,nblocksnew
!$acc parallel loop collapse(3) private(istencil,jstencil,xstencil,ystencil,zstencil)
    DO k = 1, NKnew(nbl)
    DO j = 1, NJnew(nbl)
    DO i = 1, NInew(nbl)

        !query(1) = xgridnew(i,j,k,nbl)
		!query(2) = ygridnew(i,j,k,nbl)
		!query(3) = zgridnew(i,j,k,nbl)

        !CALL search_funct(coordinates, numDimensions, numQuerys, query, results &
        !    , numResults, rootIdx, mltip, gindices, dists)

        ! Get the current index
        index = (nbl-1)*NInew(nbl)*NJnew(nbl)*NKnew(nbl)+ &
            (k-1)*NInew(nbl)*NJnew(nbl)+ &
            (j-1)*NInew(nbl)+ &
            (i-1) + 1
            
        DO indexing = 1,numResults
			convertindex = gindices((index-1)*numResults+indexing)
			dw = (dists((index-1)*numResults+indexing))**0.5
			call get_loc_index(convertindex,loc_blk, loc_i, loc_j, loc_k)
            DO indx = 1,interptsx
				istencil(indx) = Igrid(loc_i + indx - interplace, loc_j, loc_k,loc_blk)
				DO indy = 1,interptsy
					xstencil(indx,indy) = Xgrid(loc_i + indx - interplace, loc_j+ indy - interplace, loc_k,loc_blk)
					ystencil(indx,indy) = Ygrid(loc_i + indx - interplace, loc_j+ indy - interplace, loc_k,loc_blk)
					zstencil(indx,indy) = Zgrid(loc_i + indx - interplace, loc_j+ indy - interplace, loc_k,loc_blk)
					jstencil(indy) = Jgrid(loc_i, loc_j + indy - interplace, loc_k,loc_blk)
				END DO
			END DO
            iinter = Igrid(loc_i, loc_j, loc_k,loc_blk)
			jinter = Jgrid(loc_i, loc_j, loc_k,loc_blk)
            print*, zstencil
			call Coord_Interpolation(interptsx,interptsy,xstencil,ystencil,zstencil,istencil,jstencil &
                ,query((index-1)*numDimensions+1),query((index-1)*numDimensions+2),query((index-1)*numDimensions+3),iinter,jinter)
			Igridnew(i,j,k,nbl) = iinter
			Jgridnew(i,j,k,nbl) = jinter
        
			print*, istencil
			print*, jstencil
			print*, iinter, jinter
			pause
        END DO
    END DO
	END DO
	END DO
	END DO

!$acc end data
!$acc end data
!$acc end data
!$acc end data



!DO i=1, numQuerys*numResults
!    CALL get_loc_index(gindices(i), loc_blk, loc_i, loc_j, loc_k)
!    write(*,*) "xgrid(", i, ",", dists(i), ") =", xgrid(loc_i, loc_j, loc_k, loc_blk)
!    write(*,*) "ygrid(", i, ",", dists(i), ") =", ygrid(loc_i, loc_j, loc_k, loc_blk) 
!    write(*,*) "zgrid(", i, ",", dists(i), ") =", zgrid(loc_i, loc_j, loc_k, loc_blk)
!END DO



!*******************************************************END*****************************************
END PROGRAM main

!*************************************Local Index subroutine***************************************
SUBROUTINE get_loc_index(global_index,blk,iloc,jloc,kloc)
    !@@@$acc routine seq
    USE arrs

    INTEGER :: global_index, blk, iloc, jloc, kloc
    INTEGER :: rest_index
    INTEGER :: sblk,eblk
    sblk=0

    DO nbl = 1,1
        eblk=256*256*1
        !eblk=NI(nbl)*NJ(nbl)*NK(nbl)
        IF ((global_index>sblk).and.(global_index<=eblk)) THEN 
            blk = nbl
            EXIT
        END IF
        sblk=eblk
    END DO

    rest_index = global_index - sblk
    !kloc = int(rest_index/(NJ(blk)*NI(blk)))
    kloc = int(rest_index/(256*256))
    !rest_index = rest_index - kloc*NI(nbl)*NJ(nbl)
    rest_index = rest_index - kloc*256*256
    jloc = int(rest_index/(256))
    !jloc = int(rest_index/(NI(nbl)))
    !rest_index = rest_index - jloc*NI(blk)
    rest_index = rest_index - jloc*256
    iloc = int(rest_index)
    
    iloc = iloc+1
    jloc = jloc+1
    kloc = kloc+1
        
END 

!*************************************I-J prediction Subroutine********************************************
SUBROUTINE Coord_Interpolation(nxint, nyint, xgint, ygint, zgint, igint, jgint, xpin, ypin, zpin, ipred, jpred)

!@@@@$acc routine seq
USE iso_c_binding
implicit none

integer nxint, nyint, n, iter_int
real,dimension (nxint, nyint) :: xgint, ygint, zgint
real,dimension (nxint) :: igint
real,dimension (nyint) :: jgint
real, allocatable :: iweight(:), jweight(:)
integer lint, mint, llint, mmint
real ipred, jpred, itemp, jtemp
real(c_float) xpin, ypin, zpin
real xhat, yhat,zhat, xihat, yihat, xjhat, yjhat, zihat, zjhat
real jac_int(3,2), hat_coord(3), intA(2,2) , RHS_int(2)
real detA, invA(2,2), mind, consa, consb
real, allocatable :: intweights(:,:), diffweightsi(:,:),diffweightsj(:,:),weighti(:) &
    , weightj(:), weightdiffi(:), weightdiffj(:)

mind = 0.00000001d0

!allocate (xgint(nxint,nyint))
!allocate (ygint(nxint,nyint))
!allocate (zgint(nxint,nyint))
!allocate (igint(nxint))
!allocate (jgint(nyint))

allocate (iweight(nxint))
allocate (jweight(nyint))
allocate (intweights(nxint,nyint))
allocate (weightdiffi(nxint))
allocate (weightdiffj(nyint))
allocate (diffweightsi(nxint,nyint))
allocate (diffweightsj(nxint,nyint))

!***********************************Interpolation coefficients**************************

!print*, ygint


itemp = ipred + 1.d0
jtemp = jpred + 1.d0 
iter_int=0
!print*, ipred, jpred
Do while ((abs(ipred-itemp).gt.mind).or.(abs(jpred-jtemp).gt.mind))
 iter_int= iter_int+1
! print*,iter_int
itemp = ipred
jtemp = jpred	

Do lint = 1,nxint
  iweight(lint) = 1.0D+00
    do mint = 1, nxint
      if ( mint /= lint ) then
        iweight(lint) = iweight(lint) * (itemp-igint(mint)) / ( igint(lint) - igint(mint) )
      end if
    end do
Enddo

Do lint = 1,nyint
  jweight(lint) = 1.0D+00
    do mint = 1, nyint
      if ( mint /= lint ) then
        jweight(lint) = jweight(lint) * (jtemp-jgint(mint)) / ( jgint(lint) - jgint(mint) )
      end if
    end do
Enddo	

xhat = 0.d0
yhat = 0.d0
zhat = 0.d0
Do lint = 1,nxint
	Do mint = 1,nyint
		xhat = xhat + xgint(lint,mint) * iweight(lint) * jweight(mint)
		yhat = yhat + ygint(lint,mint) * iweight(lint) * jweight(mint)
		zhat = zhat + zgint(lint,mint) * iweight(lint) * jweight(mint)
	enddo
Enddo

!print*, xhat, yhat, zhat

Do lint = 1,nxint
  weightdiffi(lint) = 0.d0
    do mint = 1, nxint
      if ( mint /= lint ) then
		consa = 1.d0
		do n = 1,nxint
			if(n/=lint) then
				consa = consa/( igint(lint) - igint(n) )
			endif
		enddo
		consb = 1.d0
		do n = 1,nxint
			if((n/=lint).and.(n/=mint)) then
				consb = consb * (itemp-igint(n))
			endif
		enddo
 !       weightdiffi(lint) = weightdiffi(lint) + iweight(lint)/(itemp-igrid(mint))
		weightdiffi(lint) = weightdiffi(lint) + consb*consa
      end if
    end do
Enddo

Do lint = 1,nyint
  weightdiffj(lint) = 0.d0
    do mint = 1, nyint
      if ( mint /= lint ) then
		consa = 1.d0
		do n = 1,nyint
			if(n/=lint) then
				consa = consa/( jgint(lint) - jgint(n) )
			endif
		enddo
		consb = 1.d0
		do n = 1,nyint
			if((n/=lint).and.(n/=mint)) then
				consb = consb * (jtemp-jgint(n))
			endif
		enddo
 !       weightdiffi(lint) = weightdiffi(lint) + iweight(lint)/(itemp-igrid(mint))
		weightdiffj(lint) = weightdiffj(lint) + consb*consa
      end if
    end do
Enddo

!print*, iweight
!print*, sum(iweight)

Do lint = 1,nxint
	Do mint = 1,nyint
		diffweightsi(lint,mint) = weightdiffi(lint) * jweight(mint)
		diffweightsj(lint,mint) = iweight(lint) * weightdiffj(mint)
	enddo
Enddo

xihat = 0.d0
yihat = 0.d0
Do lint = 1,nxint
	Do mint = 1,nyint
		xihat = xihat + diffweightsi(lint,mint) * xgint(lint,mint)
		yihat = yihat + diffweightsi(lint,mint) * ygint(lint,mint)
		xjhat = xjhat + diffweightsj(lint,mint) * xgint(lint,mint)
		yjhat = yjhat + diffweightsj(lint,mint) * ygint(lint,mint)
		zihat = zihat + diffweightsi(lint,mint) * zgint(lint,mint)
		zjhat = zjhat + diffweightsj(lint,mint) * zgint(lint,mint)
	enddo
Enddo

!print*, xihat, yihat, zihat
!print*, xjhat, yjhat, zjhat


jac_int(1,1) = xihat
jac_int(1,2) = xjhat
jac_int(2,1) = yihat
jac_int(2,2) = yjhat
jac_int(3,1) = zihat
jac_int(3,2) = zjhat

hat_coord(1) = xhat-xpin
hat_coord(2) = yhat-ypin
hat_coord(3) = zhat-zpin

intA = matmul(transpose(jac_int),jac_int)
RHS_int = matmul(transpose(jac_int),hat_coord)

detA = intA(2,2) *intA(1,1) - intA(1,2)*intA(2,1)

invA(1,1) = intA(2,2)/detA
invA(1,2) = -1.d0 * intA(2,1)/detA
invA(2,1) = -1.d0 * intA(1,2)/detA
invA(2,2) = intA(1,1)/detA

ipred = itemp - (invA(1,1)*RHS_int(1) + invA(1,2)*RHS_int(2))
jpred = jtemp - (invA(2,1)*RHS_int(1) + invA(2,2)*RHS_int(2))

!print*, ipred
!print*, jpred

Enddo

END SUBROUTINE
