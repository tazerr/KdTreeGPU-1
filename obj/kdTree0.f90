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

program mainn

	use call_kd_tree
	implicit none
	
  	integer(c_int) :: numPoints, numResults, rootIdx
  	integer :: numDimensions = 3
  	integer :: sizeOfarray 
  	real(c_float), allocatable :: coordinates(:), results(:)
  	real(c_float), allocatable :: query(:)
  	integer(c_int) :: numQuerys
  	integer :: nblocks, nbl, i, j, k, index
  	integer, allocatable :: NI(:), NJ(:), NK(:)
  	real, allocatable :: xgrid(:,:,:,:), ygrid(:,:,:,:), zgrid(:,:,:,:)
  	character(100) :: filename="grid.xyz"
	real :: lx, ly, lz

  	! Open the file for reading
  	open(unit=1, file=trim(filename), form='unformatted', status='old', access='sequential')

  	! Read the number of blocks
  	read(1) nblocks
  	write(*,*) "Block:", nblocks

  	! Allocate arrays to store NI, NJ, NK
  	allocate(NI(nblocks), NJ(nblocks), NK(nblocks))

  	! Read the values of NI, NJ, NK for each block
  	do nbl = 1, nblocks
    	read(1) NI(nbl), NJ(nbl), NK(nbl)
  	end do

	NI(1) = 64
	NJ(1) = 32
	NK(1) = 32

	lx = 10.0
	ly = 10.0
	lz = 10.0

	! Allocate arrays to store xgrid, ygrid, zgrid
  	allocate(xgrid(maxval(NI), maxval(NJ), maxval(NK), nblocks))
  	allocate(ygrid(maxval(NI), maxval(NJ), maxval(NK), nblocks))
  	allocate(zgrid(maxval(NI), maxval(NJ), maxval(NK), nblocks))

	do nbl=1,nblocks
		do k=1,NK(nbl)
			do j=1,NJ(nbl)
				do i=1,NI(nbl)
					xgrid(i,j,k,nbl) = i !((i)/(NI(nbl)-1.0)*lx)*100   !((j-1)*32.0+i + k*0.128)*22/7
					!(i-1)/(NI(nbl)-1.0)*lx
					ygrid(i,j,k,nbl) = j !((j)/(NJ(nbl)-1.0)*ly)*100
					zgrid(i,j,k,nbl) = k !((k)/(NK(nbl)-1.0)*lz)*100
				end do
			end do
		end do
	end do

  	
  	! Read the xgrid, ygrid, zgrid values for each block
  	!do nbl = 1, nblocks
    !	read(1) (((xgrid(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))    &
    !		,              (((ygrid(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))    &
	!	,              (((zgrid(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))

  	!end do

  	! Close the file
 	close(1)
  

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
  
	! TILL HERE THE GGRID FILE IS READ AND THE XYZ COORDINATES ARE
  	! STORED IN THERE RESPECTIVE ARRAYS
  	! NOW CONVERT IT TO THE TYPE OF ARRAY THE CUDA CODE IS EXCEPTING
  	! I.E. {{x1,y1,z1},{x2,y2,z2}....}
  	
  	! integer(c_int), value :: numPoints
  	! Initialize the total number of points for the KDtree
  	numPoints = size(xgrid)
  	
  	! integer, parameter :: numDimensions = 3
  	! integer, parameter :: sizeOfarray = numPoints*numDimensions
  	sizeOfarray = numPoints*numDimensions
  	
  	! real(c_float), dimension(3) :: query
  	! Initialize the query point [x,y,z]
	numQuerys=3
	allocate(query(numQuerys*numResults))
  	query = [0.0, 0.0, 0.0, 100.0, 100.0, 100.0, 15.0, 0.0, 1.0]

  	
  	! Build the coordinate array which the CUDA function is expecting {{x1,y1,z1},{x2,y2,z2}....}
  	! real(c_float), dimension(sizeOfarray) :: coordinates
  	allocate(coordinates(sizeOfarray))
  	
  	do nbl = 1, nblocks
  		do k = 1, NK(nbl)
  			do j = 1, NJ(nbl)
  				do i = 1, NI(nbl)
  				
  				! Get the current index
                index = (nbl-1)*NI(nbl)*NJ(nbl)*NK(nbl)*numDimensions + &
        			(k-1)*NI(nbl)*NJ(nbl)*numDimensions + &
        			(j-1)*NI(nbl)*numDimensions + &
        			(i-1)*numDimensions + 1

                
                !store in appropriate indices
                coordinates(index) = xgrid(i,j,k,nbl)
                coordinates(index+1) = ygrid(i,j,k,nbl)
                coordinates(index+2) = zgrid(i,j,k,nbl)
                
  				end do
  			end do
  		end do
  	end do

	! Prompt the user to enter the number of points closest points required
  	write(*,*) "Enter the number of points (integer) required:"
  	read(*,*) numResults

	allocate(results(numResults*numDimensions*numQuerys))
  	
  	call custom_funct(coordinates, numPoints, numQuerys, query, results, numResults, rootIdx)
	write(*, *) "ROOT: ", rootIdx

	call search_funct(coordinates, numDimensions, numQuerys, query, results, numResults, rootIdx)

	i=1
  	j=1
  	
	write(*,*) "Closest points points found are as follows : "
  	do while (i<=size(results))
    	write(*, *) "X (", j, ") : ", results(i)
    	write(*, *) "Y (", j, ") : ", results(i+1)
    	write(*, *) "Z (", j, ") : ", results(i+2)
    	j=j+1
    	i=i+3
  	end do

  	
  	! Deallocate arrays
  	deallocate(NI, NJ, NK)
  	deallocate(xgrid, ygrid, zgrid)
  	deallocate(coordinates)
	deallocate(results)
	deallocate(query)
end program mainn
