module call_kd_tree
    use iso_c_binding
    implicit none
    interface
        subroutine custom_funct(coordinates, numPoints, searchDistance, query, results, numResults) bind(C, name="custom_funct")
            use iso_c_binding
            real(c_float), dimension(*) :: coordinates
            integer(c_int), value :: numPoints
            integer(c_int), value :: searchDistance
            real(c_float), dimension(*) :: query
			real(c_float), dimension(*) :: results
			integer(c_int), value :: numResults
        end subroutine
    end interface
end module call_kd_tree

program mainn

	use call_kd_tree
	implicit none
	! integer :: nblocks, nbl, i, j, k, i1, i2, i3, i4, i5, index
  	! integer, allocatable :: NI(:), NJ(:), NK(:)
  	! real, allocatable :: xgrid(:,:,:,:), ygrid(:,:,:,:), zgrid(:,:,:,:)
  	! character(50) :: filename = "grid.xyz"
  	integer(c_int) :: numPoints, numResults
  	integer :: numDimension = 3
  	integer :: sizeOfarray 
  	real(c_float), allocatable :: coordinates(:), results(:)
  	real(c_float), dimension(3) :: query
  	integer(c_int) :: searchDistance
  	integer :: nblocks, nbl, i, j, k, index
  	integer, allocatable :: NI(:), NJ(:), NK(:)
  	real, allocatable :: xgrid(:,:,:,:), ygrid(:,:,:,:), zgrid(:,:,:,:)
  	character(100) :: filename="grid.xyz"
	real :: lx, ly

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

	NI(1) = 32
	NJ(1) = 32
	NK(1) = 1

	lx = 10.0
	ly = 10.0

	! Allocate arrays to store xgrid, ygrid, zgrid
  	allocate(xgrid(maxval(NI), maxval(NJ), maxval(NK), nblocks))
  	allocate(ygrid(maxval(NI), maxval(NJ), maxval(NK), nblocks))
  	allocate(zgrid(maxval(NI), maxval(NJ), maxval(NK), nblocks))

	do nbl=1,nblocks
		do k=1,NK(nbl)
			do j=1,NJ(nbl)
				do i=1,NI(nbl)
					xgrid(i,j,k,nbl) = ((j-1)*32+i)*22/7
					!(i-1)/(NI(nbl)-1.0)*lx
					ygrid(i,j,k,nbl) = (j-1)/(NJ(nbl)-1.0)*ly
					zgrid(i,j,k,nbl) = 0.0
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
  	do nbl = 1, nblocks
    	write(*,*) "Block:", nbl
    	do k = 1, NK(nbl)
    	  do j = 1, NJ(nbl)
    	    do i = 1, NI(nbl)
    	      write(*,*) "xgrid(", i, j, k, nbl, ") =", xgrid(i, j, k, nbl)
    	      write(*,*) "ygrid(", i, j, k, nbl, ") =", ygrid(i, j, k, nbl)
    	      write(*,*) "zgrid(", i, j, k, nbl, ") =", zgrid(i, j, k, nbl)
    	    end do
    	  end do
    	end do
  	end do
  
	! TILL HERE THE GGRID FILE IS READ AND THE XYZ COORDINATES ARE
  	! STORED IN THERE RESPECTIVE ARRAYS
  	! NOW CONVERT IT TO THE TYPE OF ARRAY THE CUDA CODE IS EXCEPTING
  	! I.E. {{x1,y1,z1},{x2,y2,z2}....}
  	
  	! integer(c_int), value :: numPoints
  	! Initialize the total number of points for the KDtree
  	numPoints = size(xgrid)
  	
  	! integer, parameter :: numDimension = 3
  	! integer, parameter :: sizeOfarray = numPoints*numDimension
  	sizeOfarray = numPoints*numDimension
  	
  	! real(c_float), dimension(3) :: query
  	! Initialize the query point [x,y,z]
  	query = [0,0,0]
  	
  	! integer(c_int), value :: searchDistance
  	! Initialize the value of the search distance from the queryPoints
  	searchDistance = 2
  	
  	! Build the coordinate array which the CUDA function is expecting {{x1,y1,z1},{x2,y2,z2}....}
  	! real(c_float), dimension(sizeOfarray) :: coordinates
  	allocate(coordinates(sizeOfarray))
  	
  	do nbl = 1, nblocks
  		do k = 1, NK(nbl)
  			do j = 1, NJ(nbl)
  				do i = 1, NI(nbl)
  				
  				! Get the current index
                index = (nbl-1)*NI(nbl)*NJ(nbl)*NK(nbl)*numDimension + &
        			(k-1)*NI(nbl)*NJ(nbl)*numDimension + &
        			(j-1)*NI(nbl)*numDimension + &
        			(i-1)*numDimension + 1

                
                !store in appropriate indices
                coordinates(index) = xgrid(i,j,k,nbl)
                coordinates(index+1) = ygrid(i,j,k,nbl)
                coordinates(index+2) = zgrid(i,j,k,nbl)
                
  				end do
  			end do
  		end do
  	end do
  	
  	! i=1
  	! j=1
  	
  	! do while (i<=size(coordinates))
    !	write(*, *) "X (", j, ") : ", coordinates(i)
    !	write(*, *) "Y (", j, ") : ", coordinates(i+1)
    !	write(*, *) "Z (", j, ") : ", coordinates(i+2)
    !	j=j+1
    !	i=i+3
  	! end do

	! Prompt the user to enter the number of points closest points required
  	write(*,*) "Enter the number of points (integer) required:"
  	read(*,*) numResults

	allocate(results(numResults*numDimension))
  	
  	call custom_funct(coordinates, numPoints, searchDistance, query, results, numResults)
  	
  	! Deallocate arrays
  	deallocate(NI, NJ, NK)
  	deallocate(xgrid, ygrid, zgrid)
  	deallocate(coordinates)
end program mainn
