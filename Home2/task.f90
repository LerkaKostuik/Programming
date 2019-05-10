module task
   use :: mpi
   implicit none
   contains
       
      subroutine FindMaxCoordinates(A, x1, y1, x2, y2)
            
         real(8), intent(in), dimension(:,:) :: A
         integer(4), intent(out) :: x1, y1, x2, y2
         integer(4) :: n, L, R, Up, Down, m, tmp
         real(8), allocatable :: current_column(:), B(:,:)
         real(8) :: current_sum, max_sum
         logical :: transpos
         integer :: mpierr, mpiSize, mpiRank, location_of_maximum
         real(8), allocatable, dimension(:) :: most_of_max_sum 
         m = size(A, dim=1) 
         n = size(A, dim=2) 
         transpos = .FALSE.
           
         if (m < n) then 
             transpos = .TRUE.   
             B = transpose(A)
             m = size(B, dim=1) 
             n = size(B, dim=2) 
         else
             B = A     
         endif
       
         allocate(current_column(m))
        
         call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
         call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
         allocate(most_of_max_sum(mpiSize))
           
         x1=1
         y1=1
         x2=1
         y2=1
         max_sum=B(1,1)   
        
         do L=mpiRank+1, n, mpiSize
             current_column = B(:, L)            
             do R=L,n
       
                 if (R > L) then  
                     current_column = current_column + B(:, R)
                 endif 
               
                 call FindMaxInArray(current_column, current_sum, Up, Down) 

                 if (current_sum > max_sum) then
                     max_sum = current_sum
                     x1 = Up
                     x2 = Down
                     y1 = L
                     y2 = R
                 endif
                
             end do
         end do
  
         call mpi_gather(max_sum,1,MPI_REAL8,most_of_max_sum,1,MPI_REAL8,0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(most_of_max_sum,mpiSize,MPI_REAL8, 0, MPI_COMM_WORLD, mpierr)
         location_of_maximum=maxloc(most_of_max_sum,dim=1)-1
         
         
         call mpi_bcast(x1, 1, MPI_INTEGER, location_of_maximum, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(x2, 1, MPI_INTEGER, location_of_maximum, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(y1, 1, MPI_INTEGER, location_of_maximum, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(y2, 1, MPI_INTEGER, location_of_maximum, MPI_COMM_WORLD, mpierr)

         deallocate(current_column)
         deallocate(most_of_max_sum)
         if (transpos) then  
             tmp = x1
             x1 = y1
             y1 = tmp
    
             tmp = y2
             y2 = x2
             x2 = tmp
         endif

      end subroutine


      subroutine FindMaxInArray(a, Summ, Up, Down)
            
          real(8), intent(in), dimension(:) :: a
          integer(4), intent(out) :: Up, Down
          real(8), intent(out) :: Summ
          real(8) :: cur_sum
          integer(4) :: minus_pos, i

          Summ = a(1)
          Up = 1
          Down = 1
          cur_sum = 0
          minus_pos = 0
           
          do i=1, size(a)
              cur_sum = cur_sum + a(i)
              if (cur_sum > Summ) then
                  Summ = cur_sum
                  Up = minus_pos + 1
                  Down = i
              endif
         
              if (cur_sum < 0) then
                  cur_sum = 0
                  minus_pos = i
              endif

          enddo

      end subroutine FindMaxInArray

end module task
