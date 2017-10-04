subroutine differetial_evolution_f90(obj_func, x_out, fval, n_dim ,x_min, x_max, &
																		pop_size, mutation, crossover, strategy, itermaxer_max)

	implicit none
	include "mpif.h"
	external  obj_func
	integer, intent(in) :: n_dim, pop_size, strategy, iter_max
	double precision, intent(in) :: mutation, crossover
	double precision, intent(in), dimension(n_dim) :: x_min, x_max
	
	double precision, intent(out), dimension(n_dim) :: x_out
	double precision, intent(out) :: fval
	
  integer :: i, pop_size, strategy, iter, ierr, num_procs, my_id, n_slice
	double precision :: mutation, crossover
	double precision, dimension(pop_size*n_dim, n_dim) :: x_pop
	double precision, dimension(n_dim) :: d_dum1, d_dum2

	vim

.........................................................................
!             Dim_XC : Dimension of the real decision parameters.
!      XCmin(Dim_XC) : The lower bound of the real decision parameters.
!      XCmax(Dim_XC) : The upper bound of the real decision parameters.
!                VTR : The expected fitness value to reach.
!                 NP : Population size.
!            itermax : The maximum number of iteration.
!               F_XC : Mutation scaling factor for real decision parameters.
!              CR_XC : Crossover factor for real decision parameters.
!           strategy : The strategy of the mutation operations is used in HDE.
!            refresh : The intermediate output will be produced after "refresh"
!                      iterations. No intermediate output will be produced if
!                      "refresh < 1".
!             iwrite : The unit specfier for writing to an external data file.
! bestmen_XC(Dim_XC) : The best real decision parameters.
!              bestval : The best objective function.
!             nfeval : The number of function call.
!         method(1) = 0, Fixed mutation scaling factors (F_XC)
!                   = 1, Random mutation scaling factors F_XC=[0, 1]
!                   = 2, Random mutation scaling factors F_XC=[-1, 1] 
!         method(2) = 1, Random combined factor (F_CR) used for strategy = 6
!                        in the mutation operation 
!                   = other, fixed combined factor provided by the user 
!         method(3) = 1, Saving results in a data file.
!                   = other, displaying results only.


	 real(kind=8), dimension(Dim_XC), intent(in) :: XCmin, XCmax
     real(kind=8), dimension(Dim_XC), intent(inout) :: bestmem_XC 
     real(kind=8), dimension(NP,Dim_XC) :: pop_XC, bm_XC, mui_XC, mpo_XC,   &
	                                        popold_XC, rand_XC, ui_XC
     integer(kind=4), dimension(NP) :: rot, a1, a2, a3, a4, a5, rt
     integer(kind=4), dimension(4) :: ind
	 real(kind=8), dimension(NP) :: val
     real(kind=8), dimension(Dim_XC) :: bestmemit_XC
     real(kind=8), dimension(Dim_XC) :: rand_C1
	 integer(kind=4), dimension(3), intent(in) :: method

 !!-----Initialize a population --------------------------------------------!!
	dummy1 = (x_max-x_min)/pop_size
	call random_number(dummy2)
	do i=1,pop_size
		do j=1, n_dim
			x_pop(i+(j-1)*pop_size,:)=x_min +dummy1*(i-1+dummy2)
	end do
 
!!--------------------------------------------------------------------------!!

!!------Evaluate fitness functions and find the best member-----------------!!
	fval = 0.0d
	iter = 0
	call MPI_INIT(ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
	n_slice = (n_dim*pop_size / num_procs) + 1
	if ( my_id == 0 ) then 
		do i=1, num_procs-1
			! distribute population
			call MPI_SEND( x_pop( (i-1)*n_slice+1,: ), n_slice, MPI_DOUBLE_PRECISION, &
										 i, 1000, MPI_COMM_WORLD, ierr)
			pint *, ierr
		end do
		!do i=(num_procs-1)*n_slice+1, n_dim*pop_size
		!	f_pop(i)=obj_func(x_pop(i),my_id)
		!end do
	else
		call MPI_RECV( x_pop, n_slice, MPI_DOUBLE_PRECISION, 0, &
                   1000, MPI_COMM_WORLD, status, ierr)
		pint *, my_id, "th processor, err = ", ierr
		pint *, my_id, "th processor, xpop = ", x_pop
		
end subroutine differetial_evolution_f90