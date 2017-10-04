  subroutine track_scrf(new_t, new_E, old_t, old_E, rf_L, rf_V, rf_phase, rf_coeff1, rf_coeff2, n_size)
    implicit none
    double precision, intent(in) :: old_t, old_E, rf_L, rf_V, rf_phase
    double precision, intent(in), dimension(n_size)  :: rf_coeff1, rf_coeff2
    integer, intent(in) :: n_size
    double precision, intent(out) :: new_t, new_E
    integer, parameter :: n_step=10000
    double precision, parameter :: pi=3.14159265359, mass=938.27231E6, omega=4084070449.67, cLight=299792458
    double precision :: dz, Ez
    integer i, j
    
    new_t = old_t
    new_E = old_E
    dz = rf_L/n_step
    
    new_t = new_t + 1.0/sqrt(1.0-mass*mass/( (new_E+mass)*(new_E+mass)) )/cLight* 0.5*dz
    do i=1, n_step
      Ez = 0.5*rf_coeff1(1)
      do j=2, n_size
        Ez = Ez + rf_coeff1(j)*cos((j-1)*2*pi*(dz*i-0.5*dz-0.5*rf_L)/rf_L) + &
                  rf_coeff2(j)*sin((j-1)*2*pi*(dz*i-0.5*dz-0.5*rf_L)/rf_L)
      end do
      new_E = new_E - rf_V*Ez* cos(omega*new_t + rf_phase) *dz
      new_t = new_t + 1.0/sqrt(1.0-mass*mass/( (new_E+mass)*(new_E+mass)) )/cLight* dz
      
    end do
    new_t = new_t + 1.0/sqrt(1.0-mass*mass/( (new_E+mass)*(new_E+mass)) )/cLight* 0.5*dz
  end subroutine track_scrf
    
