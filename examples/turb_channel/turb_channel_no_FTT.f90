! Martin Karp 13/3-2023
! updated initial condition Philipp Schlatter 09/07/2024
module user
  use neko
  ! TorchFort
  use bc, only: bc_t
  use wall_model_bc, only: wall_model_bc_t
  use wall_model, only: wall_model_t
  use spalding, only: spalding_t
  use field_registry, only : neko_field_registry
  use math
  implicit none

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_user_ic => user_ic
    u%user_mesh_setup => user_mesh_scale
    u%user_init_modules => initialize ! TorchFort
    u%user_check => usercheck ! TorchFort
  end subroutine user_setup

  ! Initialize user variables or external objects
  subroutine initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    ! insert your initialization code here
    logical :: found

    call neko_field_registry%add_field(coef%dof, "u_old")
    found = neko_field_registry%field_exists("u_old")
    print *, "u_old field_exists: ", found

    call neko_field_registry%add_field(coef%dof, "u_older")
    found = neko_field_registry%field_exists("u_older")
    print *, "u_older field_exists: ", found

  end subroutine initialize

  ! Rescale mesh, we create a mesh with some refinement close to the wall.
  ! initial mesh: 0..4, -1..1, 0..1.5
  ! mesh size (4*pi,2*delta,4/3*pi)
  ! New mesh can easily be genreated with genmeshbox
  ! OBS refinement is not smooth and the constant values are a bit ad hoc.
  ! Stats converge close to reference DNS
  subroutine user_mesh_scale(msh)
    type(mesh_t), intent(inout) :: msh
    integer :: i, p, nvert

    real(kind=rp) :: d, y, viscous_layer, visc_el_h, el_h
    real(kind=rp) :: center_el_h, dist_from_wall
    integer :: el_in_visc_lay, el_in_y
    real(kind=rp) :: llx, llz

    ! target mesh size
    llx = 4.*pi
    llz = 4./3.*pi

    ! rescale mesh
    el_in_y = 18
    el_in_visc_lay = 2
    viscous_layer = 0.0888889
    el_h = 2.0_rp/el_in_y
    visc_el_h = viscous_layer/el_in_visc_lay
    center_el_h = (1.0_rp-viscous_layer)/(el_in_y/2-el_in_visc_lay)

    nvert = size(msh%points)
    do i = 1, nvert
       msh%points(i)%x(1) = llx/4.*msh%points(i)%x(1)
       y = msh%points(i)%x(2)
       if ((1-abs(y)) .le. (el_in_visc_lay*el_h)) then
          dist_from_wall = (1-abs(y))/el_h*visc_el_h
       else
          dist_from_wall = viscous_layer + (1-abs(y)- &
               el_in_visc_lay*el_h)/el_h*center_el_h
       end if
       if (y .gt. 0) msh%points(i)%x(2) = 1.0_rp - dist_from_wall
       if (y .lt. 0) msh%points(i)%x(2) = -1.0_rp + dist_from_wall
       msh%points(i)%x(3) = 2./3.*llz*msh%points(i)%x(3)
    end do

  end subroutine user_mesh_scale

  ! User defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i

    real(kind=rp) :: uvw(3)

    do i = 1, u%dof%size()
       uvw = channel_ic(u%dof%x(i,1,1,1),u%dof%y(i,1,1,1),u%dof%z(i,1,1,1))
       u%x(i,1,1,1) = uvw(1)
       v%x(i,1,1,1) = uvw(2)
       w%x(i,1,1,1) = uvw(3)
    end do
  end subroutine user_ic

  ! This is called at the end of every time step
  subroutine usercheck(t, tstep, u, v, w, p, coef, param)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: param

    ! insert code below
    type(field_t), pointer :: u_old, u_older, dudy_old
    integer :: i, j, k, e, n, i_wm
    class(bc_t), pointer :: bc
    type(wall_model_bc_t), pointer :: wall_bc
    type(field_t), pointer :: state_old, state_older, action_old, action_older, terminal_old, terminal_older

    u_old => neko_field_registry%get_field("u_old")
    u_older => neko_field_registry%get_field("u_older")

!     state_old => neko_field_registry%get_field("state_old")
!     state_older => neko_field_registry%get_field("state_older")

!     action_old => neko_field_registry%get_field("action_old")
!     action_older => neko_field_registry%get_field("action_older")

    terminal_old => neko_field_registry%get_field("terminal_old")
    terminal_older => neko_field_registry%get_field("terminal_older")

    print *, "***** USERCHECK *****"
!     do e = 1, 2
!       do k = 1, 2
!         do j = 1, 2
!           do i = 1, 2
!             print *, u_older%x(i,j,k,e), u_old%x(i,j,k,e), u%x(i,j,k,e)
!           end do
!         end do
!       end do
!     end do
    print *, "***** USERCHECK *****"

    call copy(u_older%x, u_old%x, size(u_old%x))
    call copy(u_old%x, u%x, size(u%x)) ! coef%dof%size()

    n = neko_simcomps%case%fluid%bcs_vel%size()
    do i_wm = 1, n
      bc => neko_simcomps%case%fluid%bcs_vel%get(i_wm)
      select type (bc)
      type is (wall_model_bc_t)
        wall_bc => bc
        select type(spald => wall_bc%wall_model)
        type is (spalding_t)
          if (allocated(spald%state)) then
!             print *, "spalding state(1,1): ", spald%state(1,1)
!             print *, "shape(spald%state): ", shape(spald%state)

!             call copy(state_older%x, state_old%x, size(state_older%x))
!             call copy(state_old%x, spald%state, size(spald%state))
            call copy(spald%sor, spald%so, size(spald%so))
            call copy(spald%so, spald%state, size(spald%state))
            print *, "size(spald%so)", size(spald%so)
            print *, "size(spald%sor)", size(spald%sor)
            print *, "No. of Zeros: ", count(spald%sor == 0.0)

!             call copy(action_older%x, action_old%x, size(action_old%x))
!             call copy(action_old%x, spald%action, size(spald%action))
            call copy(spald%aor, spald%ao, size(spald%ao))
            call copy(spald%ao, spald%action, size(spald%action))
!             print *, "Total number of wall nodes = ", spald%n_nodes
!             print *, "Current time = ", t
!             print *, "End time = ", neko_simcomps%case%time%end_time

            do i = 1, spald%n_nodes

              if (abs(t - neko_simcomps%case%time%end_time) .le. 1.0e-3_rp) then
                spald%terminal(i,1) = 1.0_rp
                if (i >= 1 .and. i <= 5) then
!                   print *, spald%terminal(1,i)
                end if
              else
                spald%terminal(i,1) = 0.0_rp
                if (i >= 1 .and. i <= 5) then
!                   print *, spald%terminal(i,1)
                end if
              end if

                if (i >= 1 .and. i <= 5) then
                    if (i == 1) then
                    write(*, '(A20, F10.4)') '1st do loop: t =', t
                    write(*, '(A6, 2X, A20, 2X, A20, 2X, A20)') &
                        'i', 'so(i,1)', 'so(i,2)', 'so(i,3)'
                    end if
                    write(*, '(I6, 2X, ES20.10, 2X, ES20.10, 2X, ES20.10)') &
                    i, spald%so(i,1), spald%so(i,2), spald%so(i,3)
                end if

                if (i >= 1 .and. i <= 5) then
                    if (i == 1) then
                    write(*, '(A20, F10.4)') '1st do loop: t =', t
                    write(*, '(A6, 2X, A20, 2X, A20, 2X, A20)') &
                        'i', 'sor(i,1)', 'sor(i,2)', 'sor(i,3)'
                    end if
                    write(*, '(I6, 2X, ES20.10, 2X, ES20.10, 2X, ES20.10)') &
                    i, spald%sor(i,1), spald%sor(i,2), spald%sor(i,3)
                end if

            end do
            call copy(terminal_older%x, terminal_old%x, size(terminal_old%x))
            call copy(terminal_old%x, spald%terminal, size(spald%terminal))
          end if
        end select
      end select
    end do

  end subroutine usercheck

  ! Kind of brute force with rather large initial disturbances
  function channel_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: ux, uy, uz, eps, Re_tau, yp, Re_b, alpha, beta
    real(kind=rp) :: C, k, kx, kz, eps1, ran

    real(kind=rp) :: llx, llz

    llx = 4.*pi
    llz = 4./3.*pi

    Re_tau = 180
    C      = 5.17
    k      = 0.41
    Re_b   = 2800

    yp = (1-y)*Re_tau
    if (y .lt. 0) yp = (1+y)*Re_tau

    ! Reichardt function
    ux  = 1/k*log(1.0+k*yp) + (C - (1.0/k)*log(k)) * &
         (1.0 - exp(-yp/11.0) - yp/11*exp(-yp/3.0))
    ux  = ux * Re_tau/Re_b

    ! actually, sometimes one may not use the turbulent profile, but
    ! rather the parabolic lamianr one
    ! ux = 1.5*(1-y**2)

    ! add perturbations to trigger turbulence
    ! base flow
    uvw(1)  = ux
    uvw(2)  = 0
    uvw(3)  = 0

    ! first, large scale perturbation
    eps = 0.05
    kx  = 3
    kz  = 4
    alpha = kx * 2*PI/llx
    beta  = kz * 2*PI/llz
    uvw(1)  = uvw(1) + eps*beta  * sin(alpha*x)*cos(beta*z)
    uvw(2)  = uvw(2) + eps       * sin(alpha*x)*sin(beta*z)
    uvw(3)  = uvw(3) -eps*alpha * cos(alpha*x)*sin(beta*z)

    ! second, small scale perturbation
    eps = 0.005
    kx  = 17
    kz  = 13
    alpha = kx * 2*PI/llx
    beta  = kz * 2*PI/llz
    uvw(1)  = uvw(1) + eps*beta  * sin(alpha*x)*cos(beta*z)
    uvw(2)  = uvw(2) + eps       * sin(alpha*x)*sin(beta*z)
    uvw(3)  = uvw(3) -eps*alpha * cos(alpha*x)*sin(beta*z)

    ! finally, random perturbations only in y
    eps1 = 0.001
    ran = sin(-20*x*z+y**3*tan(x*z**2)+100*z*y-20*sin(x*y*z)**5)
    uvw(2)  = uvw(2) + eps1*ran

  end function channel_ic

end module user


! ! Martin Karp 13/3-2023
! ! updated initial condition Philipp Schlatter 09/07/2024
! module user
!   use neko
!   ! TorchFort
!   use bc, only: bc_t
!   use wall_model_bc, only: wall_model_bc_t
!   use wall_model, only: wall_model_t
!   use spalding, only: spalding_t
!   use field_registry, only : neko_field_registry
!   use math
!   implicit none
!
! contains
!
!   ! Register user defined functions (see user_intf.f90)
!   subroutine user_setup(u)
!     type(user_t), intent(inout) :: u
!     u%fluid_user_ic => user_ic
!     u%user_mesh_setup => user_mesh_scale
!     u%user_init_modules => initialize ! TorchFort
!     u%user_check => usercheck ! TorchFort
!   end subroutine user_setup
!
!   ! Initialize user variables or external objects
!   subroutine initialize(t, u, v, w, p, coef, params)
!     real(kind=rp) :: t
!     type(field_t), intent(inout) :: u
!     type(field_t), intent(inout) :: v
!     type(field_t), intent(inout) :: w
!     type(field_t), intent(inout) :: p
!     type(coef_t), intent(inout) :: coef
!     type(json_file), intent(inout) :: params
!
!     ! insert your initialization code here
!     logical :: found
!
!     call neko_field_registry%add_field(coef%dof, "u_old")
!     found = neko_field_registry%field_exists("u_old")
!     print *, "u_old field_exists: ", found
!
!     call neko_field_registry%add_field(coef%dof, "u_older")
!     found = neko_field_registry%field_exists("u_older")
!     print *, "u_older field_exists: ", found
!
!   end subroutine initialize
!
!   ! Rescale mesh, we create a mesh with some refinement close to the wall.
!   ! initial mesh: 0..4, -1..1, 0..1.5
!   ! mesh size (4*pi,2*delta,4/3*pi)
!   ! New mesh can easily be genreated with genmeshbox
!   ! OBS refinement is not smooth and the constant values are a bit ad hoc.
!   ! Stats converge close to reference DNS
!   subroutine user_mesh_scale(msh)
!     type(mesh_t), intent(inout) :: msh
!     integer :: i, p, nvert
!
!     real(kind=rp) :: d, y, viscous_layer, visc_el_h, el_h
!     real(kind=rp) :: center_el_h, dist_from_wall
!     integer :: el_in_visc_lay, el_in_y
!     real(kind=rp) :: llx, llz
!
!     ! target mesh size
!     llx = 4.*pi
!     llz = 4./3.*pi
!
!     ! rescale mesh
!     el_in_y = 18
!     el_in_visc_lay = 2
!     viscous_layer = 0.0888889
!     el_h = 2.0_rp/el_in_y
!     visc_el_h = viscous_layer/el_in_visc_lay
!     center_el_h = (1.0_rp-viscous_layer)/(el_in_y/2-el_in_visc_lay)
!
!     nvert = size(msh%points)
!     do i = 1, nvert
!        msh%points(i)%x(1) = llx/4.*msh%points(i)%x(1)
!        y = msh%points(i)%x(2)
!        if ((1-abs(y)) .le. (el_in_visc_lay*el_h)) then
!           dist_from_wall = (1-abs(y))/el_h*visc_el_h
!        else
!           dist_from_wall = viscous_layer + (1-abs(y)- &
!                el_in_visc_lay*el_h)/el_h*center_el_h
!        end if
!        if (y .gt. 0) msh%points(i)%x(2) = 1.0_rp - dist_from_wall
!        if (y .lt. 0) msh%points(i)%x(2) = -1.0_rp + dist_from_wall
!        msh%points(i)%x(3) = 2./3.*llz*msh%points(i)%x(3)
!     end do
!
!   end subroutine user_mesh_scale
!
!   ! User defined initial condition
!   subroutine user_ic(u, v, w, p, params)
!     type(field_t), intent(inout) :: u
!     type(field_t), intent(inout) :: v
!     type(field_t), intent(inout) :: w
!     type(field_t), intent(inout) :: p
!     type(json_file), intent(inout) :: params
!     integer :: i
!
!     real(kind=rp) :: uvw(3)
!
!     do i = 1, u%dof%size()
!        uvw = channel_ic(u%dof%x(i,1,1,1),u%dof%y(i,1,1,1),u%dof%z(i,1,1,1))
!        u%x(i,1,1,1) = uvw(1)
!        v%x(i,1,1,1) = uvw(2)
!        w%x(i,1,1,1) = uvw(3)
!     end do
!   end subroutine user_ic
!
!   ! This is called at the end of every time step
!   subroutine usercheck(t, tstep, u, v, w, p, coef, param)
!     real(kind=rp), intent(in) :: t
!     integer, intent(in) :: tstep
!     type(coef_t), intent(inout) :: coef
!     type(field_t), intent(inout) :: u
!     type(field_t), intent(inout) :: v
!     type(field_t), intent(inout) :: w
!     type(field_t), intent(inout) :: p
!     type(json_file), intent(inout) :: param
!
!     ! insert code below
!     type(field_t), pointer :: u_old, u_older, dudy_old
!     integer :: i, j, k, e, n, i_node
!     class(bc_t), pointer :: bc
!     type(wall_model_bc_t), pointer :: wall_bc
!     type(field_t), pointer :: state_old, state_older, action_old, action_older, terminal_old, terminal_older
!
!     u_old => neko_field_registry%get_field("u_old")
!     u_older => neko_field_registry%get_field("u_older")
!
!     state_old => neko_field_registry%get_field("state_old")
!     state_older => neko_field_registry%get_field("state_older")
!
!     action_old => neko_field_registry%get_field("action_old")
!     action_older => neko_field_registry%get_field("action_older")
!
!     terminal_old => neko_field_registry%get_field("terminal_old")
!     terminal_older => neko_field_registry%get_field("terminal_older")
!
!     print *, "***** USERCHECK *****"
!     do e = 1, 2
!       do k = 1, 2
!         do j = 1, 2
!           do i = 1, 2
!             print *, u_older%x(i,j,k,e), u_old%x(i,j,k,e), u%x(i,j,k,e)
!           end do
!         end do
!       end do
!     end do
!      print *, "***** USERCHECK *****"
! !      call dudxyz(dudy, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
!
!     print *, "coef%dof%size() = ", coef%dof%size()
!     print *, "size(u%x) = ", size(u%x)
!     call copy(u_older%x, u_old%x, size(u_old%x))
!     call copy(u_old%x, u%x, size(u%x)) ! coef%dof%size()
! !     call copy(dudy_old%x, u%x, size(u%x)) ! coef%dof%size()
!
!     n = neko_simcomps%case%fluid%bcs_vel%size()
!     do i = 1, n
!       bc => neko_simcomps%case%fluid%bcs_vel%get(i)
!       select type (bc)
!       type is (wall_model_bc_t)
!         wall_bc => bc
!         select type(spald => wall_bc%wall_model)
!         type is (spalding_t)
!           if (allocated(spald%state)) then
!             print *, "spalding state(1,1): ", spald%state(1,1)
!             print *, "shape(spald%state): ", shape(spald%state)
!             call copy(state_older%x, state_old%x, size(state_old%x))
!             call copy(state_old%x, spald%state, size(spald%state))
!             call copy(action_older%x, action_old%x, size(action_old%x))
!             call copy(action_old%x, spald%action, size(spald%action))
!             print *, "Total number of wall nodes = ", spald%n_nodes
!             print *, "Current time = ", t
!             print *, "End time = ", neko_simcomps%case%time%end_time
!             do i_node = 1, spald%n_nodes
!               if (abs(t - neko_simcomps%case%time%end_time) .le. 1.0e-3_rp) then
!                 spald%terminal(1,i_node) = 1.0_rp
!                 if (i_node >= 1 .and. i_node <= 5) then
!                   print *, spald%terminal(1,i_node)
!                 end if
!               else
!                 spald%terminal(1,i_node) = 0.0_rp
!                 if (i_node >= 1 .and. i_node <= 5) then
!                   print *, spald%terminal(1,i_node)
!                 end if
!               end if
!             end do
!             call copy(terminal_older%x, terminal_old%x, size(terminal_old%x))
!             call copy(terminal_old%x, spald%terminal, size(spald%terminal))
!           end if
!         end select
!       end select
!     end do
!
!   end subroutine usercheck
!
!   ! Kind of brute force with rather large initial disturbances
!   function channel_ic(x, y, z) result(uvw)
!     real(kind=rp) :: x, y, z
!     real(kind=rp) :: uvw(3)
!     real(kind=rp) :: ux, uy, uz, eps, Re_tau, yp, Re_b, alpha, beta
!     real(kind=rp) :: C, k, kx, kz, eps1, ran
!
!     real(kind=rp) :: llx, llz
!
!     llx = 4.*pi
!     llz = 4./3.*pi
!
!     Re_tau = 180
!     C      = 5.17
!     k      = 0.41
!     Re_b   = 2800
!
!     yp = (1-y)*Re_tau
!     if (y .lt. 0) yp = (1+y)*Re_tau
!
!     ! Reichardt function
!     ux  = 1/k*log(1.0+k*yp) + (C - (1.0/k)*log(k)) * &
!          (1.0 - exp(-yp/11.0) - yp/11*exp(-yp/3.0))
!     ux  = ux * Re_tau/Re_b
!
!     ! actually, sometimes one may not use the turbulent profile, but
!     ! rather the parabolic lamianr one
!     ! ux = 1.5*(1-y**2)
!
!     ! add perturbations to trigger turbulence
!     ! base flow
!     uvw(1)  = ux
!     uvw(2)  = 0
!     uvw(3)  = 0
!
!     ! first, large scale perturbation
!     eps = 0.05
!     kx  = 3
!     kz  = 4
!     alpha = kx * 2*PI/llx
!     beta  = kz * 2*PI/llz
!     uvw(1)  = uvw(1) + eps*beta  * sin(alpha*x)*cos(beta*z)
!     uvw(2)  = uvw(2) + eps       * sin(alpha*x)*sin(beta*z)
!     uvw(3)  = uvw(3) -eps*alpha * cos(alpha*x)*sin(beta*z)
!
!     ! second, small scale perturbation
!     eps = 0.005
!     kx  = 17
!     kz  = 13
!     alpha = kx * 2*PI/llx
!     beta  = kz * 2*PI/llz
!     uvw(1)  = uvw(1) + eps*beta  * sin(alpha*x)*cos(beta*z)
!     uvw(2)  = uvw(2) + eps       * sin(alpha*x)*sin(beta*z)
!     uvw(3)  = uvw(3) -eps*alpha * cos(alpha*x)*sin(beta*z)
!
!     ! finally, random perturbations only in y
!     eps1 = 0.001
!     ran = sin(-20*x*z+y**3*tan(x*z**2)+100*z*y-20*sin(x*y*z)**5)
!     uvw(2)  = uvw(2) + eps1*ran
!
!   end function channel_ic
!
! end module user
