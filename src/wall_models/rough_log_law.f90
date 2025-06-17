! ! Copyright (c) 2024, The Neko Authors
! ! All rights reserved.
! !
! ! Redistribution and use in source and binary forms, with or without
! ! modification, are permitted provided that the following conditions
! ! are met:
! !
! !   * Redistributions of source code must retain the above copyright
! !     notice, this list of conditions and the following disclaimer.
! !
! !   * Redistributions in binary form must reproduce the above
! !     copyright notice, this list of conditions and the following
! !     disclaimer in the documentation and/or other materials provided
! !     with the distribution.
! !
! !   * Neither the name of the authors nor the names of its
! !     contributors may be used to endorse or promote products derived
! !     from this software without specific prior written permission.
! !
! ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! ! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! ! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! ! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! ! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! ! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! ! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! ! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! ! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! ! POSSIBILITY OF SUCH DAMAGE.
! !
! !
! !> Implements `rough_log_law_t`.
module rough_log_law
  use field, only: field_t
  use num_types, only : rp
  use json_module, only : json_file
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use wall_model, only : wall_model_t
  use field_registry, only : neko_field_registry
  use json_utils, only : json_get_or_default, json_get
  use utils, only : neko_error

  ! TorchFort
  use torchfort
  use operators, only : grad, dudxyz
  use comm, only : pe_rank, pe_size, NEKO_COMM
  use utils, only : linear_index
  use tf_module
  use iso_c_binding
  use math
  use num_types, only : sp

  implicit none
  private

  !> Wall model based on the log-law for a rough wall.
  !! The formula defining the law is \f$ u^+ = log(z/z_0)/\kappa + B \f$.
  !! Here, \f$ z \f$ is the wall-normal distance, as per tradition in
  !! atmospheric sciences, where this law is often used.
  type, public, extends(wall_model_t) :: rough_log_law_t

     !> The von Karman coefficient.
     real(kind=rp) :: kappa = 0.41_rp
     !> The log-law intercept
     real(kind=rp) :: B = 0.0_rp
     !> The roughness height
     real(kind=rp) :: z0 = 0.0_rp
     ! Array Initialization for TorchFort
     real(kind=rp), dimension(:,:,:,:), allocatable :: dudy
     real(kind=rp), dimension(:,:), allocatable :: state, action, reward, terminal, total_reward
     real(kind=rp), dimension(:,:), allocatable :: state_2d, action_2d, terminal_2d
     real(kind=rp), dimension(:,:), allocatable :: l_star, u_plus, g_plus, h_plus

   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => rough_log_law_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
       rough_log_law_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => rough_log_law_free
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => rough_log_law_compute
  end type rough_log_law_t

contains
  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param nu The molecular kinematic viscosity.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param json A dictionary with parameters.
  subroutine rough_log_law_init(this, coef, msk, facet, nu, h_index, json)

    class(rough_log_law_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    real(kind=rp), intent(in) :: nu
    integer, intent(in) :: h_index
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, B, z0

    call json_get_or_default(json, "kappa", kappa, 0.41_rp)
    call json_get(json, "B", B)
    call json_get(json, "z0", z0)

    call this%init_from_components(coef, msk, facet, nu, h_index, kappa, B, z0)

    ! Allocation of arrays
    allocate(this%dudy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(this%state(this%n_nodes,3), this%action(this%n_nodes,1), this%terminal(this%n_nodes,1))
    allocate(this%reward(this%n_nodes,1), this%total_reward(this%n_nodes,1))
    allocate(this%state_2d(this%n_nodes,3), this%action_2d(this%n_nodes,1), this%terminal_2d(this%n_nodes,1))
    allocate(this%l_star(this%n_nodes,1), this%u_plus(this%n_nodes,1), this%g_plus(this%n_nodes,1), this%h_plus(this%n_nodes,1))

    call neko_field_registry%add_field(coef%dof, "state_old", ignore_existing = .true.)
    call neko_field_registry%add_field(coef%dof, "state_older", ignore_existing = .true.)
    call neko_field_registry%add_field(coef%dof, "action_old", ignore_existing = .true.)
    call neko_field_registry%add_field(coef%dof, "action_older", ignore_existing = .true.)
    call neko_field_registry%add_field(coef%dof, "terminal_old", ignore_existing = .true.)
    call neko_field_registry%add_field(coef%dof, "terminal_older", ignore_existing = .true.)

    print *, "rough_log_law_init called!"

  end subroutine rough_log_law_init

  !> Constructor from components.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param nu The molecular kinematic viscosity.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param kappa The von Karman coefficient.
  !! @param B The log-law intercept.
  !! @param z0 The roughness height.
  subroutine rough_log_law_init_from_components(this, coef, msk, facet,&
                                                nu, h_index, kappa, B, z0)
    class(rough_log_law_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    real(kind=rp), intent(in) :: nu
    integer, intent(in) :: h_index
    real(kind=rp), intent(in) :: kappa
    real(kind=rp), intent(in) :: B
    real(kind=rp), intent(in) :: z0

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("The rough loglaw is only available on the CPU backend.")
    end if

    call this%init_base(coef, msk, facet, nu, h_index)

    this%kappa = kappa
    this%B = B
    this%z0 = z0

  end subroutine rough_log_law_init_from_components

  !> Destructor for the rough_log_law_t (base) class.
  subroutine rough_log_law_free(this)

    class(rough_log_law_t), intent(inout) :: this
    call this%free_base()

    if (allocated(this%dudy)) deallocate(this%dudy)
    if (allocated(this%state)) deallocate(this%state)
    if (allocated(this%action)) deallocate(this%action)
    if (allocated(this%reward)) deallocate(this%reward)
    if (allocated(this%total_reward)) deallocate(this%total_reward)
    if (allocated(this%terminal)) deallocate(this%terminal)

    if (allocated(this%state_2d)) deallocate(this%state_2d)
    if (allocated(this%action_2d)) deallocate(this%action_2d)
    if (allocated(this%terminal_2d)) deallocate(this%terminal_2d)

    if (allocated(this%l_star)) deallocate(this%l_star)
    if (allocated(this%u_plus)) deallocate(this%u_plus)
    if (allocated(this%g_plus)) deallocate(this%g_plus)
    if (allocated(this%h_plus)) deallocate(this%h_plus)

    print *, "spalding_free called!"

  end subroutine rough_log_law_free

  !> Compute the wall shear stress.
  !> @param t The time value.
  !> @param tstep The time iteration.
  subroutine rough_log_law_compute(this, t, tstep)

    class(rough_log_law_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu

    ! TorchFort Variable Declarations
    integer :: res, zeros
    logical :: is_ready = .false., is_end = .false.
    real(kind=sp) :: p_loss_val, q_loss_val
    type(field_t), pointer :: u_prev
    real(kind=rp) :: tau_old, tau_new, reward_mean
    real(kind=rp) :: tau_true = 1.0, error_now, error_prev, base_reward, rel_error, bonus
    type(field_t), pointer :: state_old, state_older, action_old, action_older, terminal_old

    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")

    ! Gradient Tensor
    call dudxyz(this%dudy, u%x, this%coef%drdy, this%coef%dsdy, this%coef%dtdy, this%coef)
    print *, "size(this%dudy) = ", size(this%dudy)

    do i = 1, this%n_nodes

      ! Sample the velocity
      ui = u%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))
      vi = v%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))
      wi = w%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))

      ! Project on tangential direction
      normu = ui * this%n_x%x(i) + vi * this%n_y%x(i) + wi * this%n_z%x(i)
      ui = ui - normu * this%n_x%x(i)
      vi = vi - normu * this%n_y%x(i)
      wi = wi - normu * this%n_z%x(i)

      ! Magnitude of Velocity
      magu = sqrt(ui**2 + vi**2 + wi**2)

      ! Compute the stress
      utau = (magu - this%B) * this%kappa / log(this%h%x(i) / this%z0)

      ! Normalization of inputs
      this%l_star(i,1) = this%nu / utau
      this%state(i,1) = magu / utau
      this%state(i,2) = this%dudy(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i)) / (utau / this%l_star(i,1))
      this%state(i,3) = this%h%x(i) / this%l_star(i,1)

      if (i >= 1 .and. i <= 5) then
        if (i == 1) then
          write(*, '(A20, F10.4)') '1st do loop: t =', t
          write(*, '(A6, 2X, A20, 2X, A20, 2X, A20)') &
            'i', 'state(i,1)', 'state(i,2)', 'state(i,3)'
        end if
        write(*, '(I6, 2X, ES20.10, 2X, ES20.10, 2X, ES20.10)') &
        i, this%state(i,1), this%state(i,2), this%state(i,3)
      end if

      ! Distribute according to the velocity vector
      this%tau_x(i) = -utau**2 * ui / magu
      this%tau_y(i) = -utau**2 * vi / magu
      this%tau_z(i) = -utau**2 * wi / magu

    end do

    ! Exploitation
    res = torchfort_rl_off_policy_predict(tf_key, reshape(this%state, [3,this%n_nodes]), &
    reshape(this%action, [1,this%n_nodes]))
    if (res /= TORCHFORT_RESULT_SUCCESS) stop
    print *, "result of predict_float_2d_2d: ", res

    ! Exploration
    res = torchfort_rl_off_policy_predict_explore(tf_key, reshape(this%state, [3,this%n_nodes]), &
    reshape(this%action, [1,this%n_nodes]))
    if (res /= TORCHFORT_RESULT_SUCCESS) stop
    print *, "result of predict_explore: ", res

    ! Ready for training
    res = torchfort_rl_off_policy_is_ready(tf_key, is_ready)
    if (res /= TORCHFORT_RESULT_SUCCESS) stop
    print *, "result of policy_is_ready: ", res

    ! Get fields from Neko Field Registry
    state_old => neko_field_registry%get_field("state_old")
    state_older => neko_field_registry%get_field("state_older")
    action_old => neko_field_registry%get_field("action_old")
    action_older => neko_field_registry%get_field("action_older")
    terminal_old => neko_field_registry%get_field("terminal_old")

  end subroutine rough_log_law_compute


end module rough_log_law
