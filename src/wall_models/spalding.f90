! Copyright (c) 2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!
!> Implements `spalding_t`.
module spalding
  use field, only: field_t
  use num_types, only : rp
  use json_module, only : json_file
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use wall_model, only : wall_model_t
  use field_registry, only : neko_field_registry
  use json_utils, only : json_get_or_default
  use logger, only : neko_log, NEKO_LOG_DEBUG
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

  !> Wall model based on Spalding's law of the wall.
  !! Reference: http://dx.doi.org/10.1115/1.3641728
  type, public, extends(wall_model_t) :: spalding_t

     !> The von Karman coefficient.
     real(kind=rp) :: kappa = 0.41_rp
     !> The log-law intercept.
     real(kind=rp) :: B = 5.2_rp

     ! TorchFort
     real(kind=rp), dimension(:,:,:,:), allocatable :: dudy
     real(kind=rp), dimension(:,:), allocatable :: state, action, reward, terminal, total_reward, so, sor, ao, aor
     real(kind=rp), dimension(:,:), allocatable :: state_2d, action_2d, terminal_2d
     real(kind=rp), dimension(:), allocatable :: l_star, u_plus, g_plus, h_plus
     real(kind=rp), dimension(:), allocatable :: ui_l, vi_l, wi_l, normu_l, magu_l, vg_l, utau_l, tau_old_l, tau_new_l
     real(kind=rp), dimension(:), allocatable :: error_new_l, error_old_l, base_reward_l, rel_error_l, bonus_l

   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => spalding_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
       spalding_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => spalding_free
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => spalding_compute
     !> Solve for the friction velocity
     procedure, private, pass(this) :: solve

  end type spalding_t

contains
  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param nu The molecular kinematic viscosity.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param json A dictionary with parameters.
  subroutine spalding_init(this, coef, msk, facet, nu, h_index, json)

    class(spalding_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    real(kind=rp), intent(in) :: nu
    integer, intent(in) :: h_index
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, B

    ! TorchFort
    logical :: found

    call json_get_or_default(json, "kappa", kappa, 0.41_rp)
    call json_get_or_default(json, "B", B, 5.2_rp)

    call this%init_from_components(coef, msk, facet, nu, h_index, kappa, B)

    ! Allocation of arrays
!     allocate(this%dudy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
!     allocate(this%state(this%n_nodes,3), this%action(this%n_nodes,1), this%terminal(this%n_nodes,1))
!     allocate(this%reward(this%n_nodes,1), this%total_reward(this%n_nodes,1))
!     allocate(this%state_2d(this%n_nodes,3), this%action_2d(this%n_nodes,1), this%terminal_2d(this%n_nodes,1))
!     allocate(this%l_star(this%n_nodes,1), this%u_plus(this%n_nodes,1), this%g_plus(this%n_nodes,1), this%h_plus(this%n_nodes,1))
!     allocate(this%ui_l(this%n_nodes), this%vi_l(this%n_nodes), this%wi_l(this%n_nodes), this%normu_l(this%n_nodes), this%magu_l(this%n_nodes))
!     allocate(this%vg_l(this%n_nodes), this%utau_l(this%n_nodes), this%tau_old_l(this%n_nodes), this%tau_new_l(this%n_nodes))

    allocate(this%dudy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv), &
    this%state(this%n_nodes,3), this%action(this%n_nodes,1), this%terminal(this%n_nodes,1), &
    this%reward(this%n_nodes,1), this%total_reward(this%n_nodes,1), &
    this%state_2d(this%n_nodes,3), this%action_2d(this%n_nodes,1), this%terminal_2d(this%n_nodes,1), &
    this%l_star(this%n_nodes), this%u_plus(this%n_nodes), this%g_plus(this%n_nodes), this%h_plus(this%n_nodes), &
    this%ui_l(this%n_nodes), this%vi_l(this%n_nodes), this%wi_l(this%n_nodes), this%normu_l(this%n_nodes), &
    this%magu_l(this%n_nodes), this%vg_l(this%n_nodes), this%utau_l(this%n_nodes), &
    this%tau_old_l(this%n_nodes), this%tau_new_l(this%n_nodes), this%error_new_l(this%n_nodes), this%error_old_l(this%n_nodes), &
    this%base_reward_l(this%n_nodes), this%rel_error_l(this%n_nodes), this%bonus_l(this%n_nodes), &
    this%so(this%n_nodes,3), this%sor(this%n_nodes,3), this%ao(this%n_nodes,1), this%aor(this%n_nodes,1))

    ! Adding new fields to Neko Field Registry
    call neko_field_registry%add_field(coef%dof, "state_old", ignore_existing = .true.)
    call neko_field_registry%add_field(coef%dof, "state_older", ignore_existing = .true.)
    call neko_field_registry%add_field(coef%dof, "action_old", ignore_existing = .true.)
    call neko_field_registry%add_field(coef%dof, "action_older", ignore_existing = .true.)
    call neko_field_registry%add_field(coef%dof, "terminal_old", ignore_existing = .true.)
    call neko_field_registry%add_field(coef%dof, "terminal_older", ignore_existing = .true.)
    found = neko_field_registry%field_exists("state_old")
    print *, "state_old field_exists: ", found

    print *, "spalding_init called!"

  end subroutine spalding_init

  !> Constructor from components.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param nu The molecular kinematic viscosity.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param kappa The von Karman coefficient.
  !! @param B The log-law intercept.
  subroutine spalding_init_from_components(this, coef, msk, facet, nu, h_index,&
                                           kappa, B)
    class(spalding_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    real(kind=rp), intent(in) :: nu
    real(kind=rp), intent(in) :: kappa
    real(kind=rp), intent(in) :: B

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("Spalding's law is only available on the CPU backend.")
    end if

    call this%init_base(coef, msk, facet, nu, h_index)

    this%kappa = kappa
    this%B = B

  end subroutine spalding_init_from_components


  !> Destructor for the spalding_t (base) class.
  subroutine spalding_free(this)
    class(spalding_t), intent(inout) :: this

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

    if (allocated(this%ui_l)) deallocate(this%ui_l)
    if (allocated(this%vi_l)) deallocate(this%vi_l)
    if (allocated(this%wi_l)) deallocate(this%wi_l)
    if (allocated(this%normu_l)) deallocate(this%normu_l)
    if (allocated(this%magu_l)) deallocate(this%magu_l)
    if (allocated(this%vg_l)) deallocate(this%vg_l)
    if (allocated(this%utau_l)) deallocate(this%utau_l)
    if (allocated(this%tau_old_l)) deallocate(this%tau_old_l)
    if (allocated(this%tau_new_l)) deallocate(this%tau_new_l)

    if (allocated(this%error_new_l)) deallocate(this%error_new_l)
    if (allocated(this%error_old_l)) deallocate(this%error_old_l)
    if (allocated(this%base_reward_l)) deallocate(this%base_reward_l)
    if (allocated(this%rel_error_l)) deallocate(this%rel_error_l)
    if (allocated(this%bonus_l)) deallocate(this%bonus_l)
    if (allocated(this%so)) deallocate(this%so)
    if (allocated(this%sor)) deallocate(this%sor)

    print *, "spalding_free called!"

  end subroutine spalding_free

  !> Compute the wall shear stress.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine spalding_compute(this, t, tstep)

    class(spalding_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu, guess

    ! TorchFort Variable Declarations
    integer :: res, zeros
    logical :: is_ready = .false., is_end = .false.
    real(kind=sp) :: p_loss_val, q_loss_val
    type(field_t), pointer :: u_prev
    real(kind=rp) :: tau_old, tau_new, reward_mean
    real(kind=rp) :: tau_true = 0.004132653061224489_rp, error_new, error_old, base_reward, rel_error, bonus
    type(field_t), pointer :: state_old, state_older, action_old, action_older, terminal_old

    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")

    print *, "size(u%x) = ", size(u%x)
    print *, "size(this%nx%x) = ", size(this%n_x%x)

    u_prev => neko_field_registry%get_field("u_older")

    ! Gradient Tensor
    call dudxyz(this%dudy, u%x, this%coef%drdy, this%coef%dsdy, this%coef%dtdy, this%coef)
    print *, "size(this%dudy) = ", size(this%dudy)

    print *, "** ** ** ** **"
    print *, "At t = ", t, "Rank = ", pe_rank, "n_nodes = ", this%n_nodes
    print *, "** ** ** ** **"

    do i=1, this%n_nodes

      ! Sample the velocity
      this%ui_l(i) = u%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))
      this%vi_l(i) = v%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))
      this%wi_l(i) = w%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))

      ! Project on tangential direction
      this%normu_l(i) = this%ui_l(i) * this%n_x%x(i) + this%vi_l(i) * this%n_y%x(i) + this%wi_l(i) * this%n_z%x(i)
      this%ui_l(i) = this%ui_l(i) - this%normu_l(i) * this%n_x%x(i)
      this%vi_l(i) = this%vi_l(i) - this%normu_l(i) * this%n_y%x(i)
      this%wi_l(i) = this%wi_l(i) - this%normu_l(i) * this%n_z%x(i)

      this%magu_l(i) = sqrt(this%ui_l(i)**2 + this%vi_l(i)**2 + this%wi_l(i)**2)

      if (tstep .eq. 1) then

        this%vg_l(i) = sqrt(this%magu_l(i) * this%nu / this%h%x(i))
        this%utau_l(i) =  this%solve(this%magu_l(i), this%h%x(i), this%vg_l(i))

        this%state(i,1) = this%magu_l(i)
        this%state(i,2) = this%dudy(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))
        this%state(i,3) = this%h%x(i)

        if (i >= 1 .and. i <= 5) then
            if (i == 1) then
                write(*, '(A20, F10.4)') '1st do loop: t =', t
                write(*, '(A6, 2X, A20, 2X, A20, 2X, A20)') &
                    'i', 'state(i,1)', 'state(i,2)', 'state(i,3)'
            end if
            write(*, '(I6, 2X, ES20.10, 2X, ES20.10, 2X, ES20.10)') &
            i, this%state(i,1), this%state(i,2), this%state(i,3)
        end if

      else

        this%tau_old_l(i) = sqrt(this%tau_x(i)**2 + this%tau_y(i)**2 + this%tau_z(i)**2)
        this%vg_l(i) = sqrt(this%tau_old_l(i))
        this%utau_l(i) =  this%solve(this%magu_l(i), this%h%x(i), this%vg_l(i))

        ! Normalization of inputs
        this%l_star(i) = this%nu / this%utau_l(i)
        this%u_plus(i) = this%magu_l(i) / this%utau_l(i)
        this%g_plus(i) = this%dudy(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i)) &
        / (this%utau_l(i) / this%l_star(i))
        this%h_plus(i) = this%h%x(i) / this%l_star(i)

        ! Changing to normalized states
        this%state(i,1) = this%u_plus(i)
        this%state(i,2) = this%g_plus(i)
        this%state(i,3) = this%h_plus(i)

        if (i >= 1 .and. i <= 5) then
            if (i == 1) then
            write(*, '(A20, F10.4)') '1st do loop: t =', t
            write(*, '(A6, 2X, A20, 2X, A20, 2X, A20)') &
                'i', 'state(i,1)', 'state(i,2)', 'state(i,3)'
            end if
            write(*, '(I6, 2X, ES20.10, 2X, ES20.10, 2X, ES20.10)') &
            i, this%state(i,1), this%state(i,2), this%state(i,3)
        end if

      end if

    end do ! End of 1st Do Loop

!     res = torchfort_rl_off_policy_predict_float_2d_2d(tf_key, &
!     reshape(this%state, [3, this%n_nodes]), reshape(this%action, [1, this%n_nodes]))
    res = torchfort_rl_off_policy_predict_float_2d_2d(tf_key, transpose(this%state), transpose(this%action))
    if (res /= TORCHFORT_RESULT_SUCCESS) stop
    print *, "result of predict_float_2d_2d: ", res

!     res = torchfort_rl_off_policy_predict_explore_float_2d_2d(tf_key, reshape(this%state, [3,this%n_nodes]), &
!     reshape(this%action, [1,this%n_nodes]))
!     if (res /= TORCHFORT_RESULT_SUCCESS) stop
!     print *, "result of predict_explore: ", res

    res = torchfort_rl_off_policy_is_ready(tf_key, is_ready)
    if (res /= TORCHFORT_RESULT_SUCCESS) stop
    print *, "result of policy_is_ready: ", res

    ! Get fields from Neko Field Registry
    state_old => neko_field_registry%get_field("state_old")
    state_older => neko_field_registry%get_field("state_older")
    action_old => neko_field_registry%get_field("action_old")
    action_older => neko_field_registry%get_field("action_older")
    terminal_old => neko_field_registry%get_field("terminal_old")

    print *, "====== Start of 2nd do loop ======"
    do i = 1, this%n_nodes

      if (tstep .eq. 1) then

        ! Distribute according to the velocity vector
        this%tau_x(i) = -this%utau_l(i)**2 * this%ui_l(i) / this%magu_l(i)
        this%tau_y(i) = -this%utau_l(i)**2 * this%vi_l(i) / this%magu_l(i)
        this%tau_z(i) = -this%utau_l(i)**2 * this%wi_l(i) / this%magu_l(i)

      else

        ! Distribute according to the velocity vector
        this%tau_x(i) = -this%utau_l(i)**2 * this%ui_l(i) / this%magu_l(i)
        this%tau_y(i) = -this%utau_l(i)**2 * this%vi_l(i) / this%magu_l(i)
        this%tau_z(i) = -this%utau_l(i)**2 * this%wi_l(i) / this%magu_l(i)

        this%tau_new_l(i) = sqrt(this%tau_x(i)**2 + this%tau_y(i)**2 + this%tau_z(i)**2)

        ! Base Reward
        this%error_new_l(i) = abs(tau_true - this%tau_new_l(i))
        this%error_old_l(i) = abs(tau_true - this%tau_old_l(i))
        this%base_reward_l(i)= (this%error_new_l(i) - this%error_old_l(i)) / tau_true

        ! Bonus Reward: +1 if rel_error is within 1% of true value
        this%rel_error_l(i) = this%error_new_l(i) / tau_true
        if (this%rel_error_l(i) < 0.01) then
          this%bonus_l(i) = 1.0
        else
          this%bonus_l(i) = 0.0
        end if

        ! Total Reward
        this%reward(i,1) = this%base_reward_l(i) + this%bonus_l(i)
        this%total_reward(i,1) = this%total_reward(i,1) + this%reward(i,1)

!        res = torchfort_rl_off_policy_wandb_log_double_int32step(tf_key, "total_reward", tstep, this%total_reward(i,1))

      end if

      if (i >= 1 .and. i <= 5) then
        if (i == 1) then
            write(*, '(A6, 2X, A10, 2X, A10, 2X, A15, 2X, A15, 2X, A15, 2X, A15)') &
                'i', 'tau_old', 'tau_new', 'sor(i,1)', 'action(i,1)', 'state(i,1)', 'total_reward(i,1)'
        end if
        write(*, '(I6, 2X, ES10.3, 2X, ES10.3, 2X, ES15.6, 2X, ES15.6, 2X, ES15.6, 2X, ES15.6)') &
        i, this%tau_old_l(i), this%tau_new_l(i), this%sor(i,1), this%action(i,1), this%state(i,1), &
        this%total_reward(i,1)
      end if

    end do
    print *, "====== End of 2nd do loop ======"

!     print *, "****** State ******"
!
! !     print *, "state_older%x(1,1) = ", state_older%x(1,1,1,1)
! ! !     zeros = count(state_older%x == 0.0)
! !     print *, "Zeros in state_older%x = ", zeros
!
!     this%state_2d = reshape(state_older%x, [3,this%n_nodes])
!     zeros = count(this%state_2d == 0.0)
! !     print *, "Zeros in this%state_2d = ", zeros
! !     print *, "this%state_2d(1,1) = ", this%state_2d(1,1)
!
!     print *, "****** Action ******"
!
! !     print *, "action_older%x(1,1) = ", action_older%x(1,1,1,1)
!
!     this%action_2d = reshape(action_older%x, [1,this%n_nodes])
!     zeros = count(this%action_2d == 0.0)
! !     print *, "Zeros in this%action_2d = ", zeros
! !     print *, "this%action_2d(1,1) = ", this%action_2d(1,1)
!
!     print *, "****** Terminal ******"
!
! !     print *, "Current time = ", t
! ! !     this%terminal_2d = reshape(terminal_old%x, [1,this%n_nodes])
! ! !     zeros = count(this%terminal_2d == 0.0)
! !     print *, "Zeros in this%terminal_2d = ", zeros
! !     print *, "this%terminal_2d(1,1) = ", this%terminal_2d(1,1)
!
! !     print *, "shape(state_older%x) = ", shape(state_older%x)
!
!     print *, "********************"
!
!     reward_mean = sum(this%total_reward)/this%n_nodes
!     print *, "--- Mean Reward over wall nodes ---"
!     print *, reward_mean
!     print *, "--- --- --- --- --- --- --- --- ---"
!
!     this%state = reshape(this%state, [3,this%n_nodes])
!     this%total_reward = reshape(this%total_reward, [1,this%n_nodes])
!     this%terminal_2d = reshape(this%terminal_2d, [1,this%n_nodes])
!     print *, "shape(this%state_2d) = ", shape(this%state_2d)
!     print *, "shape(this%action_2d) = ", shape(this%action_2d)
!     print *, "shape(this%state) = ", shape(this%state)
!     print *, "shape(this%total_reward) = ", shape(this%total_reward)
!     print *, "shape(this%terminal_2d) = ", shape(this%terminal_2d)

! !     res = torchfort_rl_off_policy_update_replay_buffer_multi_float_2d_2d(tf_key, &
! !     this%state_2d, this%action_2d, this%state, this%total_reward, this%terminal_2d)
! !     if (res /= TORCHFORT_RESULT_SUCCESS) stop
! !     print *, "result of update_replay_buffer_multi: ", res

    res = torchfort_rl_off_policy_update_replay_buffer_multi_float_2d_2d(tf_key, &
    transpose(this%sor), transpose(this%action), transpose(this%state), &
    transpose(this%total_reward), transpose(this%terminal))
    if (res /= TORCHFORT_RESULT_SUCCESS) stop
    print *, "result of update_replay_buffer_multi: ", res

! !     res = torchfort_rl_off_policy_update_replay_buffer_float_2d_2d(tf_key, &
! !     this%state_2d, this%action_2d, this%state, reward_mean, is_end)
! !     if (res /= TORCHFORT_RESULT_SUCCESS) stop
! !     print *, "result of update_replay_buffer: ", res

    res = torchfort_rl_off_policy_train_step_float(tf_key, p_loss_val, q_loss_val)
    print *, "p_loss_val = ", p_loss_val
    print *, "q_loss_val = ", q_loss_val
    if (res /= TORCHFORT_RESULT_SUCCESS) stop
    print *, "result of train_step_float: ", res

    res = torchfort_rl_off_policy_wandb_log_float_int32step(tf_key, "actor_loss", tstep, p_loss_val)
    res = torchfort_rl_off_policy_wandb_log_float_int32step(tf_key, "critic_loss", tstep, q_loss_val)

  end subroutine spalding_compute

  !> Newton solver for the algebraic equation defined by the law.
  !! @param u The velocity value.
  !! @param y The wall-normal distance.
  !! @param guess Initial guess.
  function solve(this, u,  y, guess) result(utau)
    class(spalding_t), intent(inout) :: this
    real(kind=rp), intent(in) :: u
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: guess
    real(kind=rp) :: yp, up, kappa, B, utau
    real(kind=rp) :: error, f, df, old
    integer :: niter, k, maxiter

    utau = guess
    kappa = this%kappa
    B = this%B

    maxiter = 100

    do k=1, maxiter
      up = u / utau
      yp = y * utau / this%nu
      niter = k
      old = utau

      ! Evaluate function and its derivative
      f = (up + exp(-kappa*B)* &
          (exp(kappa*up) - 1.0_rp - kappa*up - 0.5_rp*(kappa*up)**2 - &
           1.0_rp/6*(kappa*up)**3) - yp)

      df = (-y / this%nu - u/utau**2 - kappa*up/utau*exp(-kappa*B) * &
           (exp(kappa*up) - 1 - kappa*up - 0.5*(kappa*up)**2))

      ! Update solution
      utau = utau - f / df

      error = abs((old - utau)/old)

      if (error < 1e-3) then
        exit
      endif

    enddo

    if ((niter .eq. maxiter) .and. (neko_log%level_ .eq. NEKO_LOG_DEBUG)) then
       write(*,*) "Newton not converged", error, f, utau, old, guess
    end if
  end function solve


end module spalding
