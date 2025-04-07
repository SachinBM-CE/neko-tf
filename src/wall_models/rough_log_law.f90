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
!> Implements `rough_log_law_t`.
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

    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")

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

      magu = sqrt(ui**2 + vi**2 + wi**2)

      ! Compute the stress
      utau = (magu - this%B) * this%kappa / log(this%h%x(i) / this%z0)

      ! Distribute according to the velocity vector
      this%tau_x(i) = -utau**2 * ui / magu
      this%tau_y(i) = -utau**2 * vi / magu
      this%tau_z(i) = -utau**2 * wi / magu
    end do

  end subroutine rough_log_law_compute


end module rough_log_law

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
! module rough_log_law
!   use field, only: field_t
!   use num_types, only : rp
!   use json_module, only : json_file
!   use dofmap, only : dofmap_t
!   use coefs, only : coef_t
!   use neko_config, only : NEKO_BCKND_DEVICE
!   use wall_model, only : wall_model_t
!   use field_registry, only : neko_field_registry
!   use json_utils, only : json_get_or_default, json_get
!   use utils, only : neko_error
!   ! Reinforcement Learning (sbm)
!   use torchfort
!   use operators, only: grad
!   implicit none
!   private
!
!   !> Wall model based on the log-law for a rough wall.
!   !! The formula defining the law is \f$ u^+ = log(z/z_0)/\kappa + B \f$.
!   !! Here, \f$ z \f$ is the wall-normal distance, as per tradition in
!   !! atmospheric sciences, where this law is often used.
!   type, public, extends(wall_model_t) :: rough_log_law_t
!
!      !> The von Karman coefficient.
!      real(kind=rp) :: kappa = 0.41_rp
!      !> The log-law intercept
!      real(kind=rp) :: B = 0.0_rp
!      !> The roughness height
!      real(kind=rp) :: z0 = 0.0_rp
!      ! New member variable to store coef (sbm)
! !      type(coef_t) :: coef
!    contains
!      !> Constructor from JSON.
!      procedure, pass(this) :: init => rough_log_law_init
!      !> Constructor from components.
!      procedure, pass(this) :: init_from_components => &
!        rough_log_law_init_from_components
!      !> Destructor.
!      procedure, pass(this) :: free => rough_log_law_free
!      !> Compute the wall shear stress.
!      procedure, pass(this) :: compute => rough_log_law_compute
!   end type rough_log_law_t
!
! contains
!   !> Constructor from JSON.
!   !! @param coef SEM coefficients.
!   !! @param msk The boundary mask.
!   !! @param facet The boundary facets.
!   !! @param nu The molecular kinematic viscosity.
!   !! @param h_index The off-wall index of the sampling cell.
!   !! @param json A dictionary with parameters.
!   subroutine rough_log_law_init(this, coef, msk, facet, nu, h_index, json)
!     class(rough_log_law_t), intent(inout) :: this
!     type(coef_t), intent(in) :: coef
!     integer, intent(in) :: msk(:)
!     integer, intent(in) :: facet(:)
!     real(kind=rp), intent(in) :: nu
!     integer, intent(in) :: h_index
!     type(json_file), intent(inout) :: json
!     real(kind=rp) :: kappa, B, z0
!
!     call json_get_or_default(json, "kappa", kappa, 0.41_rp)
!     call json_get(json, "B", B)
!     call json_get(json, "z0", z0)
!
!     call this%init_from_components(coef, msk, facet, nu, h_index, kappa, B, z0)
!   end subroutine rough_log_law_init
!
!   !> Constructor from components.
!   !! @param coef SEM coefficients.
!   !! @param msk The boundary mask.
!   !! @param facet The boundary facets.
!   !! @param nu The molecular kinematic viscosity.
!   !! @param h_index The off-wall index of the sampling cell.
!   !! @param kappa The von Karman coefficient.
!   !! @param B The log-law intercept.
!   !! @param z0 The roughness height.
!   subroutine rough_log_law_init_from_components(this, coef, msk, facet,&
!                                                 nu, h_index, kappa, B, z0)
!     class(rough_log_law_t), intent(inout) :: this
!     type(coef_t), intent(in) :: coef
!     integer, intent(in) :: msk(:)
!     integer, intent(in) :: facet(:)
!     real(kind=rp), intent(in) :: nu
!     integer, intent(in) :: h_index
!     real(kind=rp), intent(in) :: kappa
!     real(kind=rp), intent(in) :: B
!     real(kind=rp), intent(in) :: z0
!
!     if (NEKO_BCKND_DEVICE .eq. 1) then
!        call neko_error("The rough loglaw is only available on the CPU backend.")
!     end if
!
!     call this%init_base(coef, msk, facet, nu, h_index)
!
!     this%kappa = kappa
!     this%B = B
!     this%z0 = z0
!     this%coef = coef  ! Store coef (sbm)
!
!   end subroutine rough_log_law_init_from_components
!
!   !> Destructor for the rough_log_law_t (base) class.
!   subroutine rough_log_law_free(this)
!     class(rough_log_law_t), intent(inout) :: this
!
!     call this%free_base()
!
!   end subroutine rough_log_law_free
!
!   !> Compute the wall shear stress.
!   !> @param t The time value.
!   !> @param tstep The time iteration.
!   subroutine rough_log_law_compute(this, t, tstep)
!     class(rough_log_law_t), intent(inout) :: this
!     real(kind=rp), intent(in) :: t
!     integer, intent(in) :: tstep
!     type(field_t), pointer :: u
!     type(field_t), pointer :: v
!     type(field_t), pointer :: w
!     integer :: i
!     real(kind=rp) :: ui, vi, wi, magu, utau, normu
!     real(kind=rp), dimension(:,:), allocatable :: dus_dx, dus_dy, dus_dz !(sbm)
!     real(kind=rp), dimension(:,:), allocatable :: dvs_dx, dvs_dy, dvs_dz !(sbm)
!     real(kind=rp), dimension(:,:), allocatable :: dws_dx, dws_dy, dws_dz !(sbm)
!     real(kind=rp), dimension(:), allocatable :: dus_dn, dvs_dn, dws_dn !(sbm)
!     integer :: ir, is, it, ie, lid, lx, ly, lxyz, nelv !(sbm)
!     lx = this%coef%Xh%lx !(sbm)
!     ly = this%coef%Xh%ly !(sbm)
!     lxyz = this%coef%Xh%lxyz !(sbm)
!     nelv = this%coef%msh%nelv !(sbm)
!
!     ! Allocate arrays for gradient components (sbm)
!     allocate(dus_dx(lxyz, nelv), dus_dy(lxyz, nelv), dus_dz(lxyz, nelv))
!     allocate(dvs_dx(lxyz, nelv), dvs_dy(lxyz, nelv), dvs_dz(lxyz, nelv))
!     allocate(dws_dx(lxyz, nelv), dws_dy(lxyz, nelv), dws_dz(lxyz, nelv))
!     allocate(dus_dn(this%n_nodes), dvs_dn(this%n_nodes), dws_dn(this%n_nodes))
!
!     u => neko_field_registry%get_field("u")
!     v => neko_field_registry%get_field("v")
!     w => neko_field_registry%get_field("w")
!
!     call grad(dus_dx, dus_dy, dus_dz, u%x, this%coef)
!     call grad(dvs_dx, dvs_dy, dvs_dz, v%x, this%coef)
!     call grad(dws_dx, dws_dy, dws_dz, w%x, this%coef)
!
!     do i = 1, this%n_nodes
!
!       ir = this%ind_r(i)
!       is = this%ind_s(i)
!       it = this%ind_t(i)
!       ie = this%ind_e(i)
!       lid = ir + lx * (is - 1 + ly * (it - 1))
!
!       ! Sample the velocity
!       ui = u%x(ir, is, it, ie)
!       vi = v%x(ir, is, it, ie)
!       wi = w%x(ir, is, it, ie)
!
!       ! Project on tangential direction
!       normu = ui * this%n_x%x(i) + vi * this%n_y%x(i) + wi * this%n_z%x(i)
!
!       ui = ui - normu * this%n_x%x(i)
!       vi = vi - normu * this%n_y%x(i)
!       wi = wi - normu * this%n_z%x(i)
!
!       dus_dn(i) = dus_dx(lid, ie)*this%n_x%x(i) + dus_dy(lid, ie)*this%n_y%x(i) + dus_dz(lid, ie)*this%n_z%x(i)
!       dvs_dn(i) = dvs_dx(lid, ie)*this%n_x%x(i) + dvs_dy(lid, ie)*this%n_y%x(i) + dvs_dz(lid, ie)*this%n_z%x(i)
!       dws_dn(i) = dws_dx(lid, ie)*this%n_x%x(i) + dws_dy(lid, ie)*this%n_y%x(i) + dws_dz(lid, ie)*this%n_z%x(i)
!
!       print *, i, dus_dn(i), dvs_dn(i), dws_dn(i)
! !       write(*, '(A,I4,A,3(ES13.5))') 'Wall-normal grads at node ', i, ': ', dus_dn(i), dvs_dn(i), dws_dn(i)
!
!       magu = sqrt(ui**2 + vi**2 + wi**2)
!
!       ! Compute the stress
!       utau = (magu - this%B) * this%kappa / log(this%h%x(i) / this%z0)
!
!       ! Distribute according to the velocity vector
!       this%tau_x(i) = -utau**2 * ui / magu
!       this%tau_y(i) = -utau**2 * vi / magu
!       this%tau_z(i) = -utau**2 * wi / magu
!
!     end do
!
!     deallocate(dus_dx, dus_dy, dus_dz, dus_dn, dvs_dn, dws_dn)
!
!   end subroutine rough_log_law_compute
!
! end module rough_log_law
