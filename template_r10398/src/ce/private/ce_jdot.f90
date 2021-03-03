
      module ce_jdot

      use ce_def
      use const_def
      use star_def
      use star_lib
      use binary_def
      use binary_lib
      
      implicit none

      contains


      ! set extra_jdot
      subroutine eval_ce_jdot(b, ierr)
         use ce_utils, only: eval_binding_energy
         type(binary_info), pointer :: b
         integer, intent(out) :: ierr
         real(dp) :: gravitational_energy, internal_energy
         real(dp) :: removed_binding_energy
         real(dp) :: alpha, separation
         real(dp) :: jorb

         ierr = 0

         b% extra_jdot = 0

         ! just be extra-sure that binary is circular
         if (b% eccentricity > 0d0) b% eccentricity = 0d0

         ! binding energy of the whole star
         call eval_binding_energy(b% s_donor, &
                         ! limits of integration
                         b% s_donor% nz, 0, &
                         gravitational_energy, internal_energy, &
                         ierr)
         if (ierr /= 0) return
         donor_bind_energy = gravitational_energy + internal_energy

         if (b% model_number <= ce_initial_model_number) return

         ! removed energy is just the difference between the actual binding energy and the previous timestep
         removed_binding_energy = donor_bind_energy_prev_step - donor_bind_energy

         if (add_accretion_on_ce) then
            write(*,'(a)') 'This feature is not ready so nothing will be calculated'
         end if

         ! check for type of ce
         if (ce_type == ce_two_stars) then
            alpha = alpha_two_stars
         else
            alpha = alpha_xrb
         end if

         separation = -(alpha * standard_cgrav * b% m(b% d_i) * b% m(b% a_i)) / &
            (2d0 * (removed_binding_energy + (alpha * orbital_energy_prev_step)))

         separation = min(separation, b% separation_old)

         jorb = b% m(b% d_i) * b% m(b% a_i) * sqrt(standard_cgrav * &
            separation / (b% m(b% d_i) + b% m(b% a_i)))

         ! In `binary_evolve`: b% jdot*b% time_step*secyer
         ! so b% extra_jdot must be in cgs units
         b% extra_jdot = jorb - b% angular_momentum_j_old
         b% extra_jdot = b% extra_jdot / (b% time_step * secyer)

      end subroutine eval_ce_jdot

      end module ce_jdot

