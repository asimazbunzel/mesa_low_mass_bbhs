
      module ce_utils

      use ce_def
      use const_def
      use star_def
      use star_lib
      use binary_def
      use binary_lib
      
      implicit none

      contains


      ! count ce events
      subroutine num_events(ce_id, ierr)
         integer, intent(in) :: ce_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b

         ierr = 0
         call binary_ptr(ce_id, b, ierr)

         if (b% d_i == 1) then
            ce_num_star_1 = ce_num_star_1 + 1
         else
            ce_num_star_2 = ce_num_star_2 + 1
         end if

         if (ce_type == ce_two_stars) then
            ce_num_two_stars = ce_num_two_stars + 1
         else
            ce_num_xrb = ce_num_xrb + 1
         end if

      end subroutine num_events


      ! evaluate binding energy of star
      subroutine eval_binding_energy(s, inner_cell, outer_cell, gravitational_energy, internal_energy, ierr)
         type(star_info), pointer :: s
         integer, intent(in) :: inner_cell, outer_cell
         integer :: k
         real(dp), intent(out) :: gravitational_energy, internal_energy
         integer, intent(out) :: ierr

         ierr = 0

         gravitational_energy = 0d0
         internal_energy = 0d0
         do k = outer_cell, inner_cell
            internal_energy = internal_energy + s% dm(k) * s% energy(k)
            gravitational_energy = gravitational_energy - standard_cgrav * s% m(k) * s% dm_bar(k) / s% r(k)
         end do

      end subroutine eval_binding_energy

      end module ce_utils
