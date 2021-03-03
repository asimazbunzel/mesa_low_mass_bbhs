
      module cc_support

      use cc_def
      use cc_mechanism
      use cc_data
      use star_def
      use star_lib
      use binary_def
      use binary_lib

      implicit none

      contains

      logical function is_binary(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         is_binary = .false.
         if (s% binary_id == 1) is_binary = .true.

      end function is_binary


      subroutine do_compact_object_formation(cc_id, cc_binary_id, ierr)
         integer, intent(in) :: cc_id
         integer, intent(in) :: cc_binary_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         include 'formats.inc'

         ierr = 0

         ! load star pointer
         call star_ptr(cc_id, s, ierr)
         if (failed('failed star_ptr', ierr)) return

         ! evaluate core-collapse
         call do_explode_star(s, ierr)
         if (failed('failed do_explode_star', ierr)) return

         ! update star parameters
         call do_star_update_after_cc(s, ierr)
         if (failed('failed do_star_update_after_cc', ierr)) return

         ! update binary parameters if needed
         if (cc_binary_id == 1) then
            if (continue_binary_evolution) then 
               call do_binary_update_after_cc(cc_binary_id, cc_id, ierr)
               if (failed('failed do_star_update_after_cc', ierr)) return
            end if
         end if
         
      end subroutine do_compact_object_formation


      logical function failed(str,ierr)
         character (len=*), intent(in) :: str
         integer, intent(in) :: ierr

         include 'formats.inc'

         failed = (ierr /= 0)

         if (failed) write(*,11) trim(str) // ' ierr', ierr

      end function failed

      end module cc_support
