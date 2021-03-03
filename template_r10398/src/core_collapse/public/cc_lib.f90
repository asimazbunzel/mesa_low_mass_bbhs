
      module cc_lib

      use cc_def
      use star_def
      use star_lib
      use binary_def
      use binary_lib

      implicit none

      contains

      ! core-collapse controls
      subroutine cc_set_controls(inlist_fname, ierr)
         use cc_ctrls_io, only: set_default_controls, read_cc_controls
         character(len=*), intent(in) :: inlist_fname
         integer, intent(out) :: ierr

         call set_default_controls

         call read_cc_controls(inlist_fname, ierr)

      end subroutine


      ! type of core-collapse model
      subroutine cc_set_explosion_mechanism(mechanism_name, ierr)
         use cc_mechanism, only: set_explosion_mechanism
         character(*), intent(in) :: mechanism_name
         integer, intent(out) :: ierr

         call set_explosion_mechanism(mechanism_name, ierr)

      end subroutine cc_set_explosion_mechanism


      ! check if in binary
      subroutine cc_is_binary(cc_id, ierr)
         use cc_support, only: is_binary
         integer, intent(in) :: cc_id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s

         logical :: val

         call star_ptr(cc_id, s, ierr)
         if (ierr /= 0) return

         val = is_binary(s, ierr)

         cc_binary_id = 0
         if (val) cc_binary_id = 1

      end subroutine cc_is_binary

      
      ! evaluate core-collapse of star
      subroutine cc_compact_object_formation(cc_id, ierr)
         use cc_support, only: do_compact_object_formation
         integer, intent(in) :: cc_id
         integer, intent(out) :: ierr
         character(len=strlen) :: star_data_fname, binary_data_fname

         call cc_set_explosion_mechanism(model_name, ierr)

         call cc_is_binary(cc_id, ierr)

         call cc_get_star_info_pre_cc(cc_id, ierr)
         if (cc_binary_id == 1) call cc_get_binary_info_pre_cc(cc_id, cc_binary_id, ierr)

         star_data_fname = trim(cc_data_directory) // '/' // trim(filename_for_star_data) // '.data'
         binary_data_fname = trim(cc_data_directory) // '/' // trim(filename_for_binary_data) // '.data'
         call cc_write_star_info_pre_cc(cc_id, star_data_fname, ierr)
         if (cc_binary_id == 1) call cc_write_binary_info_pre_cc(cc_binary_id, binary_data_fname, ierr)

         call do_compact_object_formation(cc_id, cc_binary_id, ierr)

         call cc_write_star_info_after_cc(cc_id, star_data_fname, ierr)
         if (cc_binary_id == 1 ) then
            if (continue_binary_evolution) &
               call cc_write_binary_info_after_cc(cc_binary_id, binary_data_fname, ierr)
         end if

      end subroutine cc_compact_object_formation


      ! get info pre-cc of star
      subroutine cc_get_star_info_pre_cc(cc_id, ierr)
         use cc_data, only: get_star_info_pre_cc
         integer, intent(in) :: cc_id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s

         call star_ptr(cc_id, s, ierr)
         if (ierr /= 0) return

         call get_star_info_pre_cc(s, ierr)

      end subroutine cc_get_star_info_pre_cc

      ! get info pre-cc of binary
      subroutine cc_get_binary_info_pre_cc(cc_id, cc_binary_id, ierr)
         use cc_data, only: get_binary_info_pre_cc
         integer, intent(in) :: cc_id
         integer, intent(in) :: cc_binary_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b

         call binary_ptr(cc_binary_id, b, ierr)
         if (ierr /= 0) return

         call get_binary_info_pre_cc(cc_id, b, ierr)

      end subroutine cc_get_binary_info_pre_cc


     ! write info pre-cc of star
      subroutine cc_write_star_info_pre_cc(cc_id, filename, ierr)
         use cc_data, only: write_star_info_pre_cc
         integer, intent(in) :: cc_id
         character(len=*), intent(in) :: filename
         integer, intent(out) :: ierr

         call write_star_info_pre_cc(cc_id, filename, ierr)

      end subroutine cc_write_star_info_pre_cc
     
      ! write info pre-cc of binary
      subroutine cc_write_binary_info_pre_cc(cc_binary_id, filename, ierr)
         use cc_data, only: write_binary_info_pre_cc
         integer, intent(in) :: cc_binary_id
         character(len=*), intent(in) :: filename
         integer, intent(out) :: ierr

         call write_binary_info_pre_cc(cc_binary_id, filename, ierr)

      end subroutine cc_write_binary_info_pre_cc


      ! write info after-cc of star
      subroutine cc_write_star_info_after_cc(cc_id, filename, ierr)
         use cc_data, only: write_star_info_after_cc
         integer, intent(in) :: cc_id
         character(len=*), intent(in) :: filename
         integer, intent(out) :: ierr

         call write_star_info_after_cc(cc_id, filename, ierr)

      end subroutine cc_write_star_info_after_cc

      ! write info after-cc of binary
      subroutine cc_write_binary_info_after_cc(cc_binary_id, filename, ierr)
         use cc_data, only: write_binary_info_after_cc
         integer, intent(in) :: cc_binary_id
         character(len=*), intent(in) :: filename
         integer, intent(out) :: ierr

         call write_binary_info_after_cc(cc_binary_id, filename, ierr)

      end subroutine cc_write_binary_info_after_cc

      ! Save model of pre-cc star
      !subroutine save_pre_cc_model(id, filename, ierr)
      !   use star_lib, only: star_write_model

      !   integer, intent(in) :: id
      !   character, intent(in) :: filename
      !   integer, intent(out) :: ierr

      !   call star_write_model(id, filename, ierr)
      
      !end subroutine save_pre_cc_model

      end module cc_lib

