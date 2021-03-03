
      module ce_ctrls_io

      use ce_def
      use const_def
      use star_def
      use star_lib
      use binary_def
      use binary_lib
      
      implicit none

      namelist /ce_controls/ &
         alpha_two_stars, &
         alpha_xrb, &
         edd_scaling_factor, &
         tol_two_stars, &
         tol_xrb, &
         years_in_detachment, &
         years_to_max_mdot_rlof, &
         max_mdot_rlof, &
         add_accretion_on_ce, &
         save_profile_pre_ce, &
         save_profile_after_ce, &
         save_model_pre_ce, &
         save_model_after_ce, &
         ce_data_directory, &
         filename_donor_profile_pre_ce, &
         filename_donor_profile_after_ce, &
         filename_donor_model_pre_ce, &
         filename_donor_model_after_ce, &
         filename_accretor_profile_pre_ce, &
         filename_accretor_profile_after_ce, &
         filename_accretor_model_pre_ce, &
         filename_accretor_model_after_ce, &
         max_relative_gap, &
         max_number_retries_during_ce

      contains

      subroutine set_default_controls

         alpha_two_stars = 1d0
         alpha_xrb = 2d0
         edd_scaling_factor = 1d0
         tol_two_stars = -1d-2
         tol_xrb = -1d-1
         years_in_detachment = 1d0
         years_to_max_mdot_rlof = 7.5d0
         max_mdot_rlof = 0.1d0
         add_accretion_on_ce = .false.
         save_profile_pre_ce = .true.
         save_profile_after_ce = .true.
         save_model_pre_ce = .true.
         save_model_after_ce = .true.
         ce_data_directory = 'ce_data'
         filename_donor_profile_pre_ce = 'profile_donor_pre_ce'
         filename_donor_profile_after_ce = 'profile_donor_after_ce'
         filename_donor_model_pre_ce = 'donor_pre_ce'
         filename_donor_model_after_ce = 'donor_after_ce'
         filename_accretor_profile_pre_ce = 'profile_accretor_pre_ce'
         filename_accretor_profile_after_ce = 'profile_accretor_after_ce'
         filename_accretor_model_pre_ce = 'accretor_pre_ce'
         filename_accretor_model_after_ce = 'accretor_after_ce'
         max_relative_gap = 1d2
         max_number_retries_during_ce = 200

      end subroutine set_default_controls


      subroutine read_ce_controls(inlist_fname, ierr)
         use utils_lib, only: alloc_iounit, free_iounit
         character (len=*), intent(in) :: inlist_fname
         integer, intent(out) :: ierr
         integer :: iounit

         include 'formats.inc'

         ierr = 0

         iounit = alloc_iounit(ierr)
         open(unit=iounit, file=trim(inlist_fname), action='read', delim='quote', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(inlist_fname)
            call free_iounit(iounit)
            return
         end if

         read(iounit, nml=ce_controls, iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to read ce_controls'
            return
         end if
         close(iounit)
         call free_iounit(iounit)

      end subroutine read_ce_controls


      subroutine save_ce_profile(ce_id, step_id, ierr)
         integer, intent(in) :: ce_id
         character(len=*), intent(in) :: step_id
         integer, intent(out) :: ierr
         character(len=strlen) :: fname
         logical :: do_accretor_save
         integer :: id_extra
         integer :: ce_num

         100 format(a,i0.3,a)

         ierr = 0

         ce_num = ce_num_xrb
         ! some more info on two stars case
         do_accretor_save = .false.
         if (ce_type == ce_two_stars) then
            do_accretor_save = .true.
            ce_num = ce_num_two_stars
         end if

         if (step_id == 'pre') then
            write(fname,100) trim(ce_data_directory) // '/' // trim(filename_donor_profile_pre_ce) &
               // '_event', ce_num, '.data'
            call star_write_profile_info(ce_donor_id, fname, id_extra, ierr)
            if (ierr /= 0) return

            if (do_accretor_save) then
               write(fname,100) trim(ce_data_directory) // '/' // trim(filename_accretor_profile_pre_ce) &
                  //  '_event', ce_num, '.data'
               call star_write_profile_info(ce_accretor_id, fname, id_extra, ierr)
               if (ierr /= 0) return
            end if
         else
            write(fname,100) trim(ce_data_directory) // '/' // trim(filename_donor_profile_after_ce) &
               // '_event', ce_num, '.data'
            call star_write_profile_info(ce_donor_id, fname, id_extra, ierr)
            if (ierr /= 0) return

            if (do_accretor_save) then
               write(fname,100) trim(ce_data_directory) // '/' // trim(filename_accretor_profile_after_ce) &
                  // '_event', ce_num, '.data'
               call star_write_profile_info(ce_accretor_id, fname, id_extra, ierr)
               if (ierr /= 0) return
            end if
         end if

      end subroutine save_ce_profile

            
      subroutine save_ce_model(ce_id, step_id, ierr)
         integer, intent(in) :: ce_id
         character(len=*), intent(in) :: step_id
         integer, intent(out) :: ierr
         character(len=strlen) :: fname
         logical :: do_accretor_save
         integer :: ce_num

         100 format(a, i0.3, a)

         ierr = 0

         ce_num = ce_num_xrb
         ! some more info on two stars case
         do_accretor_save = .false.
         if (ce_type == ce_two_stars) then
            do_accretor_save = .true.
            ce_num = ce_num_two_stars
         end if
         
         if (step_id == 'pre') then
            write(fname,100) trim(ce_data_directory) // '/' // trim(filename_donor_model_pre_ce) &
               // '_event', ce_num, '.mod'
            call star_write_model(ce_donor_id, fname, ierr)
            if (ierr /= 0) return

            if (do_accretor_save) then
               write(fname,100)trim(ce_data_directory) // '/' //  trim(filename_accretor_model_pre_ce) &
                  // '_event', ce_num, '.mod'
               call star_write_model(ce_accretor_id, fname, ierr)
               if (ierr /= 0) return
            end if
         else
            write(fname,100) trim(ce_data_directory) // '/' // trim(filename_donor_model_after_ce) &
               // '_event' , ce_num, '.mod'
            call star_write_model(ce_donor_id, fname, ierr)
            if (ierr /= 0) return

            if (do_accretor_save) then
               write(fname,100) trim(ce_data_directory) // '/' // trim(filename_accretor_model_after_ce) &
                  // '_event', ce_num, '.mod'
               call star_write_model(ce_accretor_id, fname, ierr)
               if (ierr /= 0) return
            end if
         end if

      end subroutine save_ce_model

      end module ce_ctrls_io

