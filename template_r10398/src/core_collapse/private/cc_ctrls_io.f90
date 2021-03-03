
      module cc_ctrls_io

      use cc_def
      use const_def

      implicit none

      namelist /cc_controls/ &
         model_name, &
         max_ns_mass, &
         cc_data_directory, &
         filename_for_star_data, &
         filename_for_binary_data, &
         continue_binary_evolution, &
         add_asymmetric_kick

      contains

      subroutine set_default_controls

         model_name = 'rapid'
         max_ns_mass = 2.5d0
         cc_data_directory = 'cc_data'
         filename_for_star_data = 'core_collapse_star'
         filename_for_binary_data = 'core_collapse_binary'
         continue_binary_evolution = .false.
         add_asymmetric_kick = .false.

      end subroutine set_default_controls

      subroutine read_cc_controls(inlist_fname, ierr)
         use utils_lib, only: alloc_iounit, free_iounit, mkdir
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

         read(iounit, nml=cc_controls, iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to read cc_controls'
            return
         end if
         close(iounit)
         call free_iounit(iounit)

         ! after reading cc_data_directory, run `mkdir -p` on it
         call mkdir(cc_data_directory)

      end subroutine read_cc_controls

      end module cc_ctrls_io
