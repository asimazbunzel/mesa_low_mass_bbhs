      program binary_run

      use bin2dco_utils
      use bin2dco_io
      use const_def, only: dp, strlen
      use utils_def, only: min_io_unit, max_io_unit
      use utils_lib, only: alloc_iounit, free_iounit
      use star_def, only: have_initialized_star_handles
      use binary_def, only: have_initialized_binary_handles
      use binary_lib, only: run1_binary
      use run_star_support, only: MESA_INLIST_RESOLVED
      use run_star_extras
      use run_binary_extras

      implicit none

      character(len=strlen) :: inlist_fname_arg
      logical :: tst
      integer :: ierr
      integer :: iounit
      integer :: i, k
      integer :: Nsim
      character(len=strlen) :: cc_model_name
      character(len=strlen) :: str_index
      character(len=strlen) :: termination_code_star_plus_star
      logical :: has_reached_cc


      MESA_INLIST_RESOLVED = .true.
      tst = .true.


      ! get bin2dco_options
      iounit = alloc_iounit(ierr)
      if (ierr /= 0) stop 'could not alloc_iounit to get bin2dco_options namelist'
      open(unit=iounit, file='bin2dco_controls', status='old', action='read', iostat=ierr)
      if (ierr /= 0) stop 'failed to open bin2dco_controls'
      read(iounit, nml=bin2dco_options, iostat=ierr)
      if (ierr /= 0) stop 'failed to read bin2dco_options'
      close(iounit)
      call free_iounit(iounit)


      ! if file .skip_star_plus_star exist, then do not make the first part of the simulation
      iounit = alloc_iounit(ierr)
      if (ierr /= 0) stop 'could not alloc_iounit to find .skip_star_plus_star file'
      open(unit=iounit, file='.skip_star_plus_star', status='old', action='read', iostat=ierr)
      if (ierr /= 0) then

         close(iounit)
         call free_iounit(iounit)

         ! Call MESA binary
         inlist_fname_arg = star_plus_star_filename
         call run1_binary(tst, &
            ! star extras
            extras_controls, &
            ! binary extras
            extras_binary_controls, &
            ierr, &
            inlist_fname_arg)

         ! if require to stop after the star + star, then do it
         if (stop_after_star_plus_star) then
            write(*,*) 'star + star simulation ended. exit program'
            call exit()
         end if

         if (do_kicks) then
            has_reached_cc = end_as_core_collapse('termination_codes/termination_code_star_plus_star')
            if (.not. has_reached_cc) then
               write(*,*) 'star + star simulation did not reach first core-collapse. exit program'
               call exit()
            end if
         end if

      else

         close(iounit)
         call free_iounit(iounit)

      end if

      ! only do natal kick study if requested by bin2dco_options
      if (do_kicks) then

         ! get masses from files produced by core_collapse module
         call read_parameter('remnant_mass', star_info_at_cc_filename, mass_of_remnant,  ierr)
         if (ierr /= 0) stop 'failed to read_parameter remnant_mass'
         call read_parameter('progenitor_mass', binary_info_at_cc_filename, mass_of_progenitor, ierr)
         if (ierr /= 0) stop 'failed to read_parameter progenitor_mass'
         call read_parameter('companion_mass', binary_info_at_cc_filename, mass_of_companion, ierr)
         if (ierr /= 0) stop 'failed to read_parameter companion_mass'
         call read_parameter('separation_pre_cc', binary_info_at_cc_filename, pre_cc_separation, ierr)
         if (ierr /= 0) stop 'failed to read_parameter separation_pre_cc'

         write(*,'(a)')
         write(*,'(a)') 'binary components at core-collapse read from file:'
         write(*,'(a)')
         write(*,'(a32, 2x, 1pes40.16e3)') 'progenitor_mass =', mass_of_progenitor
         write(*,'(a32, 2x, 1pes40.16e3)') 'remnant_mass =', mass_of_remnant
         write(*,'(a32, 2x, 1pes40.16e3)') 'companion_mass =', mass_of_companion
         write(*,'(a32, 2x, 1pes40.16e3)') 'separation_pre_cc =', pre_cc_separation
         write(*,'(a)')

         if (do_kicks_in_one_run) then

            ! get number of simulations to perform
            Nsim = number_of_kicks(natal_kicks_filename, ierr)
            if (ierr /= 0) stop 'failed in get number_of_kicks'
               
            write(*,'(a)')
            write(*,'(a, 2x, i12)') 'the number of simulations of star and a compact object to perform is', Nsim
            write(*,'(a)')

            do i = 1, Nsim

               ! remove .restart file to avoid MESA doing a restart run
               call execute_command_line('rm -f .restart')

               ! prepare_to_apply_kick: it reads info from natal-kicks filename
               ! evaluates post-cc binary parameters, checks for binary disruption
               ! and returns the post-cc parameters to use in the run and an integer
               ! that will tell the code to do it (0), or not (-1)
               call prepare_kick(natal_kicks_filename, i, &
                  mass_of_progenitor, mass_of_remnant, mass_of_companion, pre_cc_separation, &
                  name_id, after_cc_period, after_cc_eccentricity, is_disrupted, ierr)
               if (ierr /= 0) stop 'failed in prepare_kick'

               if (is_disrupted) cycle

               ! reset initilized_*_handles to false, else after 10 runs, MESA will not continue evolving binaries
               have_initialized_binary_handles = .false.
               have_initialized_star_handles = .false.

               ! call MESA binary
               inlist_fname_arg = star_plus_pm_filename
               call run1_binary(tst, &
                  ! star extras
                  extras_controls, &
                  ! binary extras
                  extras_binary_controls, &
                  ierr, &
                  inlist_fname_arg)

               ! grab end code of second part which is in file 'tmp_second_part'
               call store_termination_code(name_id)

            end do

         else

            ! do just one natal-kick run given the argument from the command-line
            write(*,'(a)')
            write(*,'(a)') 'running only one kick'
            write(*,'(a)')

            ! remove .restart file to avoid MESA doing a restart run
            call execute_command_line('rm -f .restart')


            call get_command_argument(1, str_index, status=ierr)
            if (ierr /= 0) stop 'failed to get index of ' // trim(natal_kicks_filename)
            ! convert str_index to integer
            read(str_index,*) i
            
            ! prepare_to_apply_kick: it reads info from natal-kicks filename
            ! evaluates post-cc binary parameters, checks for binary disruption
            ! and returns the post-cc parameters to use in the run and an integer
            ! that will tell the code to do it (0), or not (-1)
            call prepare_kick(natal_kicks_filename, i, &
               mass_of_progenitor, mass_of_remnant, mass_of_companion, pre_cc_separation, &
               name_id, after_cc_period, after_cc_eccentricity, is_disrupted, ierr)
            if (ierr /= 0) stop 'failed in prepare_kick'


            if (is_disrupted) then
               write(*,'(a)') 'exit program'
               call exit()
            end if

            ! reset initilized_*_handles to false, else after 10 runs, MESA will not continue evolving binaries
            have_initialized_binary_handles = .false.
            have_initialized_star_handles = .false.

            ! call MESA binary
            inlist_fname_arg = star_plus_pm_filename
            call run1_binary(tst, &
               ! star extras
               extras_controls, &
               ! binary extras
               extras_binary_controls, &
               ierr, &
               inlist_fname_arg)

            ! grab end code of second part which is in file 'tmp_second_part'
            call store_termination_code(name_id)

         end if

      else

         ! name_id needed to store the termination code in a file which name
         ! will be: 'termination_codes/termination_code_<name_id>'
         name_id = 'star_plus_pm'

         ! do star plus pm but for a different case, .e.g, UCWD binaries
         ! for that we will use the MESA_INLIST env variable set by the
         ! script that launches the runs
         inlist_fname_arg = star_plus_pm_filename
         call run1_binary(tst, &
            ! star extras
            extras_controls, &
            ! binary extras
            extras_binary_controls, &
            ierr, &
            inlist_fname_arg)

         ! grab end code of second part which is in file 'tmp_second_part'
         call store_termination_code(name_id)

      end if


      end program binary_run
