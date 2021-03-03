! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_binary_extras

      use star_def
      use binary_def
      use const_def
      use chem_def
      use ce_def
      use cc_def

      use star_lib
      use num_lib
      use crlibm_lib
      use utils_lib
      use ce_lib
      use cc_lib

      implicit none

      ! see `bin2dco_controls` for their usage
      character(len=strlen) :: star_plus_star_filename, star_plus_pm_filename
      character(len=strlen) :: cc1_inlist_filename, cc2_inlist_filename
      character(len=strlen) :: ce1_inlist_filename, ce2_inlist_filename
      logical :: stop_after_star_plus_star
      logical :: do_kicks, do_kicks_in_one_run
      character(len=strlen) :: natal_kicks_filename
      character(len=strlen) :: star_info_at_cc_filename, binary_info_at_cc_filename
      namelist /bin2dco_options/ &
         star_plus_star_filename, star_plus_pm_filename, &
         cc1_inlist_filename, cc2_inlist_filename, &
         ce1_inlist_filename, ce2_inlist_filename, &
         stop_after_star_plus_star, &
         do_kicks, do_kicks_in_one_run, &
         natal_kicks_filename, &
         star_info_at_cc_filename, binary_info_at_cc_filename

      ! only used together with natal-kicks
      real(dp) :: mass_of_progenitor
      real(dp) :: mass_of_remnant, mass_of_companion
      real(dp) :: pre_cc_separation
      character(len=strlen) :: name_id
      real(dp) :: after_cc_period, after_cc_eccentricity
      logical :: is_disrupted

      ! utilities
      integer :: num_switches
      logical :: second_collapse

      logical :: dbg = .false.

      contains
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! Set these function pinters to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
         b% warn_binary_extra =.false.
         
      end subroutine extras_binary_controls


      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id

         how_many_extra_binary_history_columns = 10

      end function how_many_extra_binary_history_columns
     

      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         use const_def, only: dp
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: zeta_1, zeta_2, zeta_rl_1, zeta_rl_2

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! common-envelope variables
         names(1) = 'ce_phase'
         names(2) = 'cumulative_binding_energy'
         names(3) = 'ce_num_star_1'
         names(4) = 'ce_num_star_2' 
         if (ce_on) then
            vals(1) = 1d0
         else
            vals(1) = 0d0
         end if
         vals(2) = cumulative_removed_binding_energy
         vals(3) = ce_num_star_1
         vals(4) = ce_num_star_2

         ! switches during almost-overcontact phase
         names(5) = 'num_switches'
         vals(5) = num_switches

         ! zetas exponents for mass-transfer stability study
         names(6) = 'zeta_1'
         names(7) = 'zeta_2'
         names(8) = 'zeta_rl_1'
         names(9) = 'zeta_rl_2'
         
         ! evaluate mass-radius exponents here using *_older copies as the *_old already have been
         ! changed to the values obtained in the last timestep
         if (b% point_mass_i /= 1) then
            zeta_1 = ((b% m(1)/Msun) / (b% r(1)/Rsun)) * &
               ((b% r(1)/Rsun) - (b% r_older(1)/Rsun)) / ((b% m(1)/Msun) - (b% m_older(1)/Msun))
         else
            zeta_1 = -99d0
         end if
         if (b% point_mass_i /= 2) then
            zeta_2 = ((b% m(2)/Msun) / (b% r(2)/Rsun)) * &
               ((b% r(2)/Rsun) - (b% r_older(2)/Rsun)) / ((b% m(2)/Msun) - (b% m_older(2)/Msun))
         else
            zeta_2 = -99d0
         end if
         zeta_rl_1 = ((b% m(1)/Msun) / (b% rl(1)/Rsun)) * &
            ((b% rl(1)/Rsun) - (b% rl_older(1)/Rsun)) / ((b% m(1)/Msun) - (b% m_older(1)/Msun))
         zeta_rl_2 = ((b% m(2)/Msun) / (b% rl(2)/Rsun)) * &
            ((b% rl(2)/Rsun) - (b% rl_older(2)/Rsun)) / ((b% m(2)/Msun) - (b% m_older(2)/Msun))

         vals(6) = zeta_1
         vals(7) = zeta_2
         vals(8) = zeta_rl_1
         vals(9) = zeta_rl_2

         ! outer lagrangian point equivalent radii
         names(10) = 'rl2_donor'
         vals(10) = eval_outer_roche_lobe(b% m(b% d_i), b% m(b% a_i), b% separation) / Rsun

      end subroutine data_for_extra_binary_history_columns


      real function eval_outer_roche_lobe(donor_mass, accretor_mass, separation) result(rl2)
         real(dp), intent(in) :: donor_mass, accretor_mass, separation
         real(dp) :: q, q13

         ! Based on Eggleton 2006 'Evolutionary processes in binary and multiple stars'
         q = donor_mass / accretor_mass
         q13 = pow_cr(q,one_third)
         if (q > 1d0) then
            rl2 = 0.49d0 * q13*q13 + 0.15d0
         else
            rl2 = 0.49d0 * q13*q13 + 0.27d0 * q - 0.12d0 * q13*q13*q13*q13
         end if

         rl2 = (separation * rl2) / (0.6d0 * q13*q13 + log1p_cr(q13))
      end function eval_outer_roche_lobe
      
      
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         integer :: iounit
         
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in binary_ptr
            return
         end if

         ! initialize variables only used inside the run_binary_extras module
         num_switches = 0
         second_collapse = .false.

         ! read bin2dco_controls
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) stop 'could not alloc_iounit to get bin2dco_options namelist'
         open(unit=iounit, file='bin2dco_controls', status='old', action='read', iostat=ierr)
         if (ierr /= 0) stop 'failed to open bin2dco_controls'
         read(iounit, nml=bin2dco_options, iostat=ierr)
         if (ierr /= 0) stop 'failed to read bin2dco_options'
         close(iounit)
         call free_iounit(iounit)

         ! set initial variables of ce
         if (b% point_mass_i == 0) then
            call ce_variables_on_startup(ce1_inlist_filename, ierr)
            if (ierr /= 0) return
         else
            call ce_variables_on_startup(ce2_inlist_filename, ierr)
            if (ierr /= 0) return
            if (do_kicks) then
               write(filename_donor_profile_pre_ce, '(a)') trim(filename_donor_profile_pre_ce) // "_" // trim(name_id)
               write(filename_donor_profile_after_ce, '(a)') trim(filename_donor_profile_after_ce) // "_" // trim(name_id)
               write(filename_donor_model_pre_ce, '(a)') trim(filename_donor_model_pre_ce) // "_" // trim(name_id)
               write(filename_donor_model_after_ce, '(a)') trim(filename_donor_model_after_ce) // "_" // trim(name_id)
            end if
         end if

         ! we want to always have the same model_number for binary & star
         if (b% model_number /= b% s_donor% model_number .or. b% binary_age /= b% s_donor% star_age) then
            b% model_number = b% s_donor% model_number
            b% model_number_old = b% model_number
            b% model_number_older = b% model_number
            b% binary_age = b% s_donor% star_age
            b% binary_age_old = b% binary_age
            b% binary_age_older = b% binary_age
         end if

         ! This is a patch that needs to be set, else we cannot 
         ! do the evolution of a binary from star + star and up to
         ! the formation of a DCO for many cases because there is
         ! an issue in the binary/private/run_binary_support.f90 file
         ! in which around line 695 a do loop is done between 1 and 2
         ! when it should be between 1 and num_stars
         if (b% point_mass_i /= 0) then
            b% star_ids(2) = -1
            b% star_extra_ids(2) = -1
         end if

         ! update binary parameters based on the core-collapse data and natal-kick
         if (b% point_mass_i /= 0 .and. do_kicks) then

            b% m(b% a_i) = mass_of_remnant * Msun  ! cgs
            b% initial_period_in_days = after_cc_period
            b% initial_eccentricity = after_cc_eccentricity
            call binary_set_period_eccentricity(binary_id, &
               b% initial_period_in_days*(24d0*60d0*60d0), b% initial_eccentricity)
            write(b% history_name, '(a)') 'binary_history_' // trim(name_id) // '.data'
            write(b% s_donor% star_history_name, '(a)') 'secondary_history_' // trim(name_id) // '.data'
            
            if (dbg) then
               write(*,'(a32, 2x, 1pes40.16e3)') 'm2 =', b% m(b% a_i) / Msun
               write(*,'(a32, 2x, 1pes40.16e3)') 'period =', b% period / (24d0*60d0*60d0)
               write(*,'(a32, 2x, 1pes40.16e3)') 'eccentricity =', b% eccentricity
               write(*,'(a32, 2x, a)') 'binary history name =', trim(b% history_name)
               write(*,'(a32, 2x, a)') 'star history name =', trim(b% s_donor% star_history_name)
            end if

         end if

         extras_binary_startup = keep_going
      
      end function  extras_binary_startup


      !Return either rety,backup,keep_going or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr, id_extra

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         extras_binary_check_model = keep_going

         ! check if ce
         if (b% r(b% d_i) > b% rl(b% d_i) .and. ce_off) then
            call ce_unstable_mt_phase(binary_id, ierr)
            if (ierr /= 0) return

            if (ce_on) call ce_init(binary_id, ierr)
            if (ierr /= 0) return
         end if

         ! check ce end
         if (ce_on) then
            call ce_check_state(binary_id, ierr)
            if (ce_merge) then
               extras_binary_check_model = terminate
               return
            end if
         end if

      end function extras_binary_check_model

      logical function switch_donor_star(binary_id) result(do_switch)
         ! HINT: Switching donor will only work by setting
         !       keep_donor_fixed = .false. in &binary_controls
         !       within "inlist"  and adding these lines in
         !       run_binary_extras.f90
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         real(dp) :: F1, q, rho, p, grav, hp, v_th, rl3, q_temp
         real(dp) :: mdot_thin_donor, mdot_thin_accretor
         
         do_switch = .false.

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         if (b% point_mass_i == 0  .and. .not. b% keep_donor_fixed .and. b% r(b% d_i) < b% rl(b% d_i)) then
            !--------------------- Optically thin MT rate --------------------------
            ! As described in H. Ritter 1988, A&A 202,93-100 and
            ! U. Kolb and H. Ritter 1990, A&A 236,385-392
            rho = b% s_donor% rho(1)  ! density at surface in g/cm^3
            p = b% s_donor% p(1)  ! pressure at surface in dynes/cm^2
            grav = b% s_donor% cgrav(1)*b% m(b% d_i)/(b% r(b% d_i))**2  ! local gravitational
            ! acceleration
            hp = p/(grav*rho)  ! pressure scale height
            v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1)))
            q = b% m(b% a_i)/b% m(b% d_i)  ! Mass ratio, as defined in Ritter 1988

            ! (Kolb & Ritter 1990 use the opposite!)
            ! consider range of validity for F1, do not extrapolate Eq. A9 of Ritter 1988
            q_temp = min(max(q,0.5d0),10d0)
            F1 = (1.23d0  + 0.5d0* log10_cr(q_temp))
            rl3 = (b% rl(b% d_i))*(b% rl(b% d_i))*(b% rl(b% d_i))
            mdot_thin_donor = (2.0d0*pi/exp_cr(0.5d0)) * v_th*v_th*v_th * &
               rl3/(b% s_donor% cgrav(1)*b% m(b% d_i)) * rho * F1

            rho = b% s_accretor% rho(1)  ! density at surface in g/cm^3
            p = b% s_accretor% p(1)  ! pressure at surface in dynes/cm^2
            grav = b% s_accretor% cgrav(1)*b% m(b% a_i)/(b% r(b% a_i))**2  ! local gravitational acceleration
            hp = p/(grav*rho) ! pressure scale height
            v_th = sqrt(kerg * b% s_accretor% T(1) / (mp * b% s_accretor% mu(1)))
            q = b% m(b% d_i)/b% m(b% a_i)  ! Mass ratio, as defined in Ritter 1988

            ! (Kolb & Ritter 1990 use the opposite!)
            ! consider range of validity for F1, do not extrapolate! Eq. A9 of Ritter 1988
            q_temp = min(max(q,0.5d0),10d0)
            F1 = (1.23d0  + 0.5d0* log10_cr(q_temp))
            rl3 = (b% rl(b% a_i))*(b% rl(b% a_i))*(b% rl(b% a_i))
            mdot_thin_accretor = (2.0d0*pi/exp_cr(0.5d0)) * v_th*v_th*v_th * &
               rl3/(b% s_accretor% cgrav(1)*b% m(b% a_i)) * rho * F1
         
            if (abs(mdot_thin_accretor) < abs(mdot_thin_donor)) do_switch = .true.
         end if

      end function switch_donor_star

      
      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         logical :: do_switch
         integer :: star_cc_id
         integer :: number_io
         real(dp), parameter :: chandra_mass = 1.4d0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         extras_binary_finish_step = keep_going

         ! if mass-transfer is really high even though no RLOF, then terminate
         if (b% r(b% d_i) < b% rl(b% d_i) .and. ce_off &
            .and. abs(b% mtransfer_rate) * secyer/Msun > max_mdot_rlof) then
            write(*,'(a)') 'reach a really high MT rate without having RLOF'
            b% s_donor% termination_code = t_xtra1
            termination_code_str(t_xtra1) = 'mdot_atmospheric > max_mdot_rlof'
            extras_binary_finish_step = terminate
            return
         end if

         ! take success step on for ce vars
         if (ce_on) then
            call ce_step_success(binary_id, ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'failed in ce_step_success'
               return
            end if
         end if


         ! first collapse
         star_cc_id = 0
         if (b% point_mass_i == 0) then
            if (b% s1% center_c12 < 1d-3 * b% s1% initial_z .and. b% s1% center_he4 < 1d-6) then
               star_cc_id = 1
               call star_write_model(2, 'companion_at_core_collapse.mod', ierr)
               b% s1% termination_code = t_xtra1
               termination_code_str(t_xtra1) = 'core-collapse'
               extras_binary_finish_step = terminate
            else if (b% s2% center_c12 < 1d-3 * b% s2% initial_z .and. b% s2% center_he4 < 1d-6) then
               star_cc_id = 2
               call star_write_model(1, 'companion_at_core_collapse.mod', ierr)
               b% s2% termination_code = t_xtra1
               termination_code_str(t_xtra1) = 'core-collapse'
               extras_binary_finish_step = terminate
            end if
         ! second collapse
         else if (b% point_mass_i == 1) then
            if (b% s2% center_c12 < 1d-3 * b% s2% initial_z .and. b% s2% center_he4 < 1d-6) then
               second_collapse = .true.
               star_cc_id = 2
               b% s2% termination_code = t_xtra1
               termination_code_str(t_xtra1) = 'core-collapse'
               extras_binary_finish_step = terminate
            end if
         else if (b% point_mass_i == 2) then
            if (b% s1% center_c12 < 1d-3 * b% s1% initial_z .and. b% s1% center_he4 < 1d-6) then
               second_collapse = .true.
               star_cc_id = 1
               b% s1% termination_code = t_xtra1
               termination_code_str(t_xtra1) = 'core-collapse'
               extras_binary_finish_step = terminate
            end if
         end if

         if (star_cc_id > 0 .and. .not. second_collapse) then
            call cc_set_controls(cc1_inlist_filename, ierr)
            if (do_kicks) then
               write(filename_for_star_data, '(a)') trim(filename_for_star_data) // "_1"
               write(filename_for_binary_data, '(a)') trim(filename_for_binary_data) // "_1"
            end if
            if (ierr /= 0) return
            call cc_compact_object_formation(star_cc_id, ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'failed in cc_compact_object_formation'
               return
            end if
            return
         else if (star_cc_id > 0 .and. second_collapse) then
            write(*,'(a)') 'calling second collapse'
            if (do_kicks) then
               call cc_set_controls(cc2_inlist_filename, ierr)
               if (ierr /= 0) return
               ! replace filenames by adding natal-kick id
               write(filename_for_star_data, '(a)') trim(filename_for_star_data) // "_" // trim(name_id)
               write(filename_for_binary_data, '(a)') trim(filename_for_binary_data) // "_" // trim(name_id)
               call cc_compact_object_formation(star_cc_id, ierr)
               if (dbg) then
                  write(*,'(a32, 2x, a)') 'filename_for_star_data =', trim(filename_for_star_data)
                  write(*,'(a32, 2x, a)') 'filename_for_binary_data =', trim(filename_for_binary_data)
               end if
            else
               call cc_set_controls(cc1_inlist_filename, ierr)
               if (do_kicks) then
                  write(filename_for_star_data, '(a)') trim(filename_for_star_data) // "_1"
                  write(filename_for_binary_data, '(a)') trim(filename_for_binary_data) // "_1"
               end if
               if (ierr /= 0) return
               call cc_compact_object_formation(star_cc_id, ierr)
            end if
            if (ierr /= 0) then
               write(*,'(a)') 'failed in cc_compact_object_formation'
               return
            end if
            return
         end if
         
         ! check if we need to switch donor star only when ce is off
         if (ce_off) then
            do_switch = switch_donor_star(binary_id)
            if (do_switch) then
               write(*,'(a)') "switching donor"
               num_switches = num_switches + 1
               if (b% d_i == 2) then
                  b% d_i = 1
                  b% d_i_old = 1
                  b% d_i_older = 1
                  b% a_i = 2
                  b% a_i_old = 2
                  b% a_i_older = 2
                  b% s_donor => b% s1
                  b% s_accretor => b% s2
               else
                  b% d_i = 2
                  b% d_i_old = 2
                  b% d_i_older = 2
                  b% a_i = 1
                  b% a_i_old = 1
                  b% a_i_older = 1
                  b% s_donor => b% s2
                  b% s_accretor => b% s1
               end if
            end if
         end if

      end function extras_binary_finish_step

      
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         integer :: ios
         character(len=strlen) :: fname, termination_code

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in binary_ptr
            return
         end if

         if (b% point_mass_i == 0) then
            if (b% s1% termination_code > 0) then
               termination_code = trim(termination_code_str(b% s1% termination_code))
            else if (b% s2% termination_code > 0) then
               termination_code = trim(termination_code_str(b% s2% termination_code))
            else
               termination_code = 'unknown'
            end if
         else
            if (b% s1% termination_code > 0) then
               termination_code = trim(termination_code_str(b% s1% termination_code))
            else if (b% s2% termination_code > 0) then
               termination_code = trim(termination_code_str(b% s2% termination_code))
            else
               termination_code = 'unknown'
            end if
         end if

         ! run mkdir -p `termination_codes`
         call mkdir('termination_codes')

         fname = ''
         if (b% point_mass_i == 0) then
            fname = 'termination_codes/termination_code_star_plus_star'
         else
            fname = 'tmp_second_part'
         end if

         if (fname /= '') then
            open(unit=22, file=trim(fname), iostat=ios)
            if (ios /= 0) stop 'error opening file ' // trim(fname) 
            write(unit=22, fmt='(a)', iostat=ios, advance='no') trim(termination_code)
            if (ios /= 0) stop 'write error in file unit 22'
            close(22)
         else
            write(*,*) 'failed to open', fname, 'to write end code'
            call execute_command_line('echo unknown >> tmp_second_part')
         end if

      end subroutine extras_binary_after_evolve     
      
      end module run_binary_extras
