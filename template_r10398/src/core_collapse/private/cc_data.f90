
      module cc_data

      use cc_def
      use const_def
      use star_def
      use binary_def

      implicit none

      public :: do_star_update_after_cc, do_binary_update_after_cc

      contains

      subroutine do_star_update_after_cc(s, ierr)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr

         include 'formats.inc'

         ierr = 0

         s% mstar = M_remnant
         s% photosphere_r = 0d0  ! point-mass approximation

      end subroutine do_star_update_after_cc

      subroutine do_binary_update_after_cc(cc_binary_id, cc_id, ierr)
         integer, intent(in) :: cc_binary_id, cc_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b

         ierr = 0
         call binary_ptr(cc_binary_id, b, ierr)
         if (ierr /= 0) return
   
         b% point_mass_i = cc_id
         if (cc_id == 1) then
            b% d_i = 2
            b% a_i = 1
            b% s_donor => b% s2
            b% s_accretor => b% s1
         else
            b% d_i = 1
            b% a_i = 2
            b% s_donor => b% s1
            b% s_accretor => b% s2
         end if

         b% m(cc_id) = M_remnant * Msun

      end subroutine do_binary_update_after_cc

      subroutine get_star_info_pre_cc(s, ierr)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr

         ierr = 0

         model_number_pre_cc = s% model_number
         age_pre_cc = s% star_age
         M_pre_cc = s% star_mass
         he_core_mass_pre_cc = s% he_core_mass
         he_core_radius_pre_cc = s% he_core_radius
         c_core_mass_pre_cc = s% c_core_mass
         c_core_radius_pre_cc = s% c_core_radius
         R_pre_cc = s% r(1) / Rsun
         Teff_pre_cc = s% Teff
         L_pre_cc = s% photosphere_L

      end subroutine get_star_info_pre_cc

      subroutine get_binary_info_pre_cc(cc_id, b, ierr)
         integer, intent(in) :: cc_id
         type(binary_info), pointer :: b
         integer, intent(out) :: ierr

         ierr = 0

         if (cc_id == 1) then
            progenitor_mass = b% m(1) / Msun
            companion_mass = b% m(2) / Msun
         else
            progenitor_mass = b% m(2) / Msun
            companion_mass = b% m(1) / Msun
         end if

         angular_momentum_pre_cc = b% angular_momentum_j
         separation_pre_cc = b% separation / Rsun
         period_pre_cc = b% period / (60d0*60d0*24d0)
         eccentricity_pre_cc = b% eccentricity
         radius_pre_cc(1) = b% r(1) / Rsun
         radius_pre_cc(2) = b% r(2) / Rsun
         rl_pre_cc(1) = b% rl(1) / Rsun
         rl_pre_cc(2) = b% rl(2) / Rsun
         mt_rate_pre_cc = b% mtransfer_rate * secyer / Msun

      end subroutine get_binary_info_pre_cc


      subroutine write_star_info_pre_cc(cc_id, filename, ierr)
         use utils_lib, only: alloc_iounit, free_iounit
         integer, intent(in) :: cc_id
         character(len=*) :: filename
         integer, intent(out) :: ierr
         integer :: iounit

         1 format(a32, 2x, 1pes40.16e3)
         2 format(a32, 2x, i40)
         3 format(a32, 3x, a8)
         4 format(a32, 19x, a40)

         ierr = 0

         iounit = alloc_iounit(ierr)
         open(unit=iounit, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(filename)
            call free_iounit(iounit)
            return
         end if

         write(iounit,2) 'model_number_pre_cc', model_number_pre_cc
         write(iounit,1) 'age_pre_cc', age_pre_cc
         write(iounit,1) 'mass_pre_cc', M_pre_cc
         write(iounit,1) 'he_core_mass_pre_cc', he_core_mass_pre_cc
         write(iounit,1) 'c_core_mass_pre_cc', c_core_mass_pre_cc
         write(iounit,1) 'he_core_radius_pre_cc', he_core_radius_pre_cc
         write(iounit,1) 'c_core_radius_pre_cc', c_core_radius_pre_cc
         write(iounit,1) 'star_radius_pre_cc', R_pre_cc
         write(iounit,1) 'Teff_pre_cc', Teff_pre_cc
         write(iounit,1) 'L_pre_cc', L_pre_cc

         close(iounit)
         call free_iounit(iounit)

      end subroutine write_star_info_pre_cc


      subroutine write_binary_info_pre_cc(cc_binary_id, filename, ierr)
         use utils_lib, only: alloc_iounit, free_iounit
         integer, intent(in) :: cc_binary_id
         character(len=*) :: filename
         integer, intent(out) :: ierr
         integer :: iounit

         1 format(a32, 2x, 1pes40.16e3)
         2 format(a32, 2x, i40)
         3 format(a32, 3x, a8)
         4 format(a32, 19x, a40)

         ierr = 0

         iounit = alloc_iounit(ierr)
         open(unit=iounit, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(filename)
            call free_iounit(iounit)
            return
         end if

         write(iounit,1) 'progenitor_mass', progenitor_mass
         write(iounit,1) 'companion_mass', companion_mass
         write(iounit,1) 'orbital_angular_momentum_pre_cc', angular_momentum_pre_cc
         write(iounit,1) 'period_pre_cc', period_pre_cc
         write(iounit,1) 'separation_pre_cc', separation_pre_cc
         write(iounit,1) 'r_1_pre_cc', radius_pre_cc(1)
         write(iounit,1) 'r_2_pre_cc', radius_pre_cc(2)
         write(iounit,1) 'rl_1_pre_cc', rl_pre_cc(1)
         write(iounit,1) 'rl_2_pre_cc', rl_pre_cc(2)
         write(iounit,1) 'mt_rate_pre_cc', mt_rate_pre_cc

         close(iounit)
         call free_iounit(iounit)

      end subroutine write_binary_info_pre_cc


      subroutine write_star_info_after_cc(cc_id, filename, ierr)
         use utils_lib, only: alloc_iounit, free_iounit
         integer, intent(in) :: cc_id
         character(len=*) :: filename
         integer, intent(out)  :: ierr
         integer :: iounit

         1 format(a32, 2x, 1pes40.16e3)
         2 format(a32, 2x, i40)
         3 format(a32, 3x, a8)
         4 format(a32, 19x, a40)

         ierr = 0
         
         iounit = alloc_iounit(ierr)
         open(unit=iounit, file=trim(filename), action='write', position='append', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(filename)
            call free_iounit(iounit)
            return
         end if

         write(iounit,4) 'sn_model', model_name
         write(iounit,1) 'max_ns_mass', max_ns_mass
         write(iounit,1) 'baryonic_mass', M_baryonic
         write(iounit,1) 'remnant_mass', M_remnant
         write(iounit,1) 'ejected_mass', M_ejected
         write(iounit,1) 'fallback_mass', M_fallback
         write(iounit,1) 'fallback_fraction', fallback_fraction

         close(iounit)
         call free_iounit(iounit)

      end subroutine write_star_info_after_cc


      subroutine write_binary_info_after_cc(cc_binary_id, filename, ierr)
         use utils_lib, only: alloc_iounit, free_iounit
         integer, intent(in) :: cc_binary_id
         character(len=*) :: filename
         integer, intent(out)  :: ierr
         integer :: iounit

         1 format(a32, 2x, 1pes40.16e3)
         2 format(a32, 2x, i40)
         3 format(a32, 3x, a8)
         4 format(a32, 19x, a40)

         ierr = 0
         
         iounit = alloc_iounit(ierr)
         open(unit=iounit, file=trim(filename), action='write', position='append', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(filename)
            call free_iounit(iounit)
            return
         end if
         
         write(iounit,1) 'orbital_angular_momentum_after_cc', angular_momentum_after_cc
         write(iounit,1) 'separation_after_cc', separation_after_cc
         write(iounit,1) 'rl_1_after_cc', rl_after_cc(1)
         write(iounit,1) 'rl_2_after_cc', rl_after_cc(2)
         write(iounit,1) 'mt_rate_after_cc', mt_rate_after_cc

         close(iounit)
         call free_iounit(iounit)

      end subroutine write_binary_info_after_cc


      end module cc_data
