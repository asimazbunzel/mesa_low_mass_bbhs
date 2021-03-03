
      module bin2dco_utils

      implicit none

      integer :: max_Nsim = 100000

      contains

      ! check version of codes
      subroutine check_version_numbers
         use const_def, only: strlen
         use utils_lib, only: alloc_iounit, free_iounit
         integer :: iounit, ierr
         character(len=strlen) :: mesa_dirname
         integer :: mesa_version_used, mesa_version_supported
         logical :: dbg = .false.

         ! get version number from $MESA_DIR
         call get_environment_variable('MESA_DIR', mesa_dirname)
         iounit = alloc_iounit(ierr)
         open(unit=iounit, file=trim(mesa_dirname) // '/data/version_number', status='old', action='read', iostat=ierr)
         read(iounit, *, iostat=ierr) mesa_version_used
         if (ierr /= 0) then
            write(*,'(a)') 'failed to open ' // trim(mesa_dirname) // '/data/version_number'
            close(iounit)
            call free_iounit(ierr)
         end if

         ! get version number from bin2dco
         open(unit=iounit, file='version_number', status='old', action='read', iostat=ierr)
         read(iounit, *, iostat=ierr) mesa_version_supported
         if (ierr /= 0) then
            write(*,'(a)') 'failed to open version_number'
            close(iounit)
            call free_iounit(ierr)
         end if

         if (dbg) then
            write(*,'(a32, 2x, i12)') 'MESA version:', mesa_version_used
            write(*,'(a32, 2x, i12)') 'bin2dco version:', mesa_version_supported
         end if

         ! check that they match
         if (mesa_version_used /= mesa_version_supported) then
            stop 'MESA version used does not match supported one'
         end if

      end subroutine check_version_numbers

      ! search for string in file and return value associated to it
      subroutine read_parameter(string, filename, val, ierr)
         use const_def, only: dp, strlen
         use utils_lib, only: alloc_iounit, free_iounit, number_iounits_allocated
         character(len=*) :: string
         character(len=*) :: filename
         real(dp), intent(out) :: val
         integer, intent(out) :: ierr
         character(len=strlen) :: buffer, trimmed_buffer, label
         integer :: pos
         logical :: dbg = .false.
         integer :: iounit

         ! allocate file unit
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) stop 'could not alloc_iounit to read_parameter'
         open(unit=iounit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(filename)
            call free_iounit(iounit)
            return
         end if

         do while (ierr == 0)
            read(iounit, '(A)', iostat=ierr) buffer
            trimmed_buffer = trim(adjustl(buffer))

            ! Find the first instance of whitespace.  Split label and data.
            pos = index(trimmed_buffer,' ')
            label = trim(adjustl(trimmed_buffer(1:pos)))
            buffer = trim(adjustl(trimmed_buffer(pos+1:)))

            if (label == trim(string)) then
               read(buffer, *, iostat=ierr) val
               if (dbg) write(*, '(a32, 2x, 1pes40.16e3)') 'Read ' // trim(label) // ': ', val
               exit
            end if
         end do

         close(iounit)
         call free_iounit(iounit)

      end subroutine read_parameter


      ! get number of natal kicks saved in the natal-kick file
      integer function number_of_kicks(filename, ierr)
         use utils_lib, only: alloc_iounit, free_iounit
         character(len=*) :: filename
         integer, intent(out) :: ierr
         integer :: iounit
         integer :: Nsim

         Nsim = 0
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) stop 'could not alloc_iounit to get number of simulations after first collapse'
         open(unit=iounit, file=trim(filename), status='old', action='read', iostat=ierr)
         if (ierr /= 0) stop 'could not open ' // trim(filename)
         do
            if (Nsim > max_Nsim) then
               write(*,'(a, 2x, i12)') 'reach max number of kicks that can be simulated of', Nsim
               stop 'modify value of max_Nsim on bin2dc0_utils.f and re-compile to have more kicks'
            end if

            read(iounit,*,iostat=ierr)
            if (ierr /= 0) exit
            Nsim = Nsim + 1
         end do
         ierr = 0
         close(iounit)
         call free_iounit(iounit)

         number_of_kicks = Nsim

      end function number_of_kicks


      ! read from file natal kick information: id, velocity, angles
      subroutine read_natal_kick(filename, line_number, name_id, kick, theta, phi, ierr)
         use const_def, only: dp, strlen
         use utils_lib, only: alloc_iounit, free_iounit
         character(len=*) :: filename
         character(len=*), intent(out) :: name_id
         integer, intent(in) :: line_number
         real(dp), intent(out) :: kick
         real(dp), intent(out) :: theta
         real(dp), intent(out) :: phi
         integer, intent(out) :: ierr

         integer :: iounit
         integer :: i, j
         character(len=strlen), dimension(4) :: str_entry
         character(len=strlen) :: buffer, trimmed_buffer, string
         logical :: dbg = .false.

         ! allocate file unit
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) stop 'could not alloc_iounit to read_natal_kick'
         open(unit=iounit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(filename)
            call free_iounit(iounit)
            return
         end if
   
         i = 1
         do while (ierr == 0)

            ! skip first (line_number - 1) lines
            do
               if (i > line_number - 1) exit
               read(iounit, *, iostat=ierr)
               if (ierr > 0) then
                  write(*,*) 'problem in reading natal-kick file'
                  return
               else if (ierr < 0) then
                  write(*,*) 'end of file reached'
                  return
               else
                  i = i + 1
               end if
            end do

            read(iounit, '(A)', iostat=ierr) buffer
            read(buffer,*) (str_entry(i), i=1,4)

            do j = 1, 4
               string = trim(adjustl(str_entry(j)))
               select case (j)
               case (1)
                  read(string, *, iostat=ierr) name_id
                  if (dbg) write(*, '(a32, 2x, a40)') 'Read: ', name_id
               case (2)
                  read(string, *, iostat=ierr) kick
                  if (dbg) write(*, '(a32, 2x, 1pes40.16e3)') 'Read: ', kick
               case (3)
                  read(string, *, iostat=ierr) theta
                  if (dbg) write(*, '(a32, 2x, 1pes40.16e3)') 'Read: ', theta
               case (4)
                  read(string, *, iostat=ierr) phi
                  if (dbg) write(*, '(a32, 2x, 1pes40.16e3)') 'Read: ', phi
               end select
            end do
            ierr = -1
         end do

         ierr = 0
         close(iounit)
         call free_iounit(iounit)

      end subroutine read_natal_kick


      ! calculate orbital parameters for a post core-collapse binary
      subroutine binary_parameters_post_cc(a_pre_cc, mass_pre_cc, mass_after_cc, companion_mass, &
            vk, theta, phi, porbf, ef)
         use const_def, only: dp, pi
         use crlibm_lib, only: cos_cr, sin_cr, pow2, pow3
         real(dp), intent(in) :: a_pre_cc
         real(dp), intent(in) :: mass_pre_cc
         real(dp), intent(in) :: mass_after_cc
         real(dp), intent(in) :: companion_mass
         real(dp), intent(in) :: vk, theta, phi
         real(dp), intent(out) :: porbf, ef

         real(dp) :: standard_cgrav, Msun, Rsun
         real(dp) :: ai, af, m1i, m1f, m2, w
         real(dp) :: v_pre
         real(dp) :: wx, wy, wz

         ! set constants
         standard_cgrav = 6.67428d-8
         Msun = 1.9892d33
         Rsun = 6.9598d10

         ! change to cgs units
         ai = a_pre_cc * Rsun
         m1i = mass_pre_cc * Msun
         m1f = mass_after_cc * Msun
         m2 = companion_mass * Msun
         w = vk * 1d5

         ! project kick to (x,y,z)
         wx = w * cos_cr(phi) * sin_cr(theta)
         wy = w * cos_cr(theta)
         wz = w * sin_cr(phi) * sin_cr(theta)

         ! velocity pre_cc in cgs (assuming circular orbit)
         v_pre = sqrt(standard_cgrav * (m1i + m2) / ai)

         ! separation after_cc in cgs
         af = standard_cgrav * (m1f + m2) / ((2d0 * standard_cgrav * (m1f + m2) / ai) &
            - pow2(w) - pow2(v_pre) - (2d0 * wy * v_pre))

         ! eccentricity
         ef = sqrt(1d0 - (pow2(wz) + pow2(wy) + pow2(v_pre) + 2d0 * wy * v_pre) * &
            pow2(ai) / (standard_cgrav * (m1f + m2) * af))

         porbf = (2*pi) * sqrt(pow3(af) / (standard_cgrav * (m1f + m2)))
         porbf = porbf / 86400d0  ! in days

      end subroutine binary_parameters_post_cc


      ! prepare run to add a natal kick
      subroutine prepare_kick(filename, k, &
            mass_of_progenitor, mass_of_remnant, mass_of_companion, pre_cc_separation, &
            name_id, after_cc_period, after_cc_eccentricity, is_disrupted, ierr)
         use const_def, only: dp
         use utils_lib, only: alloc_iounit, free_iounit
         character(len=*) :: filename
         integer, intent(in) :: k
         real(dp), intent(in) :: mass_of_progenitor, mass_of_remnant, mass_of_companion
         real(dp), intent(in) :: pre_cc_separation
         character(len=*), intent(out) :: name_id
         real(dp), intent(out) :: after_cc_period, after_cc_eccentricity
         logical, intent(out) :: is_disrupted
         integer, intent(out) :: ierr

         integer :: iounit
         real(dp) :: kick, theta, phi

         is_disrupted = .false.

         call read_natal_kick(filename, k, name_id, kick, theta, phi, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'failed in read_natal_kick'
            return
         end if

         write(*,'(a)') 'natal kicks parameters to be used in next MESA run:'
         write(*,'(a)')
         write(*,'(a32, 19x, a40)') 'name_id =', name_id
         write(*,'(a32, 2x, 1pes40.16e3)') 'kick =', kick
         write(*,'(a32, 2x, 1pes40.16e3)') 'theta =', theta
         write(*,'(a32, 2x, 1pes40.16e3)') 'phi =', phi
         write(*,'(a)')
         
         ! update of parameters after core-collapse
         call binary_parameters_post_cc(pre_cc_separation, mass_of_progenitor, mass_of_remnant, &
            mass_of_companion, kick, theta, phi, after_cc_period, after_cc_eccentricity)

         ! if binary is disrupted, then skip simulation and write name_id in file
         if (after_cc_eccentricity < 0d0 .or. after_cc_eccentricity > 1d0) then
            write(*,'(a)') 'binary disrupted after core-collapse'
            write(*,'(a)')
            write(*,'(a32, 2x, 1pes40.16e3)') 'period_after_cc =', after_cc_period
            write(*,'(a32, 2x, 1pes40.16e3)') 'eccentricity_after_cc =', after_cc_eccentricity
            write(*,'(a)')
            ! save is_disrupted now
            is_disrupted = .true.
            ! write natal-kick id into a disrupted file
            iounit = alloc_iounit(ierr)
            if (ierr /= 0) stop 'could not alloc_iounit for disrupted cases'
            open(unit=iounit, file='disrupted_ids.data', action='write', position='append', iostat=ierr)
            if (ierr /= 0) stop 'could not open disrupted_ids.data'
            write(iounit,*) name_id
            close(iounit)
            call free_iounit(iounit)
            ! cycle
         else
            write(*,'(a)') 'binary parameters after core-collapse:'
            write(*,'(a)')
            write(*,'(a32, 2x, 1pes40.16e3)') 'period_after_cc =', after_cc_period
            write(*,'(a32, 2x, 1pes40.16e3)') 'excentricity_after_cc =', after_cc_eccentricity
         end if

      end subroutine prepare_kick

      end module bin2dco_utils
         
