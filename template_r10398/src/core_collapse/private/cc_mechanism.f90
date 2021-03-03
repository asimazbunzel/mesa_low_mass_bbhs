
      module cc_mechanism

      use cc_def
      use const_def
      use utils_lib
      use crlibm_lib
      use num_lib
      use star_lib
      
      implicit none

      public :: set_explosion_mechanism, do_explode_star

      contains

      subroutine set_explosion_mechanism(mechanism_name, ierr)
         character(len=*), intent(in) :: mechanism_name
         integer, intent(out) :: ierr

         include 'formats.inc'

         ! assign name to a number as defined in public/cc_def.f
         if (mechanism_name == 'rapid') then
            model_id = cc_rapid
         else if (mechanism_name == 'delayed') then
            model_id = cc_delayed
         else if (mechanism_name == 'startrack') then
            model_id = cc_startrack
         else if (mechanism_name == 'combine') then
            model_id = cc_combine
         else
            write(*,*) 'unknown sn model'
            ierr = -1
            return
         end if

         ierr = 0

      end subroutine set_explosion_mechanism

      subroutine do_explode_star(s, ierr)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr

         real(dp) :: he_core_mass, c_core_mass
         real(dp) :: baryon2grav, baryon2grav_factor, baryonic_mass
         real(dp) :: h_envelope, he_envelope
         real(dp) :: remnant_mass
         real(dp) :: proto_mass, fallback_mass

         real(dp) :: a1, b1, a2, b2

         real(dp) :: yy1, yy2, yy3

         include 'formats.inc'

         he_core_mass = s% he_core_mass
         c_core_mass = s% c_core_mass

         if (c_core_mass <= 0d0) then
            write(*,*) 'failed in do_explode_star. carbon core mass <= 0'
            ierr = -1
            return
         end if

         ! First, mass of proto-compact object is evaluated
         select case (model_id)
         case (cc_rapid)
            proto_mass = 1d0
         case (cc_delayed)
            if (c_core_mass < 3.5d0) then
               proto_mass = 1.2d0
            else if (c_core_mass < 6d0) then
               proto_mass = 1.3d0
            else if (c_core_mass < 11d0) then
               proto_mass = 1.4d0
            else if (c_core_mass > 11d0) then
               proto_mass = 1.6d0
            else
               stop 'failed to get range of proto_mass for delayed mechanism'
            end if
         case (cc_startrack)
            if (c_core_mass < 4.82d0) then
               proto_mass = 1.5d0
            else if (c_core_mass < 6.31d0) then
               proto_mass = 2.11d0
            else if (c_core_mass < 6.75d0) then
               proto_mass = 0.69d0 * c_core_mass - 2.26d0
            else if (c_core_mass > 6.75d0) then
               proto_mass = 0.37d0 * c_core_mass - 0.07d0
            else
               stop 'failed to get range of proto_mass for startrack mechanism'
            end if
         case (cc_combine)
            if (c_core_mass > 1.37d0 .and. c_core_mass < 1.435d0) then
               proto_mass = c_core_mass
            else if (c_core_mass > 1.435d0 .and. c_core_mass < 6.5d0) then
               ! Find if star is stripped of the He envelope
               h_envelope = s% star_mass - s% he_core_mass
               if (abs(h_envelope - 0d0) < 1d-6) then
                  he_envelope = s% star_mass - s% c_core_mass
                  if (abs(he_envelope - 0d0) < 1d-6) then 
                     proto_mass = - (1d0/0.618d0) + sqrt(pow2(1d0/0.618d0) + (1.06d0/0.084d0) * pow_cr(c_core_mass,0.454d0))
                  else
                     proto_mass = 0.23d0 * c_core_mass + 0.83d0
                  end if
               else
                  proto_mass = 0.23d0 * c_core_mass + 0.83d0
               end if
            else if (c_core_mass > 6.5d0) then
               h_envelope = s% star_mass - s% he_core_mass
               he_envelope = (s% star_mass - s% c_core_mass) - h_envelope
               proto_mass = c_core_mass + 0.8d0 * he_envelope
            else
               stop 'failed to get range of proto_mass for combine mechanism'
            end if
         case default
            write(*,*) 'bad name for sn_model'
            ierr = -1
            return
         end select

         ! Second, evaluate fallback for rapid, delayed & startrack
         if (model_id /= cc_combine) then
            select case (model_id)
            case (cc_rapid)
               if (c_core_mass < 2.5d0) then
                  fallback_mass = 0.2d0
               else if (c_core_mass < 6d0) then
                  fallback_mass = 0.286d0 * c_core_mass - 0.514d0
               else if (c_core_mass < 7d0) then
                  fallback_mass = 1d0 * (s% star_mass - proto_mass)
               else if (c_core_mass < 11d0) then
                  a1 = 0.25d0 - (1.275d0 / (s% star_mass - proto_mass))
                  b1 = -11d0 * a1 + 1d0
                  fallback_mass = (a1 * c_core_mass + b1) * (s% star_mass - proto_mass)
               else if (c_core_mass > 11d0) then
                  fallback_mass = 1.0d0 * (s% star_mass - proto_mass)
               else
                  stop 'failed to get range of fallback_mass for rapid mechanism'
               end if
            case (cc_delayed)
               if (c_core_mass < 2.5d0) then
                  fallback_mass = 0.2d0
               else if (c_core_mass < 3.5d0) then
                  fallback_mass = 0.5d0 * c_core_mass - 1.05d0
               else if (c_core_mass < 11d0) then
                  a2 = 0.133d0 - (0.093d0 / (s% star_mass - proto_mass))
                  b2 = -11d0 * a2 + 1d0
                  fallback_mass = (a2 * c_core_mass + b2) * (s% star_mass - proto_mass)
               else if (c_core_mass > 11d0) then
                  fallback_mass = 1.0d0 * (s% star_mass - proto_mass)
               else
                  stop 'failed to get range of fallback_mass for rapid mechanism'
               end if
            case (cc_startrack)
               if (c_core_mass < 5d0) then
                  fallback_mass = 0d0
               else if (c_core_mass < 7.6d0) then
                  fallback_mass = (0.378d0 * c_core_mass - 1.889d0) * (s% star_mass - proto_mass)
               else if (c_core_mass > 7.6d0) then
                  fallback_mass = 1d0 * (s% star_mass - proto_mass)
               else
                  stop 'failed to get range of fallback_mass for startrack mechanism'
               end if
            case default
               write(*,*) 'bad name for sn_model'
               ierr = -1
               return
            end select

         else
            fallback_mass = 0d0
         end if

         ! eval fallback fraction. Note that this is well defined in the StarTrack, rapid and
         ! delayed mechanism but it is rather difficult to find a simple way to write it for
         ! the combine mechanism. Thus we ONLY consider reasonable values those from any
         ! mechanism that is not combine
         if (model_id /= cc_combine) then
            fallback_fraction = fallback_mass / (s% star_mass - proto_mass)
         else
            fallback_fraction = -1d0
         end if


         ! Now that fallback_mass & proto_mass are set, calculate gravitational mass of compact
         ! remnant
         baryonic_mass = proto_mass + fallback_mass
         if (baryonic_mass < max_ns_mass) then
            if (model_id == cc_combine) then
               remnant_mass = 0.9d0 * baryonic_mass
            else
               ! use eq. 13 of Fryer et al. 2012 for NS mass remnant:
               ! if x = remnant_mass; x0 = baryonic_mass, then
               !        0.075 * x**2 + x - x0 = 0
               ! and solve for x
               yy1 = aux_grav_mass(baryonic_mass, 0d0)
               yy2 = aux_grav_mass(baryonic_mass, max_ns_mass/2d0)
               yy3 = aux_grav_mass(baryonic_mass, max_ns_mass)
               remnant_mass = find0_quadratic(0d0, yy1, max_ns_mass/2d0, yy2, max_ns_mass, yy3, ierr)

               if (ierr /= 0) then
                  write(*,*) 'bad root in find0_quadratic for remnant < 2.5 Msun'
                  return
               end if
            end if
         else
            baryon2grav_factor = 0.9d0
            if (model_id == cc_combine) baryon2grav_factor = 0.8d0
            remnant_mass = baryon2grav_factor * baryonic_mass
         end if

         ! Update core-collapse info
         M_baryonic = baryonic_mass
         M_remnant = remnant_mass
         M_ejected = s% star_mass - baryonic_mass
         M_fallback = fallback_mass

      end subroutine do_explode_star

      real(dp) function aux_grav_mass(bar_mass, grav_mass)
         real(dp), intent(in) :: bar_mass, grav_mass

         aux_grav_mass = 0.075d0 * pow2(grav_mass) + grav_mass - bar_mass

      end function aux_grav_mass 

      end module cc_mechanism
