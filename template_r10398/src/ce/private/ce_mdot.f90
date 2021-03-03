
      module ce_mdot

      use ce_def
      use const_def
      use star_def
      use star_lib
      use binary_def
      use binary_lib
      use crlibm_lib
      
      implicit none

      contains

      
      ! evaluate ce mdot
      ! NOTE: rlo_mdot is the variable that must be set here, and it must be in g/s
      real(dp) function eval_ce_mdot(b, ierr) result(mdot)
         type(binary_info), pointer :: b
         integer, intent(out) :: ierr
         real(dp) :: rl_rel_limit
         real(dp) :: x0, x1, y0, y1, m

         ierr = 0

         if (ce_duration <= years_to_max_mdot_rlof .and. &
             abs(b% mtransfer_rate) * secyer/Msun < max_mdot_rlof) then
            x0 = 0d0; x1 = years_to_max_mdot_rlof
            y0 = lg_mtransfer_rate_start_ce; y1 = safe_log10_cr(max_mdot_rlof)
            m = (y1 - y0) / (x1 - x0)
            mdot = - exp10_cr(m*(ce_duration-x0)+y0) * Msun/secyer
            return
         end if

         if (ce_type == ce_two_stars) then
            rl_rel_limit = tol_two_stars
         else
            rl_rel_limit = tol_xrb
         end if

         if (b% rl_relative_gap(b% d_i) > 0d0) then
            mdot = - max_mdot_rlof * Msun/secyer
         ! If R<RL apply another linear decreasing in log10(mdot)
         !else (b% rl_relative_gap(b% d_i) > rl_rel_limit) then
         else
            write(*,*) 'reducing mdot_rlof'
            x0 = 0d0; x1 = rl_rel_limit
            y0 = safe_log10_cr(max_mdot_rlof)
            y1 = lg_mtransfer_rate_start_ce -1d0
            m = (y1 - y0) / (x1 - x0)
            mdot = - exp10_cr(m*(b% rl_relative_gap(b% d_i)-x0)+y0) * Msun/secyer
         end if
      end function eval_ce_mdot

      end module ce_mdot

