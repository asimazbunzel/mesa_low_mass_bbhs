
      module cc_def

      use const_def

      implicit none

      ! integer numbers associated with different
      ! SN engine explosions
      integer, parameter :: cc_rapid = 1
      integer, parameter :: cc_delayed = 2
      integer, parameter :: cc_startrack = 3
      integer, parameter :: cc_combine = 4

      ! id of exploding star
      integer :: cc_id

      ! name of SN model
      character (len=strlen) :: model_name
      integer :: model_id

      ! filenames
      character(len=strlen) :: cc_data_directory
      character(len=strlen) :: filename_for_binary_data
      character(len=strlen) :: filename_for_star_data

      ! max NS mass
      real(dp) :: max_ns_mass

      ! pre-cc data
      integer :: model_number_pre_cc
      real(dp) :: age_pre_cc
      real(dp) :: he_core_mass_pre_cc, c_core_mass_pre_cc
      real(dp) :: he_core_radius_pre_cc, c_core_radius_pre_cc
      real(dp) :: M_pre_cc, R_pre_cc, Teff_pre_cc, L_pre_cc
      logical :: save_model_pre_cc

      ! fraction of mass that falls back
      real(dp) :: fallback_fraction

      ! post-cc data (in Msun)
      real(dp) :: M_baryonic, M_remnant
      real(dp) :: M_ejected, M_fallback

      ! id if binary system (0 if not, 1 else)
      integer :: cc_binary_id

      ! binary pre-cc data
      real(dp) :: progenitor_mass, companion_mass
      real(dp) :: angular_momentum_pre_cc, separation_pre_cc, &
         period_pre_cc, eccentricity_pre_cc
      real(dp) :: radius_pre_cc(2), rl_pre_cc(2)
      real(dp) :: mt_rate_pre_cc
      
      ! binary post-cc data
      logical :: continue_binary_evolution
      logical :: add_asymmetric_kick
      real(dp) :: angular_momentum_after_cc, separation_after_cc, &
         rl_after_cc(2), mt_rate_after_cc

      end module cc_def

