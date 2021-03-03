
      module ce_def

      use const_def
      
      implicit none

      ! binary id copy to internal id
      integer :: ce_id

      ! id for each binary component
      integer :: ce_donor_id, ce_accretor_id

      ! type of ce:
      !     - star + star (alpha_two_stars)
      !     - star + compact object (alpha_xrb)
      integer :: ce_type

      ! integer number associated to different ce types
      integer, parameter :: ce_two_stars = 1
      integer, parameter :: ce_xrb = 2

      ! number of events the same type of ce has happened
      integer :: ce_num_two_stars, ce_num_xrb
      integer :: ce_num_star_1, ce_num_star_2

      ! efficiency for converting orbital energy into kinetic energy of the envelope
      ! for different binary types
      real(dp) :: alpha_two_stars, alpha_xrb

      ! scaling factor of the Eddington mass-loss limit to trigger ce
      real(dp) :: edd_scaling_factor

      ! tolerance to determine end of ce, depending on the binary system
      real(dp) :: tol_two_stars, tol_xrb

      ! wait at least this many years before turning ce off
      real(dp) :: years_in_detachment

      ! years to reach the maximum mass-transfer rate during ce
      real(dp) :: years_to_max_mdot_rlof

      ! maximum mass-transfer rate during ce
      real(dp) :: max_mdot_rlof

      ! consider accretion energy in the ce formalism
      logical :: add_accretion_on_ce

      ! check if  we want to save profiles before and after ce
      logical :: save_profile_pre_ce, save_profile_after_ce
      ! and the same for models of star & donor (when possible)
      logical :: save_model_pre_ce, save_model_after_ce

      ! filenames for pre and after ce
      character(len=strlen) :: ce_data_directory
      character(len=strlen) :: filename_donor_profile_pre_ce, filename_donor_profile_after_ce
      character(len=strlen) :: filename_donor_model_pre_ce, filename_donor_model_after_ce
      character(len=strlen) :: filename_accretor_profile_pre_ce, filename_accretor_profile_after_ce
      character(len=strlen) :: filename_accretor_model_pre_ce, filename_accretor_model_after_ce

      ! max value for rl_relative_gap to consider ce merge
      real(dp) :: max_relative_gap
      ! max number of retries before considering ce merge
      integer :: max_number_retries_during_ce

      ! -----------------------------
      ! parameters during ce
      logical :: ce_on, ce_off
      logical :: ce_detach, ce_merge
      real(dp) :: ce_initial_age
      integer :: ce_initial_model_number
      real(dp) :: ce_duration, ce_years_in_detach
      real(dp) :: alpha_mt_start_ce, beta_mt_start_ce, delta_mt_start_ce, gamma_mt_start_ce
      real(dp) :: fj_start_ce
      real(dp) :: accretor_overflow_terminate_start_ce
      real(dp) :: lg_mtransfer_rate_start_ce
      real(dp) :: donor_mass_start_ce, donor_radius_start_ce
      real(dp) :: accretor_mass_start_ce, accretor_radius_start_ce

      real(dp) :: donor_grav_energy_prev_step, donor_int_energy_prev_step
      real(dp) :: donor_bind_energy_prev_step, donor_bind_energy
      real(dp) :: cumulative_removed_binding_energy
      real(dp) :: orbital_energy_prev_step

      end module ce_def

