! ----------------------------- GRASSLAND CARBON MODEL (DALEC.A3.H2.M2) -------------------------------------------------------------------------------
!
!                   autotrophic    heterotrophic     loss due to      --->    manure from          
!                   respiration    respiration       grazing/cutting         grazing livestock       
!                        ^            ^                   ^                        |       
!                        |            |                   |                        V      
!                                                                                   
! PHOTOSYNTHESIS -----> [0] -------> [0] --------------> [0] <------------------> [0]     
!                 GPP         NPP            NEE                   NBE                 
!
! ---------------------------------------------------------------------------------------------------------------------------------------------------------
! < POOLS - 6 > 
! [1] Labile/Stem [2] Foliar [3] Root [4] Litter [5] SOM [6] Soil Water Content (0-30cm)
!
! < FLUXES - 41 > 
! [1] GPP [2] Temprate [3] respiration auto [4] leaf production [5] labile production [6] root production [7] aboveground production
! [8] labile consumption → leaves [9] leaffall factor [10] leaf litter production [11] woodlitter production 
! [12] rootlitter production [13] respiration het litter [14] respiration het SOM 
! [15] litter → som [16] labile release factor(leaf growth) [17] C flux due to fire 
! [19] animal_manure_to_soil [20] animal_respiration [21] animal_methane 
! [25] HARVESTextracted_labile [26] HARVESTextracted_foliage [27] HARVESTextracted_roots [28] HARVESTlitter_labile
! [29] HARVESTlitter_foliage [30] HARVESTlitter_roots
! [31] GRAZINGextracted_labile [32] GRAZINGextracted_foliage [33] GRAZINGextracted_roots 
! [34] GRAZINGlitter_labile [35] GRAZINGlitter_foliage [36] GRAZINGlitter_roots
! [46] evaporation [47] transpiration [48] soil evaporation [49] wet canopy evaporation [50] runoff [51] underflow
! [* 18 & 22-24 & 37-45 are not used / empty] 
! [* all outputs come in per day timesteps]
 
! < PARAMETERS - 37 > 
! [1] Decomposition rate [1e-5, 0.01] [2] Fraction of GPP respired [~0.54] [3] Growing Season Index sensitivity for leaf growth  [1.0, 1.025] 
! [4] NPP belowground allocation parameter [0.01, 1.00] [5] GSI max leaf turnover [1e-5, 0.2] [6] TurnOver Rate roots [0.0001, 0.01]
! [7] TOR litter [0.0001, 0.01] [8] TOR SOM [1e-7, 0.001] [9] Temp factor Q10 [0.018,  0.08] [10] GSI max labile turnover  [1e-6, 0.2]
! [11] Photosynthetic N use efficiency [12] GSI min temperature threshold (K) [225, 330] [13] GSI max temperature threshold (K) [225, 330] 
! [14] GSI min photoperiod threshold (sec) [3600, 36000] [15] LCA (g.C.leaf_m-2) [20, 60]
! [16] C labile (initial)  [1, 1000] [17] C foliar (initial) [1, 1000] [18] C roots  (initial) [1, 1000]
! [19] C litter (initial) [1, 10000] [20] GSI max photoperiod threshold (sec) [3600, 64800] [21] GSI min VPD threshold (Pa) [1, 5500] 
! [22] GSI max VPD threshold (Pa) [1, 5500] [23] C SOM (initial) [24] GSI senstivity for leaf senescence [0.96, 1.00]
! [25] GSI - growing state (>1)  [0.50, 1.5] [26] GSI - initial GSI value [1.0, 2.0] [27] DM min lim for grazing (gC/m2) [500,1500]
! [28] DM min lim for cutting (gC/m2) [1500,3000] [29] leaf-vs-stem allocation factor [0.05, 0.75] [30] critical GPP for LAI growth [1e-10, 0.30]
! [31] DM demand of animal weight (fraction) [0.015,0.035] [32] Post-grazing labile loss (fraction) [0.01,0.1] [33] Post-cut labile loss (fraction) [0.5,0.9]
! [34] Minimum removed biomass to allow grazing [gC/m2/day] [35] Initial soil water - fraction of Field Capacity [0.5,1.0]
! [36] BUCKET - coarse root biomass needed to reach 50 % (i.e. gbio/m2 not gC/m2) [10,250] [37] BUCKET - maximum rooting depth [0.35,5.0]
!
! < INPUTS: MET & MGMT - 16 >
! [1] run day [2] min T (C) (3) max T (C) [4] Radiation (MJ.m-2) [5] CO2 (ppm) [6] DOY [7] Precipitation kg.H2O.m-2.s-1
! [8] vegetation management : [8.1] LAI removed (m2.m-2) -OR- [8.2] Livestock units per ha 
! [9*] Burnt area fraction [10] 21-day avg min T (K) [11] 21-day avg photoperiod (sec) [12] 21-day avg VPD (Pa) [13*] Forest mgmt after clearing
! [14] avg daily temperature (oC) [15] avg daily wind speed (m.s-1) [16] vapour pressure deficit (Pa)
! * Not used 
! ----------------------------------------------------------------------------------------------------------------------------------------------------------

module CARBON_MODEL_MOD

  implicit none

  ! make all private
  ! private

   ! explicit publics
  public :: CARBON_MODEL     &
           ,linear_model_gradient &
           ,layer_thickness  &
           ,wSWP_time        &
           ,rSWP_time        &
           ,cica_time        &
           ,root_depth_time        &
           ,gs_demand_supply_ratio &
           ,gs_total_canopy        &
           ,gb_total_canopy        &
           ,canopy_par_MJday_time  &
           ,sw_par_fraction  &
           ,snow_storage_time&
           ,soil_frac_clay   &
           ,soil_frac_sand   &
           ,nos_root_layers  &
           ,nos_soil_layers  &
           ,dim_1,dim_2      &
           ,nos_trees        &
           ,nos_inputs       &
           ,leftDaughter     &
           ,rightDaughter    &
           ,nodestatus       &
           ,xbestsplit       &
           ,nodepred         &
           ,bestvar 
     
 
  !!!!!!!!!
  ! Parameters
  !!!!!!!!!

  ! useful technical parameters
  double precision, parameter :: vsmall = tiny(0d0)*1d3 & ! *1d3 to add a little breathing room
                                ,vlarge = huge(0d0)

  integer, parameter :: nos_root_layers = 2 &
                       ,nos_soil_layers = 3

  double precision, parameter :: pi = 3.1415927d0,  &
                               pi_1 = 0.3183099d0,  & ! pi**(-1d0)
                             two_pi = 6.283185d0,   & ! pi*2d0
                         deg_to_rad = 0.01745329d0, & ! pi/180d0
                sin_dayl_deg_to_rad = 0.3979486d0,  & ! sin( 23.45d0 * deg_to_rad )
                              boltz = 5.670400d-8,  & ! Boltzmann constant (W.m-2.K-4)
                         emissivity = 0.96d0,       &
                        emiss_boltz = 5.443584d-08, & ! emissivity * boltz
                        ppfd_to_par = 4.6d0,        & ! Conversion of umolPAR -> J
                    sw_par_fraction = 0.5d0,        & ! fraction of short-wave radiation which is PAR
                             freeze = 273.15d0,     & !
                  gs_H2Ommol_CO2mol = 0.001646259d0,& ! The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
              gs_H2Ommol_CO2mol_day = 142.2368d0,   & ! The ratio of H20:CO2 diffusion for gs, including seconds per day correction
                         gb_H2O_CO2 = 1.37d0,       & ! The ratio of H20:CO2 diffusion for gb (Jones appendix 2)
            partial_molar_vol_water = 18.05d-6,     & ! partial molar volume of water, m3 mol-1 at 20C
                     mol_to_g_water = 18d0,         & ! molecular mass of water
                   mmol_to_kg_water = 1.8d-5,       & ! milli mole conversion to kg
                         umol_to_gC = 1.2d-5,       & ! conversion of umolC -> gC
                         gC_to_umol = 83333.33d0,   & ! conversion of gC -> umolC; umol_to_gC**(-1d0)
                               Rcon = 8.3144d0,     & ! Universal gas constant (J.K-1.mol-1)
                          vonkarman = 0.41d0,       & ! von Karman's constant
                        vonkarman_1 = 2.439024d0,   & ! 1 / von Karman's constant
                              cpair = 1004.6d0        ! Specific heat capacity of air; used in energy balance J.kg-1.K-1

  ! hydraulic parameters
  double precision, parameter :: &
                         tortuosity = 2.5d0,        & ! tortuosity
                             gplant = 4d0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                        root_resist = 25d0,         & ! Root resistivity (MPa s g mmol−1 H2O)
                        root_radius = 0.00029d0,    & ! root radius (m) Bonen et al 2014 = 0.00029
                                                      ! Williams et al 1996 = 0.0001
                      root_radius_1 = root_radius**(-1d0), &
                root_cross_sec_area = pi * root_radius**2, & ! root cross sectional area (m2)
                                                             ! = pi * root_radius * root_radius
                       root_density = 0.31d6,       & ! root density (g biomass m-3 root)
                                                      ! 0.5e6 Williams et al 1996
                                                      ! 0.31e6 Bonan et al 2014
            root_mass_length_coef_1 = (root_cross_sec_area * root_density)**(-1d0), &
                 const_sfc_pressure = 101325d0,     & ! (Pa)  Atmospheric surface pressure
                               head = 0.009807d0,   & ! head of pressure (MPa/m)
                             head_1 = 101.968d0       ! inverse head of pressure (m/MPa)

  ! structural parameters
  double precision, parameter :: &
                      canopy_height = 1d0,              & ! canopy height assumed to be 1 m for grass
                       tower_height = canopy_height + 2d0, & ! tower (observation) height assumed to be 2 m above canopy
                           min_wind = 0.2d0,            & ! minimum wind speed at canopy top
                          min_layer = 0.03d0,           & ! minimum thickness of the third rooting layer (m)
                        soil_roughl = 0.00085d0,        & ! soil roughness length (m), Meier et al., (2022), https://doi.org/10.5194/gmd-15-2365-2022
                       min_drythick = soil_roughl*10d0, & ! minimum dry thickness depth (m) 0.01 WRF-SPA 
                     top_soil_depth = 0.30d0,           & ! thickness of the top soil layer (m)
                           min_root = 5d0,              & ! minimum root biomass (gBiomass.m-2)
                            min_lai = 0.01d0,           & ! minimum LAI assumed for aerodynamic conductance calculations (m2/m2)
                        min_storage = 0.1d0               ! minimum canopy water (surface) storage (mm)

  ! timing parameters
  double precision, parameter :: &
                   seconds_per_hour = 3600d0,       & ! Number of seconds per hour
                    seconds_per_day = 86400d0,      & ! Number of seconds per day
                  seconds_per_day_1 = 1.157407d-05    ! Inverse of seconds per day

  ! ACM-GPP-ET parameters
  double precision, parameter :: &
                       Vc_minT = -6.991d0,     & ! Temperature at which all photosynthetic activity is shutdown
                       Vc_coef = 0.1408d0,     & ! Temperature above Vc_minT that 50% limitation of cold shutdown occurs  
                   ! Assumption that photosythesis will be limited by Jmax temperature response
                   pn_max_temp = 57.05d0,      & ! Maximum daily max temperature for photosynthesis (oC)
                   pn_min_temp = -1d6,         & ! Minimum daily max temperature for photosynthesis (oC)
                   pn_opt_temp = 30d0,         & ! Optimum daily max temperature for photosynthesis (oC)
                   pn_kurtosis = 0.172d0,      & ! Kurtosis of Jmax temperature response
                ko_half_sat_25C = 157.46892d0,  & ! photorespiration O2 half sat(mmolO2/mol), achieved at 25oC
           ko_half_sat_gradient = 14.93643d0,   & ! photorespiration O2 half sat gradient
                kc_half_sat_25C = 319.58548d0,  & ! carboxylation CO2 half sat (umolCO2/mol), achieved at 25oC
           kc_half_sat_gradient = 24.72297d0,   & ! carboxylation CO2 half sat gradient
                co2comp_sat_25C = 36.839214d0,  & ! carboxylation CO2 compensation point(umolCO2/mol), saturation
               co2comp_gradient = 9.734371d0,   & ! carboxylation CO2 comp point, achieved at oC
                                                  ! Each of these are temperature sensitivty
                            e0 = 3.2d+00,       & ! Quantum yield (gC/MJ/m2/day PAR), SPA apparent yield
                minlwp_default =-1.808224d+00,  & ! minimum leaf water potential (MPa). NOTE: actual SPA = -2 MPa
      soil_iso_to_net_coef_LAI =-2.717467d+00,  & ! Coefficient relating soil isothermal net radiation to net.
       soil_iso_to_net_coef_SW =-3.500964d-02,  & ! Coefficient relating soil isothermal net radiation to net.
         soil_iso_to_net_const = 3.455772d+00,  & ! Constant relating soil isothermal net radiation to net
     canopy_iso_to_net_coef_SW = 1.480105d-02,  & ! Coefficient relating SW to the adjustment between isothermal and net LW
       canopy_iso_to_net_const = 3.753067d-03,  & ! Constant relating canopy isothermal net radiation to net
    canopy_iso_to_net_coef_LAI = 2.455582d+00,  & ! Coefficient relating LAI to the adjustment between isothermal and net LW
                          iWUE = 4.6875d-04       ! Intrinsic water use efficiency (umolC/mmolH2O-1/m2leaf/s-1)

  double precision :: minlwp = minlwp_default

  ! Photosynthetic Metrics
  double precision, allocatable, dimension(:) :: gs_demand_supply_ratio, & ! actual:potential stomatal conductance
                                                        gs_total_canopy, & ! stomatal conductance (mmolH2O/m2ground/s)
                                                        gb_total_canopy, & ! boundary conductance (mmolH2O/m2ground/s)
                                                  canopy_par_MJday_time    ! Absorbed PAR by canopy (MJ/m2ground/day)

  ! arrays for the emulator, just so we load them once and that is it cos they be
  ! massive
  integer ::    dim_1, & ! dimension 1 of response surface
                dim_2, & ! dimension 2 of response surface
            nos_trees, & ! number of trees in randomForest
           nos_inputs    ! number of driver inputs

  double precision, allocatable, dimension(:,:) ::     leftDaughter, & ! left daughter for forest
                                                      rightDaughter, & ! right daughter for forets
                                                         nodestatus, & ! nodestatus for forests 
                                                         xbestsplit, & ! for forest
                                                           nodepred, & ! prediction value for each tree
                                                            bestvar    ! for randomForests

  ! hydraulic model variables
  integer :: water_retention_pass, soil_layer
  double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand ! clay and soil fractions of soil
  double precision, dimension(nos_root_layers) :: uptake_fraction, & ! fraction of water uptake from each root layer
                                                           demand, & ! maximum potential canopy hydraulic demand
                                            water_flux_mmolH2Om2s    ! potential transpiration flux (mmolH2O.m-2.s-1)
  double precision, dimension(nos_soil_layers+1) :: SWP, & ! soil water potential (MPa)
                                            SWP_initial, &
                                      soil_conductivity, & ! soil conductivity
                                            waterchange, & ! net water change by specific soil layers (m)
                                         field_capacity, & ! soil field capacity (m3.m-3)
                                 field_capacity_initial, &
                                         soil_waterfrac, & ! soil water content (m3.m-3)
                                 soil_waterfrac_initial, &
                                               porosity, & ! soil layer porosity, (fraction)
                                       porosity_initial, &
                                        layer_thickness, & ! thickness of soil layers (m)
                        cond1, cond2, cond3, potA, potB    ! Saxton equation values

  double precision :: root_reach, root_biomass, &
                             fine_root_biomass, & ! root depth, coarse+fine, and fine root biomass
                              total_water_flux, & ! potential transpiration flux (kgH2O.m-2.day-1)
                                      drythick, & ! estimate of the thickness of the dry layer at soil surface (m)
                                          wSWP, & ! soil water potential weighted by canopy supply (MPa)
                                          rSWP, & ! soil water potential weighted by root presence (MPa)
                                     max_depth, & ! maximum possible root depth (m)
                                        root_k, & ! biomass to reach half max_depth
                                        runoff, & ! runoff (kgH2O.m-2.day-1)
                                     underflow, & ! drainage from the bottom of soil column (kgH2O.m-2.day-1)
                                previous_depth, & ! depth of bottom of soil profile
                                   canopy_wind, & ! wind speed (m.s-1) at canopy top
                                         ustar, & ! friction velocity (m.s-1)
                                      ustar_Uh, &
                                air_density_kg, & ! air density kg/m3
                                ET_demand_coef, & ! air_density_kg * vpd_kPa * cpair
                                        roughl, & ! roughness length (m)
                                  displacement, & ! zero plane displacement (m)
                                    max_supply, & ! maximum water supply (mmolH2O/m2/day)
                                         meant, & ! mean air temperature (oC)
                                         leafT, & ! canopy day time temperature temperature (oC)
                            canopy_swrad_MJday, & ! canopy_absorbed shortwave radiation (MJ.m-2.day-1)
                              canopy_par_MJday, & ! canopy_absorbed PAR radiation (MJ.m-2.day-1)
                              soil_swrad_MJday, & ! soil absorbed shortwave radiation (MJ.m-2.day-1)
                              canopy_lwrad_Wm2, & ! canopy absorbed longwave radiation (W.m-2)
                                soil_lwrad_Wm2, & ! soil absorbed longwave radiation (W.m-2)
                                 sky_lwrad_Wm2, & ! sky absorbed longwave radiation (W.m-2)
                          stomatal_conductance, & ! canopy scale stomatal conductance (mmolH2O.m-2.d-1)
                         potential_conductance, & ! potential stomatal conductance (mmolH2O.m-2ground.s-1)
                           minimum_conductance, & ! potential stomatal conductance (mmolH2O.m-2ground.s-1)
                       aerodynamic_conductance, & ! aerodynamic conductance at canopy top (m.s-1)
                              soil_conductance, & ! soil surface conductance (m.s-1)
                             convert_ms1_mol_1, & ! Conversion ratio for m.s-1 -> mol.m-2.s-1
                            convert_ms1_mmol_1, & ! Conversion ratio for m/s -> mmol/m2/s
                           air_vapour_pressure, & ! Vapour pressure of the air (kPa)
                                        lambda, & ! latent heat of vapourisation (J.kg-1)
                                         psych, & ! psychrometric constant (kPa K-1)
                                         slope, & ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
                        water_vapour_diffusion, & ! Water vapour diffusion coefficient in (m2/s)
                           kinematic_viscosity, & ! kinematic viscosity (m2.s-1)
                                  snow_storage, & ! snow storage (kgH2O/m2)
                                canopy_storage, & ! water storage on canopy (kgH2O.m-2)
                          intercepted_rainfall    ! intercepted rainfall rate equivalent (kgH2O.m-2.s-1)

  ! Module level variables for the Sellers (1985) 2-stream radiative transfer scheme approximation
  integer, parameter :: no_wavelength = 2 ! Number of wavelenths (order NIR, PAR)
  double precision, parameter :: & !Vc = 0.75d0   & ! Clumping factor / vegetation cover (1 = uniform, 0 totally clumped, mean = 0.75)
                                                 ! He et al., (2012) http://dx.doi.org/10.1016/j.rse.2011.12.008
               soil_nir_reflectance = 0.023d0, & ! Soil reflectance to near infrared radiation
               soil_par_reflectance = 0.033d0, & ! Soil reflectance to photosynthetically active radiation
             canopy_nir_reflectance = 0.43d0,  & ! Canopy NIR reflectance
             canopy_par_reflectance = 0.16d0,  & ! Canopty PAR reflectance
           canopy_nir_transmittance = 0.26d0,  & ! Canopy NIR reflectance
           canopy_par_transmittance = 0.16d0,  & ! Canopty PAR reflectance
                    newsnow_nir_abs = 0.27d0,  & ! NIR absorption fraction
                    newsnow_par_abs = 0.05d0,  & ! PAR absorption fraction
         leaf_distribution_deviance = 0.01d0     ! Deviation from spherical, min absolute value (0.01) required for numerical security.
                                                 ! Leaf angle distribution, quantified as the deviation from a spherical distribution.
                                                 ! The default assumption in many models, including SPA, is that leaves have a spherical distribution (=0).
                                                 ! However, =-1 would indicate vertical leaves, while =+1 are horizontal leaves.
                                                 ! See note book or references given above for the complete integral equation
  double precision ::        leaf_angle, & ! Mean leaf angle deviation from the horizontal (radians)
                              cos2theta, & ! Analytical correction for leaf angle (radians) on light scattering within the canopy
                                 Vc, Vg, & ! Define the vegetated and covered soil (i.e. by litter) fractions
                                mu_obar, & ! The average inverse diffuse optical depth per unit leaf area.
                                 O1,O2    ! Empirical coefficients related to the leaf angle distribution
  double precision, dimension(no_wavelength) :: &
                     canopy_reflectance, & ! = (/canopy_nir_reflectance,canopy_par_reflectance/), & !
                   canopy_transmittance, & ! = (/canopy_nir_transmittance,canopy_par_transmittance/), & !
                       soil_reflectance, & ! = (/soil_nir_reflectance,soil_par_reflectance/), & !
                      canopy_scattering, & ! Canopy scattering of incident light, varied by wavelength
                                     bb, & ! Downward scatting of diffuse radiation
                                     cc, & ! Upward scattering as diffuse radiation, a function of canopy_transmittance, canopy_reflectance and leaf angle.
                                     hh, & ! Extinction coefficient for diffuse radiation
                                   beta, & ! The upward scattering fraction / coefficient for diffuse radiation.
                                  beta0    ! The upward scattering fraction / coefficient for direct radiation.

  ! Module level variables for ACM_GPP_ET parameters
  double precision ::   delta_gs, & ! day length corrected gs increment mmolH2O/m2/day
                            ceff, & ! Maximum rate of carboxylation (umolC/m2/s), Vcmax_ref = avN*NUE
!                             avN, & ! average foliar N (gN/m2)
                       iWUE_step, & ! Intrinsic water use efficiency for that day (gC/m2leaf/dayl/mmolH2Ogs)
!                             NUE, & ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
!                                    ! ,unlimited by CO2, light and photoperiod (umolC/gN/m2leaf)
metabolic_limited_photosynthesis, & ! temperature, leaf area and foliar N limiterd photosynthesis (gC/m2/day)
    light_limited_photosynthesis, & ! light limited photosynthesis (gC/m2/day)
                              ci, & ! Internal CO2 concentration (ppm)
                        rb_mol_1, & ! Canopy boundary layer resistance (day/m2/molCO2)
                     o2_half_sat, & ! O2 at which photorespiration is 50 % of maximum
                    co2_half_sat, & ! CO2 at which photosynthesis is 50 % of maximum (ppm)
                  co2_comp_point    ! CO2 at which photosynthesis > 0 (ppm)

  ! Module level variables for step specific met drivers
  double precision :: mint, & ! minimum temperature (oC)
                      maxt, & ! maximum temperature (oC)
        airt_zero_fraction, & ! fraction of air temperature above freezing
                     swrad, & ! incoming short wave radiation (MJ/m2/day)
                       co2, & ! CO2 (ppm)
                       doy, & ! Day of year
                  rainfall, & ! rainfall (kgH2O/m2/s)
                  snowfall, &
                 snow_melt, & ! snow melt (kgH2O/m2/s)
                  wind_spd, & ! wind speed (m/s)
                   vpd_kPa, & ! Vapour pressure deficit (kPa)
 leaf_canopy_light_scaling, & ! approximate scaling factor from leaf to canopy for gpp and gs
  leaf_canopy_wind_scaling, & ! approximate scaling factor from leaf to canopy for gb
                       lai    ! leaf area index (m2/m2)

  ! Module level variables for step specific timing and location information
  integer :: steps_per_year, gsi_lag_remembered, two_week_lag, four_week_lag     
  double precision ::       seconds_per_step, & !
                               days_per_step, & !
                             days_per_step_1, & !
                          mean_days_per_step, &
                                dayl_seconds, & ! day length in seconds
                              dayl_seconds_1, &
                         dayl_hours_fraction, &
                                  dayl_hours, & ! day length in hours
                                    latitude, & ! latitude ()-90/90)
                            latitude_radians, & ! latitude in radians
                        sin_latitude_radians, & ! sin(latitude_radians)
                        cos_latitude_radians, & ! cos(latitude_radians)
                          sunset_solar_angle, & ! Solar angle at sunset hour
                                 declination, & ! Solar declination, function of day of year
                   cosine_solar_zenith_angle    ! Cosine zenith angle of the timestep

  double precision, dimension(:), allocatable :: deltat_1, & ! inverse of decimal days
                                             tmp_x, tmp_m, &
                                  airt_zero_fraction_time, &
                                          daylength_hours, &
                                        daylength_seconds, &
                                      daylength_seconds_1, &
                                            rainfall_time, &
                                                cica_time, & ! Internal vs ambient CO2 concentrations
                                          root_depth_time, &
                                        snow_storage_time, &
                                                rSWP_time, & ! Soil water potential weighted by access water
                                                wSWP_time    ! Soil water potential weighted by supply of water

  ! GSI and management variables
  double precision :: fol_turn_crit,lab_turn_crit &
                     ,gsi_history(22),just_grown  &
                     ,foliage_frac_res   & ! Management related variables
                     ,labile_frac_res    &
                     ,roots_frac_res     &
                     ,roots_frac_death                                              


  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine CARBON_MODEL(start,finish,deltat,lat,met,pars, &
                          nodays,nopars,nomet,nopools,nofluxes, &
                          lai_out,GPP,NEE,FLUXES,POOLS)

    ! The Data Assimilation Linked Ecosystem Carbon - Combined Deciduous
    ! Evergreen Analytical - ACMv2 - BUCKET (DALEC_CDEA_ACM2_BUCKET) model.
    ! The subroutine calls the Aggregated Canopy Model version 2 to simulate GPP and partitions
    ! between various ecosystem carbon pools. These pools are subject
    ! to turnovers / decompostion resulting in ecosystem phenology and fluxes of CO2
    ! ACMv2 simulates coupled photosynthesis-transpiration (via stomata), soil and intercepted canopy
    ! evaporation and soil water balance (4 layers).
    ! This version include grassland management specific functions and lacks a wood / structural C pool

    implicit none

    ! declare input variables
    integer, intent(in) :: start    &
                          ,finish   & 
                          ,nodays   & ! number of days in simulation
                          ,nopars   & ! number of paremeters in vector
                          ,nomet    & ! number of meteorological fields
                          ,nopools  & ! number of model pools
                          ,nofluxes   ! number of model fluxes

    double precision, intent(in) :: deltat(nodays)    & ! time step in decimal days
                                   ,lat               & ! site latitude (degrees)
                                   ,met(nomet,nodays) & ! met drivers
                                   ,pars(nopars)        ! number of parameters

    double precision, dimension(nodays), intent(out) :: lai_out & ! leaf area index
                                                             ,GPP & ! Gross primary productivity
                                                             ,NEE   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1),nopools), intent(out) :: POOLS ! vector of ecosystem pools

    double precision, dimension(nodays,nofluxes), intent(out) :: FLUXES ! vector of ecosystem fluxes    

    !f2py intent(in) :: start, finish, deltat, lat, met, pars, nodays, nopars, nomet, nopools, nofluxes   

    !f2py intent(out) :: lai_out, GPP, NEE, POOLS, FLUXES

    ! declare local variables
    double precision ::      infi &
                   ,transpiration & ! kgH2O/m2/day
                 ,soilevaporation & ! kgH2O/m2/day
                  ,wetcanopy_evap & ! kgH2O/m2/day
                 ,snowsublimation   ! kgH2O/m2/day

    ! GSI or dynamics allocation variables
    double precision :: gsi_lai_reduction &
                       ,tot_abg_exp       &
                       ,fol_frac,lab_frac &
                       ,f_root,NPP                

    integer :: f,m,n,test,gsi_lag

    ! Set some initial states for the io variables
    infi = 0d0 ; FLUXES = 0d0 ; POOLS = 0d0
    ! Reset hydrology variables
    intercepted_rainfall = 0d0 ; canopy_storage = 0d0 ; snow_storage = 0d0
    transpiration = 0d0 ; soilevaporation = 0d0 ; wetcanopy_evap = 0d0 ; snowsublimation = 0d0
    ! Reset radiation variabes
    canopy_swrad_MJday = 0d0 ; canopy_par_MJday = 0d0 ; soil_swrad_MJday = 0d0 
    canopy_lwrad_Wm2 = 0d0 ; soil_lwrad_Wm2 = 0d0 ; sky_lwrad_Wm2 = 0d0
    ! Reset conductance variables
    soil_conductance = 0d0

    ! post-removal residues and root death | 0:none 1:all
    foliage_frac_res  = 0.05d0  ! fraction of removed foliage that goes to litter
    labile_frac_res   = 0.05d0  ! fraction of removed labile that goes to litter
    roots_frac_res    = 1d0     ! fraction of roots that die which go to litter 
    roots_frac_death  = 0.01d0  ! fraction of roots that dies
    
    ! How many steps in 2 weeks
    two_week_lag = ceiling(14d0/deltat(1))
    ! How many steps in 4 weeks
    four_week_lag = ceiling(28d0/deltat(1))

    ! Generate some generic location specific variables for radiation balance
    !call calculate_radiation_commons(lat,pars(33:38))
    call calculate_radiation_commons(lat)

    ! load ACM-GPP-ET parameters
    ceff = pars(11) ! Canopy efficiency (umolC/m2/s)
                    ! This is in the full model the product of Nitrogen use efficiency (umolC/gN/m2leaf)
                    ! and average foliar nitrogen gN/m2leaf
    ! Rooting parameters
    root_k = pars(36) ; max_depth = pars(37)

    ! assigning initial conditions
    POOLS(1,1) = pars(16)
    POOLS(1,2) = pars(17)
    POOLS(1,3) = pars(18)
    POOLS(1,4) = pars(19)
    POOLS(1,5) = pars(23)

    ! Some time consuming variables we only want to set once
    if (.not.allocated(deltat_1)) then
        ! allocate variables dimension which are fixed per site only the once
        allocate(deltat_1(nodays),wSWP_time(nodays),rSWP_time(nodays),gs_demand_supply_ratio(nodays), &
                 gs_total_canopy(nodays),gb_total_canopy(nodays),canopy_par_MJday_time(nodays), &
                 daylength_hours(nodays),daylength_seconds(nodays),daylength_seconds_1(nodays), &
                 rainfall_time(nodays),cica_time(nodays),root_depth_time(nodays),snow_storage_time(nodays))

        !
        ! Timing variables which are needed first
        !

        deltat_1 = deltat**(-1d0)

        !
        ! Iteration independent variables using functions and thus need to be in a loop
        !

        ! first those linked to the time period of the analysis
        do n = 1, nodays
           ! check positive values only for rainfall input
           rainfall_time(n) = max(0d0,met(7,n))
           ! Calculate declination for the day of year
           declination = calculate_declination((met(6,n)-(deltat(n)*0.5d0)))
           ! calculate daylength in hours and seconds
           call calculate_daylength
           daylength_hours(n) = dayl_hours ; daylength_seconds(n) = dayl_seconds
        end do

        ! calculate inverse for each time step in seconds
        daylength_seconds_1 = daylength_seconds ** (-1d0)
        ! fraction of temperture period above freezing
        airt_zero_fraction_time = (met(3,:)-0d0) / (met(3,:)-met(2,:))

        ! number of time steps per year
        steps_per_year = nint(dble(nodays)/(sum(deltat)*0.002737851d0))
        ! mean days per step
        mean_days_per_step = sum(deltat) / dble(nodays)

        !
        ! Initialise the water model
        !

        ! zero variables not done elsewhere
        total_water_flux = 0d0 ; water_flux_mmolH2Om2s = 0d0
        ! initialise some time invarient parameters
        call saxton_parameters(soil_frac_clay,soil_frac_sand)
        call initialise_soils(soil_frac_clay,soil_frac_sand)
        call update_soil_initial_conditions(pars(35))
        ! save the initial conditions for later
        soil_waterfrac_initial = soil_waterfrac
        SWP_initial = SWP
        field_capacity_initial = field_capacity
        porosity_initial = porosity

    else ! deltat_1 allocated?

        !
        ! Load initial soil water conditions from memory
        !

        total_water_flux = 0d0 ; water_flux_mmolH2Om2s = 0d0
        soil_waterfrac = soil_waterfrac_initial
        SWP = SWP_initial
        field_capacity = field_capacity_initial
        porosity = porosity_initial

        ! input initial soil water fraction then
        ! update SWP and soil conductivity accordingly
        call update_soil_initial_conditions(pars(35))

    endif ! deltat_1 allocated

    ! calculate some values once as these are invarient between DALEC runs
    if (.not.allocated(tmp_x)) then
        ! 21 days is the maximum potential so we will fill the maximum potential
        ! + 1 for safety
        allocate(tmp_x(22),tmp_m(nodays))
        do f = 1, 22
           tmp_x(f) = f
        end do
        do n = 1, nodays
           ! calculate the gradient / trend of GSI
           if (sum(deltat(1:n)) < 21) then
               tmp_m(n) = n-1
           else
               ! else we will try and work out the gradient to see what is happening
               ! to the system over all. The default assumption will be to consider
               ! the averaging period of GSI model (i.e. 21 days). If this is not
               ! possible either the time step of the system is used (if step greater
               ! than 21 days) or all available steps (if n < 21).
               m = 0 ; test = 0
               do while (test < 21)
                  m = m+1 ; test = floor(sum(deltat((n-m):n)))
                  if (m > (n-1)) then 
                      test = 21 
                  endif
               end do
               tmp_m(n) = m
           endif ! for calculating gradient
        end do ! calc daily values once
        ! allocate GSI history dimension
        gsi_lag_remembered = max(2,maxval(nint(tmp_m)))
    end if ! .not.allocated(tmp_x)
    ! assign our starting value
    gsi_history = pars(26)-1d0
    just_grown = pars(25)

    ! assign climate sensitivities
    gsi_lag = gsi_lag_remembered ! added to prevent loss from memory
    fol_turn_crit = pars(24)-1d0
    lab_turn_crit = pars(3)-1d0

    ! 
    ! Begin looping through each time step
    ! 

    ! load some needed module level values
    lai = POOLS(1,2)/pars(15)
    mint = met(2,1)  ! minimum temperature (oC)
    maxt = met(3,1)  ! maximum temperature (oC)
    swrad = met(4,1) ! incoming short wave radiation (MJ/m2/day)
    co2 = met(5,1)   ! CO2 (ppm)
    doy = met(6,1)   ! Day of year
    rainfall = rainfall_time(1) ! rainfall (kgH2O/m2/s)
    wind_spd = met(15,1) ! wind speed (m/s)
    vpd_kPa = met(16,1)*1d-3 ! vapour pressure deficit (Pa->kPa)
    leafT = (maxt*0.75d0) + (mint*0.25d0)   ! initial day time canopy temperature (oC)
    seconds_per_step = deltat(1) * seconds_per_day
    days_per_step =  deltat(1)
    days_per_step_1 =  deltat_1(1)

    ! calculate some temperature dependent meteorologial properties
    call meteorological_constants(leafT,leafT+freeze,vpd_kPa)

    ! Initialise root reach based on initial coarse root biomass
    fine_root_biomass = max(min_root,POOLS(1,3)*2d0)
    root_biomass = fine_root_biomass 
    ! calculate soil depth to which roots reach - needed here to set up
    ! layer_thickness correctly!
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! Determine initial soil layer thickness
    layer_thickness(1) = top_soil_depth ; layer_thickness(2) = max(min_layer,root_reach-top_soil_depth)
    layer_thickness(3) = max_depth - sum(layer_thickness(1:2))
    layer_thickness(4) = top_soil_depth
    previous_depth = sum(layer_thickness(1:2))
    ! Needed to initialise soils
    call calculate_Rtot
    call calculate_update_soil_water(transpiration,soilevaporation,snowsublimation,&
                                     0d0,FLUXES(1,46)) ! assume no evap or rainfall
    ! Reset variable used to track ratio of water supply used to meet demand
    gs_demand_supply_ratio = 0d0

    ! Store soil water content of the surface zone (mm)
    POOLS(1,6) = 1d3 * soil_waterfrac(1) * layer_thickness(1)

    do n = start, finish  

       !!!!!!!!!!
       ! assign drivers and update some prognostic variables
       !!!!!!!!!!

       ! Incoming drivers
       mint = met(2,n)  ! minimum temperature (oC)
       maxt = met(3,n)  ! maximum temperature (oC)
       swrad = met(4,n) ! incoming short wave radiation (MJ/m2/day)
       co2 = met(5,n)   ! CO2 (ppm)
       doy = met(6,n)   ! Day of year
       rainfall = rainfall_time(n)
       meant = (mint + maxt) * 0.5d0 ! mean air temperature (oC)
       leafT = (meant + maxt) * 0.5d0 ! estimate mean daytime air temperature (oC)
       wind_spd = met(15,n) ! wind speed (m/s)
       vpd_kPa = met(16,n)*1d-3  ! Vapour pressure deficit (Pa -> kPa)

       ! calculate LAI value
       lai_out(n) = POOLS(n,2)/pars(15)
       lai = lai_out(n)

       ! extract timing related values
       dayl_hours = daylength_hours(n)
       dayl_hours_fraction = dayl_hours * 0.04166667d0 ! 1/24 = 0.04166667
       dayl_seconds = daylength_seconds(n) ; dayl_seconds_1 = daylength_seconds_1(n)
       seconds_per_step = seconds_per_day * deltat(n)
       days_per_step = deltat(n) ; days_per_step_1 = deltat_1(n)

       !!!!!!!!!!
       ! Adjust snow balance balance based on temperture
       !!!!!!!!!!

       ! snowing or not...?
       if (((mint + maxt) * 0.5d0) > 0d0) then
           ! on average above freezing so no snow
           snowfall = 0d0
       else
           ! on average below freezing, so some snow based on proportion of temperture
           ! below freezing
           snowfall = rainfall * (1d0 - airt_zero_fraction) ; rainfall = rainfall - snowfall
           ! Add rainfall to the snowpack and clear rainfall variable
           snow_storage = snow_storage + (snowfall*seconds_per_step)
       end if

       ! melting or not...?
       if (mint < 0d0 .and. maxt > 0d0) then
           ! Also melt some of the snow based on airt_zero_fraction
           ! default assumption is that snow is melting at 10 % per day hour above freezing
           snow_melt = min(snow_storage, airt_zero_fraction * snow_storage * 0.1d0 * deltat(n))
           snow_storage = snow_storage - snow_melt
           ! adjust to rate for later addition to rainfall
           snow_melt = snow_melt / seconds_per_step
       elseif (maxt < 0d0) then
           snow_melt = 0d0
           ! Add rainfall to the snowpack and clear rainfall variable
           snow_storage = snow_storage + (snowfall*seconds_per_step)
       else if (mint > 0d0 .and. snow_storage > 0d0) then
           ! otherwise we assume snow is melting at 10 % per day above hour
           snow_melt = min(snow_storage, snow_storage * 0.1d0 * deltat(n))
           snow_storage = snow_storage - snow_melt
           ! adjust to rate for later addition to rainfall
           snow_melt = snow_melt / seconds_per_step
       else
           snow_melt = 0d0
       end if
       snow_storage_time(n) = snow_storage

       !!!!!!!!!!
       ! Calculate surface exchange coefficients
       !!!!!!!!!!

       ! calculate some temperature dependent meteorologial properties
       call meteorological_constants(leafT,leafT+freeze,vpd_kPa)
       ! pass variables from memory objects
       convert_ms1_mmol_1 = convert_ms1_mol_1 * 1d3
       ! calculate aerodynamic using consistent approach with SPA
       call calculate_aerodynamic_conductance
       ! Canopy scale aerodynamic conductance (mmolH2O/m2ground/s)
       gb_total_canopy(n) = aerodynamic_conductance * convert_ms1_mmol_1 * &
                            leaf_canopy_wind_scaling

       !!!!!!!!!!
       ! Determine net shortwave and isothermal longwave energy balance
       !!!!!!!!!!

       call calculate_radiation_balance
       canopy_par_MJday_time(n) = canopy_par_MJday

       !!!!!!!!!!
       ! Calculate physically constrained evaporation and
       ! soil water potential and total hydraulic resistance
       !!!!!!!!!!

       ! Canopy intercepted rainfall evaporation (kgH2O/m2/day)
       if (lai > 0d0) then ! is this conditional needed?
           call calculate_wetcanopy_evaporation(wetcanopy_evap,canopy_storage)
       else
           ! reset pools
           intercepted_rainfall = 0d0 ; canopy_storage = 0d0 ; wetcanopy_evap = 0d0
       endif

       ! calculate the minimum soil & root hydraulic resistance based on total
       ! fine root mass ! *2*2 => *RS*C->Bio
       fine_root_biomass = max(min_root,POOLS(n,3)*2d0)
       root_biomass = fine_root_biomass 
       call calculate_Rtot
       ! Pass wSWP to output variable
       wSWP_time(n) = wSWP ; rSWP_time(n) = rSWP ; root_depth_time(n) = root_reach

       ! calculate radiation absorption and estimate stomatal conductance
       call calculate_stomatal_conductance
       ! Estimate stomatal conductance relative to its minimum / maximum, i.e. how
       ! close are we to maxing out supply (note 0.01 taken from min_gs)
       gs_demand_supply_ratio(n) = (stomatal_conductance  - minimum_conductance) &
                                 / (potential_conductance - minimum_conductance)
       ! Store the canopy level stomatal conductance (mmolH2O/m2ground/s)
       gs_total_canopy(n) = stomatal_conductance

       ! Note that soil mass balance will be calculated after phenology
       ! adjustments

       ! Reset output variable
       if (stomatal_conductance > vsmall) then
           ! Gross primary productivity (gC/m2/day)
           ! Assumes acm_gpp_stage_1 ran as part of stomatal conductance calculation
           FLUXES(n,1) = acm_gpp_stage_2(stomatal_conductance) * umol_to_gC * dayl_seconds
           cica_time(n) = ci / co2
           ! Canopy transpiration (kgH2O/m2/day)
           call calculate_transpiration(transpiration)
           ! restrict transpiration to positive only
           transpiration = max(0d0,transpiration)
       else
           ! assume zero fluxes
           FLUXES(n,1) = 0d0 ; transpiration = 0d0 ; cica_time(n) = 0d0
       endif

      ! temprate (i.e. T modified rate of metabolic activity))
      FLUXES(n,2) = exp(pars(9)*0.5d0*(met(3,n)+met(2,n)))

      ! Allocate to autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = FLUXES(n,1) * pars(2)

      !! Determine direct allocation of GPP to plant tissues   
      ! Dynamic allocation to roots vs aboveground biomass after Reyes.et.al.2017 (10.1002/2017MS001022)
      ! min/max allocation to roots as fraction of NPP
      f_root = 1d0 - exp(-1d0*pars(4)*lai_out(n))
      if (f_root < 0.1d0) f_root = 0.1d0
      if (f_root > 0.7d0) f_root = 0.7d0      
      ! allocation to roots 
      FLUXES(n,6) = (FLUXES(n,1) - FLUXES(n,3)) * f_root
      ! allocation of ABG C to leaves
      FLUXES(n,4) = (FLUXES(n,1) - FLUXES(n,3) - FLUXES(n,6)) * (1d0 - (pars(29)*(lai_out(n)/6d0)))
      ! allocation of ABG C to labile/stem using pars(29)
      ! Ostrem.et.al.2013 (10.1080/09064710.2013.819440)            
      FLUXES(n,5) = FLUXES(n,1) - FLUXES(n,3) - FLUXES(n,6) - FLUXES(n,4)

      ! Determine GSI based canopy phenology
      call gsi_phenology(nodays,n,deltat(n),gsi_lag,             &
                         met(10,n),met(11,n),met(12,n),          &
                         pars(5),pars(10),pars(15),pars(30),     &
                         pars(12),pars(13),pars(14),pars(20),    & 
                         pars(21),pars(22),                      &
                         FLUXES(n,1),FLUXES(:,18),               &
                         FLUXES(n,15),FLUXES(n,16),FLUXES(n,17), &
                         FLUXES(n,8),FLUXES(n,14),               &
                         POOLS(n,1),POOLS(n,2))

      ! FLUXES WITH TIME DEPENDENCIES

      ! labile release = P_labile * (1-(1-leafgrowth)**deltat)/deltat
      FLUXES(n,7)  = POOLS(n,1)*(1d0-(1d0-FLUXES(n,14))**deltat(n))/deltat(n)
      ! leaf litter production = P_foliar * (1-(1-leaffall)**deltat)/deltat  
      FLUXES(n,9)  = POOLS(n,2)*(1d0-(1d0-FLUXES(n,8))**deltat(n))/deltat(n)
      ! root litter production = P_root * (1-(1-rootTOR)**deltat)/deltat  
      FLUXES(n,10) = POOLS(n,3)*(1d0-(1d0-pars(6))**deltat(n))/deltat(n)

      ! FLUXES WITH TEMP AND TIME DEPENDENCIES

      ! resp het litter = P_litter * (1-(1-GPP_respired*litterTOR)**deltat)/deltat  
      FLUXES(n,11) = POOLS(n,4)*(1d0-(1d0-FLUXES(n,2)*pars(7))**deltat(n))/deltat(n)
      ! resp het som = P_som * (1-(1-GPP_respired*somTOR)**deltat)/deltat
      FLUXES(n,12) = POOLS(n,5)*(1d0-(1d0-FLUXES(n,2)*pars(8))**deltat(n))/deltat(n)
      ! litter to som = P_litter * (1-(1-dec_rate*temprate)**deltat)/deltat
      FLUXES(n,13) = POOLS(n,4)*(1d0-(1d0-pars(1)*FLUXES(n,2))**deltat(n))/deltat(n)

      ! NEE = resp_auto + resp_het_litter + resp_het_som - GPP [i.e. '-' when CO2 sink '+' when CO2 source ]
      NEE(n) = (FLUXES(n,3) + FLUXES(n,11) + FLUXES(n,12)) - FLUXES(n,1)
      ! GPP 
      GPP(n) = FLUXES(n,1)

      ! update pools for next timestep

      ! labile pool = labile_pool[†-1] + (lab_prod - lab_cons)*deltat
      POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,7))*deltat(n)
      ! foliar pool = foliar_pool[†-1] + (leaf_prod + lab_prod2 - leaf_litter_prod)*deltat
      POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,4)+FLUXES(n,7)-FLUXES(n,9))*deltat(n)
      ! root pool = root_pool[†-1] + (root_prod - root_litter_prod)*deltat
      POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6)-FLUXES(n,10))*deltat(n)
      ! litter pool = litter_pool[†-1] + (leaf_litter_prod + root_litter_prod - resp_het_litter - litter2som)*deltat
      POOLS(n+1,4) = POOLS(n,4) + (FLUXES(n,9)+FLUXES(n,10)-FLUXES(n,11)-FLUXES(n,13))*deltat(n)
      ! som pool = som_pool[†-1] + (litter2som - resp_het_som)
      POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,13)-FLUXES(n,12))*deltat(n)

       !!!!!!!!!!
       ! Update soil water balance
       !!!!!!!!!!

       ! add any snow melt to the rainfall now that we have already dealt with the canopy interception
       rainfall = rainfall + snow_melt
       ! do mass balance (i.e. is there enough water to support ET)
       call calculate_update_soil_water(transpiration,soilevaporation,snowsublimation, &
                                        ((rainfall-intercepted_rainfall)*seconds_per_day) &
                                       ,FLUXES(n,46))
       ! now that soil mass balance has been updated we can add the wet canopy
       ! evaporation (kgH2O.m-2.day-1)
       FLUXES(n,46) = FLUXES(n,46) + wetcanopy_evap
       ! store soil water content of the surface zone (mm)
       POOLS(n+1,6) = 1d3 * soil_waterfrac(1) * layer_thickness(1)
       ! Assign all water variables to output variables (kgH2O/m2/day)
       FLUXES(n,47) =  transpiration   ! transpiration
       FLUXES(n,48) =  soilevaporation ! soil evaporation
       FLUXES(n,49) =  wetcanopy_evap  ! wet canopy evaporation
       FLUXES(n,50) =  runoff
       FLUXES(n,51) =  underflow

      ! CUTTING 
      ! ------------------------------------------------------------------------------------------------------------- ! 

      if (met(8,n) * pars(15) * 0.021 >= 1d0) then ! VM edit to allow cutting consideration 
          call grass_cutting(POOLS(n+1,1),POOLS(n+1,2),POOLS(n+1,3), & 
                             POOLS(n+1,4),POOLS(n+1,5),met(6,n),met(8,n), &
                             FLUXES(:,22),FLUXES(n,25),FLUXES(n,26),FLUXES(n,27), &
                             FLUXES(n,28),FLUXES(n,29),FLUXES(n,30), &
                             n,nodays,deltat(n), & 
                             pars(28),pars(33))
      end if

      ! GRAZING 
      ! ------------------------------------------------------------------------------------------------------------- ! 

      ! Determine whether there is any remaining losses to occur given GSI driven foliar loss - was drop met or management driven 
      gsi_lai_reduction = 0d0 ! FLUXES(n,9) * deltat(n) * pars(15)
      if (met(8,n)-gsi_lai_reduction > 0d0) then
          call grass_grazing(POOLS(n+1,1),POOLS(n+1,2),POOLS(n+1,3), & 
                             POOLS(n+1,4),POOLS(n+1,5),met(6,n),met(8,n)-gsi_lai_reduction, &
                             FLUXES(n,19),FLUXES(n,20),FLUXES(n,21), &
                             FLUXES(:,22),FLUXES(:,23),FLUXES(n,31),FLUXES(n,32), &
                             FLUXES(n,33),FLUXES(n,34),FLUXES(n,35),FLUXES(n,36), &
                             n,nodays,deltat(n), & 
                             pars(15),pars(27),pars(32),pars(34))
      end if                      

    end do ! nodays loop

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  ! Subroutines below this line
  !
  !-----------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------
  !
  subroutine acm_gpp_stage_1

    ! Estimate the light and temperature limited photosynthesis components.
    ! See acm_gpp_stage_2() for estimation of CO2 supply limitation and
    ! combination of light, temperature and CO2 co-limitation

    implicit none

    ! Declare local variables
    double precision :: a, b, c, Pl_max, PAR_m2, airt_adj

    !
    ! Metabolic limited photosynthesis
    !

    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1 -> umolC/m2/s). Scaling from leaf to canopy
    ! scaled assumed to follow integral of light environment.
    metabolic_limited_photosynthesis = gC_to_umol*leaf_canopy_light_scaling*ceff*seconds_per_day_1 &
                                     * ((leafT - Vc_minT) / ((leafT - Vc_minT) + Vc_coef))         &
                                     * opt_max_scaling(pn_max_temp,pn_min_temp,pn_opt_temp,pn_kurtosis,leafT)

    !
    ! Light limited photosynthesis
    !

    ! Calculate light limted rate of photosynthesis (umolC.m-2.s-1, daylight) as a function
    ! light capture and leaf to canopy scaling on quantum yield (e0).
    !light_limited_photosynthesis = e0 * canopy_par_MJday
    light_limited_photosynthesis = e0 * canopy_par_MJday * dayl_seconds_1 * gC_to_umol

    !
    ! Stomatal conductance independent variables for diffusion limited
    ! photosynthesis
    !

    ! Canopy level boundary layer conductance unit change
    ! (m.s-1 -> mol.m-2.s-1) assuming sea surface pressure only.
    ! Note the ratio of H20:CO2 diffusion through leaf level boundary layer is
    ! 1.37 (Jones appendix 2). Note conversion to resistance for easiler merging
    ! with stomatal conductance in acm_gpp_stage_2).
    rb_mol_1 = (aerodynamic_conductance * convert_ms1_mol_1 * gb_H2O_CO2 * &
              leaf_canopy_wind_scaling) ** (-1d0)

    ! Arrhenious Temperature adjustments for Michaelis-Menten coefficients
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
    co2_half_sat   = arrhenious(kc_half_sat_25C,kc_half_sat_gradient,leafT)
    co2_comp_point = arrhenious(co2comp_sat_25C,co2comp_gradient,leafT)

    ! don't forget to return
    return

  end subroutine acm_gpp_stage_1
  !
  !------------------------------------------------------------------
  !
  double precision function acm_gpp_stage_2(gs)

    ! Combine the temperature (pn) and light (pl) limited gross primary productivity
    ! estimates with CO2 supply limited via stomatal conductance (gs).
    ! See acm_gpp_stage_1() for additional details on pn and pl calculation.

    implicit none

    ! declare input variables
    double precision, intent(in) :: gs

    ! declare local variables
    double precision :: pp, qq, mult, rc, pd

    !
    ! Combined diffusion limitation and carboxylation limited photosynthesis
    !

    ! Estimation of ci is based on the assumption that metabilic limited
    ! photosynthesis is equal to diffusion limited. For details
    ! see Williams et al, (1997), Ecological Applications,7(3), 1997, pp. 882–894

    ! Daily canopy conductance dertmined through combination of aerodynamic and
    ! stomatal conductances. Both conductances are scaled to canopy aggregate.
    ! aerodynamic conductance already in units of molCO2.m-2.s-1 (see acm_gpp_stage_1).
    ! Stomatal conductance scaled from mmolH2O to molCO2.
    ! The ratio of H20:CO2 diffusion is 1.646259 (Jones appendix 2).
    !
    ! Combining in series the stomatal and boundary layer conductances
    ! to make canopy resistence (s/m2/molCO2)
    rc = (gs*gs_H2Ommol_CO2mol) ** (-1d0) + rb_mol_1

    ! pp and qq represent limitation by metabolic (temperature & N) and
    ! diffusion (co2 supply) respectively
    pp = metabolic_limited_photosynthesis*rc ; qq = co2_comp_point-co2_half_sat
    mult = co2+qq-pp
    ! calculate internal CO2 concentration (ppm or umol/mol)
    ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(co2*qq-pp*co2_comp_point)))

    ! calculate CO2 limited rate of photosynthesis (umolC.m-2.s-1)
    ! Then scale to day light period as this is then consistent with the light
    ! capture period (1/24 = 0.04166667)
    !pd = ((co2-ci)/rc) * umol_to_gC * dayl_hours_fraction
    pd = ((co2-ci)/rc)

    !
    ! Estimate CO2 and light co-limitation
    !

    ! calculate combined light and CO2 limited photosynthesis (umolC/m2/s)
    acm_gpp_stage_2 = light_limited_photosynthesis*pd/(light_limited_photosynthesis+pd)

    !pp = acm_gpp_stage_2*rc ; mult = co2+qq-pp
    !! calculate internal CO2 concentration (ppm or umol/mol)
    !ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(co2*qq-pp*co2_comp_point)))

    ! sanity check
    if (acm_gpp_stage_2 /= acm_gpp_stage_2 .or. acm_gpp_stage_2 < 0d0) acm_gpp_stage_2 = 0d0

    ! don't forget to return
    return

  end function acm_gpp_stage_2
  !
  !------------------------------------------------------------------
  !
  subroutine average_leaf_conductance(gv_forced)

    !
    ! Subroutine calculates the forced conductance of water vapour for non-cylinder within canopy leaves (i.e. broadleaf)
    ! Free convection (i.e. that driven by energy balance) is negelected here due to the lack of an energy balance
    ! calculation in DALEC. Should a energy balance be added then this code could be expanded include free conductance
    ! Follows a simplified approach to that used in SPA (Smallman et al 2013).
    !

    implicit none

    ! arguments
    double precision, intent(out) :: gv_forced ! canopy conductance (m/s) for water vapour under forced convection

    ! local parameters
    double precision, parameter :: leaf_width_coef = 25d0, & ! (1/leaf_width) * 0.5,
                                                               ! where 0.5 accounts for one half
                                                               ! of the leaf used in water exchange
                                        leaf_width = 0.02d0    ! leaf width (m) (alternates 0.04, 0.08)
!                                                Pr = 0.72d0, & ! Prandtl number
!                                           Pr_coef = 1.05877d0 !1.18d0*(Pr**(0.33d0))
    ! local variables
    double precision :: &
              Sh_forced & ! Sherwood number under forced convection
                    ,Re   ! Reynolds number

    ! Sherwood number under forced convection. NOTE: 0.962 * Pr_coef = 1.018537
!    Sh_forced = 0.962d0*Pr_coef*(sqrt((leaf_width*canopy_wind)/kinematic_viscosity))
    Sh_forced = 1.018537d0*(sqrt((leaf_width*canopy_wind)/kinematic_viscosity))
    ! Estimate the the forced conductance of water vapour
    gv_forced = water_vapour_diffusion*Sh_forced*leaf_width_coef

  end subroutine average_leaf_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_aerodynamic_conductance

    !
    ! Calculates the aerodynamic or bulk canopy conductance (m.s-1). Here we
    ! assume neutral conditions due to the lack of an energy balance calculation
    ! in either ACM or DALEC. The equations used here are with SPA at the time
    ! of the calibration
    !

    implicit none

    ! local variables
    double precision :: local_lai, &
           mixing_length_momentum, & ! mixing length parameter for momentum (m)
            length_scale_momentum    ! length scale parameter for momentum (m)

    ! Restrict LAI used here to greater than a minium value which prevents un-realistic outputs
    local_lai = max(min_lai,lai)

    ! calculate the zero plane displacement and roughness length
    call z0_displacement(ustar_Uh,local_lai)
    ! calculate friction velocity at tower height (reference height ) (m.s-1)
    ! WARNING neutral conditions only; see WRF module_sf_sfclay.F for 'with
    ! stability versions'
    !    ustar = (wind_spd / log((tower_height-displacement)/roughl)) * vonkarman
    ustar = wind_spd * ustar_Uh

    ! both length scale and mixing length are considered to be constant within
    ! the canopy (under dense canopy conditions) calculate length scale (lc)
    ! for momentum absorption within the canopy; Harman & Finnigan (2007)
    ! and mixing length (lm) for vertical momentum within the canopy Harman & Finnigan (2008)
    length_scale_momentum = (4d0*canopy_height) / local_lai
    mixing_length_momentum = 2d0*(ustar_Uh**3)*length_scale_momentum

    ! based on Harman & Finnigan (2008); neutral conditions only
    call log_law_decay

    ! now we are interested in the within canopy wind speed,
    ! here we assume that the wind speed just inside of the canopy is most important.
    !canopy_wind = canopy_wind*exp((ustar_Uh*((canopy_height*0.5d0)-canopy_height))/mixing_length_momentum)
    ! Estimate canopy scaling factor for use with aerodynamic conductance.
    leaf_canopy_wind_scaling = exp((ustar_Uh/mixing_length_momentum)) &
                             / (ustar_Uh/mixing_length_momentum)

    ! calculate_soil_conductance
    call calculate_soil_conductance(mixing_length_momentum,local_lai)
    ! calculate leaf level conductance (m/s) for water vapour under forced convective conditions
    call average_leaf_conductance(aerodynamic_conductance)

  end subroutine calculate_aerodynamic_conductance
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_cosine_solar_zenith_angle

    implicit none

    ! Calculate some common variables needed for the Sellers (1985)
    ! parameteric appoximation

    ! Declare local parameters

    ! Solar angle for the hour of the day, assumed in daily model to be 10 am (or 2 pm)
    ! therefore, 2 = 14 - 12, where 14 = timestep of day and 12 = half the number
    ! hour-angle = ( hr - 12 ) * 15. * ( pi / 180.d0 )
    ! of time steps per day.
    double precision, parameter :: cos_hour_angle = cos(0.5235988d0)

    ! Declare local variables

    ! Now calculate the solar-zenith-angle, as per..
    !  sin(latitude)sin(declination) + cos(latitude)cos(declination)cos(hour_angle)
    ! where hour-angle = ( hr - 12 ) * 15. * ( pi / 180.d0 )
    ! and from that the solar flux.
    ! eqn 2.15 in Hartmann's Global Physical Climatology...
    ! Minimum allowed value constraint from SIB3 implementation of Sellers (1985)
    cosine_solar_zenith_angle = max(0.01747d0, sin_latitude_radians * sin(declination) + &
                                               cos_latitude_radians * cos(declination) * &
                                               cos_hour_angle)

  end subroutine calculate_cosine_solar_zenith_angle
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_daylength

    ! Subroutine uses day of year and latitude (-90 / 90 degrees) as inputs,
    ! combined with trigonomic functions to calculate day length in hours and seconds

    implicit none

    ! local variables
    double precision :: dec, sinld, cosld, aob

    !
    ! Estimate solar geometry variables needed
    !

    ! day length is estimated as the ratio of sin and cos of the product of declination an latitude in radiation
    sinld = sin_latitude_radians * sin( declination )
    cosld = cos_latitude_radians * cos( declination )
    aob = max(-1d0,min(1d0,sinld / cosld))

    ! estimate day length in hours and seconds and upload to module variables
    dayl_hours = 12d0 * ( 1d0 + 2d0 * asin( aob ) * pi_1 )
    dayl_seconds = dayl_hours * seconds_per_hour

    ! return to user
    return

  end subroutine calculate_daylength
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_diffuse_fraction(diffuse_fraction)

    implicit none

    ! Description estimate ratio of actual to extra-solar radiation
    ! for day length period following Erbs et al., (1982)

    ! Declare arguments
    double precision, intent(out) :: diffuse_fraction

    ! Declare local parameters
    double precision, parameter :: So = 117.504d0 ! solar constant (1360 Wm-2)
                                                  ! unit scaled to MJ/m2/day
    ! Declare local variables
    double precision :: clear_day_swrad, Kt

    ! Estimate sunset solar angle
    sunset_solar_angle = acos(-tan(latitude_radians)*tan(declination))
    ! Estimate clear day light solar radiation (MJ/m2/d)
    clear_day_swrad = So * pi_1 * (1d0+0.033d0*cos((360d0*doy)/365d0)) &
                    * (cos_latitude_radians*cos(declination)*sin(sunset_solar_angle) &
                      + sin_latitude_radians*sin(declination))
    ! Estimate the ratio of actual to clear day radiation (MJ/m2/d over MJ/m2/d)
    Kt = swrad / clear_day_swrad

    ! Calculate diffuse ratio
    if (sunset_solar_angle < 1.4208d0) then
        if (Kt < 0.715d0) then
            diffuse_fraction = 1d0 - 0.2727d0*Kt + 2.4495d0*Kt**2 - 11.9514d0*Kt**3 + 9.3879d0*Kt**4
        else
            diffuse_fraction = 0.143d0
        end if
    else
        if (Kt < 0.722d0) then
            diffuse_fraction = 1d0 + 0.28832d0*Kt - 2.5557d0*Kt**2 + 0.8448d0*Kt**3
        else
            diffuse_fraction = 0.175d0
        end if
    end if

    ! Sanity check
    diffuse_fraction = min(1d0, diffuse_fraction)

    ! return
    return

  end subroutine calculate_diffuse_fraction
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_field_capacity

    ! field capacity calculations for saxton eqns !

    implicit none

    ! local variables..
    integer :: i
    double precision :: x1, x2

    x1 = 0.1d0 ; x2 = 0.7d0 ! low/high guess
    do i = 1 , nos_soil_layers+1
       water_retention_pass = i
       ! field capacity is water content at which SWP = -10 kPa
       field_capacity(i) = zbrent('water_retention:water_retention_saxton_eqns', &
                                   water_retention_saxton_eqns , x1 , x2 , 0.001d0, 0d0 )
    enddo

  end subroutine calculate_field_capacity
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_longwave_isothermal(canopy_temperature,soil_temperature)

    ! Subroutine estimates the isothermal net longwave radiation (W.m-2) for
    ! the canopy and soil surface. SPA uses a complex multi-layer radiative
    ! transfer scheme including reflectance, transmittance any absorption.
    ! However, for a given canopy vertical profiles, the LAI absorption
    ! relationship is readily predicted via Michaelis-Menten or
    ! non-rectangular hyperbola as done here.

    implicit none

    ! arguments
    double precision, intent(in) :: canopy_temperature, soil_temperature ! oC

    ! local variables
    double precision ::    dT, & ! Canopy transmittance for long wave radiation
                        lwrad, & ! downward long wave radiation from sky (W.m-2)
        longwave_release_soil, & ! emission of long wave radiation from surfaces per m2
      longwave_release_canopy    ! assuming isothermal condition (W.m-2)

    ! estimate long wave radiation from atmosphere (W.m-2)
    lwrad = emiss_boltz * (maxt+freeze-20d0) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_soil = emiss_boltz * (soil_temperature+freeze) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_canopy = emiss_boltz * (canopy_temperature+freeze) ** 4
    ! Canopy transmittance for thermal radiation
    dT = 1d0-exp(-lai/Vc*mu_obar)

    !!!!!!!!!!
    ! Isothermal net long wave canopy and soil balance (W.m-2)
    !!!!!!!!!!

    ! Diffuse longwave absorbed by the canopy
    canopy_lwrad_Wm2 = (lwrad*Vc*dT) - (Vc*dT*2d0*longwave_release_canopy) + (Vc*dT*longwave_release_soil)
    ! Diffuse longwave absorbed by the soil
    soil_lwrad_Wm2 = (lwrad*(1d0-(Vc*dT))) + (Vc*dT*longwave_release_canopy) - longwave_release_soil

  end subroutine calculate_longwave_isothermal  
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_radiation_balance

    implicit none

    !
    ! Parametric approximation for 2-stream, multi-layer radiative transfer
    !

    !
    ! References
    !
    ! Sellers (1985) Canopy reflectance, photosynthesis and transpiration.
    !                International Journal of Remote Sensing, 6(8), 1335-1772, doi: 10.1080/01431168508948283
    ! Sellers et al., (1986) A simple biosphere model (SiB) for use in general circulation models.
    !                        Journal of Atmospheric Sciences, 43(6), 505-531
    ! Sellers et al., (1996) A revised land surface parameterisation (SiB2) for atmospheric GCMs. Part 1: Model formualtion.
    !                        Journal of Climate, 9(4), 676-705, doi: 10.1175/1520_0442(1996)009<0676:ARLSPF>2.0.CO;2

    !
    ! Description
    ! The purpose of these equations is the approximation of the behaviuour of a 2-stream multiple canopy layer model
    ! accounting for solar angle, within canopy interception, reflectance and transmittance.
    ! As a result the actual parametric inputs are few, but requires a large number of empiracal relationships to describe
    ! the non-linear dynamics of radiative transfer
    ! NOTE: This code provides a daily timescale linear correction on
    ! isothermal longwave balance to net based on soil surface incident shortwave
    ! radiation. This correction is drawn from the standard ACM-GPP-ET-v1 approach

    ! Declare local variables
    double precision :: delta_iso

    ! Calculate the declination
    declination = calculate_declination(doy)
    ! Calculate cosine zenith angle
    call calculate_cosine_solar_zenith_angle
    ! Estimate shortwave radiation balance
    call calculate_shortwave_balance
    ! Estimate isothermal long wave radiation balance
    call calculate_longwave_isothermal(meant,meant)
    ! Apply linear correction to soil surface isothermal->net longwave radiation
    ! balance based on absorbed shortwave radiation
    delta_iso = (soil_iso_to_net_coef_LAI * lai) + &
                (soil_iso_to_net_coef_SW * (soil_swrad_MJday * 1d6 * seconds_per_day_1)) + &
                 soil_iso_to_net_const
    ! In addition to the iso to net adjustment, SPA analysis shows that soil net never gets much below zero
    !soil_lwrad_Wm2 = max(-0.1d0,soil_lwrad_Wm2 + delta_iso)
    soil_lwrad_Wm2 = soil_lwrad_Wm2 + delta_iso
    ! Apply linear correction to canopy isothermal->net longwave radiation
    ! balance based on absorbed shortwave radiation
    delta_iso = (canopy_iso_to_net_coef_LAI * lai) + &
                (canopy_iso_to_net_coef_SW * (canopy_swrad_MJday * 1d6 * seconds_per_day_1)) + &
                canopy_iso_to_net_const
    canopy_lwrad_Wm2 = canopy_lwrad_Wm2 + delta_iso

  end subroutine calculate_radiation_balance
  !
  !----------------------------------------------------------------------
  !    
  !subroutine calculate_radiation_commons(lat,rad_pars)
  subroutine calculate_radiation_commons(lat)

    implicit none

    ! Description

    ! Declare arguments
    double precision, intent(in) :: lat!, & ! site latitude in degrees
                               !rad_pars(6) !

    ! Calculate some common variables and place into memory
    latitude = lat
    latitude_radians = lat * deg_to_rad
    sin_latitude_radians = sin(latitude_radians)
    cos_latitude_radians = cos(latitude_radians)

    ! Load canopy optical properties to their module variables
    canopy_reflectance(1)   = canopy_nir_reflectance  !rad_pars(1) ! canopy_nir_reflectance
    canopy_reflectance(2)   = canopy_par_reflectance  !rad_pars(2) ! canopy_par_reflectance
    canopy_transmittance(1) = canopy_nir_transmittance!rad_pars(3) ! canopy_nir_transmittance
    canopy_transmittance(2) = canopy_par_transmittance!rad_pars(4) ! canopy_par_transmittance
    soil_reflectance(1)     = soil_nir_reflectance    !rad_pars(5) ! soil_nir_reflectance
    soil_reflectance(2)     = soil_par_reflectance    !rad_pars(6) ! soil_par_reflectance

    ! Canopy scattering of incident light, varied by wavelength
    ! NOTE: if we want to put snow fall on the canopy in the model, then this
    ! line will need to move into shortwave_balance() or calculate_radiation_balance()
    ! to be recalculated with each update of the canopy reflectance and transmittances.
    canopy_scattering = canopy_reflectance + canopy_transmittance

    ! Two empirical functions related to the leaf_distribution_deviance
    ! used in the calcuation of several variable varying by day of year
    ! The coefficient notation used in Sellers (1985) indicate "empty set".
    ! For clarity and visual similarity, I've used the capital letter O
    O1 = 0.5d0-(0.633d0*leaf_distribution_deviance)-(0.33d0*leaf_distribution_deviance**2)
    O2 = 0.877d0*(1d0-(2d0*O1))

    ! The average inverse diffuse optical depth per unit leaf area.
    ! mu_obar indicates mu with over bar
    ! See notes or references for the original integral equation
    mu_obar = (1d0/O2) * ( 1d0-(O1/O2)*log((O1+O2)/O1) )
    ! Mean leaf angle deviation from the horizontal (radians)
    ! SPA default assumption is 30 degrees, where (pi/180) is the conversion to radians
    leaf_angle = 0.5235988 !30d0 * (pi/180d0)
    ! Analytical correction for leaf angle (radians) on light scattering within the canopy
    ! Notation is the verbal description of that used in Sellers (1985)
    cos2theta = (1d0 + cos(2d0*leaf_angle)) * 0.5d0

    ! Upward scattering as diffuse radiation, a function of canopy_transmittance, canopy_reflectance and leaf angle.
    ! Notation in Sellers (1985) is "c"
    cc = 0.5 * (canopy_reflectance + canopy_transmittance + (canopy_reflectance - canopy_transmittance) * cos2theta)
    ! The upward scattering fraction / coefficient for diffuse radiation.
    ! NOTE there is a direct radiation equivalent found below, noted as beta0
    beta = cc / canopy_scattering
    ! Downward scatting of diffuse radiation
    ! Note that notations from Sellers (1985) have been given double letters if only single character was used,
    ! or written word for greek notation
    bb = (1d0-(1d0-beta)*canopy_scattering)
    ! Extinction coefficient for diffuse radiation, related to absorption normalised by mu_obar
    ! Note that notations from Sellers (1985) have been given double letters if only single character was used,
    ! or written word for greek notation
    hh = (((bb**2 - cc**2))**(0.5d0)) / mu_obar

    ! Define the vegetated and covered soil (i.e. by litter) fractions
    ! Vc, could also be considered a canopy clumping factor
    Vc = 0.75d0 ; Vg = 0d0

    ! Return back to user
    return

  end subroutine calculate_radiation_commons
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_Rtot

    ! Purpose of this subroutine is to calculate the minimum soil-root hydraulic
    ! resistance input into ACM. The approach used here is identical to that
    ! found in SPA.

    ! local variables
    integer :: i, rooted_layer
    double precision :: transpiration_resistance,root_reach_local, &
                        slpa, mult, prev, exp_func!, root_depth_50, bonus
    double precision, dimension(2) :: Rcond_layer, &
                                    root_mass,  &
                                    root_length
    double precision, parameter :: rootdist_tol = 13.81551d0 ! log(1d0/rootdist_tol - 1d0) were rootdist_tol = 1d-6
                                  !rootdist_tol = 1d-6!, & ! Root density assessed for the max rooting depth
                                   !root_depth_frac_50 = 0.25d0 ! fractional soil depth above which 50 %
                                                               ! of the root mass is assumed to be located

    ! reset water flux
    total_water_flux = 0d0 ; water_flux_mmolH2Om2s = 0d0 ; wSWP = 0d0 ; rSWP = 0d0
    slpa = 0d0 ; root_length = 0d0 ; root_mass = 0d0 ; Rcond_layer = 0d0
    ! calculate soil depth to which roots reach
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10 (ish).
    !    transpiration_resistance = (gplant * lai)**(-1d0)
    ! Following Weiburg function, the potential conductance is reduced under increasing
    ! potential differences between the canopy and soil
    ! SOME THOUGHT NEEDED HERE ON HOW TO QUANTIFY THIS IMPACT ON GPLANT, AS STRESS VARIES DEPENDING ON THE ROOT ZONE LAYERS...
    !transpiration_resistance = gplant * exp(-(abs(SWP(1:nos_root_layers) - minlwp) / s1)**(s2))
    !transpiration_resistance = canopy_height / (transpiration_resistance * max(min_lai,lai))
    transpiration_resistance = canopy_height / (gplant * max(min_lai,lai))

    !!!!!!!!!!!
    ! calculate current steps soil hydraulic conductivity
    !!!!!!!!!!!

    ! seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
       call calculate_soil_conductivity(i,soil_waterfrac(i),soil_conductivity(i))
    end do ! soil layers

    !!!!!!!!!!!
    ! Calculate root profile
    !!!!!!!!!!!

    ! The original SPA src generates an exponential distribution which aims
    ! to maintain 50 % of root biomass in the top 25 % of the rooting depth.
    ! In a simple 3 root layer system this can be estimates more simply

!    ! top 25 % of root profile
!    root_depth_50 = root_reach * root_depth_frac_50
!    if (root_depth_50 <= layer_thickness(1)) then
!
!        ! Greater than 50 % of the fine root biomass can be found in the top
!        ! soil layer
!
!        ! Start by assigning all 50 % of root biomass to the top soil layer
!        root_mass(1) = fine_root_biomass * 0.5d0
!        ! Then quantify how much additional root is found in the top soil layer
!        ! assuming that the top 25 % depth is found somewhere within the top
!        ! layer
!        bonus = (fine_root_biomass-root_mass(1)) &
!              * (layer_thickness(1)-root_depth_50) / (root_reach - root_depth_50)
!        root_mass(1) = root_mass(1) + bonus
!        ! partition the remaining root biomass between the seconds and third
!        ! soil layers
!        if (root_reach > sum(layer_thickness(1:2))) then
!            root_mass(2) = (fine_root_biomass - root_mass(1)) &
!                         * (layer_thickness(2)/(root_reach-layer_thickness(1)))
!            root_mass(3) = fine_root_biomass - sum(root_mass(1:2))
!        else
!            root_mass(2) = fine_root_biomass - root_mass(1)
!        endif
!
!    else if (root_depth_50 > layer_thickness(1) .and. root_depth_50 <= sum(layer_thickness(1:2))) then
!
!        ! Greater than 50 % of fine root biomass found in the top two soil
!        ! layers. We will divide the root biomass uniformly based on volume,
!        ! plus bonus for the second layer (as done above)
!        root_mass(1) = fine_root_biomass * (layer_thickness(1)/root_depth_50)
!        root_mass(2) = fine_root_biomass * ((root_depth_50-layer_thickness(1))/root_depth_50)
!        root_mass(1:2) = root_mass(1:2) * 0.5d0
!
!        ! determine bonus for the seconds layer
!        bonus = (fine_root_biomass-sum(root_mass(1:2))) &
!              * ((sum(layer_thickness(1:2))-root_depth_50)/(root_reach-root_depth_50))
!        root_mass(2) = root_mass(2) + bonus
!        root_mass(3) = fine_root_biomass - sum(root_mass(1:2))
!
!    else
!
!        ! Greater than 50 % of fine root biomass stock spans across all three
!        ! layers
!        root_mass(1:2) = fine_root_biomass * 0.5d0 * (layer_thickness(1:2)/root_depth_50)
!        root_mass(3) = fine_root_biomass - sum(root_mass(1:2))
!
!    endif
!    ! now convert root mass into lengths
!    root_length = root_mass * root_mass_length_coef_1
!!    root_length = root_mass / (root_density * root_cross_sec_area)

    !!!!!!!!!!!
    ! Calculate hydraulic properties and each rooted layer
    !!!!!!!!!!!

    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! NOTE: Depth correction already accounted for in soil resistance
    ! calculations and this is the maximum potential rate of transpiration
    ! assuming saturated soil and leaves at their minimum water potential.
    demand = max(0d0, (SWP(1:nos_root_layers) - (head*canopy_height)) - minlwp )
    ! now loop through soil layers, where root is present
    rooted_layer = 1
    ! Determine the exponential coefficient needed for an exponential decay to the current root reach
    ! Exponential decay profile following:
    ! Y = 1 / (1 + exp(-B * Z)), where Y = density at Z, B = gradient, Z = depth
    ! To determine gradient for current maximum root depth assuming density reaches rootdist_tol value, rearranges to:
    ! B = ln(1/Y - 1) / Z
    !d = seq(0,2, 0.01) ; c = -2.6 ; rmax = 1 ; d50 = 0.25 ; rd = rmax / (1+ (d/d50)**c)
    !rmax = rd * (1 + (d/d50)**c)
    !(((rmax / rd) - 1)**(1/c)) * d50 = d ! Depth at which rd = 99 %
    !slpa = log(1d0/rootdist_tol - 1d0) / root_reach
    slpa = rootdist_tol / root_reach
    prev = 1d0
    do i = 1, nos_root_layers
       ! Determine the exponential function for the current cumulative depth
       exp_func = exp(-slpa * sum(layer_thickness(1:i)))
       ! Calculate the difference in the integral between depths, i.e. the proportion of root in the current volume
       mult = prev - (1d0 - (1d0/(1d0+exp_func)) + (0.5d0 * exp_func))
       ! Assign fine roo the the current layer...
       root_mass(i) = fine_root_biomass * mult
       ! and determine the associated amount of root
       root_length(i) = root_mass(i) * root_mass_length_coef_1
       prev = prev - mult
       ! If there is root in the current layer then we should calculate the resistances
       if (root_mass(i) > 0d0) then
           ! Track the deepest root layer assessed
           rooted_layer = i
           ! if there is root then there is a water flux potential...
           root_reach_local = min(root_reach,layer_thickness(i))
           ! calculate and accumulate steady state water flux in mmol.m-2.s-1
           call plant_soil_flow(i,root_length(i),root_mass(i) &
                               ,demand(i),root_reach_local &
                               ,transpiration_resistance,Rcond_layer(i))
       else
           ! ...if there is not then we wont have any below...
           exit
       end if ! root present in current layer?
    end do ! nos_root_layers
    ! Turn the output resistance into conductance
    Rcond_layer = Rcond_layer**(-1d0)

    ! if freezing then assume soil surface is frozen, therefore no water flux
    if (meant < 1d0) then
        water_flux_mmolH2Om2s(1) = 0d0
        Rcond_layer(1) = 0d0
    end if

    ! calculate sum value (mmolH2O.m-2.s-1)
    total_water_flux = sum(water_flux_mmolH2Om2s)
    ! wSWP based on the conductance due to the roots themselves.
    ! The idea being that the plant may hedge against growth based on the majority of the
    ! profile being dry while not losing leaves within some toleration.
    rSWP = sum(SWP(1:rooted_layer) * (Rcond_layer(1:rooted_layer) / sum(Rcond_layer(1:rooted_layer))))
    if (total_water_flux <= vsmall) then
        ! Set values for no water flow situation
        uptake_fraction = (layer_thickness(1:nos_root_layers) / sum(layer_thickness(1:nos_root_layers)))
        ! Estimate weighted soil water potential based on fractional extraction from soil layers
        wSWP = sum(SWP(1:nos_root_layers) * uptake_fraction(1:nos_root_layers))
        total_water_flux = 0d0
      else
        ! calculate weighted SWP and uptake fraction
        uptake_fraction(1:nos_root_layers) = water_flux_mmolH2Om2s(1:nos_root_layers) / total_water_flux
        ! Estimate weighted soil water potential based on fractional extraction from soil layers
        wSWP = sum(SWP(1:nos_root_layers) * uptake_fraction(1:nos_root_layers))
    endif

    ! and return
    return

  end subroutine calculate_Rtot
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_shortwave_balance

    ! Subroutine estimates the canopy and soil absorbed shortwave radiation
    ! (MJ/m2/day). Radiation absorption is paritioned into NIR and PAR for
    ! canopy, and NIR + PAR for soil.
    ! Follows an implementation of the Sellers (1985) approximation

    implicit none

    ! Declare local parameters
    double precision, dimension(no_wavelength), parameter :: &
                                   newsnow_reflectance = (/0.73d0,0.95d0/) ! NIR/PAR new snow reflectance fraction
    ! local variables
    double precision :: Gu, K, mu, S2, fsnow, diffuse_fraction
    ! Local variables with different values per wavelength
    double precision, dimension(no_wavelength) :: &
                      as_mu, dd, ff, sigma, u1, u2, u3, &
                      S1, p1, p2, p3, p4, D1, D2, &
                      h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, &
                      Iup, Idown, soil_albedo, &
                      canopy_absorption_fraction_diffuse, &
                      canopy_absorption_fraction_direct, &
                      soil_absorption_fraction_diffuse, &
                      soil_absorption_fraction_direct, &
                      soil_nir_par_MJday, canopy_nir_par_MJday, &
                      swrad_direct, swrad_diffuse

    ! Estimate the diffuse fraction of shortwave radiation
    call calculate_diffuse_fraction(diffuse_fraction)
    ! Estimate multiple use par and nir components (Note: units are MJ/m2/d)
    swrad_diffuse(1) = (1d0 - sw_par_fraction) * swrad * diffuse_fraction      ! NIR
    swrad_diffuse(2) = sw_par_fraction * swrad * diffuse_fraction              ! PAR
    swrad_direct(1) = (1d0 - sw_par_fraction) * swrad * (1d0-diffuse_fraction) ! NIR
    swrad_direct(2) = sw_par_fraction * swrad * (1d0-diffuse_fraction)         ! PAR

    ! Assign cosine_solar_zenith_angle to a local variable for easier readability
    mu = cosine_solar_zenith_angle

    ! Relative projected area of leaf elements in direction of the cosine_solar_zenith_angle (mu).
    ! This variable is determined as the result of two empirical functions related to the leaf_distribution_deviance
    ! Note the notation used here Gu is varied, for clarity, from the actual used in Sellers (1985) which is G(mu).
    ! Similarly, the coefficients used in calculating Gu indicate "empty set".
    ! For clarity and visual similarity, I've used the capital letter O
    Gu = O1 + (O2 * mu)

    ! Extinction coefficient for direct radiation, related to Gu and the cosine zenith angle
    K = Gu / mu
    ! Single leaf scattering albedo within the canopy, varied by mu, leaf distribution and wavelength
    as_mu = ((canopy_scattering * 0.5d0) * (Gu / (Gu+(mu*O2)))) &
          * (1d0 - (mu*(O1/(Gu+(mu*O2)))*log((Gu+(mu*O2)+(mu*O1))/(mu*O1)) ) )

    ! Upscatting coefficient for direct radiation, varied by mu and wavelength
    beta0 = ((1d0+(mu_obar*K)) / (canopy_scattering*mu_obar*K)) * as_mu

    ! Various terms, yet to have their specific functions determined
    ! Note that notations from Sellers (1985) have been given double letters if only single character was used,
    ! or written word for greek notation
    dd = canopy_scattering * mu_obar * K * beta0
    ff = canopy_scattering * mu_obar * K * (1-beta0)
    sigma = cc**2 + bb**2 + (mu_obar*K)**2
    u1 = bb - (cc/soil_reflectance) ; u2 = bb - (cc*soil_reflectance) ; u3 = ff + (cc*soil_reflectance)
    S1 = exp(-hh*lai) ; S2 = exp(-K*lai)

    ! Related to diffuse radiation
    p1 = bb + (mu_obar*hh) ; p2 = bb - (mu_obar*hh)
    ! Related to direct radiation
    p3 = bb + (mu_obar*K) ; p4 = bb - (mu_obar*K)
    !
    D1 = (p1 * (u1 - (mu_obar*hh)) * (1d0/S1)) - (p2*(u2+(mu_obar*hh))*S1)
    D2 = ((u2+(mu_obar*hh))*(1d0/S1)) - ((u2-(mu_obar*hh))*S1)

    !
    ! Direct radiation specific components
    !

    h1 = (-dd*p4) - (cc*ff)
    h2 =  (1d0/D1) * ( ((dd-((h1/sigma)*p3))*((u1-(mu_obar*hh))*(1d0/S1))) &
                     - (p2*(dd-cc-((h1/sigma)*(u1+(mu_obar*K))))*S2) )
    h3 = (-1d0/D1) * ( ((dd-((h1/sigma)*p3))*(u1+(mu_obar*hh))*S1) &
                     - (p1*(dd-cc-((h1/sigma)*(u1+(mu_obar*K))))*S2) )
    h4 = (-dd*p3) - (cc*ff) ! NOTE: "-" at the beginning is a correction identified in Sellers et al., (1996)
    h5 = (-1d0/D2) * ( ((h4/sigma)*(u2+(mu_obar*hh))*(1d0/S1)) &
                     + (u3-((h4/sigma)*(u2-(mu_obar*K))*S2)) )
    h6 = (1d0/D2) * ( ((h4/sigma)*(u2-(mu_obar*hh))*S1) &
                     + (u3-((h4/sigma)*(u2-(mu_obar*K))*S2)) )

    ! Fraction of direct radiation which leaves the canopy top as diffuse
    Iup = ((h1*exp(-K*lai))/sigma) + (h2*exp(-hh*lai)) + (h3*exp(hh*lai))
    ! Fraction of direct radiation which leaves the canopy base as diffuse
    Idown = ((h4*exp(-K*lai))/sigma) + (h5*exp(-hh*lai)) + (h6*exp(hh*lai))

    ! Update soil reflectance based on snow cover
    if (snow_storage > 0d0) then
        fsnow = 1d0 - exp( - snow_storage * 1d-2 )  ! fraction of snow cover on the ground
        soil_albedo = ((1d0 - fsnow) * soil_reflectance) + (fsnow * newsnow_reflectance) ! NIR & PAR
    else
        ! Estimate the combined soil and surface layer reflectances, assume direct and diffuse reflectances are isotropic
        ! NOTE; this variable is also used in the diffuse calculation too
        !soil_albedo = (soil_surface_reflectance * Vg) + ((1d0-Vg) * soil_reflectance)
        soil_albedo = soil_reflectance
    endif

    ! Fraction of direct radiation absorbed by the canopy
    canopy_absorption_fraction_direct = Vc * (1d0 - Iup - (Idown*(1-soil_albedo)) &
                                             - (exp(-K*lai/Vc)*(1-soil_albedo)))
    ! Fraction of direct radiation absorbed by the soil
    soil_absorption_fraction_direct = ((1d0-Vc)*(1-soil_albedo)) &
                                    + (Vc*((Idown*(1d0-soil_albedo)) + (exp(-K*lai/Vc)*(1-soil_albedo))))

    !
    ! Diffuse radiation specific components
    !

    h7  = (cc/D1) * (u1-(mu_obar*hh)) * (1d0/S1)
    h8  = (-cc/D1) * (u1+(mu_obar*hh)) * S1
    h9  = (1d0/D2) * (u2+(mu_obar*hh)) * (1d0/S1)
    h10 = (-1d0/D2) * (u2-(mu_obar*hh)) * S1

    ! Fraction of diffuse radiation which leaves the canopy top as diffuse
    Iup = (h7*exp(-hh*lai)) + (h8*exp(hh*lai))
    ! Fraction of diffuse radiation which leaves the canopy base as diffuse
    Idown = (h9*exp(-hh*lai)) + (h10*exp(hh*lai))

    ! Fraction of diffuse radiation absorbed by the canopy
    canopy_absorption_fraction_diffuse = Vc * (1d0 - Iup - (Idown*(1d0-soil_albedo)))
    ! Fraction of diffuse radiation absorbed by the soil
    soil_absorption_fraction_diffuse = ((1d0-Vc)*(1d0-soil_albedo)) + (Vc*((Idown*(1d0-soil_albedo))))

    !
    ! Combine direct and diffuse, convert into actual units of energy (MJ/m2/d)
    !

    ! Determine the combined direct and diffuse absorptions across wavelength
    soil_nir_par_MJday = (soil_absorption_fraction_diffuse * swrad_diffuse) &
                       + (soil_absorption_fraction_direct * swrad_direct)
    canopy_nir_par_MJday = (canopy_absorption_fraction_diffuse * swrad_diffuse) &
                       + (canopy_absorption_fraction_direct * swrad_direct)
    ! Assign to specific output variables for use elsewhere in the model
    canopy_par_MJday = canopy_nir_par_MJday(2)
    ! combine totals for use is soil evaporation
    soil_swrad_MJday = sum(soil_nir_par_MJday)
    ! Combine to estimate total shortwave canopy absorbed radiation
    canopy_swrad_MJday = sum(canopy_nir_par_MJday)

    ! Estimate the integral of light interception for use as a leaf to canopy
    ! scaler for photosynthesis, transpiration, and gs
    leaf_canopy_light_scaling = (1d0-S2) / K

    ! check energy balance
!    balance = swrad - canopy_par_MJday - canopy_nir_MJday - refl_par_MJday - refl_nir_MJday - soil_swrad_MJday
!    if ((balance - swrad) / swrad > 0.01) then
!        print*,"SW residual frac = ",(balance - swrad) / swrad,"SW residual = ",balance,"SW in = ",swrad
!   endif

  end subroutine calculate_shortwave_balance
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_soil_conductance(lm,local_lai)

    ! proceedsure to solve for soil surface resistance based on Monin-Obukov
    ! similarity theory stability correction momentum & heat are integrated
    ! through the under canopy space and canopy air space to the surface layer
    ! references are Nui & Yang 2004; Qin et al 2002
    ! NOTE: conversion to conductance at end

    implicit none

    ! declare arguments
    double precision, intent(in) :: lm, local_lai

    ! local variables
    double precision :: canopy_decay & ! canopy decay coefficient for soil exchange
                       ,Kh_canht       ! eddy diffusivity at canopy height (m2.s-1)

    ! parameters
    double precision, parameter :: foliage_drag = 0.2d0 ! foliage drag coefficient

    ! calculate eddy diffusivity at the top of the canopy (m2.s-1)
    ! Kaimal & Finnigan 1994; for near canopy approximation
    Kh_canht = vonkarman*ustar*(canopy_height-displacement)

    ! calculate canopy decay coefficient with stability correction
    ! NOTE this is not consistent with canopy momentum decay done by Harman &
    ! Finnigan (2008)
    canopy_decay = sqrt((foliage_drag*canopy_height*local_lai)/lm)

    ! approximation of integral for soil resistance (s/m) and conversion to
    ! conductance (m/s)
    soil_conductance = ( canopy_height/(canopy_decay*Kh_canht) &
                       * (exp(canopy_decay*(1d0-(soil_roughl/canopy_height)))- &
                          exp(canopy_decay*(1d0-((roughl+displacement)/canopy_height)))) ) ** (-1d0)

    return

  end subroutine calculate_soil_conductance  
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_soil_conductivity(soil_layer,waterfrac,conductivity)

    ! Calculate the soil conductivity (m s-1) of water based on soil
    ! characteristics and current water content

    implicit none

    ! arguments
    integer, intent(in) :: soil_layer
    double precision, intent(in) :: waterfrac
    double precision, intent(out) :: conductivity

    ! soil conductivity for the dynamic soil layers (i.e. not including core)
    conductivity = cond1(soil_layer) * exp(cond2(soil_layer)+cond3(soil_layer)/waterfrac)

    ! protection against floating point error
    if (waterfrac < 0.05d0) conductivity = 1d-30

  end subroutine calculate_soil_conductivity
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_soil_evaporation(soilevap)

    ! Estimate soil surface evaporation based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(out) :: soilevap ! kgH2O.m-2.day-1

    ! local variables
    double precision :: local_temp &
                   ,soil_radiation & ! isothermal net radiation (W/m2)
                            ,esurf & ! see code below
                             ,esat & ! soil air space saturation vapour pressure
                              ,gws   ! water vapour conductance through soil air space (m.s-1)

    ! oC -> K for local temperature value
    local_temp = maxt + freeze

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    soil_radiation = soil_lwrad_Wm2 + (soil_swrad_MJday * 1d6 * dayl_seconds_1)

    !!!!!!!!!!
    ! Calculate soil evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! calculate saturated vapour pressure (kPa), function of temperature.
    esat = 0.1d0 * exp( 1.80956664d0 + ( 17.2693882d0 * local_temp - 4717.306081d0 ) / ( local_temp - 35.86d0 ) )
    air_vapour_pressure = esat - vpd_kPa

    ! Soil conductance to water vapour diffusion (m s-1)...
    gws = porosity(1) * water_vapour_diffusion / (tortuosity*drythick)

    ! vapour pressure in soil airspace (kPa), dependent on soil water potential
    ! - Jones p.110. partial_molar_vol_water. Less vapour pressure of the air to
    ! estimate the deficit between soil and canopy air spaces
    esurf = (esat * exp( 1d6 * SWP(1) * partial_molar_vol_water / (Rcon * local_temp) )) - air_vapour_pressure

    ! Estimate potential soil evaporation flux (kgH2O.m-2.day-1)
    soilevap = ( ((slope*soil_radiation) + (air_density_kg*cpair*esurf*soil_conductance)) &
               / (lambda*(slope+(psych*(soil_conductance/gws)))) ) * dayl_seconds

    return

  end subroutine calculate_soil_evaporation
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_stomatal_conductance

    ! Determines an approximation of canopy scale stomatal conductance (gc)
    ! mmolH2O.m-2.s-1 based on potential hydraulic flow, air temperature and absorbed radiation.

    implicit none

    ! local variables
    double precision :: denom, iWUE_upper!, iWUE_lower
    double precision, parameter :: max_gs = 1000d0, &  ! mmolH2O.m-2.s-1 (leaf area)
                                   min_gs = 0.01d0, &  ! mmolH2O.m-2.s-1 (leaf area)
                                   tol_gs = 0.01d0     ! mmolH2O.m-2.s-1 (leaf area)

    !!!!!!!!!!
    ! Calculate stomatal conductance under H2O and CO2 limitations
    !!!!!!!!!!

    if (aerodynamic_conductance > vsmall .and. total_water_flux > vsmall .and. leafT > Vc_minT) then

        ! Determine potential water flow rate (mmolH2O.m-2.s-1)
        max_supply = total_water_flux

        ! Pass minimum conductance from local parameter to global value
        minimum_conductance = min_gs * leaf_canopy_light_scaling

        ! Invert Penman-Monteith equation to give gs (m.s-1) needed to meet
        ! maximum possible evaporation for the day.
        ! This will then be reduced based on CO2 limits for diffusion based
        ! photosynthesis
        denom = slope * (((canopy_swrad_MJday * 1d6 * dayl_seconds_1) + canopy_lwrad_Wm2)) &
              + (ET_demand_coef * aerodynamic_conductance * leaf_canopy_wind_scaling)
        denom = (denom / (lambda * max_supply * mmol_to_kg_water)) - slope
        potential_conductance = (aerodynamic_conductance * leaf_canopy_wind_scaling) / (denom / psych)

        ! convert m.s-1 to mmolH2O.m-2.s-1, per unit ground area, note that this
        ! is implicitly the canopy scaled value
        potential_conductance = potential_conductance * convert_ms1_mmol_1
        ! if conditions are dew forming then set conductance to maximum as we
        ! are not going to be limited by water demand
        if (potential_conductance <= 0d0 .or. potential_conductance > max_gs*leaf_canopy_light_scaling) then
            potential_conductance = max_gs*leaf_canopy_light_scaling
        end if

        ! If there is a positive demand for water then we will solve for
        ! photosynthesis limits on gs through iterative solution

        ! Determine the appropriate canopy scaled gs increment and return threshold
        delta_gs = 1d0 * leaf_canopy_light_scaling ! mmolH2O/m2leaf/s
        iWUE_step = iWUE * leaf_canopy_light_scaling ! umolC/mmolH2Ogs/s

        ! Calculate stage one acm, temperature and light limitation which
        ! are independent of stomatal conductance effects
        call acm_gpp_stage_1

!        if (do_iWUE) then
            ! Intrinsic WUE optimisation
            ! Check that the water restricted water range brackets the root solution for the bisection
            iWUE_upper = find_gs_iWUE(potential_conductance) !; iWUE_lower = find_gs_iWUE(min_gs)
            if ( iWUE_upper * find_gs_iWUE(min_gs) > 0d0 ) then
                ! Then both proposals indicate that photosynthesis
                ! would be increased by greater opening of the stomata
                ! and is therefore water is limiting!
                stomatal_conductance = potential_conductance
                ! Exception being if both are positive - therefore assume
                ! lowest
                if (iWUE_upper > 0d0) stomatal_conductance = minimum_conductance
            else

                ! In all other cases iterate
                stomatal_conductance = zbrent('calculate_gs:find_gs_iWUE', &
                                              find_gs_iWUE,minimum_conductance,potential_conductance,tol_gs*lai,iWUE_step*0.10d0)

            end if

    else

        ! if no LAI then there can be no stomatal conductance
        potential_conductance = max_gs ; minimum_conductance = vsmall
        stomatal_conductance = vsmall
        ! set minimum (computer) precision level flow
        max_supply = vsmall

    endif ! if aerodynamic conductance > vsmall

  end subroutine calculate_stomatal_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_transpiration(transpiration)

    ! Models leaf cnaopy transpiration based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(out) :: transpiration ! kgH2O.m-2.day-1

    ! local variables
    double precision :: canopy_radiation & ! isothermal net radiation (W/m2)
                                  ,gs,gb   ! stomatal and boundary layer conductance (m.s-1)

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * dayl_seconds_1)

    !!!!!!!!!!
    ! Calculate canopy conductance (to water vapour)
    !!!!!!!!!!

    ! Change units of potential stomatal conductance
    ! (mmolH2O.m-2.d-1 -> m.s-1).
    ! Note assumption of sea surface pressure only
 !   gs = (stomatal_conductance / convert_ms1_mmol_1) * dayl_seconds_1
    gs = stomatal_conductance / convert_ms1_mmol_1
    ! Scale aerodynamic conductance to canopy scale
    gb = aerodynamic_conductance * leaf_canopy_wind_scaling

    !!!!!!!!!!
    ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! Calculate numerator of Penman Montheith (kgH2O.m-2.day-1)
    ! NOTE: that restriction within water supply restriction is determined
    ! during stomatal conductance level.
    transpiration = ( ( (slope*canopy_radiation) + (ET_demand_coef*gb) ) &
                      / (lambda*(slope+(psych*(1d0+gb/gs)))) )*dayl_seconds

  end subroutine calculate_transpiration
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_update_soil_water(Eleaf,Esoil,Esnow,rainfall_in,corrected_ET)

    !
    ! Function updates the soil water status and layer thickness
    ! Soil water profile is updated in turn with evaporative losses,
    ! rainfall infiltration and gravitational drainage
    ! Root layer thickness is updated based on changes in the rooting depth from
    ! the previous step
    !

    implicit none

    ! arguments
    double precision, intent(in) :: rainfall_in   ! rainfall (kgH2O.m-2.day-1)
    double precision, intent(inout) :: Eleaf, Esoil, Esnow ! evapotranspiration estimate (kgH2O.m-2.day-1)
    double precision, intent(out) :: corrected_ET      ! water balance corrected evapotranspiration (kgH2O/m2/day)

    ! local variables
    integer :: day, a
    double precision :: depth_change, water_change, initial_soilwater, balance, mass_check, &
                        Esoil_local, Esnow_local
    double precision, dimension(2) :: avail_flux, evaporation_losses, pot_evap_losses
    !logical :: iter_soil = .true.

    ! set soil water exchanges
    Esoil = 0d0 ; Esnow = 0d0  
    underflow = 0d0 ; runoff = 0d0 ; corrected_ET = 0d0 ; evaporation_losses = 0d0 ; pot_evap_losses = 0d0
    initial_soilwater = 1d3 * sum(soil_waterfrac(1:nos_soil_layers) * layer_thickness(1:nos_soil_layers))

    !! Assume leaf transpiration is drawn from the soil based on the
    !! update_fraction estimated in calculate_Rtot
    pot_evap_losses = Eleaf * uptake_fraction
    !! Assume all soil evaporation comes from the soil surface only
    !pot_evap_losses(1) = pot_evap_losses(1) + Esoil

! Conditions under which iterative solution should not be needed...
! Scenario 1
!          (i) Soil layers at or below field capacity, therefore no drainage
!         (ii) The existing water supply and rainfall can support evaporative demanded by the canopy and soil
!     Outcome: Extract all needed water, potentially leaving soil in negetative status, followed by infilatration.
!              Allow for drainage if soil is above field capacity as a result of this proecss
! Scenario 2
!          (i) Soil layers ABOVE field capacity, therefore THERE is drainage
!         (ii) The existing water supply and rainfall can support evaporative demanded by the canopy and soil
!     Outcome: Extract allow water and add all infiltration into the soil.
!              Allow for drainage in the final instance as strongly exponential drainage flow should negate time difference.
!              NOTE: that this may bias between runoff and underflow estimation
! UPDATE (19/04/2024): Soil layer above or near field capacity do require daily iteration or other solution at a later time
!                      due to the non-linear effects on soil surface evaporation as a result of drythick development

    ! to allow for smooth water balance integration carry this out at daily time step
    ! NOTE: Should inestigate conditions where we don't need to do all the days.
    !       For example, can we check how the various terms are changing and whether 
    !       we can assume an average for the rest of the time step?
    do day = 1, nint(days_per_step)
       ! Possible conditions for avoiding looping all days
       ! 1) When drythick == min_drythick and rainfall_in > Esoil in first day
       ! 2) If initially, drythick > min_drythick, but rainfall_in > Esoil, iteration can stop once drythick == min_drythick
       ! 3) If drythick == min_drythick and rainfall_in < Esoil in first day, iteration can still be avoided if the time step multiple does not result in drythick >  min_drythick
       ! 4) 

       !!!!!!!!!!
       ! Evaporative losses
       !!!!!!!!!!

       ! Estimate drythick for the current step
       !drythick = max(min_drythick, top_soil_depth * max(0d0,(1d0 - (soil_waterfrac(1) / field_capacity(1)))))
       drythick = max(min_drythick, top_soil_depth * max(0d0,(1d0 - (soil_waterfrac(1) / porosity(1)))))
       ! Soil surface (kgH2O.m-2.day-1)
       call calculate_soil_evaporation(Esoil_local)

       ! If snow present assume that soil evaporation is sublimation of soil first
       if (Esoil_local > 0d0 .and. snow_storage > 0d0) then
           if (snow_storage > Esoil_local) then 
               Esnow_local = Esoil_local
               snow_storage = snow_storage - Esnow_local
               Esoil_local = 0d0
           else 
               Esnow_local = snow_storage
               Esoil_local = Esoil_local - Esnow_local
               snow_storage = 0d0
           end if
       else 
           Esnow_local = 0d0
       end if
       ! Accumulate overall time step soil and snow evaporation
       Esoil = Esoil + Esoil_local
       Esnow = Esnow + Esnow_local

       ! load transpiration losses from soil profile
       evaporation_losses = pot_evap_losses
       ! Update with current daily estimate of soil evaporation losses, surface only
       evaporation_losses(1) = evaporation_losses(1) + Esoil_local

       ! can not evaporate from soil more than is available (m -> mm)
       ! NOTE: This is due to the fact that both soil evaporation and transpiration
       !       are drawing from the same water supply.
       avail_flux = soil_waterfrac(1:nos_root_layers) * layer_thickness(1:nos_root_layers) * 1d3
       do a = 1, nos_root_layers ! note: timed comparison between "where" and do loop supports do loop for smaller vectors
          if (evaporation_losses(a) > avail_flux(a)) evaporation_losses(a) = avail_flux(a) * 0.9999d0
       end do
       ! this will update the ET estimate outside of the function
       ! days_per_step corrections happens outside of the loop below
       corrected_ET = corrected_ET + sum(evaporation_losses)

       ! adjust water already committed to evaporation
       ! convert kg.m-2 (or mm) -> Mg.m-2 (or m)
       soil_waterfrac(1:nos_root_layers) = soil_waterfrac(1:nos_root_layers) &
                                         + ((-evaporation_losses(1:nos_root_layers)*1d-3) / layer_thickness(1:nos_root_layers))

       ! Correct for dew formation; any water above porosity in the top layer is assumed runoff
       ! NOTE: layer_thickness * 1d3 scales between m3/m3 to kg/m2
       ! Worth investigating whether this term is actually important, dew formation in unlikely
       ! to be large. Therefore, maybe a conditional statment and calculation unrequired
       if (soil_waterfrac(1) > porosity(1)) then
           runoff = runoff + ((soil_waterfrac(1)-porosity(1)) * layer_thickness(1) * 1d3)
           soil_waterfrac(1) = porosity(1)
       endif

       !!!!!!!!!!
       ! Rainfall infiltration drainage
       !!!!!!!!!!

       ! Determine infiltration from rainfall (kgH2O/m2/day),
       ! if rainfall is probably liquid / soil surface is probably not frozen
       ! reset soil water change variable
       waterchange = 0d0
       call infiltrate(rainfall_in)
       ! update soil profiles. Convert fraction into depth specific values
       ! (rather than m3/m3) then update fluxes
       soil_waterfrac(1:nos_soil_layers) = soil_waterfrac(1:nos_soil_layers) &
                                         + (waterchange(1:nos_soil_layers) / layer_thickness(1:nos_soil_layers))
       ! soil waterchange variable reset in gravitational_drainage()

       !!!!!!!!!!
       ! Gravitational drainage
       !!!!!!!!!!

       ! Determine drainage flux between surface -> sub surface
       call gravitational_drainage(1)
!print*,day,rainfall_in,Esoil_local,Esnow_local,runoff/dble(day),underflow/dble(day)
    end do ! days_per_step

    ! apply time step correction kgH2O/m2/step -> kgH2O/m2/day
    corrected_ET = corrected_ET * days_per_step_1
    underflow = underflow * days_per_step_1
    runoff = runoff * days_per_step_1
    Esoil = Esoil * days_per_step_1
    Esnow = Esnow * days_per_step_1

    ! Based on the soil mass balance corrected_ET, make assumptions to correct Eleaf and Esnow
    balance = (corrected_ET-Eleaf) / Esoil
    if (balance < 1d0 .and. balance > 0d0) Esoil = Esoil * balance

    ! Update corrected_ET with snow sublimation
    corrected_ET = corrected_ET + Esnow

    !!!!!!!!!!
    ! Update soil layer thickness
    !!!!!!!!!!

    depth_change = (top_soil_depth+min_layer) ; water_change = 0
    ! if roots extent down into the bucket
    if (root_reach > depth_change .and. previous_depth <= depth_change) then

        !!!!!!!!!!
        ! Soil profile is within the bucket layer (layer 3)
        !!!!!!!!!!

        if (previous_depth > depth_change) then
            ! how much has root depth extended since last step?
            depth_change = root_reach - previous_depth
        else
            ! how much has root depth extended since last step?
            depth_change = root_reach - depth_change
        endif

        ! if there has been an increase
        if (depth_change > 0.05d0) then

            ! determine how much water (mm) is within the new volume of soil
            water_change = soil_waterfrac(nos_soil_layers) * depth_change
            ! now assign that new volume of water to the deep rooting layer
            soil_waterfrac(nos_root_layers) = ((soil_waterfrac(nos_root_layers)*layer_thickness(nos_root_layers))+water_change) &
                                            / (layer_thickness(nos_root_layers)+depth_change)

            ! explicitly update the soil profile if there has been rooting depth
            ! changes
            layer_thickness(1) = top_soil_depth
            layer_thickness(2) = root_reach - top_soil_depth
            layer_thickness(3) = max_depth - sum(layer_thickness(1:2))

            ! keep track of the previous rooting depth
            previous_depth = root_reach

        else if (depth_change < -0.05d0) then

            ! make positive to ensure easier calculations
            depth_change = -depth_change

            ! determine how much water is lost from the old volume of soil
            water_change = soil_waterfrac(nos_root_layers) * depth_change
            ! now assign that new volume of water to the deep rooting layer
            soil_waterfrac(nos_soil_layers) = ((soil_waterfrac(nos_soil_layers)*layer_thickness(nos_soil_layers))+water_change) &
                                            / (layer_thickness(nos_soil_layers)+depth_change)

            ! explicitly update the soil profile if there has been rooting depth
            ! changes
            layer_thickness(1) = top_soil_depth
            layer_thickness(2) = root_reach - top_soil_depth
            layer_thickness(3) = max_depth - sum(layer_thickness(1:2))

            ! keep track of the previous rooting depth
            previous_depth = root_reach

        else

            ! keep track of the previous rooting depth
            previous_depth = previous_depth

        end if ! depth change

    else if (root_reach < depth_change .and. previous_depth > depth_change) then

        !!!!!!!!!!
        ! Model has explicitly contracted from the bucket layer
        !!!!!!!!!!

        ! In this circumstance we want to return the soil profile to it's
        ! default structure with a minimum sized third layer
        depth_change = previous_depth - depth_change

        ! determine how much water is lost from the old volume of soil
        water_change = soil_waterfrac(nos_root_layers) * depth_change
        ! now assign that new volume of water to the deep rooting layer
        soil_waterfrac(nos_soil_layers) = ((soil_waterfrac(nos_soil_layers)*layer_thickness(nos_soil_layers))+water_change) &
                                        / (layer_thickness(nos_soil_layers)+depth_change)

        ! explicitly update the soil profile if there has been rooting depth
        ! changes
        layer_thickness(1) = top_soil_depth
        layer_thickness(2) = min_layer
        layer_thickness(3) = max_depth - sum(layer_thickness(1:2))

        ! keep track of the previous rooting depth
        previous_depth = min_layer

    else ! root_reach > (top_soil_depth + min_layer)

        ! if we are outside of the range when we need to consider rooting depth changes keep track in case we move into a zone when we do
        previous_depth = previous_depth

    endif ! root reach beyond top layer

    ! Update soil water potential
    call soil_water_potential

!    ! check water balance
!    balance = (rainfall_in - corrected_ET - underflow - runoff) * days_per_step
!    balance = balance &
!            - (sum(soil_waterfrac(1:nos_soil_layers) * layer_thickness(1:nos_soil_layers) * 1d3) &
!            - initial_soilwater)
!
!    if (abs(balance) > 1d-6 .or. soil_waterfrac(1) < -1d-6) then
!        print*,"Soil water miss-balance (mm)",balance
!        print*,"Initial_soilwater (mm) = ",initial_soilwater
!        print*,"Final_soilwater (mm) = ",sum(soil_waterfrac(1:nos_soil_layers) * layer_thickness(1:nos_soil_layers) * 1d3)
!        print*,"State balance = ",sum(soil_waterfrac(1:nos_soil_layers)*layer_thickness(1:nos_soil_layers)*1d3)-initial_soilwater
!        print*,"Flux balance = ",(rainfall_in - corrected_ET - underflow - runoff) * days_per_step
!        print*,"Top soilwater (fraction)",soil_waterfrac(1)
!        print*,"Rainfall (mm/step)",rainfall_in,"ET",corrected_ET,"underflow",underflow,"runoff",runoff
!        print*,"Rainfall (kgH2O/m2/s)",rainfall
!        print*,"Soil Water Fraction = ",soil_waterfrac
!    end if ! abs(balance) > 1d-10

    ! explicit return needed to ensure that function runs all needed code
    return

  end subroutine calculate_update_soil_water  
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_wetcanopy_evaporation(wetcanopy_evap,storage)

    ! Estimates evaporation of canopy intercepted rainfall based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(inout) :: storage      ! canopy water storage kgH2O/m2
    double precision, intent(out) :: wetcanopy_evap ! kgH2O.m-2.day-1

    ! local variables
    double precision :: canopy_radiation, & ! isothermal net radiation (W/m2)
                                      gb    ! stomatal and boundary layer conductance (m.s-1)

    ! Assuming there is any rainfall, currently water on the canopy or dew formation
    if (rainfall > 0d0 .or. storage > 0d0) then

        !!!!!!!!!!
        ! Calculate canopy conductance (to water vapour)
        !!!!!!!!!!

        ! Combine in series stomatal conductance with boundary layer
        gb = aerodynamic_conductance * leaf_canopy_wind_scaling

        !!!!!!!!!!
        ! Estimate energy radiation balance (W.m-2)
        !!!!!!!!!!

        ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
        canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * dayl_seconds_1)

        !!!!!!!!!!
        ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
        !!!!!!!!!!

        ! Calculate potential Penman Montheith (kgH2O.m-2.day-1)
        wetcanopy_evap = max(0d0,(((slope*canopy_radiation) + (ET_demand_coef*gb)) &
                                 / (lambda*(slope+psych))) * dayl_seconds)

        ! Update based on canopy water storage
        call canopy_interception_and_storage(wetcanopy_evap,storage)

    else

        ! there is no water movement possible
        intercepted_rainfall = 0d0 ; wetcanopy_evap = 0d0

    endif

  end subroutine calculate_wetcanopy_evaporation
  !
  !-----------------------------------------------------------------
  !
  subroutine canopy_interception_and_storage(potential_evaporation,storage)

    ! Simple daily time step integration of canopy rainfall interception, runoff
    ! and rainfall (kgH2O.m-2.s-1). NOTE: it is possible for intercepted rainfall to be
    ! negative if stored water running off into the soil is greater than
    ! rainfall (i.e. when leaves have died between steps)

    implicit none

    ! arguments
    double precision, intent(inout) :: storage, & ! canopy water storage (kgH2O/m2)
                         potential_evaporation    ! wet canopy evaporation (kgH2O.m-2.day-1),
                                                  ! enters as potential but leaves as water balance adjusted.
                                                  ! Note that this assumes a completely wet leaf surface
    ! local variables
    double precision :: a, through_fall, max_storage, max_storage_1, daily_addition, wetcanopy_evaporation &
                       ,potential_drainage_rate ,drain_rate, evap_rate, initial_canopy, co_mass_balance, dx, dz, tmp(3)
    ! local parameters
    double precision, parameter :: CanIntFrac = -0.5d0,     & ! Coefficient scaling rainfall interception fraction with LAI
                                        clump = 0.75d0,     & ! Clumping factor (1 = uniform, 0 totally clumped)
                                                              ! He et al., (2012) http://dx.doi.org/10.1016/j.rse.2011.12.008
                                  CanStorFrac = 0.2d0,      & ! Coefficient scaling canopy water storage with LAI
                                 RefDrainRate = 0.002d0,    & ! Reference drainage rate (mm/min; Rutter et al 1975)
                                  RefDrainLAI = 0.952381d0, & ! Reference drainage 1/LAI (m2/m2; Rutter et al 1975, 1/1.05)
                                 RefDrainCoef = 3.7d0,      & ! Reference drainage Coefficient (Rutter et al 1975)
                               RefDrainCoef_1 = RefDrainCoef ** (-1d0)

    ! hold initial canopy storage in memory
    initial_canopy = storage
    ! determine maximum canopy storage & through fall fraction
    through_fall = exp(CanIntFrac*lai*clump)
    ! maximum canopy storage (mm); minimum is applied to prevent errors in
    ! drainage calculation. Assume minimum capacity due to wood.
    max_storage = max(min_storage,CanStorFrac*lai)
    ! caclulate inverse for efficient calculations below
    max_storage_1 = max_storage**(-1d0)
    ! potential intercepted rainfall (kgH2O.m-2.s-1)
    intercepted_rainfall = rainfall * (1d0 - through_fall)

    ! calculate drainage coefficients (Rutter et al 1975); Corsican Pine
    ! 0.002 is canopy specific coefficient modified by 0.002*(max_storage/1.05)
    ! where max_storage is the canopy maximum capacity (mm) (LAI based) and
    ! 1.05 is the original canopy capacitance
    a = log( RefDrainRate * ( max_storage * RefDrainLAI ) ) - RefDrainCoef * max_storage

    ! average rainfall intercepted by canopy (kgH2O.m-2.day-1)
    daily_addition = intercepted_rainfall * seconds_per_day

    ! reset cumulative variables
    through_fall = 0d0 ; wetcanopy_evaporation = 0d0
    drain_rate = 0d0 ; evap_rate = 0d0

    ! add rain to the canopy and overflow as needed
    storage = storage + daily_addition

    if (storage > max_storage) then

        if (potential_evaporation > 0d0) then

            ! Assume water drainage will always occur an order magnitude above evaporation
            ! so water above canopy capacity is drained. Water below max_storage is accessable by evaporation only.

            ! Assume drainage is all water above the maximum canopy storage (kg/m2/day)
            drain_rate = storage - max_storage

            ! Estimate evaporation from remaining water (i.e. that left after
            ! initial co-access of evaporation and drainage).
            ! Assume evaporation is now restricted by:
            ! 1) energy already spent on evaporation (the -evap_rate) and
            ! 2) linear increase in surface resistance as the leaf surface
            ! dries (i.e. the 0.5).
            evap_rate = min(potential_evaporation * 0.5d0 * storage * max_storage_1, storage - drain_rate)

        else

            ! Load dew formation to the current local evap_rate variable
            evap_rate = potential_evaporation
            ! Restrict drainage the quantity above max_storage, adding dew formation too
            drain_rate = (storage - evap_rate) - max_storage

        endif

    else

        ! no drainage just apply evaporation / dew formation fluxes directly
        drain_rate = 0d0 ; evap_rate = potential_evaporation
        if (evap_rate > 0d0) then
            ! evaporation restricted by fraction of surface actually covered
            ! in water and integrated over period to bare leaf (i.e. the *0.5)
            evap_rate = evap_rate * storage * max_storage_1 * 0.5d0
            ! and the total amount of water
            evap_rate = min(evap_rate,storage)
        else
            ! then dew formation has occurred, if this pushes storage > max_storage add it to drainage
            drain_rate = max(0d0,(storage - evap_rate) - max_storage)
        endif ! evap_rate > 0

    endif ! storage > max_storage

    ! update canopy storage with water flux
    storage = storage - evap_rate - drain_rate
    wetcanopy_evaporation = wetcanopy_evaporation + evap_rate
    through_fall = through_fall + drain_rate

    ! correct intercepted rainfall rate to kgH2O.m-2.s-1
    intercepted_rainfall = intercepted_rainfall - (through_fall * seconds_per_day_1)

!    ! sanity checks; note 1e-8 prevents precision errors causing flags
!    if (intercepted_rainfall > rainfall .or. storage < -1d-8 .or. &
!       (wetcanopy_evaporation * days_per_step_1) > (1d-8 + initial_canopy + (rainfall*seconds_per_day)) ) then
!        print*,"Condition 1",intercepted_rainfall > rainfall
!        print*,"Condition 2",storage < -1d-8
!        print*,"Condition 3",(wetcanopy_evaporation * days_per_step_1) > (1d-8 + initial_canopy + (rainfall*seconds_per_day))
!        print*,"storage (kgH2O/m2)",storage,"max_storage (kgH2O/m2)",max_storage,"initial storage (kgH2O/m2)", initial_canopy
!        print*,"rainfall (kgH2O/m2/day)", rainfall*seconds_per_day, "through_fall (kgH2O/m2/day)", (through_fall * days_per_step_1)
!        print*,"through_fall_total (kgH2O/m2/step)",through_fall
!        print*,"potential_evaporation (kgH2O/m2/day)",potential_evaporation
!        print*,"actual evaporation    (kgH2O/m2/day)",wetcanopy_evaporation * days_per_step_1
!        stop
!    endif

    ! average evaporative flux to daily rate (kgH2O/m2/day)
    potential_evaporation = wetcanopy_evaporation

    ! final clearance of canopy storage of version small values at the level of system precision
    if (storage < 10d0*vsmall) storage = 0d0

  end subroutine canopy_interception_and_storage
  !
  !----------------------------------------------------------------------
  !
  double precision function find_gs_iWUE(gs_in)

    ! Calculate CO2 limited photosynthesis as a function of metabolic limited
    ! photosynthesis (pn), atmospheric CO2 concentration and stomatal
    ! conductance (gs_in). Photosynthesis is calculated twice to allow for
    ! testing of senstivity to iWUE.

    ! arguments
    double precision, intent(in) :: gs_in

    !!!!!!!!!!
    ! Optimise intrinsic water use efficiency
    !!!!!!!!!!

    ! Determine impact of gs increment on pd and how far we are from iWUE
    find_gs_iWUE = iWUE_step - (acm_gpp_stage_2(gs_in + delta_gs) - acm_gpp_stage_2(gs_in))

    ! Remember to return back to the user
    return

  end function find_gs_iWUE
  !
  !------------------------------------------------------------------
  !
  subroutine grass_cutting(labile,foliage,roots,litter,som,doy,lai_reduction, &
                           harvest,HARVESTextracted_labile,                   & 
                           HARVESTextracted_foliage,HARVESTextracted_roots,   &
                           HARVESTlitter_labile,HARVESTlitter_foliage,        &
                           HARVESTlitter_roots,                               &
                           timestep,nodays,step_length,                       &
                           cutting_threshold,post_cutting_labile_loss)
  
    !! Determine whether cutting has occured and the resulting impacts on the C-cycle
    
    implicit none
    
    ! Arguments
    integer, intent(in) :: timestep, nodays
    double precision, intent(in) :: doy, &
                            step_length, &
                          lai_reduction, &
                      cutting_threshold, &
               post_cutting_labile_loss
    
    double precision, intent(inout) :: labile, &
                                      foliage, &
                                        roots, &
                                       litter, & 
                                          som, &
                      HARVESTextracted_labile, &
                     HARVESTextracted_foliage, &
                       HARVESTextracted_roots, &
                         HARVESTlitter_labile, &
                        HARVESTlitter_foliage, &
                          HARVESTlitter_roots   
                          
    double precision, dimension(nodays), intent(inout) :: harvest
                           
    ! Local variables
    double precision :: labile_loss    & 
                       ,foliar_loss    &
                       ,roots_loss     &
                       ,labile_residue &
                       ,foliar_residue &
                       ,roots_residue  
           
    ! Determine whether cutting is plausible
    ! 1) Labile+leaf C > cutting threshold 
    ! 2) LAI > 3 (note replaced here with day of year constraints?)
    ! 3) & LAI reduction = -1 & no cut in past month 
    ! TLS: Question, conditions here are very temperate centric, can these be modifed?
    ! TLS: The timestepping of the cutting assumptions need to be dynamics in code to timestep
    if ( (labile+foliage) >= cutting_threshold   & 
          .and. doy >= 91d0 .and. doy <= 304d0 & 
          ! .and. LAI(n) >= 3 & 
          .and. sum(harvest(max(1,timestep-four_week_lag):timestep)) == 0d0 ) then
                  
        ! direct C losses
        labile_loss  = labile * post_cutting_labile_loss
        foliar_loss  = foliage * 0.95d0 ! 95% of leaves lost after cutting probably 99% lost in reality 
        roots_loss   = 0d0 ! POOLS(n+1,3) * roots_frac_death ! allocation to roots will be reduced due to reduced LAI 

        ! fraction of harvest wasted 
        labile_residue = labile_loss * labile_frac_res
        foliar_residue = foliar_loss * foliage_frac_res
        roots_residue  = roots_loss  * roots_frac_res

        ! if havest yields > 1500 kg.DM.ha-1 proceed with cut
        ! Note converted to gC/m2 equivalent assuming 47.5 % C content
        ! yields 71.25 gC/m2
        if ( ( (foliar_loss-foliar_residue)+ &
               (labile_loss-labile_residue)+ &
               (roots_loss -roots_residue ) ) >= 71.25d0 ) then
                      
            ! Assign to output the biomass extracted due to cutting
            HARVESTextracted_labile  = (labile_loss-labile_residue) / step_length
            HARVESTextracted_foliage = (foliar_loss-foliar_residue) / step_length
            HARVESTextracted_roots   = (roots_loss -roots_residue) / step_length
            ! Assign the output the biomass entering litter due to cutting
            HARVESTlitter_labile  = labile_residue / step_length
            HARVESTlitter_foliage = foliar_residue / step_length
            HARVESTlitter_roots   = roots_residue / step_length

            ! Combine to total extracted based on cutting
            harvest(timestep) = HARVESTextracted_labile  + &
                                HARVESTextracted_foliage + & 
                                HARVESTextracted_roots 

            ! update pools 
            labile  = max(0d0,labile-labile_loss)
            foliage = max(0d0,foliage-foliar_loss)
            roots   = max(0d0,roots-roots_loss)
            litter  = max(0d0,litter + (labile_residue+foliar_residue+roots_residue))
            som     = max(0d0,som)

        endif ! Determine whether the cut is actually plausible

    endif ! end cutting process   

  end subroutine grass_cutting
  !
  !------------------------------------------------------------------
  !
  subroutine grass_grazing(labile,foliage,roots,litter,som,doy,lai_reduction,       &
                           animal_manure_to_soil,animal_respiration,animal_methane, &
                           harvest,grazing,GRAZINGextracted_labile,                 & 
                           GRAZINGextracted_foliage,GRAZINGextracted_roots,         &
                           GRAZINGlitter_labile,GRAZINGlitter_foliage,              &
                           GRAZINGlitter_roots,                                     &
                           timestep,nodays,step_length,                             &
                           lca,grazing_threshold,post_grazing_labile_loss,          &
                           min_grazing_removal_threshold)
  
    !! Determine whether grazing has occured and the resulting impacts on the C-cycle
    
    implicit none
    
    ! Arguments
    integer, intent(in) :: timestep, nodays
    double precision, intent(in) :: doy, &
                            step_length, &
                          lai_reduction, &
                                    lca, &
                      grazing_threshold, &
          min_grazing_removal_threshold, &
               post_grazing_labile_loss
    
    double precision, intent(inout) :: labile, &
                                      foliage, &
                                        roots, &
                                       litter, & 
                                          som, & 
                      GRAZINGextracted_labile, &
                     GRAZINGextracted_foliage, &
                       GRAZINGextracted_roots, &
                         GRAZINGlitter_labile, &
                        GRAZINGlitter_foliage, &
                          GRAZINGlitter_roots, &
                        animal_manure_to_soil, &
                           animal_respiration, &
                               animal_methane
                       
    double precision, dimension(nodays), intent(inout) :: harvest, &
                                                          grazing

    ! Local variables
    double precision :: labile_loss    & 
                       ,foliar_loss    &
                       ,roots_loss     &
                       ,labile_residue &
                       ,foliar_residue &
                       ,roots_residue 
        
    ! Determine whether grazing is plausible
    ! 1) An LAI reduction is specified
    ! 2) Labile+leaf C > grazing threshold 
    ! 3) No cutting in the last 2 weeks
    if (labile+foliage >= grazing_threshold .and. & 
        sum(harvest(max(1,timestep-two_week_lag):timestep)) == 0d0) then
           
        ! direct C losses
        labile_loss  = labile * post_grazing_labile_loss
        foliar_loss  = max(0d0,(lai_reduction * lca) - labile_loss)  
        roots_loss   = 0d0 ! POOLS(n+1,3) * roots_frac_death

        ! fraction of harvest wasted 
        labile_residue = labile_loss * labile_frac_res
        foliar_residue = foliar_loss * foliage_frac_res
        roots_residue = roots_loss * roots_frac_res

        ! extracted C via grazing: if remaining AGB > pre-grazing limit DM & grazed biomass > pars(34) gC/m2/step
        if ( (labile+foliage)-foliar_loss-labile_loss >= grazing_threshold .and. & 
               (foliar_loss+labile_loss) >= min_grazing_removal_threshold*step_length ) then

            ! Assign to output the biomass extracted due to cutting
            GRAZINGextracted_labile = (labile_loss-labile_residue) / step_length
            GRAZINGextracted_foliage = (foliar_loss-foliar_residue) / step_length
            GRAZINGextracted_roots = (roots_loss-roots_residue) / step_length
            ! Assign the output the biomass entering litter due to cutting
            GRAZINGlitter_labile = labile_residue / step_length
            GRAZINGlitter_foliage = foliar_residue / step_length
            GRAZINGlitter_roots = roots_residue / step_length

            ! Combine to total extracted based on grazing
            grazing(timestep) = GRAZINGextracted_labile  + &
                                GRAZINGextracted_foliage + &
                                GRAZINGextracted_roots
               
            ! Constants used for animal C fluxes from Vertes.et.al.2019 (10.1016/B978-0-12-811050-8.00002-9)

            ! animal manure-C production (gC/m2)
            animal_manure_to_soil = (grazing(timestep) * 0.32d0) / step_length
            ! animal respiration CO2-C (gC/m2)
            animal_respiration = (grazing(timestep) * 0.54d0) / step_length
            ! animal CH4-C (gC/m2)
            animal_methane = (grazing(timestep) * 0.04d0) / step_length                

            ! update pools 
            labile  = max(0d0,labile-labile_loss)
            foliage = max(0d0,foliage-foliar_loss)
            roots   = max(0d0,roots-roots_loss)
            litter  = max(0d0,litter+((animal_manure_to_soil+GRAZINGlitter_labile + &
                                       GRAZINGlitter_foliage+GRAZINGlitter_roots) * step_length))
            som     = max(0d0,som)

        endif 

        ! ! Determine whether to graze, assuming grazing criterior above was not met 
        ! ! extracted C via grazing: if not done above & postgraze AGB < pre-grazing AGB 
        ! if ( (grazing(timestep) == 0d0) .and. & 
        !      ((labile+foliage)-foliar_loss-labile_loss) <= grazing_threshold ) then
              
        !     ! Bulk C losses
        !     labile_loss  = labile * post_grazing_labile_loss
        !     foliar_loss  = foliage - grazing_threshold - (labile - labile_loss)
        !     !foliar_loss  = foliage - (grazing_threshold - labile_loss)
        !     roots_loss   = 0d0 ! roots * roots_frac_death

        !     ! fraction of harvest wasted 
        !     labile_residue = labile_loss * labile_frac_res
        !     foliar_residue = foliar_loss * foliage_frac_res
        !     roots_residue = roots_loss * roots_frac_res

        !     ! proceed if simulating this grazing will remove > ~0.5 gCm-2 from AGB
        !     if ((foliar_loss+labile_loss) >= min_grazing_removal_threshold) then

        !         ! Assign to output the biomass extracted due to cutting
        !         GRAZINGextracted_labile  = (labile_loss-labile_residue) / step_length
        !         GRAZINGextracted_foliage = (foliar_loss-foliar_residue) / step_length
        !         GRAZINGextracted_roots   = (roots_loss -roots_residue) / step_length
        !         ! Assign the output the biomass entering litter due to cutting
        !         GRAZINGlitter_labile  = labile_residue / step_length
        !         GRAZINGlitter_foliage = foliar_residue / step_length
        !         GRAZINGlitter_roots   = roots_residue / step_length

        !         ! Combine to total extracted based on grazing
        !         grazing(timestep) = GRAZINGextracted_labile + &
        !                             GRAZINGextracted_foliage + &
        !                             GRAZINGextracted_roots

        !         ! animal manure-C production (gC/m2)
        !         animal_manure_to_soil = (grazing(timestep) * 0.32d0) / step_length
        !         ! animal respiration CO2-C (gC/m2)
        !         animal_respiration = (grazing(timestep) * 0.54d0) / step_length
        !         ! animal CH4-C (gC/m2)
        !         animal_methane = (grazing(timestep) * 0.04d0) / step_length                

        !         ! update pools 
        !         labile  = max(0d0,labile-labile_loss)
        !         foliage = max(0d0,foliage-foliar_loss)
        !         roots   = max(0d0,roots-roots_loss)
        !         litter  = max(0d0,litter+((animal_manure_to_soil+GRAZINGlitter_labile + &
        !                                    GRAZINGlitter_foliage+GRAZINGlitter_roots) * step_length))
        !         som     = max(0d0,som)

        !     endif ! carry out grazing?

        ! endif ! Is grazing plausible

    endif ! end grazing process

  end subroutine grass_grazing
  !
  !-----------------------------------------------------------------
  !
  subroutine gravitational_drainage(time_period_days)

    ! Integrator for soil gravitational drainage.
    ! Due to the longer time steps undertake by ACM / DALEC and the fact that
    ! drainage is a concurrent processes we assume that drainage occurs at
    ! the bottom of the column first creating space into which water can drain
    ! from the top down. Therefore we draing from the bottom first and then the top.
    ! NOTE: Assumes that any previous water movement due to infiltration and evaporation
    !       has already been updated in soil mass balance

    implicit none

    ! arguments
    integer, intent(in) :: time_period_days

    ! local variables..
    integer :: t
    double precision, dimension(nos_soil_layers) :: dx, & ! range between the start and end points of the integration
                                               halfway, & ! half way point between start and end point of integration
                                                liquid, & ! liquid water in local soil layer (m3/m3)
                                         avail_to_flow, & ! liquid content above field capacity (m3/m3)
                                               iceprop, & ! fraction of soil layer which is ice
                                          pot_drainage    ! estimats of time step potential drainage rate (m/s)
    double precision  :: tmp1,tmp2,tmp3 &
                                 ,unsat & ! unsaturated pore space in soil_layer below the current (m3/m3)
                                ,change   ! absolute volume of water drainage in current layer (m3/day)

    ! calculate soil ice proportion; at the moment
    ! assume everything liquid
    iceprop = 0d0

    ! except the surface layer in the mean daily temperature is < 0oC
    if (meant < 1d0) iceprop(1) = 1d0

    ! zero water fluxes
    waterchange = 0d0

    ! underflow is tracked in kgH2O/m2/day but estimated here in MgH2O/m2/day
    ! therefore we must convert
    underflow = underflow * 1d-3

    ! estimate potential drainage rate for the current time period
    liquid = soil_waterfrac(1:nos_soil_layers) * ( 1d0 - iceprop(1:nos_soil_layers) )
    ! estimate how much liquid is available to flow
    avail_to_flow = liquid - field_capacity(1:nos_soil_layers)
    ! trapezium rule scaler and the half-way point between current and field capacity
    dx = avail_to_flow*0.5d0 ; halfway = liquid - dx
    do t = 1, nos_soil_layers
       if (avail_to_flow(t) > 0d0) then
           ! Trapezium rule for approximating integral of drainage rate
           call calculate_soil_conductivity(t,liquid(t),tmp1)
           call calculate_soil_conductivity(t,field_capacity(t),tmp2)
           call calculate_soil_conductivity(t,halfway(t),tmp3)
           pot_drainage(t) = 0.5d0 * dx(t) * ((tmp1 + tmp2) + 2d0 * tmp3)
       else
           ! We are at field capacity currently even after rainfall has been infiltrated.
           ! Assume that the potential drainage rate is that at field capacity
           call calculate_soil_conductivity(t,field_capacity(t),pot_drainage(t))
       endif ! water above field capacity to flow?
    end do ! soil layers
    ! Scale potential drainage from per second to per day
    pot_drainage = pot_drainage * seconds_per_day

    ! Integrate drainage over each day until time period has been reached or
    ! each soil layer has reached field capacity
    t = 1
    do while (t < (time_period_days+1) .and. maxval(soil_waterfrac - field_capacity) > vsmall)

       ! Estimate liquid content and how much is available to flow / drain
       avail_to_flow = ( soil_waterfrac(1:nos_soil_layers) * (1d0 - iceprop(1:nos_soil_layers)) ) &
                     - field_capacity(1:nos_soil_layers)

       ! ...then from the top down
       do soil_layer = 1, nos_soil_layers

          ! initial conditions; i.e. is there liquid water and more water than
          ! layer can hold
          if (avail_to_flow(soil_layer) > 0d0 .and. soil_waterfrac(soil_layer+1) < porosity(soil_layer+1)) then

              ! Unsaturated volume of layer below (m3 m-2)
              unsat = ( porosity(soil_layer+1) - soil_waterfrac(soil_layer+1) ) &
                    * layer_thickness(soil_layer+1) / layer_thickness(soil_layer)
              ! Restrict potential rate calculate above for the available water
              ! and available space in the layer below.
              ! NOTE: * layer_thickness(soil_layer) converts units from m3/m2 -> (m3)
              change = min(unsat,min(pot_drainage(soil_layer),avail_to_flow(soil_layer))) * layer_thickness(soil_layer)
              ! update soil layer below with drained liquid
              waterchange( soil_layer + 1 ) = waterchange( soil_layer + 1 ) + change
              waterchange( soil_layer     ) = waterchange( soil_layer     ) - change

          end if ! some liquid water and drainage possible

       end do ! soil layers

       ! update soil water profile
       soil_waterfrac(1:nos_soil_layers) = soil_waterfrac(1:nos_soil_layers) &
                                         + (waterchange(1:nos_soil_layers)/layer_thickness(1:nos_soil_layers))
       ! estimate drainage from bottom of soil column (MgH2O/m2/day)
       ! NOTES: that underflow is reset outside of the daily soil loop
       underflow = underflow + waterchange(nos_soil_layers+1)

       ! Reset now we have moves that liquid
       waterchange = 0d0
       ! integerate through time period
       t = t + 1

    end do ! while condition

    ! convert underflow from MgH2O/m2/day -> kgH2O/m2/day
    underflow = underflow * 1d3

  end subroutine gravitational_drainage

  !
  !------------------------------------------------------------------
  !
  subroutine gsi_phenology(nodays,timestep,step_length,gsi_lag, &
                           maxt_21day,dayl_21day,vpd_21day,     &
                           potential_foliar_turnover,potential_labile_turnover, &
                           lca,critical_gpp_return, &
                           Tfac_min,Tfac_max,Pfac_min,Pfac_max, & 
                           VPDfac_min,VPDfac_max,               &
                           current_gpp,gsi,     & 
                           Tfac,Pfac,VPDfac,    &
                           foliar_turnover_fraction, &
                           labile_turnover_fraction, &
                           labile,foliage)
 
    !! Description
  
    implicit none
  
    ! Declare arguments
    integer, intent(in) :: nodays, gsi_lag, timestep
    double precision, intent(in) :: step_length, &
                                     maxt_21day, & 
                                     dayl_21day, & 
                                      vpd_21day, &
                      potential_foliar_turnover, & 
                      potential_labile_turnover, & 
                              Tfac_min,Tfac_max, &
                              Pfac_min,Pfac_max, & 
                          VPDfac_min,VPDfac_max, & 
                                            lca, &
                                    current_gpp, &
                            critical_gpp_return, &
                                         labile, & 
                                        foliage
    double precision, dimension(nodays), intent(inout) :: gsi
    double precision, intent(out) :: foliar_turnover_fraction, & 
                                     labile_turnover_fraction, &
                                                         Tfac, & 
                                                         Pfac, & 
                                                       VPDfac
  
    ! Declare local arguments
    integer :: interval
    double precision :: tmp, gradient, &
                       canopy_lw_save, &
                         soil_lw_save, &
                       canopy_sw_save, &
                      canopy_par_save, &
                         soil_sw_save, &
                              gs_save, &
                              ga_save, &
                             lai_save

    ! Store variables for restoration at the end
    ! save original values for re-allocation later
    !canopy_lw_save = canopy_lwrad_Wm2 ; soil_lw_save = soil_lwrad_Wm2
    canopy_sw_save = canopy_swrad_MJday ; canopy_par_save  = canopy_par_MJday
    soil_sw_save = soil_swrad_MJday ; gs_save = stomatal_conductance
    ga_save = aerodynamic_conductance ; lai_save = lai

    ! Calculate the Growing Season Index based on Jolly et al. 
    ! doi: 10.1111/j.1365-2486.2005.00930.x doi:10.1029/2010JG001545.
    ! It is the product of 3 limiting factors for temperature, photoperiod and
    ! vapour pressure deficit that grow linearly from 0 to 1 between a calibrated 
    ! min and max value. Photoperiod, VPD and avgTmin are direct input

    ! temperature limitation, then restrict to 0-1; correction for k-> oC
    ! Tfac = (met(10,n)-(pars(12)-273.15)) / (pars(13)-pars(12)) ! no need to K->C 
    Tfac = ( maxt_21day-(Tfac_min-273.15d0)) / (Tfac_max-Tfac_min )
    Tfac = min(1d0,max(0d0,Tfac))
    ! photoperiod limitation (seconds)
    Pfac = ( dayl_21day-Pfac_min) / (Pfac_max-Pfac_min )
    Pfac = min(1d0,max(0d0,Pfac))
    ! VPD limitation (Pa)
    VPDfac = 1d0 - ( (vpd_21day-VPDfac_min) / (VPDfac_max-VPDfac_min) )
    VPDfac = min(1d0,max(0d0,VPDfac))
    ! calculate and store the GSI index
    gsi(timestep) = Tfac * Pfac * VPDfac
  
    ! Determine GSI section to have linear regression applied to and
    ! determine the number of values, i.e. the interval
    if (timestep < gsi_lag) then
        if (timestep == 1) then
            gsi_history(2) = gsi(timestep)
            interval = 2
        else
            gsi_history(1:timestep) = gsi(1:timestep)
            interval = timestep
        endif
    else
        gsi_history(1:gsi_lag) = gsi((timestep-gsi_lag+1):timestep)
        interval = gsi_lag
    end if
 
    ! Now calculate the linear gradient
    gradient = linear_model_gradient(tmp_x(1:interval),gsi_history(1:interval),interval)      
    gsi_lag_remembered = gsi_lag

    ! now update foliage and labile conditions based on gradient calculations
    if (gradient < fol_turn_crit .or. gsi(timestep) == 0d0) then
        ! we are in a decending condition so foliar turnover
        foliar_turnover_fraction = potential_foliar_turnover*(1d0-gsi(timestep))
        just_grown = 0.5d0
    else if (gradient > lab_turn_crit) then
        ! we are in a assending condition so labile turnover
        labile_turnover_fraction = potential_labile_turnover*gsi(timestep)
        just_grown = 1.5d0
        ! check carbon return
        tmp = labile*(1d0-(1d0-labile_turnover_fraction)**step_length)/step_length
        ! Update the canopy
        lai = (foliage+tmp)/lca
        tmp = lai / lai_save
        ! And estimate the canopy conditions for the new leaf area under the
        ! current time step
        aerodynamic_conductance = aerodynamic_conductance * tmp
        stomatal_conductance = stomatal_conductance * tmp
        call calculate_shortwave_balance
        if (lai_save < vsmall) then
            call calculate_aerodynamic_conductance
            call calculate_stomatal_conductance
        endif
        ! And estimate potential GPP
        if (stomatal_conductance > vsmall) then
             call acm_gpp_stage_1
            tmp = acm_gpp_stage_2(stomatal_conductance)
        else
            tmp = 0d0
        endif
        ! determine if increase in LAI leads to an improvement in GPP greater
        ! than critical value, if not then no labile turnover allowed      
        if ( ((tmp - current_gpp)/current_gpp) < critical_gpp_return ) then
            labile_turnover_fraction = 0d0
        endif
    else
        ! probably we want nothing to happen, however if we are at the seasonal
        ! maximum we will consider further growth still
        if (just_grown >= 1d0) then
            ! we are between so definitely not losing foliage and we have
            ! previously been growing so maybe we still have a marginal return on
            ! doing so again
            labile_turnover_fraction = potential_labile_turnover*gsi(timestep)
            ! but possibly gaining some?
            ! determine if this is a good idea based on GPP increment
            tmp = labile*(1d0-(1d0-labile_turnover_fraction)**step_length)/step_length
            ! Update the canopy
            lai = (foliage+tmp)/lca
            tmp = lai / lai_save
            ! And estimate the canopy conditions for the new leaf area under the
            ! current time step
            aerodynamic_conductance = aerodynamic_conductance * tmp
            stomatal_conductance = stomatal_conductance * tmp
            call calculate_shortwave_balance
            if (lai_save < vsmall) then
                call calculate_aerodynamic_conductance
                call calculate_stomatal_conductance
            endif
            ! And estimate potential GPP
            if (stomatal_conductance > vsmall) then
                 call acm_gpp_stage_1
                tmp = acm_gpp_stage_2(stomatal_conductance)
            else
                    tmp = 0d0
            endif
            ! determine if increase in LAI leads to an improvement in GPP greater
            ! than critical value, if not then no labile turnover allowed
            if ( ((tmp - current_gpp)/current_gpp) < critical_gpp_return ) then
                labile_turnover_fraction = 0d0
            endif
        end if ! Just grown?
    endif ! gradient choice

    ! restore original value back from memory
    lai = lai_save
    !canopy_lwrad_Wm2 = canopy_lw_save ; soil_lwrad_Wm2 = soil_lw_save
    canopy_swrad_MJday = canopy_sw_save ; canopy_par_MJday = canopy_par_save
    soil_swrad_MJday = soil_sw_save ; stomatal_conductance = gs_save
    aerodynamic_conductance = ga_save
        
  end subroutine gsi_phenology
  !
  !-----------------------------------------------------------------
  !
  subroutine infiltrate(rainfall_in)

    ! Takes surface_watermm and distributes it among top !
    ! layers. Assumes total infilatration in timestep.   !
    ! NOTE: Assumes that any previous water movement due to infiltration and evaporation
    !       has already been updated in soil mass balance

    implicit none

    ! arguments
    double precision, intent(in) :: rainfall_in ! rainfall (kg.m-2.day-1)

    ! local argumemts
    integer :: i
    double precision :: add, & ! surface water available for infiltration (m)
                      wdiff    ! available space in a given soil layer for water to fill (m)

    ! convert rainfall water from mm -> m (or kgH2O.m-2.day-1 -> MgH2O.m-2.day-1)
    add = rainfall_in * 1d-3

    do i = 1 , nos_soil_layers

       ! is the input of water greater than available space
       ! if so fill and subtract from input and move on to the next
       ! layer determine the available pore space in current soil layer
       wdiff = max(0d0,(porosity(i)-soil_waterfrac(i)) * layer_thickness(i))

       if (add > wdiff) then
           ! if so fill and subtract from input and move on to the next layer
           waterchange(i) = waterchange(i) + wdiff
           add = add - wdiff
       else
           ! otherwise infiltate all in the current layer
           waterchange(i) = waterchange(i) + add
           add = 0d0 ; exit
       end if

    end do ! nos_soil_layers

    ! if after all of this we have some water left assume it is runoff (kgH2O.m-2.day-1)
    ! NOTE that runoff is reset outside of the daily soil loop
    runoff = runoff + (add * 1d3)

  end subroutine infiltrate
  !
  !---------------------------------------------------------------------
  !
  subroutine initialise_soils(soil_frac_clay,soil_frac_sand)

    !
    ! Subroutine calculate the soil layers field capacities and sets the initial
    ! soil water potential set to field capacity
    !

    implicit none

    ! arguments
    double precision, dimension(3) :: soil_frac_clay &
                                                   ,soil_frac_sand

    ! local variables
    integer :: i

    ! Include some hardcoded boundaries for the Saxton equations
    ! NOTE: do loop was found to be faster than 'where' for small vectors
    do i = 1, nos_soil_layers
       if (soil_frac_sand(i) < 5d0) soil_frac_sand(i) = 5d0
       if (soil_frac_clay(i) < 5d0) soil_frac_clay(i) = 5d0
       !if (soil_frac_sand(i) > 60d0) soil_frac_sand(i) = 60d0
       if (soil_frac_sand(i) > 70d0) soil_frac_sand(i) = 70d0
    end do
    ! calculate soil porosity (m3/m3)
    call soil_porosity(soil_frac_clay,soil_frac_sand)
    ! calculate field capacity (m3/m-3)
    call calculate_field_capacity

    ! final sanity check for porosity
    do i = 1, nos_soil_layers+1
       if (porosity(i) < (field_capacity(i)+0.05d0)) porosity(i) = field_capacity(i) + 0.05d0
    end do

  end subroutine initialise_soils  
  !
  !------------------------------------------------------------------
  !
  subroutine log_law_decay

    ! Standard log-law above canopy wind speed (m.s-1) decay under neutral
    ! conditions.
    ! See Harman & Finnigan 2008; Jones 1992 etc for details.

    implicit none

    ! log law decay, NOTE: given canopy height (9 m) the log function reduces
    ! to a constant value down to ~ 7 decimal place (0.3161471806). Therefore
    ! 1/vonkarman * 0.31 = 0.7710906
    canopy_wind = ustar * vonkarman_1 * log((canopy_height-displacement) / roughl)

    ! set minimum value for wind speed at canopy top (m.s-1)
!    canopy_wind = max(min_wind,canopy_wind)

  end subroutine log_law_decay
  !
  !------------------------------------------------------------------
  !
  subroutine meteorological_constants(input_temperature,input_temperature_K,input_vpd_kPa)

    ! Determine some multiple use constants used by a wide range of functions
    ! All variables here are linked to air temperature and thus invarient between
    ! iterations and can be stored in memory...

    implicit none

    ! arguments
    double precision, intent(in) :: input_temperature, input_temperature_K, &
                                    input_vpd_kPa

    ! local variables
    double precision :: mult, &
              dynamic_viscosity    ! dynamic viscosity (kg.m-2.s-1)

    !
    ! Used for soil, canopy evaporation and transpiration
    !

    ! Density of air (kg/m3)
    air_density_kg = 353d0/input_temperature_K
    ! Conversion ratio for m.s-1 -> mol.m-2.s-1
    convert_ms1_mol_1 = const_sfc_pressure / (input_temperature_K*Rcon)
    ! latent heat of vapourisation,
    ! function of air temperature (J.kg-1)
    lambda = 2501000d0-2364d0*input_temperature

    ! psychrometric constant (kPa K-1)
    psych = (0.0646d0*exp(0.00097d0*input_temperature))
    ! Straight line approximation of the true slope; used in determining
    ! relationship slope
    mult = input_temperature+237.3d0
    ! 2502.935945 = 0.61078*17.269*237.3
    ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    slope = (2502.935945d0*exp(17.269d0*input_temperature/mult)) / (mult*mult)

    ! estimate frequently used atmmospheric demand component
    ET_demand_coef = air_density_kg*cpair*input_vpd_kPa

    !
    ! Used for soil evaporation and leaf level conductance
    !

    ! Determine diffusion coefficient (m2.s-1), temperature dependant (pressure dependence neglected). Jones p51; appendix 2
    ! Temperature adjusted from standard 20oC (293.15 K), NOTE that 1/293.15 = 0.003411223
    ! 0.0000242 = conversion to make diffusion specific for water vapor (um2.s-1)
    water_vapour_diffusion = 0.0000242d0*((input_temperature_K/293.15d0)**1.75d0)

    !
    ! Used for calculation of leaf level conductance
    !

    ! Calculate the dynamic viscosity of air (kg.m-2.s-1)
    dynamic_viscosity = ((input_temperature_K**1.5d0)/(input_temperature_K+120d0))*1.4963d-6
    ! and kinematic viscosity (m2.s-1)
    kinematic_viscosity = dynamic_viscosity/air_density_kg

  end subroutine meteorological_constants
  !
  !------------------------------------------------------------------
  !
  subroutine plant_soil_flow(root_layer,root_length,root_mass &
                            ,demand,root_reach_in,transpiration_resistance &
                            ,Rtot_layer)

    !
    ! Calculate soil layer specific water flow form the soil to canopy (mmolH2O.m-2.s-1)
    ! Accounting for soil, root and plant resistance, and canopy demand
    !

    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! From the current soil layer given an amount of root within the soil layer.

    implicit none

    ! arguments
    integer, intent(in) :: root_layer
    double precision, intent(in) :: root_length, &
                                      root_mass, &
                                         demand, &
                                  root_reach_in, &
                       transpiration_resistance
    double precision, intent(out) :: Rtot_layer

    ! local arguments
    double precision :: soilR1, soilR2

    ! Estimate soil hydraulic resistance to water flow (MPa m2 s mmol-1)
    ! Note: 1) soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head.
    !       2) soil resistance calculation in single line to reduce assignment costs
    soilR1 = ( log(root_radius_1*(root_length*pi)**(-0.5d0)) &
               /(two_pi*root_length*root_reach_in*(soil_conductivity(root_layer)*head_1))) &
             * 1d-9 * mol_to_g_water
    ! Calculates root hydraulic resistance (MPa m2 s mmol-1) in a soil-root zone
    soilR2 = root_resist / (root_mass*root_reach_in)
    ! Estimate the total hydraulic resistance for the layer
    Rtot_layer = transpiration_resistance + soilR1 + soilR2

    ! Estimate the soil to plant flow of water mmolH2O/m2/s
    water_flux_mmolH2Om2s(root_layer) = demand/Rtot_layer

    ! return
    return

  end subroutine plant_soil_flow
  !
  !------------------------------------------------------------------
  !
  subroutine saxton_parameters(soil_frac_clay,soil_frac_sand)

    ! Calculate the key parameters of the Saxton, that is cond1,2,3 !
    ! and potA,B                                                    !

    implicit none

    ! arguments
    double precision, dimension(3) :: soil_frac_clay &
                                                   ,soil_frac_sand

    ! local variables
    double precision, parameter :: A = -4.396d0,  B = -0.0715d0,   CC = -4.880d-4, D = -4.285d-5, &
                                   E = -3.140d0,  F = -2.22d-3,     G = -3.484d-5, H = 0.332d0,   &
                                   J = -7.251d-4, K = 0.1276d0,     P = 12.012d0,  Q = -7.551d-2, &
                                   R = -3.895d0,  T = 3.671d-2,     U = -0.1103d0, V = 8.7546d-4, &
                                   mult1 = 100d0, mult2 = 2.778d-6

    ! layed out in this manor to avoid memory management issues in module
    ! variables
    potA(1:nos_soil_layers) = A + (B * soil_frac_clay) + &
                             (CC * soil_frac_sand * soil_frac_sand) + &
                              (D * soil_frac_sand * soil_frac_sand * soil_frac_clay)
    potA(1:nos_soil_layers) = exp(potA(1:nos_soil_layers))
    potA(1:nos_soil_layers) = potA(1:nos_soil_layers) * mult1

    potB(1:nos_soil_layers) = E + (F * soil_frac_clay * soil_frac_clay) + &
                                  (G * soil_frac_sand * soil_frac_sand * soil_frac_clay)

    cond1(1:nos_soil_layers) = mult2
    cond2(1:nos_soil_layers) = P + (Q * soil_frac_sand)
    cond3(1:nos_soil_layers) = R + (T * soil_frac_sand) + (U * soil_frac_clay) + &
                                   (V * soil_frac_clay * soil_frac_clay)

    ! assign bottom of soil column value to core
    potA(nos_soil_layers+1)  = potA(nos_soil_layers)
    potB(nos_soil_layers+1)  = potB(nos_soil_layers)
    cond1(nos_soil_layers+1) = mult2
    cond2(nos_soil_layers+1) = cond2(nos_soil_layers)
    cond3(nos_soil_layers+1) = cond3(nos_soil_layers)

  end subroutine saxton_parameters
  !
  !-----------------------------------------------------------------
  !
  subroutine soil_porosity(soil_frac_clay,soil_frac_sand)

    ! Porosity is estimated from Saxton equations. !

    implicit none

    ! arguments
    double precision, dimension(3) :: soil_frac_clay &
                                                   ,soil_frac_sand
    ! local variables..
    double precision, parameter :: H = 0.332d0, &
                                 J = -7.251d-4, &
                                  K = 0.1276d0

    ! loop over soil layers..
    porosity(1:nos_soil_layers) = H + J * soil_frac_sand(1:nos_soil_layers) + &
                                  K * log10(soil_frac_clay(1:nos_soil_layers))
    ! then assign same to core layer
    porosity(nos_soil_layers+1) = porosity(nos_soil_layers)

  end subroutine soil_porosity
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_water_potential

    ! Find SWP without updating waterfrac yet (we do that in !
    ! waterthermal). Waterfrac is m3 m-3, soilwp is MPa.     !

    implicit none

    integer :: i

    ! reformulation aims to remove if statement within loop to hopefully improve
    ! optimisation
    SWP(1:nos_soil_layers) = -0.001d0 * potA(1:nos_soil_layers) &
                           * soil_waterfrac(1:nos_soil_layers)**potB(1:nos_soil_layers)
    ! NOTE: profiling indiates that 'where' is slower for very short vectors
    do i = 1, nos_soil_layers
       if (SWP(i) < -20d0 .or. SWP(i) /= SWP(i)) SWP(i) = -20d0
    end do

  end subroutine soil_water_potential  
  !
  !---------------------------------------------------------------------
  !
  subroutine update_soil_initial_conditions(input_soilwater_frac)

    !
    ! Subroutine calculate the soil layers field capacities and sets the initial
    ! soil water potential set to field capacity
    !

    implicit none

    ! arguments
    double precision, intent(in) :: input_soilwater_frac ! initial soil water status as fraction of field capacity

    ! local variables
    integer :: i

    ! Default assumption to be field capacity
    soil_waterfrac = field_capacity
    SWP = SWP_initial

    ! If prior value has been given
    if (input_soilwater_frac > -9998d0) then
        ! calculate initial soil water fraction
        soil_waterfrac(1:nos_soil_layers) = input_soilwater_frac * field_capacity(1:nos_soil_layers)
        ! calculate initial soil water potential
        call soil_water_potential
    endif

    ! Seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
       call calculate_soil_conductivity(i,soil_waterfrac(i),soil_conductivity(i))
    end do ! soil layers
    ! but apply the lowest soil layer to the core as well in initial conditions
    soil_conductivity(nos_soil_layers+1) = soil_conductivity(nos_soil_layers)

  end subroutine update_soil_initial_conditions  
  !
  !------------------------------------------------------------------
  !
  subroutine z0_displacement(ustar_Uh,local_lai)

    ! dynamic calculation of roughness length and zero place displacement (m)
    ! based on canopy height and lai. Raupach (1994)

    implicit none

    ! arguments
    double precision, intent(out) :: ustar_Uh ! ratio of friction velocity over wind speed at canopy top
    double precision, intent(in) :: local_lai
    ! local variables
    double precision  sqrt_cd1_lai
    double precision, parameter :: cd1 = 7.5d0,   & ! Canopy drag parameter; fitted to data
                                    Cs = 0.003d0, & ! Substrate drag coefficient
                                    Cr = 0.3d0,   & ! Roughness element drag coefficient
                          ustar_Uh_max = 0.35d0,  &
                          ustar_Uh_min = 0.05d0,  &
                                    Cw = 2d0,     &  ! Characterises roughness sublayer depth (m)
                                 phi_h = 0.19314718056d0 ! Roughness sublayer influence function;

    ! describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law

    ! Estimate canopy drag coefficient
    sqrt_cd1_lai = sqrt(cd1 * local_lai)

    ! calculate estimate of ratio of friction velocity / canopy wind speed.
    ! NOTE: under current min LAI and fixed canopy height (9 m) this ratio is
    ! fixed at 0.3
    ustar_Uh = 0.3d0
!    ustar_Uh = max(ustar_Uh_min,min(sqrt(Cs+Cr*local_lai*0.5d0),ustar_Uh_max))
!    ustar_Uh = sqrt(Cs+Cr*local_lai*0.5d0)

    ! calculate displacement (m); assume minimum lai 1.0 or 1.5 as height is not
    ! varied
    displacement = (1d0-((1d0-exp(-sqrt_cd1_lai))/sqrt_cd1_lai))*canopy_height

    ! calculate roughness sublayer influence function;
    ! this describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law
    ! phi_h = log(Cw)-1d0+Cw**(-1d0) ! DO NOT FORGET TO UPDATE IF Cw CHANGES

    ! finally calculate roughness length, dependant on displacement, friction
    ! velocity and lai.
    roughl = ((1d0-displacement/canopy_height)*exp(-vonkarman*ustar_Uh-phi_h))*canopy_height

    ! sanity check
!    if (roughl /= roughl) then
!        write(*,*)"TLS:  ERROR roughness length calculations"
!        write(*,*)"Roughness lenght", roughl, "Displacement", displacement
!        write(*,*)"canopy height", canopy_height, "lai", lai
!    endif

  end subroutine z0_displacement  
  !
  !------------------------------------------------------------------
  !
  ! Functions below this line
  !
  !------------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------
  !
  pure function arrhenious( a , b , t )

    ! The equation is simply...                        !
    !    a * exp( b * ( t - 25.0 ) / ( t + 273.15 ) )  !
    ! However, precision in this routine matters as it !
    ! affects many others. To maximise precision, the  !
    ! calculations have been split.                    !

    implicit none

    ! arguments..
    double precision,intent(in) :: a , b , t
    double precision            :: arrhenious

    arrhenious = a * exp( b * (t - 25d0) / (t + freeze) )

  end function arrhenious
  !
  !----------------------------------------------------------------------
  !
  double precision function calculate_declination(doy)

    implicit none

     ! Declare arguments
     double precision, intent(in) :: doy

     ! Declination calculation
     ! NOTE: 0.002739726d0 = 1/365
     !    dec = - asin( sin( 23.45d0 * deg_to_rad ) * cos( 2d0 * pi * ( doy + 10d0 ) / 365d0 ) )
     !    dec = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10d0 ) / 365d0 ) )
     calculate_declination = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10d0 ) * 0.002739726d0 ) )

     ! return to user
     return

  end function calculate_declination  
  !
  !------------------------------------------------------------------
  !
  double precision function linear_model_gradient(x,y,interval)


    ! Function to calculate the gradient of a linear model for a given depentent
    ! variable (y) based on predictive variable (x). The typical use of this
    ! function will in fact be to assume that x is time.

    implicit none

    ! declare input variables
    integer :: interval
    double precision, dimension(interval) :: x,y

    ! declare local variables
    double precision :: sum_x, sum_y!, sumsq_x,sum_product_xy

    ! calculate the sum of x
    sum_x = sum(x)
    ! calculate the sum of y
    sum_y = sum(y)
    ! calculate the sum of squares of x
    !sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    !sum_product_xy = sum(x*y)
    ! calculate the gradient
    !linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) )
    !&
    !                      / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )
    ! Linear regression done as single line to reduce assignment requirements
    linear_model_gradient = ( (dble(interval)*sum(x*y)) - (sum_x*sum_y) ) &
                          / ( (dble(interval)*sum(x*x)) - (sum_x*sum_x) )

    ! for future reference here is how to calculate the intercept
!    intercept = ( (sum_y*sumsq_x) - (sum_x*sum_product_xy) ) &
!              / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )

    ! don't forget to return to the user
    return

  end function linear_model_gradient
  !
  !----------------------------------------------------------------------
  !
  double precision function opt_max_scaling( max_val, min_val , optimum , kurtosis , current )

    ! Estimates a 0-1 scaling based on a skewed guassian distribution with a
    ! given optimum, maximum and kurtosis. Minimum is assumed to be at infinity
    ! (or near enough)

    implicit none

    ! arguments..
    double precision, intent(in) :: max_val, min_val, optimum, kurtosis, current

!    ! local variables..
!    double precision, parameter :: min_val = -1d6

    ! Code with implicit assumption of min bound at infinity
!    if ( current >= max_val ) then
!         opt_max_scaling = 0d0
!    else
!         dummy     = exp( log((max_val - current) / (max_val - optimum)) * kurtosis * (max_val - optimum) )
!         opt_max_scaling = dummy * exp( kurtosis * ( current - optimum ) )
!    end if

    ! Code with explicit min bound
    opt_max_scaling = exp( kurtosis * log((max_val-current)/(max_val-optimum)) * (max_val-optimum) ) &
                    * exp( kurtosis * log((current-min_val)/(optimum-min_val)) * (optimum-min_val) )
    ! Sanity check, allows for overlapping parameter ranges
    if (opt_max_scaling /= opt_max_scaling) opt_max_scaling = 0d0

  end function opt_max_scaling  
  !
  !------------------------------------------------------------------
  !
  double precision function water_retention_saxton_eqns( xin )

    ! field capacity calculations for saxton eqns !

    implicit none

    ! arguments..
    double precision, intent(in) :: xin

    ! local variables..
    double precision :: soil_wp

    ! calculate the soil water potential (kPa)..
    soil_wp = -potA(water_retention_pass) * xin**potB(water_retention_pass)
    water_retention_saxton_eqns = soil_wp + 10d0    ! 10 kPa represents air-entry swp

    return

  end function water_retention_saxton_eqns          
  !
  !------------------------------------------------------------------
  !
double precision function zbrent(called_from, f, x1, x2, tol, toltol)

  implicit none

  ! arguments..
  character(len=*), intent(in) :: called_from ! name of procedure calling (used to pass through for errors)
  double precision, intent(in) :: tol, toltol, x1, x2

  ! Function pointer for user-provided function
  double precision, external :: f

  ! local variables..
  integer                  :: iter
  integer, parameter      :: ITMAX = 8
  double precision         :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, tol0, xm
  double precision, parameter :: EPS = 6d-8

  ! calculations...
  a = x1
  b = x2
  fa = f(a)
  fb = f(b)
  tol0 = tol * 0.5d0

  ! Check that we haven't (by fluke) already started with the root...
  if (abs(fa) < toltol) then
    zbrent = a
    return
  elseif (abs(fb) < toltol) then
    zbrent = b
    return
  end if

  c = b
  fc = fb

  do iter = 1, ITMAX

    ! If the new value (f(c)) doesn't bracket
    ! the root with f(b) then adjust it..
    if (fb * fc > 0d0) then
      c = a
      fc = fa
      d = b - a
      e = d
    end if
    if (abs(fc) < abs(fb)) then
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
    end if
    tol1 = EPS * abs(b) + tol0
    xm = 0.5d0 * (c - b)
    if ((abs(xm) <= tol1) .or. (abs(fb) < toltol)) then
      zbrent = b
      return
    end if
    if ((abs(e) >= tol1) .and. (abs(fa) > abs(fb))) then
      s = fb / fa
      if (a == c) then
        p = 2d0 * xm * s
        q = 1d0 - s
      else
        q = fa / fc
        r = fb / fc
        p = s * (2d0 * xm * q * (q - r) - (b - a) * (r - 1d0))
        q = (q - 1d0) * (r - 1d0) * (s - 1d0)
      end if
      if (p > 0d0) q = -q
      p = abs(p)
      if ((2d0 * p) < min(3d0 * xm * q - abs(tol1 * q), abs(e * q))) then
        e = d
        d = p / q
      else
        d = xm
        e = d
      end if
    else
      d = xm
      e = d
    end if
    a = b
    fa = fb
    if (abs(d) > tol1) then
      b = b + d
    else
      b = b + sign(tol1, xm)
    end if
    fb = f(b)
  enddo

  zbrent = b

end function zbrent

  !
  !------------------------------------------------------------------
  !
!
!--------------------------------------------------------------------
!
end module CARBON_MODEL_MOD