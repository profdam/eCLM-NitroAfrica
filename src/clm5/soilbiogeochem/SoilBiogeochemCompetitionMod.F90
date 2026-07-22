module SoilBiogeochemCompetitionMod
  
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Resolve plant/heterotroph competition for mineral N
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use shr_const_mod                  , only : shr_const_tkfrz
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use clm_varcon                      , only : dzsoi_decomp, &
                                               zsoi, &   ! soil depth [m]; added by fkm for NH3 volatilization
                                               zlnd, &   ! Roughness length for soil [m]; added by fkm for canopy reduction
                                               spval !mvm for 1e36 deltaph values
  use clm_varctl                      , only : use_nitrif_denitrif
  use abortutils                      , only : endrun
  use decompMod                       , only : bounds_type
  use SoilBiogeochemStateType         , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemNitrogenUptakeMod , only : SoilBiogeochemNitrogenUptake
  use ColumnType                      , only : col                
  use CNVegstateType                  , only : cnveg_state_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegnitrogenstateType          , only : cnveg_nitrogenstate_type
  use CNVegnitrogenfluxType           , only : cnveg_nitrogenflux_type
  !use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use WaterStateType                  , only : waterstate_type
  use WaterfluxType                   , only : waterflux_type
  use TemperatureType                 , only : temperature_type
  use Atm2lndType                     , only : atm2lnd_type
  use PatchType                       , only : patch             ! added for NH3 canopy capture
  use FrictionVelocityMod             , only : frictionvel_type  ! added for NH3 canopy capture
  use DryDepVelocity                  , only : drydepvel_type    ! added for NOx canopy reduction factor (dry deposition)
  use CropType                        , only : crop_type         ! crop inputs (pH adjustment hook)
  use SoilStateType                   , only : soilstate_type
  use CanopyStateType                 , only : CanopyState_type
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams
  public :: SoilBiogeochemCompetitionInit         ! Initialization
  public :: SoilBiogeochemCompetition             ! run method

  type :: params_type
     real(r8) :: bdnr              ! bulk denitrification rate (1/s)
     real(r8) :: compet_plant_no3  ! (unitless) relative compettiveness of plants for NO3
     real(r8) :: compet_plant_nh4  ! (unitless) relative compettiveness of plants for NH4
     real(r8) :: compet_decomp_no3 ! (unitless) relative competitiveness of immobilizers for NO3
     real(r8) :: compet_decomp_nh4 ! (unitless) relative competitiveness of immobilizers for NH4
     real(r8) :: compet_denit      ! (unitless) relative competitiveness of denitrifiers for NO3
     real(r8) :: compet_nit        ! (unitless) relative competitiveness of nitrifiers for NH4
  end type params_type
  !
  type(params_type), private :: params_inst  ! params_inst is populated in readParamsMod  
  !
  ! !PUBLIC DATA MEMBERS:
  character(len=* ), public, parameter :: suplnAll='ALL'       ! Supplemental Nitrogen for all PFT's
  character(len=* ), public, parameter :: suplnNon='NONE'      ! No supplemental Nitrogen
  character(len=15), public            :: suplnitro = suplnNon ! Supplemental Nitrogen mode
  !
  ! !PRIVATE DATA MEMBERS:
  real(r8) :: dt   ! decomp timestep (seconds)
  real(r8) :: bdnr ! bulk denitrification rate (1/s)

  ! Val Martin-style rain-pulse memory for soil NOx from nitrification.
  ! The pulse is computed in SoilBiogeochemCompetition because final NOx
  ! must use realised f_nit_vr and f_denit_vr after N competition.
  !
  ! Restart-safety note:
  ! prev_h2osoi_vol_vr is an in-memory helper used only to detect the
  ! moisture increment between two consecutive calls. It is not a science
  ! state variable and is not written to the restart file. On a fresh start
  ! or after a model restart it is initialised to the current H2OSOI value,
  ! which makes h2osoi_diff = 0 for the first post-restart call. This avoids
  ! a false rain-pulse trigger at restart boundaries. The persistent pulse
  ! state itself is carried by ldry_vr and pfactor_vr in the nitrogen state.
  real(r8), allocatable, save :: prev_h2osoi_vol_vr(:,:)
  logical , save :: nox_pulse_memory_initialized = .false.

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t,ncd_io

    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNAllocParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! read in parameters

    tString='bdnr'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%bdnr=tempr

    tString='compet_plant_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%compet_plant_no3=tempr

    tString='compet_plant_nh4'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%compet_plant_nh4=tempr
   
    tString='compet_decomp_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%compet_decomp_no3=tempr

    tString='compet_decomp_nh4'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%compet_decomp_nh4=tempr
   
    tString='compet_denit'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%compet_denit=tempr

    tString='compet_nit'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%compet_nit=tempr   

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemCompetitionInit ( bounds)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use clm_varcon      , only: secspday
    use clm_time_manager, only: get_step_size
    use clm_varctl      , only: iulog, cnallocate_carbon_only_set
    use shr_infnan_mod  , only: nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'SoilBiogeochemCompetitionInit'
    logical :: carbon_only
    !-----------------------------------------------------------------------

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! set space-and-time parameters from parameter file
    bdnr = params_inst%bdnr * (dt/secspday)

    ! Change namelist settings into private logical variables
    select case(suplnitro)
    case(suplnNon)
       carbon_only = .false.
    case(suplnAll)
       carbon_only = .true.
    case default
       write(iulog,*) 'Supplemental Nitrogen flag (suplnitro) can only be: ', &
            suplnNon, ' or ', suplnAll
       call endrun(msg='ERROR: supplemental Nitrogen flag is not correct'//&
            errMsg(sourcefile, __LINE__))
    end select

    call cnallocate_carbon_only_set(carbon_only)

  end subroutine SoilBiogeochemCompetitionInit

  !-----------------------------------------------------------------------
   subroutine SoilBiogeochemCompetition (bounds, num_soilc, filter_soilc,num_soilp, filter_soilp, waterstate_inst, &
                                         waterflux_inst, temperature_inst,soilstate_inst,                          &
                                         cnveg_state_inst,cnveg_carbonstate_inst,                                  &
                                         cnveg_carbonflux_inst,cnveg_nitrogenstate_inst,cnveg_nitrogenflux_inst,   &
                                         soilbiogeochem_carbonflux_inst,                                           &              
                                         soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst,             &
                                         soilbiogeochem_nitrogenflux_inst,canopystate_inst,                        &
                                         atm2lnd_inst,         & !added by fkm for NH3 volativlity
                                         drydepvel_inst, & !added by mvm for NOx fluxes
                                         crop_inst, &  ! mvm added for NH3 ph adjustemnet
                                         frictionvel_inst) ! added by fkm for canopy reduction
    !
    ! !USES:
    use clm_varctl       , only: cnallocate_carbon_only, iulog
    use clm_varpar       , only: nlevdecomp, ndecomp_cascade_transitions
    use clm_varcon       , only: nitrif_n2o_loss_frac
    use CNSharedParamsMod, only: use_fun
    use CNFUNMod         , only: CNFUN
    use subgridAveMod    , only: p2c_2d, p2c
    use perf_mod         , only : t_startf, t_stopf
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds
    integer                                 , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:)  ! filter for soil patches
    type(waterstate_type)                   , intent(in)    :: waterstate_inst
    type(waterflux_type)                    , intent(in)    :: waterflux_inst
    type(temperature_type)                  , intent(in)    :: temperature_inst
    type(atm2lnd_type)                     , intent(in)    :: atm2lnd_inst
    type(drydepvel_type)                  , intent(in)    :: drydepvel_inst
    type(crop_type)                       , intent(in)    :: crop_inst
    type(frictionvel_type)                , intent(in)    :: frictionvel_inst
    type(soilstate_type)                    , intent(in)    :: soilstate_inst
    type(cnveg_state_type)                  , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_state_type)         , intent(inout) :: soilbiogeochem_state_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(canopystate_type)                  , intent(inout) :: canopystate_inst   

    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,l,pi,j,g                                            ! indices
    integer  :: fc                                                    ! filter column index
integer  :: fp                                                    ! patch filter index (NH3 canopy capture)
real(r8) :: canopy_top, canopy_bot, wind_canopy                   ! canopy geometry + wind proxy (m, m, m/s)
real(r8) :: nh3_conc, fac_canopy_height                           ! intermediate terms
real(r8) :: fac_leaf_geometry, fac_leaf_moisture, vg_nh3          ! leaf terms + dep velocity
real(r8) :: nh3_vol_to_canopy_vr_patch                            ! patch contribution
    logical :: local_use_fun                                          ! local version of use_fun
    real(r8) :: compet_plant_no3                                      ! (unitless) relative compettiveness of plants for NO3
    real(r8) :: compet_plant_nh4                                      ! (unitless) relative compettiveness of plants for NH4
    real(r8) :: compet_decomp_no3                                     ! (unitless) relative competitiveness of immobilizers for NO3
    real(r8) :: compet_decomp_nh4                                     ! (unitless) relative competitiveness of immobilizers for NH4
    real(r8) :: compet_denit                                          ! (unitless) relative competitiveness of denitrifiers for NO3
    real(r8) :: compet_nit                                            ! (unitless) relative competitiveness of nitrifiers for NH4
    real(r8) :: fpi_no3_vr(bounds%begc:bounds%endc,1:nlevdecomp)      ! fraction of potential immobilization supplied by no3(no units)
    real(r8) :: fpi_nh4_vr(bounds%begc:bounds%endc,1:nlevdecomp)      ! fraction of potential immobilization supplied by nh4 (no units)
    real(r8) :: sum_nh4_demand(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: sum_nh4_demand_scaled(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: sum_no3_demand(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: sum_no3_demand_scaled(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: sum_ndemand_vr(bounds%begc:bounds%endc, 1:nlevdecomp) !total column N demand (gN/m3/s) at a given level
    real(r8) :: nuptake_prof(bounds%begc:bounds%endc, 1:nlevdecomp)
    real(r8) :: sminn_tot(bounds%begc:bounds%endc)
    integer  :: nlimit(bounds%begc:bounds%endc,0:nlevdecomp)          !flag for N limitation
    integer  :: nlimit_no3(bounds%begc:bounds%endc,0:nlevdecomp)      !flag for NO3 limitation
    integer  :: nlimit_nh4(bounds%begc:bounds%endc,0:nlevdecomp)      !flag for NH4 limitation
    real(r8) :: residual_sminn_vr(bounds%begc:bounds%endc, 1:nlevdecomp)
    real(r8) :: residual_sminn(bounds%begc:bounds%endc)
    real(r8) :: residual_smin_nh4_vr(bounds%begc:bounds%endc, 1:nlevdecomp)
    real(r8) :: residual_smin_no3_vr(bounds%begc:bounds%endc, 1:nlevdecomp)
    real(r8) :: residual_smin_nh4(bounds%begc:bounds%endc)
    real(r8) :: residual_smin_no3(bounds%begc:bounds%endc)
    real(r8) :: residual_plant_ndemand(bounds%begc:bounds%endc)
    ! --- NH3 volatilisation (FANv2 / Fung-style) locals ---
    real(r8) :: tC, Kw, Ka, hydrogen, hydroxide
    real(r8) :: mol_nh4, mol_nh3, cvf, f_adsorption, aq_NH3
    real(r8) :: depth_fac, wind_fac
    real(r8) :: sminn_to_plant_new(bounds%begc:bounds%endc)
    real(r8) :: temp_fac

    ! --- Val Martin-style NOx and rain-pulse locals ---
    real(r8) :: dthr
    real(r8) :: h2osoi_diff
    real(r8) :: crf_drydep_col(bounds%begc:bounds%endc)
    real(r8) :: crf_drydep_patch_local
    real(r8) :: crf_col_wtsum(bounds%begc:bounds%endc)
    real(r8) :: pulse_weighted_num(bounds%begc:bounds%endc)
    real(r8) :: pulse_weighted_den(bounds%begc:bounds%endc)
    real(r8) :: fN2Onit_tmp
    real(r8) :: nox_pulse_raw
    real(r8) :: soilph_nox
    real(r8) :: wfps_frac_nox
    real(r8) :: afps_frac_nox
    real(r8) :: dr_nox
    real(r8) :: nox_nit_moist_fac
    real(r8) :: nox_denit_moist_fac
    real(r8) :: nox_nit_diff_fac
    real(r8) :: nox_denit_diff_fac
    real(r8) :: y_no_nit
    real(r8) :: y_no_denit
    real(r8) :: nox_denit_direct
    real(r8) :: nox_nit_before_pulse
    real(r8) :: nox_nit_after_pulse
    real(r8), parameter :: nox_dry_theta_crit    = 0.175_r8
    real(r8), parameter :: nox_rewet_dtheta_crit = 0.005_r8
    real(r8), parameter :: nox_min_dry_hours     = 72._r8
    real(r8), parameter :: nox_pulse_a           = 13.01_r8
    real(r8), parameter :: nox_pulse_b           = 53.6_r8
    real(r8), parameter :: nox_pulse_decay_hr    = 0.068_r8
    real(r8), parameter :: nox_pulse_min         = 1.0e-6_r8
    real(r8), parameter :: nox_canopy_vg         = 0.05_r8

    ! Parton/DAYCENT-inspired realised-process NOx implementation.
    !
    ! This is not the Val Martin shortcut NOx = N2O * NOx:N2O.  NOx is treated
    ! as a separate gas product from the realised process fluxes:
    !
    !   F_NOx_NIT   = F_NIT   * Y_NO_NIT(WFPS, diffusivity) * PULSE_FAC
    !   F_NOx_DENIT = F_DENIT * Y_NO_DENIT(WFPS, diffusivity)
    !
    ! The yield terms are environmental functions, not fixed percentages.
    ! The maximum yields below are calibration ceilings.  The realised yield is
    ! reduced by WFPS and air-filled-pore-space diffusivity terms.  This keeps
    ! the implementation closer to the Parton/DAYCENT and Davidson logic:
    ! nitrification NO is favoured in moderately moist, aerated soil, while
    ! denitrification NO is allowed in wetter transition conditions but is
    ! suppressed when gas diffusivity is low.
    real(r8), parameter :: y_no_nit_max            = 0.070_r8
    real(r8), parameter :: y_no_denit_max          = 0.030_r8
    real(r8), parameter :: nox_nit_wfps_opt        = 0.45_r8
    real(r8), parameter :: nox_nit_wfps_width      = 0.20_r8
    real(r8), parameter :: nox_denit_wfps_opt      = 0.65_r8
    real(r8), parameter :: nox_denit_wfps_width    = 0.18_r8
    real(r8), parameter :: nox_nit_moist_fac_min   = 0.02_r8
    real(r8), parameter :: nox_denit_moist_fac_min = 0.00_r8
    real(r8), parameter :: nox_nit_dr_scale        = 0.15_r8
    real(r8), parameter :: nox_denit_dr_scale      = 0.10_r8
    !-----------------------------------------------------------------------

    associate(                                                                                           &
         fpg                          => soilbiogeochem_state_inst%fpg_col                             , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)    
         fpi                          => soilbiogeochem_state_inst%fpi_col                             , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
         fpi_vr                       => soilbiogeochem_state_inst%fpi_vr_col                          , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)
         nfixation_prof               => soilbiogeochem_state_inst%nfixation_prof_col                  , & ! Output: [real(r8) (:,:) ]                                        
         plant_ndemand                => soilbiogeochem_state_inst%plant_ndemand_col                   , & ! Input:  [real(r8) (:)   ]  column-level plant N demand

         sminn_vr                     => soilbiogeochem_nitrogenstate_inst%sminn_vr_col                , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N                
         smin_nh4_vr                  => soilbiogeochem_nitrogenstate_inst%smin_nh4_vr_col             , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral NH4              
         smin_no3_vr                  => soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col             , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral NO3              
         h2osoi_vol                    => waterstate_inst%h2osoi_vol_col                             , & ! Input:  [real(r8) (:,:) ] soil water volume [m3/m3]
         watsat                        => soilstate_inst%watsat_col                                   , & ! Input:  [real(r8) (:,:) ] saturated volumetric water content [m3/m3]
         t_soisno                      => temperature_inst%t_soisno_col                               , & ! Input:  [real(r8) (:,:) ] soil temperature [K]

         forc_wind                     => atm2lnd_inst%forc_wind_grc                                  , & ! Input:  [real(r8) (:)   ] wind speed (m/s)
         cellclay                      => soilstate_inst%cellclay_col                                  , & ! Input:  [real(r8) (:,:) ] clay content (%)
         fr_ph_nh3_vr                  => soilbiogeochem_nitrogenflux_inst%fr_ph_nh3_vr_col           , & ! Output: [real(r8) (:,:) ] pH factor for NH3
         pot_f_nh3_vol_vr              => soilbiogeochem_nitrogenflux_inst%pot_f_nh3_vol_vr_col       , & ! Output: [real(r8) (:,:) ] potential NH3 volatilisation [gN/m3/s]
         f_nh3_vol_vr                  => soilbiogeochem_nitrogenflux_inst%f_nh3_vol_vr_col           , & ! Output: [real(r8) (:,:) ] NH3 volatilisation [gN/m3/s]

         f_nh3_vol_to_atmos_vr         => soilbiogeochem_nitrogenflux_inst%f_nh3_vol_to_atmos_vr_col   , & ! Output: [real(r8) (:,:) ] NH3 volatilisation to atmosphere [gN/m3/s]
         f_nh3_vol_to_canopy_vr        => soilbiogeochem_nitrogenflux_inst%f_nh3_vol_to_canopy_vr_col  , & ! Output: [real(r8) (:,:) ] NH3 volatilisation captured by canopy [gN/m3/s]

         aq_nh3_diag_vr               => soilbiogeochem_nitrogenflux_inst%aq_nh3_vr_col              , & ! Output: [real(r8) (:,:) ] available NH3 pool [gN/m3]
         cvf_nh3_diag_vr              => soilbiogeochem_nitrogenflux_inst%cvf_nh3_vr_col             , & ! Output: [real(r8) (:,:) ] NH4/(NH4+NH3) partition factor
         fads_nh3_diag_vr             => soilbiogeochem_nitrogenflux_inst%fads_nh3_vr_col            , & ! Output: [real(r8) (:,:) ] adsorption factor
         h2osoi_nh3_diag_vr           => soilbiogeochem_nitrogenflux_inst%h2osoi_nh3_vr_col          , & ! Output: [real(r8) (:,:) ] soil water used in NH3 calc [m3/m3]
         tsoi_nh3_diag_vr             => soilbiogeochem_nitrogenflux_inst%tsoi_nh3_vr_col            , & ! Output: [real(r8) (:,:) ] soil temperature used in NH3 calc [K]
         pot_f_nit_vr                 => soilbiogeochem_nitrogenflux_inst%pot_f_nit_vr_col             , & ! Input:  [real(r8) (:,:) ]  (gN/m3/s) potential soil nitrification flux
         pot_f_denit_vr               => soilbiogeochem_nitrogenflux_inst%pot_f_denit_vr_col           , & ! Input:  [real(r8) (:,:) ]  (gN/m3/s) potential soil denitrification flux
         f_nit_vr                     => soilbiogeochem_nitrogenflux_inst%f_nit_vr_col                 , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil nitrification flux     
         f_denit_vr                   => soilbiogeochem_nitrogenflux_inst%f_denit_vr_col               , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil denitrification flux   
         potential_immob              => soilbiogeochem_nitrogenflux_inst%potential_immob_col          , & ! Output: [real(r8) (:)   ]                                          
         actual_immob                 => soilbiogeochem_nitrogenflux_inst%actual_immob_col             , & ! Output: [real(r8) (:)   ]                                          
         sminn_to_plant               => soilbiogeochem_nitrogenflux_inst%sminn_to_plant_col           , & ! Output: [real(r8) (:)   ]                                          
         sminn_to_denit_excess_vr     => soilbiogeochem_nitrogenflux_inst%sminn_to_denit_excess_vr_col , & ! Output: [real(r8) (:,:) ]                                        
         actual_immob_no3_vr          => soilbiogeochem_nitrogenflux_inst%actual_immob_no3_vr_col      , & ! Output: [real(r8) (:,:) ]                                        
         actual_immob_nh4_vr          => soilbiogeochem_nitrogenflux_inst%actual_immob_nh4_vr_col      , & ! Output: [real(r8) (:,:) ]                                        
         smin_no3_to_plant_vr         => soilbiogeochem_nitrogenflux_inst%smin_no3_to_plant_vr_col     , & ! Output: [real(r8) (:,:) ]                                        
         smin_nh4_to_plant_vr         => soilbiogeochem_nitrogenflux_inst%smin_nh4_to_plant_vr_col     , & ! Output: [real(r8) (:,:) ]                                        
         n2_n2o_ratio_denit_vr        => soilbiogeochem_nitrogenflux_inst%n2_n2o_ratio_denit_vr_col    , & ! Output: [real(r8) (:,:) ]  ratio of N2 to N2O production by denitrification [gN/gN]
         f_n2o_denit_vr               => soilbiogeochem_nitrogenflux_inst%f_n2o_denit_vr_col           , & ! Output: [real(r8) (:,:) ]  flux of N2O from denitrification [gN/m3/s]
         f_n2o_nit_vr                 => soilbiogeochem_nitrogenflux_inst%f_n2o_nit_vr_col             , & ! Output: [real(r8) (:,:) ]  flux of N2O from nitrification [gN/m3/s]
         supplement_to_sminn_vr       => soilbiogeochem_nitrogenflux_inst%supplement_to_sminn_vr_col   , & ! Output: [real(r8) (:,:) ]                                        
         sminn_to_plant_vr            => soilbiogeochem_nitrogenflux_inst%sminn_to_plant_vr_col        , & ! Output: [real(r8) (:,:) ]                                        
         potential_immob_vr           => soilbiogeochem_nitrogenflux_inst%potential_immob_vr_col       , & ! Input:  [real(r8) (:,:) ]                                        
         actual_immob_vr              => soilbiogeochem_nitrogenflux_inst%actual_immob_vr_col          , & ! Output: [real(r8) (:,:) ]                                        
         sminn_to_plant_fun_vr        => soilbiogeochem_nitrogenflux_inst%sminn_to_plant_fun_vr_col    , & ! Iutput: [real(r8) (:)   ]  Total layer soil N uptake of FUN (gN/m2/s) 
         sminn_to_plant_fun_no3_vr    => soilbiogeochem_nitrogenflux_inst%sminn_to_plant_fun_no3_vr_col, & ! Iutput: [real(r8) (:)   ]  Total layer no3 uptake of FUN (gN/m2/s)
         sminn_to_plant_fun_nh4_vr    => soilbiogeochem_nitrogenflux_inst%sminn_to_plant_fun_nh4_vr_col  , & ! Iutput: [real(r8) (:)   ]  Total layer nh4 uptake of FUN (gN/m2/s)

         ! Val Martin-style NOx fields. NOx is constructed here because this
         ! subroutine has realised f_nit_vr and f_denit_vr after competition.
         nox_n2o_ratio_vr             => soilbiogeochem_nitrogenflux_inst%nox_n2o_ratio_vr_col            , &
         fN2Onit_vr                   => soilbiogeochem_nitrogenflux_inst%fN2Onit_vr_col                  , &
         pot_f_nox_nit_vr             => soilbiogeochem_nitrogenflux_inst%pot_f_nox_nit_vr_col            , &
         pot_f_nox_denit_vr           => soilbiogeochem_nitrogenflux_inst%pot_f_nox_denit_vr_col          , &
         f_nox_nit_vr                 => soilbiogeochem_nitrogenflux_inst%f_nox_nit_vr_col                , &
         f_nox_denit_vr               => soilbiogeochem_nitrogenflux_inst%f_nox_denit_vr_col              , &
         f_nox_nit_atmos_vr           => soilbiogeochem_nitrogenflux_inst%f_nox_nit_atmos_vr_col          , &
         f_nox_denit_atmos_vr         => soilbiogeochem_nitrogenflux_inst%f_nox_denit_atmos_vr_col        , &
         pot_f_nox_nit_col            => soilbiogeochem_nitrogenflux_inst%pot_f_nox_nit_col               , &
         pot_f_nox_denit_col          => soilbiogeochem_nitrogenflux_inst%pot_f_nox_denit_col             , &
         f_nox_nit_col                => soilbiogeochem_nitrogenflux_inst%f_nox_nit_col                   , &
         f_nox_denit_col              => soilbiogeochem_nitrogenflux_inst%f_nox_denit_col                 , &
         f_nox_nit_atmos_col          => soilbiogeochem_nitrogenflux_inst%f_nox_nit_atmos_col             , &
         f_nox_denit_atmos_col        => soilbiogeochem_nitrogenflux_inst%f_nox_denit_atmos_col           , &
         soil_nox_total_col           => soilbiogeochem_nitrogenflux_inst%soil_nox_total_col              , &
         soil_nox_crop_col            => soilbiogeochem_nitrogenflux_inst%soil_nox_crop_col               , &
         pulse_fac_col                => soilbiogeochem_nitrogenflux_inst%pulse_fac_col                   , &
         pfactor_vr                   => soilbiogeochem_nitrogenstate_inst%pfactor_vr_col                  , &
         ldry_vr                      => soilbiogeochem_nitrogenstate_inst%ldry_vr_col                       &
         )

      ! calcualte nitrogen uptake profile
      ! nuptake_prof(:,:) = nan
      ! call SoilBiogelchemNitrogenUptakeProfile(bounds, &
      !     nlevdecomp, num_soilc, filter_soilc, &
      !     sminn_vr, dzsoi_decomp, nfixation_prof, nuptake_prof)

      ! column loops to resolve plant/heterotroph competition for mineral N

      sminn_to_plant_new(bounds%begc:bounds%endc)  =  0._r8

      ! Timestep in hours for the Val Martin rain-pulse decay and dry counter.
      dthr = dt / 3600._r8

      ! NOx canopy-reduction pathway.
      !
      ! The Val Martin branch stores a precomputed patch-level NOx canopy
      ! reduction factor in drydepvel_inst%crf_drydep_patch. This eCLM branch
      ! does not contain that field, so using that name causes a compile error.
      !
      ! Instead, compute the reduction directly from canopy quantities already
      ! available in this subroutine. This keeps the same conceptual order as
      ! Val Martin: final NOx is computed first, then only the above-canopy
      ! atmospheric part is reduced.
      !
      ! crf_drydep_col = 1 means no canopy loss.
      ! crf_drydep_col = 0 means full canopy removal before the free atmosphere.
      crf_drydep_col(bounds%begc:bounds%endc) = 0._r8
      crf_col_wtsum(bounds%begc:bounds%endc) = 0._r8

      do fp = 1, num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)

         canopy_top        = max(zlnd, canopystate_inst%htop_patch(p))
         canopy_bot        = max(zlnd, canopystate_inst%hbot_patch(p))
         wind_canopy       = max(0.001_r8, frictionvel_inst%fv_patch(p))
         fac_canopy_height = max(0._r8, canopy_top - canopy_bot)
         fac_leaf_geometry = max(0._r8, canopystate_inst%tlai_patch(p))
         fac_leaf_moisture = min(1._r8, max(0._r8, waterstate_inst%rh_af_patch(p)))

         ! Canopy transmission factor for soil NOx.
         ! Larger canopy height, LAI, and leaf moisture increase removal.
         ! Larger wind/friction velocity reduces canopy residence time.
         crf_drydep_patch_local = exp(-nox_canopy_vg * fac_canopy_height * &
              fac_leaf_geometry * fac_leaf_moisture / wind_canopy)
         crf_drydep_patch_local = max(0._r8, min(1._r8, crf_drydep_patch_local))

         crf_drydep_col(c) = crf_drydep_col(c) + crf_drydep_patch_local * patch%wtcol(p)
         crf_col_wtsum(c)  = crf_col_wtsum(c)  + patch%wtcol(p)
      end do

      do fc = 1, num_soilc
         c = filter_soilc(fc)
         if (crf_col_wtsum(c) > 0._r8) then
            crf_drydep_col(c) = crf_drydep_col(c) / crf_col_wtsum(c)
         else
            crf_drydep_col(c) = 1._r8
         end if
      end do

      ! Allocate previous soil-water memory once. This detects rewetting events.
      !
      ! Restart-safe behaviour:
      ! When the model starts or restarts, prev_h2osoi_vol_vr is initialised
      ! from the current soil water. Therefore the first computed moisture
      ! increment after restart is zero, preventing a spurious pulse. The
      ! longer-lived pulse state variables ldry_vr and pfactor_vr are held in
      ! soilbiogeochem_nitrogenstate_inst and should be included in restart
      ! handling with the rest of the nitrogen state.
      if (.not. nox_pulse_memory_initialized) then
         allocate(prev_h2osoi_vol_vr(bounds%begc:bounds%endc,1:nlevdecomp))
         prev_h2osoi_vol_vr(bounds%begc:bounds%endc,1:nlevdecomp) = &
              h2osoi_vol(bounds%begc:bounds%endc,1:nlevdecomp)
         nox_pulse_memory_initialized = .true.
      end if

      ! Reset column NOx diagnostics before building them from realised layer fluxes.
      f_nox_nit_col(bounds%begc:bounds%endc)         = 0._r8
      f_nox_denit_col(bounds%begc:bounds%endc)       = 0._r8
      f_nox_nit_atmos_col(bounds%begc:bounds%endc)   = 0._r8
      f_nox_denit_atmos_col(bounds%begc:bounds%endc) = 0._r8
      pot_f_nox_nit_col(bounds%begc:bounds%endc)     = 0._r8
      pot_f_nox_denit_col(bounds%begc:bounds%endc)   = 0._r8
      soil_nox_total_col(bounds%begc:bounds%endc)    = 0._r8
      soil_nox_crop_col(bounds%begc:bounds%endc)     = 0._r8
      pulse_fac_col(bounds%begc:bounds%endc)         = 1._r8
      pulse_weighted_num(bounds%begc:bounds%endc)    = 0._r8
      pulse_weighted_den(bounds%begc:bounds%endc)    = 0._r8

      local_use_fun = use_fun

      if (.not. use_nitrif_denitrif) then

         ! init sminn_tot
         do fc=1,num_soilc
            c = filter_soilc(fc)
            sminn_tot(c) = 0.
         end do

         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)
               sminn_tot(c) = sminn_tot(c) + sminn_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)      
               if (sminn_tot(c)  >  0.) then
                  nuptake_prof(c,j) = sminn_vr(c,j) / sminn_tot(c)
               else
                  nuptake_prof(c,j) = nfixation_prof(c,j)
               endif
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)      
               sum_ndemand_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j) + potential_immob_vr(c,j)
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)      
               l = col%landunit(c)
               if (sum_ndemand_vr(c,j)*dt < sminn_vr(c,j)) then

                  ! N availability is not limiting immobilization or plant
                  ! uptake, and both can proceed at their potential rates
                  nlimit(c,j) = 0
                  fpi_vr(c,j) = 1.0_r8
                  actual_immob_vr(c,j) = potential_immob_vr(c,j)
                  sminn_to_plant_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j)
               else if ( cnallocate_carbon_only()) then !.or. &
                  ! this code block controls the addition of N to sminn pool
                  ! to eliminate any N limitation, when Carbon_Only is set.  This lets the
                  ! model behave essentially as a carbon-only model, but with the
                  ! benefit of keeping track of the N additions needed to
                  ! eliminate N limitations, so there is still a diagnostic quantity
                  ! that describes the degree of N limitation at steady-state.

                  nlimit(c,j) = 1
                  fpi_vr(c,j) = 1.0_r8
                  actual_immob_vr(c,j) = potential_immob_vr(c,j)
                  sminn_to_plant_vr(c,j) =  plant_ndemand(c) * nuptake_prof(c,j)
                  supplement_to_sminn_vr(c,j) = sum_ndemand_vr(c,j) - (sminn_vr(c,j)/dt)
               else
                  ! N availability can not satisfy the sum of immobilization and
                  ! plant growth demands, so these two demands compete for available
                  ! soil mineral N resource.

                  nlimit(c,j) = 1
                  if (sum_ndemand_vr(c,j) > 0.0_r8) then
                     actual_immob_vr(c,j) = (sminn_vr(c,j)/dt)*(potential_immob_vr(c,j) / sum_ndemand_vr(c,j))
                  else
                     actual_immob_vr(c,j) = 0.0_r8
                  end if

                  if (potential_immob_vr(c,j) > 0.0_r8) then
                     fpi_vr(c,j) = actual_immob_vr(c,j) / potential_immob_vr(c,j)
                  else
                     fpi_vr(c,j) = 0.0_r8
                  end if

                  sminn_to_plant_vr(c,j) = (sminn_vr(c,j)/dt) - actual_immob_vr(c,j)
               end if
            end do
         end do

         if ( local_use_fun ) then
            call t_startf( 'CNFUN' )
            call CNFUN(bounds,num_soilc,filter_soilc,num_soilp,filter_soilp,waterstate_inst                 ,&
                      waterflux_inst,temperature_inst,soilstate_inst,cnveg_state_inst,cnveg_carbonstate_inst,&
                      cnveg_carbonflux_inst,cnveg_nitrogenstate_inst,cnveg_nitrogenflux_inst                ,&
                      soilbiogeochem_nitrogenflux_inst,soilbiogeochem_carbonflux_inst,canopystate_inst,      &
                      soilbiogeochem_nitrogenstate_inst)
            call p2c_2d(bounds, nlevdecomp, &
                      cnveg_nitrogenflux_inst%sminn_to_plant_fun_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
                      soilbiogeochem_nitrogenflux_inst%sminn_to_plant_fun_vr_col(bounds%begc:bounds%endc,1:nlevdecomp), &
                      'unity')
            call t_stopf( 'CNFUN' )
         end if

         ! sum up N fluxes to plant
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)    
               sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
               if ( local_use_fun ) then
                  if (sminn_to_plant_fun_vr(c,j).gt.sminn_to_plant_vr(c,j)) then
                      sminn_to_plant_fun_vr(c,j)  = sminn_to_plant_vr(c,j)
                  end if
               end if
            end do
         end do

         ! give plants a second pass to see if there is any mineral N left over with which to satisfy residual N demand.
         do fc=1,num_soilc
            c = filter_soilc(fc)    
            residual_sminn(c) = 0._r8
         end do

         ! sum up total N left over after initial plant and immobilization fluxes
         do fc=1,num_soilc
            c = filter_soilc(fc)    
            residual_plant_ndemand(c) = plant_ndemand(c) - sminn_to_plant(c)
         end do
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)    
               if (residual_plant_ndemand(c)  >  0._r8 ) then
                  if (nlimit(c,j) .eq. 0) then
                     residual_sminn_vr(c,j) = max(sminn_vr(c,j) - (actual_immob_vr(c,j) + sminn_to_plant_vr(c,j) ) * dt, 0._r8)
                     residual_sminn(c) = residual_sminn(c) + residual_sminn_vr(c,j) * dzsoi_decomp(j)
                  else
                     residual_sminn_vr(c,j)  = 0._r8
                  endif
               endif
            end do
         end do

         ! distribute residual N to plants
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)    
               if ( residual_plant_ndemand(c)  >  0._r8 .and. residual_sminn(c)  >  0._r8 .and. nlimit(c,j) .eq. 0) then
                  sminn_to_plant_vr(c,j) = sminn_to_plant_vr(c,j) + residual_sminn_vr(c,j) * &
                       min(( residual_plant_ndemand(c) *  dt ) / residual_sminn(c), 1._r8) / dt
               endif
            end do
         end do

         ! re-sum up N fluxes to plant
         do fc=1,num_soilc
            c = filter_soilc(fc)    
            sminn_to_plant(c) = 0._r8
         end do
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)    
               sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
               if ( .not. local_use_fun ) then
                  sum_ndemand_vr(c,j) = potential_immob_vr(c,j) + sminn_to_plant_vr(c,j)
               else
                  sminn_to_plant_new(c)  = sminn_to_plant_new(c)   + sminn_to_plant_fun_vr(c,j) * dzsoi_decomp(j)
                  sum_ndemand_vr(c,j)    = potential_immob_vr(c,j) + sminn_to_plant_fun_vr(c,j)
               end if
            end do
         end do

         ! under conditions of excess N, some proportion is assumed to
         ! be lost to denitrification, in addition to the constant
         ! proportion lost in the decomposition pathways
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)    
               if ( .not. local_use_fun ) then
                  if ((sminn_to_plant_vr(c,j) + actual_immob_vr(c,j))*dt < sminn_vr(c,j)) then
                     sminn_to_denit_excess_vr(c,j) = max(bdnr*((sminn_vr(c,j)/dt) - sum_ndemand_vr(c,j)),0._r8)
                  else
                     sminn_to_denit_excess_vr(c,j) = 0._r8
                  endif
               else
                  if ((sminn_to_plant_fun_vr(c,j)  + actual_immob_vr(c,j))*dt < sminn_vr(c,j))  then
                     sminn_to_denit_excess_vr(c,j) = max(bdnr*((sminn_vr(c,j)/dt) - sum_ndemand_vr(c,j)),0._r8)
                  else
                     sminn_to_denit_excess_vr(c,j) = 0._r8
                  endif
               end if
            end do
         end do

         ! sum up N fluxes to immobilization
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)    
               actual_immob(c) = actual_immob(c) + actual_immob_vr(c,j) * dzsoi_decomp(j)
               potential_immob(c) = potential_immob(c) + potential_immob_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         do fc=1,num_soilc
            c = filter_soilc(fc)    
            ! calculate the fraction of potential growth that can be
            ! acheived with the N available to plants      
            if (plant_ndemand(c) > 0.0_r8) then
               if ( .not. local_use_fun ) then
                  fpg(c) = sminn_to_plant(c) / plant_ndemand(c)
               else
                  fpg(c) = sminn_to_plant_new(c) / plant_ndemand(c)
               end if
            else
               fpg(c) = 1.0_r8
            end if

            ! calculate the fraction of immobilization realized (for diagnostic purposes)
            if (potential_immob(c) > 0.0_r8) then
               fpi(c) = actual_immob(c) / potential_immob(c)
            else
               fpi(c) = 1.0_r8
            end if
         end do

      else  !----------NITRIF_DENITRIF-------------!

         ! column loops to resolve plant/heterotroph/nitrifier/denitrifier competition for mineral N
         !read constants from external netcdf file
         compet_plant_no3  = params_inst%compet_plant_no3
         compet_plant_nh4  = params_inst%compet_plant_nh4
         compet_decomp_no3 = params_inst%compet_decomp_no3
         compet_decomp_nh4 = params_inst%compet_decomp_nh4
         compet_denit      = params_inst%compet_denit
         compet_nit        = params_inst%compet_nit

         ! init total mineral N pools
         do fc=1,num_soilc
            c = filter_soilc(fc)
            sminn_tot(c) = 0.
         end do

         ! sum up total mineral N pools
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)
               sminn_tot(c) = sminn_tot(c) + (smin_no3_vr(c,j) + smin_nh4_vr(c,j)) * dzsoi_decomp(j)
            end do
         end do

         ! define N uptake profile for initial vertical distribution of plant N uptake, assuming plant seeks N from where it is most abundant
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)
               if (sminn_tot(c)  >  0.) then
                  nuptake_prof(c,j) = sminn_vr(c,j) / sminn_tot(c)
               else
                  nuptake_prof(c,j) = nfixation_prof(c,j)
               endif
            end do
         end do


           ! --------------------------------------------------------------------
! NH3 volatilisation potential (Val-Martin-compatible, current-tree safe)
! --------------------------------------------------------------------
do j = 1, nlevdecomp
   do fc = 1, num_soilc
      c = filter_soilc(fc)
      g = col%gridcell(c)
      l = col%landunit(c)

      ! initialise NH3 diagnostics and fluxes
      fr_ph_nh3_vr(c,j)          = 1.0_r8
      pot_f_nh3_vol_vr(c,j)      = 0.0_r8
      f_nh3_vol_vr(c,j)          = 0.0_r8
      f_nh3_vol_to_atmos_vr(c,j) = 0.0_r8
      f_nh3_vol_to_canopy_vr(c,j)= 0.0_r8

      aq_nh3_diag_vr(c,j)        = 0.0_r8
      cvf_nh3_diag_vr(c,j)       = 0.0_r8
      fads_nh3_diag_vr(c,j)      = 0.0_r8
      h2osoi_nh3_diag_vr(c,j)    = 0.0_r8
      tsoi_nh3_diag_vr(c,j)      = 0.0_r8

      if (smin_nh4_vr(c,j) > 0.0_r8 .and. h2osoi_vol(c,j) > 0.0_r8) then

         ! adsorption factor based on clay fraction
         f_adsorption = 0.99_r8 * (7.2733_r8 * (cellclay(c,j) / 100._r8) ** 3.0_r8 - &
              11.22_r8 * (cellclay(c,j) / 100._r8) ** 2.0_r8 + &
              5.7198_r8 * (cellclay(c,j) / 100._r8) + 0.0263_r8)
         f_adsorption = max(0.01_r8, min(f_adsorption, 0.999_r8))

         ! equilibrium constants
         Kw = 10._r8 ** (0.08946_r8 + 0.03605_r8 * (t_soisno(c,j) - shr_const_tkfrz)) * 1.0e-15_r8
         Ka = (1.416_r8 + 0.01357_r8 * (t_soisno(c,j) - shr_const_tkfrz)) * 1.0e-5_r8

         ! fixed soil pH = 6.5 for current-tree-safe Dahra test
         hydrogen  = 10.0_r8 ** (-6.0_r8)
         hydroxide = Kw / max(hydrogen, 1.0e-30_r8)

         ! Val-Martin concentration conversion
         mol_nh4 = smin_nh4_vr(c,j) / 14.0_r8 * 1000.0_r8 / max(h2osoi_vol(c,j), 1.0e-12_r8)
         mol_nh3 = mol_nh4 * hydroxide / max(Ka, 1.0e-30_r8)
         cvf     = mol_nh4 / max(mol_nh4 + mol_nh3, 1.0e-30_r8)

         ! keep pH factor neutral unless you later wire explicit soil pH
         fr_ph_nh3_vr(c,j) = 0.6_r8

         aq_NH3 = fr_ph_nh3_vr(c,j) * smin_nh4_vr(c,j) * (1.0_r8 - f_adsorption) * (1.0_r8 - cvf)

         ! diagnostics
         aq_nh3_diag_vr(c,j)     = aq_NH3
         cvf_nh3_diag_vr(c,j)    = cvf
         fads_nh3_diag_vr(c,j)   = f_adsorption
         h2osoi_nh3_diag_vr(c,j) = h2osoi_vol(c,j)
         tsoi_nh3_diag_vr(c,j)   = t_soisno(c,j)

         ! Val-Martin-style volatilisation controls
         depth_fac = (zsoi(nlevdecomp) - zsoi(j)) / zsoi(nlevdecomp)
         depth_fac = max(min(depth_fac, 1.0_r8), 0.0_r8)

         wind_fac  = 1.5_r8 * max(forc_wind(g), 0.0_r8) / (1.0_r8 + max(forc_wind(g), 0.0_r8))

         temp_fac  = max(t_soisno(c,j) - shr_const_tkfrz, 0.0_r8) / &
                    (50.0_r8 + max(t_soisno(c,j) - shr_const_tkfrz, 0.0_r8))

         pot_f_nh3_vol_vr(c,j) = aq_NH3 * depth_fac * wind_fac * temp_fac

         ! convert [gN/m3] to [gN/m3/s], bounded by available NH3 and NH4
         pot_f_nh3_vol_vr(c,j) = max(min(pot_f_nh3_vol_vr(c,j), aq_NH3), 0.0_r8) / dt
         pot_f_nh3_vol_vr(c,j) = min(pot_f_nh3_vol_vr(c,j), smin_nh4_vr(c,j) / dt)

      end if
   end do
end do



         ! main column/vertical loop
         do j = 1, nlevdecomp  
            do fc=1,num_soilc
               c = filter_soilc(fc)
               l = col%landunit(c)

               !  first compete for nh4
               sum_nh4_demand(c,j) = plant_ndemand(c) * nuptake_prof(c,j) + potential_immob_vr(c,j) + pot_f_nit_vr(c,j) + pot_f_nh3_vol_vr(c,j)
               sum_nh4_demand_scaled(c,j) = plant_ndemand(c)* nuptake_prof(c,j) * compet_plant_nh4 + &
                    potential_immob_vr(c,j)*compet_decomp_nh4 + pot_f_nit_vr(c,j)*compet_nit + pot_f_nh3_vol_vr(c,j)

               if (sum_nh4_demand(c,j)*dt < smin_nh4_vr(c,j)) then

                  ! NH4 availability is not limiting immobilization or plant
                  ! uptake, and all can proceed at their potential rates
                  nlimit_nh4(c,j) = 0
                  fpi_nh4_vr(c,j) = 1.0_r8
                  actual_immob_nh4_vr(c,j) = potential_immob_vr(c,j)
                  !RF added new term. 

                  f_nit_vr(c,j) = pot_f_nit_vr(c,j)
                  
                  f_nh3_vol_vr(c,j) = pot_f_nh3_vol_vr(c,j)
                  if ( .not. local_use_fun ) then
                     smin_nh4_to_plant_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j)
                  else
                     smin_nh4_to_plant_vr(c,j) = smin_nh4_vr(c,j)/dt - actual_immob_nh4_vr(c,j) - f_nit_vr(c,j) - f_nh3_vol_vr(c,j)
                  end if

               else

                  ! NH4 availability can not satisfy the sum of immobilization, nitrification, and
                  ! plant growth demands, so these three demands compete for available
                  ! soil mineral NH4 resource.
                  nlimit_nh4(c,j) = 1
                  if (sum_nh4_demand(c,j) > 0.0_r8) then
                  ! RF microbes compete based on the hypothesised plant demand. 
                     actual_immob_nh4_vr(c,j) = min((smin_nh4_vr(c,j)/dt)*(potential_immob_vr(c,j)* &
                          compet_decomp_nh4 / sum_nh4_demand_scaled(c,j)), potential_immob_vr(c,j))

                     f_nit_vr(c,j) =  min((smin_nh4_vr(c,j)/dt)*(pot_f_nit_vr(c,j)*compet_nit / &
                          sum_nh4_demand_scaled(c,j)), pot_f_nit_vr(c,j))
                     f_nh3_vol_vr(c,j) =  min((smin_nh4_vr(c,j)/dt)*(pot_f_nh3_vol_vr(c,j) / &
                          sum_nh4_demand_scaled(c,j)), pot_f_nh3_vol_vr(c,j))

                                                 
                     if ( .not. local_use_fun ) then
                         smin_nh4_to_plant_vr(c,j) = min((smin_nh4_vr(c,j)/dt)*(plant_ndemand(c)* &
                          nuptake_prof(c,j)*compet_plant_nh4 / sum_nh4_demand_scaled(c,j)), plant_ndemand(c)*nuptake_prof(c,j))
                          
                     else
                        ! RF added new term. send rest of N to plant - which decides whether it should pay or not? 
                        smin_nh4_to_plant_vr(c,j) = smin_nh4_vr(c,j)/dt - actual_immob_nh4_vr(c,j) - f_nit_vr(c,j) - f_nh3_vol_vr(c,j)
                     end if
                    
                  else
                     actual_immob_nh4_vr(c,j) = 0.0_r8
                     smin_nh4_to_plant_vr(c,j) = 0.0_r8
                     f_nit_vr(c,j) = 0.0_r8
                     f_nh3_vol_vr(c,j) = 0.0_r8
                  end if

                  if (potential_immob_vr(c,j) > 0.0_r8) then
                     fpi_nh4_vr(c,j) = actual_immob_nh4_vr(c,j) / potential_immob_vr(c,j)
                  else
                     fpi_nh4_vr(c,j) = 0.0_r8
                  end if

               end if


               ! Split realised NH3 volatilisation between canopy and free atmosphere
          
              
              
               if(.not.local_use_fun)then
                   sum_no3_demand(c,j) = (plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j)) + &
                  (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j)) + pot_f_denit_vr(c,j)
                   sum_no3_demand_scaled(c,j) = (plant_ndemand(c)*nuptake_prof(c,j) &
                                                 -smin_nh4_to_plant_vr(c,j))*compet_plant_no3 + &
                  (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))*compet_decomp_no3 + pot_f_denit_vr(c,j)*compet_denit
               else
                  sum_no3_demand(c,j) = plant_ndemand(c)*nuptake_prof(c,j) + &
                  (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j)) + pot_f_denit_vr(c,j)
                   sum_no3_demand_scaled(c,j) = (plant_ndemand(c)*nuptake_prof(c,j))*compet_plant_no3 + &
                  (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))*compet_decomp_no3 + pot_f_denit_vr(c,j)*compet_denit
               endif
                  
          

               if (sum_no3_demand(c,j)*dt < smin_no3_vr(c,j)) then

                  ! NO3 availability is not limiting immobilization or plant
                  ! uptake, and all can proceed at their potential rates
                  nlimit_no3(c,j) = 0
                  fpi_no3_vr(c,j) = 1.0_r8 -  fpi_nh4_vr(c,j)
                  actual_immob_no3_vr(c,j) = (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))

                  f_denit_vr(c,j) = pot_f_denit_vr(c,j)

                  if(.not.local_use_fun)then
                     smin_no3_to_plant_vr(c,j) = (plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))
                  else
                     ! This restricts the N uptake of a single layer to the value determined from the total demands and the 
                     ! hypothetical uptake profile above. Which is a strange thing to do, since that is independent of FUN
                     ! do we need this at all? 
                     smin_no3_to_plant_vr(c,j) = plant_ndemand(c)*nuptake_prof(c,j)
                     ! RF added new term. send rest of N to plant - which decides whether it should pay or not? 
                     if ( local_use_fun ) then
                        smin_no3_to_plant_vr(c,j) = smin_no3_vr(c,j)/dt - actual_immob_no3_vr(c,j) - f_denit_vr(c,j)
                     end if
                  endif
                
               else 

                  ! NO3 availability can not satisfy the sum of immobilization, denitrification, and
                  ! plant growth demands, so these three demands compete for available
                  ! soil mineral NO3 resource.
                  nlimit_no3(c,j) = 1
                                  
                  if (sum_no3_demand(c,j) > 0.0_r8) then
                     if(.not.local_use_fun)then
                        actual_immob_no3_vr(c,j) = min((smin_no3_vr(c,j)/dt)*((potential_immob_vr(c,j)- &
                        actual_immob_nh4_vr(c,j))*compet_decomp_no3 / sum_no3_demand_scaled(c,j)), &
                                  potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))
        
                        smin_no3_to_plant_vr(c,j) = min((smin_no3_vr(c,j)/dt)*((plant_ndemand(c)* &
                                  nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))*compet_plant_no3 / sum_no3_demand_scaled(c,j)), &
                                  plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))
        
                        f_denit_vr(c,j) = min((smin_no3_vr(c,j)/dt)*(pot_f_denit_vr(c,j)*compet_denit / &
                                  sum_no3_demand_scaled(c,j)), pot_f_denit_vr(c,j))
                     else
                        actual_immob_no3_vr(c,j) = min((smin_no3_vr(c,j)/dt)*((potential_immob_vr(c,j)- &
                        actual_immob_nh4_vr(c,j))*compet_decomp_no3 / sum_no3_demand_scaled(c,j)), &
                                  potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))

                        f_denit_vr(c,j) = min((smin_no3_vr(c,j)/dt)*(pot_f_denit_vr(c,j)*compet_denit / &
                        sum_no3_demand_scaled(c,j)), pot_f_denit_vr(c,j))
        
                        smin_no3_to_plant_vr(c,j) = (smin_no3_vr(c,j)/dt)*((plant_ndemand(c)* &
                                  nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))*compet_plant_no3 / sum_no3_demand_scaled(c,j))
                                  
                        ! RF added new term. send rest of N to plant - which decides whether it should pay or not? 
                        smin_no3_to_plant_vr(c,j) = (smin_no3_vr(c,j) / dt) - actual_immob_no3_vr(c,j) - f_denit_vr(c,j)
                        
  
                     end if ! use_fun

                  else ! no no3 demand. no uptake fluxes.
                     actual_immob_no3_vr(c,j) = 0.0_r8
                     smin_no3_to_plant_vr(c,j) = 0.0_r8
                     f_denit_vr(c,j) = 0.0_r8

                  end if !any no3 demand?
                  
                  
                  

                  if (potential_immob_vr(c,j) > 0.0_r8) then
                     fpi_no3_vr(c,j) = actual_immob_no3_vr(c,j) / potential_immob_vr(c,j)
                  else
                     fpi_no3_vr(c,j) = 0.0_r8
                  end if

               end if

               
                    

               !----------------------------------------------------------------
               ! Gaseous N products after competition.
               ! N2O remains pH-dependent as before; nitrification NOx uses
               ! Option 1 direct yield from F_NIT below.
               !
               ! At this point f_nit_vr and f_denit_vr are realised fluxes after
               ! NH4 and NO3 competition. This is the correct place to derive
               ! final NOx, rather than using potential fluxes upstream.
               !----------------------------------------------------------------

               ! pH used for the nitrification N2O loss fraction. The original
               ! Val Martin code uses soilph(c,j), with optional crop/enhanced-
               ! weathering pH perturbation. Your current tree does not expose
               ! soilph in this subroutine, so this keeps the same fixed Dahra
               ! test pH used elsewhere. Replace this with a real soil pH field
               ! when that plumbing exists.
               soilph_nox = 6.0_r8

               ! pH-dependent nitrification N2O fraction from Val Martin.
               fN2Onit_tmp = 721.86_r8 * exp(-2.387_r8 * soilph_nox)
               fN2Onit_tmp = min(1._r8, max(0._r8, fN2Onit_tmp))
               fN2Onit_vr(c,j) = fN2Onit_tmp

               ! N2O from realised nitrification and denitrification.
               f_n2o_nit_vr(c,j)   = f_nit_vr(c,j)   * fN2Onit_vr(c,j)
               f_n2o_denit_vr(c,j) = f_denit_vr(c,j) / (1._r8 + n2_n2o_ratio_denit_vr(c,j))

               !-------------------------------------------------------------
               ! Rain pulse for nitrification-derived NOx.
               !
               ! Val Martin logic:
               !   - ldry_vr is a dry-period counter in hours.
               !   - a rewetting increment of about 0.005 m3/m3 triggers a pulse.
               !   - dry periods shorter than about 72 h do not produce a pulse.
               !   - the pulse multiplier decays exponentially with 0.068 h-1.
               !-------------------------------------------------------------

               ! Initialise state fields defensively in case this is the first
               ! timestep after a restart or if the field still contains spval/nan.
               if (ldry_vr(c,j) /= ldry_vr(c,j) .or. ldry_vr(c,j) > 1.0e30_r8) ldry_vr(c,j) = 0._r8
               if (pfactor_vr(c,j) /= pfactor_vr(c,j) .or. pfactor_vr(c,j) > 1.0e30_r8) pfactor_vr(c,j) = 1._r8
               if (pfactor_vr(c,j) < 1._r8) pfactor_vr(c,j) = 1._r8

               ! Change in volumetric soil moisture since the previous timestep.
               h2osoi_diff = h2osoi_vol(c,j) - prev_h2osoi_vol_vr(c,j)

               if (h2osoi_diff > nox_rewet_dtheta_crit) then

                  ! Rewetting event. A pulse is produced only if the preceding
                  ! dry period is long enough.
                  if (ldry_vr(c,j) > nox_min_dry_hours) then
                     nox_pulse_raw = nox_pulse_a * log(max(ldry_vr(c,j), 1._r8)) - nox_pulse_b
                     pfactor_vr(c,j) = max(1._r8, nox_pulse_raw)
                  else
                     pfactor_vr(c,j) = 1._r8
                  end if

                  ! The dry spell ends once rewetting occurs.
                  ldry_vr(c,j) = 0._r8

               else

                  ! No new rewetting event. Accumulate dry hours only when
                  ! volumetric soil water is below the Val Martin/Hudman dry threshold.
                  if (h2osoi_vol(c,j) < nox_dry_theta_crit) then
                     ldry_vr(c,j) = ldry_vr(c,j) + dthr
                  else
                     ldry_vr(c,j) = 0._r8
                  end if

                  ! Continue decay of any active pulse.
                  if (pfactor_vr(c,j) > 1._r8) then
                     pfactor_vr(c,j) = pfactor_vr(c,j) * exp(-nox_pulse_decay_hr * dthr)
                     if (pfactor_vr(c,j) <= 1._r8 + nox_pulse_min) pfactor_vr(c,j) = 1._r8
                  end if

               end if

               ! Store current soil water for the next timestep's rewetting test.
               prev_h2osoi_vol_vr(c,j) = h2osoi_vol(c,j)

               !-------------------------------------------------------------
               ! Parton/DAYCENT-inspired separate NOx gas submodel.
               !
               ! This block avoids the Val Martin failure mode in which
               ! nitrification NOx is throttled by the very small nitrification
               ! N2O yield:
               !
               !   old: F_NOx_NIT = F_N2O_NIT * NOx_N2O_RATIO * PULSE_FAC
               !
               ! Here, NOx is a separate product from the realised process fluxes,
               ! with environmental yield functions:
               !
               !   F_NOx_NIT   = F_NIT   * Y_NO_NIT(WFPS, Dr) * PULSE_FAC
               !   F_NOx_DENIT = F_DENIT * Y_NO_DENIT(WFPS, Dr)
               !
               ! WFPS sets the moisture window.  Dr is a simple relative gas
               ! diffusivity proxy from air-filled pore space.  This means that
               ! large yields occur only when the process flux and the physical
               ! environment are co-located.
               !-------------------------------------------------------------

               if (watsat(c,j) > 0._r8) then
                  wfps_frac_nox = max(0._r8, min(1._r8, h2osoi_vol(c,j) / watsat(c,j)))
               else
                  wfps_frac_nox = 0._r8
               end if

               ! Air-filled pore space and a simple relative diffusivity proxy.
               ! This follows the same broad logic as Parton/Val Martin: NO is
               ! favoured when gas diffusion is not strongly restricted.
               afps_frac_nox = max(0._r8, 1._r8 - wfps_frac_nox)
               dr_nox        = 0.209_r8 * afps_frac_nox**(4._r8/3._r8)

               ! Moisture response for nitrification NO: peak at moderately
               ! aerated soil.  A small floor preserves microsite/background
               ! production, but F_NIT still controls the magnitude.
               nox_nit_moist_fac = exp(-((wfps_frac_nox - nox_nit_wfps_opt) / &
                                          nox_nit_wfps_width)**2)
               nox_nit_moist_fac = max(nox_nit_moist_fac_min, &
                                         min(1._r8, nox_nit_moist_fac))

               ! Diffusivity response for nitrification NO.  This limits NO
               ! production under poorly aerated conditions without using N2O as
               ! the parent flux.
               if (nox_nit_dr_scale > 0._r8) then
                  nox_nit_diff_fac = max(0._r8, min(1._r8, dr_nox / nox_nit_dr_scale))
               else
                  nox_nit_diff_fac = 1._r8
               end if

               ! Realised nitrification NO yield.  y_no_nit_max is a calibration
               ! ceiling; the realised yield is environmentally reduced.
               y_no_nit = y_no_nit_max * nox_nit_moist_fac * nox_nit_diff_fac

               ! Base nitrification NOx before pulse.
               nox_nit_before_pulse = max(0._r8, f_nit_vr(c,j) * y_no_nit)

               ! Co-located pulse enhancement.  A large pfactor only matters if
               ! realised nitrification and the environmental yield are active in
               ! the same layer and timestep.
               nox_nit_after_pulse = nox_nit_before_pulse * pfactor_vr(c,j)

               pot_f_nox_nit_vr(c,j) = nox_nit_after_pulse
               f_nox_nit_vr(c,j)     = max(0._r8, nox_nit_after_pulse)

               ! Effective column pulse diagnostic.  This reports the actual
               ! NOx-weighted pulse enhancement, not the maximum pfactor in any
               ! layer.
               pulse_weighted_num(c) = pulse_weighted_num(c) + nox_nit_after_pulse  * dzsoi_decomp(j)
               pulse_weighted_den(c) = pulse_weighted_den(c) + nox_nit_before_pulse * dzsoi_decomp(j)

               if (pulse_weighted_den(c) > 0._r8) then
                  pulse_fac_col(c) = pulse_weighted_num(c) / pulse_weighted_den(c)
               else
                  pulse_fac_col(c) = 1._r8
               end if

               ! Denitrification NOx from realised denitrification using the
               ! same source philosophy as nitrification.  The moisture response
               ! peaks at wetter transition conditions, but diffusivity suppresses
               ! NO when gas diffusion becomes too restricted.
               nox_denit_moist_fac = exp(-((wfps_frac_nox - nox_denit_wfps_opt) / &
                                           nox_denit_wfps_width)**2)
               nox_denit_moist_fac = max(nox_denit_moist_fac_min, &
                                          min(1._r8, nox_denit_moist_fac))

               if (nox_denit_dr_scale > 0._r8) then
                  nox_denit_diff_fac = max(0._r8, min(1._r8, dr_nox / nox_denit_dr_scale))
               else
                  nox_denit_diff_fac = 1._r8
               end if

               ! Realised denitrification NO yield.  y_no_denit_max is a
               ! calibration ceiling; the realised yield is reduced by moisture
               ! and gas diffusivity.
               y_no_denit = y_no_denit_max * nox_denit_moist_fac * nox_denit_diff_fac

               nox_denit_direct = max(0._r8, f_denit_vr(c,j) * y_no_denit)

               pot_f_nox_denit_vr(c,j) = nox_denit_direct
               f_nox_denit_vr(c,j)     = nox_denit_direct

               ! Canopy-reduced / above-canopy NOx. This follows Val Martin:
               ! final capped NOx is multiplied by the column-mapped dry
               ! deposition canopy reduction factor.
               f_nox_nit_atmos_vr(c,j)   = f_nox_nit_vr(c,j)   * crf_drydep_col(c)
               f_nox_denit_atmos_vr(c,j) = f_nox_denit_vr(c,j) * crf_drydep_col(c)

               ! Column-integrated diagnostics [gN/m2/s].
               f_nox_nit_col(c)         = f_nox_nit_col(c)         + f_nox_nit_vr(c,j)         * dzsoi_decomp(j)
               f_nox_denit_col(c)       = f_nox_denit_col(c)       + f_nox_denit_vr(c,j)       * dzsoi_decomp(j)
               f_nox_nit_atmos_col(c)   = f_nox_nit_atmos_col(c)   + f_nox_nit_atmos_vr(c,j)   * dzsoi_decomp(j)
               f_nox_denit_atmos_col(c) = f_nox_denit_atmos_col(c) + f_nox_denit_atmos_vr(c,j) * dzsoi_decomp(j)
               pot_f_nox_nit_col(c)     = pot_f_nox_nit_col(c)     + pot_f_nox_nit_vr(c,j)     * dzsoi_decomp(j)
               pot_f_nox_denit_col(c)   = pot_f_nox_denit_col(c)   + pot_f_nox_denit_vr(c,j)   * dzsoi_decomp(j)

               ! Total soil NOx to atmosphere.
               soil_nox_total_col(c) = f_nox_nit_atmos_col(c) + f_nox_denit_atmos_col(c)
               soil_nox_crop_col(c)  = soil_nox_total_col(c)


               ! this code block controls the addition of N to sminn pool
               ! to eliminate any N limitation, when Carbon_Only is set.  This lets the
               ! model behave essentially as a carbon-only model, but with the
               ! benefit of keeping track of the N additions needed to
               ! eliminate N limitations, so there is still a diagnostic quantity
               ! that describes the degree of N limitation at steady-state.

               if ( cnallocate_carbon_only()) then !.or. &
                  if ( fpi_no3_vr(c,j) + fpi_nh4_vr(c,j) < 1._r8 ) then
                     fpi_nh4_vr(c,j) = 1.0_r8 - fpi_no3_vr(c,j)
                     supplement_to_sminn_vr(c,j) = (potential_immob_vr(c,j) &
                                                  - actual_immob_no3_vr(c,j)) - actual_immob_nh4_vr(c,j)
                     ! update to new values that satisfy demand
                     actual_immob_nh4_vr(c,j) = potential_immob_vr(c,j) -  actual_immob_no3_vr(c,j)   
                  end if
                  if ( smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j) < plant_ndemand(c)*nuptake_prof(c,j) ) then
                     supplement_to_sminn_vr(c,j) = supplement_to_sminn_vr(c,j) + &
                          (plant_ndemand(c)*nuptake_prof(c,j) - smin_no3_to_plant_vr(c,j)) - smin_nh4_to_plant_vr(c,j)  ! use old values
                     smin_nh4_to_plant_vr(c,j) = plant_ndemand(c)*nuptake_prof(c,j) - smin_no3_to_plant_vr(c,j)
                  end if
                  sminn_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)
               end if

               ! sum up no3 and nh4 fluxes
               fpi_vr(c,j) = fpi_no3_vr(c,j) + fpi_nh4_vr(c,j)
               sminn_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)
               actual_immob_vr(c,j) = actual_immob_no3_vr(c,j) + actual_immob_nh4_vr(c,j)
            end do
         end do

         if ( local_use_fun ) then
            call t_startf( 'CNFUN' )
            call CNFUN(bounds,num_soilc,filter_soilc,num_soilp,filter_soilp,waterstate_inst                 ,&
                      waterflux_inst,temperature_inst,soilstate_inst,cnveg_state_inst,cnveg_carbonstate_inst,&
                      cnveg_carbonflux_inst,cnveg_nitrogenstate_inst,cnveg_nitrogenflux_inst                ,&
                      soilbiogeochem_nitrogenflux_inst,soilbiogeochem_carbonflux_inst,canopystate_inst,      &
                      soilbiogeochem_nitrogenstate_inst)
                      
            ! sminn_to_plant_fun is output of actual N uptake from FUN
            call p2c_2d(bounds,nlevdecomp, &
                       cnveg_nitrogenflux_inst%sminn_to_plant_fun_no3_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
                       soilbiogeochem_nitrogenflux_inst%sminn_to_plant_fun_no3_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
                       'unity')

            call p2c_2d(bounds,nlevdecomp, &
                       cnveg_nitrogenflux_inst%sminn_to_plant_fun_nh4_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
                       soilbiogeochem_nitrogenflux_inst%sminn_to_plant_fun_nh4_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
                       'unity')
            call t_stopf( 'CNFUN' )
         end if



         if(.not.local_use_fun)then
            do fc=1,num_soilc
               c = filter_soilc(fc)
               ! sum up N fluxes to plant after initial competition
               sminn_to_plant(c) = 0._r8
            end do
            do j = 1, nlevdecomp  
               do fc=1,num_soilc
                  c = filter_soilc(fc)
                  sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
               end do
            end do
         else
            do fc=1,num_soilc
               c = filter_soilc(fc)
               ! sum up N fluxes to plant after initial competition
               sminn_to_plant(c) = 0._r8 !this isn't use in fun. 
               do j = 1, nlevdecomp
                  if ((sminn_to_plant_fun_no3_vr(c,j)-smin_no3_to_plant_vr(c,j)).gt.0.0000000000001_r8) then
                      write(iulog,*) 'problem with limitations on no3 uptake', &
                                 sminn_to_plant_fun_no3_vr(c,j),smin_no3_to_plant_vr(c,j)
                      call endrun("too much NO3 uptake predicted by FUN")
                  end if
!KO                  if ((sminn_to_plant_fun_nh4_vr(c,j)-smin_nh4_to_plant_vr(c,j)).gt.0.0000000000001_r8) then
!KO
                  if ((sminn_to_plant_fun_nh4_vr(c,j)-smin_nh4_to_plant_vr(c,j)).gt.0.0000001_r8) then
!KO
                      write(iulog,*) 'problem with limitations on nh4 uptake', &
                                  sminn_to_plant_fun_nh4_vr(c,j),smin_nh4_to_plant_vr(c,j)
                      call endrun("too much NH4 uptake predicted by FUN")
                  end if
               end do
            end do

         end if

         if(.not.local_use_fun)then
            ! give plants a second pass to see if there is any mineral N left over with which to satisfy residual N demand.
            ! first take frm nh4 pool; then take from no3 pool
            do fc=1,num_soilc
               c = filter_soilc(fc)
               residual_plant_ndemand(c) = plant_ndemand(c) - sminn_to_plant(c)
               residual_smin_nh4(c) = 0._r8
            end do
            do j = 1, nlevdecomp  
               do fc=1,num_soilc
                  c = filter_soilc(fc)
                  if (residual_plant_ndemand(c)  >  0._r8 ) then
                     if (nlimit_nh4(c,j) .eq. 0) then
                        residual_smin_nh4_vr(c,j) = max(smin_nh4_vr(c,j) - (actual_immob_nh4_vr(c,j) + &
                                                    smin_nh4_to_plant_vr(c,j) + f_nit_vr(c,j) ) * dt, 0._r8)

                        residual_smin_nh4(c) = residual_smin_nh4(c) + residual_smin_nh4_vr(c,j) * dzsoi_decomp(j)
                     else
                        residual_smin_nh4_vr(c,j)  = 0._r8
                     endif
   
                     if ( residual_smin_nh4(c) > 0._r8 .and. nlimit_nh4(c,j) .eq. 0 ) then
                        smin_nh4_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + residual_smin_nh4_vr(c,j) * &
                             min(( residual_plant_ndemand(c) *  dt ) / residual_smin_nh4(c), 1._r8) / dt
                     endif
                  end if
               end do
            end do

            ! re-sum up N fluxes to plant after second pass for nh4
            do fc=1,num_soilc
               c = filter_soilc(fc)
               sminn_to_plant(c) = 0._r8
            end do
            do j = 1, nlevdecomp
               do fc=1,num_soilc
                  c = filter_soilc(fc)
                  sminn_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + smin_no3_to_plant_vr(c,j)
                  sminn_to_plant(c) = sminn_to_plant(c) + (sminn_to_plant_vr(c,j)) * dzsoi_decomp(j)
               end do
            end do

            !
            ! and now do second pass for no3
            do fc=1,num_soilc
               c = filter_soilc(fc)
               residual_plant_ndemand(c) = plant_ndemand(c) - sminn_to_plant(c)
               residual_smin_no3(c) = 0._r8
            end do

            do j = 1, nlevdecomp
               do fc=1,num_soilc
                  c = filter_soilc(fc)
                  if (residual_plant_ndemand(c) > 0._r8 ) then
                     if (nlimit_no3(c,j) .eq. 0) then
                       residual_smin_no3_vr(c,j) = max(smin_no3_vr(c,j) - (actual_immob_no3_vr(c,j) + &
                                                   smin_no3_to_plant_vr(c,j) + f_denit_vr(c,j) ) * dt, 0._r8)
                        residual_smin_no3(c) = residual_smin_no3(c) + residual_smin_no3_vr(c,j) * dzsoi_decomp(j)
                     else
                        residual_smin_no3_vr(c,j)  = 0._r8
                     endif
   
                     if ( residual_smin_no3(c) > 0._r8 .and. nlimit_no3(c,j) .eq. 0) then
                        smin_no3_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + residual_smin_no3_vr(c,j) * &
                             min(( residual_plant_ndemand(c) *  dt ) / residual_smin_no3(c), 1._r8) / dt
                     endif
                  endif
               end do
            end do

            ! re-sum up N fluxes to plant after second passes of both no3 and nh4
            do fc=1,num_soilc
               c = filter_soilc(fc)
               sminn_to_plant(c) = 0._r8
            end do
            do j = 1, nlevdecomp
               do fc=1,num_soilc
                  c = filter_soilc(fc)
                  sminn_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + smin_no3_to_plant_vr(c,j)
                  sminn_to_plant(c) = sminn_to_plant(c) + (sminn_to_plant_vr(c,j)) * dzsoi_decomp(j)
               end do
            end do
   
         else !use_fun
         !calculate maximum N available to plants. 
            do fc=1,num_soilc
               c = filter_soilc(fc)
               sminn_to_plant(c) = 0._r8
            end do
            do j = 1, nlevdecomp
               do fc=1,num_soilc
                  c = filter_soilc(fc)
                  sminn_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + smin_no3_to_plant_vr(c,j)
                  sminn_to_plant(c) = sminn_to_plant(c) + (sminn_to_plant_vr(c,j)) * dzsoi_decomp(j)
               end do
            end do
   

             ! add up fun fluxes from SMINN to plant. 
             do j = 1, nlevdecomp
                do fc=1,num_soilc
                   c = filter_soilc(fc)
                   sminn_to_plant_new(c)  = sminn_to_plant_new(c) + &
                             (sminn_to_plant_fun_no3_vr(c,j) + sminn_to_plant_fun_nh4_vr(c,j)) * dzsoi_decomp(j)
                      
                end do
             end do
                             
                              
         end if !use_f
         ! sum up N fluxes to immobilization
         do fc=1,num_soilc
            c = filter_soilc(fc)
            actual_immob(c) = 0._r8
            potential_immob(c) = 0._r8
         end do
         do j = 1, nlevdecomp  
            do fc=1,num_soilc
               c = filter_soilc(fc)
               actual_immob(c) = actual_immob(c) + actual_immob_vr(c,j) * dzsoi_decomp(j)
               potential_immob(c) = potential_immob(c) + potential_immob_vr(c,j) * dzsoi_decomp(j)
            end do
         end do
        
        
     
       
         do fc=1,num_soilc
            c = filter_soilc(fc)   
            ! calculate the fraction of potential growth that can be
            ! acheived with the N available to plants
            ! calculate the fraction of immobilization realized (for diagnostic purposes)
            if(.not.local_use_fun)then !FUN has no concept of FPG.
            
               if (plant_ndemand(c) > 0.0_r8) then
                  fpg(c) = sminn_to_plant(c) / plant_ndemand(c)
               else
                  fpg(c) = 1._r8
               end if
            end if

            if (potential_immob(c) > 0.0_r8) then
               fpi(c) = actual_immob(c) / potential_immob(c)
            else
               fpi(c) = 1._r8
            end if
         end do ! end of column loops

      end if  !end of if_not_use_nitrif_denitrif


     ! ===== BEG: NH3 canopy capture (Fung/DNDC-style, patch-weighted) =====
     ! Split soil NH3 volatilisation into canopy-captured and free-air components
     ! using friction velocity, canopy geometry, LAI, and near-leaf humidity.

     ! Reset layer diagnostics before aggregating patch contributions
     do j = 1, nlevdecomp
        do fc = 1, num_soilc
           c = filter_soilc(fc)
           f_nh3_vol_to_canopy_vr(c,j) = 0.0_r8
           f_nh3_vol_to_atmos_vr(c,j)  = 0.0_r8
        end do
     end do

     do j = 1, nlevdecomp
        do fp = 1, num_soilp
           p = filter_soilp(fp)
           c = patch%column(p)

           canopy_top     = max(zlnd, canopystate_inst%htop_patch(p))
           canopy_bot     = max(zlnd, canopystate_inst%hbot_patch(p))
           wind_canopy    = max(0.001_r8, frictionvel_inst%fv_patch(p))

           nh3_conc          = f_nh3_vol_vr(c,j) * patch%wtcol(p)
           fac_canopy_height = 14.0_r8 * (canopy_top - canopy_bot)
           fac_leaf_geometry = max(0.0_r8, canopystate_inst%tlai_patch(p))
           fac_leaf_moisture = min(1.0_r8, max(0.0_r8, waterstate_inst%rh_af_patch(p)))
           vg_nh3            = 0.05_r8

           nh3_vol_to_canopy_vr_patch = nh3_conc * (fac_canopy_height * fac_leaf_geometry * fac_leaf_moisture / wind_canopy * vg_nh3)
           nh3_vol_to_canopy_vr_patch = min( f_nh3_vol_vr(c,j) * patch%wtcol(p), max(0.0_r8, nh3_vol_to_canopy_vr_patch) )

           f_nh3_vol_to_canopy_vr(c,j) = f_nh3_vol_to_canopy_vr(c,j) + nh3_vol_to_canopy_vr_patch
        end do

        do fc = 1, num_soilc
           c = filter_soilc(fc)
           f_nh3_vol_to_atmos_vr(c,j) = max(0.0_r8, f_nh3_vol_vr(c,j) - f_nh3_vol_to_canopy_vr(c,j))
        end do
     end do
     ! ===== END: NH3 canopy capture =====

    end associate

  end subroutine SoilBiogeochemCompetition

end module SoilBiogeochemCompetitionMod
