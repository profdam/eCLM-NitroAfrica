module SoilBiogeochemNitrifDenitrifMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate nitrification and denitrification rates.
  !
  ! NitroAfrica / Dahra modification in this version:
  !   - This module computes the nitrification/denitrification controls and the
  !     NOx:N2O ratio following the Val Martin implementation.
  !   - Final NOx flux construction, rain-pulse multiplication, capping and
  !     canopy reduction are handled in SoilBiogeochemCompetitionMod.F90 because
  !     Val Martin derives NOx after realised f_nit_vr and f_denit_vr are known.
  !
  ! Notes on units used in the pulse block:
  !   - h2osoi_vol is volumetric water content, m3 water / m3 soil.
  !   - dry threshold is 0.175 m3 m-3, equivalent to 17.5 percent volumetric water.
  !   - antecedent dry period is accumulated in days.
  !   - the Val Martin / Hudman pulse expression uses dry period in hours.
  !   - the decay coefficient is 0.068 h-1, equivalent to 1.632 d-1.
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use shr_const_mod                   , only : SHR_CONST_TKFRZ
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use clm_varpar                      , only : nlevdecomp
  use clm_varcon                      , only : rpi, grav, spval
  use clm_varcon                      , only : d_con_g, d_con_w, secspday
  use clm_varctl                      , only : use_lch4
  use abortutils                      , only : endrun
  use decompMod                       , only : bounds_type
  use SoilStatetype                   , only : soilstate_type
  use WaterStateType                  , only : waterstate_type
  use TemperatureType                 , only : temperature_type
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use ch4Mod                          , only : ch4_type
  use ColumnType                      , only : col                
  use CropType                        , only : crop_type

  implicit none
  private

  public :: readParams
  public :: nitrifReadNML
  public :: SoilBiogeochemNitrifDenitrif

  type, private :: params_type
     real(r8) :: k_nitr_max
     real(r8) :: surface_tension_water
     real(r8) :: rij_kro_a
     real(r8) :: rij_kro_alpha
     real(r8) :: rij_kro_beta
     real(r8) :: rij_kro_gamma
     real(r8) :: rij_kro_delta
     real(r8) :: denitrif_respiration_coefficient
     real(r8) :: denitrif_respiration_exponent
     real(r8) :: denitrif_nitrateconc_coefficient
     real(r8) :: denitrif_nitrateconc_exponent

     ! NitroAfrica NH3 volatilisation parameters.
     real(r8) :: k_nh3_vol_max
     real(r8) :: nh3_q10
     real(r8) :: nh3_tref
     real(r8) :: nh3_wfps_opt
     real(r8) :: nh3_wfps_width
     real(r8) :: nh3_pKa

     ! NitroAfrica placeholders.
     real(r8), pointer :: pot_f_nit_naf(:,:)             => null()
     real(r8), pointer :: pot_f_denit_naf(:,:)           => null()
     real(r8), pointer :: soil_n2o_naf_col(:)            => null()
     real(r8), pointer :: soil_n2o_naf_crop_col(:)       => null()
     real(r8), pointer :: f_n2o_nit_naf_atmos_patch(:)   => null()
     real(r8), pointer :: f_n2o_denit_naf_atmos_patch(:) => null()
  end type params_type

  type(params_type), private :: params_inst

  !-----------------------------------------------------------------------
  ! Persistent rain-pulse state arrays.
  !
  ! These arrays are saved between calls because the rain pulse depends on the
  ! previous soil-moisture state and on the number of dry days before rewetting.
  ! They are stored by column.
  !
  ! antecedent_dry_days(c):
  !   Number of consecutive dry days before the current timestep.
  !
  ! pulse_dry_days(c):
  !   Dry-period length saved at the exact timestep when rewetting triggers
  !   a pulse. This is important because antecedent_dry_days is reset after
  !   the pulse starts.
  !
  ! pulse_elapsed_days(c):
  !   Time since the pulse was triggered.
  !
  ! prev_theta_top(c):
  !   Top-layer volumetric water content at the previous timestep.
  !
  ! pulse_active(c):
  !   Logical flag indicating whether a pulse is currently decaying.
  !-----------------------------------------------------------------------
  real(r8), allocatable, save :: antecedent_dry_days(:)
  real(r8), allocatable, save :: pulse_dry_days(:)
  real(r8), allocatable, save :: pulse_elapsed_days(:)
  real(r8), allocatable, save :: prev_theta_top(:)
  logical , allocatable, save :: pulse_active(:)
  logical , save :: pulse_arrays_initialized = .false.

  logical, public :: no_frozen_nitrif_denitrif = .false.

  character(len=*), parameter, private :: sourcefile = __FILE__

contains

  !-----------------------------------------------------------------------  
  subroutine readParams ( ncid )

    use ncdio_pio, only: file_desc_t,ncd_io

    type(file_desc_t),intent(inout) :: ncid

    character(len=32)  :: subname = 'CNNitrifDenitrifParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv
    real(r8)           :: tempr
    character(len=100) :: tString

    tString='surface_tension_water'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%surface_tension_water=tempr

    tString='rij_kro_a'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rij_kro_a=tempr

    tString='rij_kro_alpha'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rij_kro_alpha=tempr

    tString='rij_kro_beta'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rij_kro_beta=tempr

    tString='rij_kro_gamma'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rij_kro_gamma=tempr

    tString='rij_kro_delta'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rij_kro_delta=tempr

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine nitrifReadNML( NLFilename )

    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog

    character(len=*), intent(in) :: NLFilename

    integer :: ierr
    integer :: unitn

    character(len=*), parameter :: subname = 'ReadNML'
    character(len=*), parameter :: nmlname = 'nitrif_inparm'

    real(r8) :: k_nitr_max_perday, denitrif_respiration_coefficient, &
                denitrif_respiration_exponent, denitrif_nitrateconc_coefficient, &
                denitrif_nitrateconc_exponent

    real(r8) :: k_nh3_vol_max_perday, nh3_q10, nh3_tref, nh3_wfps_opt, nh3_wfps_width, nh3_pka

    namelist /nitrif_inparm/ k_nitr_max_perday, denitrif_respiration_coefficient, &
             denitrif_respiration_exponent, denitrif_nitrateconc_coefficient, &
             denitrif_nitrateconc_exponent, k_nh3_vol_max_perday, nh3_q10, &
             nh3_tref, nh3_wfps_opt, nh3_wfps_width, nh3_pka

    denitrif_respiration_coefficient = 0.1_r8
    denitrif_respiration_exponent    = 1.3_r8
    denitrif_nitrateconc_coefficient = 1.15_r8
    denitrif_nitrateconc_exponent    = 0.57_r8
    k_nitr_max_perday                = 0.1_r8

    k_nh3_vol_max_perday = 0.05_r8
    nh3_q10              = 2.0_r8
    nh3_tref             = 298.15_r8
    nh3_wfps_opt         = 60.0_r8
    nh3_wfps_width       = 20.0_r8
    nh3_pka              = 9.25_r8

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=nitrif_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (k_nitr_max_perday                , mpicom)
    call shr_mpi_bcast (denitrif_respiration_coefficient , mpicom)
    call shr_mpi_bcast (denitrif_respiration_exponent    , mpicom)
    call shr_mpi_bcast (denitrif_nitrateconc_coefficient , mpicom)
    call shr_mpi_bcast (denitrif_nitrateconc_exponent    , mpicom)

    call shr_mpi_bcast (k_nh3_vol_max_perday, mpicom)
    call shr_mpi_bcast (nh3_q10            , mpicom)
    call shr_mpi_bcast (nh3_tref           , mpicom)
    call shr_mpi_bcast (nh3_wfps_opt       , mpicom)
    call shr_mpi_bcast (nh3_wfps_width     , mpicom)
    call shr_mpi_bcast (nh3_pka            , mpicom)

    params_inst%k_nitr_max =  k_nitr_max_perday / secspday
    params_inst%denitrif_respiration_coefficient = denitrif_respiration_coefficient
    params_inst%denitrif_respiration_exponent    = denitrif_respiration_exponent
    params_inst%denitrif_nitrateconc_coefficient = denitrif_nitrateconc_coefficient
    params_inst%denitrif_nitrateconc_exponent    = denitrif_nitrateconc_exponent

    params_inst%k_nh3_vol_max = k_nh3_vol_max_perday / secspday
    params_inst%nh3_q10       = nh3_q10
    params_inst%nh3_tref      = nh3_tref
    params_inst%nh3_wfps_opt  = nh3_wfps_opt
    params_inst%nh3_wfps_width= nh3_wfps_width
    params_inst%nh3_pKa       = nh3_pka

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=nitrif_inparm)
       write(iulog,*) ' '
    end if

  end subroutine nitrifReadNML

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemNitrifDenitrif(bounds, num_soilc, filter_soilc, &
       soilstate_inst, waterstate_inst, temperature_inst, ch4_inst, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst, &
       soilbiogeochem_nitrogenflux_inst, canopy_red_fac_col, crop_inst)

    use clm_time_manager  , only : get_step_size
    use CNSharedParamsMod , only : anoxia_wtsat, CNParamsShareInst
    use landunit_varcon   , only : istcrop
    use LandunitType      , only : lun
    use subgridAveMod     , only : p2c

    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc
    integer                                 , intent(in)    :: filter_soilc(:)
    type(soilstate_type)                    , intent(in)    :: soilstate_inst
    type(waterstate_type)                   , intent(in)    :: waterstate_inst
    type(temperature_type)                  , intent(in)    :: temperature_inst
    type(ch4_type)                          , intent(in)    :: ch4_inst
    type(soilbiogeochem_carbonflux_type)    , intent(in)    :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(in)    :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    real(r8), intent(in), optional :: canopy_red_fac_col(bounds%begc:bounds%endc)
    type(crop_type), intent(inout), optional :: crop_inst

    integer  :: c, fc, j
    real(r8) :: soil_hr_vr(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: g_per_m3__to__ug_per_gsoil
    real(r8) :: g_per_m3_sec__to__ug_per_gsoil_day
    real(r8) :: co2diff_con(2)
    real(r8) :: eps, f_a
    real(r8) :: surface_tension_water
    real(r8) :: rij_kro_a, rij_kro_alpha, rij_kro_beta, rij_kro_gamma, rij_kro_delta
    real(r8) :: rho_w  = 1.e3_r8
    real(r8) :: r_max
    real(r8) :: r_min(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: ratio_diffusivity_water_gas(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: om_frac
    real(r8) :: anaerobic_frac_sat, r_psi_sat, r_min_sat
    real(r8) :: organic_max
    character(len=32) :: subname='nitrif_denitrif'

    real(r8) :: f_T_nh3, f_wfps_nh3, f_pH_nh3
    real(r8) :: wfps_loc, pH_loc
    real(r8) :: k1_nitr
    real(r8) :: maxrate_perday, maxrate
    real(r8) :: base_flux_perday, base_flux
    real(r8) :: nreduce
    real(r8) :: min_ammonium_perday, min_ammonium
         
    real(r8) :: fr_ph_denit
    real(r8) :: dt
    real(r8) :: Dr
    real(r8) :: fno3_co2
    real(r8) :: pot_f_denit_vm, n2_n2o_ratio_denit_vm
    real(r8) :: f_n2o_nit_frac_vm
         
    real(r8) :: pot_f_nit_vm
    real(r8) :: f_n2o_nit_layer_vm
    real(r8) :: f_n2o_denit_layer_vm
         
    real(r8) :: f_n2o_nit_layer, f_n2o_denit_layer
    real(r8) :: f_nox_nit_layer, f_nox_denit_layer

    real(r8) :: pH(bounds%begc:bounds%endc)
    real(r8) :: fac
    real(r8) :: fcan

    !--------------------------------------------------------------------
    ! Rain-pulse local variables.
    !--------------------------------------------------------------------
    real(r8) :: theta_now
    real(r8) :: delta_theta
    real(r8) :: dry_days_before_update
    real(r8) :: l_dry_days
    real(r8) :: l_dry_hours
    real(r8) :: t_pulse_days
    real(r8) :: pulse_raw_peak
    real(r8) :: pulse_term
    real(r8) :: pulse_fac
    logical  :: pulse_trigger

    !--------------------------------------------------------------------
    ! Rain-pulse constants.
    !
    ! These are deliberately declared as local constants to make sensitivity
    ! tests easy. For Dahra calibration, start by changing pulse_a and pulse_b
    ! only. Keep the thresholds fixed while testing a and b.
    !--------------------------------------------------------------------
    real(r8), parameter :: pulse_a                 = 13.01_r8
    real(r8), parameter :: pulse_b                 = 53.60_r8
    real(r8), parameter :: pulse_decay_per_day     = 1.632_r8
    real(r8), parameter :: pulse_dry_theta_crit    = 0.175_r8
    real(r8), parameter :: pulse_rewet_dtheta_crit = 1.0e-3_r8
    real(r8), parameter :: pulse_min_dry_days      = 1.0_r8
    real(r8), parameter :: pulse_term_min          = 1.0e-6_r8

    ! Model timestep in seconds.
    dt = real(get_step_size(), r8)

    !--------------------------------------------------------------------
    ! Allocate persistent pulse arrays once.
    !
    ! This block is called on the first entry into the subroutine. It stores
    ! state by column and then reuses that state at later timesteps.
    !--------------------------------------------------------------------
    if (.not. pulse_arrays_initialized) then
       allocate(antecedent_dry_days(bounds%begc:bounds%endc))
       allocate(pulse_dry_days(bounds%begc:bounds%endc))
       allocate(pulse_elapsed_days(bounds%begc:bounds%endc))
       allocate(prev_theta_top(bounds%begc:bounds%endc))
       allocate(pulse_active(bounds%begc:bounds%endc))

       antecedent_dry_days(bounds%begc:bounds%endc) = 0._r8
       pulse_dry_days(bounds%begc:bounds%endc)      = 0._r8
       pulse_elapsed_days(bounds%begc:bounds%endc)  = 0._r8
       prev_theta_top(bounds%begc:bounds%endc)      = &
            waterstate_inst%h2osoi_vol_col(bounds%begc:bounds%endc,1)
       pulse_active(bounds%begc:bounds%endc)        = .false.

       pulse_arrays_initialized = .true.
    end if

    !--------------------------------------------------------------------
    ! Reset column-integrated diagnostic fluxes for the current timestep.
    !--------------------------------------------------------------------
    soilbiogeochem_nitrogenflux_inst%f_n2o_nit_vm_col        (bounds%begc:bounds%endc) = 0._r8
    soilbiogeochem_nitrogenflux_inst%f_n2o_denit_vm_col      (bounds%begc:bounds%endc) = 0._r8
    soilbiogeochem_nitrogenflux_inst%f_nox_nit_col           (bounds%begc:bounds%endc) = 0._r8
    soilbiogeochem_nitrogenflux_inst%f_nox_denit_col         (bounds%begc:bounds%endc) = 0._r8
    soilbiogeochem_nitrogenflux_inst%f_nox_nit_atmos_col     (bounds%begc:bounds%endc) = 0._r8
    soilbiogeochem_nitrogenflux_inst%f_nox_denit_atmos_col   (bounds%begc:bounds%endc) = 0._r8
    soilbiogeochem_nitrogenflux_inst%pot_f_nox_nit_col       (bounds%begc:bounds%endc) = 0._r8
    soilbiogeochem_nitrogenflux_inst%pot_f_nox_denit_col     (bounds%begc:bounds%endc) = 0._r8

    ! PULSE_FAC history variable. Default is 1, meaning no pulse enhancement.
    soilbiogeochem_nitrogenflux_inst%pulse_fac_col           (bounds%begc:bounds%endc) = 1._r8

    associate(                                                                      &
         watsat                        => soilstate_inst%watsat_col               , &
         watfc                         => soilstate_inst%watfc_col                , &
         bd                            => soilstate_inst%bd_col                   , &
         bsw                           => soilstate_inst%bsw_col                  , &
         cellorg                       => soilstate_inst%cellorg_col              , &
         sucsat                        => soilstate_inst%sucsat_col               , &
         soilpsi                       => soilstate_inst%soilpsi_col              , &
         h2osoi_vol                    => waterstate_inst%h2osoi_vol_col          , &
         h2osoi_liq                    => waterstate_inst%h2osoi_liq_col          , &
         t_soisno                      => temperature_inst%t_soisno_col           , &
         o2_decomp_depth_unsat         => ch4_inst%o2_decomp_depth_unsat_col      , &
         conc_o2_unsat                 => ch4_inst%conc_o2_unsat_col              , &
         o2_decomp_depth_sat           => ch4_inst%o2_decomp_depth_sat_col        , &
         conc_o2_sat                   => ch4_inst%conc_o2_sat_col                , &
         finundated                    => ch4_inst%finundated_col                 , &
         smin_nh4_vr                   => soilbiogeochem_nitrogenstate_inst%smin_nh4_vr_col , &
         smin_no3_vr                   => soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col , &
         phr_vr                        => soilbiogeochem_carbonflux_inst%phr_vr_col, &
         w_scalar                      => soilbiogeochem_carbonflux_inst%w_scalar_col, &
         t_scalar                      => soilbiogeochem_carbonflux_inst%t_scalar_col, &
         denit_resp_coef               => params_inst%denitrif_respiration_coefficient, &
         denit_resp_exp                => params_inst%denitrif_respiration_exponent   , &
         denit_nitrate_coef            => params_inst%denitrif_nitrateconc_coefficient, &
         denit_nitrate_exp             => params_inst%denitrif_nitrateconc_exponent   , &
         k_nitr_max                    => params_inst%k_nitr_max                      , &
         gross_nmin_vr                 => soilbiogeochem_nitrogenflux_inst%gross_nmin_vr_col, &
         net_nmin_vr                   => soilbiogeochem_nitrogenflux_inst%net_nmin_vr_col  , &
         r_psi                         => soilbiogeochem_nitrogenflux_inst%r_psi_col        , &
         anaerobic_frac                => soilbiogeochem_nitrogenflux_inst%anaerobic_frac_col, &
         smin_no3_massdens_vr          => soilbiogeochem_nitrogenflux_inst%smin_no3_massdens_vr_col, &
         k_nitr_t_vr                   => soilbiogeochem_nitrogenflux_inst%k_nitr_t_vr_col , &
         k_nitr_ph_vr                  => soilbiogeochem_nitrogenflux_inst%k_nitr_ph_vr_col, &
         k_nitr_h2o_vr                 => soilbiogeochem_nitrogenflux_inst%k_nitr_h2o_vr_col, &
         k_nitr_vr                     => soilbiogeochem_nitrogenflux_inst%k_nitr_vr_col   , &
         wfps_vr                       => soilbiogeochem_nitrogenflux_inst%wfps_vr_col     , &
         fmax_denit_carbonsubstrate_vr => soilbiogeochem_nitrogenflux_inst%fmax_denit_carbonsubstrate_vr_col, &
         fmax_denit_nitrate_vr         => soilbiogeochem_nitrogenflux_inst%fmax_denit_nitrate_vr_col, &
         f_denit_base_vr               => soilbiogeochem_nitrogenflux_inst%f_denit_base_vr_col, &
         diffus                        => soilbiogeochem_nitrogenflux_inst%diffus_col     , &
         ratio_k1                      => soilbiogeochem_nitrogenflux_inst%ratio_k1_col   , &
         ratio_no3_co2                 => soilbiogeochem_nitrogenflux_inst%ratio_no3_co2_col, &
         soil_co2_prod                 => soilbiogeochem_nitrogenflux_inst%soil_co2_prod_col, &
         fr_WFPS                       => soilbiogeochem_nitrogenflux_inst%fr_WFPS_col   , &
         fr_pH                         => soilbiogeochem_nitrogenflux_inst%fr_pH_col     , &
         soil_bulkdensity              => soilbiogeochem_nitrogenflux_inst%soil_bulkdensity_col, &
         pot_f_nit_vr                  => soilbiogeochem_nitrogenflux_inst%pot_f_nit_vr_col, &
         pot_f_denit_vr                => soilbiogeochem_nitrogenflux_inst%pot_f_denit_vr_col, &
         n2_n2o_ratio_denit_vr         => soilbiogeochem_nitrogenflux_inst%n2_n2o_ratio_denit_vr_col, &
         afps_vr                       => soilbiogeochem_nitrogenflux_inst%afps_vr_col   , &
         nox_n2o_ratio_vr              => soilbiogeochem_nitrogenflux_inst%nox_n2o_ratio_vr_col, &
         adjsoilph_vr                  => soilbiogeochem_nitrogenflux_inst%adjsoilph_vr_col )

      surface_tension_water = params_inst%surface_tension_water
      rij_kro_a             = params_inst%rij_kro_a
      rij_kro_alpha         = params_inst%rij_kro_alpha
      rij_kro_beta          = params_inst%rij_kro_beta
      rij_kro_gamma         = params_inst%rij_kro_gamma
      rij_kro_delta         = params_inst%rij_kro_delta
      organic_max           = CNParamsShareInst%organic_max

      co2diff_con(1) = 0.1325_r8
      co2diff_con(2) = 0.0009_r8

      pH(bounds%begc:bounds%endc) = 6.5_r8
      adjsoilph_vr(bounds%begc:bounds%endc,1:nlevdecomp) = 6.0_r8

      k1_nitr             = 0.2_r8
      maxrate_perday      = 0.15_r8
      maxrate             = maxrate_perday / secspday
      base_flux_perday    = 0.1_r8
      base_flux           = base_flux_perday / (secspday * 1.e4_r8)
      nreduce             = 0.6_r8
      min_ammonium_perday = 0.3_r8
      min_ammonium        = min_ammonium_perday / secspday

      !-----------------------------------------------------------------
      ! Reset layer-level NOx diagnostics before layer loop calculations.
      !-----------------------------------------------------------------
      do j = 1, nlevdecomp
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            soilbiogeochem_nitrogenflux_inst%f_nox_nit_vr_col(c,j)         = 0._r8
            soilbiogeochem_nitrogenflux_inst%f_nox_denit_vr_col(c,j)       = 0._r8
            soilbiogeochem_nitrogenflux_inst%f_nox_nit_atmos_vr_col(c,j)   = 0._r8
            soilbiogeochem_nitrogenflux_inst%f_nox_denit_atmos_vr_col(c,j) = 0._r8
            soilbiogeochem_nitrogenflux_inst%pot_f_nox_nit_vr_col(c,j)     = 0._r8
            soilbiogeochem_nitrogenflux_inst%pot_f_nox_denit_vr_col(c,j)   = 0._r8
         end do
      end do

      do j = 1, nlevdecomp
         do fc = 1, num_soilc
            c = filter_soilc(fc)

            ! Default pulse values for all layers and all columns.
            ! Only j == 1 updates the column pulse state. Deeper layers inherit
            ! no direct pulse multiplier in this implementation because the
            ! rain-pulse trigger is defined from the top soil layer.
            pulse_fac              = 1._r8
            l_dry_days             = 0._r8
            l_dry_hours            = 1._r8
            t_pulse_days           = 0._r8
            pulse_raw_peak         = 0._r8
            pulse_term             = 0._r8
            pulse_trigger          = .false.
            dry_days_before_update = 0._r8
            theta_now              = 0._r8
            delta_theta            = 0._r8

            !------------------------------------------------------------
            ! Val Martin / Hudman-style rain-pulse state machine.
            !
            ! This block is evaluated only in the top soil layer. It updates a
            ! column-level multiplier, pulse_fac_col(c), later used by the
            ! nitrification NOx flux.
            !
            ! Important implementation detail:
            ! The dry-period length is saved in pulse_dry_days(c) at the moment
            ! of rewetting. Without this, antecedent_dry_days(c) may be reset
            ! before the pulse expression is evaluated, which would suppress
            ! the pulse or make it inconsistent.
            !------------------------------------------------------------
            if (j == 1) then

               ! Current top-layer volumetric water content.
               ! Bound it between 0 and the local saturated water content.
               theta_now = max(0._r8, min(watsat(c,1), h2osoi_vol(c,1)))

               ! Change in top-layer volumetric water content since the previous timestep.
               ! Positive values indicate rewetting.
               delta_theta = theta_now - prev_theta_top(c)

               ! Save the dry-period length before updating it at this timestep.
               ! This is the dry-period memory used to decide whether a rewetting
               ! event is eligible to start a pulse.
               dry_days_before_update = antecedent_dry_days(c)

               ! A new pulse can only be started if no previous pulse is still decaying.
               if (.not. pulse_active(c)) then

                  if (theta_now < pulse_dry_theta_crit) then

                     ! Soil remains dry. Accumulate antecedent dry period.
                     antecedent_dry_days(c) = antecedent_dry_days(c) + dt / secspday

                  else

                     ! Soil is not dry at this timestep. Check whether this wetting
                     ! followed a sufficiently long dry spell and whether the moisture
                     ! increase is large enough to count as rewetting.
                     pulse_trigger = (dry_days_before_update >= pulse_min_dry_days) .and. &
                                     (delta_theta > pulse_rewet_dtheta_crit)

                     if (pulse_trigger) then

                        ! Start a new pulse.
                        pulse_active(c) = .true.

                        ! Reset the pulse decay clock to the rewetting timestep.
                        pulse_elapsed_days(c) = 0._r8

                        ! Store the dry-period length that caused the pulse.
                        pulse_dry_days(c) = dry_days_before_update

                        ! Reset antecedent dry counter after the pulse has been triggered.
                        antecedent_dry_days(c) = 0._r8

                     else

                        ! Soil is wet, but there was no valid pulse trigger.
                        ! Reset dry memory because the dry spell has ended.
                        antecedent_dry_days(c) = 0._r8
                        pulse_dry_days(c)      = 0._r8

                     end if

                  end if

               end if

               if (pulse_active(c)) then

                  ! Dry-period length attached to the active pulse, in days.
                  l_dry_days = max(pulse_dry_days(c), 0._r8)

                  ! Convert dry-period length to hours because the published
                  ! logarithmic expression uses l_dry in hours.
                  l_dry_hours = max(l_dry_days * 24._r8, 1.0_r8)

                  ! Time since pulse activation, in days.
                  t_pulse_days = pulse_elapsed_days(c)

                  ! Peak pulse amplitude before exponential decay.
                  ! The max() prevents negative enhancement for weak dry spells.
                  pulse_raw_peak = max(0._r8, pulse_a * log(l_dry_hours) - pulse_b)

                  ! Decayed pulse term. 1.632 d-1 is equivalent to 0.068 h-1.
                  pulse_term = pulse_raw_peak * exp(-pulse_decay_per_day * t_pulse_days)

                  ! Final multiplier applied to nitrification-derived NOx.
                  ! pulse_fac = 1 means no enhancement.
                  pulse_fac = 1._r8 + pulse_term

                  ! Advance pulse age by one model timestep.
                  pulse_elapsed_days(c) = pulse_elapsed_days(c) + dt / secspday

                  ! Stop tracking the pulse once its contribution is negligible.
                  ! This avoids carrying a permanently active but numerically dead pulse.
                  if (pulse_term <= pulse_term_min) then
                     pulse_active(c)       = .false.
                     pulse_elapsed_days(c) = 0._r8
                     pulse_dry_days(c)     = 0._r8
                  end if

               end if

               ! Store current moisture for rewetting calculation at the next timestep.
               prev_theta_top(c) = theta_now

               ! Write pulse multiplier to column-level history diagnostic.
               soilbiogeochem_nitrogenflux_inst%pulse_fac_col(c) = pulse_fac

            end if

            !---------------- soil anoxia state / diffusivity
            f_a = 1._r8 - watfc(c,j) / watsat(c,j)
            eps = watsat(c,j) - watfc(c,j)

            if (use_lch4) then

               if (organic_max > 0._r8) then
                  om_frac = min(cellorg(c,j) / organic_max, 1._r8)
               else
                  om_frac = 1._r8
               end if

               diffus(c,j) = (d_con_g(2,1) + d_con_g(2,2) * t_soisno(c,j)) * 1.e-4_r8 * &
                             (om_frac * f_a**(10._r8/3._r8) / watsat(c,j)**2 + &
                             (1._r8 - om_frac) * eps**2 * f_a**(3._r8 / bsw(c,j)))

               r_min(c,j) = 2._r8 * surface_tension_water / (rho_w * grav * abs(soilpsi(c,j)))
               r_max      = 2._r8 * surface_tension_water / (rho_w * grav * 0.1_r8)
               r_psi(c,j) = sqrt(r_min(c,j) * r_max)

               ratio_diffusivity_water_gas(c,j) = (d_con_g(2,1) + d_con_g(2,2)*t_soisno(c,j)) * 1.e-4_r8 / &
                    ((d_con_w(2,1) + d_con_w(2,2)*t_soisno(c,j) + d_con_w(2,3)*t_soisno(c,j)**2) * 1.e-9_r8)

               if (o2_decomp_depth_unsat(c,j) > 0._r8) then
                  anaerobic_frac(c,j) = exp(-rij_kro_a * r_psi(c,j)**(-rij_kro_alpha) * &
                       o2_decomp_depth_unsat(c,j)**(-rij_kro_beta) * &
                       conc_o2_unsat(c,j)**rij_kro_gamma * (h2osoi_vol(c,j) + ratio_diffusivity_water_gas(c,j) * &
                       watsat(c,j))**rij_kro_delta)
               else
                  anaerobic_frac(c,j) = 0._r8
               end if

               if (anoxia_wtsat) then
                  r_min_sat = 2._r8 * surface_tension_water / (rho_w * grav * abs(grav * 1.e-6_r8 * sucsat(c,j)))
                  r_psi_sat = sqrt(r_min_sat * r_max)
                  if (o2_decomp_depth_sat(c,j) > 0._r8) then
                     anaerobic_frac_sat = exp(-rij_kro_a * r_psi_sat**(-rij_kro_alpha) * &
                          o2_decomp_depth_sat(c,j)**(-rij_kro_beta) * &
                          conc_o2_sat(c,j)**rij_kro_gamma * (watsat(c,j) + ratio_diffusivity_water_gas(c,j) * &
                          watsat(c,j))**rij_kro_delta)
                  else
                     anaerobic_frac_sat = 0._r8
                  end if
                  anaerobic_frac(c,j) = (1._r8 - finundated(c)) * anaerobic_frac(c,j) + &
                                         finundated(c) * anaerobic_frac_sat
               end if

            else
               anaerobic_frac(c,j) = 0._r8
               diffus(c,j)         = 0._r8
            end if

            !---------------- nitrification
            k_nitr_t_vr(c,j)  = min(t_scalar(c,j), 1._r8)
            k_nitr_ph_vr(c,j) = 0.56_r8 + atan(rpi * 0.45_r8 * (-5._r8 + pH(c))) / rpi
            k_nitr_h2o_vr(c,j)= w_scalar(c,j)

            k_nitr_vr(c,j) = k_nitr_max * k_nitr_t_vr(c,j) * k_nitr_h2o_vr(c,j) * k_nitr_ph_vr(c,j)

            pot_f_nit_vr(c,j) = max(k1_nitr * gross_nmin_vr(c,j) + smin_nh4_vr(c,j) * k_nitr_vr(c,j), 0._r8)
            pot_f_nit_vr(c,j) = pot_f_nit_vr(c,j) * (1._r8 - anaerobic_frac(c,j))
            pot_f_nit_vm      = pot_f_nit_vr(c,j)

            if (t_soisno(c,j) <= SHR_CONST_TKFRZ .and. no_frozen_nitrif_denitrif) then
               pot_f_nit_vr(c,j) = 0._r8
               pot_f_nit_vm      = 0._r8
            end if

            !---------------- denitrification
            soil_hr_vr(c,j) = phr_vr(c,j)

            soil_bulkdensity(c,j) = bd(c,j) + h2osoi_liq(c,j) / col%dz(c,j)

            g_per_m3__to__ug_per_gsoil         = 1.e3_r8 / soil_bulkdensity(c,j)
            g_per_m3_sec__to__ug_per_gsoil_day = g_per_m3__to__ug_per_gsoil * secspday

            smin_no3_massdens_vr(c,j) = max(smin_no3_vr(c,j), 0._r8) * g_per_m3__to__ug_per_gsoil
            soil_co2_prod(c,j)        = soil_hr_vr(c,j) * g_per_m3_sec__to__ug_per_gsoil_day

            fmax_denit_carbonsubstrate_vr(c,j) = (denit_resp_coef * soil_co2_prod(c,j)**denit_resp_exp) / &
                                                 g_per_m3_sec__to__ug_per_gsoil_day

            fmax_denit_nitrate_vr(c,j) = (denit_nitrate_coef * smin_no3_massdens_vr(c,j)**denit_nitrate_exp) / &
                                         g_per_m3_sec__to__ug_per_gsoil_day

            f_denit_base_vr(c,j) = max(min(fmax_denit_carbonsubstrate_vr(c,j), fmax_denit_nitrate_vr(c,j)), 0._r8)

            if (t_soisno(c,j) <= SHR_CONST_TKFRZ .and. no_frozen_nitrif_denitrif) then
               f_denit_base_vr(c,j) = 0._r8
            end if

            fr_ph_denit    = 0.0016_r8 * exp(1.006_r8 * pH(c))
            pot_f_denit_vr(c,j) = f_denit_base_vr(c,j) * anaerobic_frac(c,j)
            pot_f_denit_vm = pot_f_denit_vr(c,j) * fr_ph_denit

            !---------------- N2O/N2 ratio
            ratio_k1(c,j) = max(1.7_r8, 38.4_r8 - 350._r8 * diffus(c,j))

            if (soil_co2_prod(c,j) > 1.0e-9_r8) then
               ratio_no3_co2(c,j) = smin_no3_massdens_vr(c,j) / soil_co2_prod(c,j)
            else
               ratio_no3_co2(c,j) = 100._r8
            end if

            wfps_vr(c,j) = max(min(h2osoi_vol(c,j)/watsat(c,j), 1._r8), 0._r8) * 100._r8
            fr_WFPS(c,j) = max(0.1_r8, 0.015_r8 * wfps_vr(c,j) - 0.32_r8)

            if (use_lch4) then
               if (anoxia_wtsat) then
                  fr_WFPS(c,j) = fr_WFPS(c,j) * (1._r8 - finundated(c)) + finundated(c) * 1.18_r8
               end if
            end if

            n2_n2o_ratio_denit_vr(c,j) = max(0.16_r8*ratio_k1(c,j), ratio_k1(c,j)*exp(-0.8_r8*ratio_no3_co2(c,j))) * fr_WFPS(c,j)

            if (pH(c) <= 4._r8) then
               fr_pH(c,j) = 0.001_r8
            elseif (pH(c) < 7._r8) then
               fr_pH(c,j) = 0.001_r8 + (pH(c) - 4._r8) / 3._r8
            else
               fr_pH(c,j) = 1._r8
            end if

            n2_n2o_ratio_denit_vm = max(1.0e-12_r8, n2_n2o_ratio_denit_vr(c,j) * fr_pH(c,j))

            !---------------- NOx/N2O ratio
            afps_vr(c,j) = 1._r8 - max(min(h2osoi_vol(c,j)/watsat(c,j), 1._r8), 0._r8)
            Dr           = 0.209_r8 * afps_vr(c,j)**(4._r8/3._r8)

            nox_n2o_ratio_vr(c,j) = 15.2_r8 + (35.5_r8 * atan(0.75_r8 * rpi * (10.0_r8 * Dr - 1.86_r8))) / rpi
            nox_n2o_ratio_vr(c,j) = max(0._r8, nox_n2o_ratio_vr(c,j))

            ! ---- Layer N2O fluxes
            f_n2o_nit_frac_vm  = 721.86_r8 * exp(-2.387_r8 * pH(c))
            f_n2o_nit_frac_vm  = min(1._r8, max(0._r8, f_n2o_nit_frac_vm))
            f_n2o_nit_layer_vm = pot_f_nit_vm * f_n2o_nit_frac_vm

            if (pot_f_denit_vm > 0._r8 .and. n2_n2o_ratio_denit_vm > 0._r8) then
               f_n2o_denit_layer_vm = pot_f_denit_vm / (1._r8 + n2_n2o_ratio_denit_vm)
            else
               f_n2o_denit_layer_vm = 0._r8
            end if

            ! ---- NOx construction intentionally disabled in this module.
            !
            ! Val Martin constructs final NOx in SoilBiogeochemCompetitionMod.F90
            ! after realised f_nit_vr and f_denit_vr are known. To avoid two
            ! competing NOx pathways, this module now keeps only the supporting
            ! diagnostics needed downstream, especially nox_n2o_ratio_vr,
            ! n2_n2o_ratio_denit_vr, f_n2o_nit_layer_vm and f_n2o_denit_layer_vm.
            !
            ! The old potential-style NOx code is retained below as comments for
            ! traceability, but it must not write final NOx history fields here.
            f_nox_nit_layer   = 0._r8
            f_nox_denit_layer = 0._r8

            ! Old potential-style NOx calculation, no longer active:
            ! f_nox_nit_layer   = nox_n2o_ratio_vr(c,j) * f_n2o_nit_layer_vm
            ! f_nox_denit_layer = nox_n2o_ratio_vr(c,j) * f_n2o_denit_layer_vm

            ! Old debug print for disabled NOx pathway, no longer active:
            ! if (c == bounds%begc .and. j == 1) then
            !    write(*,*) 'PULSEDBG', 'c=', c, 'j=', j, &
            !               'theta_now=', theta_now, &
            !               'delta_theta=', delta_theta, &
            !               'dry_days_mem=', antecedent_dry_days(c), &
            !               'pulse_dry_days=', pulse_dry_days(c), &
            !               'pulse_elapsed_days=', pulse_elapsed_days(c), &
            !               'pulse_active=', pulse_active(c), &
            !               'pulse_fac=', pulse_fac, &
            !               'f_nox_nit_layer=', f_nox_nit_layer, &
            !               'f_nox_denit_layer=', f_nox_denit_layer
            !
            !    write(*,*) 'PULSEDBG_FORM', &
            !               'a=', pulse_a, &
            !               'b=', pulse_b, &
            !               'dry_days=', l_dry_days, &
            !               'dry_hours=', l_dry_hours, &
            !               'raw_peak=', pulse_raw_peak, &
            !               'pulse_term=', pulse_term, &
            !               'pulse_fac=', pulse_fac
            ! end if

            ! Do not write NOx fluxes from this module. They are reset to zero
            ! above and are overwritten by SoilBiogeochemCompetitionMod.F90.
            ! Old writes retained as comments:
            ! soilbiogeochem_nitrogenflux_inst%f_nox_nit_vr_col(c,j)       = f_nox_nit_layer
            ! soilbiogeochem_nitrogenflux_inst%f_nox_denit_vr_col(c,j)     = f_nox_denit_layer
            ! soilbiogeochem_nitrogenflux_inst%pot_f_nox_nit_vr_col(c,j)   = nox_n2o_ratio_vr(c,j) * pot_f_nit_vr(c,j)
            ! soilbiogeochem_nitrogenflux_inst%pot_f_nox_denit_vr_col(c,j) = nox_n2o_ratio_vr(c,j) * pot_f_denit_vr(c,j)

            fcan = 1._r8
            if (present(canopy_red_fac_col)) fcan = max(0._r8, min(1._r8, canopy_red_fac_col(c)))

            ! Old atmospheric NOx writes retained as comments:
            ! soilbiogeochem_nitrogenflux_inst%f_nox_nit_atmos_vr_col(c,j)   = fcan * f_nox_nit_layer
            ! soilbiogeochem_nitrogenflux_inst%f_nox_denit_atmos_vr_col(c,j) = fcan * f_nox_denit_layer

            soilbiogeochem_nitrogenflux_inst%f_n2o_nit_vm_col(c)      = soilbiogeochem_nitrogenflux_inst%f_n2o_nit_vm_col(c)      + f_n2o_nit_layer_vm   * col%dz(c,j)
            soilbiogeochem_nitrogenflux_inst%f_n2o_denit_vm_col(c)    = soilbiogeochem_nitrogenflux_inst%f_n2o_denit_vm_col(c)    + f_n2o_denit_layer_vm * col%dz(c,j)

            ! Old column-integrated NOx writes retained as comments:
            ! soilbiogeochem_nitrogenflux_inst%f_nox_nit_col(c)         = soilbiogeochem_nitrogenflux_inst%f_nox_nit_col(c)         + f_nox_nit_layer      * col%dz(c,j)
            ! soilbiogeochem_nitrogenflux_inst%f_nox_denit_col(c)       = soilbiogeochem_nitrogenflux_inst%f_nox_denit_col(c)       + f_nox_denit_layer    * col%dz(c,j)
            ! soilbiogeochem_nitrogenflux_inst%f_nox_nit_atmos_col(c)   = soilbiogeochem_nitrogenflux_inst%f_nox_nit_atmos_col(c)   + fcan * f_nox_nit_layer   * col%dz(c,j)
            ! soilbiogeochem_nitrogenflux_inst%f_nox_denit_atmos_col(c) = soilbiogeochem_nitrogenflux_inst%f_nox_denit_atmos_col(c) + fcan * f_nox_denit_layer * col%dz(c,j)
            ! soilbiogeochem_nitrogenflux_inst%pot_f_nox_nit_col(c)     = soilbiogeochem_nitrogenflux_inst%pot_f_nox_nit_col(c)     + &
            !                                                            nox_n2o_ratio_vr(c,j) * pot_f_nit_vr(c,j)   * col%dz(c,j)
            ! soilbiogeochem_nitrogenflux_inst%pot_f_nox_denit_col(c)   = soilbiogeochem_nitrogenflux_inst%pot_f_nox_denit_col(c)   + &
            !                                                            nox_n2o_ratio_vr(c,j) * pot_f_denit_vr(c,j) * col%dz(c,j)

         end do
      end do

    end associate

  end subroutine SoilBiogeochemNitrifDenitrif

end module SoilBiogeochemNitrifDenitrifMod
