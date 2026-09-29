module metrics
  use nrtype
  implicit none
  private
  public :: get_KGE, get_KGEp, get_composite_KGE, get_KGE_CSI, get_NSE, get_MAE, get_RMSE

contains

! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Cyril Thébault, 2024
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Allow different transformation on streamflow before using metrics
! ---------------------------------------------------------------------------------------
  
  ! Compute transformation on streamflow
  subroutine apply_transformation(obs, sim, transfo)
  
    real(wp), intent(inout)      :: obs(:), sim(:)     ! observed and simulated streamflow time series
    character(len=*), intent(in) :: transfo            ! transformation applied to streamflow time series
    real(wp)                     :: eps, transfo_val   ! epsilon value to avoid infinite number, power transformation as float
    integer                      :: i                  ! loop index
  
    ! log transformation
    if (transfo == 'log') then
  
      eps = sum(obs) / (100 * size(obs)) ! use to avoid log(0)
    
      do i = 1, size(obs)
          obs(i) = log(eps + obs(i))
          sim(i) = log(eps + sim(i))
      end do
    
    ! box-cox transformation  
    else if (transfo == 'boxcox') then
      do i = 1, size(obs)
        obs(i) = (obs(i)**0.25 - 1) / 0.25
        sim(i) = (sim(i)**0.25 - 1) / 0.25
      end do
      
    ! power transformation
    else
      transfo_val = char_to_float(transfo)
    
      if (transfo_val < 0) then
    
        eps = sum(obs) / (100 * size(obs)) ! use to avoid division by 0
      
        do i = 1, size(obs)
          obs(i) = (eps + obs(i)) ** transfo_val
          sim(i) = (eps+ sim(i)) ** transfo_val
        end do
      
      else
    
        do i = 1, size(obs)
          obs(i) = obs(i) ** transfo_val
          sim(i) = sim(i) ** transfo_val
        end do
      
      end if
    
    end if
  
  end subroutine apply_transformation

! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Cyril Thébault, 2024
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Allow different metrics as objective function
! ---------------------------------------------------------------------------------------

  ! Compute KGE
  function get_KGE(obs, sim, transfo) result(kge)
  
    real(wp), intent(in)                   :: obs(:), sim(:)                                    ! observed and simulated streamflow time series
    character(len=*), intent(in), optional :: transfo                                           ! transformation applied to streamflow time series
    real(wp)                               :: kge                                               ! KGE
    real(wp), allocatable                  :: obs_t(:), sim_t(:)                                      ! transformed time series without NA values
    real(wp)                               :: sd_sim, sd_obs, m_sim, m_obs, r, alpha, beta      ! KGE variables
    integer                                :: n  
    
    n = count(.not. (isnan(obs) .or. isnan(sim)))
    
    allocate(obs_t(n), sim_t(n))
    
    obs_t = pack(obs, .not. (isnan(obs) .or. isnan(sim)))
    sim_t = pack(sim, .not. (isnan(obs) .or. isnan(sim)))
    
    call apply_transformation(obs_t, sim_t, transfo)
    
    
    ! compute mean
    m_obs  = sum(obs_t) / n
    m_sim  = sum(sim_t) / n
      
    ! compute standard deviation
    sd_obs  = standard_deviation(obs_t)
    sd_sim  = standard_deviation(sim_t)
              
    ! compute correlation
    r = correlation(sim_t, obs_t)  ! correlation
    alpha = sd_sim / sd_obs        ! variation
    beta = m_sim / m_obs           ! bias
    
    kge = 1.0 - sqrt((r-1)**2 + (alpha-1)**2 + (beta-1)**2)
    
    if (isnan(kge)) then
        kge = -1.0e6                ! bad value
    end if
    
  end function get_KGE


  ! Compute KGEp
  function get_KGEp(obs, sim, transfo) result(kgep)
  
    real(wp), intent(in)                   :: obs(:), sim(:)                                    ! observed and simulated streamflow time series
    character(len=*), intent(in), optional :: transfo                                           ! transformation applied to streamflow time series
    real(wp)                               :: kgep                                              ! Kling-Gupta Efficiency' (Kling et al. in 2012)
    real(wp), allocatable                  :: obs_t(:), sim_t(:)                                      ! transformed time series without NA values
    real(wp)                               :: sd_sim, sd_obs, m_sim, m_obs, r, alphap, beta     ! KGE' variables
    integer                                :: n                                                 ! length of the time series without NA values
    
    n = count(.not. (isnan(obs) .or. isnan(sim)))
    
    allocate(obs_t(n), sim_t(n))
    
    obs_t = pack(obs, .not. (isnan(obs) .or. isnan(sim)))
    sim_t = pack(sim, .not. (isnan(obs) .or. isnan(sim)))
    
    call apply_transformation(obs_t, sim_t, transfo)
    
    
    ! compute mean
    m_obs  = sum(obs_t) / n
    m_sim  = sum(sim_t) / n
      
    ! compute standard deviation
    sd_obs  = standard_deviation(obs_t)
    sd_sim  = standard_deviation(sim_t)
              
    ! KGEp component
    r = correlation(sim_t, obs_t)             ! correlation
    alphap = (sd_sim/m_sim) / (sd_obs/m_obs)  ! relative variation
    beta = m_sim / m_obs                      ! bias
    
    kgep = 1.0 - sqrt((r-1)**2 + (alphap-1)**2 + (beta-1)**2)
    
    if (isnan(kgep)) then
        kgep = -1.0e6                          ! bad value
    end if
        
  end function get_KGEp


  ! Compute NSE
  function get_NSE(obs, sim, transfo) result(nse)
  
    real(wp), intent(in)                   :: obs(:), sim(:)      ! observed and simulated streamflow time series
    character(len=*), intent(in), optional :: transfo             ! transformation applied to streamflow time series
    real(wp)                               :: nse                 ! Nash-Sutcliff score
    real(wp), allocatable                  :: obs_t(:), sim_t(:)  ! transformed time series without NA values
    integer                                :: n                   ! length of the time series without NA values
    
    n = count(.not. (isnan(obs) .or. isnan(sim)))
    
    allocate(obs_t(n), sim_t(n))
    
    obs_t = pack(obs, .not. (isnan(obs) .or. isnan(sim)))
    sim_t = pack(sim, .not. (isnan(obs) .or. isnan(sim)))
    
    call apply_transformation(obs_t, sim_t, transfo)
    
    nse = 1 - (sum((obs_t - sim_t)**2) / sum((obs_t - (sum(obs_t)/n))**2))
    
    if (isnan(nse)) then
        nse = -1.0e6                          ! bad value
    end if
    
  end function get_NSE


  ! Compute MAE
  function get_MAE(obs, sim, transfo) result(mae)
  
    real(wp), intent(in)                   :: obs(:), sim(:)      ! observed and simulated streamflow time series
    character(len=*), intent(in), optional :: transfo             ! transformation applied to streamflow time series
    real(wp)                               :: mae                 ! mean absolute error
    real(wp), allocatable                  :: obs_t(:), sim_t(:)  ! transformed time series without NA values
    integer                                :: n                   ! length of the time series without NA values

        
    n = count(.not. (isnan(obs) .or. isnan(sim)))
    
    allocate(obs_t(n), sim_t(n))
    
    obs_t = pack(obs, .not. (isnan(obs) .or. isnan(sim)))
    sim_t = pack(sim, .not. (isnan(obs) .or. isnan(sim)))
    
    call apply_transformation(obs_t, sim_t, transfo)
    
    mae = sum(abs(obs_t - sim_t)) / n
    
    if (isnan(mae)) then
        mae = 1.0e6                          ! bad value
    end if
    
  end function get_MAE


  ! Compute RMSE
  function get_RMSE(obs, sim, transfo) result(rmse)
  
    real(wp), intent(in)                   :: obs(:), sim(:)      ! observed and simulated streamflow time series
    character(len=*), intent(in), optional :: transfo             ! transformation applied to streamflow time series
    real(wp)                               :: rmse                ! root mean square error
    real(wp), allocatable                  :: obs_t(:), sim_t(:)  ! transformed time series without NA values
    integer                                :: n                   ! length of the time series without NA values
    
    n = count(.not. (isnan(obs) .or. isnan(sim)))
    
    allocate(obs_t(n), sim_t(n))
    
    obs_t = pack(obs, .not. (isnan(obs) .or. isnan(sim)))
    sim_t = pack(sim, .not. (isnan(obs) .or. isnan(sim)))
    
    call apply_transformation(obs_t, sim_t, transfo)
    
    rmse = sqrt(sum((obs_t - sim_t)**2) / n)
    
    if (isnan(rmse)) then
        rmse = 1.0e6                          ! bad value
    end if
    
  end function get_RMSE


  ! Compute composite KGE: mean of KGE on raw flows and KGE on inverse flows
  ! (ported from the FUSE_KGECOMP branch of CyrilThebault/fuse, commit e6e23a4)
  function get_composite_KGE(obs, sim) result(composite_kge)

    real(wp), intent(in) :: obs(:), sim(:)           ! observed and simulated streamflow time series
    real(wp)             :: composite_kge            ! composite KGE
    real(wp)             :: kge1, kge2               ! KGE values for different transformations

    ! Compute KGE for transfo = 1 (raw flows, emphasis on high flows)
    kge1 = get_KGE(obs, sim, "1")

    ! Compute KGE for transfo = -1 (inverse flows, emphasis on low flows)
    kge2 = get_KGE(obs, sim, "-1")

    ! Compute composite KGE
    composite_kge = (kge1 + kge2) / 2.0_wp

  end function get_composite_KGE


  ! Compute KGE-CSI: f0*CSI + (1-f0)*KGE_nonzero
  !   tau         = zero-flow threshold, mean of the positive observed flows / 1000
  !   f0          = fraction of observed time steps with obs <= tau
  !   CSI         = critical success index for zero flow, H / (H + M + FA), 0 if H + M + FA = 0
  !                 (H: obs <= tau and sim <= tau; M: obs > tau and sim <= tau; FA: obs <= tau and sim > tau)
  !   KGE_nonzero = KGE on the time steps with obs > tau, after the transformation TRANSFO
  ! tau, f0 and CSI always use the untransformed flows. If KGE_nonzero cannot be computed
  ! (fewer than 10 non-zero observations, constant observed or simulated flow), it is
  ! replaced by the penalty KGECSI_PENALTY so that SCE moves away from that parameter set.
  ! Follows the Python implementation used for the CAMELS calibrations (external_sceua_camels_ind.py),
  ! except that f0 is computed over the time steps with observations only.
  subroutine get_KGE_CSI(obs, sim, transfo, kgecsi, csi, kge_nz, f0, tau)

    real(wp), intent(in)          :: obs(:), sim(:)   ! observed and simulated streamflow (time steps with observations only)
    character(len=*), intent(in)  :: transfo          ! transformation applied before computing KGE_nonzero
    real(wp), intent(out)         :: kgecsi           ! KGE-CSI
    real(wp), intent(out)         :: csi              ! critical success index for zero flow
    real(wp), intent(out)         :: kge_nz           ! KGE on non-zero observed flows (-9999 if it cannot be computed)
    real(wp), intent(out)         :: f0               ! fraction of zero-flow observations
    real(wp), intent(out)         :: tau              ! zero-flow threshold
    real(wp), parameter           :: KGECSI_PENALTY = -9999._wp  ! replaces KGE_nonzero when it cannot be computed
    integer, parameter            :: MIN_NONZERO = 10            ! minimum number of non-zero observations for KGE_nonzero
    logical, allocatable          :: obs_zero(:), sim_zero(:)    ! zero-flow flags
    integer                       :: n, n_hit, n_miss, n_false, n_nz

    n = size(obs)

    ! zero-flow threshold and fraction of zero-flow observations
    if (any(obs > 0._wp)) then
      tau = max(sum(obs, mask=obs > 0._wp) / count(obs > 0._wp) / 1000._wp, 1.e-12_wp)
    else
      tau = 1.e-9_wp
    end if
    obs_zero = obs <= tau
    sim_zero = sim <= tau
    f0 = real(count(obs_zero), wp) / real(n, wp)

    ! critical success index for zero flow
    n_hit   = count(      obs_zero  .and.       sim_zero )
    n_miss  = count((.not. obs_zero) .and.       sim_zero )
    n_false = count(      obs_zero  .and. (.not. sim_zero))
    if (n_hit + n_miss + n_false > 0) then
      csi = real(n_hit, wp) / real(n_hit + n_miss + n_false, wp)
    else
      csi = 0._wp
    end if

    ! KGE on the non-zero observed flows
    n_nz = count(.not. obs_zero)
    kge_nz = KGECSI_PENALTY
    if (n_nz >= MIN_NONZERO) then
      if (standard_deviation(pack(obs, .not. obs_zero)) > 1.e-12_wp .and. &
          standard_deviation(pack(sim, .not. obs_zero)) > 0._wp) then
        kge_nz = get_KGE(pack(obs, .not. obs_zero), pack(sim, .not. obs_zero), transfo)
        if (isnan(kge_nz) .or. kge_nz <= -1.0e6_wp) kge_nz = KGECSI_PENALTY  ! get_KGE returns -1.0e6 for NaN
      end if
    end if

    kgecsi = f0 * csi + (1._wp - f0) * kge_nz

  end subroutine get_KGE_CSI


! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Cyril Thébault, 2024
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Additional functions
! ---------------------------------------------------------------------------------------

  ! Function to convert character to float
  function char_to_float(char_val) result(float_val)
  
    character(len=*), intent(in) :: char_val
    real(wp)                     :: float_val
    integer                      :: io_stat

    read(char_val, *, iostat=io_stat) float_val
    if (io_stat /= 0) then
      print *, "Error: Unable to convert '", char_val, "' to float"
      float_val = 1.0  ! set to 1 for no transformation
    end if
  end function char_to_float


  ! Function to compute standard deviation
  function standard_deviation(x) result(sd)
    real(wp), intent(in) :: x(:)
    real(wp)             :: sd
    real(wp)             :: mean_x
    integer              :: n
    
    n = size(x)
    mean_x = sum(x) / n
    sd = sqrt(sum((x - mean_x)**2) / (n - 1))
  end function standard_deviation

  ! Function to compute correlation
  function correlation(x, y) result(r)
    real(wp), intent(in) :: x(:), y(:)
    real(wp)             :: r
    real(wp)             :: mean_x, mean_y, sd_x, sd_y
    integer              :: n, i
    
    n = size(x)
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    sd_x = standard_deviation(x)
    sd_y = standard_deviation(y)
    
    r = sum((x - mean_x) * (y - mean_y)) / (n * sd_x * sd_y)
  end function correlation

  elemental function isnan(x) result(is_nan)
    real(wp), intent(in) :: x
    logical              :: is_nan
  
    is_nan = (x /= x)
  end function isnan

end module metrics