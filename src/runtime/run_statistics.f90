module run_statistics

  use nrtype,      only: i4b, wp, lgt
  use info_types,  only: fuse_info
  use work_types,  only: fuse_work
  use domain_types,only: domain_data

  use fuse_globaldata, only: isPrint
  use fuse_globaldata, only: do_mizuRoute

  use globalData, only: length_conv,time_conv ! mizuRoute shim

  implicit none
  private

  public :: accumulate_numerical_stats
  public :: finalize_numerical_stats
  public :: compute_performance_stats

contains

  SUBROUTINE compute_performance_stats(info, work, domain, ierr,message)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Nans Addor to deal with NA values in QOBS, 2016
  ! Modified by Cyril Thébault to allow different metrics as objective function, 2024
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes summary statistics from model simulations
  ! ---------------------------------------------------------------------------------------
  ! Modules Modified:
  ! -----------------
  ! MODULE multistats -- summary statistics stored in MODULE multistats
  ! ---------------------------------------------------------------------------------------
  USE nrtype                                            ! variable types, etc.
  USE fuse_fileManager,only: METRIC, TRANSFO            ! metric and transformation requested in the filemanager
  ! FUSE modules
  USE metrics                                           ! available metrics and transformations
  USE multiforce, only: NA_VALUE                        ! model forcing structure (temporally constant)
  USE multiforce, only: sim_beg, eval_beg, eval_end     ! model forcing structure (temporally constant)

  IMPLICIT NONE

  ! input/output
  type(fuse_info),   intent(in)          :: info        ! domain info
  type(fuse_work),   intent(inout)       :: work        ! work structures that depend on npar/nState
  type(domain_data), intent(inout)       :: domain      ! data for the full domain

  integer(i4b)     , intent(out)         :: ierr
  character(*)     , intent(out)         :: message

  ! internal
  INTEGER(I4B)                           :: I           ! looping
  INTEGER(I4B)                           :: NS          ! number of samples
  REAL(WP), DIMENSION(:), ALLOCATABLE    :: QOBS        ! observed runoff - whole time series
  REAL(WP), DIMENSION(:), ALLOCATABLE    :: QOBS_AVAIL  ! observed runoff - only time steps with QOBS available
  INTEGER(I4B)                           :: NUM_AVAIL   ! number of time steps with QOBS available
  LOGICAL(LGT),DIMENSION(:), ALLOCATABLE :: QOBS_MASK   ! boolean: is QOBS available?
  REAL(WP), DIMENSION(:), ALLOCATABLE    :: QSIM        ! simulated runoff - whole time series
  REAL(WP), DIMENSION(:), ALLOCATABLE    :: QSIM_AVAIL  ! simulated runoff - only time steps with QOBS available
  REAL(WP), DIMENSION(:), ALLOCATABLE    :: DOBS        ! observed runoff anomalies
  REAL(WP), DIMENSION(:), ALLOCATABLE    :: DSIM        ! simulated runoff anomalies
  REAL(WP), DIMENSION(:), ALLOCATABLE    :: RAWD        ! observed-simulated differences in flow
  REAL(WP), DIMENSION(:), ALLOCATABLE    :: LOGD        ! observed-simulated differences in LOG flow
  REAL(WP)                               :: XB_OBS      ! mean observed runoff
  REAL(WP)                               :: XB_SIM      ! mean simulated runoff
  REAL(WP)                               :: SS_OBS      ! sum of squared observed runoff anomalies
  REAL(WP)                               :: SS_SIM      ! sum of squared simulated runoff anomalies
  REAL(WP)                               :: SS_LOBS     ! sum of squared lagged differences in observed runoff
  REAL(WP)                               :: SS_LSIM     ! sum of squared lagged differences in simulated runoff
  REAL(WP)                               :: SS_RAW      ! sum of squared differences in observed - simulated
  REAL(WP)                               :: SS_LOG      ! sum of squared differences in LOG observed - LOG simulated
  REAL(WP)                               :: NO_ZERO     ! avoid divide by zero
  integer(i4b)                           :: iSeg        ! index of stream segment
  integer(i4b)                           :: iSpat1      ! index of the 1st spatial dimension
  integer(i4b)                           :: iSpat2      ! index of the 2nd spatial dimension
  integer(i4b)                           :: ixStart     ! start index for statistics
  integer(i4b)                           :: ixEnd       ! end index for statistics
  real(wp)                               :: upsarea     ! upstream area (m2)
  character(len=256)                     :: cmessage    ! error message of downwind routine

  ierr = 0
  message = 'compute_performance_stats/' 

  ! ---------------------------------------------------------------------------------------
  ! (1) PRELIMINARIES
  ! ---------------------------------------------------------------------------------------
  ! define sample size
  NS =  eval_end-eval_beg+1
  
  ! allocate space for observed and simulated runoff
  ALLOCATE(QOBS(NS),QOBS_MASK(NS),QSIM(NS),STAT=IERR)
  if (ierr /= 0)then
    message=trim(message)//'PROBLEM ALLOCATING SPACE'
    return
  endif
 
  ! define start/end indices
  ! note that sim_beg, eval_beg, sim_end, eval_end are all with respect to julian_day_input
  ixStart = eval_beg-sim_beg+1
  ixEnd   = eval_end-sim_beg+1

  ! ---- extract SIM for evaluation period ------------------------------------------------

  if ( do_mizuRoute ) then

    upsarea    = domain%reach%totArea(info%ntopo%ixSegOut)

    QSIM(1:NS) = domain%river_network%method(1)%streamflow( &
                 info%ntopo%ixSegOut, ixStart:ixEnd)     &
                 / (time_conv * length_conv * upsarea)

  ! no routing -- compute average weighted by overlap area
  else

    QSIM(:) = 0._wp

    do iSpat2 = 1, size(domain%route,2)
      do iSpat1 = 1, size(domain%route,1)
        QSIM(1:NS) = QSIM(1:NS)  + domain%olap_area(iSpat1,iSpat2) * &
                                   domain%route(iSpat1,iSpat2,ixStart:ixEnd)%Q_ROUTED
      end do 
    end do
    
    QSIM(:) = QSIM(:) / sum(domain%olap_area)
 
  endif

  ! ---------------------------------------------------------------------------------------
  
  ! currently assuming one observation (checked earlier)
  QOBS(:) = work%obs%q(1, ixStart:ixEnd) 
  
  ! check for missing QOBS values
  QOBS_MASK = QOBS.ne.REAL(NA_VALUE, KIND(WP)) ! find the time steps for which QOBS is available
  NUM_AVAIL = COUNT(QOBS_MASK) ! number of time steps for which QOBS is available
  
  if(isPrint)then
    PRINT *, 'Number of time steps in evaluation period (EP) = ', NS
    PRINT *, 'Number of time steps with observed streamflow in EP = ', NUM_AVAIL
  endif
  
  IF (NUM_AVAIL.EQ.0) THEN
  
    PRINT *, 'Skiping computation of error statistics because no observed streamflow data'
    work%run%stats%NASH_SUTT=-9999
    work%run%stats%RAW_RMSE=-9999
    work%run%stats%KGE=-9999
    work%run%stats%KGECOMP=-9999
    work%run%stats%METRIC_VAL=-9999
  
  ELSE
  
    ! extract elements from QOBS and QSIM for time steps with QOBS available
    ALLOCATE(QOBS_AVAIL(NUM_AVAIL),QSIM_AVAIL(NUM_AVAIL),DOBS(NUM_AVAIL),DSIM(NUM_AVAIL),RAWD(NUM_AVAIL),LOGD(NUM_AVAIL),STAT=IERR)
    if (ierr /= 0)then
      message=trim(message)//'PROBLEM ALLOCATING SPACE FOR AVAILABLE DATA'
      return
    endif
  
    QOBS_AVAIL=PACK(QOBS,QOBS_MASK,QOBS_AVAIL)  ! moves QOBS time steps indicated by QOBS_MASK to QOBS_AVAIL,
                                                ! if no values is missing (i.e. NS = NUM_AVAIL) then QOBS_AVAIL
                                                ! should be a copy of QOBS
    QSIM_AVAIL=PACK(QSIM,QOBS_MASK,QSIM_AVAIL)  ! moves QSIM time steps indicated by QOBS_MASK to QSIM_AVAIL
                                                ! if no values is missing (i.e. NS = NUM_AVAIL) then QSIM_AVAIL
                                                ! should be a copy of QSIM
                                                                   
    ! compute mean
    XB_OBS  = SUM(QOBS_AVAIL(:)) / INT(NUM_AVAIL, KIND(WP))
    XB_SIM  = SUM(QSIM_AVAIL(:)) / INT(NUM_AVAIL, KIND(WP))
    
    ! define NO_ZERO as 1% of the observed mean flow
    NO_ZERO = XB_OBS/100
  
    ! compute the sum of squares of simulated and observed vectors
    DOBS(:) = QOBS_AVAIL(:) - XB_OBS
    DSIM(:) = QSIM_AVAIL(:) - XB_SIM
    SS_OBS  = DOT_PRODUCT(DOBS,DOBS)  ! = SUM( DOBS(:)*DOBS(:) )
    SS_SIM  = DOT_PRODUCT(DSIM,DSIM)  ! = SUM( DSIM(:)*DSIM(:) )
  
    ! compute the sum of squares of lagged differences
    SS_LOBS = DOT_PRODUCT(DOBS(2:NUM_AVAIL),DOBS(1:NUM_AVAIL-1))
    SS_LSIM = DOT_PRODUCT(DSIM(2:NUM_AVAIL),DSIM(1:NUM_AVAIL-1))
    ! compute sum of squared differences between model and observations
    RAWD(:) = QSIM_AVAIL(:) - QOBS_AVAIL(:)
    LOGD(:) = LOG(QSIM_AVAIL(:)+NO_ZERO) - LOG(QOBS_AVAIL(:)+NO_ZERO)
    SS_RAW  = DOT_PRODUCT(RAWD,RAWD)  ! = SUM( RAWD(:)*RAWD(:) )
    SS_LOG  = DOT_PRODUCT(LOGD,LOGD)  ! = SUM( LOGD(:)*LOGD(:) )
    
    
  
    ! ---------------------------------------------------------------------------------------
    ! (2) COMPUTE ERROR STATISTICS
    ! ---------------------------------------------------------------------------------------
    ! compute the mean
    work%run%stats%QOBS_MEAN = XB_OBS
    work%run%stats%QSIM_MEAN = XB_SIM
    ! compute the coefficient of variation
    work%run%stats%QOBS_CVAR = SQRT( SS_OBS / INT(NUM_AVAIL-1, KIND(WP)) ) / (XB_OBS+NO_ZERO)
    work%run%stats%QSIM_CVAR = SQRT( SS_SIM / INT(NUM_AVAIL-1, KIND(WP)) ) / (XB_SIM+NO_ZERO)
    ! compute the lag-1 correlation coefficient
    work%run%stats%QOBS_LAG1 = SS_LOBS / (SQRT(SS_OBS*SS_OBS)+NO_ZERO)
    work%run%stats%QSIM_LAG1 = SS_LSIM / (SQRT(SS_SIM*SS_SIM)+NO_ZERO)
    
    ! Compute RMSE using the metrics module
    work%run%stats%RAW_RMSE = get_RMSE(QOBS_AVAIL, QSIM_AVAIL, '1')   ! No transformation
      
    ! Compute log RMSE using the metrics module
    work%run%stats%LOG_RMSE = get_RMSE(QOBS_AVAIL, QSIM_AVAIL, 'log') ! Log transformation
  
    ! Compute NSE using the metrics module
    work%run%stats%NASH_SUTT = get_NSE(QOBS_AVAIL, QSIM_AVAIL, '1')   ! No transformation
    
    ! Compute KGE using the metrics module
    work%run%stats%KGE = get_KGE(QOBS_AVAIL, QSIM_AVAIL, '1')         ! No transformation
  
    ! Compute KGEp using the metrics module
    work%run%stats%KGEP = get_KGEP(QOBS_AVAIL, QSIM_AVAIL, '1')       ! No transformation
    
    ! Compute MAE using the metrics module
    work%run%stats%MAE = get_MAE(QOBS_AVAIL, QSIM_AVAIL, '1')         ! No transformation

    ! Compute composite KGE using the metrics module
    work%run%stats%KGECOMP = get_composite_KGE(QOBS_AVAIL, QSIM_AVAIL) ! Transformations 1 and -1
    
    ! Savethe metric chosen as objective function
    ! Note: the objective function uses the transformation (TRANSFO) requested in the control file
    select case(METRIC)
      case ("KGE");     work%run%stats%METRIC_VAL = get_KGE(QOBS_AVAIL, QSIM_AVAIL, TRANSFO)
      case ("KGEP");    work%run%stats%METRIC_VAL = get_KGEP(QOBS_AVAIL, QSIM_AVAIL, TRANSFO)
      case ("NSE");     work%run%stats%METRIC_VAL = get_NSE(QOBS_AVAIL, QSIM_AVAIL, TRANSFO)
      case ("RMSE");    work%run%stats%METRIC_VAL = get_RMSE(QOBS_AVAIL, QSIM_AVAIL, TRANSFO)
      case ("MAE");     work%run%stats%METRIC_VAL = get_MAE(QOBS_AVAIL, QSIM_AVAIL, TRANSFO)
      case ("KGECOMP"); work%run%stats%METRIC_VAL = work%run%stats%KGECOMP  ! uses fixed transformations 1 and -1
      case default
        message=trim(message)//'The requested metric is not available for calibration'
        ierr=10; return
    end select

    ! ---------------------------------------------------------------------------------------
    DEALLOCATE(QOBS,QOBS_AVAIL,QSIM,QSIM_AVAIL,DOBS,DSIM,RAWD,LOGD,STAT=IERR)
    if (ierr /= 0)then
      message=trim(message)//'PROBLEM DEALLOCATING SPACE'
      return
    endif
  
  END IF
  
  if(isPrint)then
    PRINT *, 'NSE = ',          work%run%stats%NASH_SUTT
    PRINT *, 'KGE = ',          work%run%stats%KGE
    PRINT *, 'KGEP = ',         work%run%stats%KGEP
    PRINT *, 'KGECOMP = ',      work%run%stats%KGECOMP
    PRINT *, 'MAE = ',          work%run%stats%MAE
    PRINT *, 'RAW_RMSE = ',     work%run%stats%RAW_RMSE
    PRINT *, 'LOG_RMSE = ',     work%run%stats%LOG_RMSE
    PRINT *, 'METRIC_VAL [Metric:',METRIC,' / Transfo:',TRANSFO,'] =',   work%run%stats%METRIC_VAL
  endif
 
  END SUBROUTINE compute_performance_stats
 
  ! ----------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------------

  subroutine accumulate_numerical_stats(work)
  USE model_numerix, only: NUM_FUNCS, NUM_JACOBIAN, NUMSUB_ACCEPT, NUMSUB_REJECT, NUMSUB_NOCONV
  USE model_numerix, only: MAXNUM_ITERNS, ORD_NSUBS, PRB_NSUBS
  IMPLICIT NONE
  type(fuse_work), intent(inout)    :: work  ! structures that depend on nState/nPar
  ! ----------------------------------------------------------------------------------------
  ! compute numerical stats
  work%run%stats%NUM_FUNCS     = work%run%stats%NUM_FUNCS     + REAL(NUM_FUNCS, KIND(WP))     ! number of function calls
  work%run%stats%NUM_JACOBIAN  = work%run%stats%NUM_JACOBIAN  + REAL(NUM_JACOBIAN, KIND(WP))  ! number of times Jacobian is calculated
  work%run%stats%NUMSUB_ACCEPT = work%run%stats%NUMSUB_ACCEPT + REAL(NUMSUB_ACCEPT, KIND(WP)) ! number of sub-steps accepted (taken)
  work%run%stats%NUMSUB_REJECT = work%run%stats%NUMSUB_REJECT + REAL(NUMSUB_REJECT, KIND(WP)) ! number of sub-steps tried but rejected
  work%run%stats%NUMSUB_NOCONV = work%run%stats%NUMSUB_NOCONV + REAL(NUMSUB_NOCONV, KIND(WP)) ! number of sub-steps tried that did not converge
  ! compute maximum number of iterations
  IF (MAXNUM_ITERNS > work%run%stats%MAXNUM_ITERNS) work%run%stats%MAXNUM_ITERNS = MAXNUM_ITERNS
  ! compute probability distributions
  WHERE(ORD_NSUBS.GE.NUMSUB_ACCEPT) PRB_NSUBS = PRB_NSUBS + 1
  ! ----------------------------------------------------------------------------------------
  end subroutine accumulate_numerical_stats

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  subroutine finalize_numerical_stats(work)
  USE multiforce, only: NUMTIM_SIM                      ! model forcing structure (temporally constant)
  USE model_numerix, only: PRB_NSUBS                    ! model numerix parameters and data
  type(fuse_work),   intent(inout)       :: work        ! work structures that depend on npar/nState
  ! ---------------------------------------------------------------------------------------
  ! COMPUTE STATISTICS ON NUMERICAL ACCURACY AND EFFICIENCY
  ! ---------------------------------------------------------------------------------------
  ! compute RMSE between "more accurate" and "less accurate" solutions
  work%run%stats%NUM_FUNCS     = work%run%stats%NUM_FUNCS     / REAL(NUMTIM_SIM, KIND(WP)) ! number of function calls
  work%run%stats%NUM_JACOBIAN  = work%run%stats%NUM_JACOBIAN  / REAL(NUMTIM_SIM, KIND(WP)) ! number of times Jacobian is calculated
  work%run%stats%NUMSUB_ACCEPT = work%run%stats%NUMSUB_ACCEPT / REAL(NUMTIM_SIM, KIND(WP)) ! number of sub-steps accepted (taken)
  work%run%stats%NUMSUB_REJECT = work%run%stats%NUMSUB_REJECT / REAL(NUMTIM_SIM, KIND(WP)) ! number of sub-steps tried but rejected
  work%run%stats%NUMSUB_NOCONV = work%run%stats%NUMSUB_NOCONV / REAL(NUMTIM_SIM, KIND(WP)) ! number of sub-steps tried that did not converge
  ! compute cumulative probability distributions
  work%run%stats%NUMSUB_PROB   = REAL(PRB_NSUBS(:), KIND(WP)) / REAL(NUMTIM_SIM, KIND(WP))
  ! ---------------------------------------------------------------------------------------
 end subroutine finalize_numerical_stats

end module run_statistics
