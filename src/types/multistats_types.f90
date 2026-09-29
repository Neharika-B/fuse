MODULE multistats_types

 USE nrtype

 implicit none
 private

 public :: SUMMARY
 
 ! --------------------------------------------------------------------------------------

 TYPE SUMMARY

  ! comparisons between model output and observations
  REAL(WP)                             :: QOBS_MEAN     ! mean observed runoff (mm day-1)
  REAL(WP)                             :: QSIM_MEAN     ! mean simulated runoff (mm day-1)
  REAL(WP)                             :: QOBS_CVAR     ! coefficient of variation of observed runoff (-)
  REAL(WP)                             :: QSIM_CVAR     ! coefficient of variation of simulated runoff (-)
  REAL(WP)                             :: QOBS_LAG1     ! lag-1 correlation of observed runoff (-)
  REAL(WP)                             :: QSIM_LAG1     ! lag-1 correlation of simulated runoff (-)
  REAL(WP)                             :: RAW_RMSE      ! root-mean-squared-error of flow (mm day-1)
  REAL(WP)                             :: LOG_RMSE      ! root-mean-squared-error of LOG flow (mm day-1)
  REAL(WP)                             :: NASH_SUTT     ! Nash-Sutcliffe score
  REAL(WP)                             :: KGE           ! Kling-Gupta Efficiency score
  REAL(WP)                             :: KGEP          ! Kling-Gupta Efficiency' score
  REAL(WP)                             :: KGECOMP       ! Kling-Gupta Efficiency composite score
  REAL(WP)                             :: KGECSI        ! KGE-CSI score: F0*CSI + (1-F0)*KGE_NZ
  REAL(WP)                             :: CSI           ! critical success index for zero flow (part of KGE-CSI)
  REAL(WP)                             :: KGE_NZ        ! KGE on non-zero observed flows (part of KGE-CSI)
  REAL(WP)                             :: F0            ! fraction of zero-flow observations (part of KGE-CSI)
  REAL(WP)                             :: ZERO_THR      ! zero-flow threshold (part of KGE-CSI)
  REAL(WP)                             :: MAE           ! Mean absolute error
  REAL(WP)                             :: METRIC_VAL    ! value of the metric chosen as objective function
 
  ! attributes of model output
  REAL(WP)                             :: NUM_RMSE      ! error of the approximate solution
  REAL(WP)                             :: NUM_FUNCS     ! number of function calls
  REAL(WP)                             :: NUM_JACOBIAN  ! number of times Jacobian is calculated
  REAL(WP)                             :: NUMSUB_ACCEPT ! number of sub-steps taken
  REAL(WP)                             :: NUMSUB_REJECT ! number of sub-steps taken
  REAL(WP)                             :: NUMSUB_NOCONV ! number of sub-steps tried that did not converge
  INTEGER(I4B)                         :: MAXNUM_ITERNS ! maximum number of iterations in implicit scheme
  REAL(WP), DIMENSION(20)              :: NUMSUB_PROB   ! probability distribution for number of sub-steps
 
  ! error checking
  CHARACTER(LEN=1024)                  :: ERR_MESSAGE   ! error message
 
 ENDTYPE SUMMARY

END MODULE multistats_types
