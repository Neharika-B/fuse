MODULE SUMEXTRACT_MODULE
IMPLICIT NONE
CONTAINS
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
PURE FUNCTION SUMEXTRACT(stats, STATNAME)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! Modified by Cyril Thébault to allow different metrics as objective function, 2024
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Extracts variable "VNAME(IVAR)" from relevant data structures
! ---------------------------------------------------------------------------------------
USE nrtype
USE multistats_types, only: SUMMARY
IMPLICIT NONE
! input
type(SUMMARY), intent(in)              :: stats       ! structures that depend on nState/nPar
CHARACTER(*), INTENT(IN)               :: STATNAME    ! variable name
! internal
REAL(WP)                               :: XVAR        ! variable
! output
REAL(WP)                               :: SUMEXTRACT  ! FUNCTION name
! ---------------------------------------------------------------------------------------
! initialize XVAR
XVAR=-9999._wp
! extract summary statistics
IF (TRIM(STATNAME).EQ.'qobs_mean')   XVAR = stats%QOBS_MEAN
IF (TRIM(STATNAME).EQ.'qsim_mean')   XVAR = stats%QSIM_MEAN 
IF (TRIM(STATNAME).EQ.'qobs_cvar')   XVAR = stats%QOBS_CVAR
IF (TRIM(STATNAME).EQ.'qsim_cvar')   XVAR = stats%QSIM_CVAR
IF (TRIM(STATNAME).EQ.'qobs_lag1')   XVAR = stats%QOBS_LAG1
IF (TRIM(STATNAME).EQ.'qsim_lag1')   XVAR = stats%QSIM_LAG1
IF (TRIM(STATNAME).EQ.'raw_rmse')    XVAR = stats%RAW_RMSE
IF (TRIM(STATNAME).EQ.'log_rmse')    XVAR = stats%LOG_RMSE
IF (TRIM(STATNAME).EQ.'nash_sutt')   XVAR = stats%NASH_SUTT
IF (TRIM(STATNAME).EQ.'kge')         XVAR = stats%KGE
IF (TRIM(STATNAME).EQ.'kgep')        XVAR = stats%KGEP
IF (TRIM(STATNAME).EQ.'kgecomp')     XVAR = stats%KGECOMP
IF (TRIM(STATNAME).EQ.'mae')         XVAR = stats%MAE
IF (TRIM(STATNAME).EQ.'metric_val')  XVAR = stats%METRIC_VAL
! extract numerix stats
IF (TRIM(STATNAME).EQ.'numerx_rmse') XVAR = stats%NUM_RMSE      
IF (TRIM(STATNAME).EQ.'mean_nfuncs') XVAR = stats%NUM_FUNCS
IF (TRIM(STATNAME).EQ.'mean_njacob') XVAR = stats%NUM_JACOBIAN
IF (TRIM(STATNAME).EQ.'mean_accept') XVAR = stats%NUMSUB_ACCEPT
IF (TRIM(STATNAME).EQ.'mean_reject') XVAR = stats%NUMSUB_REJECT
IF (TRIM(STATNAME).EQ.'mean_noconv') XVAR = stats%NUMSUB_NOCONV
IF (TRIM(STATNAME).EQ.'maxnum_iter') XVAR = REAL(stats%MAXNUM_ITERNS, KIND(WP))
! and, save the output
SUMEXTRACT = XVAR
! ---------------------------------------------------------------------------------------
END FUNCTION SUMEXTRACT
END MODULE SUMEXTRACT_MODULE
