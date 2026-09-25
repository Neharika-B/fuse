FUNCTION FUNCTN(NOPT,A)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! Modified by Cyril Thébault to allow different metrics as objective function, 2024
! Modified by Cyril Thébault to allow parameter transformations, 7/2026
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Wrapper for SCE (used to compute the objective function)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE sce_callback_context, only: ctx                   ! access FUSE data structures
USE parameter_transform_module, only: vector_to_physical_space
USE fuse_evaluate_module, only: fuse_evaluate         ! run model and compute the metric chosen as objective function
USE multiforce, only: ncid_hydromet                   ! NetCDF hydromet file ID
USE fuse_fileManager,only:METRIC, TRANSFO             ! metric and transformation requested in the filemanager
USE fuse_globaldata, only: nFUSE_eval                      ! # fuse evaluations

IMPLICIT NONE
! input
INTEGER(I4B)                           :: NOPT        ! number of parameters
REAL(MSP), DIMENSION(100), INTENT(IN)  :: A            ! model parameter set - can be bumped up to 100 elements

! internal
REAL(MSP), DIMENSION(NOPT)             :: SCE_PAR_MSP ! physical parameters in SCE precision
REAL(WP), DIMENSION(NOPT)              :: SCE_PAR     ! physical parameters in FUSE precision
INTEGER(I4B)                           :: IERR        ! transformation error code
CHARACTER(LEN=256)                     :: MESSAGE     ! parameter transformation error message
LOGICAL(LGT)                           :: OUTPUT_FLAG ! .TRUE. = write model time series
REAL(WP)                               :: METRIC_VAL  ! value of the metric chosen as objective function

! output
REAL(MSP)                              :: FUNCTN      ! objective function value

! ---------------------------------------------------------------------------------------

nFUSE_eval = nFUSE_eval + 1

! Convert the optimizer vector from search space back to physical space.
IF (.NOT. ALLOCATED(ctx%transform_codes)) THEN
  STOP 'SCE parameter transformations are not initialized'
END IF

IF (SIZE(ctx%transform_codes) /= NOPT) THEN
  STOP 'Incorrect number of SCE parameter transformation codes'
END IF

CALL vector_to_physical_space(             &
  A(1:NOPT),                               &
  ctx%transform_codes,                     &
  SCE_PAR_MSP,                             &
  IERR,                                    &
  MESSAGE)

IF (IERR /= 0) THEN
  WRITE(*,'(A)') TRIM(MESSAGE)
  STOP 'Unable to transform SCE parameters to physical space'
END IF

! Convert from MSP used by SCE to WP used by FUSE.
SCE_PAR = SCE_PAR_MSP

OUTPUT_FLAG = .FALSE.  ! do not produce runs.nc files during calibration

CALL FUSE_evaluate(SCE_PAR, ctx%info, ctx%work, ctx%domain, OUTPUT_FLAG, METRIC_VAL, &
                   ierr, message)
if (ierr/=0) stop trim(message)

! save objective function value: SCE is a minimization algorithm
select case(metric)
 case ("KGE", "KGEP", "NSE", "KGECOMP"); FUNCTN = -METRIC_VAL
 case ("MAE", "RMSE");        FUNCTN = METRIC_VAL
 case default
  STOP 'The requested metric is not available in metrics module'
end select

! ---------------------------------------------------------------------------------------
END FUNCTION FUNCTN
