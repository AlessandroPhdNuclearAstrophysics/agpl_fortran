!! \author Alessandro
!! \date 2025

PROGRAM EXAMPLE_QUANTUM_NUMBERS
  USE QUANTUM_NUMBERS
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: LMAX = 6, JMAX = 10, TZ = 0
  TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: CHANNELS(:)
  INTEGER :: ICH, NCH

  WRITE(*,*) "=========================================="
  WRITE(*,*) "   QUANTUM NUMBERS MODULE EXAMPLE"
  WRITE(*,*) "=========================================="
  WRITE(*,*)
  WRITE(*,*) "Generating physical scattering channels..."
  WRITE(*,*) "LMAX =", LMAX, ", JMAX =", JMAX, ", TZ =", TZ
  WRITE(*,*)

  CALL PREPARE_CHANNELS(LMAX, JMAX, TZ, CHANNELS)
  NCH = SIZE(CHANNELS)
  
  WRITE(*,*) "Generated", NCH, "physical scattering channels:"
  WRITE(*,*)
  
  DO ICH = 1, NCH
    WRITE(*,*) "Channel", ICH, ":"
    CALL PRINT_SCATTERING_CHANNEL(CHANNELS(ICH), 6)
    WRITE(*,*)
  END DO

  WRITE(*,*) "=========================================="
  WRITE(*,*) "   EXAMPLE COMPLETED SUCCESSFULLY!"
  WRITE(*,*) "=========================================="

END PROGRAM EXAMPLE_QUANTUM_NUMBERS