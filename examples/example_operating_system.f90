!> \file example_operating_system.f90
!! \brief Example program demonstrating the OPERATING_SYSTEM_LINUX module functionalities
!!
!! This program shows how to use all the available functions in the operating_system module:
!! - Directory management (create, delete)
!! - File operations (check existence, remove)
!! - Directory listing and file searching
!! - System information queries
!! - String utilities
!!
!! \author Alessandro
!! \date 2025

PROGRAM EXAMPLE_OPERATING_SYSTEM
  USE OPERATING_SYSTEM_LINUX
  IMPLICIT NONE

  ! Variables for demonstration
  CHARACTER(LEN=256) :: current_dir
  CHARACTER(LEN=256) :: test_dir = "test_example_dir"
  CHARACTER(LEN=256) :: test_file = "example_test.txt"
  CHARACTER(LEN=256) :: os_name
  CHARACTER(LEN=256) :: found_file
  CHARACTER(LEN=256), ALLOCATABLE :: files(:)
  CHARACTER(LEN=100) :: test_string = "Hello World EXAMPLE"
  CHARACTER(LEN=100) :: pattern = "EXAMPLE"
  CHARACTER(LEN=100), ALLOCATABLE :: matches(:)
  INTEGER :: n_matches, count, ios, i

  WRITE(*,*) "=========================================="
  WRITE(*,*) "   OPERATING SYSTEM MODULE EXAMPLE"
  WRITE(*,*) "=========================================="
  WRITE(*,*)

  ! 1. Get system information
  WRITE(*,*) "1. SYSTEM INFORMATION:"
  WRITE(*,*) "----------------------"
  
  ! Get operating system name
  os_name = GET_OS_NAME()
  WRITE(*,*) "Operating System: ", TRIM(os_name)
  
  ! Get current working directory (two different methods)
  CALL GET_CURRENT_DIRECTORY(current_dir)
  WRITE(*,*) "Current Directory (method 1): ", TRIM(current_dir)
  
  current_dir = GET_CURRENT_WORKING_DIRECTORY()
  WRITE(*,*) "Current Directory (method 2): ", TRIM(current_dir)
  WRITE(*,*)

  ! 2. Directory operations
  WRITE(*,*) "2. DIRECTORY OPERATIONS:"
  WRITE(*,*) "------------------------"
  
  ! Create a test directory
  WRITE(*,*) "Creating test directory: ", TRIM(test_dir)
  CALL CREATE_DIRECTORY(test_dir)
  WRITE(*,*) "Directory created successfully!"
  
  ! Create a subdirectory
  CALL CREATE_DIRECTORY(TRIM(test_dir)//"/subdir")
  WRITE(*,*) "Subdirectory created: ", TRIM(test_dir)//"/subdir"
  WRITE(*,*)

  ! 3. File operations
  WRITE(*,*) "3. FILE OPERATIONS:"
  WRITE(*,*) "-------------------"
  
  ! Create a test file
  OPEN(UNIT=10, FILE=TRIM(test_dir)//"/"//TRIM(test_file), STATUS='NEW', ACTION='WRITE', IOSTAT=ios)
  IF (ios == 0) THEN
    WRITE(10,*) "This is a test file created by the operating system example."
    WRITE(10,*) "Line 2: Testing file operations"
    WRITE(10,*) "Line 3: Example content"
    CLOSE(10)
    WRITE(*,*) "Test file created: ", TRIM(test_dir)//"/"//TRIM(test_file)
  ELSE
    WRITE(*,*) "Error creating test file"
  END IF
  
  ! Check if file exists
  IF (FILE_EXISTS(TRIM(test_dir)//"/"//TRIM(test_file))) THEN
    WRITE(*,*) "File existence check: PASSED"
  ELSE
    WRITE(*,*) "File existence check: FAILED"
  END IF
  
  ! Create another test file
  OPEN(UNIT=11, FILE=TRIM(test_dir)//"/another_file.dat", STATUS='NEW', ACTION='WRITE', IOSTAT=ios)
  IF (ios == 0) THEN
    WRITE(11,*) "Another test file"
    CLOSE(11)
    WRITE(*,*) "Additional file created: ", TRIM(test_dir)//"/another_file.dat"
  END IF
  WRITE(*,*)

  ! 4. Directory listing
  WRITE(*,*) "4. DIRECTORY LISTING:"
  WRITE(*,*) "---------------------"
  
  ! List all files in the test directory
  ALLOCATE(files(100))  ! Allocate space for up to 100 files
  CALL LIST_FILES_IN_DIRECTORY(test_dir, files, count)
  
  WRITE(*,*) "Files in directory '", TRIM(test_dir), "':"
  DO i = 1, count
    WRITE(*,*) "  - ", TRIM(files(i))
  END DO
  WRITE(*,*) "Total files found: ", count
  WRITE(*,*)
  
  ! List only .txt files
  CALL LIST_FILES_IN_DIRECTORY(test_dir, files, count, "txt")
  WRITE(*,*) "Only .txt files:"
  DO i = 1, count
    WRITE(*,*) "  - ", TRIM(files(i))
  END DO
  WRITE(*,*)

  ! 5. File searching
  WRITE(*,*) "5. FILE SEARCHING:"
  WRITE(*,*) "------------------"
  
  ! Find a specific file
  found_file = FIND_FILE(test_file, test_dir)
  IF (TRIM(found_file) /= "") THEN
    WRITE(*,*) "File found: ", TRIM(found_file)
  ELSE
    WRITE(*,*) "File not found"
  END IF
  WRITE(*,*)

  ! 6. String utilities
  WRITE(*,*) "6. STRING UTILITIES:"
  WRITE(*,*) "--------------------"
  
  ! Case-insensitive string search
  ALLOCATE(matches(10))
  CALL FIND_STRING_CASE_INSENSITIVE(test_string, pattern, matches, n_matches)
  
  WRITE(*,*) "Original string: '", TRIM(test_string), "'"
  WRITE(*,*) "Search pattern: '", TRIM(pattern), "'"
  WRITE(*,*) "Number of matches: ", n_matches
  DO i = 1, n_matches
    WRITE(*,*) "  Match found: '", TRIM(matches(i)), "'"
  END DO
  
  ! Convert to lowercase
  CALL TO_LOWERCASE(test_string)
  WRITE(*,*) "Lowercase version: '", TRIM(test_string), "'"
  WRITE(*,*)

  ! 7. Cleanup operations
  WRITE(*,*) "7. CLEANUP OPERATIONS:"
  WRITE(*,*) "----------------------"
  
  ! Remove individual file
  CALL REMOVE_FILE(TRIM(test_dir)//"/another_file.dat")
  WRITE(*,*) "Removed file: another_file.dat"
  
  ! Check if file still exists
  IF (.NOT. FILE_EXISTS(TRIM(test_dir)//"/another_file.dat")) THEN
    WRITE(*,*) "File removal confirmed"
  END IF
  
  ! Delete the entire test directory
  CALL DELETE_DIRECTORY(test_dir)
  WRITE(*,*) "Deleted test directory: ", TRIM(test_dir)
  WRITE(*,*)

  ! 8. Final verification
  WRITE(*,*) "8. FINAL VERIFICATION:"
  WRITE(*,*) "----------------------"
  
  ! Check if directory was deleted
  IF (.NOT. FILE_EXISTS(test_dir)) THEN
    WRITE(*,*) "Directory cleanup confirmed - test directory no longer exists"
  ELSE
    WRITE(*,*) "Warning: Test directory still exists"
  END IF

  WRITE(*,*)
  WRITE(*,*) "=========================================="
  WRITE(*,*) "   EXAMPLE COMPLETED SUCCESSFULLY!"
  WRITE(*,*) "=========================================="

  ! Cleanup allocated arrays
  DEALLOCATE(files)
  DEALLOCATE(matches)

END PROGRAM EXAMPLE_OPERATING_SYSTEM