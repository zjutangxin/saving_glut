! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - Module with shared variables.
! Author:
!     Xin Tang @ IMF, Summer 2019
!  
! Record of Revisions:
!         Date:                 Description of Changes
!     ===========        =================================
!      11/11/2019:          Original Code: No paralleling
!
! Compiling Environment:
!   GNU gfortran on Ubuntu 16.04
!
! Library Used:
!   - MINPACK: by source code
!   - LAPACK, ATLAS, and BLAS: by binary library
! 
! Shared by:
!   - debt_main.f90
! =============================================================================

module global

    use parameter
    implicit none

    real(dp), dimension(nb) :: bvec
    real(dp), dimension(nz) :: zvec, pzvec
    real(dp), dimension(maxgrid) :: DebChoiceVec

    real(dp), dimension(nb) :: v1wMx, v1eMx, priMx1, DebPol1

end module global