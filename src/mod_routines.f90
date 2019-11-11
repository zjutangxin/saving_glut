! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - Module with auxilliary routines.
! Author:
!     Xin Tang @ IMF, Summer 2019
!  
! Record of Revisions:
!         Date:                 Description of Changes
!     ===========        =================================
!      07/23/2019:          Original Code: No paralleling
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

module routines
    implicit none
    integer, parameter, private :: dp = kind(1.0d0)
    
    contains

    function FFindInt(xx,xxVec,yyVec,nn)
    ! one-dimension interpolation
        implicit none
        integer, intent(in) :: nn

        real(dp), intent(in) :: xx
        real(dp), dimension(nn) :: xxvec, yyvec
        real(dp) :: FFindInt

        integer :: ind
        real(dp) :: slope

        if (xx <= xxVec(1)) then
            slope    = (yyVec(2)-yyVec(1))/(xxVec(2)-xxVec(1))
            FFindInt = yyVec(1)+slope*(xx-xxVec(1))
        elseif (xx >= xxVec(nn)) then
            slope    = (yyVec(nn)-yyVec(nn-1))/(xxVec(nn)-xxVec(nn-1))
            FFindInt = yyVec(nn)+slope*(xx-xxVec(nn))
        else
            ind = 1
            do while (xx > xxVec(ind))
                ind = ind+1
            end do
            slope    = (yyVec(ind)-yyVec(ind-1))/(xxVec(ind)-xxVec(ind-1))
            FFindInt = yyVec(ind-1)+slope*(xx-xxVec(ind-1))
        endif
    end function FFindInt

    function Afun(zz,pp)
        ! returns to capital
            use parameter, only: theta, zbar
            implicit none
    
            real(dp), intent(in):: zz, pp
            real(dp) :: afun
    
            afun = theta*zz/(zbar**(1-theta))+(zz-zbar)*pp/zbar
        end function Afun
    
        function phifun(pppr,bbpr)
        ! risky share
            use parameter, only: nz
            use global, only: zvec, pzvec
            implicit none
    
            real(dp), intent(in) :: pppr, bbpr
            real(dp) :: phifun
                
            integer :: indz
            real(dp) :: zzpr
    
            phifun = 0.0_dp
            do indz = 1, nz
                zzpr = zVec(indz)
                phifun = phifun + &
                    ((Afun(zzpr,pppr)+pppr)/(Afun(zzpr,pppr)+pppr+bbpr)) &
                    *pzVec(indz)
            end do
        end function phifun

end module routines