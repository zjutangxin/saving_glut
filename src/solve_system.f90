! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - Main subroutine implementing Property 1
! Author:
!     Xin Tang @ IMF, Winter 2019
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
! =============================================================================

subroutine SolveSystem(t,b1,b1pr,retval)

    use parameter
    use global
    use routines
    implicit none

    integer, intent(in) :: t
    real(dp), intent(in) :: b1,b1pr
    real(dp), dimension(7), intent(out) :: retval
    real(dp), parameter :: x_min = 1e-9

    real(dp) :: be1,be1pr,p1,p1pr,phi1,a1,r1,c1_e,EU1_e,c1_w,U1_w
    real(dp) :: v1wpr, v1epr
    real(dp) :: alfa_t, eta_t
    integer :: indz

    ! TODO: here the last period should be treated separately
    be1 = wgt1*b1 - d_foreign
    be1pr = wgt1*b1pr - d_foreign

    alfa_t = (1.0_dp-bbeta**(nT-t+1))/(1.0_dp-bbeta)
    eta_t  = (bbeta-bbeta**(nT-t+1))/(1.0_dp-bbeta**(nT-t+1))

    p1pr = FFindInt(b1pr,bVec,PriMx1,Nb)
    phi1 = phifun(p1pr,be1pr)
    p1   = eta_t*phi1*(Afun(zbar,p1pr)+be1)/(1.0-eta_t*phi1)
    
    if (t == nT) then
        R1    = 1.0_dp
        c1_e  = Afun(zbar,p1) + be1
        EU1_e = log(1.0_dp-eta_t)
        do indz=1, Nz
            a1 = max(Afun(zVec(indz),p1)+p1+be1,x_min)
            EU1_e = EU1_e + log(a1)*pzVec(indz)
        end do
    else
        R1    = (1.0_dp-eta_t*phi1)*be1pr/(eta_t*(1-phi1)*(Afun(zbar,p1)+be1))
        c1_e  = ((1.0_dp-eta_t)/eta_t)*(p1+be1pr/R1)
        EU1_e = log(1-eta_t)+(alfa_t-1)*log(eta_t*phi1/p1)
        do indz=1, Nz
            a1 = max(Afun(zVec(indz),p1)+p1+be1,x_min)
            EU1_e = EU1_e+alfa_t*log(a1)*pzVec(indz)
        end do
    endif

    c1_w = wbar+wgt1*(b1pr/R1-b1)
    U1_w = log(max(c1_w,x_min))
    
    V1wpr = FFindInt(b1pr,bVec,v1wMx,Nb)
    V1epr = FFindInt(b1pr,bVec,v1eMx,Nb)    

    !----CHECK IF PRICES ARE POSITIVE----!
    if (p1 < 0) then
        EU1_e = -10000.0_dp
        U1_w  = -10000.0_dp
    endif

    ! send the equilibrium values back to main
    retval(1)  = c1_e
    retval(2)  = c1_w
    retval(3)  = p1
    retval(4)  = R1
    retval(5)  = U1_w + delta*V1wpr
    retval(6)  = EU1_e + bbeta*V1epr
    retval(7)  = p1pr

end subroutine SolveSystem