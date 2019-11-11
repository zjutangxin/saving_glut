! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - Main function for one country with fixed foreign demand.
!     - Saving_glut Project
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

program saving_glut_gov_main

    use global
    implicit none

    real(dp) :: bmin, bmax

    real(dp), dimension(nt+1,nb) :: v1w_seqa, v1e_seqa, Pri1_seqa
    real(dp), dimension(nt,nb) :: DebPol1_seqa
    real(dp), dimension(nb) :: vv1wmx,vv1emx,pprimx1,val1mx
    real(dp) :: b1,b1pr,p1,v1e,v1w
    real(dp), dimension(7) :: retval

    ! Auxilliary variables
    integer :: indz, indb, indbp, indt, indb1

!=======================================================================!
!                          INITIALIZATION                               !
!=======================================================================!    

    do indz = 1,nz
        zvec(indz) = zmin+(zmax-zmin)*(real(indz-1)/real(nz-1))        
    end do
    pzvec = 1.0_dp/real(nz)
    pzvec = pzvec/sum(pzvec)

    bmin = d_foreign/wgt1
    bmax = bmin+0.8_dp

    do indb = 1, nb
        bvec(indb) = bmin+(bmax-bmin)*(real(indb-1)/real(nb-1))
    end do

    do indbp = 1, maxGrid
        DebChoiceVec(indbp) = &
            bmin+(bmax-bmin)*(real(indbp-1)/real(maxGrid-1))
    end do    
    
    v1wMx = 0.0_dp
    v1eMx = 0.0_dp
    priMx1 = 0.0_dp

    v1w_seqa(nt+1,:) = v1wMx
    v1e_seqa(nt+1,:) = v1eMx
    Pri1_seqa(nt+1,:) = PriMx1

!=======================================================================!
!                          TERMINAL PERIOD                              !
!=======================================================================!

    indt = nt

    DebPol1 = 0.0_dp
    DebPol1_seqa(indt,:) = DebPol1

    do indb1 = 1,nb
        b1 = bvec(indb1)
        b1pr = DebPol1(indb1)
        call SolveSystem(indt,b1,b1pr,retval)
        p1 = RetVal(3)
        v1w = RetVal(5)
        v1e = RetVal(6)
        vv1wMx(indb1) = v1w
        vv1eMx(indb1) = v1e
        PPriMx1(indb1) = p1
    end do

    val1mx = wgt1*vv1wMx + (1-wgt1)*vv1eMx
    v1w_seqa(indt,:) = vv1wMx
    v1e_seqa(indt,:) = vv1eMx
    Pri1_seqa(indt,:) = PPriMx1

    write (*,*) 'Export Results'
    open(1,file='./results/val1.txt',form='formatted')
    do indb = 1,nb
        write(1,'(ES14.6)') val1mx(indb)
    end do
    close(1)    

    open(1,file='./results/bvec.txt',form='formatted')
    do indb = 1,nb
        write(1,'(ES14.6)') bVec(indb)
    end do
    close(1)        

end program saving_glut_gov_main