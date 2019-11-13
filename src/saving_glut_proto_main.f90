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

    real(dp), dimension(nt+1,nb) :: v1w_seqa, v1e_seqa, Pri1_seqa,val1_seqa
    real(dp), dimension(nt,nb) :: DebPol1_seqa
    real(dp), dimension(nb) :: vv1wmx,vv1emx,pprimx1
    real(dp), dimension(nb) :: DDebPol1
    real(dp), dimension(maxgrid) :: val1mx
    real(dp) :: b1,b1pr,p1,v1e,v1w
    real(dp), dimension(7) :: retval

    ! Auxilliary variables
    integer :: indz, indb, indbp, indt, indb1, indbp1
    integer, dimension(1) :: maxind

!=======================================================================!
!                          INITIALIZATION                               !
!=======================================================================!    

    do indz = 1,nz
        zvec(indz) = zmin+(zmax-zmin)*(real(indz-1)/real(nz-1))        
    end do
    pzvec = 1.0_dp/real(nz)
    pzvec = pzvec/sum(pzvec)

    ! bmin = 0.0_dp
    bmin = d_foreign/wgt1
    bmax = bmin+0.8_dp
    ! bmax = min(bmax,0.99_dp)

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

    v1w_seqa(indt,:) = vv1wMx
    v1e_seqa(indt,:) = vv1eMx
    Pri1_seqa(indt,:) = PPriMx1

    ! write (*,*) 'Export Results'
    ! open(1,file='./results/val1.txt',form='formatted')
    ! do indb = 1,nb
    !     write(1,'(ES14.6)') val1mx(indb)
    ! end do
    ! close(1)    

    ! open(1,file='./results/bvec.txt',form='formatted')
    ! do indb = 1,nb
    !     write(1,'(ES14.6)') bVec(indb)
    ! end do
    ! close(1)        

!=======================================================================!
!                       ITERATION STARTING AT T-1                       !
!=======================================================================!

    time: do indt = nt-1,1,-1

        ! solve optimal policy 
        do indb1 = 1,nb
            b1 = bvec(indb1)
            val1mx = 0.0_dp
            do indbp1 = 1,maxgrid
                ! TODO: here needs to be greater than Dbar as well
                b1pr = DebChoiceVec(indbp1)
                call SolveSystem(indt,b1,b1pr,RetVal)
                v1w = RetVal(5)
                v1e = RetVal(6)
                val1Mx(indbp1) = wgt1*v1w+(1-wgt1)*v1e
            end do

            ! Optimal Policy
            MaxInd = MAXLOC(val1Mx)
            DDebPol1(indb1) = DebChoiceVec(MaxInd(1))
        end do

        write (*,'(a12,i6,f17.11)') "POLICY ITER", nT-indt, &
            sum(abs(DDebPol1-DebPol1))

        DebPol1 = DDebPol1
        DebPol1_seqa(indt,:) = DebPol1

        ! compute the equilibrium at optimal policy
        do indb1 = 1,nb
            b1 = bVec(indb1)
            b1pr = DebPol1(indb1)
            call SolveSystem(indt,b1,b1pr,RetVal)
            p1 = RetVal(3)
            v1w = RetVal(5)
            v1e = RetVal(6)

            vv1wMx(indb1) = v1w
            vv1eMx(indb1) = v1e
            PPriMx1(indb1) = p1
        end do

        v1wMx = vv1wMx
        v1eMx = vv1eMx
        PriMx1 = PPriMx1

        v1w_seqa(indt,:) = v1wMx
        v1e_seqa(indt,:) = v1eMx
        Pri1_seqa(indt,:) = PriMx1
    end do time

    val1_seqa = wgt1*v1w_seqa+(1-wgt1)*v1e_seqa

    write (*,*) 'Export Results'
    open(1,file='./results/val1_seqa.txt',form='formatted')
    do indt = 1,nt+1
        write(1,'(200ES14.6)') val1_seqa(indt,:)
    end do
    close(1)    

    open(1,file='./results/DebPol1_seqa.txt',form='formatted')
    do indt = 1,nt+1
        write(1,'(200ES14.6)') DebPol1_seqa(indt,:)
    end do
    close(1)

    open(1,file='./results/bvec.txt',form='formatted')
    do indb = 1,nb
        write(1,'(ES14.6)') bVec(indb)
    end do
    close(1)        

end program saving_glut_gov_main