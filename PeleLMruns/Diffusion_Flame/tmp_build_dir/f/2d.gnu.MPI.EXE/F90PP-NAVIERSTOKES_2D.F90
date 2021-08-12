













  

















































module navierstokes_2d_module
  
  implicit none

  private 

  public fort_maxval, cen2edg, FORT_AVERAGE_EDGE_STATES
  
contains

!c :: ----------------------------------------------------------
!c :: MAXVAL
!c ::             maxval = max{ rho(i,j) }
!c ::
!c :: ----------------------------------------------------------
!c ::
     subroutine fort_maxval(rho,rho_l1, rho_l2, rho_h1, rho_h2,grid_l1, grid_l2, grid_h1, grid_h2,mxval) &
          bind(C,name="fort_maxval")

       implicit none
       integer rho_l1, rho_l2, rho_h1, rho_h2
       integer grid_l1, grid_l2, grid_h1, grid_h2
       DOUBLE PRECISION  rho(rho_l1:rho_h1, rho_l2:rho_h2)
       DOUBLE PRECISION  mxval

       integer i,j

       mxval = -Huge(0.0D0)

       do j = grid_l2, grid_h2
          do i = grid_l1, grid_h1
             mxval = max(mxval, rho(i,j))
          end do
       end do

     end subroutine fort_maxval

!c ::
!c :: ----------------------------------------------------------
!c :: This routine fills an edge-centered fab from a cell-centered
!c :: fab using simple linear interpolation.
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  lo,hi      => index limits of the of the cell-centered fab
!c ::  cfab_l1, cfab_l2, cfab_h1, cfab_h2 => index limits of the cell-centered fab
!c ::  cfab       => cell-centered data
!c ::  efab_l1, efab_l2, efab_h1, efab_h2 => index limits of the edge-centered fab
!c ::  efab       => edge-centered fab to fill
!c ::  n!c         => Number of components in the fab to fill
!c ::  dir        => direction data needs to be shifted to get to edges
!c :: ----------------------------------------------------------
!c ::
      subroutine cen2edg(lo, hi, &
          cfab_l1, cfab_l2, cfab_h1, cfab_h2, cfab,&
          efab_l1, efab_l2, efab_h1, efab_h2, efab, nc, dir,&
          isharm) bind(C,name="cen2edg")
      implicit none
      integer lo(2), hi(2), nc, dir, isharm
      integer cfab_l1, cfab_l2, cfab_h1, cfab_h2
      integer efab_l1, efab_l2, efab_h1, efab_h2
      DOUBLE PRECISION  cfab(cfab_l1:cfab_h1, cfab_l2:cfab_h2, nc)
      DOUBLE PRECISION  efab(efab_l1:efab_h1, efab_l2:efab_h2, nc)

      integer i,j,n

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,n) = 0.5D0*(cfab(i,j,n) + cfab(i-1,j,n))
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,n) = 0.5D0*(cfab(i,j,n) + cfab(i,j-1,n))
                  end do
               end do
            end do
         end if
      else
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if((cfab(i,j,n) * cfab(i-1,j,n)).gt.0.0D0)then
                        efab(i,j,n)&
                            = 2*(cfab(i,j,n) * cfab(i-1,j,n))/&
                            (cfab(i,j,n) + cfab(i-1,j,n))
                     else
                        efab(i,j,n)=0.0D0
                     endif
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if((cfab(i,j,n) * cfab(i,j-1,n)).gt.0.0D0)then
                        efab(i,j,n)&
                            = 2*(cfab(i,j,n) * cfab(i,j-1,n))/&
                            (cfab(i,j,n) + cfab(i,j-1,n))
                     else
                        efab(i,j,n)=0.0D0
                     endif
                  end do
               end do
            end do
         end if
      end if
    end subroutine cen2edg

!c
!c
!c ::: -----------------------------------------------------------
!c
!c     This routine averages the mac face velocities for makeforce at 0.5D0 time

   subroutine FORT_AVERAGE_EDGE_STATES( vel, v_lo, v_hi,&
                                        umacx, ux_lo, ux_hi,&
                                        umacy, uy_lo, uy_hi,&
                                        getForceVerbose)&
                                        bind(C, name="FORT_AVERAGE_EDGE_STATES")

      implicit none

      integer :: v_lo(3), v_hi(3)
      integer :: ux_lo(3), ux_hi(3)
      integer :: uy_lo(3), uy_hi(3)
      integer :: getForceVerbose
      DOUBLE PRECISION, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3), 2) :: vel
      DOUBLE PRECISION, dimension(ux_lo(1):ux_hi(1),ux_lo(2):ux_hi(2),ux_lo(3):ux_hi(3)) :: umacx
      DOUBLE PRECISION, dimension(uy_lo(1):uy_hi(1),uy_lo(2):uy_hi(2),uy_lo(3):uy_hi(3)) :: umacy

      DOUBLE PRECISION  :: velmin(3)
      DOUBLE PRECISION  :: velmax(3)
      integer :: isioproc

      integer :: i, j, k, n

      do n = 1, 2
         velmin(n) = 1.d234
         velmax(n) = -1.d234
      enddo

      do k = v_lo(3), v_hi(3)
         do j = v_lo(2), v_hi(2)
            do i = v_lo(1), v_hi(1)
               vel(i,j,k,1) = 0.5D0*(umacx(i,j,k)+umacx(i+1,j,k))
               vel(i,j,k,2) = 0.5D0*(umacy(i,j,k)+umacy(i,j+1,k))
               do n = 1, 2
                  velmin(n) = min(velmin(n),vel(i,j,k,n))
                  velmax(n) = max(velmax(n),vel(i,j,k,n))
               enddo
            enddo
         enddo
      enddo

      if (getForceVerbose.gt.0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc.eq.1) then
            do n = 1, 2
               write (6,*) "mac velmin (",n,") = ",velmin(n)
               write (6,*) "mac velmax (",n,") = ",velmax(n)
            enddo
         endif
      endif

   end subroutine FORT_AVERAGE_EDGE_STATES

  end module navierstokes_2d_module
