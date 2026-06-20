subroutine dm_heating_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid,info
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call dmheatfine1(ind_grid,ngrid,ilevel)
  end do

111 format('   Entering dm_heating_fine for level',i2)

end subroutine dm_heating_fine
!###########################################################
!###########################################################


subroutine dmheatfine1(ind_grid,ngrid,ilevel)
   use amr_commons
   use hydro_commons
   use poisson_commons, ONLY: rho_dm
   implicit none

   integer::ilevel,ngrid
   integer,dimension(1:nvector)::ind_grid
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   integer::i,ind,iskip,idim,nleaf
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   real(kind=8)::dtcool,nISM,nCOM,gamma0,eps0,scale_eps
   integer,dimension(1:nvector),save::ind_cell,ind_leaf
   real(kind=8),dimension(1:nvector),save::nH,T2,delta_T2,ekk,T2min,Zsolar

   real(kind=8)::dx,dx_loc,scale,alpha_dx2
   real(kind=8),dimension(1:3)::skip_loc

   ! Conversion factor from user units to cgs units
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

#ifndef DM_HEATING_GAMMA
#define DM_HEATING_GAMMA 1d-27
#endif
   gamma0 = DM_HEATING_GAMMA * scale_t
   scale_eps = (scale_l/scale_t)**2
   eps0 = 9e20 / scale_eps


   ! Loop over cells
   do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do i=1,ngrid
         ind_cell(i)=iskip+ind_grid(i)
      end do

      ! Gather leaf cells
      nleaf=0
      do i=1,ngrid
         if(son(ind_cell(i))==0)then
            nleaf=nleaf+1
            ind_leaf(nleaf)=ind_cell(i)
         end if
      end do

      ! Compute dark matter density
      do i=1,nleaf
         nH(i)=MAX(rho_dm(ind_leaf(i)), 0.0d0)
      end do

      ! Compute net energy sink
         do i=1,nleaf
            delta_T2(i) = nH(i)*eps0 * gamma0 * dtnew(ilevel)
         end do

      ! Update total fluid energy
         do i=1,nleaf
            T2(i) = uold(ind_leaf(i),ndim+2)
         end do
         do i=1,nleaf
            T2(i) = T2(i)+delta_T2(i)
         end do

         do i=1,nleaf
            uold(ind_leaf(i),ndim+2) = T2(i)
         end do


   end do
   ! End loop over cells

end subroutine dmheatfine1


