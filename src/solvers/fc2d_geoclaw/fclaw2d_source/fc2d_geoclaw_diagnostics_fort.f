c    # ----------------------------------------------------------------------------------
c    # Output and diagnostics
c    # ----------------------------------------------------------------------------------
      subroutine fc2d_geoclaw_fort_conservation_check(mx,my,mbc,meqn,
     &      dx,dy,area,q,sum,c_kahan)
      implicit none

      integer mx,my,mbc,meqn
      double precision dx, dy, dxdy
      double precision sum(meqn), c_kahan
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision t, y

      !!include 'fclaw2d_metric_terms.i'
      double precision  area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j,m
      integer*8 cont, fclaw_map_get_context
      logical fclaw_map_is_used

      cont = fclaw_map_get_context()

      dxdy = dx*dy
      do m = 1,meqn
         if (fclaw_map_is_used(cont)) then
            do j = 1,my
               do i = 1,mx
                  y = q(i,j,m)*area(i,j) - c_kahan
                  t = sum(m) + y
                  c_kahan = (t-sum(m)) - y
                  sum(m) = t
c                  sum(m) = sum(m) + q(m,i,j)*area(i,j)
               enddo
            enddo
         else
            do j = 1,my
               do i = 1,mx
                  sum(m) = sum(m) + q(m,i,j)*dx*dy
               enddo
            enddo
         endif
      enddo

      end


c     # Compute area of a patch
      double precision function
     &      fc2d_geoclaw_fort_compute_patch_area(mx,my,
     &      mbc,dx,dy,area)
      implicit none

      integer mx,my, mbc
      double precision dx, dy
      double precision sum

      double precision  area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j
      integer*8 cont, fclaw_map_get_context
      logical fclaw_map_is_used

      cont = fclaw_map_get_context()

      if (fclaw_map_is_used(cont)) then
         sum = 0
         do j = 1,my
            do i = 1,mx
               sum = sum + area(i,j)
            enddo
         enddo
      else
         sum = dx*dy*mx*my
      endif

      fc2d_geoclaw_fort_compute_patch_area = sum

      end


      subroutine fc2d_geoclaw_fort_compute_error_norm(blockno, 
     &      mx,my,mbc,meqn, dx,dy,area,error,error_norm)
      implicit none

      integer blockno, mx,my,mbc,meqn
      double precision dx, dy, dxdy, eij
      double precision error_norm(meqn,3)
      double precision error(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      !!include 'fclaw2d_metric_terms.i'
      double precision  area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j,m
      integer*8 cont, fclaw_map_get_context
      logical fclaw_map_is_used

      cont = fclaw_map_get_context()

c     # error_norm(:) comes in with values;  do not initialize it here!
      dxdy = dx*dy
      do m = 1,meqn
         if (fclaw_map_is_used(cont)) then
            do j = 1,my
               do i = 1,mx
                  eij = abs(error(m,i,j))
                  error_norm(m,1) = error_norm(m,1) +
     &                  eij*area(i,j)
                  error_norm(m,2) = error_norm(m,2) +
     &                  eij**2*area(i,j)
                  error_norm(m,3) = max(eij,error_norm(m,3))
               enddo
            enddo
         else
            do j = 1,my
               do i = 1,mx
                  eij = abs(error(m,i,j))
                  error_norm(m,1) = error_norm(m,1) +
     &                  eij*dxdy
                  error_norm(m,2) = error_norm(m,2) +
     &                  eij**2*dxdy
                  error_norm(m,3) = max(eij,error_norm(m,3))
               enddo
            enddo
         endif
      enddo

      end
