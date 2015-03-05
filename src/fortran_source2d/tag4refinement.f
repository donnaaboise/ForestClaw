c     # Template function for setting refinement criteria.  The
c     # user can copy this file to their directory, and then set the
c     # vt.fort_tag4refinement = &tag4refinement.

      subroutine tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq,m
      double precision xc,yc, qmin, qmax
      double precision dq, dqi, dqj

      tag_patch = 0

c     # Refine based only on first variable in system.
      do mq = 1,meqn
         qmin = q(1,1,mq)
         qmax = q(1,1,mq)
         do j = 1,my
            do i = 1,mx
               qmin = min(q(i,j,mq),qmin)
               qmax = max(q(i,j,mq),qmax)
               if (qmax - qmin .gt. tag_threshold) then
                  tag_patch = 1
                  return
               endif
            enddo
         enddo
      enddo

      end
