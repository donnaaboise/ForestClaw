c     # --------------------------------------------
c     # Default routines
c     #
c     # fclaw2d_fort_tag4refinement
c     # fclaw2d_fort_tag4coarsening
c     # fclaw2d_fort_interpolate2fine
c     # fclaw2d_fort_average2coarse
c     # --------------------------------------------

c     # We tag for coarsening if this coarsened patch isn't tagged for refinement
      subroutine tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      integer blockno
      double precision xlower(0:3), ylower(0:3), dx, dy
      double precision coarsen_threshold
      double precision q0(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q1(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q2(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q3(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax

c     # Assume that we will coarsen a family unless we find a grid
c     # that doesn't pass the coarsening test.
      tag_patch = 1

      mq = 1
      qmin = q0(1,1,mq)
      qmax = q0(1,1,mq)

c     # If we find that (qmax-qmin > coarsen_threshold) on any
c     # grid, we return immediately, since the family will then
c     # not be coarsened.

      call periodic_get_minmax(mx,my,mbc,meqn,mq,blockno,q0,qmin,qmax,
     &      dx,dy,xlower(0),ylower(0),coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call periodic_get_minmax(mx,my,mbc,meqn,mq,blockno, q1,qmin,qmax,
     &      dx,dy,xlower(1),ylower(1),coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call periodic_get_minmax(mx,my,mbc,meqn,mq,blockno,q2,qmin,qmax,
     &      dx,dy,xlower(2),ylower(2),coarsen_threshold,tag_patch)
      if (tag_patch == 0) return

      call periodic_get_minmax(mx,my,mbc,meqn,mq,blockno,q3,qmin,qmax,
     &      dx,dy,xlower(3),ylower(3),coarsen_threshold,tag_patch)

      end

      subroutine periodic_get_minmax(mx,my,mbc,meqn,mq,blockno, q,
     &      qmin,qmax,dx,dy,xlower,ylower, 
     &      coarsen_threshold,tag_patch)

      implicit none
      integer mx,my,mbc,meqn,mq,tag_patch, blockno
      double precision coarsen_threshold
      double precision dx,dy,xlower,ylower
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer example
      common /example_comm/ example  

      double precision qmin,qmax
      double precision qx, qy, xc,yc
      integer i,j, init_flag
      logical refine

      init_flag = 0
      refine = .false.

      call tag4refinement(mx,my,mbc,meqn,xlower,ylower,dx,dy,
     &      blockno, q, coarsen_threshold, init_flag, refine)

      if (refine) then
          tag_patch = 0
      else
          tag_patch = 1
      endif

      end
