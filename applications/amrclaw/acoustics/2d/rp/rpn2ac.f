
c
c
c     =====================================================
      subroutine clawpack46_rpn2(ixy,maxm,meqn,mwaves,mbc,mx,
     &      ql,qr,auxl,auxr,wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann solver for the acoustics equations in 2d,
c
c     # Note that although there are 3 eigenvectors, the second eigenvalue
c     # is always zero and so we only need to compute 2 waves.
c     #
c     # solve Riemann problems along one slice of data.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # This data is along a slice in the x-direction if ixy=1
c     #                            or the y-direction if ixy=2.
c     # On output, wave contains the waves,
c     #            s the speeds,
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c     # aux arrays not used in this solver.
c
      implicit double precision (a-h,o-z)
c
      dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(1-mbc:maxm+mbc, meqn)
      dimension   qr(1-mbc:maxm+mbc, meqn)
      dimension apdq(1-mbc:maxm+mbc, meqn)
      dimension amdq(1-mbc:maxm+mbc, meqn)
c
c     local arrays
c     ------------
      dimension delta(3)
c
c     # density, bulk modulus, and sound speed, and impedence of medium:
c     # (should be set in setprob.f)
      common /cparam/ rho,bulk,cc,zz

c
c
c     # set mu to point to  the component of the system that corresponds
c     # to velocity in the direction of this slice, mv to the orthogonal
c     # velocity:
c
      if (ixy.eq.1) then
          mu = 2
          mv = 3
        else
          mu = 3
          mv = 2
        endif
c
c     # note that notation for u and v reflects assumption that the
c     # Riemann problems are in the x-direction with u in the normal
c     # direciton and v in the orthogonal direcion, but with the above
c     # definitions of mu and mv the routine also works with ixy=2
c     # in which case waves come from the
c     # Riemann problems u_t + g(u)_y = 0 in the y-direction.
c
c
c     # split the jump in q at each interface into waves
c
c     # find a1 and a2, the coefficients of the 2 eigenvectors:
      do 20 i = 2-mbc, mx+mbc
         delta(1) = ql(i,1) - qr(i-1,1)
         delta(2) = ql(i,mu) - qr(i-1,mu)
         a1 = (-delta(1) + zz*delta(2)) / (2.d0*zz)
         a2 = (delta(1) + zz*delta(2)) / (2.d0*zz)
c
c        # Compute the waves.
c
         wave(i,1,1) = -a1*zz
         wave(i,mu,1) = a1
         wave(i,mv,1) = 0.d0
         s(i,1) = -cc
c
         wave(i,1,2) = a2*zz
         wave(i,mu,2) = a2
         wave(i,mv,2) = 0.d0
         s(i,2) = cc
c
   20    continue
c
c
c     # compute the leftgoing and rightgoing flux differences:
c     # Note s(i,1) < 0   and   s(i,2) > 0.
c
      do 220 m=1,meqn
         do 220 i = 2-mbc, mx+mbc
            amdq(i,m) = s(i,1)*wave(i,m,1)
            apdq(i,m) = s(i,2)*wave(i,m,2)
  220       continue
c
      return
      end
