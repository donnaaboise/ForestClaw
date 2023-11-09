subroutine clawpack46_rpt2(ixy,maxm,meqn,mwaves,mbc,mx, &
    ql,qr,aux1,aux2,aux3,imp,asdq,bmasdq,bpasdq)

    implicit none

    !! # Riemann solver in the transverse direction for the shallow water
    !! # equations  on a quadrilateral grid.

    !! # Split asdq (= A^* \Delta q, where * = + or -)
    !! # into down-going flux difference bmasdq (= B^- A^* \Delta q)
    !! #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)

    !! # Uses Roe averages and other quantities which were
    !! # computed in rpn2sh and stored in the common block comroe.

    integer ixy,maxm,meqn,mwaves,mbc,mx, imp

    double precision     ql(1-mbc:maxm+mbc, meqn)
    double precision     qr(1-mbc:maxm+mbc, meqn)
    double precision   asdq(1-mbc:maxm+mbc, meqn)
    double precision bmasdq(1-mbc:maxm+mbc, meqn)
    double precision bpasdq(1-mbc:maxm+mbc, meqn)
    double precision   aux1(1-mbc:maxm+mbc, *)
    double precision   aux2(1-mbc:maxm+mbc, *)
    double precision   aux3(1-mbc:maxm+mbc, *)

    double precision u(-1:maxm),v(-1:maxm),a(-1:maxm),h(-1:maxm)
    double precision wave(-1:maxm, 4, 3)
    double precision    s(-1:maxm, 3)
    double precision enx(-1:maxm), eny(-1:maxm), enz(-1:maxm)
    double precision etx(-1:maxm), ety(-1:maxm), etz(-1:maxm)
    double precision gamma(-1:maxm)

    double precision grav
    common /swe_parms/  grav

    double precision dtcom, dxcom, dycom, tcom
    integer icom,jcom
    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

    double precision delta(4)
    double precision a1, a2, a3, dx, dy, bn
    double precision erx, ery, erz 
    integer i, j, m, i1, mw, ioff, ix1, ixm1

    if (ixy .eq. 1) then
        dx = dxcom
    else
        dx = dycom
    endif


    if (ixy .eq. 1) then
        ioff = 7
    else
        ioff = 1
    endif


    !! # imp is used to flag whether wave is going to left or right,
    !! # since states and grid are different on each side
    if (imp .eq. 1) then
        !! # asdq = amdq, moving to left
        ix1 = 2-mbc
        ixm1 = mx+mbc
    else
        !! # asdq = apdq, moving to right
        ix1 = 1-mbc
        ixm1 = mx+mbc
    endif

    !!  --------------
    !!  # up-going:
    !!  --------------

    !! # determine rotation matrix for interface above cell, using
    !! # aux3

    !! [ alf  beta ]
    !! [-beta  alf ]
    
    do i = ix1,ixm1

        if (imp.eq.1) then
            i1 = i-1
        else
            i1 = i
        endif

        enx(i) =   aux3(i1,ioff+1)
        eny(i) =   aux3(i1,ioff+2)
        enz(i) =   aux3(i1,ioff+3)

        etx(i) =   aux3(i1,ioff+4)
        ety(i) =   aux3(i1,ioff+5)
        etz(i) =   aux3(i1,ioff+6)

        gamma(i) = dsqrt(etx(i)**2 + ety(i)**2 + etz(i)**2)
        etx(i) =   etx(i) / gamma(i)
        ety(i) =   ety(i) / gamma(i)
        etz(i) =   etz(i) / gamma(i)

        h(i) = ql(i1,1)
        u(i) = (enx(i)*ql(i1,2)+eny(i)*ql(i1,3)+enz(i)*ql(i1,4)) & 
              / h(i)
        v(i) = (etx(i)*ql(i1,2)+ety(i)*ql(i1,3)+etz(i)*ql(i1,4)) &
              / h(i)
        a(i) = dsqrt(grav*h(i))
    end do



    !! # now split asdq into waves:

    !! # find a1 thru a3, the coefficients of the 3 eigenvectors:
    do  i = ix1,ixm1
        delta(1) = asdq(i,1)
        delta(2) = enx(i)*asdq(i,2)+eny(i)*asdq(i,3)+enz(i)*asdq(i,4)
        delta(3) = etx(i)*asdq(i,2)+ety(i)*asdq(i,3)+etz(i)*asdq(i,4)
        a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
        a2 = -v(i)*delta(1) + delta(3)
        a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))

        !! # Compute the waves.

        wave(i,1,1) = a1
        wave(i,2,1) = enx(i)*a1*(u(i)-a(i)) + etx(i)*a1*v(i)
        wave(i,3,1) = eny(i)*a1*(u(i)-a(i)) + ety(i)*a1*v(i)
        wave(i,4,1) = enz(i)*a1*(u(i)-a(i)) + etz(i)*a1*v(i)
        s(i,1) = (u(i)-a(i))*gamma(i)/dx

        wave(i,1,2) = 0.0d0
        wave(i,2,2) = etx(i)*a2
        wave(i,3,2) = ety(i)*a2
        wave(i,4,2) = etz(i)*a2
        s(i,2) = u(i) * gamma(i)/dx

        wave(i,1,3) = a3
        wave(i,2,3) = enx(i)*a3*(u(i)+a(i)) + etx(i)*a3*v(i)
        wave(i,3,3) = eny(i)*a3*(u(i)+a(i)) + ety(i)*a3*v(i)
        wave(i,4,3) = enz(i)*a3*(u(i)+a(i)) + etz(i)*a3*v(i)
        s(i,3) = (u(i)+a(i)) * gamma(i)/dx
    end do

    !! --------------------------------
    !! # compute flux difference bpasdq
    !! --------------------------------

    do m = 1,meqn
        do i = ix1,ixm1
            bpasdq(i,m) = 0.d0
            do mw = 1,mwaves
                bpasdq(i,m) = bpasdq(i,m) & 
                    + dmax1(s(i,mw),0.d0)*wave(i,m,mw)
            end do
        end do
    end do


    !! # project momentum component of bpasdq to tangent plane:
    do i = ix1, ixm1
        if (imp.eq.1) then
            i1 = i-1
        else
            i1 = i
        endif
        erx = aux3(i1,14)
        ery = aux3(i1,15)
        erz = aux3(i1,16)
        bn = erx*bpasdq(i,2)+ery*bpasdq(i,3)+erz*bpasdq(i,4)
        bpasdq(i,2) = bpasdq(i,2) - bn*erx
        bpasdq(i,3) = bpasdq(i,3) - bn*ery
        bpasdq(i,4) = bpasdq(i,4) - bn*erz
    end do

    !! --------------
    !! # down-going:
    !! --------------


    !! # determine rotation matrix for interface below cell, using
    !! aux2
    !! [ alf  beta ]
    !! [-beta  alf ]

    do i=ix1,ixm1

        if (imp.eq.1) then
            i1 = i-1
        else
            i1 = i
        endif

        enx(i) =   aux2(i1,ioff+1)
        eny(i) =   aux2(i1,ioff+2)
        enz(i) =   aux2(i1,ioff+3)
        etx(i) =   aux2(i1,ioff+4)
        ety(i) =   aux2(i1,ioff+5)
        etz(i) =   aux2(i1,ioff+6)
        gamma(i) = dsqrt(etx(i)**2 + ety(i)**2 + etz(i)**2)
        etx(i) =   etx(i) / gamma(i)
        ety(i) =   ety(i) / gamma(i)
        etz(i) =   etz(i) / gamma(i)
        u(i) = (enx(i)*ql(i1,2)+eny(i)*ql(i1,3)+enz(i)*ql(i1,4)) & 
                 / h(i)
        v(i) = (etx(i)*ql(i1,2)+ety(i)*ql(i1,3)+etz(i)*ql(i1,4)) &
                 / h(i)
    end do



    !!# now split asdq into waves:

    !!# find a1 thru a3, the coefficients of the 3 eigenvectors:
    do i = ix1,ixm1
        delta(1) = asdq(i,1)
        delta(2) = enx(i)*asdq(i,2)+eny(i)*asdq(i,3)+enz(i)*asdq(i,4)
        delta(3) = etx(i)*asdq(i,2)+ety(i)*asdq(i,3)+etz(i)*asdq(i,4)

        a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
        a2 = -v(i)*delta(1) + delta(3)
        a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))

        !! # Compute the waves.

        wave(i,1,1) = a1
        wave(i,2,1) = enx(i)*a1*(u(i)-a(i)) + etx(i)*a1*v(i)
        wave(i,3,1) = eny(i)*a1*(u(i)-a(i)) + ety(i)*a1*v(i)
        wave(i,4,1) = enz(i)*a1*(u(i)-a(i)) + etz(i)*a1*v(i)
        s(i,1) = (u(i)-a(i)) * gamma(i)/dx

        wave(i,1,2) = 0.0
        wave(i,2,2) = etx(i)*a2
        wave(i,3,2) = ety(i)*a2
        wave(i,4,2) = etz(i)*a2
        s(i,2) = u(i) * gamma(i)/dx

        wave(i,1,3) = a3
        wave(i,2,3) = enx(i)*a3*(u(i)+a(i)) + etx(i)*a3*v(i)
        wave(i,3,3) = eny(i)*a3*(u(i)+a(i)) + ety(i)*a3*v(i)
        wave(i,4,3) = enz(i)*a3*(u(i)+a(i)) + etz(i)*a3*v(i)
        s(i,3) = (u(i)+a(i)) * gamma(i)/dx
    end do


    !! # compute flux difference bmasdq
    !! --------------------------------

    do  m=1,meqn
        do  i = ix1,ixm1
            bmasdq(i,m) = 0.d0
            do  mw = 1,mwaves
                bmasdq(i,m) = bmasdq(i,m) &
                    + dmin1(s(i,mw), 0.d0)*wave(i,m,mw)
            end do
        end do
    end do

    !! # project momentum component of bmasdq to tangent plane:
    do i = ix1, ixm1
        if (imp.eq.1) then
            i1 = i-1
        else
            i1 = i
        endif
        erx = aux1(i1,14)
        ery = aux1(i1,15)
        erz = aux1(i1,16)
        bn = erx*bmasdq(i,2)+ery*bmasdq(i,3)+erz*bmasdq(i,4)
        bmasdq(i,2) = bmasdq(i,2) - bn*erx
        bmasdq(i,3) = bmasdq(i,3) - bn*ery
        bmasdq(i,4) = bmasdq(i,4) - bn*erz
    end do

    return
end
