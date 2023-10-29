subroutine magic_div_tau_mapped(blockno, mx,my,mz, mbc, meqn, & 
    maux, dx,dy, dz, dt, xp,yp,zp, xd,yd,zd, & 
    volumes, xrot,yrot,zrot, faceareas, &
    aux, q)

    implicit none

    INTEGER mx,my,mz,mbc, meqn, maux, blockno
    double precision dt

    !! delta on computational grid
    double precision dx, dy, dz 

    double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    !! Stores mu, lambda somewhere
    double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+1,maux)

    double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)

    double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)
    double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)
    double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)

    double precision volumes(-mbc:mx+mbc+1,-mbc:my+mbc+1, -mbc:mz+mbc+1)
    double precision faceareas(-mbc:mx+mbc+1,-mbc:my+mbc+1, -mbc:mz+mbc+2,3)

    double precision xrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)
    double precision yrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)
    double precision zrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)

    double precision g1(3), g2(3), g3(3), g1_up(3), g2_up(3), g3_up(3)

    !! Dummy arrays
    double precision up(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+1)
    double precision vp(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+1)
    double precision wp(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+1)

    double precision ud(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+1)
    double precision vd(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+1)
    double precision wd(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+1)

    double precision tau_n(1:mx+1,1:my+1,1:mz+1,3,3)
    double precision dUdx(3), dUdy(3), dUdz(3)

    double precision gxn(3), gradxn(3), gradxTn(3)
    double precision nvec(3), hex(0:1,0:1,0:1,3), quad(0:1,0:1,3)

    INTEGER i,j,k, m, dim, ii, jj, kk
    double precision div_tau, divu, lmbda, mu
    double precision du_xi(3), du_eta(3)    

    lmbda = 1
    mu = 1

    write(6,*) 'Calling div_tau'


    !! Store Cartesian components (u,v) extracted from 
    !! state vectors q. 
    do k = 1-mbc,mz+mbc
        do i = 1-mbc,mx+mbc
            do j = 1-mbc,my+mbc
                up(i,j,k) = q(i,j,k,2)/q(i,j,k,1)
                vp(i,j,k) = q(i,j,k,3)/q(i,j,k,1)
                wp(i,j,k) = q(i,j,k,4)/q(i,j,k,1)
            end do
        end do
    end do

    do k = 0,mz+1
        do i = 0,mx+1
            do j = 0,my+1
                ud(i,j,k) = (up(i-1,j,k)   + up(i,j,k)   + up(i,j-1,k)   + up(i-1,j-1,k) + & 
                             up(i-1,j,k-1) + up(i,j,k-1) + up(i,j-1,k-1) + up(i-1,j-1,k-1))/8.
                vd(i,j,k) = (vp(i-1,j,k)   + vp(i,j,k)   + vp(i,j-1,k)   + vp(i-1,j-1,k) + & 
                             vp(i-1,j,k-1) + vp(i,j,k-1) + vp(i,j-1,k-1) + vp(i-1,j-1,k-1))/8.
                wd(i,j,k) = (wp(i-1,j,k)   + wp(i,j,k)   + wp(i,j-1,k)   + wp(i-1,j-1,k) + & 
                             wp(i-1,j,k-1) + wp(i,j,k-1) + wp(i,j-1,k-1) + wp(i-1,j-1,k-1))/8.
            end do
        end do
    end do


    !! Step 1 : Compute grad v at x-faces

    do dim = 1,3
        do k = 1,mz+1
            do j = 1,my+1
                do i = 1,mx+1
                    do ii = 0,1
                        do jj = 0,1
                            do kk = 0,1
                                hex(ii,jj,kk,1) = xd(i+ii,j+jj,k+kk)
                                hex(ii,jj,kk,2) = yd(i+ii,j+jj,k+kk)
                                hex(ii,jj,kk,3) = zd(i+ii,j+jj,k+kk)
                            end do
                        end do
                    end do

                    !! We need to include the harmonic average to get the correct
                    !! corvariant vectors in the non-tangent vectors. 
                    call magic_hex_compute_vectors(hex,g1,g2,g3, g1_up, g2_up, g3_up)

                    !! Compute covariant vectors
                    if (dim == 1) then 
                        !! x-face (left)
                        dUdx(1) = (up(i,j,k) - up(i-1,j,k))/dx
                        dUdx(2) = (vp(i,j,k) - vp(i-1,j,k))/dx
                        dUdx(3) = (wp(i,j,k) - wp(i-1,j,k))/dx

                        !! Compute other two vectors
                        do kk = 0,1
                            do jj = 0,1
                                quad(jj,kk,1) = ud(i,j+jj,k+kk)
                                quad(jj,kk,2) = vd(i,j+jj,k+kk)
                                quad(jj,kk,3) = wd(i,j+jj,k+kk)
                            end do
                        end do
                        call magic_u_derivatives(quad, du_xi, du_eta)
                        do ii = 1,3
                            dUdy(ii) = du_xi(ii)
                            dUdz(ii) = du_eta(ii)
                        enddo 

                        do ii = 1,3
                            nvec(ii) = xrot(i,j,k,1,ii)
                        end do

                    elseif (dim == 2) then
                        !! y-face (bottom)
                        dUdy(1) = (up(i,j,k) - up(i,j-1,k))/dy
                        dUdy(2) = (vp(i,j,k) - vp(i,j-1,k))/dy
                        dUdy(3) = (wp(i,j,k) - wp(i,j-1,k))/dy

                        !! Compute other two vectors
                        do ii = 0,1
                            do kk = 0,1
                                quad(ii,kk,1) = ud(i+ii,j,k+kk)
                                quad(ii,kk,2) = vd(i+ii,j,k+kk)
                                quad(ii,kk,3) = wd(i+ii,j,k+kk)
                            end do
                        end do
                        call magic_u_derivatives(quad, du_xi, du_eta)
                        do ii = 1,3
                            dUdz(ii) = du_xi(ii)
                            dUdx(ii) = du_eta(ii)
                        enddo 

                        do ii = 1,3
                            nvec(ii) = yrot(i,j,k,1,ii)
                        end do
                    elseif (dim == 3) then
                        !! y-face (bottom)
                        dUdz(1) = (up(i,j,k) - up(i,j,k-1))/dz
                        dUdz(2) = (vp(i,j,k) - vp(i,j,k-1))/dz
                        dUdz(3) = (wp(i,j,k) - wp(i,j,k-1))/dz

                        !! Compute other two vectors
                        do ii = 0,1
                            do jj = 0,1
                                quad(ii,jj,1) = ud(i+ii,j+jj,k)
                                quad(ii,jj,2) = vd(i+ii,j+jj,k)
                                quad(ii,jj,3) = wd(i+ii,j+jj,k)
                            end do
                        end do
                        call magic_u_derivatives(quad, du_xi, du_eta)
                        do ii = 1,3
                            dUdx(ii) = du_xi(ii)
                            dUdy(ii) = du_eta(ii)
                        end do 

                        do ii = 1,3
                            nvec(ii) = zrot(i,j,k,1,ii)
                        end do
                    endif

                    !! ------------------------
                    !! Compute (grad x) n
                    !! ------------------------

                    !! g1_up dot n
                    gxn(1) = g1_up(1)*nvec(1) + g1_up(2)*nvec(2) + g1_up(3)*nvec(3)

                    !! g2_up dot n
                    gxn(2) = g2_up(1)*nvec(1) + g2_up(2)*nvec(2) + g2_up(3)*nvec(3)

                    !! g3_up dot n
                    gxn(3) = g3_up(1)*nvec(1) + g3_up(2)*nvec(2) + g3_up(3)*nvec(3)

                    !! Compute (grad x) normal
                    do ii = 1,3
                        gradxn(ii) = gxn(1)*dUdx(ii) + gxn(2)*dUdy(ii) + gxn(3)*dUdz(ii)
                    end do

                    !! ------------------------
                    !! Compute (grad x)^T n
                    !! ------------------------
                    gxn(1) = dUdx(1)*nvec(1) + dUdx(2)*nvec(2) + dUdx(3)*nvec(3)
                    gxn(2) = dUdy(1)*nvec(1) + dUdy(2)*nvec(2) + dUdy(3)*nvec(3)
                    gxn(3) = dUdz(1)*nvec(1) + dUdz(2)*nvec(2) + dUdz(3)*nvec(3)

                    !! Compute (grad x)^T normal
                    do ii = 1,3
                        gradxTn(ii) = gxn(1)*g1_up(ii) + gxn(2)*g2_up(ii) + gxn(3)*g3_up(ii)
                    end do

                    !! -----------------------------------------------------
                    !! compute div(u) = tr(D) when D = 1/2*(gradx + gradx^T)
                    !! -----------------------------------------------------

                    divu = (dUdx(1)*g1_up(1) + dUdx(2)*g1_up(2) + dUdx(3)*g1_up(3)) + & 
                        (dUdy(1)*g2_up(1) + dUdy(2)*g2_up(2) + dUdy(3)*g2_up(3)) + &
                        (dUdz(1)*g3_up(1) + dUdz(2)*g3_up(3) + dUdz(3)*g3_up(3))

                    do ii = 1,3
                        tau_n(i,j,k,dim,ii) = mu*(gradxn(ii) + gradxTn(ii)) + & 
                            lmbda*divu*nvec(ii)
                    enddo
                end do
            end do
        end do
    end do

    do k = 1,mz
        do i = 1,mx
            do j = 1,my

                do ii = 1,3
                    div_tau = faceareas(i+1,j,k,1) *tau_n(i+1,j,k,1,ii) -  &
                          faceareas(i,j,k,1)   *tau_n(i,j,k,1,ii) + &                          
                          faceareas(i,j+1,k,2) *tau_n(i,j+1,k,2,ii) - &
                          faceareas(i,j,k,2)   *tau_n(i,j,k,2,ii)  +    &                  
                          faceareas(i,j,k+1,3) *tau_n(i,j,k+1,3,ii) - &
                          faceareas(i,j,k,3)  *tau_n(i,j,k,3,ii)
                    q(i,j,k,1+ii) = q(i,j,k,1+ii) + dt*div_tau/volumes(i,j,k)
                end do
            end do
        end do
    end do



end

!! Won't get called
subroutine magic_compute_basis(mx, my, mz, mbc, xd, yd, zd, & 
          xrot, yrot, zrot)
    implicit none

    integer mx,my,mz, mbc

    double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)
    double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)
    double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)

    double precision xrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)
    double precision yrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)
    double precision zrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)

    integer i,j, k, ii,jj, kk
    !! double precision xcorner, ycorner, zcorner
    !! double precision xe, ye,ze, xp1, yp1, zp1
    double precision hex(0:1,0:1,0:1,3)

    double precision rot(3,3,3)
    double precision g1(3), g2(3), g3(3), g1_up(3), g2_up(3), g3_up(3)

!!    integer*8 map_context_ptr, fclaw_map_get_context
!!
!!    map_context_ptr = fclaw_map_get_context()

    do j = -mbc,my+mbc+1
        do i = -mbc,mx+mbc+1
            do k = -mbc,mz+mbc+1
                do ii = 0,1
                    do jj = 0,1
                        do kk = 0,1
                            hex(ii,jj,kk,1) = xd(i+ii,j+jj,k+kk)
                            hex(ii,jj,kk,2) = yd(i+ii,j+jj,k+kk)
                            hex(ii,jj,kk,3) = zd(i+ii,j+jj,k+kk)
                        end do
                    end do
                end do
                call magic_hex_compute_vectors(hex,g1,g2,g3,g1_up, g2_up, g3_up)
            end do
        end do
    end do
end subroutine magic_compute_basis

subroutine magic_hex_compute_vectors(hex,g1,g2,g3,g1_up, g2_up, g3_up)
    implicit none

    double precision hex(0:1,0:1,0:1,3)
    double precision g1(3), g2(3), g3(3), g1_up(3), g2_up(3), g3_up(3)

    double precision Jb(3,3), Jinv(3,3), xcv(3), ycv(3), zcv(3)

    double precision z000(3),z100(3),z010(3),z001(3), z110(3), & 
           z101(3), z011(3), z111(3)

    double precision a000(3),a100(3),a010(3),a001(3), a110(3), & 
           a101(3), a011(3), a111(3)

    logical is_id
    integer i,j,k, n, ni, info
    double precision rhs(3,3)


!!    do i = 1,3
!!        do j = 1,3
!!            if (i == j) then
!!                eye(i,j) = 1.d0
!!            else
!!                eye(i,j) = 0.d0
!!            endif
!!        end do
!!    end do


    !! Get centers of faces so we can find basis vectors at each face.
    do i = 1,3
        xcv(i) = 0.5d0
        ycv(i) = 0.5d0
        zcv(i) = 0.5d0
    end do

    xcv(1) = 0.d0
    ycv(2) = 0.d0
    zcv(3) = 0.d0

    !! # Make notation easy..
    do i = 1,3
        z000(i) = hex(0,0,0,i)
        z100(i) = hex(1,0,0,i)
        z010(i) = hex(0,1,0,i)
        z001(i) = hex(0,0,1,i)
        z110(i) = hex(1,1,0,i)
        z101(i) = hex(1,0,1,i)
        z011(i) = hex(0,1,1,i)
        z111(i) = hex(1,1,1,i)
    end do

    !! # Get coefficients do trilinear map.
    do i = 1,3
        a000(i) = z000(i)
        a001(i) = z001(i) - z000(i)
        a010(i) = z010(i) - z000(i)
        a100(i) = z100(i) - z000(i)

        a111(i) = z111(i) - z110(i) - z101(i) - & 
              z011(i) + z100(i) + z010(i) + z001(i) - z000(i)
        a011(i) = z011(i) - z010(i) - z001(i) + z000(i)
        a101(i) = z101(i) - z100(i) - z001(i) + z000(i)
        a110(i) = z110(i) - z100(i) - z010(i) + z000(i)
    end do

    !! # Start computing basis vectors.
    do n = 1,3  !! Loop over each face.

        do i = 1,3 !! Get three rows of Jacobian
            Jb(i,1) = a100(i) + a110(i)*ycv(n)  + a101(i)*zcv(n) + & 
                 a111(i)*ycv(n)*zcv(n)
            Jb(i,2) = a010(i) + a110(i)*xcv(n)  + a011(i)*zcv(n) + & 
                 a111(i)*xcv(n)*zcv(n)
            Jb(i,3) = a001(i) + a101(i)*xcv(n)  + a011(i)*ycv(n) + & 
                 a111(i)*xcv(n)*ycv(n)
        end do        

        do i = 1,3
            if (n .eq. 1) then
                g1(i) = Jb(i,1)
            else if (n .eq. 2) then
                g2(i) = Jb(i,2)
            else if (n .eq. 3) then
                g3(i) = Jb(i,3)
            endif 
        end do


        !! # Compute inverse of Jacobian.
        !!  call dgesv(3,3,Jb,3,IPIV,rhs,3,info)
        call magic_hex_compute_Jinv(Jb,rhs,info)

        do i = 1,3
            if (n .eq. 1) then
                g1_up(i) = rhs(1,1)*g1(i) + rhs(1,2)*g2(i) + rhs(1,3)*g3(i)
            else if (n .eq. 2) then
                g2_up(i) = rhs(2,1)*g1(i) + rhs(2,2)*g2(i) + rhs(2,3)*g3(i)
            else if (n .eq. 3) then
                g3_up(i) = rhs(3,1)*g1(i) + rhs(3,2)*g2(i) + rhs(3,3)*g3(i)
            endif 
        end do
    end do

end subroutine magic_hex_compute_vectors


!! # This computes the inverse of a 3x3 matrix using Cramer's rule
!! # Note : Only if det(Jac) == 0 (exactly) will it report that
!! # the matrix is singular (i.e. info == 1).   It will not detect
!! # potential ill-conditioning.
subroutine magic_hex_compute_Jinv(Jac,Jinv,info)
    implicit none

    double precision Jac(3,3), Jinv(3,3)
    double precision hex_dot_cross, detJ, s
    integer i,j, k(3,2), info

    data k /2, 1, 1, 3, 3, 2/

    info = 0

    !! # Compute determinant of Jacobian
    detJ  = hex_dot_cross(jac(1,1),jac(1,2),jac(1,3))
    if (detJ .eq. 0) then
        info = 1
        return
    endif

    !! # Apply Cramer's rule to get inverse.
    s = 1.d0/detJ
    do j = 1,3
        do i = 1,3
            Jinv(i,j) = s*(Jac(k(j,1),k(i,1))*Jac(k(j,2),k(i,2)) - & 
                           Jac(k(j,1),k(i,2))*Jac(k(j,2),k(i,1)))
            s = -s
        end do
    end do
end subroutine magic_hex_compute_Jinv


double precision function magic_hex_dot_cross(u,v,w)
    implicit none
    double precision u(3), v(3), w(3)

    magic_hex_dot_cross = u(1)*(v(2)*w(3) - v(3)*w(2)) - & 
                u(2)*(v(1)*w(3) - v(3)*w(1)) + & 
                u(3)*(v(1)*w(2) - v(2)*w(1))

    return
end function magic_hex_dot_cross


subroutine magic_u_derivatives(quad, du_xi, du_eta)
    implicit none

    double precision quad(0:1,0:1,3), du_xi(3), du_eta(3)

    double precision a00(3), a01(3),a10(3), a11(3)
    double precision xi,eta

    integer m

    !! # Coefficients do bilinear approximation to surface
    do m = 1,3
        a00(m) = quad(0,0,m)
        a01(m) = quad(1,0,m) - quad(0,0,m)
        a10(m) = quad(0,1,m) - quad(0,0,m)
        a11(m) = quad(1,1,m) - quad(1,0,m) - & 
              quad(0,1,m) + quad(0,0,m)
    end do

    eta = 0.5d0
    xi  = 0.5d0

    !! # Mesh square is approximated by
    !! #       U(xi,eta) = a00 + a01*xi + a10*eta + a11*xi*eta
    !! # Differentiate U to get dU/dxi and dU/deta
    do m = 1,3
        du_xi(m)  = (a01(m) + a11(m)*eta)
        du_eta(m) = (a10(m) + a11(m)*xi)
    end do

end subroutine magic_u_derivatives




