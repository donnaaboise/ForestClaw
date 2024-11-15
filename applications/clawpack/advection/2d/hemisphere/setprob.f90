subroutine setprob()
    implicit none

    double precision pi, pi2
    common /compi/ pi, pi2

    integer example, mapping
    common /sphere_int/ example, mapping

    double precision alpha
    common /sphere_map/ alpha

    double precision revs_per_second
    common /spherecomm/ revs_per_second

    pi = 4.d0*atan(1.d0)
    pi2 = 2*pi

    open(10,file='setprob.data')
    read(10,*) example
    read(10,*) mapping
    read(10,*) alpha
    read(10,*) revs_per_second

    close(10)

end
