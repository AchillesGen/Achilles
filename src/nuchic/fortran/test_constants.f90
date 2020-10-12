program test_constants
    use libutilities
    use libvectors
    use libpartinfo
    use liblogging
    use libparticle
    use libinterpolate
    implicit none

    type(fourvector) :: vec1
    type(threevector) :: vec2
    type(pinfo) :: info
    type(particle) :: part1
    type(interp1d) :: interp
    type(interp2d) :: interp2
    double precision, dimension(10), parameter :: x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    double precision, dimension(10), parameter :: y = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    double precision, dimension(3), parameter :: x1 = [1, 2, 3]
    double precision, dimension(3), parameter :: x2 = [1, 2, 3]
    double precision, dimension(9), parameter :: z = [11, 12, 13, 21, 22, 23, 31, 32, 33]

    print*, constants%c, constants%hbarc, constants%hbarc2, constants%mp, constants%mn, constants%mqe
    call init(constants)
    print*, constants%c, constants%hbarc, constants%hbarc2, constants%mp, constants%mn, constants%mqe

    vec1 = fourvector(1d0, 2d0, 3d0, 4d0)
    vec1 = vec1+vec1
    vec1 = 3d0*vec1 + vec1*2d0 - vec1/2d0
    vec2 = threevector(1d0, 2d0, 3d0)
    call vec1%print()
    call vec2%print()

    print*, vec1*vec1 

    info = pinfo(2212)
    call logger%info(info%name())
    call logger%warn("test")

    part1 = particle(2212, vec1, vec2, 0)

    vec1 = part1%momentum()
    call vec1%print()

    call part1%set_momentum(fourvector(1d0, 2d0, 3d0, 4d0))
    vec1 = part1%momentum()
    call vec1%print()

    interp = interp1d(x, y, 10)

    print*, interp%min()
    print*, interp%max()
    print*, interp%call(3.3d0)

    interp2 = interp2d(x, y, z, 3, 3)

    print*, interp2%call(2.3d0, 2.3d0)

end program
