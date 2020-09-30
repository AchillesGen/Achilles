program test_constants
    use libutilities
    use libvectors
    use libpartinfo
    use liblogging
    use libparticle
    implicit none

    type(fourvector) :: vec1
    type(threevector) :: vec2
    type(pinfo) :: info
    type(particle) :: part1

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

end program
