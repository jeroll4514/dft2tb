! A = lattice vectors in each column
! B = reciprocal lattice vectors in each column
! dists = radii between atoms in supercell
! distsU = unique interactions (radius, atomic species, zeta)
! eH/eS = stores tight-binding recreation of sparse elements
! has_s, p, d = logical variable if basis set includes s, p, d orbitals
! hsparse = sparse hamiltonian matrix
! indx = finds atom and orbital for each row/column
! indxuo = number of equivalent orbitals in unit cell
! ind_store = stores indices for given interactions
! line = used to search for lattice vectors within input.fdf
! lmn = normalized version of dists (for SK)
! L, M, N = components of vector from lmn
! no_u = number of orbitals in the unit cell
! no_s = number of orbitals in the super cell
! numh = number of nonzero elements per row of H
! numInters = number of unique interactions
! nSpin = number of spin components
! numDim = dimensionality of sparse matrices
! orb = index which corresponds to possible orbitals
! phase = exponential component for matrix elements (folding)
! reord = lists indices for py, pz, px -> px, py, pz for more concise coding
! signs = used to take into account opposite pz orientation in SIESTA
! soe = stores systems of equations to determine pp and dd parameters
! ssparse = sparse overlap matrix
! xij = stores locations of atoms within supercell

! ii, jj, kk, LL, t, t2 = dummy variables for loops and index storing
! tmp = dummy variable for reading orbital indices

PROGRAM MAIN
    IMPLICIT NONE

    INTEGER ii, jj, kk, LL, mm, no_u, no_s, nSpin, numDim, t, t2, t3, numInters
    INTEGER reord(4), signs(4)
    REAL distsU(1000, 5), soe(9), A(3, 3), B(3, 3), vs(3), ve(3)
    REAL L, M, N, eH, eS, s3, pi, temp1, temp2
    COMPLEX e, phase
    CHARACTER(LEN=2) tmp
    CHARACTER(LEN=3) file_status
    CHARACTER(LEN=1000) line
    LOGICAL tc, has_s, has_p, has_d, file_exists

    INTEGER, ALLOCATABLE :: indxuo(:), numh(:), listh(:), listhptr(:)
    INTEGER, ALLOCATABLE :: indx(:, :), ind_store(:, :), kden(:), xlabInd(:)
    REAL, ALLOCATABLE :: hsparse(:), ssparse(:), xij(:, :), dists(:), lmn(:, :)
    REAL, ALLOCATABLE :: kpE(:, :), kp(:, :), mags(:), xout(:), Vonsite(:, :, :)
    REAL, ALLOCATABLE :: Vss(:, :), Vsp(:, :), Vps(:, :), VppSig(:, :)
    REAL, ALLOCATABLE :: VppPi(:, :), Vsd(:, :), Vds(:, :), VpdSig(:, :)
    REAL, ALLOCATABLE :: VpdPi(:, :), VdpSig(:, :), VdpPi(:, :), VddSig(:, :)
    REAL, ALLOCATABLE :: VddPi(:, :), VddDelt(:, :), energy(:), rwork(:)
    COMPLEX, ALLOCATABLE :: H(:, :), S(:, :), work(:)
    CHARACTER, ALLOCATABLE :: xlab(:), xlabF(:, :)

    !------------------------------------------------------!
    !---------------- Reads all input files ---------------!
    !------------------------------------------------------!

    !-------------First reads data from HSX file------------!
    inquire(file='HSX_out', exist=file_exists) ! checks if exists

    if (file_exists .eqv. .FALSE.) then
        write(*, '(a38)') 'Error: HSX_out file does not exist'
        call exit(0)
    end if

    open(unit=1, file='HSX_out', status='old') ! opens file
    read(1, *) no_u, no_s, nSpin, numDim
    read(1, *) ! looks at if gamma only calculation (not our case)

    ! reads indxuo
    allocate(indxuo(no_s))
    do ii = 0, no_s/16-1 ! through last line for indxuo
        read(1, *) indxuo(16*ii+1:16*ii+16)
    enddo ! ii
    if (16*(no_s/16) /= no_s) then
        read(1, *) indxuo(16*(no_s/16)+1:no_s)
    endif

    ! reads numh
    allocate(numh(no_u))
    do ii = 0, no_u/4-1 ! through last line for numh
        read(1, *) numh(4*ii+1:4*ii+4)
    enddo
    if (4*(no_u/4) /= no_u) then
        read(1, *) numh(4*(no_u/4)+1:no_u)
    endif

    ! makes pointers for sparse matrices
    allocate(listhptr(no_u))
    do ii=1, no_u
        listhptr(ii) = sum(numh(1:ii-1))
    enddo !ii

    ! reads listh
    allocate(listh(numDim))
    do ii = 1, no_u
        t = listhptr(ii)
        do jj = 0, numh(ii)/10-1
            read(1, *) listh(t+10*jj+1:t+10*jj+10)
        enddo ! jj
        if (10*(numh(ii)/10) /= numh(ii)) then
            read(1, *) listh(t+10*(numh(ii)/10)+1:t+numh(ii))
        endif
    enddo ! ii

    ! reads sparse hamiltonian matrix
    allocate(hsparse(numDim))
    read(1, *)
    do ii = 1, no_u
        t = listhptr(ii)
        do jj = 0, numh(ii)/10-1
            read(1, *) hsparse(t+10*jj+1:t+10*jj+10)
        enddo ! jj
        if (10*(numh(ii)/10) /= numh(ii)) then
            read(1, *) hsparse(t+10*(numh(ii)/10)+1:t+numh(ii))
        endif
    enddo ! ii
    hsparse = hsparse * 13.605698066 ! Ry -> eV

    ! reads sparse ovarlap matrix
    allocate(ssparse(numDim))
    read(1, *)
    do ii = 1, no_u
        t = listhptr(ii)
        do jj = 0, numh(ii)/10-1
            read(1, *) ssparse(t+10*jj+1:t+10*jj+10)
        enddo ! jj
        if (10*(numh(ii)/10) /= numh(ii)) then
            read(1, *) ssparse(t+10*(numh(ii)/10)+1:t+numh(ii))
        endif
    enddo ! ii

    !!reads TE overlap-x matrix
    !read(1, *)
    !do ii = 1, no_u
    !    t = listhptr(ii)
    !    do jj = 0, numh(ii)/10-1
    !        read(1, *) !ssparse(t+10*jj+1:t+10*jj+10)
    !    enddo !jj
    !    if (10*(numh(ii)/10) /= numh(ii)) then
    !        read(1, *) !ssparse(t+10*(numh(ii)/10)+1:t+numh(ii))
    !    endif
    !enddo !ii

    !!reads TE overlap-y matrix
    !read(1, *)
    !do ii = 1, no_u
    !    t = listhptr(ii)
    !    do jj = 0, numh(ii)/10-1
    !        read(1, *) !ssparse(t+10*jj+1:t+10*jj+10)
    !    enddo !jj
    !    if (10*(numh(ii)/10) /= numh(ii)) then
    !        read(1, *) !ssparse(t+10*(numh(ii)/10)+1:t+numh(ii))
    !    endif
    !enddo !ii

    !!reads TE overlap-z matrix
    !read(1, *)
    !do ii = 1, no_u
    !    t = listhptr(ii)
    !    do jj = 0, numh(ii)/10-1
    !        read(1, *) !ssparse(t+10*jj+1:t+10*jj+10)
    !    enddo !jj
    !    if (10*(numh(ii)/10) /= numh(ii)) then
    !        read(1, *) !ssparse(t+10*(numh(ii)/10)+1:t+numh(ii))
    !    endif
    !enddo !ii

    ! reads xij
    allocate(xij(numDim, 3))
    read(1, *)
    do ii = 1, numDim
        read(1, *) xij(ii, :)
    enddo ! ii
    xij = xij * 0.529177249 ! Bohr -> ang
    close(1)

    !-------------Reads data from ORB_INDX file-------------!
    ! This sorts orbitals for rows/cols
    has_s = .FALSE. ! initialize
    has_p = .FALSE. ! initialize
    has_d = .FALSE. ! initialize
    allocate(indx(no_u, 3))

    inquire(file='SnS.ORB_INDX', exist=file_exists) ! checks if exists

    if (file_exists .eqv. .FALSE.) then
        write(*, '(a38)') 'Error: .ORB_INDX file does not exist'
        call exit(0)
    end if
    open(unit=2, file='SnS.ORB_INDX', status='old') ! opens file
    read(2, *)
    read(2, *)
    read(2, *)
    do ii = 1, no_u ! cycles through all rows/columns
        read(2, *) t, t, indx(ii, 1), tmp, t, t, jj, kk, indx(ii, 3)
        select case(jj)
        case(0) ! s orbital
            has_s = .TRUE.
            indx(ii, 2) = 1
        case(1) ! p orbital
            has_p = .TRUE.
            if (kk == -1) then ! py 
                indx(ii, 2) = 2
            elseif (kk == 0) then ! pz
                indx(ii, 2) = 3
            elseif (kk == 1) then ! px
                indx(ii, 2) = 4
            endif
        case(2) ! d orbital
            has_d = .TRUE.
            if (kk == -2) then ! dxy
                indx(ii, 2) = 5
            elseif (kk == -1) then ! dyz
                indx(ii, 2) = 6
            elseif (kk == 0) then ! dz2
                indx(ii, 2) = 7
            elseif (kk == 1) then ! dxz
                indx(ii, 2) = 8
            elseif (kk == 2) then ! dx2-y2
                indx(ii, 2) = 9
            end if
        end select
    enddo
    close(2)
    ! This sets indx(atomic species, orbital, zeta).  
    ! The species follow our input file and the orbitals follow above.

    !-----------Makes helpful arrays for analysis-----------!
    ! Finds radii between atoms in the supercell
    allocate(dists(numDim), lmn(numDim, 3))
    do ii = 1, numDim
        dists(ii) = norm2(xij(ii, :))
        lmn(ii, :) = xij(ii, :)/dists(ii)
    enddo ! ii

    !------Makes ordered list of radii and interactions------!
    ! First makes un-ordered list of radii and interactions
    distsU = 0 ! initialize
    numInters = 0 ! will index unique interactions
    do ii = 1, no_u ! cycles through unit cell orbitals
        do t = 1, numh(ii)
            kk = listhptr(ii)+t ! cycles through all sparse elements
            jj = indxuo(listh(kk)) ! cycles through unit cell orbitals
            tc = .FALSE. ! stores if radius is found in distsU

            ! skips onsite terms and lower triangle
            if (dists(kk) < 1e-4 .or. jj < ii) then
                cycle
            endif

            do LL=1, numInters 
                ! First checks if radius is in array
                if (abs(dists(kk)-distsU(LL, 1)) < 1e-4) then 
                    ! Second checks for atomic interactions
                    if (indx(ii, 1) == nint(distsU(LL, 2)) .and. &
                        indx(jj, 1) == nint(distsU(LL, 3))) then
                    ! Third checks for zeta interactions
                        if (indx(ii, 3) == nint(distsU(LL, 4)) .and. &
                            indx(jj, 3) == nint(distsU(LL, 5))) then
                            tc = .TRUE. ! interaction is within list
                        endif
                endif
                endif
            enddo ! LL

            if (tc .eqv. .FALSE.) then
                numInters = numInters + 1
                distsU(numInters, 1) = dists(kk)  ! stores radius
                distsU(numInters, 2) = indx(ii, 1) ! stores species 1
                distsU(numInters, 3) = indx(jj, 1) ! stores species 2
                distsU(numInters, 4) = indx(ii, 3) ! stores zeta 1
                distsU(numInters, 5) = indx(jj, 3) ! stores zeta 2
            endif

        enddo ! t
    enddo ! ii
    
    ! Orders list of radii and interactions
    ! First orders radii
    do jj=1, numInters
        do ii=jj, numInters ! hits all non-zero elements of distsU
            if (distsU(ii, 1) < distsU(jj, 1)) then
                distsU(numInters+1, :) = distsU(jj, :) ! larger rad to temp
                distsU(jj, :) = distsU(ii, :) ! less radius to start
                distsU(ii, :) = distsU(numInters+1, :) ! larger rad to end
            endif
        enddo ! ii
    enddo ! jj
    ! Next cycles interactions (species then zeta)
    do LL=2, 5 ! cycles through elements of distsU
    do jj=1, numInters
    do ii=jj, numInters ! hits all non-zero elements of distsU
        if ((abs(distsU(ii, 1)-distsU(jj, 1)) < 1e-4) .and. &
            (norm2(distsU(ii, 2:LL-1)-distsU(jj, 2:LL-1)) < 1e-4) .and. &
            (nint(distsU(ii, LL))<nint(distsU(jj, LL)))) then
            distsU(numInters+1, :) = distsU(jj, :) ! larger radius to temp
            distsU(jj, :) = distsU(ii, :) ! less radius to start
            distsU(ii, :) = distsU(numInters+1, :) ! larger radius to end
        endif
    enddo ! ii
    enddo ! jj
    enddo ! LL

    reord(:) = (/0, 2, 3, 1/) ! reorders p orbitals for easier coding
    signs(:) = (/0, 1, -1, 1/) ! corrects for opposite sign on pz orbitals

    !------------------------------------------------------!
    !----------- Finds tight binding parameters -----------!
    !------------------------------------------------------!

    ! Allocates storage arrays depending on basis set
    if (has_s .eqv. .TRUE.) then ! if we have s orbitals
        allocate(Vss(numInters, 2))
        if (has_p .eqv. .TRUE.) then ! if we have s and p orbitals
            allocate(Vsp(numInters, 2), Vps(numInters, 2))
        endif
        if (has_d .eqv. .TRUE.) then ! if we have s and d orbitals
            allocate(Vsd(numInters, 2), Vds(numInters, 2))
        endif
    endif
    if (has_p .eqv. .TRUE.) then ! if we have p orbitals
        allocate(VppSig(numInters, 2), VppPi(numInters, 2))
        if (has_d .eqv. .TRUE.) then ! if we have p and d orbitals
            allocate(VpdSig(numInters, 2), VpdPi(numInters, 2))
            allocate(VdpSig(numInters, 2), VdpPi(numInters, 2))
        endif
    endif
    if (has_d .eqv. .TRUE.) then ! if we have d orbitals
        allocate(VddSig(numInters, 2), VddPi(numInters, 2))
        allocate(VddDelt(numInters, 2))
        VddSig  = 0 ! will store tight-binding parameters
        VddPi   = 0 ! will store tight-binding parameters
        VddDelt = 0 ! will store tight-binding parameters
    endif

    allocate(Vonsite(no_u, no_u, 2), ind_store(numDim, 7))
    ! ind_store(sparse indx, orb type1-2, species1-2, zeta1-2)

    ! Finds onsite terms
    Vonsite = 0 ! will store onsite terms
    do ii = 1, no_u ! cycles through unit cell orbitals ! row
        do t = 1, numh(ii)
            kk = listhptr(ii)+t ! cycles through all sparse elements
            jj = indxuo(listh(kk)) ! cycles through unit cell orbs ! col

            if (dists(kk) < 1e-4) then
                Vonsite(ii, jj, 1) = hsparse(kk)
                Vonsite(ii, jj, 2) = ssparse(kk)
            endif
        enddo ! t
    enddo ! ii

    ! Finds any s-s interactions
    if (has_s .eqv. .TRUE.) then 
    ind_store = 0 ! will store indices for data points
    Vss = 0 ! will store tight-binding parameters
    LL = 0 ! will store number of data points to analyse
    do ii = 1, no_u ! cycles through unit cell orbitals ! row
        do t = 1, numh(ii)
            kk = listhptr(ii)+t ! cycles through all sparse elements
            jj = indxuo(listh(kk)) ! cycles through unit cell orbs ! col

            if (indx(ii, 2) == 1) then ! row is s-orbital
            if (indx(jj, 2) == 1) then ! col is s-orbital
                LL = LL + 1 ! count of elements

                ind_store(LL, 1) = kk ! stores index for data point
                ind_store(LL, 2) = indx(ii, 2) ! stores orbital for row
                ind_store(LL, 3) = indx(jj, 2) ! stores orbital for col
                ind_store(LL, 4:5) = (/indx(ii, 1), indx(jj, 1)/) ! species
                ind_store(LL, 6:7) = (/indx(ii, 3), indx(jj, 3)/) ! zetas
            endif
            endif
        enddo ! t
    enddo ! ii
    
    do ii = 1, numInters ! cycles through all possible interactions
        t = 0 ! will store number of data points per interaction
        do jj = 1, LL ! cycles through indices for data point
            distsU(numInters+1, 1) = dists(ind_store(jj, 1))
            distsU(numInters+1, 2:5) = ind_store(jj, 4:7)
            ! runs if distances, atomic species, and zetas match
            if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then
                t = t + 1 ! stores number of points
                Vss(ii, 1) = Vss(ii, 1) + hsparse(ind_store(jj, 1))
                Vss(ii, 2) = Vss(ii, 2) + ssparse(ind_store(jj, 1))

            endif
        enddo ! jj
        if (t == 0) then ! if no interaction is calculated
            t = 1
        endif

        Vss(ii, :) = Vss(ii, :)/t ! takes average
    enddo ! ii
    endif
    
    ! Finds any s-p interactions
    if ((has_s .eqv. .TRUE.) .and. (has_p .eqv. .TRUE.)) then
    ind_store = 0 ! will store indices for data points
    Vsp = 0 ! will store tight-binding parameters
    LL = 0 ! will store number of data points to analyse
    do ii = 1, no_u ! cycles through unit cell orbitals ! row
        do t = 1, numh(ii)
            kk = listhptr(ii)+t ! cycles through all sparse elements
            jj = indxuo(listh(kk)) ! cycles through unit cell orbs ! col

            if (indx(ii, 2) == 1) then ! row is s-orbital
            if (2 <= indx(jj, 2).and.indx(jj, 2) <= 4) then ! col is p-orb
                LL = LL + 1 ! count of elements

                ind_store(LL, 1) = kk ! stores index for data point
                ind_store(LL, 2) = indx(ii, 2) ! stores orbital for row
                ind_store(LL, 3) = indx(jj, 2) ! stores orbital for col
                ind_store(LL, 4:5) = (/indx(ii, 1), indx(jj, 1)/) ! species
                ind_store(LL, 6:7) = (/indx(ii, 3), indx(jj, 3)/) ! zetas
            endif
            endif
        enddo ! t
    enddo ! ii
    do ii = 1, numInters ! cycles through all possible s-p interactions
        t = 0 ! will store number of data points per interaction
        do jj = 1, LL ! cycles through indices for data point
            distsU(numInters+1, 1) = dists(ind_store(jj, 1))
            distsU(numInters+1, 2:5) = ind_store(jj, 4:7)
            ! runs if distances, atomic species, and zetas match
            if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then
                select case(ind_store(jj, 3)) ! handles each p orbital
                case(2) ! py orbital
                    if (abs(lmn(ind_store(jj, 1), 2)) < 1e-4) then
                        cycle ! if y-coordinate is zero for vector
                    endif

                    t = t + 1 ! stores number of points
                    Vsp(ii, 1) = Vsp(ii, 1) + &
                        hsparse(ind_store(jj, 1))/lmn(ind_store(jj, 1), 2)
                    Vsp(ii, 2) = Vsp(ii, 2) + &
                        ssparse(ind_store(jj, 1))/lmn(ind_store(jj, 1), 2)

                case(3) ! pz orbital
                    if (abs(lmn(ind_store(jj, 1), 3)) < 1e-4) then
                        cycle ! if z-coordinate is zero for vector
                    endif

                    t = t + 1 ! stores number of points
                    Vsp(ii, 1) = Vsp(ii, 1) - &
                        hsparse(ind_store(jj, 1))/lmn(ind_store(jj, 1), 3)
                    Vsp(ii, 2) = Vsp(ii, 2) - &
                        ssparse(ind_store(jj, 1))/lmn(ind_store(jj, 1), 3)

                case(4) ! px orbital
                    if (abs(lmn(ind_store(jj, 1), 1)) < 1e-4) then
                        cycle ! if x-coordinate is zero for vector
                    endif

                    t = t + 1 ! stores number of points
                    Vsp(ii, 1) = Vsp(ii, 1) + &
                        hsparse(ind_store(jj, 1))/lmn(ind_store(jj, 1), 1)
                    Vsp(ii, 2) = Vsp(ii, 2) + &
                        ssparse(ind_store(jj, 1))/lmn(ind_store(jj, 1), 1)

                end select
            endif
        enddo ! jj
        if (t == 0) then ! if no interaction is calculated
            t = 1
        endif

        Vsp(ii, :) = Vsp(ii, :)/t ! takes average
    enddo ! ii
    
    ! Finds any p-s interactions
    ind_store = 0 ! will store indices for data points
    Vps = 0 ! will store tight-binding parameters
    LL = 0 ! will store number of data points to analyse
    do ii = 1, no_u ! cycles through unit cell orbitals ! row
        do t = 1, numh(ii)
            kk = listhptr(ii)+t ! cycles through all sparse elements
            jj = indxuo(listh(kk)) ! cycles through unit cell orbs ! col

            if (2 <= indx(ii, 2).and.indx(ii, 2) <= 4) then ! row is p-orb
            if (indx(jj, 2) == 1) then ! col is s-orbital
                LL = LL + 1 ! count of elements

                ind_store(LL, 1) = kk ! stores index for data point
                ind_store(LL, 2) = indx(ii, 2) ! stores orbital for row
                ind_store(LL, 3) = indx(jj, 2) ! stores orbital for col
                ind_store(LL, 4:5) = (/indx(ii, 1), indx(jj, 1)/) ! species
                ind_store(LL, 6:7) = (/indx(ii, 3), indx(jj, 3)/) ! zetas
            endif
            endif
        enddo ! t
    enddo ! ii
    do ii = 1, numInters ! cycles through all possible s-p interactions
        t = 0 ! will store number of data points per interaction
        do jj = 1, LL ! cycles through indices for data point
            distsU(numInters+1, 1) = dists(ind_store(jj, 1))
            distsU(numInters+1, 2:5) = ind_store(jj, 4:7)
            ! runs if distances, atomic species, and zetas match
            if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then
                select case(ind_store(jj, 2)) ! handles each p orbital
                case(2) ! py orbital
                    if (abs(lmn(ind_store(jj, 1), 2)) < 1e-4) then
                        cycle ! if y-coordinate is zero for vector
                    endif

                    t = t + 1 ! stores number of points
                    Vps(ii, 1) = Vps(ii, 1) + &
                        hsparse(ind_store(jj, 1))/lmn(ind_store(jj, 1), 2)
                    Vps(ii, 2) = Vps(ii, 2) + &
                        ssparse(ind_store(jj, 1))/lmn(ind_store(jj, 1), 2)
                case(3) ! pz orbital
                    if (abs(lmn(ind_store(jj, 1), 3)) < 1e-4) then
                        cycle ! if z-coordinate is zero for vector
                    endif

                    t = t + 1 ! stores number of points
                    Vps(ii, 1) = Vps(ii, 1) - &
                        hsparse(ind_store(jj, 1))/lmn(ind_store(jj, 1), 3)
                    Vps(ii, 2) = Vps(ii, 2) - &
                        ssparse(ind_store(jj, 1))/lmn(ind_store(jj, 1), 3)
                case(4) ! px orbital
                    if (abs(lmn(ind_store(jj, 1), 1)) < 1e-4) then
                        cycle ! if x-coordinate is zero for vector
                    endif

                    t = t + 1 ! stores number of points
                    Vps(ii, 1) = Vps(ii, 1) + &
                        hsparse(ind_store(jj, 1))/lmn(ind_store(jj, 1), 1)
                    Vps(ii, 2) = Vps(ii, 2) + &
                        ssparse(ind_store(jj, 1))/lmn(ind_store(jj, 1), 1)
                end select
            endif
        enddo ! jj
        if (t == 0) then ! if no interaction is calculated
            t = 1
        endif

        Vps(ii, :) = Vps(ii, :)/t ! takes average
    enddo ! ii
    endif
    
    ! Finds any p-p interactions
    if (has_p .eqv. .TRUE.) then
    ind_store = 0 ! will store indices for data points
    VppPi  = 0 ! will store tight-binding parameters
    VppSig = 0 ! will store tight-binding parameters
    LL = 0 ! will store number of data points to analyse
    do ii = 1, no_u ! cycles through unit cell orbitals ! row
        do t = 1, numh(ii)
            kk = listhptr(ii)+t ! cycles through all sparse elements
            jj = indxuo(listh(kk)) ! cycles through unit cell orbs ! col

            if (2 <= indx(ii, 2).and.indx(ii, 2) <= 4) then ! row is p-orb
            if (2 <= indx(jj, 2).and.indx(jj, 2) <= 4) then ! col is p-orb
                LL = LL + 1 ! count of elements

                ind_store(LL, 1) = kk ! stores index for data point
                ind_store(LL, 2) = indx(ii, 2) ! stores orbital for row
                ind_store(LL, 3) = indx(jj, 2) ! stores orbital for col
                ind_store(LL, 4:5) = (/indx(ii, 1), indx(jj, 1)/) ! species
                ind_store(LL, 6:7) = (/indx(ii, 3), indx(jj, 3)/) ! zetas
            endif
            endif
        enddo ! t
    enddo ! ii
    do ii = 1, numInters ! cycles through all possible interactions
        t  = 0 ! will store number of data points for averaging
        do jj = 1, LL-1 ! index of first row in system of equations
            distsU(numInters+1, 1) = dists(ind_store(jj, 1))
            distsU(numInters+1, 2:5) = ind_store(jj, 4:7)
        ! runs if first row matches correct interaction
        if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then

            ! At this point, we will solve our system of equations by
            ! analytically inverting a 2x2 matrix, (a b ; c d), where
            ! a, b are projection cosines from jj; soe(1:2)
            ! c, d are projection cosines from kk; soe(3:4)

            ! VppSig term
            soe(1) = signs(ind_store(jj, 2)) * &
                        lmn(ind_store(jj, 1), reord(ind_store(jj, 2))) *&
                        signs(ind_store(jj, 3)) * &
                        lmn(ind_store(jj, 1), reord(ind_store(jj, 3)))
            ! VppPi term
            if (ind_store(jj, 2)==ind_store(jj, 3)) then ! jj term is pi-pi
                soe(2)=1-lmn(ind_store(jj, 1), reord(ind_store(jj, 2)))**2
            else ! jj term is pi-pj
                soe(2) = -1*soe(1)
            endif

        do kk = jj+1, LL ! index of second row in system of equations
            distsU(numInters+1, 1) = dists(ind_store(kk, 1))
            distsU(numInters+1, 2:5) = ind_store(kk, 4:7)
        ! runs if second row matches correct interaction
        if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then

            ! VppSig term
            soe(3) = signs(ind_store(kk, 2)) * &
                        lmn(ind_store(kk, 1), reord(ind_store(kk, 2))) *&
                        signs(ind_store(kk, 3)) * &
                        lmn(ind_store(kk, 1), reord(ind_store(kk, 3)))

            ! VppPi term
            if (ind_store(kk, 2)==ind_store(kk, 3)) then ! kk term is pi-pi
                soe(4)=1-lmn(ind_store(kk, 1), reord(ind_store(kk, 2)))**2
            else ! kk term is pi-pj
                soe(4) = -1*soe(3)
            endif

            if (abs(soe(1)*soe(4)-soe(2)*soe(3)) < 1e-4) then
                cycle ! runs if no solution (singular system)
            endif
            
            ! We now have a guaranteed solution
            t = t + 1 

            ! Solves 2x2 system of equations analytically
            VppSig(ii, 1) = VppSig(ii, 1) + &
                (soe(4)*hsparse(ind_store(jj, 1)) - &
                 soe(2)*hsparse(ind_store(kk, 1))) / &
                (soe(1)*soe(4)-soe(2)*soe(3))
            VppSig(ii, 2) = VppSig(ii, 2) + &
                 (soe(4)*ssparse(ind_store(jj, 1)) - &
                  soe(2)*ssparse(ind_store(kk, 1))) / &
                 (soe(1)*soe(4)-soe(2)*soe(3))
            VppPi(ii, 1) = VppPi(ii, 1) + &
                  (soe(1)*hsparse(ind_store(kk, 1)) - &
                   soe(3)*hsparse(ind_store(jj, 1))) / &
                   (soe(1)*soe(4)-soe(2)*soe(3))
            VppPi(ii, 2) = VppPi(ii, 2) + &
                   (soe(1)*ssparse(ind_store(kk, 1)) - &
                    soe(3)*ssparse(ind_store(jj, 1))) / &
                   (soe(1)*soe(4)-soe(2)*soe(3))

            !debug block!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (ii == 1) then
            temp1 = (soe(4)*hsparse(ind_store(jj, 1)) - &
                     soe(2)*hsparse(ind_store(kk, 1))) / &
                    (soe(1)*soe(4)-soe(2)*soe(3))
            temp2 = (soe(4)*hsparse(ind_store(jj, 1)) - &
                    soe(2)*hsparse(ind_store(kk, 1))) / &
                    (soe(1)*soe(4)-soe(2)*soe(3))

            write(*, '(3f12.6, 8i3)') &
            distsU(ii, 1), temp1, temp2, ind_store(jj, 4:5), &
            ind_store(kk, 4:5), ind_store(jj, 2:3), ind_store(kk, 2:3)
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        endif ! second row
        enddo ! kk
        endif ! first row
        enddo ! jj
        if (t == 0) then ! if no interaction is calculated
            t = 1
        endif

        VppSig(ii, :) = VppSig(ii, :)/t ! takes average
        VppPi(ii, :) = VppPi(ii, :)/t  ! takes average
    enddo ! ii
    endif
    
    s3 = sqrt(3.0)
    ! Finds any s-d interactions
    if ((has_s .eqv. .TRUE.) .and. (has_d .eqv. .TRUE.)) then
    ind_store = 0 ! will store indices for data points
    Vsd = 0 ! will store tight-binding parameters
    LL = 0 ! will store number of data points to analyse
    do ii = 1, no_u ! cycles through unit cell orbitals ! row
        do t = 1, numh(ii)
            kk = listhptr(ii)+t ! cycles through all sparse elements
            jj = indxuo(listh(kk)) ! cycles through unit cell orbs ! col

            if (indx(ii, 2) == 1) then ! row is s-orbital
            if (5 <= indx(jj, 2)) then ! col is d-orb
                LL = LL + 1 ! count of elements

                ind_store(LL, 1) = kk ! stores index for data point
                ind_store(LL, 2) = indx(ii, 2) ! stores orbital for row
                ind_store(LL, 3) = indx(jj, 2) ! stores orbital for col
                ind_store(LL, 4:5) = (/indx(ii, 1), indx(jj, 1)/) ! species
                ind_store(LL, 6:7) = (/indx(ii, 3), indx(jj, 3)/) ! zetas
            endif
            endif
        enddo ! t
    enddo ! ii
    do ii = 1, numInters ! cycles through all possible s-d interactions
        t = 0 ! will store number of data points per interaction
        do jj = 1, LL ! cycles through indices for data point
            distsU(numInters+1, 1) = dists(ind_store(jj, 1))
            distsU(numInters+1, 2:5) = ind_store(jj, 4:7)
            ! runs if distances, atomic species, and zetas match
            if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then
                select case(ind_store(jj, 3)) ! handles each d orbital
                case(5) ! dxy
                    if (abs(lmn(ind_store(jj, 1), 1)) < 1e-4 .or. &
                        abs(lmn(ind_store(jj, 1), 2)) < 1e-4) then
                        cycle ! if x or y-coordinate is zero for vector
                    endif

                    t = t + 1 ! stores number of points
                    Vsd(ii, 1) = Vsd(ii, 1) + hsparse(ind_store(jj, 1))/&
                        s3/lmn(ind_store(jj, 1), 1)/lmn(ind_store(jj, 1), 2)
                    Vsd(ii, 2) = Vsd(ii, 2) + ssparse(ind_store(jj, 1))/&
                        s3/lmn(ind_store(jj, 1), 1)/lmn(ind_store(jj, 1), 2)
                case(6) ! dyz
                    if (abs(lmn(ind_store(jj, 1), 2)) < 1e-4 .or. &
                        abs(lmn(ind_store(jj, 1), 3)) < 1e-4) then
                        cycle ! if y or z-coordinate is zero for vector
                    endif

                    t = t + 1 ! stores number of points
                    Vsd(ii, 1) = Vsd(ii, 1) - hsparse(ind_store(jj, 1))/&
                        s3/lmn(ind_store(jj, 1), 2)/lmn(ind_store(jj, 1), 3)
                    Vsd(ii, 2) = Vsd(ii, 2) - ssparse(ind_store(jj, 1))/&
                        s3/lmn(ind_store(jj, 1), 2)/lmn(ind_store(jj, 1), 3)
                case(7) ! dz^2
                    L = lmn(ind_store(jj, 1), 3)**2 - &
                    (lmn(ind_store(jj, 1), 1)**2 + &
                    lmn(ind_store(jj, 1), 2)**2)/2 ! storage variable
                    if (abs(L) < 1e-4) then
                            cycle
                    endif

                    t = t + 1
                    Vsd(ii, 1) = Vsd(ii, 1) + hsparse(ind_store(jj, 1))/L
                    Vsd(ii, 2) = Vsd(ii, 2) + ssparse(ind_store(jj, 1))/L
                case(8) ! dxz
                    if (abs(lmn(ind_store(jj, 1), 1)) < 1e-4 .or. &
                        abs(lmn(ind_store(jj, 1), 3)) < 1e-4) then
                        cycle ! if x or z-coordinate is zero for vector
                    endif

                    t = t + 1 ! stores number of points
                    Vsd(ii, 1) = Vsd(ii, 1) - hsparse(ind_store(jj, 1))/&
                        s3/lmn(ind_store(jj, 1), 1)/lmn(ind_store(jj, 1), 3)
                    Vsd(ii, 2) = Vsd(ii, 2) - ssparse(ind_store(jj, 1))/&
                        s3/lmn(ind_store(jj, 1), 1)/lmn(ind_store(jj, 1), 3)
                case(9) ! dx^2-y^2
                    L = s3*(lmn(ind_store(jj, 1), 1)**2 - &
                            lmn(ind_store(jj, 1), 2)**2)/2 ! storage
                    if (abs(L) < 1e-4) then
                        cycle
                    endif

                    t = t + 1
                    Vsd(ii, 1) = Vsd(ii, 1) + hsparse(ind_store(jj, 1))/L
                    Vsd(ii, 2) = Vsd(ii, 2) + ssparse(ind_store(jj, 1))/L
                end select
            endif
        enddo ! jj
        if (t == 0) then ! if no interaction is calculated
            t = 1
        endif

        Vsd(ii, :) = Vsd(ii, :)/t ! takes average
    enddo ! ii

    ! Finds any d-s interactions
    ind_store = 0 ! will store indices for data points
    Vds = 0 ! will store tight-binding parameters
    LL = 0 ! will store number of data points to analyse
    do ii = 1, no_u ! cycles through unit cell orbitals ! row
        do t = 1, numh(ii)
            kk = listhptr(ii)+t ! cycles through all sparse elements
            jj = indxuo(listh(kk)) ! cycles through unit cell orbs ! col

            if (5 <= indx(ii, 2)) then ! row is d-orbital
            if (indx(jj, 2) == 1) then ! col is s-orbital
                LL = LL + 1 ! count of elements

                ind_store(LL, 1) = kk ! stores index for data point
                ind_store(LL, 2) = indx(ii, 2) ! stores orbital for row
                ind_store(LL, 3) = indx(jj, 2) ! stores orbital for col
                ind_store(LL, 4:5) = (/indx(ii, 1), indx(jj, 1)/) ! species
                ind_store(LL, 6:7) = (/indx(ii, 3), indx(jj, 3)/) ! zetas
            endif
            endif
        enddo ! t
    enddo ! ii
    do ii = 1, numInters ! cycles through all possible d-s interactions
        t = 0 ! will store number of data points per interaction
        do jj = 1, LL ! cycles through indices for data point
            distsU(numInters+1, 1) = dists(ind_store(jj, 1))
            distsU(numInters+1, 2:5) = ind_store(jj, 4:7)
            ! runs if distances, atomic species, and zetas match
            if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then
                select case(ind_store(jj, 2)) ! handles each d orbital
                case(5) ! dxy
                    if (abs(lmn(ind_store(jj, 1), 1)) < 1e-4 .or. &
                        abs(lmn(ind_store(jj, 1), 2)) < 1e-4) then
                        cycle ! if x or y-coordinate is zero for vector
                    endif

                    t = t + 1 ! stores number of points
                    Vds(ii, 1) = Vds(ii, 1) + hsparse(ind_store(jj, 1))/&
                        s3/lmn(ind_store(jj, 1), 1)/lmn(ind_store(jj, 1), 2)
                    Vds(ii, 2) = Vds(ii, 2) + ssparse(ind_store(jj, 1))/&
                        s3/lmn(ind_store(jj, 1), 1)/lmn(ind_store(jj, 1), 2)
                case(6) ! dyz
                    if (abs(lmn(ind_store(jj, 1), 2)) < 1e-4 .or. &
                        abs(lmn(ind_store(jj, 1), 3)) < 1e-4) then
                        cycle ! if y or z-coordinate is zero for vector
                    endif

                    t = t + 1 ! stores number of points
                    Vds(ii, 1) = Vds(ii, 1) - hsparse(ind_store(jj, 1))/&
                        s3/lmn(ind_store(jj, 1), 2)/lmn(ind_store(jj, 1), 3)
                    Vds(ii, 2) = Vds(ii, 2) - ssparse(ind_store(jj, 1))/&
                        s3/lmn(ind_store(jj, 1), 2)/lmn(ind_store(jj, 1), 3)
                case(7) ! dz^2
                    L = lmn(ind_store(jj, 1), 3)**2 - &
                    (lmn(ind_store(jj, 1), 1)**2 + &
                    lmn(ind_store(jj, 1), 2)**2)/2 ! storage variable
                    if (abs(L) < 1e-4) then
                            cycle
                    endif

                    t = t + 1
                    Vds(ii, 1) = Vds(ii, 1) + hsparse(ind_store(jj, 1))/L
                    Vds(ii, 2) = Vds(ii, 2) + ssparse(ind_store(jj, 1))/L
                case(8) ! dxz
                    if (abs(lmn(ind_store(jj, 1), 1)) < 1e-4 .or. &
                        abs(lmn(ind_store(jj, 1), 3)) < 1e-4) then
                        cycle ! if x or z-coordinate is zero for vector
                    endif

                    t = t + 1 ! stores number of points
                    Vds(ii, 1) = Vds(ii, 1) - hsparse(ind_store(jj, 1))/&
                        s3/lmn(ind_store(jj, 1), 1)/lmn(ind_store(jj, 1), 3)
                    Vds(ii, 2) = Vds(ii, 2) - ssparse(ind_store(jj, 1))/&
                        s3/lmn(ind_store(jj, 1), 1)/lmn(ind_store(jj, 1), 3)
                case(9) ! dx^2-y^2
                    L = s3*(lmn(ind_store(jj, 1), 1)**2 - &
                            lmn(ind_store(jj, 1), 2)**2)/2 ! storage
                    if (abs(L) < 1e-4) then
                        cycle
                    endif

                    t = t + 1
                    Vds(ii, 1) = Vds(ii, 1) + hsparse(ind_store(jj, 1))/L
                    Vds(ii, 2) = Vds(ii, 2) + ssparse(ind_store(jj, 1))/L
                end select
            endif
        enddo ! jj
        if (t == 0) then ! if no interaction is calculated
            t = 1
        endif

        Vds(ii, :) = Vds(ii, :)/t ! takes average
    enddo ! ii
    endif
    
    ! Finds any p-d interactions
    if ((has_p .eqv. .TRUE.) .and. (has_d .eqv. .TRUE.)) then
    ind_store = 0 ! will store indices for data points
    VpdPi  = 0 ! will store tight-binding parameters
    VpdSig = 0 ! will store tight-binding parameters
    LL = 0 ! will store number of data points to analyse
    do ii = 1, no_u ! cycles through unit cell orbitals ! row
        do t = 1, numh(ii)
            kk = listhptr(ii)+t ! cycles through all sparse elements
            jj = indxuo(listh(kk)) ! cycles through unit cell orbs ! col

            if (2 <= indx(ii, 2).and.indx(ii, 2) <= 4) then ! row is p-orb
            if (5 <= indx(jj, 2)) then ! col is d-orb
                LL = LL + 1 ! count of elements

                ind_store(LL, 1) = kk ! stores index for data point
                ind_store(LL, 2) = indx(ii, 2) ! stores orbital for row
                ind_store(LL, 3) = indx(jj, 2) ! stores orbital for col
                ind_store(LL, 4:5) = (/indx(ii, 1), indx(jj, 1)/) ! species
                ind_store(LL, 6:7) = (/indx(ii, 3), indx(jj, 3)/) ! zetas
            endif
            endif
        enddo ! t
    enddo ! ii
    do ii = 1, numInters ! cycles through all possible interactions
        t  = 0 ! will store number of data points for averaging
        do jj = 1, LL-1 ! index of first row in system of equations
            distsU(numInters+1, 1) = dists(ind_store(jj, 1))
            distsU(numInters+1, 2:5) = ind_store(jj, 4:7)
        ! runs if first row matches correct interaction
        if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then

            ! At this point, we will solve our system of equations by
            ! analytically inverting a 2x2 matrix, (a b ; c d), where
            ! a, b are projection cosines from jj; soe(1:2)
            ! c, d are projection cosines from kk; soe(3:4)

            L = lmn(ind_store(jj, 1), 1)
            M = lmn(ind_store(jj, 1), 2) ! redefining these saves space
            N = lmn(ind_store(jj, 1), 3)

            if (ind_store(jj, 3) == 5) then ! dxy
                select case(ind_store(jj, 2))
                case(2) ! py
                    soe(1) = s3*L*M*M
                    soe(2) = L*(1-2*M*M)
                case(3) ! pz
                    soe(1) = s3*L*M*N*(-1)
                    soe(2) = -2*L*M*N*(-1)
                case(4) ! px
                    soe(1) = s3*L*L*M
                    soe(2) = M*(1-2*L*L)
                end select
            elseif (ind_store(jj, 3) == 6) then ! dyz
                select case(ind_store(jj, 2))
                case(2) ! py
                    soe(1) = s3*M*M*N*(-1)
                    soe(2) = N*(1-2*M*M)*(-1)
                case(3) ! pz
                    soe(1) = s3*M*N*N 
                    soe(2) = M*(1-2*N*N)
                case(4) ! px
                    soe(1) = s3*L*M*N*(-1)
                    soe(2) = -2*L*M*N*(-1)
                end select
            elseif (ind_store(jj, 3) == 7) then ! dz^2
                select case(ind_store(jj, 2))
                case(2) ! py
                    soe(1) = M*(N*N - 0.5*(L*L+M*M))
                    soe(2) = -1*s3*M*N*N
                case(3) ! pz
                    soe(1) = N*(N*N - 0.5*(L*L+M*M))*(-1)
                    soe(2) = s3*N*(L*L+M*M)*(-1)
                case(4) ! px
                    soe(1) = L*(N*N - 0.5*(L*L+M*M))
                    soe(2) = -1*s3*L*N*N
                end select
            elseif (ind_store(jj, 3) == 8) then ! dxz
                select case(ind_store(jj, 2))
                case(2) ! py
                    soe(1) = s3*L*M*N*(-1)
                    soe(2) = -2*L*M*N*(-1)
                case(3) ! pz
                    soe(1) = s3*L*N*N 
                    soe(2) = L*(1-2*N*N)
                case(4) ! px
                    soe(1) = s3*L*L*N*(-1)
                    soe(2) = N*(1-2*L*L)*(-1)
                end select
            elseif (ind_store(jj, 3) == 9) then ! dx2-y2
                select case(ind_store(jj, 2))
                case(2) ! py
                    soe(1) = s3*M*(L*L-M*M)/2
                    soe(2) = -1*M*(1+L*L-M*M)
                case(3) ! pz
                    soe(1) = s3*N*(L*L-M*M)/2*(-1)
                    soe(2) = -1*N*(L*L-M*M)*(-1)
                case(4) ! px
                    soe(1) = s3*L*(L*L-M*M)/2
                    soe(2) = L*(1-L*L+M*M)
                end select
            else 
                write(*, '(a38)')'Error: unexpected orbital in pd system'
                call exit(0)
            endif

        do kk = jj+1, LL ! index of second row in system of equations
            distsU(numInters+1, 1) = dists(ind_store(kk, 1))
            distsU(numInters+1, 2:5) = ind_store(kk, 4:7)
        ! runs if second row matches correct interaction
        if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then

            L = lmn(ind_store(kk, 1), 1)
            M = lmn(ind_store(kk, 1), 2) ! redefining these saves space
            N = lmn(ind_store(kk, 1), 3)

            if (ind_store(kk, 3) == 5) then ! dxy
                select case(ind_store(kk, 2))
                case(2) ! py
                    soe(3) = s3*L*M*M
                    soe(4) = L*(1-2*M*M)
                case(3) ! pz
                    soe(3) = s3*L*M*N*(-1)
                    soe(4) = -2*L*M*N*(-1)
                case(4) ! px
                    soe(3) = s3*L*L*M
                    soe(4) = M*(1-2*L*L)
                end select
            elseif (ind_store(kk, 3) == 6) then ! dyz
                select case(ind_store(kk, 2))
                case(2) ! py
                    soe(3) = s3*M*M*N*(-1)
                    soe(4) = N*(1-2*M*M)*(-1)
                case(3) ! pz
                    soe(3) = s3*M*N*N 
                    soe(4) = M*(1-2*N*N)
                case(4) ! px
                    soe(3) = s3*L*M*N*(-1)
                    soe(4) = -2*L*M*N*(-1)
                end select
            elseif (ind_store(kk, 3) == 7) then ! dz^2
                select case(ind_store(kk, 2))
                case(2) ! py
                    soe(3) = M*(N*N - 0.5*(L*L+M*M))
                    soe(4) = -1*s3*M*N*N
                case(3) ! pz
                    soe(3) = N*(N*N - 0.5*(L*L+M*M))*(-1)
                    soe(4) = s3*N*(L*L+M*M)*(-1)
                case(4) ! px
                    soe(3) = L*(N*N - 0.5*(L*L+M*M))
                    soe(4) = -1*s3*L*N*N
                end select
            elseif (ind_store(kk, 3) == 8) then ! dxz
                select case(ind_store(kk, 2))
                case(2) ! py
                    soe(3) = s3*L*M*N*(-1)
                    soe(4) = -2*L*M*N*(-1)
                case(3) ! pz
                    soe(3) = s3*L*N*N 
                    soe(4) = L*(1-2*N*N)
                case(4) ! px
                    soe(3) = s3*L*L*N*(-1)
                    soe(4) = N*(1-2*L*L)*(-1)
                end select
            elseif (ind_store(kk, 3) == 9) then ! dx2-y2
                select case(ind_store(kk, 2))
                case(2) ! py
                    soe(3) = s3*M*(L*L-M*M)/2
                    soe(4) = -1*M*(1+L*L-M*M)
                case(3) ! pz
                    soe(3) = s3*N*(L*L-M*M)/2*(-1)
                    soe(4) = -1*N*(L*L-M*M)*(-1)
                case(4) ! px
                    soe(3) = s3*L*(L*L-M*M)/2
                    soe(4) = L*(1-L*L+M*M)
                end select
            else 
                write(*, '(a38)')'Error: unexpected orbital in pd system'
                call exit(0)
            endif

            ! finds determinant for system of equations
            L = soe(1)*soe(4)-soe(2)*soe(3)
            if (abs(L) < 1e-2) then
                cycle ! skips if no solution (singular system)
            endif
            
            ! We now have a guaranteed solution
            t = t + 1 
            t2 = t2 + 1

            ! Solves 2x2 system of equations analytically
            VpdSig(ii, 1) = VpdSig(ii, 1) + &
                (soe(4)*hsparse(ind_store(jj, 1)) - &
                 soe(2)*hsparse(ind_store(kk, 1))) / L
            VpdSig(ii, 2) = VpdSig(ii, 2) + &
                 (soe(4)*ssparse(ind_store(jj, 1)) - &
                  soe(2)*ssparse(ind_store(kk, 1))) / L
            VpdPi(ii, 1) = VpdPi(ii, 1) + &
                  (soe(1)*hsparse(ind_store(kk, 1)) - &
                   soe(3)*hsparse(ind_store(jj, 1))) / L
            VpdPi(ii, 2) = VpdPi(ii, 2) + &
                   (soe(1)*ssparse(ind_store(kk, 1)) - &
                    soe(3)*ssparse(ind_store(jj, 1))) / L

        endif ! second row
        enddo ! kk
        endif ! first row
        enddo ! jj
        if (t == 0) then ! if no interaction is calculated
            t = 1
        endif

        VpdSig(ii, :) = VpdSig(ii, :)/t ! takes average
        VpdPi(ii, :) = VpdPi(ii, :)/t   ! takes average
    enddo ! ii
    
    ! Finds any d-p interactions
    ind_store = 0 ! will store indices for data points
    VdpPi  = 0 ! will store tight-binding parameters
    VdpSig = 0 ! will store tight-binding parameters
    LL = 0 ! will store number of data points to analyse
    do ii = 1, no_u ! cycles through unit cell orbitals ! row
        do t = 1, numh(ii)
            kk = listhptr(ii)+t ! cycles through all sparse elements
            jj = indxuo(listh(kk)) ! cycles through unit cell orbs ! col

            if (5 <= indx(ii, 2)) then ! row is d-orb
            if (2 <= indx(jj, 2).and.indx(jj, 2) <= 4) then ! col is p-orb
                LL = LL + 1 ! count of elements

                ind_store(LL, 1) = kk ! stores index for data point
                ind_store(LL, 2) = indx(ii, 2) ! stores orbital for row
                ind_store(LL, 3) = indx(jj, 2) ! stores orbital for col
                ind_store(LL, 4:5) = (/indx(ii, 1), indx(jj, 1)/) ! species
                ind_store(LL, 6:7) = (/indx(ii, 3), indx(jj, 3)/) ! zetas
            endif
            endif
        enddo ! t
    enddo ! ii
    do ii = 1, numInters ! cycles through all possible interactions
        t  = 0 ! will store number of data points for averaging
        do jj = 1, LL-1 ! index of first row in system of equations
            distsU(numInters+1, 1) = dists(ind_store(jj, 1))
            distsU(numInters+1, 2:5) = ind_store(jj, 4:7)
        ! runs if first row matches correct interaction
        if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then

            ! At this point, we will solve our system of equations by
            ! analytically inverting a 2x2 matrix, (a b ; c d), where
            ! a, b are projection cosines from jj; soe(1:2)
            ! c, d are projection cosines from kk; soe(3:4)

            L = lmn(ind_store(jj, 1), 1)
            M = lmn(ind_store(jj, 1), 2) ! redefining these saves space
            N = lmn(ind_store(jj, 1), 3)

            if (ind_store(jj, 2) == 5) then ! dxy
                select case(ind_store(jj, 3))
                case(2) ! py
                    soe(1) = s3*L*M*M
                    soe(2) = L*(1-2*M*M)
                case(3) ! pz
                    soe(1) = s3*L*M*N*(-1)
                    soe(2) = -2*L*M*N*(-1)
                case(4) ! px
                    soe(1) = s3*L*L*M
                    soe(2) = M*(1-2*L*L)
                end select
            elseif (ind_store(jj, 2) == 6) then ! dyz
                select case(ind_store(jj, 3))
                case(2) ! py
                    soe(1) = s3*M*M*N*(-1)
                    soe(2) = N*(1-2*M*M)*(-1)
                case(3) ! pz
                    soe(1) = s3*M*N*N 
                    soe(2) = M*(1-2*N*N)
                case(4) ! px
                    soe(1) = s3*L*M*N*(-1)
                    soe(2) = -2*L*M*N*(-1)
                end select
            elseif (ind_store(jj, 2) == 7) then ! dz^2
                select case(ind_store(jj, 3))
                case(2) ! py
                    soe(1) = M*(N*N - 0.5*(L*L+M*M))
                    soe(2) = -1*s3*M*N*N
                case(3) ! pz
                    soe(1) = N*(N*N - 0.5*(L*L+M*M))*(-1)
                    soe(2) = s3*N*(L*L+M*M)*(-1)
                case(4) ! px
                    soe(1) = L*(N*N - 0.5*(L*L+M*M))
                    soe(2) = -1*s3*L*N*N
                end select
            elseif (ind_store(jj, 2) == 8) then ! dxz
                select case(ind_store(jj, 3))
                case(2) ! py
                    soe(1) = s3*L*M*N*(-1)
                    soe(2) = -2*L*M*N*(-1)
                case(3) ! pz
                    soe(1) = s3*L*N*N 
                    soe(2) = L*(1-2*N*N)
                case(4) ! px
                    soe(1) = s3*L*L*N*(-1)
                    soe(2) = N*(1-2*L*L)*(-1)
                end select
            elseif (ind_store(jj, 2) == 9) then ! dx2-y2
                select case(ind_store(jj, 3))
                case(2) ! py
                    soe(1) = s3*M*(L*L-M*M)/2
                    soe(2) = -1*M*(1+L*L-M*M)
                case(3) ! pz
                    soe(1) = s3*N*(L*L-M*M)/2*(-1)
                    soe(2) = -1*N*(L*L-M*M)*(-1)
                case(4) ! px
                    soe(1) = s3*L*(L*L-M*M)/2
                    soe(2) = L*(1-L*L+M*M)
                end select
            else 
                write(*, '(a38)')'Error: unexpected orbital in dp system'
                call exit(0)
            endif

        do kk = jj+1, LL ! index of second row in system of equations
            distsU(numInters+1, 1) = dists(ind_store(kk, 1))
            distsU(numInters+1, 2:5) = ind_store(kk, 4:7)
        ! runs if second row matches correct interaction
        if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then

            L = lmn(ind_store(kk, 1), 1)
            M = lmn(ind_store(kk, 1), 2) ! redefining these saves space
            N = lmn(ind_store(kk, 1), 3)

            if (ind_store(kk, 2) == 5) then ! dxy
                select case(ind_store(kk, 3))
                case(2) ! py
                    soe(3) = s3*L*M*M
                    soe(4) = L*(1-2*M*M)
                case(3) ! pz
                    soe(3) = s3*L*M*N*(-1)
                    soe(4) = -2*L*M*N*(-1)
                case(4) ! px
                    soe(3) = s3*L*L*M
                    soe(4) = M*(1-2*L*L)
                end select
            elseif (ind_store(kk, 2) == 6) then ! dyz
                select case(ind_store(kk, 3))
                case(2) ! py
                    soe(3) = s3*M*M*N*(-1)
                    soe(4) = N*(1-2*M*M)*(-1)
                case(3) ! pz
                    soe(3) = s3*M*N*N 
                    soe(4) = M*(1-2*N*N)
                case(4) ! px
                    soe(3) = s3*L*M*N*(-1)
                    soe(4) = -2*L*M*N*(-1)
                end select
            elseif (ind_store(kk, 2) == 7) then ! dz^2
                select case(ind_store(kk, 3))
                case(2) ! py
                    soe(3) = M*(N*N - 0.5*(L*L+M*M))
                    soe(4) = -1*s3*M*N*N
                case(3) ! pz
                    soe(3) = N*(N*N - 0.5*(L*L+M*M))*(-1)
                    soe(4) = s3*N*(L*L+M*M)*(-1)
                case(4) ! px
                    soe(3) = L*(N*N - 0.5*(L*L+M*M))
                    soe(4) = -1*s3*L*N*N
                end select
            elseif (ind_store(kk, 2) == 8) then ! dxz
                select case(ind_store(kk, 3))
                case(2) ! py
                    soe(3) = s3*L*M*N*(-1)
                    soe(4) = -2*L*M*N*(-1)
                case(3) ! pz
                    soe(3) = s3*L*N*N 
                    soe(4) = L*(1-2*N*N)
                case(4) ! px
                    soe(3) = s3*L*L*N*(-1)
                    soe(4) = N*(1-2*L*L)*(-1)
                end select
            elseif (ind_store(kk, 2) == 9) then ! dx2-y2
                select case(ind_store(kk, 3))
                case(2) ! py
                    soe(3) = s3*M*(L*L-M*M)/2
                    soe(4) = -1*M*(1+L*L-M*M)
                case(3) ! pz
                    soe(3) = s3*N*(L*L-M*M)/2*(-1)
                    soe(4) = -1*N*(L*L-M*M)*(-1)
                case(4) ! px
                    soe(3) = s3*L*(L*L-M*M)/2
                    soe(4) = L*(1-L*L+M*M)
                end select
            else 
                write(*, '(a38)')'Error: unexpected orbital in dp system'
                call exit(0)
            endif

            ! finds determinant for system of equations
            L = soe(1)*soe(4)-soe(2)*soe(3)
            if (abs(L) < 1e-2) then ! larger tol avoids wrong rounding
                cycle ! skips if no solution (singular system)
            endif
            
            ! We now have a guaranteed solution
            t = t + 1 

            ! Solves 2x2 system of equations analytically
            VdpSig(ii, 1) = VdpSig(ii, 1) + &
                (soe(4)*hsparse(ind_store(jj, 1)) - &
                 soe(2)*hsparse(ind_store(kk, 1))) / L
            VdpSig(ii, 2) = VdpSig(ii, 2) + &
                 (soe(4)*ssparse(ind_store(jj, 1)) - &
                  soe(2)*ssparse(ind_store(kk, 1))) / L
            VdpPi(ii, 1) = VdpPi(ii, 1) + &
                  (soe(1)*hsparse(ind_store(kk, 1)) - &
                   soe(3)*hsparse(ind_store(jj, 1))) / L
            VdpPi(ii, 2) = VdpPi(ii, 2) + &
                   (soe(1)*ssparse(ind_store(kk, 1)) - &
                    soe(3)*ssparse(ind_store(jj, 1))) / L

        endif ! second row
        enddo ! kk
        endif ! first row
        enddo ! jj
        if (t == 0) then ! if no interaction is calculated
            t = 1
        endif

        VdpSig(ii, :) = VdpSig(ii, :)/t ! takes average
        VdpPi(ii, :) = VdpPi(ii, :)/t   ! takes average
    enddo ! ii
    endif

    ! Finds any d-d interactions
    if (has_d .eqv. .TRUE.) then
    ind_store = 0 ! will store indices for data points
    VddPi  = 0  ! will store tight-binding parameters
    VddSig = 0  ! will store tight-binding parameters
    VddDelt = 0 ! will store tight-binding parameters
    LL = 0 ! will store number of data points to analyse
    do ii = 1, no_u ! cycles through unit cell orbitals ! row
        do t = 1, numh(ii)
            kk = listhptr(ii)+t ! cycles through all sparse elements
            jj = indxuo(listh(kk)) ! cycles through unit cell orbs ! col

            if (5 <= indx(ii, 2)) then ! row is d-orb
            if (5 <= indx(jj, 2)) then ! col is d-orb
                LL = LL + 1 ! count of elements

                ind_store(LL, 1) = kk ! stores index for data point
                ind_store(LL, 2) = indx(ii, 2) ! stores orbital for row
                ind_store(LL, 3) = indx(jj, 2) ! stores orbital for col
                ind_store(LL, 4:5) = (/indx(ii, 1), indx(jj, 1)/) ! species
                ind_store(LL, 6:7) = (/indx(ii, 3), indx(jj, 3)/) ! zetas
            endif
            endif
        enddo ! t
    enddo ! ii
    do ii = 1, numInters ! cycles through all possible interactions
        t  = 0 ! will store number of data points for averaging
        do jj = 1, LL-2 ! index of first row in system of equations
            distsU(numInters+1, 1) = dists(ind_store(jj, 1))
            distsU(numInters+1, 2:5) = ind_store(jj, 4:7)
        ! runs if first row matches correct interaction
        if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then

            ! At this point, we will solve our system of equations by
            ! analytically inverting a 3x3 matrix, 
            ! (a b c ; d e f ; g h i), where
            ! a, b, c are projection cosines from jj; soe(1:3)
            ! d, e, f are projection cosines from kk; soe(4:6)
            ! g, h, i are projection cosines from mm; soe(7:9)

            L = lmn(ind_store(jj, 1), 1)
            M = lmn(ind_store(jj, 1), 2) ! redefining these saves space
            N = lmn(ind_store(jj, 1), 3)

            if (ind_store(jj, 2) == 5) then ! dxy
                select case(ind_store(jj, 3))
                case(5) ! dxy
                    soe(1) = 3*L*L*M*M
                    soe(2) = L*L + M*M - 4*L*L*M*M
                    soe(3) = N*N + L*L*M*M
                case(6) ! dyz
                    soe(1) = 3*L*M*M*N*(-1)
                    soe(2) = L*N*(1 - 4*M*M)*(-1)
                    soe(3) = L*N*(M*M - 1)*(-1)
                case(7) ! dz^2
                    soe(1) = s3*L*M*(N*N-(L*L+M*M)/2)
                    soe(2) = -2*s3*L*M*N*N
                    soe(3) = s3*L*M*(1+N*N)/2
                case(8) ! dxz
                    soe(1) = 3*L*L*M*N*(-1)
                    soe(2) = M*N*(1 - 4*L*L)*(-1)
                    soe(3) = M*N*(L*L - 1)*(-1)
                case(9) ! dx2-y2
                    soe(1) = 3*L*M*(L*L - M*M)/2
                    soe(2) = 2*L*M*(M*M - L*L)
                    soe(3) = L*M*(L*L - M*M)/2
                end select
            elseif (ind_store(jj, 2) == 6) then ! dyz
                select case(ind_store(jj, 3))
                case(5) ! dxy
                    soe(1) = 3*L*M*M*N*(-1)
                    soe(2) = L*N*(1 - 4*M*M)*(-1)
                    soe(3) = L*N*(M*M - 1)*(-1) 
                case(6) ! dyz
                    soe(1) = 3*M*M*N*N
                    soe(2) = M*M + N*N - 4*M*M*N*N
                    soe(3) = L*L + M*M*N*N
                case(7) ! dz^2
                    soe(1) = s3*M*N*(N*N-(L*L+M*M)/2)*(-1)
                    soe(2) = s3*M*N*(L*L+M*M-N*N)*(-1)
                    soe(3) = s3*M*N*(L*L+M*M)/2*(-1)
                case(8) ! dxz
                    soe(1) = 3*L*M*N*N
                    soe(2) = L*M*(1-4*N*N)
                    soe(3) = L*M*(N*N-1)
                case(9) ! dx2-y2
                    soe(1) = 1.5*M*N*(L*L-M*M)*(-1)
                    soe(2) = -M*N*(1+2*(L*L-M*M))*(-1)
                    soe(3) = M*N*(1+0.5*(L*L-M*M))*(-1)
                end select
            elseif (ind_store(jj, 2) == 7) then ! dz^2
                select case(ind_store(jj, 3))
                case(5) ! dxy
                    soe(1) = s3*L*M*(N*N-0.5*(L*L+M*M))
                    soe(2) = -2*s3*L*M*N*N
                    soe(3) = (s3/2)*L*M*(1+N*N)
                case(6) ! dyz
                    soe(1) = s3*M*N*(N*N-(L*L+M*M)/2)*(-1)
                    soe(2) = s3*M*N*(L*L+M*M-N*N)*(-1)
                    soe(3) = s3*M*N*(L*L+M*M)/2*(-1)
                case(7) ! dz^2
                    soe(1) = (N*N-0.5*(L*L+M*M))**2
                    soe(2) = 3*N*N*(L*L+M*M)
                    soe(3) = 0.75*(L*L+M*M)**2
                case(8) ! dxz
                    soe(1) = s3*L*N*(N*N-0.5*(L*L+M*M))*(-1)
                    soe(2) = s3*L*N*(L*L+M*M-N*N)*(-1)
                    soe(3) = -(s3/2)*L*N*(L*L+M*M)*(-1)
                case(9) ! dx2-y2
                    soe(1) = (s3/2)*(L*L-M*M)*(N*N-0.5*(L*L+M*M))
                    soe(2) = s3*N*N*(M*M-L*L)
                    soe(3) = (s3/4)*(1+N*N)*(L*L-M*M)
                end select
            elseif (ind_store(jj, 2) == 8) then ! dxz
                select case(ind_store(jj, 3))
                case(5) ! dxy
                    soe(1) = 3*L*L*M*N*(-1)
                    soe(2) = M*N*(1 - 4*L*L)*(-1)
                    soe(3) = M*N*(L*L - 1)*(-1)
                case(6) ! dyz
                    soe(1) = 3*L*M*N*N
                    soe(2) = L*M*(1-4*N*N)
                    soe(3) = L*M*(N*N-1)
                case(7) ! dz^2
                    soe(1) = s3*L*N*(N*N-0.5*(L*L+M*M))*(-1)
                    soe(2) = s3*L*N*(L*L+M*M-N*N)*(-1)
                    soe(3) = -(s3/2)*L*N*(L*L+M*M)*(-1)
                case(8) ! dxz
                    soe(1) = 3*L*L*N*N
                    soe(2) = L*L + N*N - 4*L*L*N*N
                    soe(3) = M*M + L*L*N*N
                case(9) ! dx2-y2
                    soe(1) = 1.5*N*L*(L*L-M*M)*(-1)
                    soe(2) = N*L*(1-2*(L*L-M*M))*(-1)
                    soe(3) = -N*L*(1-0.5*(L*L-M*M))*(-1)
                end select
            elseif (ind_store(jj, 2) == 9) then ! dx2-y2
                select case(ind_store(jj, 3))
                case(5) ! dxy
                    soe(1) = 3*L*M*(L*L - M*M)/2
                    soe(2) = 2*L*M*(M*M - L*L)
                    soe(3) = L*M*(L*L - M*M)/2
                case(6) ! dyz
                    soe(1) = 1.5*M*N*(L*L-M*M)*(-1)
                    soe(2) = -M*N*(1+2*(L*L-M*M))*(-1)
                    soe(3) = M*N*(1+0.5*(L*L-M*M))*(-1)
                case(7) ! dz^2
                    soe(1) = (s3/2)*(L*L-M*M)*(N*N-0.5*(L*L+M*M))
                    soe(2) = s3*N*N*(M*M-L*L)
                    soe(3) = (s3/4)*(1+N*N)*(L*L-M*M)
                case(8) ! dxz
                    soe(1) = 1.5*N*L*(L*L-M*M)*(-1)
                    soe(2) = N*L*(1-2*(L*L-M*M))*(-1)
                    soe(3) = -N*L*(1-0.5*(L*L-M*M))*(-1)
                case(9) ! dx2-y2
                    soe(1) = 0.75*(L*L-M*M)**2
                    soe(2) = (L*L+M*M - (L*L-M*M)**2)
                    soe(3) = (N*N+0.25*(L*L- M*M)**2)
                end select 
            else 
                write(*, '(a35)') 'Error: wrong orbital in dd system 1'
                call exit(0)
            endif

        do kk = jj+1, LL-1 ! index of second row in system of equations
            distsU(numInters+1, 1) = dists(ind_store(kk, 1))
            distsU(numInters+1, 2:5) = ind_store(kk, 4:7)
        ! runs if second row matches correct interaction
        if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then

            L = lmn(ind_store(kk, 1), 1)
            M = lmn(ind_store(kk, 1), 2) ! redefining these saves space
            N = lmn(ind_store(kk, 1), 3)

            if (ind_store(kk, 2) == 5) then ! dxy
                select case(ind_store(kk, 3))
                case(5) ! dxy
                    soe(4) = 3*L*L*M*M
                    soe(5) = L*L + M*M - 4*L*L*M*M
                    soe(6) = N*N + L*L*M*M
                case(6) ! dyz
                    soe(4) = 3*L*M*M*N*(-1)
                    soe(5) = L*N*(1 - 4*M*M)*(-1)
                    soe(6) = L*N*(M*M - 1)*(-1)
                case(7) ! dz^2
                    soe(4) = s3*L*M*(N*N-(L*L+M*M)/2)
                    soe(5) = -2*s3*L*M*N*N
                    soe(6) = s3*L*M*(1+N*N)/2
                case(8) ! dxz
                    soe(4) = 3*L*L*M*N*(-1)
                    soe(5) = M*N*(1 - 4*L*L)*(-1)
                    soe(6) = M*N*(L*L - 1)*(-1)
                case(9) ! dx2-y2
                    soe(4) = 3*L*M*(L*L - M*M)/2
                    soe(5) = 2*L*M*(M*M - L*L)
                    soe(6) = L*M*(L*L - M*M)/2
                end select
            elseif (ind_store(kk, 2) == 6) then ! dyz
                select case(ind_store(kk, 3))
                case(5) ! dxy
                    soe(4) = 3*L*M*M*N*(-1)
                    soe(5) = L*N*(1 - 4*M*M)*(-1)
                    soe(6) = L*N*(M*M - 1)*(-1) 
                case(6) ! dyz
                    soe(4) = 3*M*M*N*N
                    soe(5) = M*M + N*N - 4*M*M*N*N
                    soe(6) = L*L + M*M*N*N
                case(7) ! dz^2
                    soe(4) = s3*M*N*(N*N-(L*L+M*M)/2)*(-1)
                    soe(5) = s3*M*N*(L*L+M*M-N*N)*(-1)
                    soe(6) = s3*M*N*(L*L+M*M)/2*(-1)
                case(8) ! dxz
                    soe(4) = 3*L*M*N*N
                    soe(5) = L*M*(1-4*N*N)
                    soe(6) = L*M*(N*N-1)
                case(9) ! dx2-y2
                    soe(4) = 1.5*M*N*(L*L-M*M)*(-1)
                    soe(5) = -M*N*(1+2*(L*L-M*M))*(-1)
                    soe(6) = M*N*(1+0.5*(L*L-M*M))*(-1)
                end select
            elseif (ind_store(kk, 2) == 7) then ! dz^2
                select case(ind_store(kk, 3))
                case(5) ! dxy
                    soe(4) = s3*L*M*(N*N-0.5*(L*L+M*M))
                    soe(5) = -2*s3*L*M*N*N
                    soe(6) = (s3/2)*L*M*(1+N*N)
                case(6) ! dyz
                    soe(4) = s3*M*N*(N*N-(L*L+M*M)/2)*(-1)
                    soe(5) = s3*M*N*(L*L+M*M-N*N)*(-1)
                    soe(6) = s3*M*N*(L*L+M*M)/2*(-1)
                case(7) ! dz^2
                    soe(4) = (N*N-0.5*(L*L+M*M))**2
                    soe(5) = 3*N*N*(L*L+M*M)
                    soe(6) = 0.75*(L*L+M*M)**2
                case(8) ! dxz
                    soe(4) = s3*L*N*(N*N-0.5*(L*L+M*M))*(-1)
                    soe(5) = s3*L*N*(L*L+M*M-N*N)*(-1)
                    soe(6) = -(s3/2)*L*N*(L*L+M*M)*(-1)
                case(9) ! dx2-y2
                    soe(4) = (s3/2)*(L*L-M*M)*(N*N-0.5*(L*L+M*M))
                    soe(5) = s3*N*N*(M*M-L*L)
                    soe(6) = (s3/4)*(1+N*N)*(L*L-M*M)
                end select
            elseif (ind_store(kk, 2) == 8) then ! dxz
                select case(ind_store(kk, 3))
                case(5) ! dxy
                    soe(4) = 3*L*L*M*N*(-1)
                    soe(5) = M*N*(1 - 4*L*L)*(-1)
                    soe(6) = M*N*(L*L - 1)*(-1)
                case(6) ! dyz
                    soe(4) = 3*L*M*N*N
                    soe(5) = L*M*(1-4*N*N)
                    soe(6) = L*M*(N*N-1)
                case(7) ! dz^2
                    soe(4) = s3*L*N*(N*N-0.5*(L*L+M*M))*(-1)
                    soe(5) = s3*L*N*(L*L+M*M-N*N)*(-1)
                    soe(6) = -(s3/2)*L*N*(L*L+M*M)*(-1)
                case(8) ! dxz
                    soe(4) = 3*L*L*N*N
                    soe(5) = L*L + N*N - 4*L*L*N*N
                    soe(6) = M*M + L*L*N*N
                case(9) ! dx2-y2
                    soe(4) = 1.5*N*L*(L*L-M*M)*(-1)
                    soe(5) = N*L*(1-2*(L*L-M*M))*(-1)
                    soe(6) = -N*L*(1-0.5*(L*L-M*M))*(-1)
                end select
            elseif (ind_store(kk, 2) == 9) then ! dx2-y2
                select case(ind_store(kk, 3))
                case(5) ! dxy
                    soe(4) = 3*L*M*(L*L - M*M)/2
                    soe(5) = 2*L*M*(M*M - L*L)
                    soe(6) = L*M*(L*L - M*M)/2
                case(6) ! dyz
                    soe(4) = 1.5*M*N*(L*L-M*M)*(-1)
                    soe(5) = -M*N*(1+2*(L*L-M*M))*(-1)
                    soe(6) = M*N*(1+0.5*(L*L-M*M))*(-1)
                case(7) ! dz^2
                    soe(4) = (s3/2)*(L*L-M*M)*(N*N-0.5*(L*L+M*M))
                    soe(5) = s3*N*N*(M*M-L*L)
                    soe(6) = (s3/4)*(1+N*N)*(L*L-M*M)
                case(8) ! dxz
                    soe(4) = 1.5*N*L*(L*L-M*M)*(-1)
                    soe(5) = N*L*(1-2*(L*L-M*M))*(-1)
                    soe(6) = -N*L*(1-0.5*(L*L-M*M))*(-1)
                case(9) ! dx2-y2
                    soe(4) = 0.75*(L*L-M*M)**2
                    soe(5) = (L*L+M*M - (L*L-M*M)**2)
                    soe(6) = (N*N+0.25*(L*L- M*M)**2)
                end select 
            else
                write(*, '(a35)') 'Error: wrong orbital in dd system 2'
                call exit(0)
            endif

        do mm = kk+1, LL
            distsU(numInters+1, 1) = dists(ind_store(mm, 1))
            distsU(numInters+1, 2:5) = ind_store(mm, 4:7)
        ! runs if third row matches correct interaction
        if (norm2(distsU(numInters+1, :)-distsU(ii, :)) < 1e-4) then

            L = lmn(ind_store(mm, 1), 1)
            M = lmn(ind_store(mm, 1), 2) ! redefining these saves space
            N = lmn(ind_store(mm, 1), 3)

            if (ind_store(mm, 2) == 5) then ! dxy
                select case(ind_store(mm, 3))
                case(5) ! dxy
                    soe(7) = 3*L*L*M*M
                    soe(8) = L*L + M*M - 4*L*L*M*M
                    soe(9) = N*N + L*L*M*M
                case(6) ! dyz
                    soe(7) = 3*L*M*M*N*(-1)
                    soe(8) = L*N*(1 - 4*M*M)*(-1)
                    soe(9) = L*N*(M*M - 1)*(-1)
                case(7) ! dz^2
                    soe(7) = s3*L*M*(N*N-(L*L+M*M)/2)
                    soe(8) = -2*s3*L*M*N*N
                    soe(9) = s3*L*M*(1+N*N)/2
                case(8) ! dxz
                    soe(7) = 3*L*L*M*N*(-1)
                    soe(8) = M*N*(1 - 4*L*L)*(-1)
                    soe(9) = M*N*(L*L - 1)*(-1)
                case(9) ! dx2-y2
                    soe(7) = 3*L*M*(L*L - M*M)/2
                    soe(8) = 2*L*M*(M*M - L*L)
                    soe(9) = L*M*(L*L - M*M)/2
                end select
            elseif (ind_store(mm, 2) == 6) then ! dyz
                select case(ind_store(mm, 3))
                case(5) ! dxy
                    soe(7) = 3*L*M*M*N*(-1)
                    soe(8) = L*N*(1 - 4*M*M)*(-1)
                    soe(9) = L*N*(M*M - 1)*(-1) 
                case(6) ! dyz
                    soe(7) = 3*M*M*N*N
                    soe(8) = M*M + N*N - 4*M*M*N*N
                    soe(9) = L*L + M*M*N*N
                case(7) ! dz^2
                    soe(7) = s3*M*N*(N*N-(L*L+M*M)/2)*(-1)
                    soe(8) = s3*M*N*(L*L+M*M-N*N)*(-1)
                    soe(9) = s3*M*N*(L*L+M*M)/2*(-1)
                case(8) ! dxz
                    soe(7) = 3*L*M*N*N
                    soe(8) = L*M*(1-4*N*N)
                    soe(9) = L*M*(N*N-1)
                case(9) ! dx2-y2
                    soe(7) = 1.5*M*N*(L*L-M*M)*(-1)
                    soe(8) = -M*N*(1+2*(L*L-M*M))*(-1)
                    soe(9) = M*N*(1+0.5*(L*L-M*M))*(-1)
                end select
            elseif (ind_store(mm, 2) == 7) then ! dz^2
                select case(ind_store(mm, 3))
                case(5) ! dxy
                    soe(7) = s3*L*M*(N*N-0.5*(L*L+M*M))
                    soe(8) = -2*s3*L*M*N*N
                    soe(9) = (s3/2)*L*M*(1+N*N)
                case(6) ! dyz
                    soe(7) = s3*M*N*(N*N-(L*L+M*M)/2)*(-1)
                    soe(8) = s3*M*N*(L*L+M*M-N*N)*(-1)
                    soe(9) = s3*M*N*(L*L+M*M)/2*(-1)
                case(7) ! dz^2
                    soe(7) = (N*N-0.5*(L*L+M*M))**2
                    soe(8) = 3*N*N*(L*L+M*M)
                    soe(9) = 0.75*(L*L+M*M)**2
                case(8) ! dxz
                    soe(7) = s3*L*N*(N*N-0.5*(L*L+M*M))*(-1)
                    soe(8) = s3*L*N*(L*L+M*M-N*N)*(-1)
                    soe(9) = -(s3/2)*L*N*(L*L+M*M)*(-1)
                case(9) ! dx2-y2
                    soe(7) = (s3/2)*(L*L-M*M)*(N*N-0.5*(L*L+M*M))
                    soe(8) = s3*N*N*(M*M-L*L)
                    soe(9) = (s3/4)*(1+N*N)*(L*L-M*M)
                end select
            elseif (ind_store(mm, 2) == 8) then ! dxz
                select case(ind_store(mm, 3))
                case(5) ! dxy
                    soe(7) = 3*L*L*M*N*(-1)
                    soe(8) = M*N*(1 - 4*L*L)*(-1)
                    soe(9) = M*N*(L*L - 1)*(-1)
                case(6) ! dyz
                    soe(7) = 3*L*M*N*N
                    soe(8) = L*M*(1-4*N*N)
                    soe(9) = L*M*(N*N-1)
                case(7) ! dz^2
                    soe(7) = s3*L*N*(N*N-0.5*(L*L+M*M))*(-1)
                    soe(8) = s3*L*N*(L*L+M*M-N*N)*(-1)
                    soe(9) = -(s3/2)*L*N*(L*L+M*M)*(-1)
                case(8) ! dxz
                    soe(7) = 3*L*L*N*N
                    soe(8) = L*L + N*N - 4*L*L*N*N
                    soe(9) = M*M + L*L*N*N
                case(9) ! dx2-y2
                    soe(7) = 1.5*N*L*(L*L-M*M)*(-1)
                    soe(8) = N*L*(1-2*(L*L-M*M))*(-1)
                    soe(9) = -N*L*(1-0.5*(L*L-M*M))*(-1)
                end select
            elseif (ind_store(mm, 2) == 9) then ! dx2-y2
                select case(ind_store(mm, 3))
                case(5) ! dxy
                    soe(7) = 3*L*M*(L*L - M*M)/2
                    soe(8) = 2*L*M*(M*M - L*L)
                    soe(9) = L*M*(L*L - M*M)/2
                case(6) ! dyz
                    soe(7) = 1.5*M*N*(L*L-M*M)*(-1)
                    soe(8) = -M*N*(1+2*(L*L-M*M))*(-1)
                    soe(9) = M*N*(1+0.5*(L*L-M*M))*(-1)
                case(7) ! dz^2
                    soe(7) = (s3/2)*(L*L-M*M)*(N*N-0.5*(L*L+M*M))
                    soe(8) = s3*N*N*(M*M-L*L)
                    soe(9) = (s3/4)*(1+N*N)*(L*L-M*M)
                case(8) ! dxz
                    soe(7) = 1.5*N*L*(L*L-M*M)*(-1)
                    soe(8) = N*L*(1-2*(L*L-M*M))*(-1)
                    soe(9) = -N*L*(1-0.5*(L*L-M*M))*(-1)
                case(9) ! dx2-y2
                    soe(7) = 0.75*(L*L-M*M)**2
                    soe(8) = (L*L+M*M - (L*L-M*M)**2)
                    soe(9) = (N*N+0.25*(L*L- M*M)**2)
                end select
            else
                write(*, '(a35)') 'Error: wrong orbital in dd system 3'
                call exit(0)
            endif

            ! finds determinant for system of equations
            L =  soe(1)*(soe(5)*soe(9)-soe(6)*soe(8)) &
                -soe(2)*(soe(4)*soe(9)-soe(6)*soe(7)) &
                +soe(3)*(soe(4)*soe(8)-soe(5)*soe(7))
            if (abs(L) < 1e-2) then 
                cycle ! skips if no solution (singular system)
            endif
            
            ! We now have a guaranteed solution
            t = t + 1 

            ! Solves 3x3 system of equations analytically
            VddSig(ii, 1) = VddSig(ii, 1) + (1/L)*(&
                (soe(5)*soe(9)-soe(6)*soe(8))*hsparse(ind_store(jj, 1))-&
                (soe(2)*soe(9)-soe(3)*soe(8))*hsparse(ind_store(kk, 1))+&
                (soe(2)*soe(6)-soe(3)*soe(5))*hsparse(ind_store(mm, 1)))
            VddSig(ii, 2) = VddSig(ii, 2) + (1/L)*(&
                (soe(5)*soe(9)-soe(6)*soe(8))*ssparse(ind_store(jj, 1))-&
                (soe(2)*soe(9)-soe(3)*soe(8))*ssparse(ind_store(kk, 1))+&
                (soe(2)*soe(6)-soe(3)*soe(5))*ssparse(ind_store(mm, 1)))
            VddPi(ii, 1) = VddPi(ii, 1) + (1/L)*(-&
                (soe(4)*soe(9)-soe(6)*soe(7))*hsparse(ind_store(jj, 1))+&
                (soe(1)*soe(9)-soe(3)*soe(7))*hsparse(ind_store(kk, 1))-&
                (soe(1)*soe(6)-soe(3)*soe(4))*hsparse(ind_store(mm, 1)))
            VddPi(ii, 2) = VddPi(ii, 2) + (1/L)*(-&
                (soe(4)*soe(9)-soe(6)*soe(7))*ssparse(ind_store(jj, 1))+&
                (soe(1)*soe(9)-soe(3)*soe(7))*ssparse(ind_store(kk, 1))-&
                (soe(1)*soe(6)-soe(3)*soe(4))*ssparse(ind_store(mm, 1)))
            VddDelt(ii, 1) = VddDelt(ii, 1) + (1/L)*(&
                (soe(4)*soe(8)-soe(5)*soe(7))*hsparse(ind_store(jj, 1))-&
                (soe(1)*soe(8)-soe(2)*soe(7))*hsparse(ind_store(kk, 1))+&
                (soe(1)*soe(5)-soe(2)*soe(4))*hsparse(ind_store(mm, 1)))
            VddDelt(ii, 2) = VddDelt(ii, 2) + (1/L)*(&
                (soe(4)*soe(8)-soe(5)*soe(7))*ssparse(ind_store(jj, 1))-&
                (soe(1)*soe(8)-soe(2)*soe(7))*ssparse(ind_store(kk, 1))+&
                (soe(1)*soe(5)-soe(2)*soe(4))*ssparse(ind_store(mm, 1)))

        endif ! third row
        enddo ! mm
        endif ! second row
        enddo ! kk
        endif ! first row
        enddo ! jj
        if (t == 0) then ! if no interaction is calculated
            t = 1
        endif

        VddSig(ii, :) = VddSig(ii, :)/t ! takes average
        VddPi(ii, :) = VddPi(ii, :)/t ! takes average
        VddDelt(ii, :) = VddDelt(ii, :)/t ! takes average
    enddo ! ii
    endif

    !------------------------------------------------------!
    !---------- Outputs tight binding parameters ----------!
    !------------------------------------------------------!

    inquire(file='TB_PARAMS', exist=file_exists) ! checks if it exists

    if (file_exists .eqv. .True.) then
            file_status = 'old' ! if it exists, then override it
    else
            file_status = 'new' ! if it does not exist, create it
    end if

    open(unit=2, file='TB_PARAMS', status=file_status)
    write(2, *)
    write(2, '(a31)') 'Distance (angstroms) and H (eV)' ! writes units
    write(2, *)
    write(2, '(i3, a41)') numInters, &
                    ' ! number of unique distances/interactions'
    write(2, '(i3, a37)') no_u, ' ! number of atomic orbitals in system'
    write(2, *)

    ! outputs hamiltonian section headers
    LL = 0 ! will count number of columns
    write(2, '(a12)', advance='no') '   Distance '
    if (has_s .eqv. .TRUE.) then
        LL = LL + 1
        write(2, '(a12)', advance='no') '      Hss   '
        if (has_p .eqv. .TRUE.) then
            LL = LL + 2
            write(2, '(a12)', advance='no') '      Hsp   '
            write(2, '(a12)', advance='no') '      Hps   '
        endif
    endif
    if (has_p .eqv. .TRUE.) then
        LL = LL + 2
        write(2, '(a12)', advance='no') '    HppSig  '
        write(2, '(a12)', advance='no') '     HppPi  '
    endif
    write(2, '(a12)', advance='no') '   Species  '
    write(2, '(a5)') 'Zetas'
    do ii = 1, LL+2
        write(2, '(a12)', advance='no') '------------'
    enddo ! ii
    write(2, '(a6)') '------'

    ! outputs hamiltonian tight-binding parameters
    do ii = 1, numInters
        write(2, '(f12.6)', advance='no') distsU(ii, 1)
        if (has_s .eqv. .TRUE.) then ! For s orbitals
            write(2, '(f12.6)', advance='no') Vss(ii, 1) ! Hss
            if (has_p .eqv. .TRUE.) then ! For s and p orbitals
                write(2, '(2f12.6)', advance='no') Vsp(ii, 1), Vps(ii, 1)
            endif
        endif
        if (has_p .eqv. .TRUE.) then ! For p orbitals
            write(2, '(2f12.6)', advance='no') VppSig(ii, 1), VppPi(ii, 1)
        endif
        write(2, '(a4, 2i2, a4)', advance='no') '    ', &
                                            nint(distsU(ii, 2:3)), '    '
        write(2, '(2i2)', advance='no') nint(distsU(ii, 4:5))
        write(2, *) ! end of line
    enddo ! ii
    write(2, *)

    ! outputs hamiltonian section headers
    LL = 0 ! will count number of columns
    if (has_d .eqv. .TRUE.) then
        write(2, '(a12)', advance='no') '   Distance '
        if (has_s .eqv. .TRUE.) then
            LL = LL + 1
            write(2, '(a12)', advance='no') '      Hsd   '
            write(2, '(a12)', advance='no') '      Hds   '
        endif
        if (has_p .eqv. .TRUE.) then
            LL = LL + 4
            write(2, '(a12)', advance='no') '    HpdSig  '
            write(2, '(a12)', advance='no') '     HpdPi  '
            write(2, '(a12)', advance='no') '    HdpSig  '
            write(2, '(a12)', advance='no') '     HdpPi  '
        endif
        LL = LL + 3
        write(2, '(a12)', advance='no') '    HddSig  '
        write(2, '(a12)', advance='no') '     HddPi  '
        write(2, '(a12)', advance='no') '    HddDelt '
        write(2, '(a12)', advance='no') '   Species  '
        write(2, '(a5)') 'Zetas'
        do ii = 1, LL+3
            write(2, '(a12)', advance='no') '------------'
        enddo ! ii
        write(2, '(a6)') '------'

        ! outputs hamiltonian tight-binding parameters
        do ii = 1, numInters
            write(2, '(f12.6)', advance='no') distsU(ii, 1)
            if (has_s .eqv. .TRUE.) then ! For s orbitals
                write(2, '(f12.6)', advance='no') Vsd(ii, 1) ! Hsd
                write(2, '(f12.6)', advance='no') Vds(ii, 1) ! Hds
            endif
            if (has_p .eqv. .TRUE.) then ! For p orbitals
                write(2, '(2f12.6)', advance='no')VpdSig(ii, 1), VpdPi(ii, 1)
                write(2, '(2f12.6)', advance='no')VdpSig(ii, 1), VdpPi(ii, 1)
            endif
            write(2, '(2f12.6)', advance='no') VddSig(ii, 1), VddPi(ii, 1)
            write(2, '(f12.6)', advance='no') VddDelt(ii, 1)
            write(2, '(a4, 2i2, a4)', advance='no') '    ', &
                                            nint(distsU(ii, 2:3)), '    '
            write(2, '(2i2)', advance='no') nint(distsU(ii, 4:5))
            write(2, *) ! end of line
        enddo ! ii
        write(2, *)
    endif

    ! outputs overlap section headers
    LL = 0 ! will count number of columns
    write(2, '(a12)', advance='no') '   Distance '
    if (has_s .eqv. .TRUE.) then
        LL = LL + 1
        write(2, '(a12)', advance='no') '      Sss   '
        if (has_p .eqv. .TRUE.) then
            LL = LL + 2
            write(2, '(a12)', advance='no') '      Ssp   '
            write(2, '(a12)', advance='no') '      Sps   '
        endif
    endif
    if (has_p .eqv. .TRUE.) then
        LL = LL + 2
        write(2, '(a12)', advance='no') '    SppSig  '
        write(2, '(a12)', advance='no') '     SppPi  '
    endif
    write(2, '(a12)', advance='no') '   Species  '
    write(2, '(a5)') 'Zetas'
    do ii = 1, LL+2
        write(2, '(a12)', advance='no') '------------'
    enddo ! ii
    write(2, '(a6)') '------'

    ! outputs overlap tight-binding parameters
    do ii = 1, numInters
        write(2, '(f12.6)', advance='no') distsU(ii, 1)
        if (has_s .eqv. .TRUE.) then ! For s orbitals
            write(2, '(f12.6)', advance='no') Vss(ii, 2) ! Hss
            if (has_p .eqv. .TRUE.) then ! For s and p orbitals
                write(2, '(2f12.6)', advance='no') Vsp(ii, 2), Vps(ii, 2)
            endif
        endif
        if (has_p .eqv. .TRUE.) then ! For p orbitals
            write(2, '(2f12.6)', advance='no') VppSig(ii, 2), VppPi(ii, 2)
        endif
        write(2, '(a4, 2i2, a4)', advance='no') '    ', &
                                            nint(distsU(ii, 2:3)), '    '
        write(2, '(2i2)', advance='no') nint(distsU(ii, 4:5))
        write(2, *) ! end of line
    enddo ! ii
    write(2, *)

    ! outputs overlap section headers
    LL = 0 ! will count number of columns
    if (has_d .eqv. .TRUE.) then
        write(2, '(a12)', advance='no') '   Distance '
        if (has_s .eqv. .TRUE.) then
            LL = LL + 1
            write(2, '(a12)', advance='no') '      Ssd   '
            write(2, '(a12)', advance='no') '      Sds   '
        endif
        if (has_p .eqv. .TRUE.) then
            LL = LL + 4
            write(2, '(a12)', advance='no') '    SpdSig  '
            write(2, '(a12)', advance='no') '     SpdPi  '
            write(2, '(a12)', advance='no') '    SdpSig  '
            write(2, '(a12)', advance='no') '     SdpPi  '
        endif
        LL = LL + 3
        write(2, '(a12)', advance='no') '    SddSig  '
        write(2, '(a12)', advance='no') '     SddPi  '
        write(2, '(a12)', advance='no') '    SddDelt '
        write(2, '(a12)', advance='no') '   Species  '
        write(2, '(a5)') 'Zetas'
        do ii = 1, LL+3
            write(2, '(a12)', advance='no') '------------'
        enddo ! ii
        write(2, '(a6)') '------'

        ! outputs overlap tight-binding parameters
        do ii = 1, numInters
            write(2, '(f12.6)', advance='no') distsU(ii, 1)
            if (has_s .eqv. .TRUE.) then ! For s orbitals
                write(2, '(f12.6)', advance='no') Vsd(ii, 2)
                write(2, '(f12.6)', advance='no') Vds(ii, 2)
            endif
            if (has_p .eqv. .TRUE.) then ! For p orbitals
                write(2, '(2f12.6)', advance='no')VpdSig(ii, 2), VpdPi(ii, 2)
                write(2, '(2f12.6)', advance='no')VdpSig(ii, 2), VdpPi(ii, 2)
            endif
            write(2, '(3f12.6)', advance='no') VddSig(ii, 2), VddPi(ii, 2)
            write(2, '(f12.6)', advance='no') VddDelt(ii, 2)
            write(2, '(a4, 2i2, a4)', advance='no') '    ', &
                                            nint(distsU(ii, 2:3)), '    '
            write(2, '(2i2)', advance='no') nint(distsU(ii, 4:5))
            write(2, *) ! end of line
        enddo ! ii
        write(2, *)
    endif

    ! outputs onsite matrices
    write(2, '(a25)') 'Onsite hamiltonian matrix'
    do ii = 1, no_u
        do jj = 1, no_u
            write(2, '(f12.6)', advance='no') Vonsite(ii, jj, 1)
        enddo ! jj
        write(2, *)
    enddo ! ii
    write(2, *)
    write(2, '(a21)') 'Onsite overlap matrix'
    do ii = 1, no_u
        do jj = 1, no_u
            write(2, '(f12.6)', advance='no') Vonsite(ii, jj, 2)
        enddo ! jj
        write(2, *)
    enddo ! ii
    write(2, *)

    close(2)
    
    !------------------------------------------------------!
    !--------- Inputs lattice vectors and kpoints ---------!
    !------------------------------------------------------!
    ! Inputs lattice vectors
    inquire(file='input.fdf', exist=file_exists) ! checks if it exists

    if (file_exists .eqv. .FALSE.) then
        write(*, '(a40)') 'Error: input.fdf file does not exist'
        call exit(0)
    end if

    tc = .FALSE.
    open(unit=3, file='input.fdf', status='old')
    do while (tc .eqv. .FALSE.)
        read(3, '(A)', end=99) line
        if (trim(line)=='%block LatticeVectors') then
            read(3, *) A(:, 1)
            read(3, *) A(:, 2)
            read(3, *) A(:, 3)
            tc = .TRUE.
        endif
    enddo
    99 continue
    close(3)
    if (tc .eqv. .FALSE.) then
        write(*, '(a45)') 'Error: lattice vectors not found in input.fdf'
    endif

    pi = 4*atan(1.0)
    ! calculates unit cell volume
    vs(1) = A(2, 2)*A(3, 3) - A(3, 2)*A(2, 3)
    vs(2) = A(3, 2)*A(1, 3) - A(3, 3)*A(1, 2)
    vs(3) = A(1, 2)*A(2, 3) - A(2, 2)*A(1, 3)
    L = sum(A(:, 1)*vs(:)) ! stores unit cell volume

    ! makes each reciprocal lattice vector
    B(:, 1) = vs(:)
    B(1, 2) = A(2, 3)*A(3, 1) - A(3, 3)*A(2, 1)
    B(2, 2) = A(3, 3)*A(1, 1) - A(1, 3)*A(3, 1)
    B(3, 2) = A(1, 3)*A(2, 1) - A(2, 3)*A(1, 1)
    B(1, 3) = A(2, 1)*A(3, 2) - A(3, 1)*A(2, 2)
    B(2, 3) = A(3, 1)*A(1, 2) - A(1, 1)*A(3, 2)
    B(3, 3) = A(1, 1)*A(2, 2) - A(2, 1)*A(1, 2)
    B = B*2*pi/L

    ! inputs kpoints
    inquire(file='kpoints.dat', exist=file_exists) ! checks if it exists

    if (file_exists .eqv. .FALSE.) then
        write(*, '(a38)') 'Error: kpoints.dat file does not exist'
        call exit(0)
    end if

    tc = .FALSE.
    open(unit=4, file='kpoints.dat', status='old')
    read(4, *) LL ! inputs number of kpoint-paths within file
    allocate(kden(LL), kpE(2*LL, 3), mags(LL))
    allocate(xlab(2*LL), xlabF(LL+1, 3), xlabInd(LL+1))

    do ii=1, LL
        read(4, *) xlab(2*ii-1), kpE(2*ii-1, :), xlab(2*ii), &
                    kpE(2*ii, :), kden(ii)
    enddo ! ii
    close(4)

    ! Checks for discontinuities within k-point path
    xlabF(1, :) = (/xlab(1), '', ''/)
    xlabInd(1) = 1
    do ii=1, LL-1
        if (xlab(2*ii)==xlab(2*ii+1)) then
            xlabF(ii+1, :)=(/xlab(2*ii), '', ''/)
            xlabInd(ii+1) = 1
        else
            xlabF(ii+1, :) = (/ xlab(2*ii), '|', xlab(2*ii+1) /)
            xlabInd(ii+1) = 3
        endif
    enddo
    xlabF(LL+1, :) = (/xlab(2*LL), '', ''/)
    xlabInd(LL+1) = 1

    allocate(kp(sum(kden)+LL, 3)) ! (1, 2, 3) is (x, y, z) coords.
    allocate(xout(sum(kden)+LL))

    ! makes kpoint paths
    do ii=1, LL ! iterates through each path
        vs = matmul(B(:, :), kpE(2*ii-1, :)) ! start coord in cartesian
        ve = matmul(B(:, :), kpE(2*ii, :)) ! end coord in cartesian
        mags(ii) = norm2(ve-vs) ! magnitude of k-point path
        do jj=1, kden(ii)+1 ! iterates over k-point density
            kp(jj+sum(kden(1:ii-1))+ii-1, 1) = vs(1)+((jj-1.0)/kden(ii))&
                                                *(ve(1)-vs(1))
            kp(jj+sum(kden(1:ii-1))+ii-1, 2) = vs(2)+((jj-1.0)/kden(ii))&
                                                *(ve(2)-vs(2))
            kp(jj+sum(kden(1:ii-1))+ii-1, 3) = vs(3)+((jj-1.0)/kden(ii))&
                                                *(ve(3)-vs(3))
            xout(jj+sum(kden(1:ii-1))+ii-1) = sum(mags(1:ii-1))+&
                                            ((jj-1.0)/kden(ii))*mags(ii)
        enddo ! jj
    enddo ! ii

    !------------------------------------------------------!
    !--------- Makes hamiltonian/overlap matrices ---------!
    !------------------------------------------------------!

    e = complex(0, 1) ! defines sqrt(-1)
    allocate(H(no_u, no_u), S(no_u, no_u)) ! creates hamiltonian/overlap
    allocate(energy(no_u), work(2*no_u-1), rwork(3*no_u-2)) ! for eigs

    ! opens output file for tight-binding bands
    inquire(file='TB_EIGS', exist=file_exists) ! checks if it exists

    if (file_exists .eqv. .True.) then
            file_status = 'old' ! if it exists, then override it
    else
            file_status = 'new' ! if it does not exist, create it
    end if

    open(unit=5, file='TB_EIGS', status=file_status)
    write(5, '(i4, a29)') sum(kden)+LL, ' ! number of kpoints evaluated'
    write(5, '(i4, a20)') no_u, ' ! number of orbitals'
    write(5, *)

    do t2 = 1, sum(kden)+LL ! cycles through k-points
    vs(:) = (/ kp(t2, 1), kp(t2, 2), kp(t2, 3) /) ! defines current k-point 

    ! resets hamiltonian/overlap matrices for each k-point
    H = Vonsite(:, :, 1)
    S = Vonsite(:, :, 2)

    do ii = 1, no_u ! cycles through unit cell orbitals ! row
    do t = 1, numh(ii)
        kk = listhptr(ii)+t ! cycles through all sparse elements
        jj = indxuo(listh(kk)) ! cycles through unit cell orbs ! col

        ! skips onsite terms and lower triangle
        if (dists(kk) < 1e-4 .or. jj < ii) then
            cycle
        endif

        phase = exp(e*sum(vs(:)*xij(kk, :)))
        L = lmn(kk, 1) ! redefining these significantly reduces
        M = lmn(kk, 2) ! space needed for coding next part
        N = lmn(kk, 3)

        ! finds index for current interaction
        distsU(numInters+1, 1) = dists(kk) ! distance
        distsU(numInters+1, 2:3) = (/indx(ii, 1), indx(jj, 1)/) ! species
        distsU(numInters+1, 4:5) = (/indx(ii, 3), indx(jj, 3)/) ! zetas
        do t3 = 1, numInters
            if ((norm2(distsU(numInters+1, :)-distsU(t3, :)) < 1e-4)) then
                exit ! stores index within t3 variable
            endif
        enddo ! t3

        if (indx(ii, 2) == 1) then
            select case(indx(jj, 2))
                case(1) ! s-s
                    eH = Vss(t3, 1)
                    eS = Vss(t3, 2)
                case(2) ! s-py
                    eH = M*Vsp(t3, 1)
                    eS = M*Vsp(t3, 2)
                case(3) ! s-pz
                    eH = N*Vsp(t3, 1)*(-1)
                    eS = N*Vsp(t3, 2)*(-1)
                case(4) ! s-px
                    eH = L*Vsp(t3, 1)
                    eS = L*Vsp(t3, 2)
                case(5) ! s-dxy
                    eH = s3*L*M*Vsd(t3, 1)
                    eS = s3*L*M*Vsd(t3, 2)
                case(6) ! s-dyz
                    eH = s3*M*N*Vsd(t3, 1)*(-1)
                    eS = s3*M*N*Vsd(t3, 2)*(-1)
                case(7) ! s-dz^2r^2
                    eH = (N*N-0.5*(L*L+M*M))*Vsd(t3, 1)
                    eS = (N*N-0.5*(L*L+M*M))*Vsd(t3, 2)
                case(8) ! s-dxz
                    eH = s3*L*N*Vsd(t3, 1)*(-1)
                    eS = s3*L*N*Vsd(t3, 2)*(-1)
                case(9) ! s-dx^2y^2
                    eH = (s3/2)*(L*L-M*M)*Vsd(t3, 1)
                    eS = (s3/2)*(L*L-M*M)*Vsd(t3, 2)
            end select
        elseif (indx(ii, 2) == 2) then
            select case(indx(jj, 2))
                case(1) ! py-s
                    eH = M*Vps(t3, 1)
                    eS = M*Vps(t3, 2)
                case(2) ! py-py
                    eH = M*M*VppSig(t3, 1) + (1-M*M)*VppPi(t3, 1)
                    eS = M*M*VppSig(t3, 2) + (1-M*M)*VppPi(t3, 2)
                case(3) ! py-pz
                    eH = M*N*(VppSig(t3, 1)-VppPi(t3, 1))*(-1)
                    eS = M*N*(VppSig(t3, 2)-VppPi(t3, 2))*(-1)
                case(4) ! py-px
                    eH = M*L*(VppSig(t3, 1)-VppPi(t3, 1))
                    eS = M*L*(VppSig(t3, 2)-VppPi(t3, 2))
                case(5) ! py-dxy
                    eH = s3*M*M*L*VpdSig(t3, 1)+L*(1-2*M*M)*VpdPi(t3, 1)
                    eS = s3*M*M*L*VpdSig(t3, 2)+L*(1-2*M*M)*VpdPi(t3, 2)
                case(6) ! py-dyz
                    eH = (s3*M*M*N*VpdSig(t3, 1)+N*(1-2*M*M)*&
                         VpdPi(t3, 1))*(-1)
                    eS = (s3*M*M*N*VpdSig(t3, 2)+N*(1-2*M*M)*&
                         VpdPi(t3, 2))*(-1)
                case(7) ! py-dz^2r^2
                    eH = M*(N*N-0.5*(L*L+M*M))*VpdSig(t3, 1) - s3*M*N*N*&
                         VpdPi(t3, 1)
                    eS = M*(N*N-0.5*(L*L+M*M))*VpdSig(t3, 2) - s3*M*N*N*&
                         VpdPi(t3, 2)
                case(8) ! py-dxz
                    eH = (s3*L*M*N*VpdSig(t3, 1)-2*L*M*N*VpdPi(t3, 1))&
                         *(-1)
                    eS = (s3*L*M*N*VpdSig(t3, 2)-2*L*M*N*VpdPi(t3, 2))&
                         *(-1)
                case(9) ! py-dx^2y^2
                    eH = (s3/2)*M*(L*L-M*M)*VpdSig(t3, 1) - M*(1+L*L-&
                         M*M)*VpdPi(t3, 1)
                    eS = (s3/2)*M*(L*L-M*M)*VpdSig(t3, 2) - M*(1+L*L-&
                         M*M)*VpdPi(t3, 2)
            end select
        elseif (indx(ii, 2) == 3) then
            select case(indx(jj, 2))
                case(1) ! pz-s
                    eH = N*Vps(t3, 1)*(-1)
                    eS = N*Vps(t3, 2)*(-1)
                case(2) ! pz-py
                    eH = N*M*(VppSig(t3, 1)-VppPi(t3, 1))*(-1)
                    eS = N*M*(VppSig(t3, 2)-VppPi(t3, 2))*(-1)
                case(3) ! pz-pz
                    eH = N*N*VppSig(t3, 1) + (1-N*N)*VppPi(t3, 1)
                    eS = N*N*VppSig(t3, 2) + (1-N*N)*VppPi(t3, 2)
                case(4) ! pz-px
                    eH = L*N*(VppSig(t3, 1)-VppPi(t3, 1))*(-1)
                    eS = L*N*(VppSig(t3, 2)-VppPi(t3, 2))*(-1)
                case(5) ! pz-dxy
                    eH = (s3*L*M*N*VpdSig(t3, 1)-2*L*M*N*VpdPi(t3, 1))&
                         *(-1)
                    eS = (s3*L*M*N*VpdSig(t3, 2)-2*L*M*N*VpdPi(t3, 2))&
                         *(-1)
                case(6) ! pz-dyz
                    eH = s3*N*N*M*VpdSig(t3, 1)+M*(1-2*N*N)*VpdPi(t3, 1)
                    eS = s3*N*N*M*VpdSig(t3, 2)+M*(1-2*N*N)*VpdPi(t3, 2)
                case(7) ! pz-dz^2r^2
                    eH = (N*(N*N-0.5*(L*L+M*M))*VpdSig(t3, 1)+s3*N*&
                          (L*L+M*M)*VpdPi(t3, 1))*(-1)
                    eS = (N*(N*N-0.5*(L*L+M*M))*VpdSig(t3, 2)+s3*N*&
                          (L*L+M*M)*VpdPi(t3, 2))*(-1)
                case(8) ! pz-dxz
                    eH = s3*N*N*L*VpdSig(t3, 1)+L*(1-2*N*N)*VpdPi(t3, 1)
                    eS = s3*N*N*L*VpdSig(t3, 2)+L*(1-2*N*N)*VpdPi(t3, 2)
                case(9) ! pz-dx^2y^2
                    eH = ((s3/2)*N*(L*L-M*M)*VpdSig(t3, 1)-N*(L*L-M*M)*&
                         VpdPi(t3, 1))*(-1)
                    eS = ((s3/2)*N*(L*L-M*M)*VpdSig(t3, 2)-N*(L*L-M*M)*&
                         VpdPi(t3, 2))*(-1)
            end select
        elseif (indx(ii, 2) == 4) then
            select case(indx(jj, 2))
                case(1) ! px-s
                    eH = L*Vps(t3, 1)
                    eS = L*Vps(t3, 2)
                case(2) ! px-py
                    eH = L*M*(VppSig(t3, 1)-VppPi(t3, 1))
                    eS = L*M*(VppSig(t3, 2)-VppPi(t3, 2))
                case(3) ! px-pz
                    eH = L*N*(VppSig(t3, 1)-VppPi(t3, 1))*(-1)
                    eS = L*N*(VppSig(t3, 2)-VppPi(t3, 2))*(-1)
                case(4) ! px-px
                    eH = L*L*VppSig(t3, 1) + (1-L*L)*VppPi(t3, 1)
                    eS = L*L*VppSig(t3, 2) + (1-L*L)*VppPi(t3, 2)
                case(5) ! px-dxy
                    eH = s3*L*L*M*VpdSig(t3, 1)+M*(1-2*L*L)*VpdPi(t3, 1)
                    eS = s3*L*L*M*VpdSig(t3, 2)+M*(1-2*L*L)*VpdPi(t3, 2)
                case(6) ! px-dyz
                    eH = (s3*L*M*N*VpdSig(t3, 1)-2*L*M*N*VpdPi(t3, 1))&
                         *(-1)
                    eS = (s3*L*M*N*VpdSig(t3, 2)-2*L*M*N*VpdPi(t3, 2))&
                         *(-1)
                case(7) ! px-dz^2r^2
                    eH = L*(N*N-0.5*(L*L+M*M))*VpdSig(t3, 1)-s3*L*N*N*&
                         VpdPi(t3, 1)
                    eS = L*(N*N-0.5*(L*L+M*M))*VpdSig(t3, 2)-s3*L*N*N*&
                         VpdPi(t3, 2)
                case(8) ! px-dxz
                    eH = (s3*L*L*N*VpdSig(t3, 1)+N*(1-2*L*L)*&
                         VpdPi(t3, 1))*(-1)
                    eS = (s3*L*L*N*VpdSig(t3, 2)+N*(1-2*L*L)*&
                         VpdPi(t3, 2))*(-1)
                case(9) ! px-dx^2y^2
                    eH = (s3/2)*L*(L*L-M*M)*VpdSig(t3, 1)+L*(1-L*L+M*M)*&
                         VpdPi(t3, 1)
                    eS = (s3/2)*L*(L*L-M*M)*VpdSig(t3, 2)+L*(1-L*L+M*M)*&
                         VpdPi(t3, 2)
            end select
        elseif (indx(ii, 2) == 5) then
            select case(indx(jj, 2))
                case(1) ! dxy-s
                    eH = s3*L*M*Vds(t3, 1)
                    eS = s3*L*M*Vds(t3, 2)
                case(2) ! dxy-py
                    eH = s3*M*M*L*VdpSig(t3, 1)+L*(1-2*M*M)*VdpPi(t3, 1)
                    eS = s3*M*M*L*VdpSig(t3, 2)+L*(1-2*M*M)*VdpPi(t3, 2)
                case(3) ! dxy-pz
                    eH = (s3*L*M*N*VdpSig(t3, 1)-2*L*M*N*VdpPi(t3, 1))&
                         *(-1)
                    eS = (s3*L*M*N*VdpSig(t3, 2)-2*L*M*N*VdpPi(t3, 2))&
                         *(-1)
                case(4) ! dxy-px
                    eH = s3*L*L*M*VdpSig(t3, 1)+M*(1-2*L*L)*VdpPi(t3, 1)
                    eS = s3*L*L*M*VdpSig(t3, 2)+M*(1-2*L*L)*VdpPi(t3, 2)
                case(5) ! dxy-dxy
                    eH = 3*L*L*M*M*VddSig(t3, 1)+(L*L+M*M-4*L*L*M*M)*&
                         VddPi(t3, 1) + (N*N+L*L*M*M)*VddDelt(t3, 1)
                    eS = 3*L*L*M*M*VddSig(t3, 2)+(L*L+M*M-4*L*L*M*M)*&
                         VddPi(t3, 2) + (N*N+L*L*M*M)*VddDelt(t3, 2)
                case(6) ! dxy-dyz
                    eH = (3*L*M*M*N*VddSig(t3, 1)+L*N*(1-4*M*M)*&
                         VddPi(t3, 1)+L*N*(M*M-1)*VddDelt(t3, 1))*(-1)
                    eS = (3*L*M*M*N*VddSig(t3, 2)+L*N*(1-4*M*M)*&
                         VddPi(t3, 2)+L*N*(M*M-1)*VddDelt(t3, 2))*(-1)
                case(7) ! dxy-dz^2r^2
                    eH = s3*L*M*(N*N-0.5*(L*L+M*M))*VddSig(t3, 1)-2*s3*&
                         L*M*N*N*VddPi(t3, 1)+(s3/2)*L*M*(1+N*N)*&
                         VddDelt(t3, 1)
                    eS = s3*L*M*(N*N-0.5*(L*L+M*M))*VddSig(t3, 2)-2*s3*&
                         L*M*N*N*VddPi(t3, 2)+(s3/2)*L*M*(1+N*N)*&
                         VddDelt(t3, 2)
                case(8) ! dxy-dxz
                    eH = (3*L*L*M*N*VddSig(t3, 1)+M*N*(1-4*L*L)*&
                         VddPi(t3, 1)+M*N*(L*L-1)*VddDelt(t3, 1))*(-1)
                    eS = (3*L*L*M*N*VddSig(t3, 2)+M*N*(1-4*L*L)*&
                         VddPi(t3, 2)+M*N*(L*L-1)*VddDelt(t3, 2))*(-1)
                case(9) ! dxy-dx^2y^2
                    eH = 1.5*L*M*(L*L-M*M)*VddSig(t3, 1)+2*L*M*(M*M-L*L)&
                         *VddPi(t3, 1)+0.5*L*M*(L*L-M*M)*VddDelt(t3, 1)
                    eS = 1.5*L*M*(L*L-M*M)*VddSig(t3, 2)+2*L*M*(M*M-L*L)&
                         *VddPi(t3, 2)+0.5*L*M*(L*L-M*M)*VddDelt(t3, 2)
            end select
        elseif (indx(ii, 2) == 6) then
            select case(indx(jj, 2))
                case(1) ! dyz-s
                    eH = s3*M*N*Vds(t3, 1)*(-1)
                    eS = s3*M*N*Vds(t3, 2)*(-1)
                case(2) ! dyz-py
                    eH = (s3*M*M*N*VdpSig(t3, 1)+N*(1-2*M*M)*&
                         VdpPi(t3, 1))*(-1)
                    eS = (s3*M*M*N*VdpSig(t3, 2)+N*(1-2*M*M)*&
                         VdpPi(t3, 2))*(-1)
                case(3) ! dyz-pz
                    eH = s3*N*N*M*VdpSig(t3, 1)+M*(1-2*N*N)*VdpPi(t3, 1)
                    eS = s3*N*N*M*VdpSig(t3, 2)+M*(1-2*N*N)*VdpPi(t3, 2)
                case(4) ! dyz-px
                    eH = (s3*L*M*N*VdpSig(t3, 1)-2*L*M*N*VdpPi(t3, 1))*&
                         (-1)
                    eS = (s3*L*M*N*VdpSig(t3, 2)-2*L*M*N*VdpPi(t3, 2))*&
                         (-1)
                case(5) ! dyz-dxy
                    eH = (3*L*M*M*N*VddSig(t3, 1)+L*N*(1-4*M*M)*&
                         VddPi(t3, 1)+L*N*(M*M-1)*VddDelt(t3, 1))*(-1)
                    eS = (3*L*M*M*N*VddSig(t3, 2)+L*N*(1-4*M*M)*&
                         VddPi(t3, 2)+L*N*(M*M-1)*VddDelt(t3, 2))*(-1)
                case(6) ! dyz-dyz
                    eH = 3*M*M*N*N*VddSig(t3, 1)+(M*M+N*N-4*M*M*N*N)*&
                         VddPi(t3, 1) + (L*L+M*M*N*N)*VddDelt(t3, 1)
                    eS = 3*M*M*N*N*VddSig(t3, 2)+(M*M+N*N-4*M*M*N*N)*&
                         VddPi(t3, 2) + (L*L+M*M*N*N)*VddDelt(t3, 2)
                case(7) ! dyz-dz^2r^2
                    eH = (s3*M*N*(N*N-0.5*(L*L+M*M))*VddSig(t3, 1) + &
                         s3*M*N*(L*L+M*M-N*N)*VddPi(t3, 1)-(s3/2)*M*N*&
                         (L*L+M*M)*VddDelt(t3, 1))*(-1)
                    eS = (s3*M*N*(N*N-0.5*(L*L+M*M))*VddSig(t3, 2) + &
                         s3*M*N*(L*L+M*M-N*N)*VddPi(t3, 2)-(s3/2)*M*N*&
                         (L*L+M*M)*VddDelt(t3, 2))*(-1)
                case(8) ! dyz-dxz
                    eH = 3*L*M*N*N*VddSig(t3, 1)+L*M*(1-4*N*N)*&
                         VddPi(t3, 1)+L*M*(N*N-1)*VddDelt(t3, 1)
                    eS = 3*L*M*N*N*VddSig(t3, 2)+L*M*(1-4*N*N)*&
                         VddPi(t3, 2)+L*M*(N*N-1)*VddDelt(t3, 2)
                case(9) ! dyz-dx^2y^2
                    eH = (1.5*M*N*(L*L-M*M)*VddSig(t3, 1)-M*N*&
                         (1+2*(L*L-M*M))*VddPi(t3, 1)+M*N*(1+0.5*&
                         (L*L-M*M))*VddDelt(t3, 1))*(-1)
                    eS = (1.5*M*N*(L*L-M*M)*VddSig(t3, 2)-M*N*&
                         (1+2*(L*L-M*M))*VddPi(t3, 2)+M*N*(1+0.5*&
                         (L*L-M*M))*VddDelt(t3, 2))*(-1)
            end select
        elseif (indx(ii, 2) == 7) then
            select case(indx(jj, 2))
                case(1) ! dz^2r^2-s
                    eH = (N*N-0.5*(L*L+M*M))*Vds(t3, 1)
                    eS = (N*N-0.5*(L*L+M*M))*Vds(t3, 2)
                case(2) ! dz^2r^2-py
                    eH = M*(N*N-0.5*(L*L+M*M))*VdpSig(t3, 1) - s3*M*N*N*&
                         VdpPi(t3, 1)
                    eS = M*(N*N-0.5*(L*L+M*M))*VdpSig(t3, 2) - s3*M*N*N*&
                         VdpPi(t3, 2)
                case(3) ! dz^2r^2-pz
                    eH = (N*(N*N-0.5*(L*L+M*M))*VdpSig(t3, 1)+s3*N*&
                          (L*L+M*M)*VdpPi(t3, 1))*(-1)
                    eS = (N*(N*N-0.5*(L*L+M*M))*VdpSig(t3, 2)+s3*N*&
                          (L*L+M*M)*VdpPi(t3, 2))*(-1)
                case(4) ! dz^2r^2-px
                    eH = L*(N*N-0.5*(L*L+M*M))*VdpSig(t3, 1)-s3*L*N*N*&
                         VdpPi(t3, 1)
                    eS = L*(N*N-0.5*(L*L+M*M))*VdpSig(t3, 2)-s3*L*N*N*&
                         VdpPi(t3, 2)
                case(5) ! dz^2r^2-dxy
                    eH = s3*L*M*(N*N-0.5*(L*L+M*M))*VddSig(t3, 1)-2*s3*&
                         L*M*N*N*VddPi(t3, 1)+(s3/2)*L*M*(1+N*N)*&
                         VddDelt(t3, 1)
                    eS = s3*L*M*(N*N-0.5*(L*L+M*M))*VddSig(t3, 2)-2*s3*&
                         L*M*N*N*VddPi(t3, 2)+(s3/2)*L*M*(1+N*N)*&
                         VddDelt(t3, 2)
                case(6) ! dz^2r^2-dyz
                    eH = (s3*M*N*(N*N-0.5*(L*L+M*M))*VddSig(t3, 1) + &
                         s3*M*N*(L*L+M*M-N*N)*VddPi(t3, 1)-(s3/2)*M*N*&
                         (L*L+M*M)*VddDelt(t3, 1))*(-1)
                    eS = (s3*M*N*(N*N-0.5*(L*L+M*M))*VddSig(t3, 2) + &
                         s3*M*N*(L*L+M*M-N*N)*VddPi(t3, 2)-(s3/2)*M*N*&
                         (L*L+M*M)*VddDelt(t3, 2))*(-1)
                case(7) ! dz^2r^2-dz^2r^2
                    eH = (N*N-0.5*(L*L+M*M))**2*VddSig(t3, 1) + &
                         3*N*N*(L*L+M*M)*VddPi(t3, 1)+0.75*(L*L+M*M)**2&
                         *VddDelt(t3, 1)
                    eS = (N*N-0.5*(L*L+M*M))**2*VddSig(t3, 2) + &
                         3*N*N*(L*L+M*M)*VddPi(t3, 2)+0.75*(L*L+M*M)**2&
                         *VddDelt(t3, 2)
                case(8) ! dz^2r^2-dxz
                    eH = (s3*L*N*(N*N-0.5*(L*L+M*M))*VddSig(t3, 1)+s3*L*&
                         N*(L*L+M*M-N*N)*VddPi(t3, 1)-(s3/2)*L*N*&
                         (L*L+M*M)*VddDelt(t3, 1))*(-1)
                    eS = (s3*L*N*(N*N-0.5*(L*L+M*M))*VddSig(t3, 2)+s3*L*&
                         N*(L*L+M*M-N*N)*VddPi(t3, 2)-(s3/2)*L*N*&
                         (L*L+M*M)*VddDelt(t3, 2))*(-1)
                case(9) ! dz^2r^2-dx^2y^2
                    eH = (s3/2)*(L*L-M*M)*(N*N-0.5*(L*L+M*M))*&
                         VddSig(t3, 1)+s3*N*N*(M*M-L*L)*VddPi(t3, 1)+&
                         (s3/4)*(1+N*N)*(L*L-M*M)*VddDelt(t3, 1)
                    eS = (s3/2)*(L*L-M*M)*(N*N-0.5*(L*L+M*M))*&
                         VddSig(t3, 2)+s3*N*N*(M*M-L*L)*VddPi(t3, 2)+&
                         (s3/4)*(1+N*N)*(L*L-M*M)*VddDelt(t3, 2)
            end select
        elseif (indx(ii, 2) == 8) then
            select case(indx(jj, 2))
                case(1) ! dxz-s
                    eH = s3*L*N*Vds(t3, 1)*(-1)
                    eS = s3*L*N*Vds(t3, 2)*(-1)
                case(2) ! dxz-py
                    eH = (s3*L*M*N*VdpSig(t3, 1)-2*L*M*N*VdpPi(t3, 1))*&
                         (-1)
                    eS = (s3*L*M*N*VdpSig(t3, 2)-2*L*M*N*VdpPi(t3, 2))*&
                         (-1)
                case(3) ! dxz-pz
                    eH = s3*N*N*L*VdpSig(t3, 1)+L*(1-2*N*N)*VdpPi(t3, 1)
                    eS = s3*N*N*L*VdpSig(t3, 2)+L*(1-2*N*N)*VdpPi(t3, 2)
                case(4) ! dxz-px
                    eH = (s3*L*L*N*VdpSig(t3, 1)+N*(1-2*L*L)*&
                         VdpPi(t3, 1))*(-1)
                    eS = (s3*L*L*N*VdpSig(t3, 2)+N*(1-2*L*L)*&
                         VdpPi(t3, 2))*(-1)
                case(5) ! dxz-dxy
                    eH = (3*L*L*M*N*VddSig(t3, 1)+M*N*(1-4*L*L)*&
                         VddPi(t3, 1)+M*N*(L*L-1)*VddDelt(t3, 1))*(-1)
                    eS = (3*L*L*M*N*VddSig(t3, 2)+M*N*(1-4*L*L)*&
                         VddPi(t3, 2)+M*N*(L*L-1)*VddDelt(t3, 2))*(-1)
                case(6) ! dxz-dyz
                    eH = 3*L*M*N*N*VddSig(t3, 1)+L*M*(1-4*N*N)*&
                         VddPi(t3, 1)+L*M*(N*N-1)*VddDelt(t3, 1)
                    eS = 3*L*M*N*N*VddSig(t3, 2)+L*M*(1-4*N*N)*&
                         VddPi(t3, 2)+L*M*(N*N-1)*VddDelt(t3, 2)
                case(7) ! dxz-dz^2r^2
                    eH = (s3*L*N*(N*N-0.5*(L*L+M*M))*VddSig(t3, 1)+s3*L*&
                         N*(L*L+M*M-N*N)*VddPi(t3, 1)-(s3/2)*L*N*&
                         (L*L+M*M)*VddDelt(t3, 1))*(-1)
                    eS = (s3*L*N*(N*N-0.5*(L*L+M*M))*VddSig(t3, 2)+s3*L*&
                         N*(L*L+M*M-N*N)*VddPi(t3, 2)-(s3/2)*L*N*&
                         (L*L+M*M)*VddDelt(t3, 2))*(-1)
                case(8) ! dxz-dxz
                    eH = 3*L*L*N*N*VddSig(t3, 1)+(L*L+N*N-4*L*L*N*N)*&
                         VddPi(t3, 1) + (M*M+L*L*N*N)*VddDelt(t3, 1)
                    eS = 3*L*L*N*N*VddSig(t3, 2)+(L*L+N*N-4*L*L*N*N)*&
                         VddPi(t3, 2) + (M*M+L*L*N*N)*VddDelt(t3, 2)
                case(9) ! dxz-dx^2y^2
                    eH = (1.5*N*L*(L*L-M*M)*VddSig(t3, 1)+N*L*(1-2*&
                         (L*L-M*M))*VddPi(t3, 1)-N*L*(1-0.5*(L*L-M*M))*&
                         VddDelt(t3, 1))*(-1)
                    eS = (1.5*N*L*(L*L-M*M)*VddSig(t3, 2)+N*L*(1-2*&
                         (L*L-M*M))*VddPi(t3, 2)-N*L*(1-0.5*(L*L-M*M))*&
                         VddDelt(t3, 2))*(-1)
            end select
        elseif (indx(ii, 2) == 9) then
            select case(indx(jj, 2))
                case(1) ! dx^2y^2-s
                    eH = (s3/2)*(L*L-M*M)*Vds(t3, 1)
                    eS = (s3/2)*(L*L-M*M)*Vds(t3, 2)
                case(2) ! dx^2y^2-py
                    eH = (s3/2)*M*(L*L-M*M)*VdpSig(t3, 1) - M*(1+L*L-&
                         M*M)*VdpPi(t3, 1)
                    eS = (s3/2)*M*(L*L-M*M)*VdpSig(t3, 2) - M*(1+L*L-&
                         M*M)*VdpPi(t3, 2)
                case(3) ! dx^2y^2-pz
                    eH = ((s3/2)*N*(L*L-M*M)*VdpSig(t3, 1)-N*(L*L-M*M)*&
                         VdpPi(t3, 1))*(-1)
                    eS = ((s3/2)*N*(L*L-M*M)*VdpSig(t3, 2)-N*(L*L-M*M)*&
                         VdpPi(t3, 2))*(-1)
                case(4) ! dx^2y^2-px
                    eH = (s3/2)*L*(L*L-M*M)*VdpSig(t3, 1)+L*(1-L*L+M*M)*&
                         VdpPi(t3, 1)
                    eS = (s3/2)*L*(L*L-M*M)*VdpSig(t3, 2)+L*(1-L*L+M*M)*&
                         VdpPi(t3, 2)
                case(5) ! dx^2y^2-dxy
                    eH = 1.5*L*M*(L*L-M*M)*VddSig(t3, 1)+2*L*M*(M*M-L*L)&
                         *VddPi(t3, 1)+0.5*L*M*(L*L-M*M)*VddDelt(t3, 1)
                    eS = 1.5*L*M*(L*L-M*M)*VddSig(t3, 2)+2*L*M*(M*M-L*L)&
                         *VddPi(t3, 2)+0.5*L*M*(L*L-M*M)*VddDelt(t3, 2)
                case(6) ! dx^2y^2-dyz
                    eH = (1.5*M*N*(L*L-M*M)*VddSig(t3, 1)-M*N*(1+2*&
                         (L*L-M*M))*VddPi(t3, 1)+M*N*(1+0.5*(L*L-M*M))*&
                         VddDelt(t3, 1))*(-1)
                    eS = (1.5*M*N*(L*L-M*M)*VddSig(t3, 2)-M*N*(1+2*&
                         (L*L-M*M))*VddPi(t3, 2)+M*N*(1+0.5*(L*L-M*M))*&
                         VddDelt(t3, 2))*(-1)
                case(7) ! dx^2y^2-dz^2r^2
                    eH = (s3/2)*(L*L-M*M)*(N*N-0.5*(L*L+M*M))*&
                         VddSig(t3, 1)+s3*N*N*(M*M-L*L)*VddPi(t3, 1)+&
                         (s3/4)*(1+N*N)*(L*L-M*M)*VddDelt(t3, 1)
                    eS = (s3/2)*(L*L-M*M)*(N*N-0.5*(L*L+M*M))*&
                         VddSig(t3, 2)+s3*N*N*(M*M-L*L)*VddPi(t3, 2)+&
                         (s3/4)*(1+N*N)*(L*L-M*M)*VddDelt(t3, 2)
                case(8) ! dx^2y^2-dxz
                    eH = (1.5*N*L*(L*L-M*M)*VddSig(t3, 1)+N*L*(1-2*&
                         (L*L-M*M))*VddPi(t3, 1)-N*L*(1-0.5*(L*L-M*M))*&
                         VddDelt(t3, 1))*(-1)
                    eS = (1.5*N*L*(L*L-M*M)*VddSig(t3, 2)+N*L*(1-2*&
                         (L*L-M*M))*VddPi(t3, 2)-N*L*(1-0.5*(L*L-M*M))*&
                         VddDelt(t3, 2))*(-1)
                case(9) ! dx^2y^2-dx^2y^2
                    eH = 0.75*(L*L-M*M)**2*VddSig(t3, 1)+(L*L+M*M - &
                         (L*L-M*M)**2)*VddPi(t3, 1)+(N*N+0.25*(L*L- &
                         M*M)**2)*VddDelt(t3, 1)
                    eS = 0.75*(L*L-M*M)**2*VddSig(t3, 2)+(L*L+M*M - &
                         (L*L-M*M)**2)*VddPi(t3, 2)+(N*N+0.25*(L*L- &
                         M*M)**2)*VddDelt(t3, 2)
            end select
        else 
            write(*, '(a42)')'Error: contains extra orbital interactions'
        endif

        ! This is to make TB hamiltonian
        H(ii, jj) = H(ii, jj) + eH*phase
        S(ii, jj) = S(ii, jj) + eS*phase

    enddo ! t  - cols/indices
    enddo ! ii - rows

    ! We now have hamiltonian/overlap matrix for given k-point
    ! Finds energy eigenvalues
    call chegv(1, 'N', 'U', no_u, H, no_u, S, no_u, energy, &
                work, 2*no_u-1, rwork, ii)
    if (ii /= 0) then
        write(*, '(a30, i5)') 'Error: chegv unsuccessful exit', ii
        write(*, '(a25, 3f12.6)') 'Error occured at k-point:', vs(:)
        call exit(0)
    endif

    ! ouputs energy data
    write(5, '(f12.6)', advance='no') xout(t2)
    do ii = 1, no_u
        write(5, '(f12.6)', advance='no') energy(ii)
    enddo ! ii
    write(5, *) ! end line

    enddo ! t2 - kpoints

    ! outputs x labels/values
    write(5, *)
    do ii=1, LL+1
        if (xlabInd(ii)==1) then
            write(5, '(a1, f12.6)') xlabF(ii, 1), sum(mags(1:ii-1))
        else
            write(5, '(a3, f12.6)') xlabF(ii, :), sum(mags(1:ii-1))
        endif
    enddo ! ii

    close(5) ! closes TB energy output file

END PROGRAM MAIN
