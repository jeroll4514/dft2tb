
! B           = matrix of reciprocal lattice vectors (one per column)
! colData     = column for given matrix entry (array of col indices)
! distInd     = given distance index
! distIndData = index of distance/interaction with respect to ordering within
!               SKPARAM file (array of distance indices)
! eH          = given entry (without phase) for hamiltonian matrix
! file_exists = used for checking/writing output files
! file_status = used for checking/writing output files
! H           = hamiltonian matrix
! HMonsiteA/B = onsite matrices for Group IV metal atoms A and B
! Hsparce     = sparce hamiltonian entries (from HSX file)
! HXonsiteA/B = onsite matrices for chalcogen atoms A and B
! kden        = array with density of given k-point path
! kp          = vectors defining k-point paths
! kpE         = vectors defining endpoints for k-point paths
! kvec        = vector for given k-point to analyze within matrices
! L/M/N       = x/y/z-coordinate of normalized xij vector
! mags        = magnitudes of k-point paths
! nbands      = the number of bands
!               16 -> s and p orbitals
!               36 -> s, p, and d orbitals
! npaths      = number of k-point paths to plot
! numElem     = number of entries within Slater-Koster matrices
! phase       = phase factor between interacting orbitals
! rowData     = row for given matrix entry (array of row indices)
! S           = overlap matrix
! Ssparce     = sparce overlap matrices (from HSX file)
! t/tr        = temporary integer/real variable
! vs/ve       = start/end vector respectively for k-point paths
! xij         = vectors between interacting orbitals
! xijNorm     = normalized xij vectors.  Used for Slater-Koster integrals
! xlab(F)     = (final) array of labels for k-point paths
! xout        = x-values for output plot

PROGRAM MAIN
    IMPLICIT NONE

    CHARACTER(LEN=3)   file_status
    LOGICAL            file_exists

    INTEGER nbands, numElem, row, col, distInd, t, ii, jj, kk
    INTEGER npaths, info
    REAL tr, L, M, N, B(3, 3), kvec(3), eH, eS, vs(3), ve(3), s3, pi
    COMPLEX phase, e

    CHARACTER, ALLOCATABLE :: xlab(:), xlabF(:, :)
    INTEGER, ALLOCATABLE :: xlabInd(:)
    INTEGER, ALLOCATABLE :: rowData(:), colData(:), distIndData(:), kden(:)
    REAL, ALLOCATABLE :: HMonsiteA(:, :), HMonsiteB(:, :)
    REAL, ALLOCATABLE :: HXonsiteA(:, :), HXonsiteB(:, :)
    REAL, ALLOCATABLE :: Hss(:), Hsp(:), Hps(:), HppSig(:), HppPi(:)
    REAL, ALLOCATABLE :: Hsd(:), Hds(:), HpdSig(:), HpdPi(:), HdpSig(:), HdpPi(:)
    REAL, ALLOCATABLE :: HddSig(:), HddPi(:), HddDelt(:)
    REAL, ALLOCATABLE :: Sss(:), Ssp(:), Sps(:), SppSig(:), SppPi(:)
    REAL, ALLOCATABLE :: Ssd(:), Sds(:), SpdSig(:), SpdPi(:), SdpSig(:), SdpPi(:)
    REAL, ALLOCATABLE :: SddSig(:), SddPi(:), SddDelt(:), Ssparce(:), Hsparce(:)
    REAL, ALLOCATABLE :: xij(:, :), xijNorm(:, :), kpE(:, :), kp(:, :), mags(:), xout(:)
    REAL, ALLOCATABLE :: energy(:), rwork(:), eVals(:, :)
    COMPLEX, ALLOCATABLE :: work(:), H(:, :), S(:, :)

    !-------------------------------------------------------------!
    ! This section defines constants for use later in code        !
    !-------------------------------------------------------------!

    s3 = sqrt(3.0)
    pi = 4.0*atan(1.0)
    e  = complex(0, 1)

    !-------------------------------------------------------------!
    ! This section inputs all data from xij.dat file.  This input !
    ! file is obtained through the hsxREADERv4-0.py or newer code !
    !-------------------------------------------------------------!
        
    open(unit=1, file='xij.dat', status='old') ! opens file for input

    read(1, *) ! skips comment
    read(1, *) B(:, 1) ! reads reciprocal lattice vectors (for making kpoint-paths)
    read(1, *) B(:, 2)
    read(1, *) B(:, 3)
  
    read(1, *) nbands ! inputs number of bands (16:s, p)(36:s, p, d)
    read(1, *) ! skips comment

    ! defines size of onsite matrices -> (nbands/4 , nbands/4)
    allocate(HMonsiteA(nbands/4, nbands/4), HMonsiteB(nbands/4, nbands/4))
    allocate(HXonsiteA(nbands/4, nbands/4), HXonsiteB(nbands/4, nbands/4))
    allocate(H(nbands, nbands), S(nbands, nbands))
    allocate(energy(nbands), work(2*nbands-1), rwork(3*nbands-2)) ! for eigvals

    ! reads onsite matrices
    do ii=1, 4
        do jj=1, nbands/4
            select case(ii)
            case(1)
                read(1, *) HMonsiteA(jj, :)
            case(2)
                read(1, *) HMonsiteB(jj, :)
            case(3)
                read(1, *) HXonsiteA(jj, :)
            case(4)
                read(1, *) HXonsiteB(jj, :)
            end select
        enddo ! jj
        if (ii == 4) cycle ! does not skip line on last iteration
        read(1, *) ! skips comment
    enddo ! ii

    read(1, *) numElem ! stores number of entries within matrices
    t = numElem ! for the sake of space
    allocate(rowData(t), colData(t), distIndData(t))
    allocate(xij(t, 3), xijNorm(t, 3), Hsparce(t), Ssparce(t))
    read(1, *) ! skips comment
    do ii=1, t
        read(1, *) rowData(ii), colData(ii), distIndData(ii), xij(ii, :), xijNorm(ii, :)&
                  , Hsparce(ii), Ssparce(ii)
    enddo ! ii

    ! Note that WE DO NOT READ THE DISTANCE/INTERACTION.  That is just to 
    ! help with possible error checking in the future.

    close(1)

    !-------------------------------------------------------------!
    ! This section inputs all data from SKPARAM.dat file. This    !
    ! input file is obtained through the hsxREADERv4-0.py or      !
    ! newer code                                                  !
    !-------------------------------------------------------------!

    open(unit=2, file='SKPARAM.dat', status='old') ! opens file for input

    read(2, *)
    read(2, *)
    read(2, *)

    read(2, *) t ! reads number of unique dists/inters
    allocate(Hss(t), Hsp(t), Hps(t), HppSig(t), HppPi(t))
    allocate(Sss(t), Ssp(t), Sps(t), SppSig(t), SppPi(t))

    if (nbands == 36) then
        allocate(Hsd(t), Hds(t), HpdSig(t), HpdPi(t), HdpSig(t), HdpPi(t))
        allocate(HddSig(t), HddPi(t), HddDelt(t))
        allocate(Ssd(t), Sds(t), SpdSig(t), SpdPi(t), SdpSig(t), SdpPi(t))
        allocate(SddSig(t), SddPi(t), SddDelt(t))
    endif

    read(2, *)
    read(2, *)
    read(2, *)
    do ii=1, t ! reads SK Hamiltonian parameters for 16+ bands
        read(2, *) tr, Hss(ii), Hsp(ii), Hps(ii), HppSig(ii), HppPi(ii)
    enddo ! ii

    if (nbands == 36) then
        read(2, *)
        read(2, *)
        read(2, *)
        do ii=1, t ! reads SK Hamiltonian parameters for 36 bands
            read(2, *) tr, Hsd(ii), Hds(ii), HpdSig(ii), HpdPi(ii), HdpSig(ii), HdpPi(ii), &
                      HddSig(ii), HddPi(ii), HddDelt(ii)
        enddo ! ii
    endif

    read(2, *)
    read(2, *)
    read(2, *)

    do ii=1, t ! reads SK Overlap parameters for 16+ bands
        read(2, *) tr, Sss(ii), Ssp(ii), Sps(ii), SppSig(ii), SppPi(ii)
    enddo ! ii

    if (nbands == 36) then
        read(2, *)
        read(2, *)
        read(2, *)
        do ii=1, t ! reads SK Overlap paramters for 36 bands
            read(2, *) tr, Ssd(ii), Sds(ii), SpdSig(ii), SpdPi(ii), SdpSig(ii), SdpPi(ii), &
                      SddSig(ii), SddPi(ii), SddDelt(ii)
        enddo ! ii
    endif

    close(2)
    
    !-------------------------------------------------------------!
    ! This section inputs all data from kpoints.dat file. This    !
    ! file must be made by hand                                   !
    !-------------------------------------------------------------!

    open(unit=3, file='kpoints.dat', status='old') ! opens file for input

    read(3, *) npaths ! inputs number of kpoint-paths within file
    allocate(kden(npaths), kpE(2*npaths, 3), mags(npaths))
    allocate(xlab(2*npaths), xlabF(npaths+1, 3), xlabInd(npaths+1))

    do ii=1, npaths
        read(3, *) xlab(2*ii-1), kpE(2*ii-1, :), xlab(2*ii), kpE(2*ii, :), kden(ii)
    enddo ! ii

    ! Checks for discontinuities within k-point path
    xlabF(1, :) = (/xlab(1), '', ''/)
    xlabInd(1) = 1
    do ii=1, npaths-1
        if (xlab(2*ii)==xlab(2*ii+1)) then
            xlabF(ii+1, :)=(/xlab(2*ii), '', ''/)
            xlabInd(ii+1) = 1
        else
            xlabF(ii+1, :) = (/ xlab(2*ii), '|', xlab(2*ii+1) /)
            xlabInd(ii+1) = 3
        endif
    enddo
    xlabF(npaths+1, :) = (/xlab(2*npaths), '', ''/)
    xlabInd(npaths+1) = 1

    allocate(kp(3, sum(kden)+npaths)) ! (1) is x coord. (2) is y coord. (3) is z coord
    allocate(xout(sum(kden)+npaths))

    close(3)

    !-------------------------------------------------------------!
    ! This section creates the k-point paths                      !
    !-------------------------------------------------------------!

    do ii=1, npaths ! iterates through each path
        vs = matmul(B(:, :), kpE(2*ii-1, :)) ! start coord in cartesian coords
        ve = matmul(B(:, :), kpE(2*ii, :)) ! end coord in cartesian coords
        mags(ii) = norm2(ve-vs) ! magnitude of k-point path
        do jj=1, kden(ii)+1 ! iterates over k-point density
            kp(1, jj+sum(kden(1:ii-1))+ii-1) = vs(1)+((jj-1.0)/kden(ii))*(ve(1)-vs(1))
            kp(2, jj+sum(kden(1:ii-1))+ii-1) = vs(2)+((jj-1.0)/kden(ii))*(ve(2)-vs(2))
            kp(3, jj+sum(kden(1:ii-1))+ii-1) = vs(3)+((jj-1.0)/kden(ii))*(ve(3)-vs(3))
            xout(jj+sum(kden(1:ii-1))+ii-1) = sum(mags(1:ii-1))+((jj-1.0)/kden(ii))*mags(ii)
        enddo ! jj
    enddo ! ii
    
    allocate(eVals(sum(kden)+npaths, nbands))

    !-------------------------------------------------------------!
    ! This section creates the hamiltonian/overlap matrices.      !
    !-------------------------------------------------------------!

    do ii=1, sum(kden)+npaths

    ! Resets matrices for each k-point
    H = 0.0
    S = 0.0

    ! Enters onsite Hamiltonian terms
    ! Comment this block if making HSX bands
    do jj=1, nBands/4
        do kk=jj, nBands/4
            H(jj, kk) = HMonsiteA(jj, kk)
            H(jj+nBands/4, kk+nBands/4) = HMonsiteB(jj, kk)
            H(jj+nBands/2, kk+nBands/2) = HXonsiteA(jj, kk)
            H(jj+3*nBands/4, kk+3*nBands/4) = HXonsiteB(jj, kk)
        enddo ! kk
    enddo ! jj

    kvec(:) = (/ kp(1, ii), kp(2, ii), kp(3, ii) /) ! defines current k-point

    do jj=1, numElem
        
        row = rowData(jj)
        col = colData(jj)
        distInd = distIndData(jj)
        phase = exp(e*sum(xij(jj, :)*kvec(:)))

        if (distInd > 0) then ! does not find SK for onsite-terms

        L = xijNorm(jj, 1)
        M = xijNorm(jj, 2)
        N = xijNorm(jj, 3)

        if (row==1 .or. row==10 .or. row==19 .or. row==28) then
            select case(col)
                case(1, 10, 19, 28) ! s-s
                    eH = Hss(distInd)
                case(2, 11, 20, 29) ! s-py
                    eH = M*Hsp(distInd)
                case(3, 12, 21, 30) ! s-pz
                    eH = N*Hsp(distInd)*(-1)
                case(4, 13, 22, 31) ! s-px
                    eH = L*Hsp(distInd)
                case(5, 14, 23, 32) ! s-dxy
                    eH = s3*L*M*Hsd(distInd)
                case(6, 15, 24, 33) ! s-dyz
                    eH = s3*M*N*Hsd(distInd)*(-1)
                case(7, 16, 25, 34) ! s-dz^2r^2
                    eH = (N*N-0.5*(L*L+M*M))*Hsd(distInd)
                case(8, 17, 26, 35) ! s-dxz
                    eH = s3*L*N*Hsd(distInd)*(-1)
                case(9, 18, 27, 36) ! s-dx^2y^2
                    eH = (s3/2)*(L*L-M*M)*Hsd(distInd)
            end select
        elseif (row==2 .or. row==11 .or. row==20 .or. row==29) then
            select case(col)
                case(10, 19, 28) ! py-s
                    eH = M*Hps(distInd)
                case(2, 11, 20, 29) ! py-py
                    eH = M*M*HppSig(distInd) + (1-M*M)*HppPi(distInd)
                case(3, 12, 21, 30) ! py-pz
                    eH = M*N*(HppSig(distInd)-HppPi(distInd))*(-1)
                case(4, 13, 22, 31) ! py-px
                    eH = M*L*(HppSig(distInd)-HppPi(distInd))
                case(5, 14, 23, 32) ! py-dxy
                    eH = s3*M*M*L*HpdSig(distInd)+L*(1-2*M*M)*HpdPi(distInd)
                case(6, 15, 24, 33) ! py-dyz
                    eH = (s3*M*M*N*HpdSig(distInd)+N*(1-2*M*M)*HpdPi(distInd))&
                         *(-1)
                case(7, 16, 25, 34) ! py-dz^2r^2
                    eH = M*(N*N-0.5*(L*L+M*M))*HpdSig(distInd) - s3*M*N*N*&
                         HpdPi(distInd)
                case(8, 17, 26, 35) ! py-dxz
                    eH = (s3*L*M*N*HpdSig(distInd)-2*L*M*N*HpdPi(distInd))*(-1)
                case(9, 18, 27, 36) ! py-dx^2y^2
                    eH = (s3/2)*M*(L*L-M*M)*HpdSig(distInd) - M*(1+L*L-&
                         M*M)*HpdPi(distInd)
            end select
        elseif (row==3 .or. row==12 .or. row==21 .or. row==30) then
            select case(col)
                case(10, 19, 28) ! pz-s
                    eH = N*Hps(distInd)*(-1)
                case(11, 20, 29) ! pz-py
                    eH = N*M*(HppSig(distInd)-HppPi(distInd))*(-1)
                case(3, 12, 21, 30) ! pz-pz
                    eH = N*N*HppSig(distInd) + (1-N*N)*HppPi(distInd)
                case(4, 13, 22, 31) ! pz-px
                    eH = L*N*(HppSig(distInd)-HppPi(distInd))*(-1)
                case(5, 14, 23, 32) ! pz-dxy
                    eH = (s3*L*M*N*HpdSig(distInd)-2*L*M*N*HpdPi(distInd))*(-1)
                case(6, 15, 24, 33) ! pz-dyz
                    eH = s3*N*N*M*HpdSig(distInd)+M*(1-2*N*N)*HpdPi(distInd)
                case(7, 16, 25, 34) ! pz-dz^2r^2
                    eH = (N*(N*N-0.5*(L*L+M*M))*HpdSig(distInd)+s3*N*&
                          (L*L+M*M)*HpdPi(distInd))*(-1)
                case(8, 17, 26, 35) ! pz-dxz
                    eH = s3*N*N*L*HpdSig(distInd)+L*(1-2*N*N)*HpdPi(distInd)
                case(9, 18, 27, 36) ! pz-dx^2y^2
                    eH = ((s3/2)*N*(L*L-M*M)*HpdSig(distInd) - N*(L*L-M*M)*&
                         HpdPi(distInd))*(-1)
            end select
        elseif (row==4 .or. row==13 .or. row==22 .or. row==31) then
            select case(col)
                case(10, 19, 28) ! px-s
                    eH = L*Hps(distInd)
                case(11, 20, 29) ! px-py
                    eH = L*M*(HppSig(distInd)-HppPi(distInd))
                case(12, 21, 30) ! px-pz
                    eH = L*N*(HppSig(distInd)-HppPi(distInd))*(-1)
                case(4, 13, 22, 31) ! px-px
                    eH = L*L*HppSig(distInd) + (1-L*L)*HppPi(distInd)
                case(5, 14, 23, 32) ! px-dxy
                    eH = s3*L*L*M*HpdSig(distInd)+M*(1-2*L*L)*HpdPi(distInd)
                case(6, 15, 24, 33) ! px-dyz
                    eH = (s3*L*M*N*HpdSig(distInd)-2*L*M*N*HpdPi(distInd))*(-1)
                case(7, 16, 25, 34) ! px-dz^2r^2
                    eH = L*(N*N-0.5*(L*L+M*M))*HpdSig(distInd)-s3*L*N*N*&
                         HpdPi(distInd)
                case(8, 17, 26, 35) ! px-dxz
                    eH = (s3*L*L*N*HpdSig(distInd)+N*(1-2*L*L)*HpdPi(distInd))&
                         *(-1)
                case(9, 18, 27, 36) ! px-dx^2y^2
                    eH = (s3/2)*L*(L*L-M*M)*HpdSig(distInd)+L*(1-L*L+M*M)*&
                         HpdPi(distInd)
            end select
        elseif (row==5 .or. row==14 .or. row==23 .or. row==32) then
            select case(col)
                case(10, 19, 28) ! dxy-s
                    eH = s3*L*M*Hds(distInd)
                case(11, 20, 29) ! dxy-py
                    eH = s3*M*M*L*HdpSig(distInd)+L*(1-2*M*M)*HdpPi(distInd)
                case(12, 21, 30) ! dxy-pz
                    eH = (s3*L*M*N*HdpSig(distInd)-2*L*M*N*HdpPi(distInd))*(-1)
                case(13, 22, 31) ! dxy-px
                    eH = s3*L*L*M*HdpSig(distInd)+M*(1-2*L*L)*HdpPi(distInd)
                case(5, 14, 23, 32) ! dxy-dxy
                    eH = 3*L*L*M*M*HddSig(distInd)+(L*L+M*M-4*L*L*M*M)*&
                         HddPi(distInd) + (N*N+L*L*M*M)*HddDelt(distInd)
                case(6, 15, 24, 33) ! dxy-dyz
                    eH = (3*L*M*M*N*HddSig(distInd)+L*N*(1-4*M*M)*HddPi(distInd)&
                         +L*N*(M*M-1)*HddDelt(distInd))*(-1)
                case(7, 16, 25, 34) ! dxy-dz^2r^2
                    eH = s3*L*M*(N*N-0.5*(L*L+M*M))*HddSig(distInd)-2*s3*L*&
                         M*N*N*HddPi(distInd)+(s3/2)*L*M*(1+N*N)*&
                         HddDelt(distInd)
                case(8, 17, 26, 35) ! dxy-dxz
                    eH = (3*L*L*M*N*HddSig(distInd)+M*N*(1-4*L*L)*HddPi(distInd)&
                         +M*N*(L*L-1)*HddDelt(distInd))*(-1)
                case(9, 18, 27, 36) ! dxy-dx^2y^2
                    eH = 1.5*L*M*(L*L-M*M)*HddSig(distInd)+2*L*M*(M*M-L*L)*&
                         HddPi(distInd)+0.5*L*M*(L*L-M*M)*HddDelt(distInd)
            end select
        elseif (row==6 .or. row==15 .or. row==24 .or. row==33) then
            select case(col)
                case(10, 19, 28) ! dyz-s
                    eH = s3*M*N*Hds(distInd)*(-1)
                case(11, 20, 29) ! dyz-py
                    eH = (s3*M*M*N*HdpSig(distInd)+N*(1-2*M*M)*HdpPi(distInd))&
                         *(-1)
                case(12, 21, 30) ! dyz-pz
                    eH = s3*N*N*M*HdpSig(distInd)+M*(1-2*N*N)*HdpPi(distInd)
                case(13, 22, 31) ! dyz-px
                    eH = (s3*L*M*N*HdpSig(distInd)-2*L*M*N*HdpPi(distInd))*(-1)
                case(14, 23, 32) ! dyz-dxy
                    eH = (3*L*M*M*N*HddSig(distInd)+L*N*(1-4*M*M)*HddPi(distInd)&
                         +L*N*(M*M-1)*HddDelt(distInd))*(-1)
                case(6, 15, 24, 33) ! dyz-dyz
                    eH = 3*M*M*N*N*HddSig(distInd)+(M*M+N*N-4*M*M*N*N)*&
                         HddPi(distInd) + (L*L+M*M*N*N)*HddDelt(distInd)
                case(7, 16, 25, 34) ! dyz-dz^2r^2
                    eH = (s3*M*N*(N*N-0.5*(L*L+M*M))*HddSig(distInd) + &
                         s3*M*N*(L*L+M*M-N*N)*HddPi(distInd)-(s3/2)*M*N*&
                         (L*L+M*M)*HddDelt(distInd))*(-1)
                case(8, 17, 26, 35) ! dyz-dxz
                    eH = 3*L*M*N*N*HddSig(distInd)+L*M*(1-4*N*N)*HddPi(distInd)&
                         +L*M*(N*N-1)*HddDelt(distInd)
                case(9, 18, 27, 36) ! dyz-dx^2y^2
                    eH = (1.5*M*N*(L*L-M*M)*HddSig(distInd)-M*N*(1+2*(L*L-M*M))*&
                         HddPi(distInd)+M*N*(1+0.5*(L*L-M*M))*HddDelt(distInd))*&
                         (-1)
            end select
        elseif (row==7 .or. row==16 .or. row==25 .or. row==34) then
            select case(col)
                case(10, 19, 28) ! dz^2r^2-s
                    eH = (N*N-0.5*(L*L+M*M))*Hds(distInd)
                case(11, 20, 29) ! dz^2r^2-py
                    eH = M*(N*N-0.5*(L*L+M*M))*HdpSig(distInd) - s3*M*N*N*&
                         HdpPi(distInd)
                case(12, 21, 30) ! dz^2r^2-pz
                    eH = (N*(N*N-0.5*(L*L+M*M))*HdpSig(distInd)+s3*N*&
                          (L*L+M*M)*HdpPi(distInd))*(-1)
                case(13, 22, 31) ! dz^2r^2-px
                    eH = L*(N*N-0.5*(L*L+M*M))*HdpSig(distInd)-s3*L*N*N*&
                         HdpPi(distInd)
                case(14, 23, 32) ! dz^2r^2-dxy
                    eH = s3*L*M*(N*N-0.5*(L*L+M*M))*HddSig(distInd)-2*s3*L*&
                         M*N*N*HddPi(distInd)+(s3/2)*L*M*(1+N*N)*&
                         HddDelt(distInd)
                case(15, 24, 33) ! dz^2r^2-dyz
                    eH = (s3*M*N*(N*N-0.5*(L*L+M*M))*HddSig(distInd) + &
                         s3*M*N*(L*L+M*M-N*N)*HddPi(distInd)-(s3/2)*M*N*&
                         (L*L+M*M)*HddDelt(distInd))*(-1)
                case(7, 16, 25, 34) ! dz^2r^2-dz^2r^2
                    eH = (N*N-0.5*(L*L+M*M))**2*HddSig(distInd) + &
                         3*N*N*(L*L+M*M)*HddPi(distInd)+0.75*(L*L+M*M)**2&
                         *HddDelt(distInd)
                case(8, 17, 26, 35) ! dz^2r^2-dxz
                    eH = (s3*L*N*(N*N-0.5*(L*L+M*M))*HddSig(distInd)+s3*L*N&
                         *(L*L+M*M-N*N)*HddPi(distInd)-(s3/2)*L*N*&
                         (L*L+M*M)*HddDelt(distInd))*(-1)
                case(9, 18, 27, 36) ! dz^2r^2-dx^2y^2
                    eH = (s3/2)*(L*L-M*M)*(N*N-0.5*(L*L+M*M))*HddSig(distInd)&
                         +s3*N*N*(M*M-L*L)*HddPi(distInd)+(s3/4)*(1+N*N)*&
                         (L*L-M*M)*HddDelt(distInd)
            end select
        elseif (row==8 .or. row==17 .or. row==26 .or. row==35) then
            select case(col)
                case(10, 19, 28) ! dxz-s
                    eH = s3*L*N*Hds(distInd)*(-1)
                case(11, 20, 29) ! dxz-py
                    eH = (s3*L*M*N*HdpSig(distInd)-2*L*M*N*HdpPi(distInd))*(-1)
                case(12, 21, 30) ! dxz-pz
                    eH = s3*N*N*L*HdpSig(distInd)+L*(1-2*N*N)*HdpPi(distInd)
                case(13, 22, 31) ! dxz-px
                    eH = (s3*L*L*N*HdpSig(distInd)+N*(1-2*L*L)*HdpPi(distInd))&
                         *(-1)
                case(14, 23, 32) ! dxz-dxy
                    eH = (3*L*L*M*N*HddSig(distInd)+M*N*(1-4*L*L)*HddPi(distInd)&
                         +M*N*(L*L-1)*HddDelt(distInd))*(-1)
                case(15, 24, 33) ! dxz-dyz
                    eH = 3*L*M*N*N*HddSig(distInd)+L*M*(1-4*N*N)*HddPi(distInd)&
                         +L*M*(N*N-1)*HddDelt(distInd)
                case(16, 25, 34) ! dxz-dz^2r^2
                    eH = (s3*L*N*(N*N-0.5*(L*L+M*M))*HddSig(distInd)+s3*L*N&
                         *(L*L+M*M-N*N)*HddPi(distInd)-(s3/2)*L*N*&
                         (L*L+M*M)*HddDelt(distInd))*(-1)
                case(8, 17, 26, 35) ! dxz-dxz
                    eH = 3*L*L*N*N*HddSig(distInd)+(L*L+N*N-4*L*L*N*N)*&
                         HddPi(distInd) + (M*M+L*L*N*N)*HddDelt(distInd)
                case(9, 18, 27, 36) ! dxz-dx^2y^2
                    eH = (1.5*N*L*(L*L-M*M)*HddSig(distInd)+N*L*(1-2*&
                         (L*L-M*M))*HddPi(distInd)-N*L*(1-0.5*(L*L-M*M))*&
                         HddDelt(distInd))*(-1)
            end select
        elseif (row==9 .or. row==18 .or. row==27 .or. row==36) then
            select case(col)
                case(10, 19, 28) ! dx^2y^2-s
                    eH = (s3/2)*(L*L-M*M)*Hds(distInd)
                case(11, 20, 29) ! dx^2y^2-py
                    eH = (s3/2)*M*(L*L-M*M)*HdpSig(distInd) - M*(1+L*L-&
                         M*M)*HdpPi(distInd)
                case(12, 21, 30) ! dx^2y^2-pz
                    eH = ((s3/2)*N*(L*L-M*M)*HdpSig(distInd) - N*(L*L-M*M)*&
                         HdpPi(distInd))*(-1)
                case(13, 22, 31) ! dx^2y^2-px
                    eH = (s3/2)*L*(L*L-M*M)*HdpSig(distInd)+L*(1-L*L+M*M)*&
                         HdpPi(distInd)
                case(14, 23, 32) ! dx^2y^2-dxy
                    eH = 1.5*L*M*(L*L-M*M)*HddSig(distInd)+2*L*M*(M*M-L*L)*&
                         HddPi(distInd)+0.5*L*M*(L*L-M*M)*HddDelt(distInd)
                case(15, 24, 33) ! dx^2y^2-dyz
                    eH = (1.5*M*N*(L*L-M*M)*HddSig(distInd)-M*N*(1+2*(L*L-M*M))*&
                         HddPi(distInd)+M*N*(1+0.5*(L*L-M*M))*HddDelt(distInd))*&
                         (-1)
                case(16, 25, 34) ! dx^2y^2-dz^2r^2
                    eH = (s3/2)*(L*L-M*M)*(N*N-0.5*(L*L+M*M))*HddSig(distInd)&
                         +s3*N*N*(M*M-L*L)*HddPi(distInd)+(s3/4)*(1+N*N)*&
                         (L*L-M*M)*HddDelt(distInd)
                case(17, 26, 35) ! dx^2y^2-dxz
                    eH = (1.5*N*L*(L*L-M*M)*HddSig(distInd)+N*L*(1-2*&
                         (L*L-M*M))*HddPi(distInd)-N*L*(1-0.5*(L*L-M*M))*&
                         HddDelt(distInd))*(-1)
                case(9, 18, 27, 36) ! dx^2y^2-dx^2y^2
                    eH = 0.75*(L*L-M*M)**2*HddSig(distInd)+(L*L+M*M - &
                         (L*L-M*M)**2)*HddPi(distInd)+(N*N+0.25*(L*L- &
                         M*M)**2)*HddDelt(distInd)
            end select
        else 
            write(*, *) 'Error: contains extra orbital interactions'
        endif

        ! This is to make TB hamiltonian
        H(row, col) = H(row, col) + eH*phase

        endif

        ! This is to make HSX hamiltonian/overlap matrices
!         H(row, col) = H(row, col) + Hsparce(jj)*phase
        S(row, col) = S(row, col) + Ssparce(jj)*phase
 
    enddo ! jj
    
    call chegv(1, 'N', 'U', nbands, H, 36, S, 36, energy, work, 2*nbands-1, rwork, info)
    if (info /= 0) then
        write(*, *) 'Error: chegv unsuccessful exit'
        write(*, *) info
    endif

    ! Puts energy eigenvalues within correct arrays for plotting
    do jj=1, nbands
        eVals(ii, jj) = energy(jj)
    enddo ! jj

    enddo ! ii

    !-------------------------------------------------------------!
    ! This section outputs energy eigenvalues for plotting        !
    !-------------------------------------------------------------!

    inquire(file='band_plot.dat', exist=file_exists) ! checks if output file exists

    if (file_exists .eqv. .True.) then
        file_status = 'old' ! if it exists, then override it
    else
        file_status = 'new' ! if it does not exist, create it
    end if

    open (unit=10, file='band_plot.dat', status=file_status)
    do ii=1, nbands
        do jj=1, sum(kden)+npaths
            write(10, *) xout(jj), eVals(jj, ii)
        enddo ! jj
        write(10, *)
    enddo ! ii
    close(10)

    !-------------------------------------------------------------!
    ! This section outputs gnuplot file for plotting              !
    !-------------------------------------------------------------!

    inquire(file='band_plot.gnu', exist=file_exists) ! checks if output file exists

    if (file_exists .eqv. .True.) then
        file_status = 'old' ! if it exists, then override it
    else
        file_status = 'new' ! if it does not exist, create it
    end if

    open (unit=11, file='band_plot.gnu', status=file_status)
    write(11, '(a9)') 'set nokey'
    write(11, '(a24)') 'set ylabel "Energy (eV)"'
    write(11, '(a14, f10.5, a1)') 'set xrange [0:', xout(sum(kden)+npaths), ']'
    write(11, '(a12, f10.5, a1, f10.5, a1)') 'set yrange [', minval(eVals(1, :))*1.1, ':', &
                 maxval(eVals(nbands, :))*1.1, ']'
    do ii=1, npaths-1
        write(11, '(a15, f10.5, a1, f10.5, a4, f10.5, a1, f10.5, a7)') &
                    'set arrow from ', sum(mags(1:ii)), ', ', &
                    minval(eVals(1, :))*1.1, ' to ', sum(mags(1:ii)), &
                    ', ', maxval(eVals(nbands, :))*1.1, ' nohead'
    enddo ! ii
    write(11, '(a11)', advance='no') 'set xtics ('
    do ii=1, npaths+1
        if (xlabInd(ii)==1) then
            write(11, '(a1, a1, a1, f10.5)', advance='no') '"', xlabF(ii, 1), '"', &
                 sum(mags(1:ii-1))
        else
            write(11, '(a1, a3, a1, f10.5)', advance='no') '"', xlabF(ii, :), '"', &
                 sum(mags(1:ii-1))
        endif
        if (ii<npaths+1) then
            write(11, '(a1)', advance='no') ', '
        endif
    enddo ! ii
    write(11, '(a1)') ')'
    write(11, '(a47)') 'set style line 1 lt 1 linecolor rgb "royalblue"'
    write(11, '(a29)') 'plot "band_plot.dat" ls 1 w l'
    close(11)

END PROGRAM MAIN