!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!
!                                                            !
!                                                            !
!------------------------------------------------------------!

module w90_mcae
  !! Routine for calculating magnetocrystalline anisotropy energy(MCAE).
  !! Actually what this does is summing the band energies.
  !! The band energy difference between two magnetic directions should
  !! be a good approximation to MCAE.
  !!
  !! written by Junfeng Qiao
  !! Beihang University (China)
  !! Feb, 2018

  use w90_constants
  use w90_parameters, only    : num_wann, recip_lattice, real_lattice, &
          timing_level, mcae_kmesh, mcae_adpt_kmesh, mcae_adpt_kmesh_thresh, &
          fermi_energy, mcae_adpt_smr, mcae_adpt_smr_fac, mcae_adpt_smr_max, &
          mcae_smr_fixed_en_width, mcae_smr_index, mcae_num_elec, mcae_no_smr
  use w90_io, only            : io_error,stdout,io_stopwatch,io_file_unit,seedname
  use w90_comms
  use w90_io, only            : io_date

  implicit none

  private 
  public                                          :: mcae_main

  integer                                         :: num_kpts
  real(kind=dp), dimension(:), allocatable        :: globalsum
  real(kind=dp), dimension(:,:), allocatable      :: kpoints, localkpoints
  real(kind=dp)                                   :: ef
  character(len=50)                               :: file_name
  integer                                         :: file_unit

contains
    
  subroutine internal_pre_info()
    !! This routine sums the band energies

    implicit none

    integer            :: ikpt

    !%%%%%%%%%%% variable declaration ends %%%%%%%%%%%%

    ! write header
    if (on_root .and. (timing_level>0)) call io_stopwatch('mcae_main',1)

    if (on_root) then
       write(stdout,*) 
       write(stdout,'(1x,a)') '*---------------------------------------------------------------------------*'
       write(stdout,'(1x,a)') '|                                MCAE routines                              |'
       write(stdout,'(1x,a)') '*---------------------------------------------------------------------------*'

    end if

    ! write information
    if(on_root) then

        write(stdout,'(1x,a)') ' '

        write(stdout, '(1x,a28,3(i5,1x))')&
                'Regular interpolation grid: ', mcae_kmesh(1:3)
        write(stdout, '(1x,a28,3(i5,1x))') 'Adaptive refinement grid: ',&
                mcae_adpt_kmesh
        write(stdout, '(1x,a28,a60,f6.2,a)') 'Refinement threshold: ',&
                        'Difference of MCAE between neighbour k-points >',&
                        mcae_adpt_kmesh_thresh,' meV'
        !write(stdout,'(1x,a30,i5,a,f5.2,a)')&
        !        ' Points triggering refinement: ',&
        !        adpt_counter,'(',&
        !        100*real(adpt_counter,dp)/product(berry_kmesh),'%)'
        write(stdout, '(1x,a30,f10.6)') 'Fermi energy: ', fermi_energy

        write(stdout, '(1x,a30,l2)') 'Adaptive smearing: ', mcae_adpt_smr
        write(stdout, '(1x,a30,l2)') 'No smearing: ', mcae_no_smr
        write(stdout, '(1x,a30,f10.6)') 'Adaptive_smr_fac: ', mcae_adpt_smr_fac
        write(stdout, '(1x,a30,f10.6)') 'Adaptive_smr_max: ', mcae_adpt_smr_max
        write(stdout, '(1x,a30,i4)') 'smr_index: ', mcae_smr_index
        write(stdout, '(1x,a30,f10.6)') 'smr_fixed_en_width: ', mcae_smr_fixed_en_width

    end if !on_root

    return

  end subroutine internal_pre_info

  subroutine internal_post_info()
    !! This routine sums the band energies
    implicit none

    integer            :: ikpt

    !%%%%%%%%%%% variable declaration ends %%%%%%%%%%%%
    if (on_root) then
        ! write k-resolved mcae to file
        write(stdout, '(1x,a30,f18.14)') 'fermi_energy: ', ef
        write(stdout,'(/,1x,a)')&
                '---------------------------------'
        write(stdout,'(1x,a)')&
                'Output data files related to k-resolved MCAE:'
        write(stdout,'(1x,a)')&
                '---------------------------------'
        file_name=trim(seedname)//'-mcae_k.dat'
        write(stdout,'(/,3x,a)') '* '//file_name

        if (on_root .and. (timing_level>0)) call io_stopwatch('mcae_main: write file',1)
        file_unit=io_file_unit()
        open(file_unit,FILE=file_name,STATUS='UNKNOWN',FORM='FORMATTED')
        write(file_unit,'(/,1x,a18,F18.10)') 'Total MCAE (meV):',&
                sum(globalsum)
        write(file_unit,'(/,1x,a)')&
                'index      kx          ky          kz          mcae_k'
        do ikpt=1,num_kpts

          ! First calculate the absolute coordinates for printing
          !frac = 0._dp
          !do j=1,3
          !  frac(j)=recip_lattice(1,j)*kpt(1) + recip_lattice(2,j)*kpt(2) + recip_lattice(3,j)*kpt(3)
          !end do

          write(file_unit,'(1x,i6,3(f10.7,1x),3(f18.10,1x),/)') &
                    ikpt, kpoints(:,ikpt), globalsum(ikpt)
        enddo
        close(file_unit)
        if (on_root .and. (timing_level>0)) call io_stopwatch('mcae_main: write file',2)

        write(stdout,'(1x,a)') '*---------------------------------------------------------------------------*'
        write(stdout,'(1x,a)') '|                                All done.                                  |'
        write(stdout,'(1x,a)') '*---------------------------------------------------------------------------*'

    end if

    return

  end subroutine internal_post_info


  subroutine mcae_main()

    use w90_get_oper, only      : get_HH_R, HH_R

    !!!!!!!!!!!
    !integer, dimension(:), allocatable              :: kpointidx, localkpointidx

    integer            :: ierr, enidx
    character(len=50)  :: cdum
    !real(kind=dp), dimension(3)                     :: kpt, frac
    real(kind=dp), dimension(:,:), allocatable      :: localeig
    ! first indice is band, second is kpoint
    !real(kind=dp), dimension(:,:), allocatable      :: globaleig
    ! first indice is band, second is kpoint
    real(kind=dp), dimension(:,:,:), allocatable    :: localdel_eig
    ! first indice is band, second is xyz direction, third is kpoint
    real(kind=dp), dimension(:,:), allocatable      :: locallevelspacing

    integer, dimension(0:num_nodes-1)               :: counts
    integer, dimension(0:num_nodes-1)               :: displs

    real(kind=dp), dimension(:), allocatable      :: localsum

    !! adaptive kmesh
    real(kind=dp), allocatable    :: adkpt(:,:)
    integer           :: adpt_counter

    real(kind=dp)     :: kweight,kweight_adpt,kpt(3),kpt_ad(3),&
            db1,db2,db3,fac,rdum,vdum(3)

    !! tmp variables
    integer           :: i, j, k, ikpt
    integer           :: loop_x,loop_y,loop_z,loop_xyz,loop_adpt
    character(len=200)       :: outfile_name_node

    !%%%%%%%%%%% variable declaration ends %%%%%%%%%%%%

    call internal_pre_info()

    if (on_root) then
      num_kpts = PRODUCT(mcae_kmesh)
    end if
    call comms_bcast(num_kpts,1)

    ! I call once the routine to calculate the Hamiltonian in real-space <0n|H|Rm>
    if (on_root .and. (timing_level>0)) call io_stopwatch('mcae_main: get_HH_R',1)
    call get_HH_R
    if (on_root .and. (timing_level>0)) call io_stopwatch('mcae_main: get_HH_R',2)
    if (on_root) then
        write(stdout, '(a)') 'get_HH_R done'
    end if

    if (on_root) then
       !allocate(kpointidx(num_kpts),stat=ierr)
       !if (ierr/=0) call io_error('Error allocating kpointidx in mcae_interp.')
       allocate(kpoints(3,num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error allocating kpoints in mcae_interp.')
       allocate(globalsum(num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error allocating globalsum in mcae_interp.')
    else
       ! On the other nodes, I still allocate them with size 1 to avoid 
       ! that some compilers still try to access the memory
       !allocate(kpointidx(1),stat=ierr)
       !if (ierr/=0) call io_error('Error allocating kpointidx in mcae_interp.')
       allocate(kpoints(1,1),stat=ierr)
       if (ierr/=0) call io_error('Error allocating kpoints in mcae_interp.')
       allocate(globalsum(1),stat=ierr)
       if (ierr/=0) call io_error('Error allocating globalsum in mcae_interp.')
    end if
       
    ! I precalculate how to split on different nodes
    if (on_root .and. (timing_level>0)) call io_stopwatch('mcae_main: comms_array_split',1)
    call comms_array_split(num_kpts,counts,displs)
    !write(*, '(a,i4)') 'num_nodes = ',num_nodes
    !write(*, '(a,i4,a,i4)') 'counts(',my_node_id,') = ',counts(my_node_id)
    !write(*, '(a,i4,a,i4)') 'displs(',my_node_id,') = ',displs(my_node_id)

    if (on_root .and. (timing_level>0)) call io_stopwatch('mcae_main: comms_array_split',2)
    if (on_root) then
        write(stdout, '(a)') 'comms_array_split done'
    end if

    allocate(localkpoints(3,max(1,counts(my_node_id))),stat=ierr)
    if (ierr/=0) call io_error('Error allocating localkpoints in mcae_interp.')
    
    allocate(localeig(num_wann,max(1,counts(my_node_id))),stat=ierr)
    if (ierr/=0) call io_error('Error allocating localeig in mcae_interp.')

    if (mcae_adpt_smr) then
      allocate(localdel_eig(num_wann,3,max(1,counts(my_node_id))),stat=ierr)
      if (ierr/=0) call io_error('Error allocating localdel_eig in mcae_interp.')

      allocate(locallevelspacing(num_wann,max(1,counts(my_node_id))),stat=ierr)
      if (ierr/=0) call io_error('Error allocating locallevelspacing in mcae_interp.')
    end if

    allocate(localsum(max(1,counts(my_node_id))),stat=ierr)
    if (ierr/=0) call io_error('Error allocating localsum in mcae_interp.')


    if (on_root .and. (timing_level>0)) call io_stopwatch('mcae_main: kpoints',1)
    ! On root, get all the needed kpoints
    if (on_root) then

       ! Mesh spacing in reduced coordinates
       !
       db1=1.0_dp/real(mcae_kmesh(1),dp)
       db2=1.0_dp/real(mcae_kmesh(2),dp)
       db3=1.0_dp/real(mcae_kmesh(3),dp)
       ! Set up adaptive refinement mesh
       !
       !allocate(adkpt(3,mcae_adpt_kmesh**3),stat=ierr)
       !if (ierr/=0) call io_error('Error in allocating adkpt in mcae')
       !ikpt=0
       !do i=0,mcae_adpt_kmesh-1
       !    do j=0,mcae_adpt_kmesh-1
       !        do k=0,mcae_adpt_kmesh-1
       !            ikpt=ikpt+1
       !            adkpt(1,ikpt)=db1*((i+0.5_dp)/mcae_adpt_kmesh-0.5_dp)
       !            adkpt(2,ikpt)=db2*((j+0.5_dp)/mcae_adpt_kmesh-0.5_dp)
       !            adkpt(3,ikpt)=db3*((k+0.5_dp)/mcae_adpt_kmesh-0.5_dp)
       !        end do
       !    end do
       !end do

       ! Loop over a regular grid in the full BZ
       kweight=db1*db2*db3
       !kweight_adpt=kweight/mcae_adpt_kmesh**3
       !call comms_bcast(kweight_adpt,1)

       ! (in crystallographic coordinates relative to the reciprocal lattice vectors)
       do i=1,num_kpts
          !kpointidx(i) = i

          !! modified from berry.F90
          loop_xyz = i
          loop_x = loop_xyz/(mcae_kmesh(2)*mcae_kmesh(3))
          loop_y = (loop_xyz-loop_x*(mcae_kmesh(2)&
                   *mcae_kmesh(3)))/mcae_kmesh(3)
          loop_z = loop_xyz-loop_x*(mcae_kmesh(2)*mcae_kmesh(3))&
                   -loop_y*mcae_kmesh(3)

          ! Internally, I need the relative (fractional) coordinates in units of the reciprocal-lattice vectors
          kpoints(1,i) = loop_x*db1
          kpoints(2,i) = loop_y*db2
          kpoints(3,i) = loop_z*db3

       end do

    end if
    call comms_bcast(kweight,1)

    if (on_root) then
        write(stdout, '(a)') 'kpoints done'
    end if
    if (on_root .and. (timing_level>0)) call io_stopwatch('mcae_main: kpoints',2)

    ! Now, I distribute the kpoints; 3* because I send kx, ky, kz
    call comms_scatterv(localkpoints,3*counts(my_node_id),kpoints,3*counts, 3*displs)
    ! Allocate at least one entry, even if we don't use it
    !allocate(localkpointidx(max(1,counts(my_node_id))),stat=ierr)
    !if (ierr/=0) call io_error('Error allocating localkpointidx in mcae_main.')
    !call comms_scatterv(localkpointidx,counts(my_node_id),kpointidx,counts, displs)

    if (on_root .and. (timing_level>0)) call io_stopwatch('mcae_main: sum_k',1)
    ! And now, each node calculates its own k points
    write(outfile_name_node, '(a,i3)') 'eigs',102+my_node_id
    open(unit=102+my_node_id, file=outfile_name_node, form='formatted')

    ! first we interpolate eigenvalues
    call interp_eig()
    if (on_root) then
      write(stdout, '(a)') 'interp_eig done'
    end if

    ! then get fermi_energy
    ef = get_ef(localeig, num_wann, counts(my_node_id), mcae_num_elec, mcae_smr_fixed_en_width, mcae_smr_index)
    if (on_root) then
      write(stdout, '(a)') 'get_ef done'
    end if

    ! finally get the eband
    do i=1, counts(my_node_id)
       localsum(i) = sum_k(i)
    end do

    close(102+my_node_id)
    if (on_root .and. (timing_level>0)) call io_stopwatch('mcae_main: sum_k',2)

    ! Now, I get the results from the different nodes
    call comms_gatherv(localsum,counts(my_node_id),globalsum, &
         counts, displs)
    if (on_root) then
        write(stdout, '(a)') 'comms_gatherv done'
    end if


    ! interpolations on adaptive mesh should be done when all the
    ! interpolations on initial kmesh have finished
    !if(rdum>mcae_adpt_kmesh_thresh) then
    !    adpt_counter=adpt_counter+1
    !    do loop_adpt=1,mcae_adpt_kmesh**3
    !        ! Using imf_k_list here would corrupt values for other
    !        ! frequencies, hence dummy. Only if-th element is used
    !        call berry_get_imf_klist(kpt(:)+adkpt(:,loop_adpt),&
    !                imf_k_list_dummy)
    !        imf_list(:,:,if)=imf_list(:,:,if)&
    !                +imf_k_list_dummy(:,:,if)*kweight_adpt
    !    end do
    !endif



    ! Collect contributions from all nodes
    !
    !call comms_reduce(imf_list(1,1,1),3*3*nfermi,'SUM')
    !call comms_reduce(adpt_counter_list(1),nfermi,'SUM')

    ! All k points processed: Final processing/deallocations

    call internal_post_info()

    if (allocated(kpoints)) deallocate(kpoints)
    if (allocated(globalsum)) deallocate(globalsum)

    !if (allocated(kpointidx)) deallocate(kpointidx)
    !if (allocated(localkpointidx)) deallocate(localkpointidx)
    if (allocated(localkpoints)) deallocate(localkpoints)
    !if (allocated(delHH)) deallocate(delHH)
    if (allocated(localeig)) deallocate(localeig)
    if (allocated(localdel_eig)) deallocate(localdel_eig)
    if (allocated(locallevelspacing)) deallocate(locallevelspacing)
    !if (allocated(globaleig)) deallocate(globaleig)
    !if (allocated(globaldeleig)) deallocate(globaldeleig)
    if (allocated(localsum)) deallocate(localsum)

    if(on_root .and. (timing_level>0)) call io_stopwatch('mcae_main',2)

    return

  contains

    subroutine interp_eig()
    ! interpolate eigenvalues on this proccess
      use w90_utility, only        : utility_diagonalize
      use w90_postw90_common, only : pw90common_fourier_R_to_k
      use w90_wan_ham, only        : wham_get_eig_deleig
      use w90_dos, only            : dos_get_levelspacing

      implicit none

      complex(kind=dp), dimension(:,:), allocatable   :: HH
      complex(kind=dp), dimension(:,:), allocatable   :: UU
      complex(kind=dp), dimension(:,:,:), allocatable :: delHH
      integer                                         :: i
      real(kind=dp), dimension(3)                     :: kpt
      !
      allocate(HH(num_wann,num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating HH in mcae_interp')
      allocate(UU(num_wann,num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating UU in mcae_interp')
      allocate(delHH(num_wann,num_wann,3),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating delHH in dos')

      do i=1, counts(my_node_id)
        kpt = localkpoints(:,i)
        if (mcae_adpt_smr) then
          call wham_get_eig_deleig(kpt,localeig(:,i),localdel_eig(:,:,i),HH,delHH,UU)
          if (on_root) then
            write(stdout, '(a)') 'wham_get_eig_deleig done'
          end if
          call dos_get_levelspacing(localdel_eig(:,:,i),mcae_kmesh,locallevelspacing(:,i))
          if (on_root) then
            write(stdout, '(a)') 'dos_get_levelspacing done'
          end if
        else
          ! Here I get the band energies
          call pw90common_fourier_R_to_k(kpt,HH_R,HH,0)
          call utility_diagonalize(HH,num_wann,localeig(:,i),UU)
        end if
      end do

      if (allocated(HH)) deallocate(HH)
      if (allocated(UU)) deallocate(UU)
      if (allocated(delHH)) deallocate(delHH)

      return

    end subroutine interp_eig

    function sum_k(ikpt)
    ! compute band energy on the kpoint ikpt
      implicit none

      integer, intent(in)                             :: ikpt
      real(kind=dp), dimension(3)                     :: kpt
      real(kind=dp)                                   :: sum_k
      real(kind=dp), dimension(:), allocatable        :: occ

      !
      kpt = localkpoints(:,ikpt)

      allocate(occ(num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error allocating occ in mcae_interp::sum_k.')

      write(*, '(a, i4, a, 3(f10.7,1x))') 'on node ', my_node_id, ', kpt ', kpt

      ! get occupancies
      occ = get_occ(ikpt)
      !call pw90common_get_occ(localeig(:,i),occ,fermi_energy)
      write(102+my_node_id, '(3(f10.7))') kpt(1:3)
      write(102+my_node_id,'((f14.7))') localeig(:,ikpt)
      write(102+my_node_id,'((f14.7))') occ(:)

      ! local sum
      sum_k = dot_product(occ, localeig(:,ikpt)) * kweight

      if (allocated(occ)) deallocate(occ)

      write(*, '(a, i4, a, 3(f10.7,1x), a)') 'on node ', my_node_id, ', kpt ', &
              kpt, ' done'

      return

    end function sum_k

    function get_occ(ikpt)
    ! get occupancy for an array of eigenvalues
      use w90_utility, only        : utility_wgauss
      use w90_postw90_common, only : pw90common_get_occ

      implicit none

      real(kind=dp) :: get_occ(num_wann)
      real(kind=dp) :: arg, eta_smr
      integer       :: ikpt, i

      if (mcae_no_smr) then
        call pw90common_get_occ(localeig(:,ikpt),get_occ,ef)
      else
        do i=1,num_wann
          if (mcae_adpt_smr) then
            ! Eq.(35) YWVS07
            eta_smr=min(locallevelspacing(i,ikpt)*mcae_adpt_smr_fac,mcae_adpt_smr_max)
          else
            eta_smr=mcae_smr_fixed_en_width
          endif
          arg = (ef - localeig(i,ikpt))/eta_smr
          get_occ(i)=utility_wgauss(arg,mcae_smr_index)
        end do
      end if

      return

    end function get_occ

    function get_ef(eig, nbnd, nks, nelec, degauss, ngauss)
      !--------------------------------------------------------------------
      !
      !     Finds the Fermi energy - Gaussian Broadening
      !     (see Methfessel and Paxton, PRB 40, 3616 (1989 )
      !

      implicit none
      !  I/O variables
      integer, intent(in) :: nks, nbnd, ngauss
      real(dp), intent(in) :: eig(nbnd, nks), degauss, nelec
      real(dp) :: get_ef
      !
      real(DP), parameter :: eps = 1.0d-10
      integer, parameter :: maxiter = 300
      ! internal variables
      real(dp) :: Ef, Eup, Elw, sumkup, sumklw, sumkmid
      integer :: i, kpoint
      !
      !      find bounds for the Fermi energy. Very safe choice!
      !
      Elw = eig(1, 1)
      Eup = eig(nbnd, 1)
      do kpoint = 2, nks
          Elw = min(Elw, eig(1, kpoint) )
          Eup = max(Eup, eig(nbnd, kpoint) )
      enddo
      Eup = Eup + 2 * degauss
      Elw = Elw - 2 * degauss
      !
      ! find min and max across procs
      !
      call comms_allreduce(Eup, 1, 'MAX')
      call comms_allreduce(Elw, 1, 'MIN')
      !
      !      Bisection method
      !
      sumkup = sum_num_states(eig, nbnd, nks, degauss, ngauss, Eup)
      sumklw = sum_num_states(eig, nbnd, nks, degauss, ngauss, Elw)

      if (on_root) then
        if ( (sumkup - nelec) < -eps .or. (sumklw - nelec) > eps )  &
          call io_error('get_ef: internal error, cannot bracket Ef')
      end if

      do i = 1, maxiter
        Ef = (Eup + Elw) / 2.d0
        sumkmid = sum_num_states(eig, nbnd, nks, degauss, ngauss, Ef)
        if (abs (sumkmid-nelec) < eps) then
          get_ef = Ef
          return
        elseif ( (sumkmid-nelec) < -eps) then
          Elw = Ef
        else
          Eup = Ef
        endif
      enddo

      if (on_root) then
        write(stdout, '(5x,"Warning: too many iterations in bisection"/ &
              &      5x,"Ef = ",f10.6," sumk = ",f10.6," electrons")' ) &
              Ef, sumkmid
      end if
      !
      get_ef = Ef

      return

    end function get_ef

    function sum_num_states(eig, nbnd, nks, degauss, ngauss, e)
      !-----------------------------------------------------------------------
      !
      !     This function computes the number of states under a given energy e
      !
      !
      ! function which compute the smearing
      use w90_utility, only        : utility_wgauss
      use w90_postw90_common, only : pw90common_get_occ

      implicit none

      ! Output variable
      real(dp) :: sum_num_states
      ! Input variables
      integer, intent(in) :: nks, nbnd, ngauss
      ! input: the total number of K points
      ! input: the number of bands
      ! input: the type of smearing
      real(DP), intent(in) :: eig(nbnd, nks), degauss, e
      ! input: the weight of the k points
      ! input: the energy eigenvalues
      ! input: gaussian broadening
      ! input: the energy to check
      !
      ! local variables
      !
      real(DP) :: sum1, eta_smr
      integer  :: ik, ibnd
      real(dp) :: occ(nbnd)
      ! counter on k points
      ! counter on the band energy
      !
      sum_num_states = 0.d0
      do ik = 1, nks
        sum1 = 0.d0
        if (mcae_no_smr) then
          call pw90common_get_occ(eig(:,ik), occ, e)
          sum1 = sum(occ)
        else
          do ibnd = 1, nbnd
            if(mcae_adpt_smr) then
              ! Eq.(35) YWVS07
              eta_smr=min(locallevelspacing(ibnd,ik)*mcae_adpt_smr_fac,mcae_adpt_smr_max)
            else
              eta_smr=degauss
            endif
            sum1 = sum1 + utility_wgauss( (e-eig(ibnd, ik) ) / eta_smr, ngauss)
          enddo
        end if
        sum_num_states = sum_num_states + kweight * sum1
      enddo

      ! result stored on all procs
      call comms_allreduce(sum_num_states, 1, 'SUM')

      return

    end function sum_num_states

  end subroutine mcae_main

end module w90_mcae
