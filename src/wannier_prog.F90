!-*- mode: F90 -*-!
!------------------------------------------------------------!
!                                                            !
!                       WANNIER90                            !
!                                                            !
!          The Maximally-Localised Generalised               !
!                 Wannier Functions Code                     !
!                                                            !
! Please cite                                                !
!                                                            !
!  [ref] "Wannier90 as a community code:                     !
!        new features and applications",                     !
!        G. Pizzi et al.,  J. Phys. Cond. Matt. 32,          !
!        165902 (2020).                                      !
!        http://doi.org/10.1088/1361-648X/ab51ff             !
!                                                            !
! in any publications arising from the use of this code.     !
!                                                            !
! Wannier90 is based on Wannier77, written by N. Marzari,    !
! I. Souza and D. Vanderbilt. For the method please cite     !
!                                                            !
! [ref] N. Marzari and D. Vanderbilt,                        !
!       Phys. Rev. B 56 12847 (1997)                         !
!       http://dx.doi.org/10.1103/PhysRevB.56.12847          !
!                                                            !
! [ref] I. Souza, N. Marzari and D. Vanderbilt,              !
!       Phys. Rev. B 65 035109 (2001)                        !
!       http://dx.doi.org/10.1103/PhysRevB.65.035109         !
!                                                            !
! [ref] N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza,    !
!       D. Vanderbilt, "Maximally localized Wannier          !
!       functions: theory and applications",                 !
!       Rev. Mod. Phys. 84, 1419 (2012)                      !
!       http://dx.doi.org/10.1103/RevModPhys.84.1419         !
!                                                            !
! For a full list of authors and contributors, please        !
! see the README file in the root directory of the           !
! distribution.                                              !
!                                                            !
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

program wannier
  !! The main Wannier90 program

  use w90_constants
  use w90_parameters
  use w90_io
  use w90_hamiltonian
  use w90_kmesh
  use w90_disentangle
  use w90_overlap
  use w90_wannierise
  use w90_plot
  use w90_transport
  use w90_comms, only: on_root, num_nodes, comms_setup, comms_end, comms_bcast, my_node_id
  use w90_sitesym !YN:

  implicit none

  real(kind=dp) time0, time1, time2
  character(len=9) :: stat, pos, cdate, ctime
  logical :: wout_found, dryrun
  integer :: len_seedname
  character(len=50) :: prog

  call comms_setup

  library = .false.

  time0 = io_time()

  if (on_root) then
    prog = 'wannier90'
    call io_commandline(prog, dryrun)
    len_seedname = len(seedname)
  end if
  call comms_bcast(len_seedname, 1)
  call comms_bcast(seedname, len_seedname)
  call comms_bcast(dryrun, 1)

  if (on_root) then
    stdout = io_file_unit()
    open (unit=stdout, file=trim(seedname)//'.werr')
    call io_date(cdate, ctime)
    write (stdout, *) 'Wannier90: Execution started on ', cdate, ' at ', ctime

    call param_read
    close (stdout, status='delete')

    if (restart .eq. ' ') then
      stat = 'replace'
      pos = 'rewind'
    else
      inquire (file=trim(seedname)//'.wout', exist=wout_found)
      if (wout_found) then
        stat = 'old'
      else
        stat = 'replace'
      endif
      pos = 'append'
    endif

    stdout = io_file_unit()
    open (unit=stdout, file=trim(seedname)//'.wout', status=trim(stat), position=trim(pos))
    call param_write_header()
    if (num_nodes == 1) then
#ifdef MPI
      write (stdout, '(/,1x,a)') 'Running in serial (with parallel executable)'
#else
      write (stdout, '(/,1x,a)') 'Running in serial (with serial executable)'
#endif
    else
      write (stdout, '(/,1x,a,i3,a/)') &
        'Running in parallel on ', num_nodes, ' CPUs'
    endif
    call param_write()

    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'

    if (.not. explicit_nnkpts) call kmesh_get
    time2 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') &
      'Time to get kmesh        ', time2 - time1, ' (sec)'

    call param_memory_estimate
  end if

  if (dryrun) then
    if (on_root) then
      write (stdout, *) ' '
      write (stdout, *) '                       ==============================='
      write (stdout, *) '                                   DRYRUN             '
      write (stdout, *) '                       No problems found with win file'
      write (stdout, *) '                       ==============================='
    endif
    stop
  endif

  ! We now distribute the parameters to the other nodes
  call param_dist
  if (gamma_only .and. num_nodes > 1) &
    call io_error('Gamma point branch is serial only at the moment')

  if (transport .and. tran_read_ht) goto 3003

  ! Sort out restarts
  if (restart .eq. ' ') then  ! start a fresh calculation
    if (on_root) write (stdout, '(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else                      ! restart a previous calculation
    if (on_root) call param_read_chkpt()
    call param_chkpt_dist
    if (lsitesymmetry) call sitesym_read()   ! update this to read on root and bcast - JRY

    select case (restart)
    case ('default')    ! continue from where last checkpoint was written
      if (on_root) write (stdout, '(/1x,a)', advance='no') 'Resuming a previous Wannier90 calculation '
      if (checkpoint .eq. 'postdis') then
        if (on_root) write (stdout, '(a/)') 'from wannierisation ...'
        goto 1001         ! go to wann_main
      elseif (checkpoint .eq. 'postwann') then
        if (on_root) write (stdout, '(a/)') 'from plotting ...'
        goto 2002         ! go to plot_main
      else
        if (on_root) write (stdout, '(/a/)')
        call io_error('Value of checkpoint not recognised in wann_prog')
      endif
    case ('wannierise') ! continue from wann_main irrespective of value of last checkpoint
      if (on_root) write (stdout, '(1x,a/)') 'Restarting Wannier90 from wannierisation ...'
      goto 1001
    case ('plot')       ! continue from plot_main irrespective of value of last checkpoint
      if (on_root) write (stdout, '(1x,a/)') 'Restarting Wannier90 from plotting routines ...'
      goto 2002
    case ('transport')   ! continue from tran_main irrespective of value of last checkpoint
      if (on_root) write (stdout, '(1x,a/)') 'Restarting Wannier90 from transport routines ...'
      goto 3003
    case default        ! for completeness... (it is already trapped in param_read)
      call io_error('Value of restart not recognised in wann_prog')
    end select
  endif

  if (postproc_setup) then
    if (on_root) call kmesh_write()
    call kmesh_dealloc()
    call param_dealloc()
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'
    if (on_root) write (stdout, '(/a)') ' Exiting... '//trim(seedname)//'.nnkp written.'
    call comms_end
    stop
  endif

  if (lsitesymmetry) call sitesym_read()   ! update this to read on root and bcast - JRY
  call overlap_allocate()
  call overlap_read()

  time1 = io_time()
  if (on_root) write (stdout, '(/1x,a25,f11.3,a)') 'Time to read overlaps    ', time1 - time2, ' (sec)'

  have_disentangled = .false.

  if (disentanglement) then
    call dis_main()
    have_disentangled = .true.
    time2 = io_time()
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to disentangle bands', time2 - time1, ' (sec)'
  endif

  if (on_root) call param_write_chkpt('postdis')
!~  call param_write_um

1001 time2 = io_time()

  if (.not. gamma_only) then
    call wann_main()
  else
    call wann_main_gamma()
  end if

  time1 = io_time()
  if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time for wannierise      ', time1 - time2, ' (sec)'

  if (on_root) call param_write_chkpt('postwann')

2002 continue
!!!!!!!!!!!!!!!!!!!
call print_spread()
!!!!!!!!!!!!!!!!!!!
  if (on_root) then
    ! I call the routine always; the if statements to decide if/what
    ! to plot are inside the function
    time2 = io_time()
  endif
  call plot_main()
  if (on_root) then
    time1 = io_time()
    ! Now time is always printed, even if no plotting is done/required, but
    ! it shouldn't be a problem.
    write (stdout, '(1x,a25,f11.3,a)') 'Time for plotting        ', time1 - time2, ' (sec)'
  endif

3003 continue
  if (on_root) then
    time2 = io_time()
    if (transport) then
      call tran_main()
      time1 = io_time()
      write (stdout, '(1x,a25,f11.3,a)') 'Time for transport       ', time1 - time2, ' (sec)'
      if (tran_read_ht) goto 4004
    end if
  endif

  call tran_dealloc()
  call hamiltonian_dealloc()
  call overlap_dealloc()
  call kmesh_dealloc()
  call param_dealloc()
  if (lsitesymmetry) call sitesym_dealloc() !YN:

4004 continue

  if (on_root) then
    write (stdout, '(1x,a25,f11.3,a)') 'Total Execution Time     ', io_time(), ' (sec)'

    if (timing_level > 0) call io_print_timings()

    write (stdout, *)
    write (stdout, '(1x,a)') 'All done: wannier90 exiting'

    close (stdout)
  endif

  call comms_end

contains
  subroutine print_spread
    use w90_constants, only: dp, cmplx_1, cmplx_0, eps2, eps5, eps8
    use w90_io, only: stdout, io_error, io_wallclocktime, io_stopwatch &
      , io_file_unit
    use w90_parameters, only: num_wann, num_cg_steps, num_iter, nnlist, &
      nntot, wbtot, u_matrix, m_matrix, num_kpts, iprint, num_print_cycles, &
      num_dump_cycles, omega_invariant, param_write_chkpt, length_unit, &
      lenconfac, proj_site, real_lattice, write_r2mn, guiding_centres, &
      num_guide_cycles, num_no_guide_iter, timing_level, trial_step, precond, spinors, &
      fixed_step, lfixstep, write_proj, have_disentangled, conv_tol, num_proj, &
      conv_window, conv_noise_amp, conv_noise_num, wannier_centres, write_xyz, &
      wannier_spreads, omega_total, omega_tilde, optimisation, write_vdw_data, &
      write_hr_diag, kpt_latt, bk, ccentres_cart, slwf_num, selective_loc, &
      slwf_constrain, slwf_lambda
    use w90_utility, only: utility_frac_to_cart, utility_zgemm
    use w90_parameters, only: lsitesymmetry                !RS:
    use w90_sitesym, only: sitesym_symmetrize_gradient  !RS:

    implicit none

    ! guiding centres
    real(kind=dp), allocatable :: rguide(:, :)
    integer :: irguide

    ! local arrays used and passed in subroutines
    complex(kind=dp), allocatable :: csheet(:, :, :)
    complex(kind=dp), allocatable :: cdodq(:, :, :)
    complex(kind=dp), allocatable :: cdodq_r(:, :, :)
    complex(kind=dp), allocatable :: k_to_r(:, :)
    complex(kind=dp), allocatable :: cdodq_precond(:, :, :)
    complex(kind=dp), allocatable :: cdodq_precond_loc(:, :, :)
    real(kind=dp), allocatable :: sheet(:, :, :)
    real(kind=dp), allocatable :: rave(:, :), r2ave(:), rave2(:)
    real(kind=dp), dimension(3) :: rvec_cart

    real(kind=dp) :: lambda_loc
    type(localisation_vars) :: old_spread
    type(localisation_vars) :: wann_spread
    logical       :: lquad
    integer       :: i, n, iter, ind, ierr, iw, ncg, info, nkp, nkp_loc, nn
    !
    irguide = 0
    if (guiding_centres .and. (num_no_guide_iter .le. 0)) then
      call wann_phases(csheet, sheet, rguide, irguide)
      irguide = 1
    endif

    ! constrained centres part
    lambda_loc = 0.0_dp
    if (selective_loc .and. slwf_constrain) then
      lambda_loc = slwf_lambda
    end if

    ! calculate initial centers and spread
    call wann_omega(csheet, sheet, rave, r2ave, rave2, wann_spread)

    ! public variables
    if (.not. selective_loc) then
      omega_total = wann_spread%om_tot
      omega_invariant = wann_spread%om_i
      omega_tilde = wann_spread%om_d + wann_spread%om_od
    else
      omega_total = wann_spread%om_tot
      ! omega_invariant = wann_spread%om_iod
      ! omega_tilde = wann_spread%om_d + wann_spread%om_nu
    end if

    ! public arrays of Wannier centres and spreads
    wannier_centres = rave
    wannier_spreads = r2ave - rave2

    if (lfixstep) lquad = .false.
    ncg = 0
    iter = 0
    old_spread%om_tot = 0.0_dp

    ! print initial state
    if (on_root) then
      write (stdout, '(1x,a78)') repeat('-', 78)
      write (stdout, '(1x,a)') 'Initial State'
      do iw = 1, num_wann
        write (stdout, 1000) iw, (rave(ind, iw)*lenconfac, ind=1, 3), &
          (r2ave(iw) - rave2(iw))*lenconfac**2
      end do
      write (stdout, 1001) (sum(rave(ind, :))*lenconfac, ind=1, 3), (sum(r2ave) - sum(rave2))*lenconfac**2
      write (stdout, *)
      if (selective_loc .and. slwf_constrain) then
        write (stdout, '(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)') &
          iter, (wann_spread%om_tot - old_spread%om_tot)*lenconfac**2, sqrt(abs(gcnorm1))*lenconfac, &
          wann_spread%om_tot*lenconfac**2, io_wallclocktime(), '<-- CONV'
        write (stdout, '(7x,a,F15.7,a,F15.7,a,F15.7,a,F15.7,a)') &
          'O_D=', wann_spread%om_d*lenconfac**2, &
          ' O_IOD=', (wann_spread%om_iod + wann_spread%om_nu)*lenconfac**2, &
          ' O_TOT=', wann_spread%om_tot*lenconfac**2, ' <-- SPRD'
        write (stdout, '(1x,a78)') repeat('-', 78)
      elseif (selective_loc .and. .not. slwf_constrain) then
        write (stdout, '(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)') &
          iter, (wann_spread%om_tot - old_spread%om_tot)*lenconfac**2, sqrt(abs(gcnorm1))*lenconfac, &
          wann_spread%om_tot*lenconfac**2, io_wallclocktime(), '<-- CONV'
        write (stdout, '(7x,a,F15.7,a,F15.7,a,F15.7,a)') &
          'O_D=', wann_spread%om_d*lenconfac**2, &
          ' O_IOD=', wann_spread%om_iod*lenconfac**2, &
          ' O_TOT=', wann_spread%om_tot*lenconfac**2, ' <-- SPRD'
        write (stdout, '(1x,a78)') repeat('-', 78)
      else
        write (stdout, '(1x,i6,2x,E12.3,2x,F15.10,2x,F18.10,3x,F8.2,2x,a)') &
          iter, (wann_spread%om_tot - old_spread%om_tot)*lenconfac**2, sqrt(abs(gcnorm1))*lenconfac, &
          wann_spread%om_tot*lenconfac**2, io_wallclocktime(), '<-- CONV'
        write (stdout, '(8x,a,F15.7,a,F15.7,a,F15.7,a)') &
          'O_D=', wann_spread%om_d*lenconfac**2, ' O_OD=', wann_spread%om_od*lenconfac**2, &
          ' O_TOT=', wann_spread%om_tot*lenconfac**2, ' <-- SPRD'
        write (stdout, '(1x,a78)') repeat('-', 78)
      end if
    endif
  end subroutine print_spread
end program wannier
