!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!------------------------------------------------------------!

!===============================================================!
!                                                               !
!  w90_mrwf : manifold-remixed Wannier functions                !
!                                                               !
!  Split a valence+conduction Wannierisation into isolated      !
!  energy manifolds by diagonalising the Wannier-gauge          !
!  Hamiltonian, then parallel-transport each manifold to a      !
!  smooth gauge. Writes per-manifold mmn/eig/amn/win files so   !
!  that each manifold becomes a standalone Wannier problem.     !
!                                                               !
!  Serial port of the Wannier.jl mrwf/split implementation.     !
!                                                               !
!===============================================================!

module w90_mrwf_mod

#ifdef MPI08
  use mpi_f08
#endif
#ifdef MPI90
  use mpi
#endif

  use w90_constants, only: dp, cmplx_0, cmplx_1
  use w90_types, only: kmesh_info_type, atom_data_type, dis_manifold_type, &
                       ws_region_type, print_output_type, wannier_data_type, timer_list_type
  use w90_error_base, only: w90_error_type
  use w90_error, only: set_error_input, set_error_fatal, set_error_alloc, set_error_file
  use w90_comms, only: w90_comm_type, mpirank
  use w90_utility, only: utility_zgemm_new, utility_diagonalize
  use w90_parallel_transport_mod, only: w90_parallel_transport

#ifdef MPIH
  include 'mpif.h'
#endif

  implicit none

  private

  public :: w90_mrwf

contains

  !> Split the Wannierisation into the isolated manifolds `manifolds(1:2, g)`
  !> (first:last band index into the num_wann WFs), parallel-transport each,
  !> and write the per-manifold files into `outdirs(g)`.
  subroutine w90_mrwf(kmesh_info, u_matrix, u_matrix_opt, m_matrix, eigval, &
                      kpt_latt, real_lattice, atom_data, mp_grid, &
                      num_bands, num_wann, num_kpts, manifolds, outdirs, &
                      log_interp, run_maxloc, write_unk, num_iter, dis_manifold, &
                      have_disentangled, wvfn_formatted, wvfn_spin, &
                      seedname, stdout, error, comm)
    type(kmesh_info_type), intent(in) :: kmesh_info
    complex(kind=dp), intent(in) :: u_matrix(:, :, :)      ! (nw, nw, nk)
    complex(kind=dp), intent(in) :: u_matrix_opt(:, :, :)  ! (nb, nw, nk)
    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)   ! Wannier gauge (nw, nw, nn, nk)
    real(kind=dp), intent(in) :: eigval(:, :)              ! (nb, nk)
    real(kind=dp), intent(in) :: kpt_latt(:, :)            ! (3, nk)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    type(atom_data_type), intent(in) :: atom_data
    integer, intent(in) :: mp_grid(3), num_bands, num_wann, num_kpts, stdout
    integer, intent(in) :: manifolds(:, :)                 ! (2, n_manifold)
    character(len=*), intent(in) :: outdirs(:)
    logical, intent(in) :: log_interp, run_maxloc, write_unk, have_disentangled, wvfn_formatted
    integer, intent(in) :: num_iter, wvfn_spin
    type(dis_manifold_type), intent(in) :: dis_manifold
    character(len=*), intent(in) :: seedname
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    complex(kind=dp), allocatable :: utot(:, :, :), hw(:, :), vmat(:, :, :), tmp(:, :)
    complex(kind=dp), allocatable :: mg(:, :, :, :), upt(:, :, :), vg(:, :), vg2(:, :)
    complex(kind=dp), allocatable :: gsplit(:, :, :), tt(:, :), bg(:, :, :)
    real(kind=dp), allocatable :: dw(:, :), evk(:)
    integer :: rank, ik, ik2, inn, ig, a, b, ng, nman, i, m1
    character(len=256) :: gseed

    rank = mpirank(comm)
    nman = size(manifolds, 2)

    ! --- Wannier-gauge Hamiltonian H^W(k) = U_tot^H diag(E) U_tot, diagonalised ---
    allocate (utot(num_bands, num_wann, num_kpts), hw(num_wann, num_wann), &
              vmat(num_wann, num_wann, num_kpts), dw(num_wann, num_kpts), &
              tmp(num_bands, num_wann), evk(num_wann))
    do ik = 1, num_kpts
      ! U_tot = u_opt . u   (nb x nw)
      call utility_zgemm_new(u_matrix_opt(:, :, ik), u_matrix(:, :, ik), utot(:, :, ik), 'N', 'N')
      ! tmp = diag(E) . U_tot
      do m1 = 1, num_bands
        tmp(m1, :) = eigval(m1, ik)*utot(m1, :, ik)
      end do
      ! H^W = U_tot^H . tmp
      call utility_zgemm_new(utot(:, :, ik), tmp, hw, 'C', 'N')
      call utility_diagonalize(hw, num_wann, evk, vmat(:, :, ik), error, comm)
      if (allocated(error)) return
      dw(:, ik) = evk
    end do
    deallocate (hw, tmp, evk)

    ! --- per manifold: subgroup overlaps, parallel transport, assemble, write ---
    do ig = 1, nman
      a = manifolds(1, ig)
      b = manifolds(2, ig)
      ng = b - a + 1
      if (rank == 0 .and. stdout > 0) then
        write (stdout, '(/,1x,a,i0,a,i0,a,i0,a)') 'mrwf: manifold ', ig, ' = bands ', a, ':', b, &
          ' -> parallel transport'
      end if

      allocate (mg(ng, ng, kmesh_info%nntot, num_kpts), upt(ng, ng, num_kpts), &
                vg(num_wann, ng), vg2(num_wann, ng), tt(ng, num_wann))
      ! M^g(k, nn) = V_g(k)^H . M_wann(:,:,nn,k) . V_g(k2)
      do ik = 1, num_kpts
        vg = vmat(:, a:b, ik)
        do inn = 1, kmesh_info%nntot
          ik2 = kmesh_info%nnlist(ik, inn)
          vg2 = vmat(:, a:b, ik2)
          call utility_zgemm_new(vg, m_matrix(:, :, inn, ik), tt, 'C', 'N') ! tt = V_g^H M
          call utility_zgemm_new(tt, vg2, mg(:, :, inn, ik), 'N', 'N')      ! M^g = tt V_g2
        end do
      end do

      ! seed identity gauge, parallel transport the manifold
      upt = cmplx_0
      do ik = 1, num_kpts
        do i = 1, ng
          upt(i, i, ik) = cmplx_1
        end do
      end do
      call w90_parallel_transport(kmesh_info, upt, mg, kpt_latt, mp_grid, ng, num_kpts, &
                                  .false., log_interp, stdout, error, comm)
      if (allocated(error)) return

      ! optionally run a full maximal localisation on the manifold; upt (gauge) and
      ! mg (overlaps) are updated in place, so the writes below pick up the MLWF gauge.
      if (run_maxloc) then
        if (rank == 0 .and. stdout > 0) write (stdout, '(1x,a,i0,a)') &
          'mrwf: manifold ', ig, ' -> maximal localisation'
        call mrwf_maxloc(kmesh_info, mg, upt, kpt_latt, real_lattice, mp_grid, &
                         num_kpts, ng, num_iter, stdout, error, comm)
        if (allocated(error)) return
      end if

      ! assembled split gauge G_g(k) = U_tot(k) . V_g(k) . U_pt_g(k)   (nb x ng)
      allocate (gsplit(num_bands, ng, num_kpts))
      do ik = 1, num_kpts
        block
          complex(kind=dp) :: uv(num_bands, ng)
          call utility_zgemm_new(utot(:, :, ik), vmat(:, a:b, ik), uv, 'N', 'N') ! U_tot V_g
          call utility_zgemm_new(uv, upt(:, :, ik), gsplit(:, :, ik), 'N', 'N')  ! . U_pt
        end block
      end do

      ! --- write per-manifold files (root only) ---
      if (rank == 0) then
        gseed = trim(outdirs(ig))//'/'//trim(seedname)
        call execute_command_line('mkdir -p '//trim(outdirs(ig)))
        call mrwf_write_mmn(trim(gseed)//'.mmn', mg, kmesh_info, ng, num_kpts)
        call mrwf_write_eig(trim(gseed)//'.eig', dw(a:b, :), ng, num_kpts)
        call mrwf_write_amn(trim(gseed)//'.amn', upt, ng, ng, num_kpts)
        call mrwf_write_amn(trim(gseed)//'_split.amn', gsplit, num_bands, ng, num_kpts)
        call mrwf_write_win(trim(gseed)//'.win', ng, mp_grid, real_lattice, atom_data, &
                            kpt_latt, num_kpts)
        write (stdout, '(1x,a,i0,a)') 'mrwf: manifold ', ig, ' written to '//trim(outdirs(ig))

        ! optionally rotate the UNK files into this manifold for WF plotting
        if (write_unk) then
          ! post-u_opt gauge B_g(k) = u(k) . V_g(k) . U_final_g(k)   (nw x ng)
          allocate (bg(num_wann, ng, num_kpts))
          do ik = 1, num_kpts
            block
              complex(kind=dp) :: uvg(num_wann, ng)
              call utility_zgemm_new(u_matrix(:, :, ik), vmat(:, a:b, ik), uvg, 'N', 'N')
              call utility_zgemm_new(uvg, upt(:, :, ik), bg(:, :, ik), 'N', 'N')
            end block
          end do
          call mrwf_write_unk_files(trim(outdirs(ig)), bg, u_matrix_opt, dis_manifold, &
                                    have_disentangled, num_bands, num_wann, ng, num_kpts, &
                                    wvfn_formatted, wvfn_spin, stdout, error, comm)
          deallocate (bg)
          if (allocated(error)) return
          write (stdout, '(1x,a,i0)') 'mrwf: UNK files rotated for manifold ', ig
        end if
      end if

      deallocate (mg, upt, vg, vg2, tt, gsplit)
    end do

    deallocate (utot, vmat, dw)
  end subroutine w90_mrwf

  !> Write overlaps in .mmn format.
  subroutine mrwf_write_mmn(fname, mg, kmesh_info, ng, num_kpts)
    character(len=*), intent(in) :: fname
    complex(kind=dp), intent(in) :: mg(:, :, :, :)
    type(kmesh_info_type), intent(in) :: kmesh_info
    integer, intent(in) :: ng, num_kpts
    integer :: iun, ik, inn, i, j
    open (newunit=iun, file=fname, form='formatted', action='write')
    write (iun, '(a)') 'mrwf manifold overlaps'
    write (iun, '(3i12)') ng, num_kpts, kmesh_info%nntot
    do ik = 1, num_kpts
      do inn = 1, kmesh_info%nntot
        write (iun, '(2i8,3i5)') ik, kmesh_info%nnlist(ik, inn), kmesh_info%nncell(:, ik, inn)
        do j = 1, ng
          do i = 1, ng
            write (iun, '(2f18.12)') real(mg(i, j, inn, ik), dp), aimag(mg(i, j, inn, ik))
          end do
        end do
      end do
    end do
    close (iun)
  end subroutine mrwf_write_mmn

  !> Write eigenvalues in .eig format.
  subroutine mrwf_write_eig(fname, dg, ng, num_kpts)
    character(len=*), intent(in) :: fname
    real(kind=dp), intent(in) :: dg(:, :)
    integer, intent(in) :: ng, num_kpts
    integer :: iun, ik, i
    open (newunit=iun, file=fname, form='formatted', action='write')
    do ik = 1, num_kpts
      do i = 1, ng
        write (iun, '(2i5,f18.12)') i, ik, dg(i, ik)
      end do
    end do
    close (iun)
  end subroutine mrwf_write_eig

  !> Write a gauge matrix (nrow x ncol x nk) in .amn format.
  subroutine mrwf_write_amn(fname, a, nrow, ncol, num_kpts)
    character(len=*), intent(in) :: fname
    complex(kind=dp), intent(in) :: a(:, :, :)
    integer, intent(in) :: nrow, ncol, num_kpts
    integer :: iun, ik, i, j
    open (newunit=iun, file=fname, form='formatted', action='write')
    write (iun, '(a)') 'mrwf gauge'
    write (iun, '(3i12)') nrow, num_kpts, ncol
    do ik = 1, num_kpts
      do j = 1, ncol
        do i = 1, nrow
          write (iun, '(3i5,2f18.12)') i, j, ik, real(a(i, j, ik), dp), aimag(a(i, j, ik))
        end do
      end do
    end do
    close (iun)
  end subroutine mrwf_write_amn

  !> Write a minimal, runnable .win for a manifold.
  subroutine mrwf_write_win(fname, ng, mp_grid, real_lattice, atom_data, kpt_latt, num_kpts)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: ng, mp_grid(3), num_kpts
    real(kind=dp), intent(in) :: real_lattice(3, 3), kpt_latt(:, :)
    type(atom_data_type), intent(in) :: atom_data
    integer :: iun, i, is, ia, ik
    open (newunit=iun, file=fname, form='formatted', action='write')
    write (iun, '(a,i0)') 'num_wann = ', ng
    write (iun, '(a)') 'num_iter = 2000'
    write (iun, '(a)') 'auto_projections = .true.'
    write (iun, '(a)') ''
    write (iun, '(a)') 'begin unit_cell_cart'
    write (iun, '(a)') 'ang'
    do i = 1, 3
      write (iun, '(3f18.12)') real_lattice(:, i)
    end do
    write (iun, '(a)') 'end unit_cell_cart'
    write (iun, '(a)') ''
    write (iun, '(a)') 'begin atoms_cart'
    write (iun, '(a)') 'ang'
    do is = 1, atom_data%num_species
      do ia = 1, atom_data%species_num(is)
        write (iun, '(a,3x,3f18.12)') trim(atom_data%symbol(is)), atom_data%pos_cart(:, ia, is)
      end do
    end do
    write (iun, '(a)') 'end atoms_cart'
    write (iun, '(a)') ''
    write (iun, '(a,3i4)') 'mp_grid = ', mp_grid
    write (iun, '(a)') ''
    write (iun, '(a)') 'begin kpoints'
    do ik = 1, num_kpts
      write (iun, '(3f18.12)') kpt_latt(:, ik)
    end do
    write (iun, '(a)') 'end kpoints'
    close (iun)
  end subroutine mrwf_write_win

  !> Run a full maximal localisation (wann_main) on one manifold. `mg` (overlaps)
  !> and `u` (gauge, seeded with the parallel-transport gauge) are updated in place
  !> to the MLWF gauge. Runs redundantly on every rank using a serial communicator
  !> (MPI_COMM_SELF), so wann_main's allreduce does not multiply the identical data.
  subroutine mrwf_maxloc(kmesh_info, mg, u, kpt_latt, real_lattice, mp_grid, &
                         num_kpts, ng, num_iter, stdout, error, comm)
    use w90_wannierise_mod, only: wann_main
    use w90_wannier90_types, only: wann_control_type, wann_omega_type, sitesym_type, &
                                   w90_calculation_type, ham_logical_type
    type(kmesh_info_type), intent(in) :: kmesh_info
    complex(kind=dp), intent(in) :: mg(:, :, :, :)     ! M^g in the V_g reference basis
    complex(kind=dp), intent(inout) :: u(:, :, :)      ! seed U_pt in, U_final out
    real(kind=dp), intent(in) :: kpt_latt(:, :), real_lattice(3, 3)
    integer, intent(in) :: mp_grid(3), num_kpts, ng, num_iter, stdout
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    type(ham_logical_type) :: hl
    type(wann_control_type) :: wctl
    type(wann_omega_type) :: om
    type(sitesym_type) :: ss
    type(print_output_type) :: po
    type(ws_region_type) :: wsr
    type(w90_calculation_type) :: wcalc
    type(wannier_data_type) :: wd
    type(timer_list_type) :: tmr
    type(w90_comm_type) :: selfcomm
    complex(kind=dp), allocatable :: ham_k(:, :, :), ham_r(:, :, :)
    integer, allocatable :: irvec(:, :), ndegen(:), distk(:)
    real(kind=dp), allocatable :: wct(:, :)
    integer :: nrpts, rpt_origin, ik, inn, ik2
    complex(kind=dp), allocatable :: rt(:, :), mgc(:, :, :, :)

    ! wann_main expects m_matrix in the CURRENT WF gauge, i.e. <w|w_kb> with the
    ! seed gauge `u` already applied, and rotates it in place. Work on a copy so
    ! the caller's mg stays in the V_g reference basis for the .mmn output:
    !   mgc(k,b) = u(k)^H . M^g(k,b) . u(k2).
    allocate (mgc(ng, ng, kmesh_info%nntot, num_kpts), rt(ng, ng))
    do ik = 1, num_kpts
      do inn = 1, kmesh_info%nntot
        ik2 = kmesh_info%nnlist(ik, inn)
        call utility_zgemm_new(u(:, :, ik), mg(:, :, inn, ik), rt, 'C', 'N') ! u^H . M
        call utility_zgemm_new(rt, u(:, :, ik2), mgc(:, :, inn, ik), 'N', 'N')
      end do
    end do
    deallocate (rt)

    wctl%num_iter = num_iter
    allocate (wd%centres(3, ng), wd%spreads(ng))
    wd%centres = 0.0_dp
    wd%spreads = 0.0_dp
    allocate (distk(num_kpts))
    distk = 0  ! all kpoints local to rank 0 of the serial communicator
    nrpts = 0
    rpt_origin = 0
#ifdef MPI
    selfcomm%comm = MPI_COMM_SELF
#else
    selfcomm = comm
#endif

    call wann_main(hl, kmesh_info, kpt_latt, wctl, om, ss, po, wd, wsr, wcalc, &
                   ham_k, ham_r, mgc, u, real_lattice, wct, irvec, mp_grid, ndegen, &
                   nrpts, num_kpts, ng, ng, 3, rpt_origin, 's-k', 'bulk', .false., &
                   stdout, tmr, distk, error, selfcomm)
    deallocate (mgc)
  end subroutine mrwf_maxloc

  !> Rotate the input UNK files into one manifold and write them to `outdir`.
  !> Follows plot_wannier's two-stage convention: stage 1 collapses the
  !> `num_inc` in-window Bloch bands to `num_wann` states via `u_opt`, stage 2
  !> applies the post-u_opt manifold gauge `bg` (nw x ng).
  subroutine mrwf_write_unk_files(outdir, bg, u_opt, dis_manifold, have_dis, &
                                  num_bands, num_wann, ng, num_kpts, formatted, spin, &
                                  stdout, error, comm)
    character(len=*), intent(in) :: outdir
    complex(kind=dp), intent(in) :: bg(:, :, :)      ! (nw, ng, nk)
    complex(kind=dp), intent(in) :: u_opt(:, :, :)   ! (nb, nw, nk)
    type(dis_manifold_type), intent(in) :: dis_manifold
    logical, intent(in) :: have_dis, formatted
    integer, intent(in) :: num_bands, num_wann, ng, num_kpts, spin, stdout
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    complex(kind=dp), allocatable :: wtmp(:, :), rwv(:, :), cwv(:, :), buf(:)
    logical, allocatable :: inc(:)
    character(len=60) :: fin, fout
    integer :: iun, oun, ik, ib, iw, iwg, num_inc, cnt, ngx, ngy, ngz, nkk, nbnd, ip, ngpts, ierr
    real(kind=dp) :: rr, ci

    do ik = 1, num_kpts
      write (fin, '(a,i5.5,a,i1)') 'UNK', ik, '.', spin
      write (fout, '(a,a,i5.5,a,i1)') trim(outdir)//'/', 'UNK', ik, '.', spin
      if (formatted) then
        open (newunit=iun, file=fin, form='formatted', status='old', iostat=ierr)
      else
        open (newunit=iun, file=fin, form='unformatted', status='old', iostat=ierr)
      end if
      if (ierr /= 0) then
        call set_error_file(error, 'mrwf_write_unk: cannot open '//trim(fin), comm)
        return
      end if
      if (formatted) then
        read (iun, *) ngx, ngy, ngz, nkk, nbnd
      else
        read (iun) ngx, ngy, ngz, nkk, nbnd
      end if
      ngpts = ngx*ngy*ngz

      allocate (inc(nbnd))
      if (have_dis) then
        inc = dis_manifold%lwindow(1:nbnd, ik)
        num_inc = dis_manifold%ndimwin(ik)
      else
        inc = .true.
        num_inc = num_bands
      end if
      allocate (wtmp(ngpts, num_inc), rwv(ngpts, num_wann), cwv(ngpts, ng), buf(ngpts))

      cnt = 0
      do ib = 1, nbnd
        if (formatted) then
          do ip = 1, ngpts
            read (iun, *) rr, ci
            buf(ip) = cmplx(rr, ci, dp)
          end do
        else
          read (iun) (buf(ip), ip=1, ngpts)
        end if
        if (inc(ib)) then
          cnt = cnt + 1
          wtmp(:, cnt) = buf
        end if
      end do
      close (iun)

      ! stage 1: collapse in-window bands to num_wann states via u_opt
      rwv = cmplx_0
      do iw = 1, num_wann
        do ib = 1, num_inc
          rwv(:, iw) = rwv(:, iw) + u_opt(ib, iw, ik)*wtmp(:, ib)
        end do
      end do
      ! stage 2: apply the post-u_opt manifold gauge
      cwv = cmplx_0
      do iwg = 1, ng
        do iw = 1, num_wann
          cwv(:, iwg) = cwv(:, iwg) + bg(iw, iwg, ik)*rwv(:, iw)
        end do
      end do

      if (formatted) then
        open (newunit=oun, file=fout, form='formatted', action='write')
        write (oun, '(5i8)') ngx, ngy, ngz, ik, ng
        do iwg = 1, ng
          do ip = 1, ngpts
            write (oun, '(2f20.12)') real(cwv(ip, iwg), dp), aimag(cwv(ip, iwg))
          end do
        end do
      else
        open (newunit=oun, file=fout, form='unformatted', action='write')
        write (oun) ngx, ngy, ngz, ik, ng
        do iwg = 1, ng
          write (oun) (cwv(ip, iwg), ip=1, ngpts)
        end do
      end if
      close (oun)

      deallocate (inc, wtmp, rwv, cwv, buf)
    end do
  end subroutine mrwf_write_unk_files

end module w90_mrwf_mod
