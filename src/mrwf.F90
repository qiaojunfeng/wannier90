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

  use w90_constants, only: dp, cmplx_0, cmplx_1
  use w90_types, only: kmesh_info_type, atom_data_type
  use w90_error_base, only: w90_error_type
  use w90_error, only: set_error_input, set_error_fatal, set_error_alloc
  use w90_comms, only: w90_comm_type, mpirank
  use w90_utility, only: utility_zgemm_new, utility_diagonalize
  use w90_parallel_transport_mod, only: w90_parallel_transport

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
                      log_interp, seedname, stdout, error, comm)
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
    logical, intent(in) :: log_interp
    character(len=*), intent(in) :: seedname
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    complex(kind=dp), allocatable :: utot(:, :, :), hw(:, :), vmat(:, :, :), tmp(:, :)
    complex(kind=dp), allocatable :: mg(:, :, :, :), upt(:, :, :), vg(:, :), vg2(:, :)
    complex(kind=dp), allocatable :: gsplit(:, :, :), tt(:, :)
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

end module w90_mrwf_mod
