!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!------------------------------------------------------------!

!===============================================================!
!                                                               !
!  w90_parallel_transport                                       !
!                                                               !
!  Construct a smooth Bloch gauge by parallel transport across  !
!  the Monkhorst-Pack k-mesh, following the homotopy algorithm  !
!  of Gontier, Levitt and Siraj-dine (J. Math. Phys. 60, 2019). !
!                                                               !
!  This is a serial port of the Wannier.jl implementation; it   !
!  operates on full k-mesh arrays on every rank (the arrays are !
!  small: num_wann x num_wann per k-point).                     !
!                                                               !
!===============================================================!

module w90_parallel_transport_mod

  use w90_constants, only: dp, cmplx_0, cmplx_1, cmplx_i, pi, twopi
  use w90_types, only: kmesh_info_type
  use w90_error_base, only: w90_error_type
  use w90_error, only: set_error_input, set_error_fatal, set_error_alloc
  use w90_comms, only: w90_comm_type, mpirank
  use w90_utility, only: utility_zgemm_new

  implicit none

  private

  public :: w90_parallel_transport

contains

  !================================================================!
  !                      Linear-algebra helpers                    !
  !================================================================!

  !> Dummy SELECT for ZGEES (not referenced when sort = 'N').
  logical function pt_select_dummy(w)
    complex(kind=dp), intent(in) :: w
    pt_select_dummy = .false.
    if (real(w, dp) > huge(1.0_dp)) pt_select_dummy = .true. ! silence unused warning
  end function pt_select_dummy

  !> Loewdin orthonormalisation: nearest (semi-)unitary to `a`, O = U V^H
  !> from the SVD  a = U S V^H.
  function pt_lowdin(a) result(o)
    complex(kind=dp), intent(in) :: a(:, :)
    complex(kind=dp), allocatable :: o(:, :)
    complex(kind=dp), allocatable :: acopy(:, :), umat(:, :), vt(:, :), work(:)
    real(kind=dp), allocatable :: svals(:), rwork(:)
    integer :: n, info, lwork

    n = size(a, 1)
    allocate (o(n, n), acopy(n, n), umat(n, n), vt(n, n), svals(n), rwork(5*n))
    lwork = 4*n
    allocate (work(lwork))
    acopy = a
    call zgesvd('A', 'A', n, n, acopy, n, svals, umat, n, vt, n, work, lwork, rwork, info)
    call utility_zgemm_new(umat, vt, o, 'N', 'N') ! O = U . V^H
    deallocate (acopy, umat, vt, svals, rwork, work)
  end function pt_lowdin

  !> Eigendecomposition-based logarithm of a (near-)unitary matrix `o`.
  !> Returns eigenvectors `v` (unitary) and `logd` such that
  !> log(o) = v . diag(logd) . v^H. The branch of each eigenphase is the
  !> principal value, but phases within 0.01 of -pi are lifted by +2pi to keep
  !> a cluster of eigenvalues near -1 on the same side of the cut.
  subroutine pt_eig_log(o, v, logd)
    complex(kind=dp), intent(in) :: o(:, :)
    complex(kind=dp), intent(out) :: v(:, :)
    complex(kind=dp), intent(out) :: logd(:)
    complex(kind=dp), allocatable :: acopy(:, :), wvals(:), work(:)
    real(kind=dp), allocatable :: rwork(:)
    logical, allocatable :: bwork(:)
    real(kind=dp) :: theta
    integer :: n, info, lwork, sdim, i

    n = size(o, 1)
    lwork = 10*n
    allocate (acopy(n, n), wvals(n), work(lwork), rwork(n), bwork(n))
    acopy = o
    call zgees('V', 'N', pt_select_dummy, n, acopy, n, sdim, wvals, v, n, work, lwork, &
               rwork, bwork, info)
    do i = 1, n
      theta = atan2(aimag(wvals(i)), real(wvals(i), dp))
      if (theta < -pi + 0.01_dp) theta = theta + twopi
      logd(i) = cmplx_i*theta
    end do
    deallocate (acopy, wvals, work, rwork, bwork)
  end subroutine pt_eig_log

  !> Right-multiply factor exp(t.log O) = V . diag(exp(t.logd)) . V^H.
  function pt_pullback_matrix(v, logd, t) result(res)
    complex(kind=dp), intent(in) :: v(:, :), logd(:)
    real(kind=dp), intent(in) :: t
    complex(kind=dp), allocatable :: res(:, :), tmp(:, :)
    integer :: n, i

    n = size(v, 1)
    allocate (res(n, n), tmp(n, n))
    do i = 1, n
      tmp(:, i) = v(:, i)*exp(t*logd(i)) ! V . diag(exp(t logd))
    end do
    call utility_zgemm_new(tmp, v, res, 'N', 'C') ! . V^H
    deallocate (tmp)
  end function pt_pullback_matrix

  !> Matrix power of a (near-)unitary matrix: O**p = V . diag(d**p) . V^H.
  function pt_powm(o, p) result(res)
    complex(kind=dp), intent(in) :: o(:, :)
    real(kind=dp), intent(in) :: p
    complex(kind=dp), allocatable :: res(:, :), v(:, :), acopy(:, :), wvals(:), work(:), tmp(:, :)
    real(kind=dp), allocatable :: rwork(:)
    logical, allocatable :: bwork(:)
    real(kind=dp) :: theta
    integer :: n, info, lwork, sdim, i

    n = size(o, 1)
    lwork = 10*n
    allocate (res(n, n), v(n, n), acopy(n, n), wvals(n), work(lwork), rwork(n), bwork(n), tmp(n, n))
    acopy = o
    call zgees('V', 'N', pt_select_dummy, n, acopy, n, sdim, wvals, v, n, work, lwork, &
               rwork, bwork, info)
    do i = 1, n
      theta = atan2(aimag(wvals(i)), real(wvals(i), dp)) ! d = exp(i theta) => d**p = exp(i p theta)
      tmp(:, i) = v(:, i)*exp(cmplx_i*p*theta)
    end do
    call utility_zgemm_new(tmp, v, res, 'N', 'C')
    deallocate (v, acopy, wvals, work, rwork, bwork, tmp)
  end function pt_powm

  !> Determinant of a complex matrix via LU factorisation.
  function pt_zdet(a) result(det)
    complex(kind=dp), intent(in) :: a(:, :)
    complex(kind=dp) :: det
    complex(kind=dp), allocatable :: acopy(:, :)
    integer, allocatable :: ipiv(:)
    integer :: n, info, i

    n = size(a, 1)
    allocate (acopy(n, n), ipiv(n))
    acopy = a
    call zgetrf(n, n, acopy, n, ipiv, info)
    det = cmplx_1
    do i = 1, n
      det = det*acopy(i, i)
      if (ipiv(i) /= i) det = -det
    end do
    deallocate (acopy, ipiv)
  end function pt_zdet

  !================================================================!
  !                   Grid / neighbour helpers                     !
  !================================================================!

  !> Build the map from grid position (i,j,k) (1-based) to global k-index, and
  !> the +/- cardinal neighbour tables nn_plus/nn_minus(ik, dir), dir=1,2,3.
  subroutine pt_build_maps(kmesh_info, kpt_latt, mp_grid, num_kpts, xyz_k, nn_plus, nn_minus, &
                           error, comm)
    type(kmesh_info_type), intent(in) :: kmesh_info
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    integer, intent(in) :: mp_grid(3), num_kpts
    integer, allocatable, intent(out) :: xyz_k(:, :, :)
    integer, allocatable, intent(out) :: nn_plus(:, :), nn_minus(:, :)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer :: ik, nn, dir, gi, gj, gk, ik2
    real(kind=dp) :: disp(3), target_plus(3)
    real(kind=dp), parameter :: tol = 1.0e-4_dp

    allocate (xyz_k(mp_grid(1), mp_grid(2), mp_grid(3)))
    allocate (nn_plus(num_kpts, 3), nn_minus(num_kpts, 3))
    xyz_k = 0
    nn_plus = 0
    nn_minus = 0

    do ik = 1, num_kpts
      gi = modulo(nint(kpt_latt(1, ik)*mp_grid(1)), mp_grid(1)) + 1
      gj = modulo(nint(kpt_latt(2, ik)*mp_grid(2)), mp_grid(2)) + 1
      gk = modulo(nint(kpt_latt(3, ik)*mp_grid(3)), mp_grid(3)) + 1
      xyz_k(gi, gj, gk) = ik
    end do

    do ik = 1, num_kpts
      do nn = 1, kmesh_info%nntot
        ik2 = kmesh_info%nnlist(ik, nn)
        disp(:) = kpt_latt(:, ik2) + real(kmesh_info%nncell(:, ik, nn), dp) - kpt_latt(:, ik)
        do dir = 1, 3
          target_plus = 0.0_dp
          target_plus(dir) = 1.0_dp/real(mp_grid(dir), dp)
          if (all(abs(disp - target_plus) < tol)) nn_plus(ik, dir) = nn
          if (all(abs(disp + target_plus) < tol)) nn_minus(ik, dir) = nn
        end do
      end do
    end do

    if (any(nn_plus == 0) .or. any(nn_minus == 0)) then
      call set_error_input(error, &
                           'parallel_transport: the k-mesh does not contain the 6 cardinal neighbours &
                           &(+/-kx,+/-ky,+/-kz); a full Monkhorst-Pack finite-difference stencil is required', &
                           comm)
      return
    end if
  end subroutine pt_build_maps

  !================================================================!
  !                    Obstruction / propagation                   !
  !================================================================!

  !> Loewdin-orthonormalised obstruction between k1 and its +dir neighbour k2,
  !> N = U(k1)^H . M_raw(:,:,nn+,k1) . U(k2).
  function pt_obstruction(u, m, nn_plus, k1, dir, num_wann, nnlist) result(o)
    complex(kind=dp), intent(in) :: u(:, :, :), m(:, :, :, :)
    integer, intent(in) :: nn_plus(:, :), k1, dir, num_wann, nnlist(:, :)
    complex(kind=dp), allocatable :: o(:, :), tmp(:, :), n(:, :)
    integer :: nn, k2

    nn = nn_plus(k1, dir)
    k2 = nnlist(k1, nn)
    allocate (tmp(num_wann, num_wann), n(num_wann, num_wann))
    ! tmp = U(k1)^H . M(:,:,nn,k1)
    call utility_zgemm_new(u(:, :, k1), m(:, :, nn, k1), tmp, 'C', 'N')
    ! n = tmp . U(k2)
    call utility_zgemm_new(tmp, u(:, :, k2), n, 'N', 'N')
    o = pt_lowdin(n)
    deallocate (tmp, n)
  end function pt_obstruction

  !> Propagate the gauge along a line of consecutive kpoints `kpts` (in +dir),
  !> starting from the first (already fixed). Uses the -dir overlaps.
  subroutine pt_propagate(u, m, nn_minus, kpts, dir, num_wann)
    complex(kind=dp), intent(inout) :: u(:, :, :)
    complex(kind=dp), intent(in) :: m(:, :, :, :)
    integer, intent(in) :: nn_minus(:, :), kpts(:), dir, num_wann
    complex(kind=dp), allocatable :: tmp(:, :)
    integer :: i, ik, ik0, nn

    allocate (tmp(num_wann, num_wann))
    do i = 2, size(kpts)
      ik = kpts(i)
      ik0 = kpts(i - 1)
      nn = nn_minus(ik, dir) ! neighbour of ik in -dir is ik0
      call utility_zgemm_new(m(:, :, nn, ik), u(:, :, ik0), tmp, 'N', 'N')
      u(:, :, ik) = pt_lowdin(tmp)
    end do
    deallocate (tmp)
  end subroutine pt_propagate

  !> Factor the U(1) determinant winding out of an edge obstruction path so the
  !> SU(N) contraction only handles the special-unitary part. The winding is
  !> spread equally over all num_wann bands. Returns the normalised path and the
  !> continuous (unwrapped) determinant phase logD.
  subroutine pt_factor_det_winding(o_path, num_wann, logd)
    complex(kind=dp), intent(inout) :: o_path(:, :, :)
    integer, intent(in) :: num_wann
    real(kind=dp), intent(out) :: logd(:)
    integer :: nk, i, kbest, k
    real(kind=dp) :: best, trial

    nk = size(o_path, 3)
    do i = 1, nk
      logd(i) = aimag(log(pt_zdet(o_path(:, :, i))))
    end do
    do i = 2, nk
      kbest = 0
      best = huge(1.0_dp)
      do k = -1, 1
        trial = abs(logd(i) + twopi*k - logd(i - 1))
        if (trial < best) then
          best = trial
          kbest = k
        end if
      end do
      logd(i) = logd(i) + twopi*kbest
    end do
    do i = 1, nk
      o_path(:, :, i) = exp(-cmplx_i*logd(i)/real(num_wann, dp))*o_path(:, :, i)
    end do
  end subroutine pt_factor_det_winding

  !> Scalar phase reattaching the determinant winding, spread over num_wann.
  complex(kind=dp) function pt_det_winding_phase(logd_i, t, num_wann)
    real(kind=dp), intent(in) :: logd_i, t
    integer, intent(in) :: num_wann
    pt_det_winding_phase = exp(cmplx_i*logd_i*t/real(num_wann, dp))
  end function pt_det_winding_phase

  !================================================================!
  !            SU(N) contraction (matrix_transport)                !
  !================================================================!

  !> Choose a pole far from the column path (GLS2019). Appends the next column
  !> index to `columns` and returns the pole.
  subroutine pt_choose_pole(matrix_path, columns, ncol_used, prev_poles, pole)
    complex(kind=dp), intent(in) :: matrix_path(:, :, :)
    integer, intent(inout), allocatable :: columns(:)
    integer, intent(inout) :: ncol_used
    complex(kind=dp), intent(in) :: prev_poles(:, :)
    complex(kind=dp), intent(out) :: pole(:)
    integer, parameter :: n_iter = 100
    integer :: n_col, n_k, col, i, j, m
    real(kind=dp) :: diam, d, max_dist, dist
    complex(kind=dp), allocatable :: vec_path(:, :), max_pole(:), proj(:, :), tmpp(:), bary(:)
    real(kind=dp), allocatable :: rr(:), ri(:)
    integer, allocatable :: newcols(:)

    n_col = size(matrix_path, 2)
    n_k = size(matrix_path, 3)

    if (ncol_used == 0) then
      allocate (columns(1)); columns(1) = 1; ncol_used = 1
    else
      allocate (newcols(ncol_used + 1))
      newcols(1:ncol_used) = columns(1:ncol_used)
      newcols(ncol_used + 1) = ncol_used + 1
      call move_alloc(newcols, columns)
      ncol_used = ncol_used + 1
    end if
    col = columns(ncol_used)

    allocate (vec_path(n_col, n_k), max_pole(n_col), proj(n_col, n_col), tmpp(n_col), bary(n_col))
    allocate (rr(n_col), ri(n_col))
    vec_path(:, :) = matrix_path(:, col, :)

    diam = 0.0_dp
    do i = 1, n_k
      do j = 1, n_k
        d = sqrt(real(dot_product(vec_path(:, i) - vec_path(:, j), vec_path(:, i) - vec_path(:, j)), dp))
        if (d > diam) diam = d
      end do
    end do

    pole = cmplx_0
    pole(col) = cmplx_1

    if (diam > 1.5_dp) then
      max_pole = pole
      max_dist = 0.0_dp
      m = 1
      do while (m < n_iter)
        dist = huge(1.0_dp)
        do i = 1, n_k
          d = sqrt(real(dot_product(pole + vec_path(:, i), pole + vec_path(:, i)), dp))
          if (d < dist) dist = d
        end do
        if (dist > max_dist) then
          max_pole = pole
          max_dist = dist
        end if
        if (m <= n_col) then
          pole = cmplx_0; pole(m) = cmplx_1
        else if (m <= 2*n_col) then
          pole = cmplx_0; pole(m - n_col) = -cmplx_1
        else
          call random_number(rr); call random_number(ri)
          pole = cmplx(2.0_dp*rr - 1.0_dp, 2.0_dp*ri - 1.0_dp, dp)
        end if
        ! P = I - prev_poles prev_poles^H ; pole = P.pole
        call utility_zgemm_new(prev_poles, prev_poles, proj, 'N', 'C')
        tmpp = pole
        pole = tmpp - matmul(proj, tmpp)
        pole = pole/sqrt(real(dot_product(pole, pole), dp))
        m = m + 1
      end do
      pole = max_pole
    else
      bary = cmplx_0
      do i = 1, n_k
        bary = bary + vec_path(:, i)
      end do
      pole = bary/sqrt(real(dot_product(bary, bary), dp))
    end if

    deallocate (vec_path, max_pole, proj, tmpp, bary, rr, ri)
  end subroutine pt_choose_pole

  !> Parallel-transport the columns listed in `columns` along the frame path
  !> (backwards, from t = n_t to t = 1). Returns U(:, :, ik, it).
  subroutine pt_matrix_parallel_transport(frame_path, matrix_path, columns, ncol, u)
    complex(kind=dp), intent(in) :: frame_path(:, :, :, :), matrix_path(:, :, :)
    integer, intent(in) :: columns(:), ncol
    complex(kind=dp), intent(out) :: u(:, :, :, :)
    integer :: n_k, n_t, n_col, it, ik, ic, iit, idx, jj, ncnt
    complex(kind=dp), allocatable :: proj(:, :), src(:)
    integer, allocatable :: not_columns(:)
    real(kind=dp) :: nrm
    logical :: in_cols

    n_col = size(matrix_path, 2)
    n_k = size(frame_path, 3)
    n_t = size(frame_path, 4)
    allocate (proj(n_col, n_col), src(n_col))

    ! list of columns NOT being transported
    ncnt = 0
    allocate (not_columns(n_col))
    do jj = 1, n_col
      in_cols = .false.
      do idx = 1, ncol
        if (columns(idx) == jj) in_cols = .true.
      end do
      if (.not. in_cols) then
        ncnt = ncnt + 1
        not_columns(ncnt) = jj
      end if
    end do

    u = cmplx_0 ! match the reference (fresh zeros each call)
    do it = 1, n_t
      do ik = 1, n_k
        iit = n_t - it + 1
        ! P = I - F F^H, with F the full frame matrix at (ik, iit): projects out
        ! the span of the already-fixed frame columns (matches the reference).
        call utility_zgemm_new(frame_path(:, :, ik, iit), frame_path(:, :, ik, iit), proj, 'N', 'C')
        proj = -proj
        do jj = 1, n_col
          proj(jj, jj) = proj(jj, jj) + cmplx_1
        end do
        do ic = 1, n_col
          in_cols = .false.
          idx = 0
          do jj = 1, ncol
            if (columns(jj) == ic) then
              in_cols = .true.
              idx = jj
            end if
          end do

          if (in_cols) then
            nrm = sqrt(real(dot_product(frame_path(:, idx, ik, iit), frame_path(:, idx, ik, iit)), dp))
            if (abs(nrm - 1.0_dp) < 1.0e-2_dp) then
              u(:, ic, ik, iit) = frame_path(:, idx, ik, iit)
            else
              if (it > 1) then
                src = u(:, ic, ik, iit + 1)
              else
                src = matrix_path(:, ic, ik)
              end if
              u(:, ic, ik, iit) = matmul(proj, src)
            end if
          else
            if (it > 1) then
              src = u(:, ic, ik, iit + 1)
            else
              src = matrix_path(:, ic, ik)
            end if
            u(:, ic, ik, iit) = matmul(proj, src)
          end if
          nrm = sqrt(real(dot_product(u(:, ic, ik, iit), u(:, ic, ik, iit)), dp))
          u(:, ic, ik, iit) = u(:, ic, ik, iit)/nrm

          ! Loewdin-orthonormalise the remaining columns at time slice `it`
          ! (load-bearing: matches the reference implementation exactly; inside
          ! the ic loop, as in the reference).
          if (ncnt > 0) then
            call pt_lowdin_cols(u, not_columns(1:ncnt), ik, it)
          end if
        end do
      end do
    end do

    deallocate (proj, src, not_columns)
  end subroutine pt_matrix_parallel_transport

  !> Loewdin-orthonormalise the sub-block u(:, cols, ik, it) in place.
  subroutine pt_lowdin_cols(u, cols, ik, it)
    complex(kind=dp), intent(inout) :: u(:, :, :, :)
    integer, intent(in) :: cols(:), ik, it
    complex(kind=dp), allocatable :: block(:, :), ob(:, :)
    integer :: n, nc, j

    n = size(u, 1)
    nc = size(cols)
    allocate (block(n, nc))
    do j = 1, nc
      block(:, j) = u(:, cols(j), ik, it)
    end do
    ob = pt_lowdin_rect(block)
    do j = 1, nc
      u(:, cols(j), ik, it) = ob(:, j)
    end do
    deallocate (block)
  end subroutine pt_lowdin_cols

  !> Loewdin orthonormalisation of a (possibly rectangular, n x nc, n >= nc)
  !> column block: O = U(:,1:nc) V^H from the thin SVD.
  function pt_lowdin_rect(a) result(o)
    complex(kind=dp), intent(in) :: a(:, :)
    complex(kind=dp), allocatable :: o(:, :)
    complex(kind=dp), allocatable :: acopy(:, :), umat(:, :), vt(:, :), work(:)
    real(kind=dp), allocatable :: svals(:), rwork(:)
    integer :: n, nc, info, lwork

    n = size(a, 1)
    nc = size(a, 2)
    allocate (o(n, nc), acopy(n, nc), umat(n, nc), vt(nc, nc), svals(nc), rwork(5*nc))
    lwork = 2*max(1, 2*nc + n)
    allocate (work(lwork))
    acopy = a
    call zgesvd('S', 'A', n, nc, acopy, n, svals, umat, n, vt, nc, work, lwork, rwork, info)
    call utility_zgemm_new(umat, vt, o, 'N', 'N')
    deallocate (acopy, umat, vt, svals, rwork, work)
  end function pt_lowdin_rect

  !> Contract the unitary matrix path (n_wann x n_wann x n_k) to a constant
  !> point using parallel transport, returning the interpolating frames
  !> U(:, :, ik, it) with it indexing the interpolation parameter t.
  subroutine pt_matrix_transport(matrix_path, t, u, stdout)
    complex(kind=dp), intent(in) :: matrix_path(:, :, :)
    real(kind=dp), intent(in) :: t(:)
    complex(kind=dp), intent(out) :: u(:, :, :, :)
    integer, intent(in) :: stdout
    integer :: n_col, n_k, n_t, col, ncol_used, ik, j, k, i, kbest, kk
    complex(kind=dp), allocatable :: prev_poles(:, :), frame_path(:, :, :, :), pole(:)
    complex(kind=dp), allocatable :: cc(:, :, :), e1(:), fp(:), ovec(:, :), oinv(:, :), w(:, :)
    integer, allocatable :: columns(:)
    real(kind=dp), allocatable :: phi(:)
    real(kind=dp) :: nrm, best, trial, mchern
    integer :: n_c

    n_col = size(matrix_path, 1)
    n_k = size(matrix_path, 3)
    n_t = size(t)

    if (n_t == 1) then
      u(:, :, :, 1) = matrix_path(:, :, :)
      return
    end if

    allocate (prev_poles(n_col, n_col), frame_path(n_col, n_col, n_k, n_t), pole(n_col))
    prev_poles = cmplx_0
    frame_path = cmplx_0
    ncol_used = 0

    do col = 1, n_col
      call pt_choose_pole(matrix_path, columns, ncol_used, prev_poles, pole)
      prev_poles(:, col) = pole

      call pt_matrix_parallel_transport(frame_path, matrix_path, columns, ncol_used, u)

      if (col < n_col) then
        n_c = n_col - col + 1
        allocate (cc(n_c, n_k, n_t), e1(n_c), fp(n_col))
        e1 = cmplx_0; e1(1) = cmplx_1
        do ik = 1, n_k
          do j = 1, n_c
            cc(j, ik, 1) = dot_product(u(:, col + j - 1, ik, n_t), pole) ! <U_j | pole>
          end do
          do j = 1, n_t
            cc(:, ik, j) = (1.0_dp - t(j))*cc(:, ik, 1)
            cc(1, ik, j) = cc(1, ik, j) + t(j)*e1(1)
            nrm = sqrt(real(dot_product(cc(:, ik, j), cc(:, ik, j)), dp))
            cc(:, ik, j) = cc(:, ik, j)/nrm
          end do
          do j = 1, n_t
            fp = cmplx_0
            do k = 1, n_c
              fp = fp + u(:, col + k - 1, ik, j)*cc(k, ik, j)
            end do
            nrm = sqrt(real(dot_product(fp, fp), dp))
            frame_path(:, col, ik, j) = fp/nrm
          end do
        end do
        deallocate (cc, e1, fp)
      else
        allocate (phi(n_k))
        do ik = 1, n_k
          phi(ik) = aimag(log(dot_product(pole, u(:, col, ik, 1)))) ! pole^H . U
          if (ik > 1) then
            kbest = 0; best = huge(1.0_dp)
            do kk = -1, 1
              trial = abs(phi(ik) + twopi*kk - phi(ik - 1))
              if (trial < best) then; best = trial; kbest = kk; end if
            end do
            phi(ik) = phi(ik) + twopi*kbest
          end if
        end do
        mchern = (phi(n_k) - phi(1))/twopi
        if (stdout > 0) write (stdout, '(3x,a,f12.6)') 'parallel_transport: Chern number = ', mchern
        do i = 1, n_k
          do j = 1, n_t
            u(:, col, i, j) = u(:, col, i, j)*exp(-cmplx_i*(1.0_dp - t(j))*phi(i))
          end do
        end do
        deallocate (phi)
      end if
    end do

    ! Bring the contraction point from Obs to I
    allocate (ovec(n_col, n_col), oinv(n_col, n_col), w(n_col, n_col))
    ovec = pt_lowdin(u(:, :, 1, 1))
    ! O^H raised to (1 - t)
    do i = 1, n_k
      do j = 1, n_t
        w = pt_powm(transpose(conjg(ovec)), 1.0_dp - t(j))
        u(:, :, i, j) = matmul(w, u(:, :, i, j))
      end do
    end do

    deallocate (prev_poles, frame_path, pole, ovec, oinv, w)
    if (allocated(columns)) deallocate (columns)
  end subroutine pt_matrix_transport

  !================================================================!
  !                          Error metric                          !
  !================================================================!

  !> Smoothness error: sum over k and the 3 cardinal directions of
  !> || lowdin(U(k1)^H M(:,:,nn+,k1) U(k2)) - I ||^2, sqrt-ed and normalised.
  function pt_compute_error(u, m, nn_plus, nnlist, num_wann, num_kpts) result(eps)
    complex(kind=dp), intent(in) :: u(:, :, :), m(:, :, :, :)
    integer, intent(in) :: nn_plus(:, :), nnlist(:, :), num_wann, num_kpts
    real(kind=dp) :: eps
    complex(kind=dp), allocatable :: o(:, :)
    integer :: ik, dir, i

    eps = 0.0_dp
    do ik = 1, num_kpts
      do dir = 1, 3
        o = pt_obstruction(u, m, nn_plus, ik, dir, num_wann, nnlist)
        do i = 1, num_wann
          o(i, i) = o(i, i) - cmplx_1
        end do
        eps = eps + real(sum(conjg(o)*o), dp)
      end do
    end do
    eps = sqrt(eps)/real(num_kpts, dp)
  end function pt_compute_error

  !================================================================!
  !                          Main driver                           !
  !================================================================!

  !> Parallel-transport the gauge from the first kpoint to all others.
  !>
  !> `m_matrix` are the raw overlaps <u_mk|u_n,k+b> (num_wann x num_wann x nntot
  !> x num_kpts) and `u_matrix` (num_wann x num_wann x num_kpts) is the current
  !> gauge, used as the seed when `use_gauge` is true (otherwise the identity).
  !> On exit `u_matrix` holds the parallel-transport gauge.
  subroutine w90_parallel_transport(kmesh_info, u_matrix, m_matrix, kpt_latt, mp_grid, &
                                    num_wann, num_kpts, use_gauge, log_interp, stdout, error, comm)
    type(kmesh_info_type), intent(in) :: kmesh_info
    complex(kind=dp), intent(inout) :: u_matrix(:, :, :)
    complex(kind=dp), intent(in) :: m_matrix(:, :, :, :)
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    integer, intent(in) :: mp_grid(3), num_wann, num_kpts, stdout
    logical, intent(in) :: use_gauge, log_interp
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, allocatable :: xyz_k(:, :, :), nn_plus(:, :), nn_minus(:, :), line(:)
    complex(kind=dp), allocatable :: u0(:, :, :), u(:, :, :)
    complex(kind=dp), allocatable :: vmat(:, :), logd(:), o(:, :)
    complex(kind=dp), allocatable :: oxy(:, :, :), uxy(:, :, :, :)
    complex(kind=dp), allocatable :: oedge(:, :, :), uedge(:, :, :, :)
    real(kind=dp), allocatable :: tx(:), ty(:), tz(:), logdw(:)
    integer :: nkx, nky, nkz, i, j, k, ik, dir, rank
    real(kind=dp) :: eps0, eps1

    rank = mpirank(comm)
    nkx = mp_grid(1); nky = mp_grid(2); nkz = mp_grid(3)

    call pt_build_maps(kmesh_info, kpt_latt, mp_grid, num_kpts, xyz_k, nn_plus, nn_minus, error, comm)
    if (allocated(error)) return

    allocate (tx(nkx), ty(nky), tz(nkz))
    do i = 1, nkx; tx(i) = real(i - 1, dp)/real(nkx, dp); end do
    do i = 1, nky; ty(i) = real(i - 1, dp)/real(nky, dp); end do
    do i = 1, nkz; tz(i) = real(i - 1, dp)/real(nkz, dp); end do

    ! initialise the RNG deterministically so all ranks agree
    call pt_seed_rng()

    ! keep the original gauge for the initial-error report
    allocate (u0(num_wann, num_wann, num_kpts), u(num_wann, num_wann, num_kpts))
    u0 = u_matrix
    if (use_gauge) then
      u = u_matrix
    else
      u = cmplx_0
      do ik = 1, num_kpts
        do i = 1, num_wann
          u(i, i, ik) = cmplx_1
        end do
      end do
    end if

    allocate (vmat(num_wann, num_wann), logd(num_wann), line(max(nkx, max(nky, nkz))))

    ! 1. propagate along kx
    do i = 1, nkx; line(i) = xyz_k(i, 1, 1); end do
    call pt_propagate(u, m_matrix, nn_minus, line(1:nkx), 1, num_wann)
    ! corner obstruction O1
    o = pt_obstruction(u, m_matrix, nn_plus, xyz_k(nkx, 1, 1), 1, num_wann, kmesh_info%nnlist)
    call pt_eig_log(o, vmat, logd)
    do i = 1, nkx
      ik = xyz_k(i, 1, 1)
      u(:, :, ik) = matmul(u(:, :, ik), pt_pullback_matrix(vmat, logd, tx(i)))
    end do

    ! 2. propagate along ky
    do i = 1, nkx
      do j = 1, nky; line(j) = xyz_k(i, j, 1); end do
      call pt_propagate(u, m_matrix, nn_minus, line(1:nky), 2, num_wann)
    end do
    ! corner O2
    o = pt_obstruction(u, m_matrix, nn_plus, xyz_k(1, nky, 1), 2, num_wann, kmesh_info%nnlist)
    call pt_eig_log(o, vmat, logd)
    do i = 1, nkx
      do j = 1, nky
        ik = xyz_k(i, j, 1)
        u(:, :, ik) = matmul(u(:, :, ik), pt_pullback_matrix(vmat, logd, ty(j)))
      end do
    end do

    ! line obstruction Oxy, at ky = end along kx
    allocate (oxy(num_wann, num_wann, nkx), logdw(nkx), uxy(num_wann, num_wann, nkx, nky))
    do i = 1, nkx
      oxy(:, :, i) = pt_obstruction(u, m_matrix, nn_plus, xyz_k(i, nky, 1), 2, num_wann, kmesh_info%nnlist)
    end do
    call pt_factor_det_winding(oxy, num_wann, logdw)
    if (log_interp) then
      do i = 1, nkx
        do j = 1, nky
          uxy(:, :, i, j) = pt_powm(oxy(:, :, i), ty(j))
        end do
      end do
    else
      call pt_matrix_transport(oxy, ty, uxy, merge(stdout, -1, rank == 0))
    end if
    do i = 1, nkx
      do j = 1, nky
        ik = xyz_k(i, j, 1)
        u(:, :, ik) = matmul(u(:, :, ik), pt_det_winding_phase(logdw(i), ty(j), num_wann)*uxy(:, :, i, j))
      end do
    end do
    deallocate (oxy, logdw, uxy)

    ! 3. propagate along kz
    do i = 1, nkx
      do j = 1, nky
        do k = 1, nkz; line(k) = xyz_k(i, j, k); end do
        call pt_propagate(u, m_matrix, nn_minus, line(1:nkz), 3, num_wann)
      end do
    end do
    ! corner O4
    o = pt_obstruction(u, m_matrix, nn_plus, xyz_k(1, 1, nkz), 3, num_wann, kmesh_info%nnlist)
    call pt_eig_log(o, vmat, logd)
    do k = 1, nkz
      do i = 1, nkx
        do j = 1, nky
          ik = xyz_k(i, j, k)
          u(:, :, ik) = matmul(u(:, :, ik), pt_pullback_matrix(vmat, logd, tz(k)))
        end do
      end do
    end do

    ! edge Oxz, at kz = end along kx
    allocate (oedge(num_wann, num_wann, nkx), logdw(nkx), uedge(num_wann, num_wann, nkx, nkz))
    do i = 1, nkx
      oedge(:, :, i) = pt_obstruction(u, m_matrix, nn_plus, xyz_k(i, 1, nkz), 3, num_wann, kmesh_info%nnlist)
    end do
    call pt_factor_det_winding(oedge, num_wann, logdw)
    if (log_interp) then
      do i = 1, nkx
        do k = 1, nkz
          uedge(:, :, i, k) = pt_powm(oedge(:, :, i), tz(k))
        end do
      end do
    else
      call pt_matrix_transport(oedge, tz, uedge, merge(stdout, -1, rank == 0))
    end if
    do i = 1, nkx
      do k = 1, nkz
        do j = 1, nky
          ik = xyz_k(i, j, k)
          u(:, :, ik) = matmul(u(:, :, ik), pt_det_winding_phase(logdw(i), tz(k), num_wann)*uedge(:, :, i, k))
        end do
      end do
    end do
    deallocate (oedge, logdw, uedge)

    ! edge Oyz, at kz = end along ky
    allocate (oedge(num_wann, num_wann, nky), logdw(nky), uedge(num_wann, num_wann, nky, nkz))
    do j = 1, nky
      oedge(:, :, j) = pt_obstruction(u, m_matrix, nn_plus, xyz_k(1, j, nkz), 3, num_wann, kmesh_info%nnlist)
    end do
    call pt_factor_det_winding(oedge, num_wann, logdw)
    if (log_interp) then
      do j = 1, nky
        do k = 1, nkz
          uedge(:, :, j, k) = pt_powm(oedge(:, :, j), tz(k))
        end do
      end do
    else
      call pt_matrix_transport(oedge, tz, uedge, merge(stdout, -1, rank == 0))
    end if
    do j = 1, nky
      do k = 1, nkz
        do i = 1, nkx
          ik = xyz_k(i, j, k)
          u(:, :, ik) = matmul(u(:, :, ik), pt_det_winding_phase(logdw(j), tz(k), num_wann)*uedge(:, :, j, k))
        end do
      end do
    end do
    deallocate (oedge, logdw, uedge)

    ! surface: full kz fix
    do i = 1, nkx
      do j = 1, nky
        o = pt_obstruction(u, m_matrix, nn_plus, xyz_k(i, j, nkz), 3, num_wann, kmesh_info%nnlist)
        do k = 1, nkz
          ik = xyz_k(i, j, k)
          u(:, :, ik) = matmul(u(:, :, ik), pt_powm(o, tz(k)))
        end do
      end do
    end do

    ! error report
    eps0 = pt_compute_error(u0, m_matrix, nn_plus, kmesh_info%nnlist, num_wann, num_kpts)
    eps1 = pt_compute_error(u, m_matrix, nn_plus, kmesh_info%nnlist, num_wann, num_kpts)
    if (rank == 0 .and. stdout > 0) then
      write (stdout, '(3x,a,f12.6)') 'parallel_transport: initial smoothness error = ', eps0
      write (stdout, '(3x,a,f12.6)') 'parallel_transport: final   smoothness error = ', eps1
    end if

    u_matrix = u

    deallocate (xyz_k, nn_plus, nn_minus, tx, ty, tz, u0, u, vmat, logd, line)
    if (allocated(o)) deallocate (o)
  end subroutine w90_parallel_transport

  !> Deterministic RNG seed (so all MPI ranks produce identical results).
  subroutine pt_seed_rng()
    integer :: nseed, i
    integer, allocatable :: seed(:)
    call random_seed(size=nseed)
    allocate (seed(nseed))
    do i = 1, nseed
      seed(i) = 1234567 + 37*i
    end do
    call random_seed(put=seed)
    deallocate (seed)
  end subroutine pt_seed_rng

end module w90_parallel_transport_mod
