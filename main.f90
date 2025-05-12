program portfolio_opt
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use omp_lib
  implicit none
  ! Declarando variáveis
  integer, parameter :: N = 30, K = 25, N_SIMS = 1000
  real(dp), parameter :: W_MAX = 0.2_dp
  integer, parameter :: MAXOBS = 2000, MAXCOMB = 142506
  character(len=256) :: fname
  integer :: unit, ierr, obs, n_obs
  character(len=2000) :: line
  character(len=32) :: dummy_date
  real(dp), dimension(MAXOBS, N) :: prices
  real(dp), dimension(MAXOBS-1, N) :: ret_daily
  real(dp), dimension(N) :: mu, mean_ret
  real(dp), dimension(N, N) :: covar
  integer, dimension(K) :: idx, bestCombo
  integer, dimension(K, MAXCOMB) :: combos
  real(dp), dimension(K) :: bestWeights
  real(dp) :: bestSR

  ! Variáveis paraleliacao 
  integer :: i, j, c, kk, pos, sims
  real(dp) :: sumw, u, sr, localSR, media, stddev, mean_daily
  integer, dimension(K) :: iidx
  real(dp), dimension(K) :: w
  real(dp), pointer :: ret_carteira(:)
  integer :: ncomb

  ! Leitura do arquivo CSV
  fname = 'precos_dowjones_2024H2.csv'
  unit = 10

  print *, '>>> Lendo arquivo CSV...'
  open(unit, file=fname, status='old', action='read', iostat=ierr)
  if (ierr /= 0) stop 'Erro abrindo arquivo CSV'
  read(unit,'(A)',iostat=ierr) line
  if (ierr /= 0) stop 'Erro header errado'

  obs = 0
  do
    read(unit,*,iostat=ierr) dummy_date, (prices(obs+1,j), j=1,N)
    if (ierr /= 0) exit
    obs = obs + 1
    if (obs >= MAXOBS) exit
  end do
  close(unit)
  n_obs = obs - 1
  if (n_obs < 1) stop 'Csv formato errado'

  ! Pegando retorno diario
  do j = 1, N
    do i = 1, n_obs
      ret_daily(i,j) = (prices(i+1,j) - prices(i,j)) / prices(i,j)
    end do
  end do
 ! media de retorno
  do j = 1, N
    mu(j) = sum(ret_daily(:,j)) / n_obs * 252
    mean_ret(j) = sum(ret_daily(:,j)) / n_obs
  end do
    !calculando a matriz de covariancia
  do i = 1, N
    do j = 1, N
      covar(i,j) = sum((ret_daily(:,i) - mean_ret(i)) * (ret_daily(:,j) - mean_ret(j))) / real(n_obs-1, dp) * 252
    end do
  end do

  ! Geração de combinações 
  print *, '>>> Gerando combinacoes '
  ncomb = 0
  do kk = 1, K
    idx(kk) = kk
  end do

  do
    if (ncomb >= MAXCOMB) exit
    combos(:,ncomb+1) = idx
    ncomb = ncomb + 1
    pos = K
    do while (pos > 0 .and. idx(pos) == N-K+pos)
      pos = pos - 1
    end do
    if (pos == 0) exit
    idx(pos) = idx(pos) + 1
    do j = pos+1, K
      idx(j) = idx(j-1) + 1
    end do
  end do

  print *, '>>> Combinacoes feitas: ', ncomb


  bestSR = -1.0e9_dp
  call random_seed()

!$omp parallel default(shared) private(c,i,j,w,sumw,u,sims,media,stddev,sr,iidx,ret_carteira,localSR,mean_daily)
!$omp do schedule(dynamic)
!loop da paralizacao simula e gera melhor carteira
  do c = 1, ncomb
    iidx = combos(:,c)
    localSR = -1.0e9_dp
    allocate(ret_carteira(n_obs))
    sims = 0

    do while (sims < N_SIMS)
      sumw = 0.0_dp
      do i = 1, K
        call random_number(u)
        w(i) = -log(u)
        sumw = sumw + w(i)
      end do
      w = w / sumw
      if (any(w > W_MAX)) cycle
      sims = sims + 1

      ret_carteira = 0.0_dp
      do i = 1, n_obs
        do j = 1, K
          ret_carteira(i) = ret_carteira(i) + w(j) * ret_daily(i, iidx(j))
        end do
      end do

      mean_daily = sum(ret_carteira) / n_obs
      media = mean_daily * 252.0_dp
      stddev = sqrt(sum((ret_carteira - mean_daily)**2) / real(n_obs-1)) * sqrt(252.0_dp)
      sr = media / stddev
      if (sr > localSR) localSR = sr
    end do

    deallocate(ret_carteira)

    !$omp critical
    if (localSR > bestSR) then
      bestSR = localSR
      bestWeights = w
      bestCombo = iidx
    end if
    !$omp end critical
  end do
!$omp end do
!$omp end parallel

  ! === Resultado final ===
  print *, '>>> Simulações concluídas'
  print *, 'Melhor Sharpe encontrado: ', bestSR
  print *, 'Carteira vencedora (25 ativos):'
  do i = 1, K
    write(*,'(A,I2,A,F8.5)') 'Ativo ', bestCombo(i), ': peso = ', bestWeights(i)
  end do

end program portfolio_opt
