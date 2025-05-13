program simula_carteira
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: N = 25, MAXOBS = 2000
  real(dp), dimension(N, N) :: covar
  real(dp), dimension(N) :: pesos

  real(dp), dimension(MAXOBS, N) :: prices
  real(dp), dimension(MAXOBS-1, N) :: ret_ativos
  real(dp), allocatable :: ret_diario(:)
  integer :: i, j, obs, n_obs, ierr, unit
  character(len=256) :: line
  character(len=32) :: dummy_date
  character(len=256) :: fname
  real(dp) :: media, sigma_p, sharpe, mean_daily

  ! Pesos da carteira
  ! filepath: c:\Users\henri\Downloads\Project-Optimization-\first_semestersim.f90
   ! Pesos da carteira
    pesos = [ &
  0.02791_dp, 0.00805_dp, 0.00318_dp, 0.00484_dp, &
  0.03859_dp, 0.03515_dp, 0.02045_dp, 0.01396_dp, &
  0.00244_dp, 0.01858_dp, 0.03036_dp, 0.08884_dp, 0.15928_dp, &
  0.08434_dp, 0.02247_dp, 0.02588_dp, 0.02332_dp, &
  0.01480_dp, 0.00851_dp, 0.03282_dp, 0.04303_dp, &
  0.02980_dp, 0.07578_dp, 0.13917_dp, 0.04843_dp ]

  ! Lendo o CSV
  fname = 'precos_dowjones_2025Q1.csv'
  unit = 10
  open(unit, file=fname, status='old', action='read', iostat=ierr)
  if (ierr /= 0) stop 'Erro abrindo arquivo de preços'

  read(unit,'(A)',iostat=ierr) line
  if (ierr /= 0) stop 'Erro lendo header do CSV'

  obs = 0
  do
    read(unit,*,iostat=ierr) dummy_date, (prices(obs+1,j), j=1,N)
    if (ierr /= 0) exit
    obs = obs + 1
    if (obs >= MAXOBS) exit
  end do
  close(unit)

  n_obs = obs - 1
  if (n_obs < 1) stop 'Poucos dados'

  allocate(ret_diario(n_obs))

  ! Retornos individuais dos ativos
  do j = 1, N
    do i = 1, n_obs
      ret_ativos(i,j) = (prices(i+1,j) - prices(i,j)) / prices(i,j)
    end do
  end do

  ! Retorno diário da carteira
  do i = 1, n_obs
    ret_diario(i) = sum(ret_ativos(i,:) * pesos(:))
  end do

  ! Média diária
  mean_daily = sum(ret_diario) / n_obs
  media = mean_daily * 252.0_dp

  ! Matriz de covariância
  do i = 1, N
    do j = 1, N
      covar(i,j) = sum((ret_ativos(:,i) - sum(ret_ativos(:,i))/n_obs) * &
                       (ret_ativos(:,j) - sum(ret_ativos(:,j))/n_obs)) / real(n_obs - 1, dp)
    end do
  end do

  ! Volatilidade da carteira
  sigma_p = sqrt( dot_product(pesos, matmul(covar, pesos)) ) * sqrt(252.0_dp)


  ! Sharpe ratio (sem taxa livre de risco)
  sharpe = media / sigma_p

  print *, '>>> Simulação da carteira 2025 concluída'
  print *, 'Retorno anualizado da carteira: ', media
  print *, 'Volatilidade anualizada:        ', sigma_p
  print *, 'Índice de Sharpe:               ', sharpe

  deallocate(ret_diario)

end program simula_carteira
