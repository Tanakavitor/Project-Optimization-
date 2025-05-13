# Portfolio Optimization

## Visão Geral  
Este projeto implementa uma otimização de portfólio em Fortran , utilizando paralelização via OpenMP para acelerar o processo



## Pré-requisitos  
- GNU Fortran (gfortran) ≥ 6.0  
- Suporte a OpenMP (já incluído no gfortran)  
- Ambiente Unix-like

## Compilação e Execução  

```bash
# Compila com otimizações e OpenMP
gfortran -O3 -fopenmp main.f90 -o portfolio_opt

# Remove limite de stack
ulimit -s unlimited

# Define número de threads (ex.: 4)
export OMP_NUM_THREADS=4

# Executa o programa
./portfolio_opt
```



## Parâmetros e Opções  
- `OMP_NUM_THREADS` — define quantos núcleos usar  
- `main.f90` — código-fonte principal  
- `precos_dowjones_2024H2.csv` — arquivo de preços de entrada (gerado via Python)

...existing content...

## Download dos Dados  
Para coletar os preços históricos de cada ativo, podemos utilizar o script Python `get_data.py`.  

```bash
# Executar na raiz do projeto
python get_data.py
```

## Estrutura do Repositório  
```
Project-Optimization-/
├── main.f90                     # Código Fortran
├── first_semestersim.f90        # Simulação da carteira no 1º trimestre de 2025
├── get_data.py                  # Script Python para baixar dados
├── precos_dowjones_2024H2.csv   # Dados de preços (2º semestre de 2024)
├── precos_dowjones_2025Q1.csv   # Dados de preços (1º trimestre de 2025)
└── readme.MD                    # Documentação do projeto          
```

## Melhor carteira encontrada
**Sharpe Ratio**: `0.16134131891971898`

| Ação  | Peso (%) |
|-------|----------|
| AAPL  | 2.79%    |
| AMZN  | 0.81%    |
| AXP   | 0.32%    |
| BA    | 0.48%    |
| CAT   | 3.86%    |
| CRM   | 3.52%    |
| CSCO  | 2.05%    |
| CVX   | 1.40%    |
| DIS   | 0.24%    |
| GS    | 1.86%    |
| HD    | 3.04%    |
| HON   | 8.88%    |
| IBM   | 15.93%   |
| JNJ   | 8.43%    |
| JPM   | 2.25%    |
| KO    | 2.59%    |
| MCD   | 2.33%    |
| MMM   | 1.48%    |
| MRK   | 0.85%    |
| MSFT  | 3.28%    |
| NKE   | 4.30%    |
| NVDA  | 2.98%    |
| PG    | 7.58%    |
| TRV   | 13.92%   |
| WMT   | 4.84%    |


## Simulação — Primeiro Trimestre de 2025
```
 Retorno anualizado da carteira:    6.5415774472862565E-002
 Volatilidade anualizada:          0.13840742477890805
 Índice de Sharpe:                 0.47263197460221290
```


## Ferramentas de IA  
Foram utilizadas ferramentas de inteligência artificial no desenvolvimento, correção e otimização deste projeto.  