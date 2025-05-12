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
├── main.f90               # Codigo Fortran
├── get_data.py            # Script Python para baixar dados
├── precos_dowjones_2024H2.csv  # Dados de preços
└── readme.MD              
```

## Melhor carteira encontrada
**Sharpe Ratio**: `3.58097`
| Ação | Peso (%) |
|------|----------|
| AAPL | 2.53% |
| AXP  | 1.16% |
| BA   | 9.88% |
| CAT  | 1.37% |
| CRM  | 10.49% |
| CSCO | 10.75% |
| CVX  | 10.95% |
| DIS  | 4.43% |
| GS   | 8.25% |
| HD   | 1.08% |
| HON  | 0.86% |
| IBM  | 1.65% |
| JNJ  | 2.67% |
| JPM  | 7.63% |
| KO   | 3.50% |
| MCD  | 0.47% |
| MMM  | 4.38% |
| MRK  | 1.11% |
| MSFT | 0.56% |
| NVDA | 1.58% |
| TRV  | 0.09% |
| UNH  | 2.93% |
| V    | 11.08% |
| VZ   | 0.33% |
| WMT  | 0.29% |

## Ferramentas de IA  
Foram utilizadas ferramentas de inteligência artificial no desenvolvimento, correção e otimização deste projeto.  