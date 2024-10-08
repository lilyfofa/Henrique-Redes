import numpy as np
from funcoes import FluxoDeCarga
from time import time


def pso(n_particulas, n_iteracoes, n_pgs):
    # Parâmetros iniciais
    posicoes = np.random.uniform(Vmax, Pgmax, (n_particulas, n_pgs))
    velocidades = np.zeros((n_particulas,n_pgs))
    # Fitness e melhores posicoes
    fitness = np.array([FluxoDeCarga(posicao, 1)[0] for posicao in posicoes])
    pbest = np.copy(posicoes)
    valor_pbest = np.copy(fitness)
    posicao_gbest = np.zeros(n_pgs)
    valor_gbest = np.inf

    # Inicialização de pbest e gbest
    for i in range(n_particulas):
        if fitness[i] < valor_gbest:
            valor_gbest = fitness[i]
            posicao_gbest = np.copy(posicoes[i])  # Usar np.copy para evitar referências

    print(f'Iteração 0 concluída! Perda mínima: {valor_gbest:.4f} pu')

    for i in range(n_iteracoes):
        # Definição de velocidade
        r1 = np.random.uniform(0, 1, (n_particulas, n_pgs))
        r2 = np.random.uniform(0, 1, (n_particulas, n_pgs))
        velocidades = w * velocidades + c1 * r1 * (pbest - posicoes) + c2 * r2 * (posicao_gbest - posicoes)

        # Atualização da posição
        posicoes += velocidades

        # Aplicação dos limites de tensão e potência
        for j in range(n_particulas):
            posicoes[j, :] = np.clip(posicoes[j, :], Pgmin, Pgmax)

        # Atualização do fitness
        fitness = np.array([FluxoDeCarga(posicao, 1)[0] for posicao in posicoes])

        # Atualização dos melhores valores
        for j in range(n_particulas):
            if fitness[j] < valor_pbest[j]:
                valor_pbest[j] = fitness[j]
                pbest[j] = np.copy(posicoes[j])

        # Atualização de gbest
        melhor_fitness_idx = np.argmin(fitness)
        if fitness[melhor_fitness_idx] < valor_gbest:
            valor_gbest = fitness[melhor_fitness_idx]
            posicao_gbest = np.copy(posicoes[melhor_fitness_idx])
        print(f'Iteração {i + 1} concluída! Perda mínima: {valor_gbest:.4f} pu')
    return valor_gbest, posicao_gbest


# Parâmetros do PSO
n_particulas = 20
n_iteracoes = 200
n_barras = 6
n_pgs = 2

# Coeficientes do algoritmo
w = 0.9
c1 = 1.25
c2 = 1.25

# Limites das variáveis de estado
Vmin = 0.94
Vmax = 1.06
Pgmin = 0
Pgmax = 1.5

# Execução do PSO
print('PSO para minimização das perdas do sistema')
print('-' * 50)

t0 = time()
perda_minima, valor_vpg = pso(n_particulas, n_iteracoes, n_pgs)
t1 = time()

# Resultados
print('-' * 50)
print(f'Resultado final:')
variaveis = ['Pg1', 'Pg2']
for i in range(len(variaveis)):
    print(f'{variaveis[i]} = {valor_vpg[i]:.4f} pu')
print(f'Perdas totais = {perda_minima:.4f} pu')
print(f'Tempo de execução: {t1 - t0:.2f} s')
print('-' * 50)
print('(c) 2024 - EmpelTec Jr - Todos os direitos reservados')
