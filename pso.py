import numpy as np
from funcao import FluxoDePotencia


def pso(n_particulas, n_iteracoes, n_tensoes, n_pgs):
    # Parâmetros iniciais
    tensoes = np.random.uniform(Vmin, Vmax, (n_particulas, n_tensoes))
    geracoes = np.random.uniform(Pgmin, Pgmax, (n_particulas, n_pgs))
    posicoes = np.block([tensoes, geracoes])
    velocidades = np.zeros((n_particulas, n_tensoes + n_pgs))
    # Fitness e melhores posicoes
    fitness = np.array([FluxoDePotencia(6, 100, 0.0001, np.concatenate(
        (np.zeros(1), posicao[n_tensoes:n_barras], np.zeros(n_barras - n_pgs - 1))).tolist(),
                                        [0, 0, 0, 0.9, 1, 0.9], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0.6, 0.7, 0.5],
                                        np.concatenate((posicao[0:n_tensoes], np.zeros(n_barras - n_tensoes))).tolist(),
                                        [0, -1, -1, -1, -1, -1],
                                        dados_linha, 1)[8] for posicao in posicoes])
    pbest = np.copy(posicoes)
    valor_gbest = np.min(fitness)
    posicao_gbest = pbest[np.argmin(fitness)]

    for i in range(n_iteracoes):
        # Definição de velocidade
        r1 = np.random.uniform(0, 1, (n_particulas, n_tensoes + n_pgs))
        r2 = np.random.uniform(0, 1, (n_particulas, n_tensoes + n_pgs))
        velocidades = w * velocidades + c1 * r1 * (pbest - posicoes) + c2 * r2 * (posicao_gbest - posicoes)
        # Atualização da posição
        posicoes += velocidades

        # Aplicação dos limites de tensão e potência
        for j in range(0, n_particulas):
            for k in range(0, n_tensoes):
                if posicoes[j][k] > Vmax:
                    posicoes[j][k] = Vmax
                elif posicoes[j][k] < Vmin:
                    posicoes[j][k] = Vmin
            for k in range(n_tensoes, n_tensoes + n_pgs):
                if posicoes[j][k] > Pgmax:
                    posicoes[j][k] = Pgmax
                elif posicoes[j][k] < Pgmin:
                    posicoes[j][k] = Pgmin

        # Atualização do fitness
        valor_funcao = np.array([FluxoDePotencia(6, 100, 0.0001, np.concatenate(
            (np.zeros(1), posicao[n_tensoes:n_barras], np.zeros(n_barras - n_pgs - 1))).tolist(),
                                                 [0, 0, 0, 0.9, 1, 0.9], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0.6, 0.7, 0.5],
                                                 np.concatenate((posicao[0:n_tensoes], np.zeros(n_barras - n_tensoes))).tolist(),
                                                 [0, -1, -1, -1, -1, -1],
                                                 dados_linha, 1)[8] for posicao in posicoes])

        # Atualização dos melhores valores
        improved_indices = np.where(valor_funcao < fitness)
        pbest[improved_indices] = posicoes[improved_indices]
        fitness[improved_indices] = valor_funcao[improved_indices]
        if np.min(valor_funcao) < valor_gbest:
            posicao_gbest = posicoes[np.argmin(valor_funcao)]
            valor_gbest = np.min(valor_funcao)
        print(f'Iteração {i+1} concluída! Perda mínima: {valor_gbest:.4f} pu')
    return valor_gbest, posicao_gbest


n_particulas = 5
n_iteracoes = 50
n_barras = 6
n_tensoes = 3
n_pgs = 2

w = 0.5
c1 = 1
c2 = 2

Vmin = 0.94
Vmax = 1.06

Pgmin = 0  # Limite mínimo?
Pgmax = 1.5  # Limite máximo?

dados_linha = [(0.1 + 1j * 0.2, 4, 1, 2),
               (0.05 + 1j * 0.2, 4, 1, 4),
               (0.08 + 1j * 0.3, 6, 1, 5),
               (0.05 + 1j * 0.25, 6, 2, 3),
               (0.05 + 1j * 0.1, 2, 2, 4),
               (0.1 + 1j * 0.3, 4, 2, 5),
               (0.07 + 1j * 0.2, 5, 2, 6),
               (0.12 + 1j * 0.26, 5, 3, 5),
               (0.02 + 1j * 0.10, 2, 3, 6),
               (0.20 + 1j * 0.40, 8, 4, 5),
               (0.10 + 1j * 0.30, 6, 5, 6)]

perda_minima, valor_vpg = pso(n_particulas=n_particulas, n_iteracoes=n_iteracoes, n_tensoes=n_tensoes, n_pgs=n_pgs)

print(valor_vpg)
