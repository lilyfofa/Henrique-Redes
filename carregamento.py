from funcoes import FluxoDePotencia
import matplotlib.pyplot as plt
from numpy import arange

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

step = 0.01
start = 1
stop = 3.0138
alpha = arange(start, stop, step)
y = []

print('Gráfico de carregamento')
print('-'*100)
print(f'Começo: {start}')
print(f'Fim: {stop}')
print(f'Passo: {step}')
print('-'*100)

for i in alpha:
    print(f"Calculando solução para lambda = {i:.2f}.")
    resultado = FluxoDePotencia(6, 100, 0.0001, [0, 1.4, 0.6, 0, 0, 0],
                        [0, 0, 0, 0.9, 1, 0.9], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0.6, 0.7, 0.5],
                                [1.05, 1.06, 1.05, 0, 0, 0], [0, -1, -1, -1, -1, -1], dados_linha, i)
    Tensoes = resultado[0]
    soma = 0
    for i in range(0, len(Tensoes)):
        soma += Tensoes[i]
    media = soma/len(Tensoes)
    y.append([tensao for tensao in Tensoes])

saida = [[], [], [], [], [], []]

for i in range(0, len(y)):
    for j in range(0, len(y[0])):
        saida[j].append(y[i][j])

print('-'*100)
print('Gerando gráfico...')

plt.figure(1)
plt.plot(alpha, saida[0], label='Barra 1')
plt.plot(alpha, saida[1], label='Barra 2')
plt.plot(alpha, saida[2], label='Barra 3')
plt.plot(alpha, saida[3], label='Barra 4')
plt.plot(alpha, saida[4], label='Barra 5')
plt.plot(alpha, saida[5], label='Barra 6')
plt.title("Tensão média do sistema em função do carregamento")
plt.xlabel("Carregamento")
plt.ylabel("Tensão média [pu]")
plt.xlim(1,  3.1)
plt.xticks(arange(1, 3.2, 0.1))
plt.grid()
plt.legend()
plt.show()
