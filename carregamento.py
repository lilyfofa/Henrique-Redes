from funcao import FluxoDePotencia
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

for i in alpha:
    print(f"{i:.2f}")
    resultado = FluxoDePotencia(6, 100, 0.0001, [0, 1.4, 0.6, -0.9, -1, -0.9],
                        [0, 0, 0, -0.6, -0.7, -0.5], [1.05, 1.06, 1.05, 0, 0, 0],
                        [0, -1, -1, -1, -1, -1], dados_linha, i)
    Tensoes = resultado[0]
    soma = 0
    for i in range(0, len(Tensoes)):
        soma += Tensoes[i]
    media = soma/len(Tensoes)
    y.append(media)

plt.figure(1)
plt.plot(alpha, y)
plt.title("Tensão média do sistema em função do carregamento")
plt.xlabel("Carregamento")
plt.ylabel("Tensão média [pu]")
plt.xlim(1,  3.1)
plt.xticks(arange(1, 3.2, 0.1))
plt.grid()
plt.show()
