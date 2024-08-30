from sympy import *

# Corinthians
print("Corinthians")

# Definição das variáveis do sistema
n_barras = 6  # Número de barras
S_base = 100  # Potência base em MVA
erro = 0.0001  # Tolerância para o método de Newton

# Definição dos dados de barra
Pesp = [0, 1.4, 0.6, -0.9, -1, -0.9]  # Potência ativa líquida nas barras
Qesp = [0, 0, 0, -0.6, -0.7, -0.5]  # Potência ativa líquida nas barras
Tensao = [1.05, 1.06, 1.05, 0, 0, 0]  # Magnitude de tensão de cada barra (0 se desconhecido)
Fase = [0, -1, -1, -1, -1, -1]  # Fase da tensão de cada barra (-1 se desconhecido)

# Definição dos dados de linha
# Valor da impedância em pu (R + jX), valor em MVar dos elementos shunt, primeira barra e segunda barra
admitancias = [(0.1 + 1j * 0.2, 4, 1, 2),
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

# Montagem da matriz Ybus
ybus = zeros(n_barras, n_barras)
for item in admitancias:
    bsh = item[1] / (2 * S_base)
    y = 1 / item[0]
    barra1 = item[2] - 1
    barra2 = item[3] - 1
    ybus[barra1, barra1] += y + 1j * bsh
    ybus[barra2, barra2] += y + 1j * bsh
    ybus[barra1, barra2] -= y
    ybus[barra2, barra1] -= y

# Display da matriz Ybus
print("-" * 30)
print("Ybus")
for i in range(0, n_barras):
    print("-" * 30)
    print(f"Linha {i + 1}")
    for j in range(0, n_barras):
        if N(im(ybus[i, j])) >= 0:
            print(f"Y{i+1}{j+1} = {float(N(re(ybus[i, j]))):.4f} + j{float(N(im(ybus[i, j]))):.4f} [pu]")
        else:
            print(f"Y{i + 1}{j + 1} = {float(N(re(ybus[i, j]))):.4f} - j{float(abs(N(im(ybus[i, j])))):.4f} [pu]")
print("-" * 30)

# Cópia das variáveis para adequação ao código antigo
V = Tensao
Theta = Fase

# Matrizes Gbus e Bbus (Ybus = Gbus + jBbus)
gbus = re(ybus)
bbus = im(ybus)

# Listas para as equações, as variáveis, os chutes e os nomes
equacoes = []
variaveis = []
chute_variaveis = []
nomes_equacoes = []

# Acréscimo das variáveis Theta
for i in range(0, len(Tensao)):
    if Theta[i] == -1:
        Theta[i] = symbols(f"Theta{i + 1}")
        variaveis.append(symbols(f"Theta{i + 1}"))
        chute_variaveis.append((symbols(f"Theta{i + 1}"), 0))

# Acréscimo das variáveis V
for i in range(0, len(Tensao)):
    if Tensao[i] == 0:
        V[i] = symbols(f"V{i + 1}")
        variaveis.append(symbols(f"V{i + 1}"))
        chute_variaveis.append((symbols(f"V{i + 1}"), 1))

# Construção das equações deltaP
for i in range(0, n_barras):
    if Pesp[i] != 0:
        expressao = Pesp[i]
        for j in range(0, n_barras):
            expressao -= V[i] * V[j] * (gbus[i, j] * cos(Theta[i] - Theta[j]) + bbus[i, j] * sin(Theta[i] - Theta[j]))
        equacoes.append(expressao)
        nomes_equacoes.append(f"deltaP{i + 1}")

# Construção das equações deltaQ
for i in range(0, n_barras):
    if Qesp[i] != 0:
        expressao = Qesp[i]
        for j in range(0, n_barras):
            expressao -= V[i] * V[j] * (gbus[i, j] * sin(Theta[i] - Theta[j]) - bbus[i, j] * cos(Theta[i] - Theta[j]))
        equacoes.append(expressao)
        nomes_equacoes.append(f"deltaQ{i + 1}")

# Construção do jacobiano
jacobiano = zeros(len(equacoes), len(variaveis))
for i in range(0, len(equacoes)):
    for j in range(0, len(variaveis)):
        jacobiano[i, j] = diff(equacoes[i], variaveis[j])

# Início do método iterativo (flatstart)
print("Início do método iterativo")
print("-" * 30)
k = 0
valores_calculados = zeros(len(equacoes), 1)

# Cálculo das equações para o primeiro chute
for i in range(0, len(valores_calculados)):
    valores_calculados[i] = [equacoes[i].subs(chute_variaveis)]

# Display dos valores encontrados
print(f"Iteração {k}")
print("-" * 30)
print("Valores das variáveis")
for i in chute_variaveis:
    if str(i[0])[0] == 'V':
        print(f"{i[0]}: {i[1]:.4f} pu")
    else:
        print(f"{i[0]}: {i[1]:.4f} rad = {float(N(i[1] * 180 / pi)):.4f} [graus]")
print("Valores das equações")
for i in range(0, len(valores_calculados)):
    print(f"{nomes_equacoes[i]} = {valores_calculados[i]:.4f} pu")

# Laço do método iterativo
while max(abs(valores_calculados)) > erro:
    jac = jacobiano.subs(chute_variaveis)  # Cálculo do valor númerico do jacobiano
    delta = -1 * jac.inv() * valores_calculados  # Encontrando os ajustes
    novos_valores = []  # Matriz auxiliar pois uma tupla não pode ser alterada diretamente
    for i in range(0, len(chute_variaveis)):
        novos_valores.append((variaveis[i], chute_variaveis[i][1] + delta[i])) # Ajustes dos chutes das variáveis
    chute_variaveis = novos_valores # Substituição da tupla antiga
    for i in range(0, len(valores_calculados)):
        valores_calculados[i] = [equacoes[i].subs(chute_variaveis)]  # Novos valores para as equações
    # Display dos ajustes
    print("Valores dos ajustes")
    for i in range(0, len(chute_variaveis)):
        if str(variaveis[i])[0] == 'V':
            print(f"delta{variaveis[i]}: {delta[i]:.4f} pu")
        else:
            print(f"delta{variaveis[i]}: {delta[i]:.4f} rad")
    k += 1
    # Display dos novos valores encontrados para as variáveis e os valores das equações
    print("-" * 30)
    print(f"Iteração {k}")
    print("-" * 30)
    print("Valores das variáveis")
    for i in chute_variaveis:
        if str(i[0])[0] == 'V':
            print(f"{i[0]}: {i[1]:.4f} pu")
        else:
            print(f"{i[0]}: {i[1]} rad = {float(N(i[1] * 180 / pi)):.4f} [graus]")
    print("Valores das equações")
    for i in range(0, len(valores_calculados)):
        print(f"{nomes_equacoes[i]} = {valores_calculados[i]:.4f} pu")
# Display do resultado final
print("-" * 30)
print(f"Resultado final")
print("-" * 30)
print("Valores das variáveis")
for i in chute_variaveis:
    if str(i[0])[0] == 'V':
        print(f"{i[0]}: {i[1]:.4f} pu")
    else:
        print(f"{i[0]}: {i[1]:.4f} rad = {float(N(i[1] * 180 / pi)):.4f} [graus]")
print("Valores das equações")
for i in range(0, len(valores_calculados)):
    print(f"{nomes_equacoes[i]} = {valores_calculados[i]:.4f} pu")
print("-" * 30)
# EmpelTec Jr.
print("(c) 2024 - EmpelTec Jr. - Todos os direitos reservados")
