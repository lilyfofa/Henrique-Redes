from sympy import *


def FluxoDePotencia(n_barras, S_base, erro, Pg, Pl, Qg, Ql, Tensao, Fase, admitancias, alpha):
    for i in range(0, n_barras):
        Pg[i] *= alpha
        Qg[i] *= alpha
        Pl[i] *= alpha
        Ql[i] *= alpha
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
    V = Tensao
    Theta = Fase
    gbus = re(ybus)
    bbus = im(ybus)
    equacoes = []
    variaveis = []
    chute_variaveis = []
    for i in range(0, len(Tensao)):
        if Theta[i] == -1:
            Theta[i] = symbols(f"Theta{i + 1}")
            variaveis.append(symbols(f"Theta{i + 1}"))
            chute_variaveis.append((symbols(f"Theta{i + 1}"), 0))
    for i in range(0, len(Tensao)):
        if Tensao[i] == 0:
            V[i] = symbols(f"V{i + 1}")
            variaveis.append(symbols(f"V{i + 1}"))
            chute_variaveis.append((symbols(f"V{i + 1}"), 1))
    for i in range(0, n_barras):
        Pn = Pg[i] - Pl[i]
        if Pn != 0:
            expressao = Pn
            for j in range(0, n_barras):
                expressao -= V[i] * V[j] * (
                        gbus[i, j] * cos(Theta[i] - Theta[j]) + bbus[i, j] * sin(Theta[i] - Theta[j]))
            equacoes.append(expressao)
    for i in range(0, n_barras):
        Qn = Qg[i] - Ql[i]
        if Qn != 0:
            expressao = Qn
            for j in range(0, n_barras):
                expressao -= V[i] * V[j] * (
                        gbus[i, j] * sin(Theta[i] - Theta[j]) - bbus[i, j] * cos(Theta[i] - Theta[j]))
            equacoes.append(expressao)
    jacobiano = zeros(len(equacoes), len(variaveis))
    for i in range(0, len(equacoes)):
        for j in range(0, len(variaveis)):
            jacobiano[i, j] = diff(equacoes[i], variaveis[j])
    k = 0
    valores_calculados = zeros(len(equacoes), 1)
    for i in range(0, len(valores_calculados)):
        valores_calculados[i] = [equacoes[i].subs(chute_variaveis)]
    while max(abs(valores_calculados)) > erro:
        jac = jacobiano.subs(chute_variaveis)  # Cálculo do valor númerico do jacobiano
        delta = -1 * jac.inv() * valores_calculados  # Encontrando os ajustes
        novos_valores = []  # Matriz auxiliar pois uma tupla não pode ser alterada diretamente
        for i in range(0, len(chute_variaveis)):
            novos_valores.append((variaveis[i], chute_variaveis[i][1] + delta[i]))  # Ajustes dos chutes das variáveis
        chute_variaveis = novos_valores  # Substituição da tupla antiga
        for i in range(0, len(valores_calculados)):
            valores_calculados[i] = [equacoes[i].subs(chute_variaveis)]  # Novos valores para as equações
        k += 1
    V_resultados = V
    Theta_resultados = Theta
    for i in range(0, len(V_resultados)):
        for j in range(0, len(chute_variaveis)):
            if V_resultados[i] == chute_variaveis[j][0]:
                V_resultados[i] = chute_variaveis[j][1]
    for i in range(0, len(Theta_resultados)):
        for j in range(0, len(chute_variaveis)):
            if Theta_resultados[i] == chute_variaveis[j][0]:
                Theta_resultados[i] = chute_variaveis[j][1]
    Pk = []
    Qk = []
    for i in range(0, n_barras):
        expressao1 = 0
        expressao2 = 0
        for j in range(0, n_barras):
            expressao1 += V_resultados[i] * V_resultados[j] * (
                    gbus[i, j] * cos(Theta_resultados[i] - Theta_resultados[j]) +
                    bbus[i, j] * sin(Theta_resultados[i] - Theta_resultados[j]))
            expressao2 += V_resultados[i] * V_resultados[j] * (
                    gbus[i, j] * sin(Theta_resultados[i] - Theta_resultados[j])
                    - bbus[i, j] * cos(Theta_resultados[i] - Theta_resultados[j]))
        Pk.append(expressao1)
        Qk.append(expressao2)
    Vk, Vm, Thetak, Thetam, ykm, bshkm = symbols('Vk Vm Thetak Thetam ykm bshkm')
    Pkm = Vk * Vk * re(ykm) - Vk * Vm * re(ykm) * cos(Thetak - Thetam) - Vk * Vm * im(ykm) * sin(Thetak - Thetam)
    Qkm = -Vk * Vk * (im(ykm) + bshkm) + Vk * Vm * im(ykm) * cos(Thetak - Thetam) - Vk * Vm * re(ykm) * sin(
        Thetak - Thetam)
    Plinha = []
    Qlinha = []
    Slinha = []
    for carga in admitancias:
        valor1 = Pkm.subs(
            [(Vk, V[carga[2] - 1]), (Vm, V[carga[3] - 1]), (ykm, 1 / carga[0]), (Thetak, Theta[carga[2] - 1]),
             (Thetam, Theta[carga[3] - 1])])
        valor2 = Qkm.subs(
            [(Vk, V[carga[2] - 1]), (Vm, V[carga[3] - 1]), (ykm, 1 / carga[0]), (Thetak, Theta[carga[2] - 1]),
             (Thetam, Theta[carga[3] - 1]), (bshkm, carga[1] / 200)])
        valor3 = sqrt(valor1 ** 2 + valor2 ** 2)
        Plinha.append(valor1)
        Qlinha.append(valor2)
        Slinha.append(valor3)
    soma = 0
    Perdas = []
    for item in admitancias:
        barra1 = item[2] - 1
        barra2 = item[3] - 1
        carga = item[0]
        V1 = V[barra1] * cos(Theta[barra1]) + I * V[barra1] * sin(Theta[barra1])
        V2 = V[barra2] * cos(Theta[barra2]) + I * V[barra2] * sin(Theta[barra2])
        queda = abs(V1 - V2)
        P = re(queda ** 2 / carga)
        soma += P
        Perdas.append(P)
    return V_resultados, Theta_resultados, Pk, Qk, Plinha, Qlinha, Slinha, Perdas, soma
