from sympy import *
from time import time
import numpy as np


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
    jacobiano = zeros(len(equacoes), len(equacoes))
    print(variaveis, len(variaveis), len(equacoes))
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


def FluxoDeCarga(posicao):
    V1 = posicao[0]
    V2 = posicao[1]
    V3 = posicao[2]
    Pg2 = posicao[3]
    Pg3 = posicao[4]
    x = np.array([1, 1, 1, 0, 0, 0, 0, 0])
    ΔP2 = Pg2 - V1*V2*(4.0*np.sin(x[3]) - 2.0*np.cos(x[3])) - 9.3282508137742*V2**2 - V2*V3*(3.84615384615385*np.sin(x[3] - x[4]) - 0.769230769230769*np.cos(x[3] - x[4])) - V2*x[0]*(8.0*np.sin(x[3] - x[5]) - 4.0*np.cos(x[3] - x[5])) - V2*x[1]*(3.0*np.sin(x[3] - x[6]) - 1.0*np.cos(x[3] - x[6])) - V2*x[2]*(4.4543429844098*np.sin(x[3] - x[7]) - 1.55902004454343*np.cos(x[3] - x[7]))
    ΔP3 = Pg3 - V2*V3*(-3.84615384615385*np.sin(x[3] - x[4]) - 0.769230769230769*np.cos(x[3] - x[4])) - 4.15572232645403*V3**2 - V3*x[1]*(3.17073170731707*np.sin(x[4] - x[6]) - 1.46341463414634*np.cos(x[4] - x[6])) - V3*x[2]*(9.61538461538461*np.sin(x[4] - x[7]) - 1.92307692307692*np.cos(x[4] - x[7]))
    ΔP4 = -V1*x[0]*(4.70588235294118*np.sin(x[5]) - 1.17647058823529*np.cos(x[5])) - V2*x[0]*(-8.0*np.sin(x[3] - x[5]) - 4.0*np.cos(x[3] - x[5])) - 6.17647058823529*x[0]**2 - x[0]*x[1]*(2.0*np.sin(x[5] - x[6]) - 1.0*np.cos(x[5] - x[6])) - 0.9
    ΔP5 = -V1*x[1]*(3.11203319502075*np.sin(x[6]) - 0.829875518672199*np.cos(x[6])) - V2*x[1]*(-3.0*np.sin(x[3] - x[6]) - 1.0*np.cos(x[3] - x[6])) - V3*x[1]*(-3.17073170731707*np.sin(x[4] - x[6]) - 1.46341463414634*np.cos(x[4] - x[6])) - x[0]*x[1]*(-2.0*np.sin(x[5] - x[6]) - 1.0*np.cos(x[5] - x[6])) - 5.29329015281854*x[1]**2 - x[1]*x[2]*(3.0*np.sin(x[6] - x[7]) - 1.0*np.cos(x[6] - x[7])) - 1
    ΔP6 = -V2*x[2]*(-4.4543429844098*np.sin(x[3] - x[7]) - 1.55902004454343*np.cos(x[3] - x[7])) - V3*x[2]*(-9.61538461538461*np.sin(x[4] - x[7]) - 1.92307692307692*np.cos(x[4] - x[7])) - x[1]*x[2]*(-3.0*np.sin(x[6] - x[7]) - 1.0*np.cos(x[6] - x[7])) - 4.48209696762035*x[2]**2 - 0.9
    ΔQ4 = -V1*x[0]*(-1.17647058823529*np.sin(x[5]) - 4.70588235294118*np.cos(x[5])) - V2*x[0]*(4.0*np.sin(x[3] - x[5]) - 8.0*np.cos(x[3] - x[5])) - 14.6358823529412*x[0]**2 - x[0]*x[1]*(-1.0*np.sin(x[5] - x[6]) - 2.0*np.cos(x[5] - x[6])) - 0.6
    ΔQ5 = -V1*x[1]*(-0.829875518672199*np.sin(x[6]) - 3.11203319502075*np.cos(x[6])) - V2*x[1]*(1.0*np.sin(x[3] - x[6]) - 3.0*np.cos(x[3] - x[6])) - V3*x[1]*(1.46341463414634*np.sin(x[4] - x[6]) - 3.17073170731707*np.cos(x[4] - x[6])) - x[0]*x[1]*(1.0*np.sin(x[5] - x[6]) - 2.0*np.cos(x[5] - x[6])) - 14.1377649023378*x[1]**2 - x[1]*x[2]*(-1.0*np.sin(x[6] - x[7]) - 3.0*np.cos(x[6] - x[7])) - 0.7
    ΔQ6 = -V2*x[2]*(1.55902004454343*np.sin(x[3] - x[7]) - 4.4543429844098*np.cos(x[3] - x[7])) - V3*x[2]*(1.92307692307692*np.sin(x[4] - x[7]) - 9.61538461538461*np.cos(x[4] - x[7])) - x[1]*x[2]*(1.0*np.sin(x[6] - x[7]) - 3.0*np.cos(x[6] - x[7])) - 17.0047275997944*x[2]**2 - 0.5
    B = np.array([ΔP2, ΔP3, ΔP4, ΔP5, ΔP6, ΔQ4, ΔQ5, ΔQ6])
    while abs(np.max(B)) > 0.0001:
        A = np.array([[-V2 * (8.0 * np.sin(x[3] - x[5]) - 4.0 * np.cos(x[3] - x[5])),
                       -V2 * (3.0 * np.sin(x[3] - x[6]) - 1.0 * np.cos(x[3] - x[6])),
                       -V2 * (4.4543429844098 * np.sin(x[3] - x[7]) - 1.55902004454343 * np.cos(x[3] - x[7])),
                       -V1 * V2 * (2.0 * np.sin(x[3]) + 4.0 * np.cos(x[3])) - V2 * V3 * (
                               0.769230769230769 * np.sin(x[3] - x[4]) + 3.84615384615385 * np.cos(
                           x[3] - x[4])) - V2 * x[0] * (
                               4.0 * np.sin(x[3] - x[5]) + 8.0 * np.cos(x[3] - x[5])) - V2 * x[1] * (
                               1.0 * np.sin(x[3] - x[6]) + 3.0 * np.cos(x[3] - x[6])) - V2 * x[2] * (
                               1.55902004454343 * np.sin(x[3] - x[7]) + 4.4543429844098 * np.cos(x[3] - x[7])),
                       -V2 * V3 * (-0.769230769230769 * np.sin(x[3] - x[4]) - 3.84615384615385 * np.cos(x[3] - x[4])),
                       -V2 * x[0] * (-4.0 * np.sin(x[3] - x[5]) - 8.0 * np.cos(x[3] - x[5])),
                       -V2 * x[1] * (-1.0 * np.sin(x[3] - x[6]) - 3.0 * np.cos(x[3] - x[6])),
                       -V2 * x[2] * (-1.55902004454343 * np.sin(x[3] - x[7]) - 4.4543429844098 * np.cos(x[3] - x[7]))],
                      [0, -V3 * (3.17073170731707 * np.sin(x[4] - x[6]) - 1.46341463414634 * np.cos(x[4] - x[6])),
                       -V3 * (9.61538461538461 * np.sin(x[4] - x[7]) - 1.92307692307692 * np.cos(x[4] - x[7])),
                       -V2 * V3 * (0.769230769230769 * np.sin(x[3] - x[4]) - 3.84615384615385 * np.cos(x[3] - x[4])),
                       -V2 * V3 * (-0.769230769230769 * np.sin(x[3] - x[4]) + 3.84615384615385 * np.cos(
                           x[3] - x[4])) - V3 * x[1] * (
                                   1.46341463414634 * np.sin(x[4] - x[6]) + 3.17073170731707 * np.cos(
                               x[4] - x[6])) - V3 * x[2] * (
                               1.92307692307692 * np.sin(x[4] - x[7]) + 9.61538461538461 * np.cos(x[4] - x[7])), 0,
                       -V3 * x[1] * (-1.46341463414634 * np.sin(x[4] - x[6]) - 3.17073170731707 * np.cos(x[4] - x[6])),
                       -V3 * x[2] * (-1.92307692307692 * np.sin(x[4] - x[7]) - 9.61538461538461 * np.cos(x[4] - x[7]))],
                      [
                          -V1 * (4.70588235294118 * np.sin(x[5]) - 1.17647058823529 * np.cos(x[5])) - V2 * (
                                  -8.0 * np.sin(x[3] - x[5]) - 4.0 * np.cos(
                              x[3] - x[5])) - 12.3529411764706 * x[0] - x[1] * (
                                  2.0 * np.sin(x[5] - x[6]) - 1.0 * np.cos(x[5] - x[6])),
                          -x[0] * (2.0 * np.sin(x[5] - x[6]) - 1.0 * np.cos(x[5] - x[6])), 0,
                          -V2 * x[0] * (4.0 * np.sin(x[3] - x[5]) - 8.0 * np.cos(x[3] - x[5])), 0,
                          -V1 * x[0] * (1.17647058823529 * np.sin(x[5]) + 4.70588235294118 * np.cos(x[5])) - V2 * x[
                              0] * (
                                  -4.0 * np.sin(x[3] - x[5]) + 8.0 * np.cos(x[3] - x[5])) - x[0] * x[1] * (
                                  1.0 * np.sin(x[5] - x[6]) + 2.0 * np.cos(x[5] - x[6])),
                          -x[0] * x[1] * (-1.0 * np.sin(x[5] - x[6]) - 2.0 * np.cos(x[5] - x[6])), 0],
                      [-x[1] * (-2.0 * np.sin(x[5] - x[6]) - 1.0 * np.cos(x[5] - x[6])),
                       -V1 * (3.11203319502075 * np.sin(x[6]) - 0.829875518672199 * np.cos(x[6])) - V2 * (
                               -3.0 * np.sin(x[3] - x[6]) - 1.0 * np.cos(x[3] - x[6])) - V3 * (
                               -3.17073170731707 * np.sin(x[4] - x[6]) - 1.46341463414634 * np.cos(
                           x[4] - x[6])) - x[0] * (-2.0 * np.sin(x[5] - x[6]) - 1.0 * np.cos(
                           x[5] - x[6])) - 10.5865803056371 * x[1] - x[2] * (
                               3.0 * np.sin(x[6] - x[7]) - 1.0 * np.cos(x[6] - x[7])),
                       -x[1] * (3.0 * np.sin(x[6] - x[7]) - 1.0 * np.cos(x[6] - x[7])),
                       -V2 * x[1] * (1.0 * np.sin(x[3] - x[6]) - 3.0 * np.cos(x[3] - x[6])),
                       -V3 * x[1] * (1.46341463414634 * np.sin(x[4] - x[6]) - 3.17073170731707 * np.cos(x[4] - x[6])),
                       -x[0] * x[1] * (1.0 * np.sin(x[5] - x[6]) - 2.0 * np.cos(x[5] - x[6])),
                       -V1 * x[1] * (0.829875518672199 * np.sin(x[6]) + 3.11203319502075 * np.cos(x[6])) - V2 * x[1] * (
                               -1.0 * np.sin(x[3] - x[6]) + 3.0 * np.cos(x[3] - x[6])) - V3 * x[1] * (
                               -1.46341463414634 * np.sin(x[4] - x[6]) + 3.17073170731707 * np.cos(
                           x[4] - x[6])) - x[0] * x[1] * (
                               -1.0 * np.sin(x[5] - x[6]) + 2.0 * np.cos(x[5] - x[6])) - x[1] * x[2] * (
                               1.0 * np.sin(x[6] - x[7]) + 3.0 * np.cos(x[6] - x[7])),
                       -x[1] * x[2] * (-1.0 * np.sin(x[6] - x[7]) - 3.0 * np.cos(x[6] - x[7]))],
                      [0, -x[2] * (-3.0 * np.sin(x[6] - x[7]) - 1.0 * np.cos(x[6] - x[7])),
                       -V2 * (-4.4543429844098 * np.sin(x[3] - x[7]) - 1.55902004454343 * np.cos(x[3] - x[7])) - V3 * (
                               -9.61538461538461 * np.sin(x[4] - x[7]) - 1.92307692307692 * np.cos(
                           x[4] - x[7])) - x[1] * (
                               -3.0 * np.sin(x[6] - x[7]) - 1.0 * np.cos(x[6] - x[7])) - 8.96419393524071 * x[2],
                       -V2 * x[2] * (1.55902004454343 * np.sin(x[3] - x[7]) - 4.4543429844098 * np.cos(x[3] - x[7])),
                       -V3 * x[2] * (1.92307692307692 * np.sin(x[4] - x[7]) - 9.61538461538461 * np.cos(x[4] - x[7])),
                       0,
                       -x[1] * x[2] * (1.0 * np.sin(x[6] - x[7]) - 3.0 * np.cos(x[6] - x[7])), -V2 * x[2] * (
                               -1.55902004454343 * np.sin(x[3] - x[7]) + 4.4543429844098 * np.cos(
                           x[3] - x[7])) - V3 * x[2] * (
                               -1.92307692307692 * np.sin(x[4] - x[7]) + 9.61538461538461 * np.cos(
                           x[4] - x[7])) - x[1] * x[2] * (-1.0 * np.sin(x[6] - x[7]) + 3.0 * np.cos(x[6] - x[7]))], [
                          -V1 * (-1.17647058823529 * np.sin(x[5]) - 4.70588235294118 * np.cos(x[5])) - V2 * (
                                  4.0 * np.sin(x[3] - x[5]) - 8.0 * np.cos(
                              x[3] - x[5])) - 29.2717647058824 * x[0] - x[1] * (
                                  -1.0 * np.sin(x[5] - x[6]) - 2.0 * np.cos(x[5] - x[6])),
                          -x[0] * (-1.0 * np.sin(x[5] - x[6]) - 2.0 * np.cos(x[5] - x[6])), 0,
                          -V2 * x[0] * (8.0 * np.sin(x[3] - x[5]) + 4.0 * np.cos(x[3] - x[5])), 0,
                          -V1 * x[0] * (4.70588235294118 * np.sin(x[5]) - 1.17647058823529 * np.cos(x[5])) - V2 * x[
                              0] * (
                                  -8.0 * np.sin(x[3] - x[5]) - 4.0 * np.cos(x[3] - x[5])) - x[0] * x[1] * (
                                  2.0 * np.sin(x[5] - x[6]) - 1.0 * np.cos(x[5] - x[6])),
                          -x[0] * x[1] * (-2.0 * np.sin(x[5] - x[6]) + 1.0 * np.cos(x[5] - x[6])), 0],
                      [-x[1] * (1.0 * np.sin(x[5] - x[6]) - 2.0 * np.cos(x[5] - x[6])),
                       -V1 * (-0.829875518672199 * np.sin(x[6]) - 3.11203319502075 * np.cos(x[6])) - V2 * (
                               1.0 * np.sin(x[3] - x[6]) - 3.0 * np.cos(x[3] - x[6])) - V3 * (
                               1.46341463414634 * np.sin(x[4] - x[6]) - 3.17073170731707 * np.cos(
                           x[4] - x[6])) - x[0] * (
                               1.0 * np.sin(x[5] - x[6]) - 2.0 * np.cos(x[5] - x[6])) - 28.2755298046756 * x[1] - x[
                           2] * (
                               -1.0 * np.sin(x[6] - x[7]) - 3.0 * np.cos(x[6] - x[7])),
                       -x[1] * (-1.0 * np.sin(x[6] - x[7]) - 3.0 * np.cos(x[6] - x[7])),
                       -V2 * x[1] * (3.0 * np.sin(x[3] - x[6]) + 1.0 * np.cos(x[3] - x[6])),
                       -V3 * x[1] * (3.17073170731707 * np.sin(x[4] - x[6]) + 1.46341463414634 * np.cos(x[4] - x[6])),
                       -x[0] * x[1] * (2.0 * np.sin(x[5] - x[6]) + 1.0 * np.cos(x[5] - x[6])),
                       -V1 * x[1] * (3.11203319502075 * np.sin(x[6]) - 0.829875518672199 * np.cos(x[6])) - V2 * x[1] * (
                               -3.0 * np.sin(x[3] - x[6]) - 1.0 * np.cos(x[3] - x[6])) - V3 * x[1] * (
                               -3.17073170731707 * np.sin(x[4] - x[6]) - 1.46341463414634 * np.cos(
                           x[4] - x[6])) - x[0] * x[1] * (
                               -2.0 * np.sin(x[5] - x[6]) - 1.0 * np.cos(x[5] - x[6])) - x[1] * x[2] * (
                               3.0 * np.sin(x[6] - x[7]) - 1.0 * np.cos(x[6] - x[7])),
                       -x[1] * x[2] * (-3.0 * np.sin(x[6] - x[7]) + 1.0 * np.cos(x[6] - x[7]))],
                      [0, -x[2] * (1.0 * np.sin(x[6] - x[7]) - 3.0 * np.cos(x[6] - x[7])),
                       -V2 * (1.55902004454343 * np.sin(x[3] - x[7]) - 4.4543429844098 * np.cos(x[3] - x[7])) - V3 * (
                               1.92307692307692 * np.sin(x[4] - x[7]) - 9.61538461538461 * np.cos(
                           x[4] - x[7])) - x[1] * (
                               1.0 * np.sin(x[6] - x[7]) - 3.0 * np.cos(x[6] - x[7])) - 34.0094551995888 * x[2],
                       -V2 * x[2] * (4.4543429844098 * np.sin(x[3] - x[7]) + 1.55902004454343 * np.cos(x[3] - x[7])),
                       -V3 * x[2] * (9.61538461538461 * np.sin(x[4] - x[7]) + 1.92307692307692 * np.cos(x[4] - x[7])),
                       0,
                       -x[1] * x[2] * (3.0 * np.sin(x[6] - x[7]) + 1.0 * np.cos(x[6] - x[7])), -V2 * x[2] * (
                               -4.4543429844098 * np.sin(x[3] - x[7]) - 1.55902004454343 * np.cos(
                           x[3] - x[7])) - V3 * x[2] * (
                               -9.61538461538461 * np.sin(x[4] - x[7]) - 1.92307692307692 * np.cos(
                           x[4] - x[7])) - x[1] * x[2] * (-3.0 * np.sin(x[6] - x[7]) - 1.0 * np.cos(x[6] - x[7]))]])
        Δx = np.linalg.solve(A, -B)
        x = x + Δx
        ΔP2 = Pg2 - V1 * V2 * (4.0 * np.sin(x[3]) - 2.0 * np.cos(x[3])) - 9.3282508137742 * V2 ** 2 - V2 * V3 * (
                    3.84615384615385 * np.sin(x[3] - x[4]) - 0.769230769230769 * np.cos(x[3] - x[4])) - V2 * x[0] * (
                          8.0 * np.sin(x[3] - x[5]) - 4.0 * np.cos(x[3] - x[5])) - V2 * x[1] * (
                          3.0 * np.sin(x[3] - x[6]) - 1.0 * np.cos(x[3] - x[6])) - V2 * x[2] * (
                          4.4543429844098 * np.sin(x[3] - x[7]) - 1.55902004454343 * np.cos(x[3] - x[7]))
        ΔP3 = Pg3 - V2 * V3 * (-3.84615384615385 * np.sin(x[3] - x[4]) - 0.769230769230769 * np.cos(
            x[3] - x[4])) - 4.15572232645403 * V3 ** 2 - V3 * x[1] * (
                          3.17073170731707 * np.sin(x[4] - x[6]) - 1.46341463414634 * np.cos(x[4] - x[6])) - V3 * x[
                  2] * (9.61538461538461 * np.sin(x[4] - x[7]) - 1.92307692307692 * np.cos(x[4] - x[7]))
        ΔP4 = -V1 * x[0] * (4.70588235294118 * np.sin(x[5]) - 1.17647058823529 * np.cos(x[5])) - V2 * x[0] * (
                    -8.0 * np.sin(x[3] - x[5]) - 4.0 * np.cos(x[3] - x[5])) - 6.17647058823529 * x[0] ** 2 - x[0] * x[
                  1] * (2.0 * np.sin(x[5] - x[6]) - 1.0 * np.cos(x[5] - x[6])) - 0.9
        ΔP5 = -V1 * x[1] * (3.11203319502075 * np.sin(x[6]) - 0.829875518672199 * np.cos(x[6])) - V2 * x[1] * (
                    -3.0 * np.sin(x[3] - x[6]) - 1.0 * np.cos(x[3] - x[6])) - V3 * x[1] * (
                          -3.17073170731707 * np.sin(x[4] - x[6]) - 1.46341463414634 * np.cos(x[4] - x[6])) - x[0] * x[
                  1] * (-2.0 * np.sin(x[5] - x[6]) - 1.0 * np.cos(x[5] - x[6])) - 5.29329015281854 * x[1] ** 2 - x[1] * \
              x[2] * (3.0 * np.sin(x[6] - x[7]) - 1.0 * np.cos(x[6] - x[7])) - 1
        ΔP6 = -V2 * x[2] * (-4.4543429844098 * np.sin(x[3] - x[7]) - 1.55902004454343 * np.cos(x[3] - x[7])) - V3 * x[
            2] * (-9.61538461538461 * np.sin(x[4] - x[7]) - 1.92307692307692 * np.cos(x[4] - x[7])) - x[1] * x[2] * (
                          -3.0 * np.sin(x[6] - x[7]) - 1.0 * np.cos(x[6] - x[7])) - 4.48209696762035 * x[2] ** 2 - 0.9
        ΔQ4 = -V1 * x[0] * (-1.17647058823529 * np.sin(x[5]) - 4.70588235294118 * np.cos(x[5])) - V2 * x[0] * (
                    4.0 * np.sin(x[3] - x[5]) - 8.0 * np.cos(x[3] - x[5])) - 14.6358823529412 * x[0] ** 2 - x[0] * x[
                  1] * (-1.0 * np.sin(x[5] - x[6]) - 2.0 * np.cos(x[5] - x[6])) - 0.6
        ΔQ5 = -V1 * x[1] * (-0.829875518672199 * np.sin(x[6]) - 3.11203319502075 * np.cos(x[6])) - V2 * x[1] * (
                    1.0 * np.sin(x[3] - x[6]) - 3.0 * np.cos(x[3] - x[6])) - V3 * x[1] * (
                          1.46341463414634 * np.sin(x[4] - x[6]) - 3.17073170731707 * np.cos(x[4] - x[6])) - x[0] * x[
                  1] * (1.0 * np.sin(x[5] - x[6]) - 2.0 * np.cos(x[5] - x[6])) - 14.1377649023378 * x[1] ** 2 - x[1] * \
              x[2] * (-1.0 * np.sin(x[6] - x[7]) - 3.0 * np.cos(x[6] - x[7])) - 0.7
        ΔQ6 = -V2 * x[2] * (1.55902004454343 * np.sin(x[3] - x[7]) - 4.4543429844098 * np.cos(x[3] - x[7])) - V3 * x[
            2] * (1.92307692307692 * np.sin(x[4] - x[7]) - 9.61538461538461 * np.cos(x[4] - x[7])) - x[1] * x[2] * (
                          1.0 * np.sin(x[6] - x[7]) - 3.0 * np.cos(x[6] - x[7])) - 17.0047275997944 * x[2] ** 2 - 0.5
        B = np.array([ΔP2, ΔP3, ΔP4, ΔP5, ΔP6, ΔQ4, ΔQ5, ΔQ6])
    dlin = [(0.1 + 1j * 0.2, 4, 1, 2),
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
    V = np.array([V1, V2, V3, x[0], x[1], x[2]])
    θ = np.array([0, x[3], x[4], x[5], x[6], x[7]])
    soma = 0
    for impedancia, _, barra1, barra2 in dlin:
        Va = complex(V[barra1 - 1] * np.cos(θ[barra1 - 1]), V[barra1 - 1] * np.sin(θ[barra1 - 1]))
        Vb = complex(V[barra2 - 1] * np.cos(θ[barra2 - 1]), V[barra2 - 1] * np.sin(θ[barra2 - 1]))
        queda = np.absolute(Va - Vb)
        soma += np.real(queda ** 2 / impedancia)
    return soma
