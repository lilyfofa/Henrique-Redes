import numpy as np

V1 = 1.06
V2 = 1.0489
V3 = 1.06
Pg2 = 1.0230
Pg3 = 1.0075

x = np.array([1, 1, 1, 0, 0, 0, 0, 0])

# V4 = x[0]
# V5 = x[1]
# V6 = x[2]
# θ2 = x[3]
# θ3 = x[4]
# θ4 = x[5]
# θ5 = x[6]
# θ6 = x[7]

ΔP2 = Pg2 - V1*V2*(4.0*np.sin(x[3]) - 2.0*np.cos(x[3])) - 9.3282508137742*V2**2 - V2*V3*(3.84615384615385*np.sin(x[3] - x[4]) - 0.769230769230769*np.cos(x[3] - x[4])) - V2*x[0]*(8.0*np.sin(x[3] - x[5]) - 4.0*np.cos(x[3] - x[5])) - V2*x[1]*(3.0*np.sin(x[3] - x[6]) - 1.0*np.cos(x[3] - x[6])) - V2*x[2]*(4.4543429844098*np.sin(x[3] - x[7]) - 1.55902004454343*np.cos(x[3] - x[7]))
ΔP3 = Pg3 - V2*V3*(-3.84615384615385*np.sin(x[3] - x[4]) - 0.769230769230769*np.cos(x[3] - x[4])) - 4.15572232645403*V3**2 - V3*x[1]*(3.17073170731707*np.sin(x[4] - x[6]) - 1.46341463414634*np.cos(x[4] - x[6])) - V3*x[2]*(9.61538461538461*np.sin(x[4] - x[7]) - 1.92307692307692*np.cos(x[4] - x[7]))
ΔP4 = -V1*x[0]*(4.70588235294118*np.sin(x[5]) - 1.17647058823529*np.cos(x[5])) - V2*x[0]*(-8.0*np.sin(x[3] - x[5]) - 4.0*np.cos(x[3] - x[5])) - 6.17647058823529*x[0]**2 - x[0]*x[1]*(2.0*np.sin(x[5] - x[6]) - 1.0*np.cos(x[5] - x[6])) - 0.9
ΔP5 = -V1*x[1]*(3.11203319502075*np.sin(x[6]) - 0.829875518672199*np.cos(x[6])) - V2*x[1]*(-3.0*np.sin(x[3] - x[6]) - 1.0*np.cos(x[3] - x[6])) - V3*x[1]*(-3.17073170731707*np.sin(x[4] - x[6]) - 1.46341463414634*np.cos(x[4] - x[6])) - x[0]*x[1]*(-2.0*np.sin(x[5] - x[6]) - 1.0*np.cos(x[5] - x[6])) - 5.29329015281854*x[1]**2 - x[1]*x[2]*(3.0*np.sin(x[6] - x[7]) - 1.0*np.cos(x[6] - x[7])) - 1
ΔP6 = -V2*x[2]*(-4.4543429844098*np.sin(x[3] - x[7]) - 1.55902004454343*np.cos(x[3] - x[7])) - V3*x[2]*(-9.61538461538461*np.sin(x[4] - x[7]) - 1.92307692307692*np.cos(x[4] - x[7])) - x[1]*x[2]*(-3.0*np.sin(x[6] - x[7]) - 1.0*np.cos(x[6] - x[7])) - 4.48209696762035*x[2]**2 - 0.9
ΔQ4 = -V1*x[0]*(-1.17647058823529*np.sin(x[5]) - 4.70588235294118*np.cos(x[5])) - V2*x[0]*(4.0*np.sin(x[3] - x[5]) - 8.0*np.cos(x[3] - x[5])) - 14.6358823529412*x[0]**2 - x[0]*x[1]*(-1.0*np.sin(x[5] - x[6]) - 2.0*np.cos(x[5] - x[6])) - 0.6
ΔQ5 = -V1*x[1]*(-0.829875518672199*np.sin(x[6]) - 3.11203319502075*np.cos(x[6])) - V2*x[1]*(1.0*np.sin(x[3] - x[6]) - 3.0*np.cos(x[3] - x[6])) - V3*x[1]*(1.46341463414634*np.sin(x[4] - x[6]) - 3.17073170731707*np.cos(x[4] - x[6])) - x[0]*x[1]*(1.0*np.sin(x[5] - x[6]) - 2.0*np.cos(x[5] - x[6])) - 14.1377649023378*x[1]**2 - x[1]*x[2]*(-1.0*np.sin(x[6] - x[7]) - 3.0*np.cos(x[6] - x[7])) - 0.7
ΔQ6 = -V2*x[2]*(1.55902004454343*np.sin(x[3] - x[7]) - 4.4543429844098*np.cos(x[3] - x[7])) - V3*x[2]*(1.92307692307692*np.sin(x[4] - x[7]) - 9.61538461538461*np.cos(x[4] - x[7])) - x[1]*x[2]*(1.0*np.sin(x[6] - x[7]) - 3.0*np.cos(x[6] - x[7])) - 17.0047275997944*x[2]**2 - 0.5

B = np.array([ΔP2, ΔP3, ΔP4, ΔP5, ΔP6, ΔQ4, ΔQ5, ΔQ6])

while np.max(abs(B)) > 0.0001:
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
                       x[3] - x[4])) - V3 * x[1] * (1.46341463414634 * np.sin(x[4] - x[6]) + 3.17073170731707 * np.cos(
                       x[4] - x[6])) - V3 * x[2] * (
                           1.92307692307692 * np.sin(x[4] - x[7]) + 9.61538461538461 * np.cos(x[4] - x[7])), 0,
                   -V3 * x[1] * (-1.46341463414634 * np.sin(x[4] - x[6]) - 3.17073170731707 * np.cos(x[4] - x[6])),
                   -V3 * x[2] * (-1.92307692307692 * np.sin(x[4] - x[7]) - 9.61538461538461 * np.cos(x[4] - x[7]))], [
                      -V1 * (4.70588235294118 * np.sin(x[5]) - 1.17647058823529 * np.cos(x[5])) - V2 * (
                              -8.0 * np.sin(x[3] - x[5]) - 4.0 * np.cos(
                          x[3] - x[5])) - 12.3529411764706 * x[0] - x[1] * (
                              2.0 * np.sin(x[5] - x[6]) - 1.0 * np.cos(x[5] - x[6])),
                      -x[0] * (2.0 * np.sin(x[5] - x[6]) - 1.0 * np.cos(x[5] - x[6])), 0,
                      -V2 * x[0] * (4.0 * np.sin(x[3] - x[5]) - 8.0 * np.cos(x[3] - x[5])), 0,
                      -V1 * x[0] * (1.17647058823529 * np.sin(x[5]) + 4.70588235294118 * np.cos(x[5])) - V2 * x[0] * (
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
                   -V3 * x[2] * (1.92307692307692 * np.sin(x[4] - x[7]) - 9.61538461538461 * np.cos(x[4] - x[7])), 0,
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
                      -V1 * x[0] * (4.70588235294118 * np.sin(x[5]) - 1.17647058823529 * np.cos(x[5])) - V2 * x[0] * (
                              -8.0 * np.sin(x[3] - x[5]) - 4.0 * np.cos(x[3] - x[5])) - x[0] * x[1] * (
                              2.0 * np.sin(x[5] - x[6]) - 1.0 * np.cos(x[5] - x[6])),
                      -x[0] * x[1] * (-2.0 * np.sin(x[5] - x[6]) + 1.0 * np.cos(x[5] - x[6])), 0],
                  [-x[1] * (1.0 * np.sin(x[5] - x[6]) - 2.0 * np.cos(x[5] - x[6])),
                   -V1 * (-0.829875518672199 * np.sin(x[6]) - 3.11203319502075 * np.cos(x[6])) - V2 * (
                           1.0 * np.sin(x[3] - x[6]) - 3.0 * np.cos(x[3] - x[6])) - V3 * (
                           1.46341463414634 * np.sin(x[4] - x[6]) - 3.17073170731707 * np.cos(
                       x[4] - x[6])) - x[0] * (
                           1.0 * np.sin(x[5] - x[6]) - 2.0 * np.cos(x[5] - x[6])) - 28.2755298046756 * x[1] - x[2] * (
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
                   -V3 * x[2] * (9.61538461538461 * np.sin(x[4] - x[7]) + 1.92307692307692 * np.cos(x[4] - x[7])), 0,
                   -x[1] * x[2] * (3.0 * np.sin(x[6] - x[7]) + 1.0 * np.cos(x[6] - x[7])), -V2 * x[2] * (
                           -4.4543429844098 * np.sin(x[3] - x[7]) - 1.55902004454343 * np.cos(
                       x[3] - x[7])) - V3 * x[2] * (
                           -9.61538461538461 * np.sin(x[4] - x[7]) - 1.92307692307692 * np.cos(
                       x[4] - x[7])) - x[1] * x[2] * (-3.0 * np.sin(x[6] - x[7]) - 1.0 * np.cos(x[6] - x[7]))]])
    Δx = np.linalg.solve(A, -B)
    x = x + Δx
    ΔP2 = Pg2 - V1*V2*(4.0*np.sin(x[3]) - 2.0*np.cos(x[3])) - 9.3282508137742*V2**2 - V2*V3*(3.84615384615385*np.sin(x[3] - x[4]) - 0.769230769230769*np.cos(x[3] - x[4])) - V2*x[0]*(8.0*np.sin(x[3] - x[5]) - 4.0*np.cos(x[3] - x[5])) - V2*x[1]*(3.0*np.sin(x[3] - x[6]) - 1.0*np.cos(x[3] - x[6])) - V2*x[2]*(4.4543429844098*np.sin(x[3] - x[7]) - 1.55902004454343*np.cos(x[3] - x[7]))
    ΔP3 = Pg3 - V2*V3*(-3.84615384615385*np.sin(x[3] - x[4]) - 0.769230769230769*np.cos(x[3] - x[4])) - 4.15572232645403*V3**2 - V3*x[1]*(3.17073170731707*np.sin(x[4] - x[6]) - 1.46341463414634*np.cos(x[4] - x[6])) - V3*x[2]*(9.61538461538461*np.sin(x[4] - x[7]) - 1.92307692307692*np.cos(x[4] - x[7]))
    ΔP4 = -V1*x[0]*(4.70588235294118*np.sin(x[5]) - 1.17647058823529*np.cos(x[5])) - V2*x[0]*(-8.0*np.sin(x[3] - x[5]) - 4.0*np.cos(x[3] - x[5])) - 6.17647058823529*x[0]**2 - x[0]*x[1]*(2.0*np.sin(x[5] - x[6]) - 1.0*np.cos(x[5] - x[6])) - 0.9
    ΔP5 = -V1*x[1]*(3.11203319502075*np.sin(x[6]) - 0.829875518672199*np.cos(x[6])) - V2*x[1]*(-3.0*np.sin(x[3] - x[6]) - 1.0*np.cos(x[3] - x[6])) - V3*x[1]*(-3.17073170731707*np.sin(x[4] - x[6]) - 1.46341463414634*np.cos(x[4] - x[6])) - x[0]*x[1]*(-2.0*np.sin(x[5] - x[6]) - 1.0*np.cos(x[5] - x[6])) - 5.29329015281854*x[1]**2 - x[1]*x[2]*(3.0*np.sin(x[6] - x[7]) - 1.0*np.cos(x[6] - x[7])) - 1
    ΔP6 = -V2*x[2]*(-4.4543429844098*np.sin(x[3] - x[7]) - 1.55902004454343*np.cos(x[3] - x[7])) - V3*x[2]*(-9.61538461538461*np.sin(x[4] - x[7]) - 1.92307692307692*np.cos(x[4] - x[7])) - x[1]*x[2]*(-3.0*np.sin(x[6] - x[7]) - 1.0*np.cos(x[6] - x[7])) - 4.48209696762035*x[2]**2 - 0.9
    ΔQ4 = -V1*x[0]*(-1.17647058823529*np.sin(x[5]) - 4.70588235294118*np.cos(x[5])) - V2*x[0]*(4.0*np.sin(x[3] - x[5]) - 8.0*np.cos(x[3] - x[5])) - 14.6358823529412*x[0]**2 - x[0]*x[1]*(-1.0*np.sin(x[5] - x[6]) - 2.0*np.cos(x[5] - x[6])) - 0.6
    ΔQ5 = -V1*x[1]*(-0.829875518672199*np.sin(x[6]) - 3.11203319502075*np.cos(x[6])) - V2*x[1]*(1.0*np.sin(x[3] - x[6]) - 3.0*np.cos(x[3] - x[6])) - V3*x[1]*(1.46341463414634*np.sin(x[4] - x[6]) - 3.17073170731707*np.cos(x[4] - x[6])) - x[0]*x[1]*(1.0*np.sin(x[5] - x[6]) - 2.0*np.cos(x[5] - x[6])) - 14.1377649023378*x[1]**2 - x[1]*x[2]*(-1.0*np.sin(x[6] - x[7]) - 3.0*np.cos(x[6] - x[7])) - 0.7
    ΔQ6 = -V2*x[2]*(1.55902004454343*np.sin(x[3] - x[7]) - 4.4543429844098*np.cos(x[3] - x[7])) - V3*x[2]*(1.92307692307692*np.sin(x[4] - x[7]) - 9.61538461538461*np.cos(x[4] - x[7])) - x[1]*x[2]*(1.0*np.sin(x[6] - x[7]) - 3.0*np.cos(x[6] - x[7])) - 17.0047275997944*x[2]**2 - 0.5
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

V = [V1, V2, V3, x[0], x[1], x[2]]
θ = [0, x[3], x[4], x[5], x[6], x[7]]
soma = 0

for impedancia, _, barra1, barra2 in dlin:
    Va = complex(V[barra1-1] * np.cos(θ[barra1-1]), V[barra1-1] * np.sin(θ[barra1-1]))
    Vb = complex(V[barra2-1] * np.cos(θ[barra2-1]), V[barra2-1] * np.sin(θ[barra2-1]))
    queda = abs(Va - Vb)
    soma += np.real(queda ** 2 / impedancia)
print(f'Resultado final:')
nomes = ['V4', 'V5', 'V6', 'θ2', 'θ3', 'θ4', 'θ5', 'θ6']
unidade= ['pu', 'pu', 'pu', 'rad', 'rad', 'rad', 'rad', 'rad']
for i in range(len(x)):
    print(f'{nomes[i]} = {x[i]:.4f} {unidade[i]}')
print(f'Perdas: {soma:.4f} pu')
print('(c) 2024 - EmpelTec Jr - Todos os direitos reservados')