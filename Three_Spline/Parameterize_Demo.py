
import numpy as np
import matplotlib.pyplot as plt
from Function import three_spline_coe, three_spline_curve, parameterize

if __name__ == "__main__" :
    ControlPoints = np.array([(1, 2), (2, 4), (3, 7), (4, 3), (5, 4), (6, 8)])
    a = 0
    b = 1
    t = np.linspace(a, b, num=101, endpoint=True)
    x = np.zeros(101)
    y = np.zeros(101)
    ControlPoints_x, ControlPoints_y = parameterize(ControlPoints, 0)
    '''计算系数'''
    x_2 = three_spline_coe(ControlPoints_x, 0, 0, 2)
    y_2 = three_spline_coe(ControlPoints_y, 0, 0, 2)
    '''计算x(t), y(t)'''
    for i in range(101):
        x[i] = three_spline_curve(ControlPoints_x, x_2, t[i])
        y[i] = three_spline_curve(ControlPoints_y, y_2, t[i])
    '''作图'''
    plt.figure()

    plt.plot(ControlPoints[:, 0], ControlPoints[:, 1], 'ro')
    plt.plot(x, y, 'k')

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Figure")

    plt.grid(True)
    plt.show()