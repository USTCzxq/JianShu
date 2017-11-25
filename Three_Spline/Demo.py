
import numpy as np
import matplotlib.pyplot as plt
from Function import three_spline_coe, three_spline_curve

if __name__ == "__main__" :
    ControlPoints = np.array([(1, 2), (2, 4), (3, 7), (4, 3), (5, 4), (6, 8)])
    a = np.max(ControlPoints[:, 0])
    b = np.min(ControlPoints[:, 0])
    t = np.linspace(a, b, num=101, endpoint=True)
    y = np.zeros(101)
    '''计算y(x)的系数'''
    y_2 = three_spline_coe(ControlPoints, 0, 0, 2)
    '''计算y(x)'''
    for i in range(101):
        y[i] = three_spline_curve(ControlPoints, y_2, t[i])
    '''作图'''
    plt.figure()

    plt.plot(ControlPoints[:, 0], ControlPoints[:, 1], 'ro')
    plt.plot(t, y, 'k')

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Figure")

    plt.grid(True)
    plt.show()