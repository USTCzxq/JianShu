import numpy as np


def phi0(x):
    y = 1 - x
    return y


def phi1(x):
    y = x
    return y


def phi2(x):
    y = -np.power(x, 3) / 6 + np.power(x, 2) / 2 - x / 3
    return y


def phi3(x):
    y = np.power(x, 3) / 6 - x / 6
    return y


def parameterize(ControlPoints, index):
    x = ControlPoints[:, 0]
    y = ControlPoints[:, 1]
    num = len(x)
    n = num - 1
    t = np.zeros(num)
    ControlPoints_x = np.zeros((num, 2))
    ControlPoints_y = np.zeros((num, 2))
    if index == 0:
        # 均匀参数化
        for i in range(num):
            t[i] = i / n
        ControlPoints_x[:, 0] = t
        ControlPoints_x[:, 1] = x
        ControlPoints_y[:, 0] = t
        ControlPoints_y[:, 1] = y
    if index == 1:
        # 弦长参数化
        d = np.zeros(n)
        for i in range(n):
            d[i] = np.sqrt(np.power(x[i+1] - x[i], 2) + np.power(y[i+1] - y[i], 2))
        for i in range(1, num):
            t[i] = np.sum(d[0:i]) / np.sum(d)
        ControlPoints_x[:, 0] = t
        ControlPoints_x[:, 1] = x
        ControlPoints_y[:, 0] = t
        ControlPoints_y[:, 1] = y

    return ControlPoints_x, ControlPoints_y


def thomas(a, b, c, d):
    """
    求解三对角矩阵
    :param a:
    :param b:
    :param c:
    :param d:
    :return:
    """
    n = len(d)
    c_star = np.array([])
    d_star = np.array([])
    x = np.zeros(n)
    for i in range(n-1):
        if i == 0:
            ele = c[i] / b[i]
            c_star = np.append(c_star, ele)
        else:
            ele = c[i] / (b[i] - a[i] * c_star[i-1])
            c_star = np.append(c_star, ele)

    for i in range(n):
        if i == 0:
            ele = d[i] / b[i]
            d_star = np.append(d_star, ele)
        else:
            ele = (d[i] - a[i] * d_star[i-1]) / (b[i] - a[i] * c_star[i-1])
            d_star = np.append(d_star, ele)

    for i in range(n, 0, -1):
        if (i-1) == (n-1):
            x[i-1] = d_star[n-1]
        else:
            x[i - 1] = d_star[i-1] - c_star[i-1] * x[i]

    return x


def three_spline_coe(ControlPoints, bound0, boundn, index):
    """
    输入控制顶点ControlPoints，和边界条件，求出曲线y(x)中未知的参数, 其中index代表三种边界条件的索引
    :param ControlPoints:
    :param bound0:
    :param boundn:
    :param index:
    :return:
    """
    x = ControlPoints[:, 0]
    y = ControlPoints[:, 1]
    num = len(x)  # num代表控制点的个数
    a = np.zeros(num)
    b = np.zeros(num)
    c = np.zeros(num)
    d = np.zeros(num)
    n = num - 1
    h = np.zeros(n)
    m = np.zeros(n)

    for i in range(n):
        h[i] = x[i+1] - x[i]
    for i in range(n):
        m[i] = (y[i+1] - y[i]) / h[i]

    if index == 0:
        '''自由端点'''
        b[0] = 1
        b[n] = 1
        d[0] = 0
        d[n] = 0
    if index == 1:
        '''抛物条件'''
        b[0] = 1
        c[0] = -1
        a[n] = 1
        b[n] = -1
        d[0] = 0
        d[n] = 0
    if index == 2:
        b[0] = 2 * h[0]
        c[0] = h[0]
        a[n] = h[n-1]
        b[n] = 2 * h[n-1]
        d[0] = m[0] - bound0
        d[n] = boundn - m[n-1]

    for i in range(1, n):
        a[i] = h[i-1]
        b[i] = 2 * (h[i-1] + h[i])
        c[i] = h[i]
        d[i] = m[i] - m[i-1]

    d = d * 6

    y_2 = thomas(a, b, c, d)

    return y_2


def three_spline_curve(ControlPoints, coe, t):
    """
    输入控制顶点，边界条件，t；返回插值曲线在t处的函数值
    :param ControlPoints:
    :param coe:
    :param t:
    :return:
    """
    x = ControlPoints[:, 0]
    y = ControlPoints[:, 1]
    num = len(x)  # num代表控制点的个数
    n = num - 1
    h = np.zeros(n)
    value_y = 0

    for i in range(n):
        h[i] = x[i+1] - x[i]

    y_2 = coe
    for i in range(n):
        if (t >= x[i]) and (t <= x[i+1]):
            value_y = y[i] * phi0((t - x[i]) / h[i]) + y[i+1] * phi1((t - x[i]) / h[i]) + \
                    y_2[i] * np.power(h[i], 2) * phi2((t - x[i]) / h[i]) + y_2[i+1] * np.power(h[i], 2) * phi3((t - x[i]) / h[i])
            break

    return value_y
