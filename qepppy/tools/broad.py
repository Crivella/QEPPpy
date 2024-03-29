import numpy as np


def _gaussian(x, broad):
    dx = x[1] - x[0]
    sigma = broad
    twoSigmaSq = 2.0 * sigma**2
    norm = np.sqrt(2 * np.pi) * sigma / dx
    #Divided for dx to take into account implementation as sum instead of integral (volume element)

    g = np.arange(-10*sigma, 10*sigma, dx)
    g = np.exp(-(g**2) / twoSigmaSq)
    g /= norm
    return g

def _lorentz(x, broad):
    dx = x[1] - x[0]
    gamma = broad
    gammaSq = gamma**2
    norm = np.pi / gamma / dx
    #Divided for dx to take into account implementation as sum instead of integral (volume element)

    g = np.arange(-15*gamma, 15*gamma, dx)
    g = 1/(g**2 + gammaSq)
    g /= norm
    return g

def _check_and_pad_(data,deg):
    x = data[:,0]
    y = data[:,1:]

    dx = x[1] - x[0]
    if np.any((x[1:
        ] - x[:-1]) - dx > 1E-4):
        raise ValueError('X-axis data must be uniformely distributed.')

    padding = int(15*deg/dx)
    x = np.pad(x, (padding,padding), 'linear_ramp', end_values=(x[0]-15*deg, x[-1] + 15*deg))
    y = np.pad(y, ((padding,padding),(0,0)), 'edge')

    return x,y

def _get_conv_function_(t):
    if t == 'gauss':
        return _gaussian
    if t == 'lorentz':
        return _lorentz
    else:
        raise NotImplemented(f'Invalid convolution function {t}.')


def broad(data, t='gauss', deg=0.1, axis=0):
    if axis:
        data = data.T

    res = np.zeros(data.shape)
    res[:,0] = data[:,0]

    x,y = _check_and_pad_(data,deg)

    conv = _get_conv_function_(t)(x,deg)
    w = np.where((res[0,0] <= x) & (x <= res[-1,0]))
    for n,yi in enumerate(y.T):
        new = np.convolve(yi, conv, mode='same')
        res[:,n+1] = new[w]

    if axis:
        res = res.T

    return res

def main():
    import sys
    argc = len(sys.argv)
    if not 2<=argc<=5 or sys.argv[1] == 'help':
        print('Incorrect use. Pleas pass arguments:'
            "\n\t'fname'\t (comma separated),"
            "\n\t'type\t(gauss/lorentz, default=gauss)',"
            "\n\t'broad\t (default=0.1)'"
            "\n\t'oname\t(optional) (comma separated)'")
        exit()

    t = 'gauss'
    if argc >= 3:
        t = str(sys.argv[2])
    deg = 0.1
    if argc >= 4:
        deg = float(sys.argv[3])
    oname = None
    if argc >= 5:
        oname = str(sys.argv[4]).split(',')
    for n, name in enumerate(str(sys.argv[1]).split(',')):
        print(f"Applying '{t.upper()}' broadening of '{deg}' to '{name}'")
        data = np.loadtxt(name, usecols=None)

        b = broad(data, t, deg)

        if oname:
            n_name = oname[n]
        else:
            n_name = f'{name}_{deg}'
        print(f"Saving broadened data to '{n_name}'")
        np.savetxt(n_name, b)

if __name__ == '__main__':
    main()
