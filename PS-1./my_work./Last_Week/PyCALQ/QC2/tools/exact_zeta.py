import os
import math
import cmath
import numpy as np
import scipy.special
import scipy.integrate
import warnings   # add this


def B(q, gamma, l, precision, verbose):
    if l != 0:
        return 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        a = 2. * sph_harm(0, 0, 0.0, 0.0) * gamma * math.pow(math.pi, 1.5)
        dawson = -1j * np.sqrt(math.pi) * np.exp(-q) * scipy.special.erf(1j * cmath.sqrt(q)) / 2.
        b = q * 2. * np.exp(q) * dawson / cmath.sqrt(q)
        c = math.exp(q)
        res = a * (b - c)
        if verbose:
            print('Term B:', res)
        return res
def Z(q2, gamma=1.0, l=0, m=0, d=np.array([0.,0.,0.]),
      m_split=1, precision=1e-6, verbose=0):
    if abs(gamma - 1.0) < 1e-12:
        gamma = 1.0
    if gamma < 1.0:
        print('Gamma must be >= 1.0')
        exit(0)
    n = np.load(os.path.join(os.path.dirname(__file__), "momenta.npy"), allow_pickle=True)
    res = A(q2, gamma, l, m, d, precision, verbose, m_split, n) + \
          B(q2, gamma, l, precision, verbose) + \
          C(q2, gamma, l, m, d, precision, verbose, m_split, n)
    if verbose:
        print('Luescher Zeta function:', res)
    return res

def appendSpherical_np(xyz):
    ptsnew = np.zeros(xyz.shape)
    xy = xyz[:,0]**2 + xyz[:,1]**2
    ptsnew[:,0] = np.sqrt(xy + xyz[:,2]**2)
    ptsnew[:,1] = np.arctan2(np.sqrt(xy), xyz[:,2]) 
    ptsnew[:,2] = np.arctan2(xyz[:,1], xyz[:,0])
    return ptsnew

def compute_r_in_spherical_coordinates(a, d, gamma, m_split):
    out = np.zeros(a.shape)
    d_norm = np.linalg.norm(d)
    if d_norm == 0.0:
        for i, r in enumerate(a):
            out[i,:] = r / gamma
    else:
        for i, r in enumerate(a):
            r_p = np.dot(r, d) / np.dot(d,d) * d
            r_o = r - r_p
            out[i,:] = (r_p - 0.5*m_split*d)/gamma + r_o
    return appendSpherical_np(out)

def sph_harm(m=0, l=0, phi=0, theta=0):
    if l == 0 and m == 0:
        return 0.28209479177387814
    elif l == 1 and m == -1:
        return 0.3454941494713355 * np.sin(theta) * np.exp(-1j*phi)
    elif l == 1 and m == 0:
        return 0.48860251190291992 * np.cos(theta)
    elif l == 1 and m == 1:
        return -0.3454941494713355 * np.sin(theta) * np.exp(1j*phi)
    elif l == 2 and m == -2:
        return 0.38627420202318957 * np.sin(theta)**2 * np.exp(-2j*phi)
    elif l == 2 and m == 0:
        return 0.31539156525252005 * (3*np.cos(theta)**2 - 1)
    elif l == 2 and m == 2:
        return 0.38627420202318957 * np.sin(theta)**2 * np.exp(2j*phi)
    else:
        return scipy.special.sph_harm_y(m, l, phi, theta)

def return_momentum_array(p, n):
    while p < len(n) and len(n[p,0]) == 0:
        p += 1
    out = n[p, 0]
    p += 1
    return out, p

def compute_summands_A(a_sph, q, l, m):
    inter = []
    counter = 0
    for r in a_sph:
        breaker = 0
        if counter > 0:
            for i in range(counter):
                if abs(float(r[0]) - float(a_sph[i,0])) < 1e-8:
                    inter.append(inter[i])
                    breaker = 1
                    break
        if breaker == 0:
            denom = r[0]**2 - q
            if abs(denom) < 1e-12:
                inter.append(0.0)
            else:
                inter.append(np.exp(-(r[0]**2 - q)) * r[0]**l / denom)
        counter += 1
    result = 0.0
    for i, r in enumerate(a_sph):
        result += inter[i] * sph_harm(m, l, r[2], r[1])
    return result

def A(q, gamma, l, m, d, precision, verbose, m_split, n):
    i = 0
    r, i = return_momentum_array(i, n)
    r_sph = compute_r_in_spherical_coordinates(r, d, gamma, m_split)
    result = compute_summands_A(r_sph, q, l, m)
    if verbose:
        print('convergence in term A:', i-1, result)
    eps = 1.0
    while eps > precision:
        r, i = return_momentum_array(i, n)
        r_sph = compute_r_in_spherical_coordinates(r, d, gamma, m_split)
        res_h = compute_summands_A(r_sph, q, l, m)
        if result != 0.0:
            eps = abs(res_h / result)
        result += res_h
        if verbose:
            print('\t', i-1, result, eps)
        if result == 0.0 and i > 4:
            break
    return result


def compute_gamma_w_in_spherical_coordinates(a, d, gamma):
    out = np.zeros(a.shape)
    d_norm = np.linalg.norm(d)
    if d_norm == 0.0:
        for i, r in enumerate(a):
            out[i,:] = r * gamma
    else:
        for i, r in enumerate(a):
            r_p = np.dot(r, d) / np.dot(d,d) * d
            r_o = r - r_p
            out[i,:] = r_p * gamma + r_o
    return appendSpherical_np(out)

def compute_summands_C(w_sph, w, q, gamma, l, m, d, m_split, precision):
    part1 = (-1j)**l * gamma * (np.abs(w_sph[:,0])**l) * \
            np.exp(-1j * m_split * math.pi * np.dot(w, d)) * \
            sph_harm(m, l, w_sph[:,2], w_sph[:,1])
    part2 = []
    counter = 0
    for ww in w_sph:
        breaker = 0
        if counter > 0:
            for i in range(counter):
                if abs(float(ww[0]) - float(w_sph[i,0])) < 1e-8:
                    part2.append(part2[i])
                    breaker = 1
                    break
        if breaker == 0:
            integrand = lambda t, qq, ll, ww2: ((math.pi/t)**(1.5 + ll)) * np.exp(qq*t - ww2/t)
            res, _ = scipy.integrate.quad(integrand, 0., 1.,
                                          args=(q, l, (math.pi*ww[0])**2),
                                          epsabs=precision*0.1, epsrel=precision*0.1, limit=1000)
            part2.append(res)
        counter += 1
    part2 = np.asarray(part2, dtype=float)
    return np.dot(part1, part2)

def C(q, gamma, l, m, d, precision, verbose, m_split, n):
    i = 1
    w, i = return_momentum_array(i, n)
    w_sph = compute_gamma_w_in_spherical_coordinates(w, d, gamma)
    result = compute_summands_C(w_sph, w, q, gamma, l, m, d, m_split, precision)
    if verbose:
        print('convergence in term C:', i-1, result)
    eps = 1.0
    while eps > precision:
        w, i = return_momentum_array(i, n)
        w_sph = compute_gamma_w_in_spherical_coordinates(w, d, gamma)
        res_h = compute_summands_C(w_sph, w, q, gamma, l, m, d, m_split, precision)
        if result != 0.0:
            eps = abs(res_h / result)
        result += res_h
        if verbose:
            print('\t', i-1, result, eps)
        if result == 0.0 and i > 4:
            break
    return result

def test():
    print('\nTest in cms:')
    Pcm = np.array([0.,0.,0.])
    q = 0.1207 * 24 / (2.*math.pi)
    gamma = 1.0
    zeta = Z(q*q, gamma, d=Pcm).real
    delta = np.arctan(math.pi**1.5 * q / zeta) * 180./math.pi
    if delta < 0: delta += 180
    print('delta:', delta, '(expected 136.6527)')

    print('\nTest in mv1:')
    Pcm = np.array([0.,0.,1.])
    L = 32
    q = 0.161 * L / (2.*math.pi)
    E = 0.440
    Ecm = 0.396
    gamma = E / Ecm
    Z00 = Z(q*q, gamma, d=Pcm).real
    Z20 = Z(q*q, gamma, d=Pcm, l=2).real
    delta = np.arctan(gamma * math.pi**1.5 * q / (Z00 + 2./(q*q*math.sqrt(5))*Z20)) * 180./math.pi
    if delta < 0: delta += 180
    print('delta:', delta, '(expected 115.7653)')

    print('\nTest in mv2:')
    Pcm = np.array([1.,1.,0.])
    L = 32
    q = 0.167 * L / (2.*math.pi)
    E = 0.490
    Ecm = 0.407
    gamma = E / Ecm
    Z00 = Z(q*q, gamma, d=Pcm).real
    Z20 = Z(q*q, gamma, d=Pcm, l=2).real
    Z22 = Z(q*q, gamma, d=Pcm, l=2, m=2).imag
    Z2_2 = Z(q*q, gamma, d=Pcm, l=2, m=-2).imag
    denominator = Z00 - (1./(q*q*math.sqrt(5)))*Z20 + (math.sqrt(3./10.)/(q*q))*(Z22 - Z2_2)
    delta = np.arctan(gamma * math.pi**1.5 * q / denominator) * 180./math.pi
    if delta < 0: delta += 180
    print('delta:', delta, '(expected 127.9930)')

if __name__ == "__main__":
    test()
