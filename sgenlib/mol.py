# Author: Samuel Genheden, samuel.genheden@gmail.com

import numpy as np

def density_scaling(xvals, area) :
    dx = (xvals[1] - xvals[0])
    lenz = xvals[-1]+dx-xvals[0]
    return len(xvals) / (area * lenz)

def density_intercept(dens1, dens2) :

    n1 = np.sum(dens1)
    n2 = np.sum(dens2)

    fi = 0
    while  dens1[fi] / n1 == 0.0 or dens2[fi] / n2 < dens1[fi] / n1 :
        fi += 1

    li = len(dens1) - 1
    while  dens1[li] / n1 == 0.0 or dens2[li] / n2 < dens1[li] / n1 :
        li -= 1

    return fi, li

def _make_samples(density, xvals) :

    density2 = density * len(xvals)
    xvals2 = list(xvals)
    samples = []
    for x, d in zip(xvals2, density2) :
        n = int(d)
        samples.extend([x]*n)
    xvals2.append(xvals2[-1]+(xvals2[1] - xvals2[0]))
    return np.asarray(samples), xvals2

def bootstrap_density(density, xvals, nboots=1000, calc_max=False) :

    samples, bins = _make_samples(density, xvals)
    hist_scale = 1 / float(len(xvals))

    dens_boots = np.zeros([density.shape[0],nboots])
    if calc_max :
        max_boots = np.zeros(nboots)
    for i in range(nboots) :
        sel = np.random.randint(0 ,samples.shape[0], samples.shape[0])
        dens_boots[:, i], bins = np.histogram(samples[sel], bins=bins)
        if calc_max :
            max_boots[i] = bins[np.argmax(dens_boots[:,i])]

    dens_std = dens_boots.std(axis=1) * hist_scale
    if calc_max :
        return dens_std, max_boots.std()
    else :
        return dens_std


def bootstrap_density_intercept(dens1, dens2, xvals, nboots=100, user_func=None, **kwargs) :

    samples1, bins = _make_samples(dens1, xvals)
    new_dens1, dummy = np.histogram(samples1, bins=bins)
    samples2, bins = _make_samples(dens2, xvals)
    new_dens2, dummy = np.histogram(samples2, bins=bins)
    hist_scale = 1 / float(len(xvals))

    intercepts = np.zeros([2, nboots])
    if user_func is not None :
        user_boots = np.zeros(nboots)
    for i in range(nboots) :
        sel1 = np.random.randint(0 ,samples1.shape[0], samples1.shape[0])
        sel2 = np.random.randint(0 ,samples2.shape[0], samples2.shape[0])
        new_dens1, dummy = np.histogram(samples1[sel1], bins=bins)
        new_dens2, dummy = np.histogram(samples2[sel2], bins=bins)
        fi, li = density_intercept(new_dens1*hist_scale, new_dens2*hist_scale)
        intercepts[:, i] = xvals[fi], xvals[li]
        if user_func is not None :
            user_boots[i] = user_func(new_dens1*hist_scale, new_dens2*hist_scale,
                                        fi, li, xvals[fi], xvals[li], **kwargs)

    if user_func is not None :
        return intercepts.std(axis=1), user_boots.std()
    else :
        return intercepts.std(axis=1)
