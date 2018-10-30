import numpy as np
from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity
import pickle
import pandas as pd

def generate(pdf, x_grid, nevents=10000):
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    values = np.random.rand(int(nevents))
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = x_grid[value_bins]
    return random_from_cdf

def kde(x, x_grid, bandwidth=0.2, **kwargs):
#    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
#    return kde.evaluate(x_grid)
    kde = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde.fit(x[:, np.newaxis])
    log_pdf = kde.score_samples(x_grid[:, np.newaxis])
    return np.exp(log_pdf)

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    data = pd.read_csv('../flux/t2k_flux.csv',index_col=0)
    x_grid = data['low'].values
    integral = np.sum(data['numu_flux'].values)
    random_from_cdf = generate(data['numu_flux'].values/integral, x_grid,1E7)

    try:
        kdepdf = pickle.load(open('kde.pickle','rb'))
    except (OSError, IOError) as e:
        kdepdf = kde(random_from_cdf,x_grid,bandwidth=0.1)
        pickle.dump(kdepdf, open('kde.pickle','wb'))
    random_from_kde = generate(kdepdf, x_grid, 1E7)

    plt.subplot(121)
    plt.plot(data['low'].values,data['numu_flux'].values/integral/0.05,ls='steps', label='hist')
    plt.plot(x_grid, kdepdf, color='r')
    plt.yscale('log')
    plt.ylim(1E-8,5)
    plt.subplot(122)
    plt.hist(random_from_cdf, len(x_grid))
    plt.hist(random_from_kde, len(x_grid))
    plt.yscale('log')
    plt.show()
