import numpy as np
import scipy.signal
import matplotlib.pyplot as plt


def ssa_decomposition(x, L, fe, do_plot=False):
    """
    Decompose a signal using Singular Spectral Analysis

    :param x: input matrix
    :param L: window length
    :param fe: sampling frequency
    :param do_plot: flag to plot results
    :return: fc : components frequencies in Hz
             sig : variances of the components
             Y : the components of the SSA
    """

    N = len(x)
    # build up of the signal matrix
    X = np.zeros((N-L+1, L))
    for k in np.arange(0, N-L+1):
        X[k, :] = np.transpose(x[k:k+L])

    U, S, V = np.linalg.svd(X, full_matrices=False)
    sig = S
    Y = np.zeros((N, L))
    fc = np.zeros((L, 1))
    y = np.zeros(N)

    for k in np.arange(0, L):

        # computation s_i*u_i*v_i^T
        u_k = U[:, k].reshape(N-L+1, 1)
        v_k = V[k, :].reshape(1, L)
        Z = np.fliplr(S[k]*u_k*v_k)

        # averaging the anti-diagonals
        for p in np.arange(-N+L+1, L+1):
            y[L-p] = np.mean(np.diag(Z, p-1))

        # component k is stored
        Y[:, k] = y

        # Estimate of component frequency
        fr, p = scipy.signal.welch(y-np.mean(y), window='hamming',
                                   nperseg=round(N/4),
                                   noverlap=None,
                                   nfft=round(N/4), fs=fe)
        ind = np.where(p == max(p))
        fc[k] = fr[ind]
        if do_plot == 1:
            t = np.arange(1/fe, N/fe, 1/fe)
            plt.plot(t, x)
            plt.plot(t, y, 'r', 'LineWidth', 2)
            plt.title('blue: signal, red: component #' + str(k))
            plt.xlabel('seconds')
            print('frequency ' + str(fc(k)) + ' Hz')

    sig_2 = sig*sig
    if do_plot == 2:
        fcs = np.sort(fc.flatten())
        sort_sig2 = list(x for _, x in sorted(zip(fc.flatten(), sig_2)))
        plt.figure()
        plt.plot(fcs, sort_sig2, marker='o', linestyle='dashed', linewidth=2, markersize=5)
        plt.xlabel('component frequencies')
        plt.ylabel('component variances')

    return Y, fc, sig_2
