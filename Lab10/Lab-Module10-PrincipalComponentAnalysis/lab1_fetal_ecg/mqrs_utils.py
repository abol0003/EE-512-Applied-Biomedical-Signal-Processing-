import numpy as np
from scipy.linalg import svd
from scipy.interpolate import interp1d


def cancel_mqrs(fs, x, mqrs, n_eigenvectors_to_keep=3):
    """Cancel maternal QRS peaks as described in section 2.6 of
    M Varanini et al 2014 Physiol. Meas. 35 1607.

    Args:
        fs: sampling frequency (in Hz)
        x: multi-channel ECG signals of size nSamples x nChannels.
        mqrs: indices to previously detected mQRS locations.        
        n_eigenvectors_to_keep: number of eigenvectors to keep for QRS template.
        
    Returns:
        residual signal of `x` with cancelled mQRS peaks.
        mQRS signals based on templates
    """

    # number of samples in window around mQRS
    qrs_window = (-np.fix(0.2 * fs), np.fix(0.5 * fs))
    # prepare indices for matrix containing all mQRS complexes
    from scipy.signal.windows import tukey
    tmp = tukey(M=101, alpha=0.5, sym=True)
    f = interp1d([qrs_window[0], 0, qrs_window[-1]], [0, (len(tmp)-1)/2, (len(tmp)-1)])
    idxs = f(np.arange(qrs_window[0], qrs_window[-1]))
    weighting_win = interp1d(np.arange(len(tmp)), tmp)(idxs) * 0.8 + 0.2
    weighting_win = weighting_win.reshape(-1, 1).repeat(len(mqrs), axis=1)
    lookup_matrix = (np.arange(*qrs_window).reshape(-1, 1).repeat(
        len(mqrs), axis=1) + mqrs.reshape(1, -1).repeat(
        int(np.diff(qrs_window)), axis=0))
    lookup_matrix = lookup_matrix.astype(int)

    # perform mQRS cancelling for each channel/lead individually
    x_mqrs = np.full(x.shape, np.nan)

    for ic in range(x.shape[1]):
        # create matrix with all mQRS complexes for this channel/lead
        #   matrix of size nSamplesInWindow x nMaternalQRSComplexes
        A = x[lookup_matrix.astype(int), ic].squeeze() * weighting_win

        # perform SVD on all mQRS complexes
        u, s, vh = svd(A, full_matrices=False)
        s_reduced = np.zeros(s.shape)
        s_reduced[:n_eigenvectors_to_keep] = s[:n_eigenvectors_to_keep]
        ar = u @ np.diag(s_reduced) @ vh

        # approximate each individual mQRS using SVD
        unweighted_mqrs = ar / weighting_win
        for iqrs in range(A.shape[1]):
            x_mqrs[lookup_matrix[:, iqrs], ic] = unweighted_mqrs[:, iqrs]

        # create smooth connections between successive mQRS
        to_interpolate = np.isnan(x_mqrs[:, ic])
        idxs = np.arange(len(to_interpolate))
        f = interp1d(idxs[~to_interpolate], x_mqrs[~to_interpolate, ic],
                     kind='linear', fill_value=0, bounds_error=False)
        x_mqrs[to_interpolate, ic] = f(idxs[to_interpolate])

    # remove mQRS to get fQRS in residual signal
    x_residual = x - x_mqrs

    return x_residual, x

