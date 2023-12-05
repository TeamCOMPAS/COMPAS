from .gpu_utils import xp


def bin_2d_data(data2d: xp, data1d: xp, bins: xp, axis=0):
    """Bin data2d and data1d along the given axis using the given 1d bins"""

    num_bins = len(bins)
    bin_data1d_falls_in = xp.digitize(data1d, bins, right=True)

    assert num_bins <= data2d.shape[axis], "More bins than data-rows!"
    assert len(bin_data1d_falls_in) == len(data1d), "Something went wrong with binning"
    assert data2d.shape[axis] == len(data1d), "Data2d and bins do not match"

    # bin data
    binned_data_shape = list(data2d.shape)
    binned_data_shape[axis] = num_bins
    binned_data = xp.zeros(binned_data_shape)

    for i in range(num_bins):
        if axis == 0:
             binned_data[i, :] = xp.sum(data2d[bin_data1d_falls_in == i, :], axis=axis)
        else:
            binned_data[:, i] = xp.sum(data2d[:, bin_data1d_falls_in == i], axis=axis)
    return binned_data
