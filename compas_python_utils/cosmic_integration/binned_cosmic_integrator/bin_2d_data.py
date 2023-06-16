from .gpu_utils import xp


def bin_2d_data(data2d: xp, data1d: xp, bins: xp, axis=0):
    """Bin data2d and data1d along the given axis using the given 1d bins"""

    num_bins = len(bins)
    assert num_bins != data2d.shape[axis], "More bins than data-rows!"

    bin_data1d_falls_in = xp.digitize(data1d, bins)
    assert len(bin_data1d_falls_in) == len(data1d), "Something went wrong with binning"
    assert data2d.shape[axis] == len(data1d), "Data2d and bins do not match"

    # bin data
    binned_data_shape = list(data2d.shape)
    binned_data_shape[axis] = num_bins - 1
    binned_data = xp.zeros(binned_data_shape)
    for bii in range(1, num_bins):
        mask = bin_data1d_falls_in == bii
        masked_data = data2d.take(indices=xp.where(mask)[0], axis=axis)
        if axis == 0:
            binned_data[bii - 1, :] = xp.sum(masked_data, axis=axis)
        else:
            binned_data[:, bii - 1] = xp.sum(masked_data, axis=axis)
    return binned_data
