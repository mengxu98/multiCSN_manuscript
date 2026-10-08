def decode_hdf5_strings(values):
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def read_sparse_group(handle, group):
    import scipy.sparse as sp

    source = handle[group]
    matrix = sp.csc_matrix(
        (source["data"][:], source["indices"][:], source["indptr"][:]),
        shape=tuple(int(value) for value in source["shape"][:]),
    )
    return matrix, decode_hdf5_strings(source["features"][:]), decode_hdf5_strings(source["barcodes"][:])
