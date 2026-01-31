/**
 * HDF5 API Compatibility Shim
 *
 * This file provides backward compatibility for CGNS libraries built against
 * older HDF5 versions that use H5Literate with H5L_iterate1_t callback.
 *
 * HDF5 1.12+ renamed H5Literate to H5Literate1, with H5Literate2 as default.
 * HDF5 1.14+ further tightens the macro mapping.
 *
 * We only provide the shim for HDF5 1.12.x and 1.13.x where H5Literate1
 * exists but the old H5Literate symbol may be needed by CGNS.
 * For HDF5 >= 1.14, the native API handles everything correctly.
 */

#include <hdf5.h>

/*
 * Only define for HDF5 1.12.x-1.13.x.
 * For < 1.12, H5Literate already exists natively.
 * For >= 1.14, the macro system handles it and redefining causes conflicts.
 */
#if H5_VERSION_GE(1,12,0) && !H5_VERSION_GE(1,14,0)

herr_t H5Literate_compat(hid_t grp_id, H5_index_t idx_type, H5_iter_order_t order,
                          hsize_t *idx, H5L_iterate1_t op, void *op_data)
{
    return H5Literate1(grp_id, idx_type, order, idx, op, op_data);
}

#endif
