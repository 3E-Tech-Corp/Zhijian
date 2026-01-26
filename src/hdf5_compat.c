/**
 * HDF5 API Compatibility Shim
 *
 * This file provides backward compatibility for CGNS libraries built against
 * older HDF5 versions (< 1.12) that use H5Literate, which was renamed to
 * H5Literate1 in HDF5 1.12+ and H5Literate2 became the new default.
 *
 * This shim forwards H5Literate calls to H5Literate1.
 */

#include <hdf5.h>

/* Only define if H5Literate is not already available (HDF5 1.12+) */
#if H5_VERSION_GE(1,12,0)

herr_t H5Literate(hid_t grp_id, H5_index_t idx_type, H5_iter_order_t order,
                  hsize_t *idx, H5L_iterate1_t op, void *op_data)
{
    return H5Literate1(grp_id, idx_type, order, idx, op, op_data);
}

#endif
