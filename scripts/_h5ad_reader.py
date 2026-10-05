"""
_h5ad_reader.py -- read the parts of an .h5ad file that the domain comparison needs, without h5py
===============================================================================================

The analysis machine used for v2.2 had the HDF5 C library (libhdf5) but no h5py, anndata or Scanpy,
and no package index. This helper talks to libhdf5 through ctypes and returns NumPy arrays. It
supports exactly what SPARTA's ``*.scored.h5ad`` files contain:

  * dense or CSR/CSC ``X`` (anndata >= 0.8 on-disk format: a group with data/indices/indptr and a
    ``shape`` attribute);
  * ``obsm/<key>`` dense arrays;
  * ``obs``/``var`` index (``_index`` attribute) and plain numeric / boolean / string columns;
  * scalar string attributes such as ``encoding-type``.

It is read-only and is not a general HDF5 reader. With h5py installed, ``read_h5ad`` uses h5py
instead and returns the same arrays, so results do not depend on which path ran.
"""
from __future__ import annotations

import ctypes as C
import ctypes.util
from pathlib import Path

import numpy as np

_LIB = None
hid_t = C.c_int64
herr_t = C.c_int
hsize_t = C.c_uint64
H5P_DEFAULT = 0
H5F_ACC_RDONLY = 0
H5S_ALL = 0
H5T_VARIABLE = C.c_size_t(-1).value
# H5T_class_t
H5T_INTEGER, H5T_FLOAT, H5T_STRING, H5T_ENUM = 0, 1, 3, 8
# H5O_type_t (from H5Oget_info / H5Iget_type)
H5I_GROUP, H5I_DATASET = 2, 5


def _lib():
    global _LIB
    if _LIB is not None:
        return _LIB
    cand = [ctypes.util.find_library("hdf5"), ctypes.util.find_library("hdf5_serial"),
            "libhdf5_serial.so.103", "libhdf5.so.103", "libhdf5.so.200", "libhdf5.so.310", "libhdf5.so"]
    for name in [c for c in cand if c]:
        try:
            lib = C.CDLL(name)
            break
        except OSError:
            continue
    else:
        raise OSError("libhdf5 not found")
    lib.H5open()
    sig = {
        "H5Fopen": (hid_t, [C.c_char_p, C.c_uint, hid_t]),
        "H5Fclose": (herr_t, [hid_t]),
        "H5Oopen": (hid_t, [hid_t, C.c_char_p, hid_t]),
        "H5Oclose": (herr_t, [hid_t]),
        "H5Iget_type": (C.c_int, [hid_t]),
        "H5Lexists": (C.c_int, [hid_t, C.c_char_p, hid_t]),
        "H5Gget_info": (herr_t, [hid_t, C.c_void_p]),
        "H5Lget_name_by_idx": (C.c_ssize_t, [hid_t, C.c_char_p, C.c_int, C.c_int, hsize_t, C.c_char_p,
                                             C.c_size_t, hid_t]),
        "H5Dopen2": (hid_t, [hid_t, C.c_char_p, hid_t]),
        "H5Dclose": (herr_t, [hid_t]),
        "H5Dget_space": (hid_t, [hid_t]),
        "H5Dget_type": (hid_t, [hid_t]),
        "H5Dread": (herr_t, [hid_t, hid_t, hid_t, hid_t, hid_t, C.c_void_p]),
        "H5Dvlen_reclaim": (herr_t, [hid_t, hid_t, hid_t, C.c_void_p]),
        "H5Sget_simple_extent_ndims": (C.c_int, [hid_t]),
        "H5Sget_simple_extent_dims": (C.c_int, [hid_t, C.POINTER(hsize_t), C.POINTER(hsize_t)]),
        "H5Sget_simple_extent_npoints": (C.c_ssize_t, [hid_t]),
        "H5Sclose": (herr_t, [hid_t]),
        "H5Tget_class": (C.c_int, [hid_t]),
        "H5Tget_size": (C.c_size_t, [hid_t]),
        "H5Tget_sign": (C.c_int, [hid_t]),
        "H5Tget_super": (hid_t, [hid_t]),
        "H5Tis_variable_str": (C.c_int, [hid_t]),
        "H5Tget_native_type": (hid_t, [hid_t, C.c_int]),
        "H5Tcopy": (hid_t, [hid_t]),
        "H5Tset_size": (herr_t, [hid_t, C.c_size_t]),
        "H5Tset_cset": (herr_t, [hid_t, C.c_int]),
        "H5Tclose": (herr_t, [hid_t]),
        "H5Aexists": (C.c_int, [hid_t, C.c_char_p]),
        "H5Aopen": (hid_t, [hid_t, C.c_char_p, hid_t]),
        "H5Aclose": (herr_t, [hid_t]),
        "H5Aget_type": (hid_t, [hid_t]),
        "H5Aget_space": (hid_t, [hid_t]),
        "H5Aread": (herr_t, [hid_t, hid_t, C.c_void_p]),
    }
    for name, (res, args) in sig.items():
        f = getattr(lib, name)
        f.restype, f.argtypes = res, args
    _LIB = lib
    return lib


class _GInfo(C.Structure):
    _fields_ = [("storage_type", C.c_int), ("nlinks", hsize_t), ("max_corder", C.c_int64),
                ("mounted", C.c_bool)]


def _chk(v, what):
    if v < 0:
        raise OSError(f"HDF5 call failed: {what}")
    return v


def _numpy_dtype(lib, tid):
    cls, size = lib.H5Tget_class(tid), lib.H5Tget_size(tid)
    if cls == H5T_FLOAT:
        return np.dtype({4: np.float32, 8: np.float64, 2: np.float16}[size])
    if cls == H5T_INTEGER:
        signed = lib.H5Tget_sign(tid) == 1
        return np.dtype(f"{'i' if signed else 'u'}{size}")
    if cls == H5T_ENUM:  # h5py stores booleans as an int8 enum
        sup = lib.H5Tget_super(tid)
        try:
            return _numpy_dtype(lib, sup)
        finally:
            lib.H5Tclose(sup)
    raise TypeError(f"unsupported HDF5 type class {cls}")


def _shape(lib, sid):
    nd = lib.H5Sget_simple_extent_ndims(sid)
    if nd <= 0:
        return ()
    dims = (hsize_t * nd)()
    lib.H5Sget_simple_extent_dims(sid, dims, None)
    return tuple(int(d) for d in dims)


def _read_strings(lib, obj_id, tid, sid, is_attr):
    n = max(int(lib.H5Sget_simple_extent_npoints(sid)), 1)
    shape = _shape(lib, sid)
    if lib.H5Tis_variable_str(tid) > 0:
        mem = lib.H5Tcopy(tid)
        buf = (C.c_char_p * n)()
        if is_attr:
            _chk(lib.H5Aread(obj_id, mem, buf), "H5Aread vlen str")
        else:
            _chk(lib.H5Dread(obj_id, mem, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf), "H5Dread vlen str")
        out = [b.decode("utf-8") if b is not None else "" for b in buf]
        lib.H5Dvlen_reclaim(mem, sid, H5P_DEFAULT, buf)
        lib.H5Tclose(mem)
    else:
        size = lib.H5Tget_size(tid)
        raw = C.create_string_buffer(size * n)
        if is_attr:
            _chk(lib.H5Aread(obj_id, tid, raw), "H5Aread str")
        else:
            _chk(lib.H5Dread(obj_id, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw), "H5Dread str")
        out = [raw.raw[i * size:(i + 1) * size].split(b"\0", 1)[0].decode("utf-8") for i in range(n)]
    return out[0] if shape == () else np.array(out, dtype=object).reshape(shape)


class H5File:
    """Minimal read-only view of an HDF5 file."""

    def __init__(self, path):
        self.lib = _lib()
        self.fid = _chk(self.lib.H5Fopen(str(path).encode(), H5F_ACC_RDONLY, H5P_DEFAULT), f"open {path}")

    def close(self):
        if self.fid is not None:
            self.lib.H5Fclose(self.fid)
            self.fid = None

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    # -- structure ----------------------------------------------------------------------------
    def exists(self, path):
        parts = [p for p in path.strip("/").split("/") if p]
        cur = ""
        for p in parts:
            cur = f"{cur}/{p}"
            if self.lib.H5Lexists(self.fid, cur.encode(), H5P_DEFAULT) <= 0:
                return False
        return True

    def kind(self, path):
        oid = _chk(self.lib.H5Oopen(self.fid, path.encode(), H5P_DEFAULT), f"open {path}")
        try:
            t = self.lib.H5Iget_type(oid)
        finally:
            self.lib.H5Oclose(oid)
        return {H5I_GROUP: "group", H5I_DATASET: "dataset"}.get(t, str(t))

    def keys(self, path="/"):
        oid = _chk(self.lib.H5Oopen(self.fid, path.encode(), H5P_DEFAULT), f"open {path}")
        try:
            info = _GInfo()
            _chk(self.lib.H5Gget_info(oid, C.byref(info)), "H5Gget_info")
            names = []
            for i in range(info.nlinks):
                n = self.lib.H5Lget_name_by_idx(oid, b".", 0, 0, i, None, 0, H5P_DEFAULT)
                buf = C.create_string_buffer(n + 1)
                self.lib.H5Lget_name_by_idx(oid, b".", 0, 0, i, buf, n + 1, H5P_DEFAULT)
                names.append(buf.value.decode())
            return names
        finally:
            self.lib.H5Oclose(oid)

    # -- data ---------------------------------------------------------------------------------
    def attr(self, path, name, default=None):
        lib = self.lib
        oid = _chk(lib.H5Oopen(self.fid, path.encode(), H5P_DEFAULT), f"open {path}")
        try:
            if lib.H5Aexists(oid, name.encode()) <= 0:
                return default
            aid = _chk(lib.H5Aopen(oid, name.encode(), H5P_DEFAULT), f"attr {name}")
            tid, sid = lib.H5Aget_type(aid), lib.H5Aget_space(aid)
            try:
                if lib.H5Tget_class(tid) == H5T_STRING:
                    return _read_strings(lib, aid, tid, sid, True)
                dt = _numpy_dtype(lib, tid)
                shape = _shape(lib, sid)
                arr = np.empty(shape if shape else (1,), dtype=dt)
                nat = lib.H5Tget_native_type(tid, 1)
                _chk(lib.H5Aread(aid, nat, arr.ctypes.data_as(C.c_void_p)), "H5Aread")
                lib.H5Tclose(nat)
                return arr if shape else arr[0]
            finally:
                lib.H5Tclose(tid)
                lib.H5Sclose(sid)
                lib.H5Aclose(aid)
        finally:
            lib.H5Oclose(oid)

    def dataset(self, path):
        lib = self.lib
        did = _chk(lib.H5Dopen2(self.fid, path.encode(), H5P_DEFAULT), f"dataset {path}")
        tid, sid = lib.H5Dget_type(did), lib.H5Dget_space(did)
        try:
            if lib.H5Tget_class(tid) == H5T_STRING:
                return _read_strings(lib, did, tid, sid, False)
            dt = _numpy_dtype(lib, tid)
            shape = _shape(lib, sid)
            arr = np.empty(shape if shape else (1,), dtype=dt)
            nat = lib.H5Tget_native_type(tid, 1)
            _chk(lib.H5Dread(did, nat, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr.ctypes.data_as(C.c_void_p)), "H5Dread")
            lib.H5Tclose(nat)
            return arr if shape else arr[0]
        finally:
            lib.H5Tclose(tid)
            lib.H5Sclose(sid)
            lib.H5Dclose(did)


def _read_matrix(f: H5File, path: str):
    """Dense dataset or anndata sparse group -> numpy array or scipy sparse matrix."""
    if f.kind(path) == "dataset":
        return f.dataset(path)
    import scipy.sparse as sp
    enc = f.attr(path, "encoding-type", "")
    shape = tuple(int(x) for x in f.attr(path, "shape"))
    data, indices, indptr = (f.dataset(f"{path}/{k}") for k in ("data", "indices", "indptr"))
    if enc == "csc_matrix":
        return sp.csc_matrix((data, indices, indptr), shape=shape)
    return sp.csr_matrix((data, indices, indptr), shape=shape)


def _read_column(f: H5File, path: str):
    """Dataset, or an anndata >= 0.12 nullable array group (values + mask; masked entries -> None/NaN)."""
    if f.kind(path) == "dataset":
        return f.dataset(path)
    enc = f.attr(path, "encoding-type", "")
    if enc.startswith("nullable"):
        vals = f.dataset(f"{path}/values")
        if f.exists(f"{path}/mask"):
            mask = f.dataset(f"{path}/mask").astype(bool)
            if mask.any():
                vals = vals.astype(object if vals.dtype == object else float)
                vals[mask] = None if vals.dtype == object else np.nan
        return vals
    if enc == "categorical":
        cats, codes = f.dataset(f"{path}/categories"), f.dataset(f"{path}/codes")
        out = np.asarray(cats, dtype=object)[np.clip(codes, 0, None)]
        out[codes < 0] = None
        return out
    raise TypeError(f"unsupported column encoding {enc!r} at {path}")


def _read_index(f: H5File, grp: str):
    key = f.attr(grp, "_index", "_index")
    return _read_column(f, f"{grp}/{key}")


def read_h5ad(path, *, obsm_keys=None, var_columns=None, obs_columns=None):
    """Return dict(X, obs_names, var_names, obsm={...}, var={...}, obs={...}, uns_log1p_base)."""
    try:  # prefer h5py when it exists
        import h5py  # noqa: F401
        return _read_h5ad_h5py(path, obsm_keys, var_columns, obs_columns)
    except ImportError:
        pass
    with H5File(Path(path)) as f:
        out = dict(X=_read_matrix(f, "/X"), obs_names=_read_index(f, "/obs"), var_names=_read_index(f, "/var"),
                   obsm={}, var={}, obs={})
        if f.exists("/obsm"):
            for k in f.keys("/obsm"):
                if obsm_keys is None or k in obsm_keys:
                    out["obsm"][k] = _read_matrix(f, f"/obsm/{k}")
        for grp, cols, dst in (("/var", var_columns, out["var"]), ("/obs", obs_columns, out["obs"])):
            for c in (cols or []):
                if f.exists(f"{grp}/{c}"):
                    dst[c] = _read_column(f, f"{grp}/{c}")
        base = None
        if f.exists("/uns/log1p/base") and f.attr("/uns/log1p/base", "encoding-type", "") != "null":
            base = f.dataset("/uns/log1p/base")
        out["uns_log1p_base"] = base
        out["uns_log1p"] = f.exists("/uns/log1p")
        return out


def _read_h5ad_h5py(path, obsm_keys, var_columns, obs_columns):
    import h5py
    import scipy.sparse as sp

    def mat(g):
        if isinstance(g, h5py.Dataset):
            return g[()]
        shape = tuple(int(x) for x in g.attrs["shape"])
        cls = sp.csc_matrix if g.attrs.get("encoding-type", "") == "csc_matrix" else sp.csr_matrix
        return cls((g["data"][()], g["indices"][()], g["indptr"][()]), shape=shape)

    def idx(g):
        k = g.attrs.get("_index", "_index")
        d = g[k] if isinstance(g[k], h5py.Dataset) else g[k]["values"]
        return np.array([x.decode() if isinstance(x, bytes) else x for x in d[()]], dtype=object)

    with h5py.File(path, "r") as f:
        out = dict(X=mat(f["X"]), obs_names=idx(f["obs"]), var_names=idx(f["var"]), obsm={}, var={}, obs={})
        for k in f.get("obsm", {}):
            if obsm_keys is None or k in obsm_keys:
                out["obsm"][k] = mat(f["obsm"][k])
        for grp, cols, dst in (("var", var_columns, out["var"]), ("obs", obs_columns, out["obs"])):
            for c in (cols or []):
                if c in f[grp] and isinstance(f[grp][c], h5py.Dataset):
                    dst[c] = f[grp][c][()]
        b = f.get("uns/log1p/base")
        out["uns_log1p_base"] = None if b is None or b.attrs.get("encoding-type", "") == "null" else b[()]
        out["uns_log1p"] = "uns/log1p" in f
        return out
