import numpy as np
import scipy.sparse as sp
import scipy.linalg as la


class BdGOperator:
    def __init__(
        self,
        con_func_a,
        con_func_b,
        mu_a,
        mu_b,
        vort_a,
        vort_b,
        azi_number_pert,
        frac_a,
        frac_b,
        ga,
        gb,
        gab,
    ):
        self.dtype = np.complex128
        self.tht_pts = con_func_a.size
        self.dtht = pi / (self.tht_pts - 1)
        self.tht = np.linspace(0, pi, self.tht_pts)
        self.shape = (4 * (self.tht_pts - 2), 4 * (self.tht_pts - 2))
        self.vort_a = vort_a
        self.vort_b = vort_b
        self.azi_number_pert = azi_number_pert
        self.mu_a = mu_a
        self.mu_b = mu_b
        self.nabla_factor = -0.5
        self.inter = np.array([ga, gb, gab])
        self.frac_a = frac_a
        self.frac_b = frab_b
        self.con_func_a = con_func_a
        self.con_func_b = con_func_b

    def sparse_matrix(self):
        dtht = self.dtht
        nabl = self.nabla_factor
        sin_th = np.sin(self.tht)
        cos_th = np.cos(self.tht)
        m = self.azi_number_pert
        sa, sb = self.vort_a, self.vort_b
        ga, gb, gab = self.inter[0], self.inter[1], self.inter[2]
        fa = self.con_func_a
        fb = self.con_func_b
        sign = -1
        block_stride = self.tht_pts - 2
        val = np.empty(4 * (6 * block_stride - 2), dtype=np.complex128)
        col = np.empty(val.size, dtype=np.int32)
        row = np.empty(4 * block_stride + 1, dtype=np.int32)
        row[0] = 0
        k = 0
        sqrt_ratio = sqrt(self.frac_b / self.frac_a)
        # block 1 ============================================================
        for i in range(1, self.tht_pts - 1):
            # Block diagonal setup -------------------------------------------
            j = i - 1
            if j > 0:
                val[k] = nabl * (
                    1.0 / dtht / dtht - 0.5 * cos_th[i] / sin_th[i] / dtht
                )
                col[k] = j - 1
                k = k + 1
            val[k] = (
                nabl * (-2.0 / dtht / dtht - (sa + m) ** 2 / sin_th[i] ** 2)
                + 2 * ga * abs(fa[i]) ** 2
                + sqrt_ratio * gab * abs(fb[i]) ** 2
                - mu_a
            )
            col[k] = j
            k = k + 1
            if j < block_stride - 1:
                val[k] = nabl * (
                    1.0 / dtht / dtht + 0.5 * cos_th[i] / sin_th[i] / dtht
                )
                col[k] = j + 1
                k = k + 1
            # finish block diagonal setup ------------------------------------
            val[k] = ga * fa[i] ** 2
            col[k] = j + block_stride
            k = k + 1
            val[k] = sqrt_ratio * gab * fb[i].conj() * fa[i]
            col[k] = j + 2 * block_stride
            k = k + 1
            val[k] = sqrt_ratio * gab * fa[i] * fb[i]
            col[k] = j + 3 * block_stride
            k = k + 1
            row[i] = k
        # block 2 ============================================================
        for i in range(1, self.tht_pts - 1):
            j = i - 1
            val[k] = sign * (ga * fa[i] ** 2).conj()
            col[k] = j
            k = k + 1
            # block diagonal setup -------------------------------------------
            if j > 0:
                val[k] = sign * (
                    nabl
                    * (1.0 / dtht / dtht - 0.5 * cos_th[i] / sin_th[i] / dtht)
                )
                col[k] = j - 1 + block_stride
                k = k + 1
            val[k] = sign * (
                nabl * (-2.0 / dtht / dtht - (sa - m) ** 2 / sin_th[i] ** 2)
                + 2 * ga * abs(fa[i]) ** 2
                + sqrt_ratio * gab * abs(fb[i]) ** 2
                - mu_a
            )
            col[k] = j + block_stride
            k = k + 1
            if j < block_stride - 1:
                val[k] = sign * (
                    nabl
                    * (1.0 / dtht / dtht + 0.5 * cos_th[i] / sin_th[i] / dtht)
                )
                col[k] = j + 1 + block_stride
                k = k + 1
            # finish block diagonal setup ------------------------------------
            val[k] = sign * sqrt_ratio * gab * (fb[i] * fa[i]).conj()
            col[k] = j + 2 * block_stride
            k = k + 1
            val[k] = sign * sqrt_ratio * gab * fa[i].conj() * fb[i]
            col[k] = j + 3 * block_stride
            k = k + 1
            row[i + block_stride] = k
        # block 3 ============================================================
        sqrt_ratio = sqrt(self.frac_a / self.frac_b)
        for i in range(1, self.tht_pts - 1):
            j = i - 1
            val[k] = sqrt_ratio * gab * fb[i] * fa[i].conj()
            col[k] = j
            k = k + 1
            val[k] = sqrt_ratio * gab * fa[i] * fb[i]
            col[k] = j + block_stride
            k = k + 1
            # block diagonal setup -------------------------------------------
            if j > 0:
                val[k] = nabl * (
                    1.0 / dtht / dtht - 0.5 * cos_th[i] / sin_th[i] / dtht
                )
                col[k] = j - 1 + 2 * block_stride
                k = k + 1
            val[k] = (
                nabl * (-2.0 / dtht / dtht - (sb + m) ** 2 / sin_th[i] ** 2)
                + 2 * gb * abs(fb[i]) ** 2
                + sqrt_ratio * gab * abs(fa[i]) ** 2
                - mu_b
            )
            col[k] = j + 2 * block_stride
            k = k + 1
            if j < block_stride - 1:
                val[k] = nabl * (
                    1.0 / dtht / dtht + 0.5 * cos_th[i] / sin_th[i] / dtht
                )
                col[k] = j + 1 + 2 * block_stride
                k = k + 1
            # finish block diagonal setup ------------------------------------
            val[k] = gb * fb[i] ** 2
            col[k] = j + 3 * block_stride
            k = k + 1
            row[i + 2 * block_stride] = k
        # block 4 ============================================================
        for i in range(1, self.tht_pts - 1):
            j = i - 1
            val[k] = sign * sqrt_ratio * gab * (fb[i] * fa[i]).conj()
            col[k] = j
            k = k + 1
            val[k] = sign * sqrt_ratio * gab * fa[i] * fb[i].conj()
            col[k] = j + block_stride
            k = k + 1
            val[k] = sign * (gb * fb[i] ** 2).conj()
            col[k] = j + 2 * block_stride
            k = k + 1
            # block diagonal part -------------------------------------------
            if j > 0:
                val[k] = sign * (
                    nabl
                    * (1.0 / dtht / dtht - 0.5 * cos_th[i] / sin_th[i] / dtht)
                )
                col[k] = j - 1 + 3 * block_stride
                k = k + 1
            val[k] = sign * (
                nabl * (-2.0 / dtht / dtht - (sb - m) ** 2 / sin_th[i] ** 2)
                + 2 * gb * abs(fb[i]) ** 2
                + sqrt_ratio * gab * abs(fa[i]) ** 2
                - mu_b
            )
            col[k] = j + 3 * block_stride
            k = k + 1
            if j < block_stride - 1:
                val[k] = sign * (
                    nabl
                    * (1.0 / dtht / dtht + 0.5 * cos_th[i] / sin_th[i] / dtht)
                )
                col[k] = j + 1 + 3 * block_stride
                k = k + 1
            # finish block diagonal setup ------------------------------------
            row[i + 3 * block_stride] = k
        return sp.csr_matrix(
            (val, col, row),
            shape=(row.size - 1, row.size - 1),
            dtype=np.complex128,
        )
