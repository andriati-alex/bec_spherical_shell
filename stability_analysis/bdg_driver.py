from math import pi, sqrt
import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
from scipy.linalg import eig
from scipy.integrate import simps


class BdGOperator:
    def __init__(
        self,
        con_func_a,
        con_func_b,
        vort_a,
        vort_b,
        frac_a,
        frac_b,
    ):
        self.tht_pts = con_func_a.size
        self.dtht = pi / (self.tht_pts - 1)
        self.tht = np.linspace(0, pi, self.tht_pts)
        self.shape = (4 * (self.tht_pts - 2), 4 * (self.tht_pts - 2))
        self.vort_a = vort_a
        self.vort_b = vort_b
        self.nabla_factor = -0.5
        self.frac_a = frac_a
        self.frac_b = frac_b
        self.con_func_a = np.sqrt(abs(con_func_a) ** 2)
        self.con_func_b = np.sqrt(abs(con_func_b) ** 2)

    def __chem(self, vort, f, f_ext, frac, frac_ext, g_self, g_inter):
        npts = self.tht_pts
        nabl = self.nabla_factor
        dtht = self.dtht
        sin_th = np.sin(self.tht)
        cos_th = np.cos(self.tht)
        sqrt_ratio = sqrt(frac_ext / frac)
        twice_der = np.empty(npts)
        single_der = np.empty(npts)

        twice_der[1 : npts - 1] = (
            f[2:npts] - 2 * f[1 : npts - 1] + f[: npts - 2]
        ) / dtht ** 2
        twice_der[0] = (
            f[1] - 2 * f[0] + (1 - 2 * (abs(vort) % 2)) * f[1]
        ) / dtht ** 2
        twice_der[npts - 1] = (
            (1 - 2 * (abs(vort) % 2)) * f[npts - 2]
            - 2 * f[npts - 1]
            + f[npts - 2]
        ) / dtht ** 2
        single_der[1 : npts - 1] = (f[2:npts] - f[: npts - 2]) / (2 * dtht)
        single_der[0] = (f[1] - (1 - 2 * (abs(vort) % 2)) * f[1]) / (2 * dtht)
        single_der[npts - 1] = (
            (1 - 2 * (abs(vort) % 2)) * f[npts - 2] - f[npts - 2]
        ) / (2 * dtht)
        phi_der = np.zeros(npts)
        phi_der[1 : npts - 1] = (
            -(vort ** 2) * f[1 : npts - 1] / sin_th[1 : npts - 1] ** 2
        )
        integ = (
            nabl
            * f
            * (sin_th * twice_der + cos_th * single_der + sin_th * phi_der)
            + sin_th * g_self * f ** 4
            + sin_th * sqrt_ratio * g_inter * f ** 2 * f_ext ** 2
        )
        return simps(integ, dx=dtht)

    def chem_a(self, ga, gab, fa=None, fb=None):
        fa = fa or self.con_func_a
        fb = fb or self.con_func_b
        return self.__chem(
            self.vort_a,
            fa,
            fb,
            self.frac_a,
            self.frac_b,
            ga,
            gab,
        )

    def chem_b(self, gb, gab, fa=None, fb=None):
        fa = fa or self.con_func_a
        fb = fb or self.con_func_b
        return self.__chem(
            self.vort_b,
            fb,
            fa,
            self.frac_b,
            self.frac_a,
            gb,
            gab,
        )

    def lowlying_eig(
        self,
        m,
        mu_a,
        mu_b,
        ga,
        gb,
        gab,
        fa=None,
        fb=None,
        neigs=80,
        sigma=0.67 + 0.33j,
    ):
        mat = self.sparse_diag2(m, mu_a, mu_b, ga, gb, gab, fa, fb)
        eigvals, eigvecs = sp.linalg.eigs(mat, neigs, which="LM", sigma=sigma)
        return np.sort(eigvals)[::-1]

    def all_eig(
        self,
        m,
        mu_a,
        mu_b,
        ga,
        gb,
        gab,
        fa=None,
        fb=None,
    ):
        mat = self.sparse_diag2(m, mu_a, mu_b, ga, gb, gab, fa, fb)
        eigvals, eigvecs = eig(mat.toarray())
        return np.sort(eigvals)[::-1]

    def sparse_diag2(self, m, mu_a, mu_b, ga, gb, gab, fa=None, fb=None):
        if isinstance(fa, type(None)):
            fa = self.con_func_a
        if isinstance(fb, type(None)):
            fb = self.con_func_b
        dtht = self.dtht
        npts = self.tht_pts
        nabl = self.nabla_factor
        sin_th = np.sin(self.tht)
        cos_th = np.cos(self.tht)
        sa, sb = self.vort_a, self.vort_b
        fa = self.con_func_a
        fb = self.con_func_b
        sqrt_ratio_a = sqrt(self.frac_b / self.frac_a)
        sqrt_ratio_b = sqrt(self.frac_a / self.frac_b)
        sign = -1

        diag_main1 = np.zeros(npts)
        diag_main1[1 : npts - 1] = (
            nabl
            * (-2.0 / dtht ** 2 - (sa + m) ** 2 / sin_th[1 : npts - 1] ** 2)
            + 2 * ga * fa[1 : npts - 1] ** 2
            + sqrt_ratio_a * gab * fb[1 : npts - 1] ** 2
            - mu_a
        )
        if sa + m == 0:
            diag_main1[0] = (
                nabl * (-2.0 / dtht / dtht)
                + 2 * ga * fa[0] ** 2
                + sqrt_ratio_a * gab * fb[0] ** 2
                - mu_a
            )
            diag_main1[npts - 1] = (
                nabl * (-2.0 / dtht / dtht)
                + 2 * ga * fa[npts - 1] ** 2
                + sqrt_ratio_a * gab * fb[npts - 1] ** 2
                - mu_a
            )

        diag_main2 = np.zeros(npts)
        diag_main2[1 : npts - 1] = sign * (
            nabl
            * (-2.0 / dtht / dtht - (sa - m) ** 2 / sin_th[1 : npts - 1] ** 2)
            + 2 * ga * fa[1 : npts - 1] ** 2
            + sqrt_ratio_a * gab * fb[1 : npts - 1] ** 2
            - mu_a
        )
        if sa - m == 0:
            diag_main2[0] = sign * (
                nabl * (-2.0 / dtht / dtht)
                + 2 * ga * fa[0] ** 2
                + sqrt_ratio_a * gab * fb[0] ** 2
                - mu_a
            )
            diag_main2[npts - 1] = sign * (
                nabl * (-2.0 / dtht / dtht)
                + 2 * ga * fa[npts - 1] ** 2
                + sqrt_ratio_a * gab * fb[npts - 1] ** 2
                - mu_a
            )

        diag_main3 = np.zeros(npts)
        diag_main3[1 : npts - 1] = (
            nabl
            * (-2.0 / dtht / dtht - (sb + m) ** 2 / sin_th[1 : npts - 1] ** 2)
            + 2 * gb * fb[1 : npts - 1] ** 2
            + sqrt_ratio_b * gab * fa[1 : npts - 1] ** 2
            - mu_b
        )
        if sb + m == 0:
            diag_main3[0] = (
                nabl * (-2.0 / dtht / dtht)
                + 2 * gb * fb[0] ** 2
                + sqrt_ratio_b * gab * fa[0] ** 2
                - mu_b
            )
            diag_main3[npts - 1] = (
                nabl * (-2.0 / dtht / dtht)
                + 2 * gb * fb[npts - 1] ** 2
                + sqrt_ratio_b * gab * fa[npts - 1] ** 2
                - mu_b
            )

        diag_main4 = np.zeros(npts)
        diag_main4[1 : npts - 1] = sign * (
            nabl
            * (-2.0 / dtht / dtht - (sb - m) ** 2 / sin_th[1 : npts - 1] ** 2)
            + 2 * gb * fb[1 : npts - 1] ** 2
            + sqrt_ratio_b * gab * fa[1 : npts - 1] ** 2
            - mu_b
        )
        if sb - m == 0:
            diag_main4[0] = sign * (
                nabl * (-2.0 / dtht / dtht)
                + 2 * gb * fb[0] ** 2
                + sqrt_ratio_b * gab * fa[0] ** 2
                - mu_b
            )
            diag_main4[npts - 1] = sign * (
                nabl * (-2.0 / dtht / dtht)
                + 2 * gb * fb[npts - 1] ** 2
                + sqrt_ratio_b * gab * fa[npts - 1] ** 2
                - mu_b
            )

        diag_main = np.concatenate(
            [diag_main1, diag_main2, diag_main3, diag_main4]
        )
        # =====================================================================

        diag_upp1 = np.zeros(npts)
        diag_upp1[1 : npts - 2] = nabl * (
            1.0 / dtht / dtht
            + 0.5 * cos_th[1 : npts - 2] / sin_th[1 : npts - 2] / dtht
        )
        if sa + m == 0:
            diag_upp1[0] = 2 * nabl * (1.0 / dtht / dtht)
            diag_upp1[npts - 2] = nabl * (
                1.0 / dtht / dtht
                + 0.5 * cos_th[npts - 2] / sin_th[npts - 2] / dtht
            )

        diag_upp2 = np.zeros(npts)
        diag_upp2[1 : npts - 2] = sign * (
            nabl
            * (
                1.0 / dtht / dtht
                + 0.5 * cos_th[1 : npts - 2] / sin_th[1 : npts - 2] / dtht
            )
        )
        if sa - m == 0:
            diag_upp2[0] = 2 * sign * nabl * (1.0 / dtht / dtht)
            diag_upp2[npts - 2] = sign * (
                nabl
                * (
                    1.0 / dtht / dtht
                    + 0.5 * cos_th[npts - 2] / sin_th[npts - 2] / dtht
                )
            )

        diag_upp3 = np.zeros(npts)
        diag_upp3[1 : npts - 2] = nabl * (
            1.0 / dtht / dtht
            + 0.5 * cos_th[1 : npts - 2] / sin_th[1 : npts - 2] / dtht
        )
        if sb + m == 0:
            diag_upp3[0] = 2 * nabl * (1.0 / dtht / dtht)
            diag_upp3[npts - 2] = nabl * (
                1.0 / dtht / dtht
                + 0.5 * cos_th[npts - 2] / sin_th[npts - 2] / dtht
            )

        diag_upp4 = np.zeros(npts - 1)
        diag_upp4[1 : npts - 2] = sign * (
            nabl
            * (
                1.0 / dtht / dtht
                + 0.5 * cos_th[1 : npts - 2] / sin_th[1 : npts - 2] / dtht
            )
        )
        if sb - m == 0:
            diag_upp4[0] = 2 * sign * nabl * (1.0 / dtht / dtht)
            diag_upp4[npts - 2] = sign * (
                nabl
                * (
                    1.0 / dtht / dtht
                    + 0.5 * cos_th[npts - 2] / sin_th[npts - 2] / dtht
                )
            )

        diag_upp = np.concatenate([diag_upp1, diag_upp2, diag_upp3, diag_upp4])
        # =====================================================================

        diag_low1 = np.zeros(npts)
        diag_low1[1 : npts - 2] = nabl * (
            1.0 / dtht / dtht
            - 0.5 * cos_th[2 : npts - 1] / sin_th[2 : npts - 1] / dtht
        )
        if sa + m == 0:
            diag_low1[0] = nabl * (
                1.0 / dtht / dtht - 0.5 * cos_th[1] / sin_th[1] / dtht
            )
            diag_low1[npts - 2] = 2 * nabl * (1.0 / dtht / dtht)

        diag_low2 = np.zeros(npts)
        diag_low2[1 : npts - 2] = sign * (
            nabl
            * (
                1.0 / dtht / dtht
                - 0.5 * cos_th[2 : npts - 1] / sin_th[2 : npts - 1] / dtht
            )
        )
        if sa - m == 0:
            diag_low2[0] = sign * (
                nabl * (1.0 / dtht / dtht - 0.5 * cos_th[1] / sin_th[1] / dtht)
            )
            diag_low2[npts - 2] = 2 * sign * nabl * (1.0 / dtht / dtht)

        diag_low3 = np.zeros(npts)
        diag_low3[1 : npts - 2] = nabl * (
            1.0 / dtht / dtht
            - 0.5 * cos_th[2 : npts - 1] / sin_th[2 : npts - 1] / dtht
        )
        if sb + m == 0:
            diag_low3[0] = nabl * (
                1.0 / dtht / dtht - 0.5 * cos_th[1] / sin_th[1] / dtht
            )
            diag_low3[npts - 2] = 2 * nabl * (1.0 / dtht / dtht)

        diag_low4 = np.zeros(npts - 1)
        diag_low4[1 : npts - 2] = sign * (
            nabl
            * (
                1.0 / dtht / dtht
                - 0.5 * cos_th[2 : npts - 1] / sin_th[2 : npts - 1] / dtht
            )
        )
        if sb - m == 0:
            diag_low4[0] = sign * (
                nabl * (1.0 / dtht / dtht - 0.5 * cos_th[1] / sin_th[1] / dtht)
            )
            diag_low4[npts - 2] = 2 * sign * nabl * (1.0 / dtht / dtht)

        diag_low = np.concatenate([diag_low1, diag_low2, diag_low3, diag_low4])

        # ---------------------------------------------------------------------
        d1 = np.zeros(npts)
        if sa - m != 0:
            d1[1 : npts - 1] = ga * fa[1 : npts - 1] ** 2
        else:
            d1 = ga * fa ** 2
        d2 = np.zeros(npts)
        if sb + m != 0:
            d2[1 : npts - 1] = (
                sign * gab * sqrt_ratio_a * fa[1 : npts - 1] * fb[1 : npts - 1]
            )
        else:
            d2 = sign * gab * sqrt_ratio_a * fa * fb
        d3 = np.zeros(npts)
        if sb - m != 0:
            d3[1 : npts - 1] = gb * fb[1 : npts - 1] ** 2
        else:
            d3 = gb * fb ** 2
        block_diag_upp1 = np.concatenate([d1, d2, d3])

        d1 = np.zeros(npts)
        if sb + m != 0:
            d1[1 : npts - 1] = (
                gab * sqrt_ratio_a * fa[1 : npts - 1] * fb[1 : npts - 1]
            )
        else:
            d1 = gab * sqrt_ratio_a * fa * fb
        d2 = np.zeros(npts)
        if sb - m != 0:
            d2[1 : npts - 1] = (
                sign * gab * sqrt_ratio_a * fa[1 : npts - 1] * fb[1 : npts - 1]
            )
        else:
            d2 = sign * gab * sqrt_ratio_a * fa * fb
        block_diag_upp2 = np.concatenate([d1, d2])

        block_diag_upp3 = np.zeros(npts)
        if sb - m != 0:
            block_diag_upp3[1 : npts - 1] = (
                sqrt_ratio_a * gab * fa[1 : npts - 1] * fb[1 : npts - 1]
            )
        else:
            block_diag_upp3 = sqrt_ratio_a * gab * fa * fb

        # ---------------------------------------------------------------------
        d1 = np.zeros(npts)
        if sa + m != 0:
            d1[1 : npts - 1] = sign * ga * fa[1 : npts - 1] ** 2
        else:
            d1 = sign * ga * fa ** 2
        d2 = np.zeros(npts)
        if sa - m != 0:
            d2[1 : npts - 1] = (
                gab * sqrt_ratio_b * fa[1 : npts - 1] * fb[1 : npts - 1]
            )
        else:
            d2 = gab * sqrt_ratio_b * fa * fb
        d3 = np.zeros(npts)
        if sb + m != 0:
            d3[1 : npts - 1] = sign * gb * fb[1 : npts - 1] ** 2
        else:
            d3 = sign * gb * fb ** 2
        block_diag_low1 = np.concatenate([d1, d2, d3])

        d1 = np.zeros(npts)
        if sa + m != 0:
            d1[1 : npts - 1] = (
                gab * sqrt_ratio_b * fa[1 : npts - 1] * fb[1 : npts - 1]
            )
        else:
            d1 = gab * sqrt_ratio_b * fa * fb
        d2 = np.zeros(npts)
        if sa - m != 0:
            d2[1 : npts - 1] = (
                sign * gab * sqrt_ratio_b * fa[1 : npts - 1] * fb[1 : npts - 1]
            )
        else:
            d2 = sign * gab * sqrt_ratio_b * fa * fb
        block_diag_low2 = np.concatenate([d1, d2])

        block_diag_low3 = np.zeros(npts)
        if sa + m == 0:
            block_diag_low3[1 : npts - 1] = sign * (
                sqrt_ratio_b * gab * fa[1 : npts - 1] * fb[1 : npts - 1]
            )
        else:
            block_diag_low3 = sign * sqrt_ratio_b * gab * fa * fb

        all_diags = [
            block_diag_low3,
            block_diag_low2,
            block_diag_low1,
            diag_low,
            diag_main,
            diag_upp,
            block_diag_upp1,
            block_diag_upp2,
            block_diag_upp3,
        ]
        offsets = [
            -3 * npts,
            -2 * npts,
            -npts,
            -1,
            0,
            1,
            npts,
            2 * npts,
            3 * npts,
        ]
        return sp.diags(
            all_diags,
            offsets,
            dtype=np.float64,
            shape=(4 * npts, 4 * npts),
        )

    def sparse_diag(self, m, mu_a, mu_b, ga, gb, gab):
        dtht = self.dtht
        npts = self.tht_pts
        nabl = self.nabla_factor
        sin_th = np.sin(self.tht)
        cos_th = np.cos(self.tht)
        sa, sb = self.vort_a, self.vort_b
        fa = self.con_func_a
        fb = self.con_func_b
        sqrt_ratio_a = sqrt(self.frac_b / self.frac_a)
        sqrt_ratio_b = sqrt(self.frac_a / self.frac_b)
        sign = -1
        diag_main1 = (
            nabl
            * (-2.0 / dtht / dtht - (sa + m) ** 2 / sin_th[1 : npts - 1] ** 2)
            + 2 * ga * fa[1 : npts - 1] ** 2
            + sqrt_ratio_a * gab * fb[1 : npts - 1] ** 2
            - mu_a
        )
        diag_main2 = sign * (
            nabl
            * (-2.0 / dtht / dtht - (sa - m) ** 2 / sin_th[1 : npts - 1] ** 2)
            + 2 * ga * fa[1 : npts - 1] ** 2
            + sqrt_ratio_a * gab * fb[1 : npts - 1] ** 2
            - mu_a
        )
        diag_main3 = (
            nabl
            * (-2.0 / dtht / dtht - (sb + m) ** 2 / sin_th[1 : npts - 1] ** 2)
            + 2 * gb * fb[1 : npts - 1] ** 2
            + sqrt_ratio_b * gab * fa[1 : npts - 1] ** 2
            - mu_b
        )
        diag_main4 = sign * (
            nabl
            * (-2.0 / dtht / dtht - (sb - m) ** 2 / sin_th[1 : npts - 1] ** 2)
            + 2 * gb * fb[1 : npts - 1] ** 2
            + sqrt_ratio_b * gab * fa[1 : npts - 1] ** 2
            - mu_b
        )
        diag_main = np.concatenate(
            [diag_main1, diag_main2, diag_main3, diag_main4]
        )
        diag_upp1 = np.zeros(npts - 2)
        diag_upp1[: npts - 3] = nabl * (
            1.0 / dtht / dtht
            + 0.5 * cos_th[1 : npts - 2] / sin_th[1 : npts - 2] / dtht
        )
        diag_upp2 = np.zeros(npts - 2)
        diag_upp2[: npts - 3] = sign * (
            nabl
            * (
                1.0 / dtht / dtht
                + 0.5 * cos_th[1 : npts - 2] / sin_th[1 : npts - 2] / dtht
            )
        )
        diag_upp3 = np.zeros(npts - 2)
        diag_upp3[: npts - 3] = nabl * (
            1.0 / dtht / dtht
            + 0.5 * cos_th[1 : npts - 2] / sin_th[1 : npts - 2] / dtht
        )
        diag_upp4 = sign * (
            nabl
            * (
                1.0 / dtht / dtht
                + 0.5 * cos_th[1 : npts - 2] / sin_th[1 : npts - 2] / dtht
            )
        )
        diag_upp = np.concatenate([diag_upp1, diag_upp2, diag_upp3, diag_upp4])
        diag_low1 = np.zeros(npts - 2)
        diag_low1[: npts - 3] = nabl * (
            1.0 / dtht / dtht
            - 0.5 * cos_th[2 : npts - 1] / sin_th[2 : npts - 1] / dtht
        )
        diag_low2 = np.zeros(npts - 2)
        diag_low2[: npts - 3] = sign * (
            nabl
            * (
                1.0 / dtht / dtht
                - 0.5 * cos_th[2 : npts - 1] / sin_th[2 : npts - 1] / dtht
            )
        )
        diag_low3 = np.zeros(npts - 2)
        diag_low3[: npts - 3] = nabl * (
            1.0 / dtht / dtht
            - 0.5 * cos_th[2 : npts - 1] / sin_th[2 : npts - 1] / dtht
        )
        diag_low4 = sign * (
            nabl
            * (
                1.0 / dtht / dtht
                - 0.5 * cos_th[2 : npts - 1] / sin_th[2 : npts - 1] / dtht
            )
        )
        diag_low = np.concatenate([diag_low1, diag_low2, diag_low3, diag_low4])
        # ---------------------------------------------------------------------
        d1 = ga * fa[1 : npts - 1] ** 2
        d2 = sign * gab * sqrt_ratio_a * fa[1 : npts - 1] * fb[1 : npts - 1]
        d3 = gb * fb[1 : npts - 1] ** 2
        block_diag_upp1 = np.concatenate([d1, d2, d3])
        d1 = gab * sqrt_ratio_a * fa[1 : npts - 1] * fb[1 : npts - 1]
        d2 = sign * gab * sqrt_ratio_a * fa[1 : npts - 1] * fb[1 : npts - 1]
        block_diag_upp2 = np.concatenate([d1, d2])
        block_diag_upp3 = (
            sqrt_ratio_a * gab * fa[1 : npts - 1] * fb[1 : npts - 1]
        )
        # ---------------------------------------------------------------------
        d1 = sign * ga * fa[1 : npts - 1] ** 2
        d2 = gab * sqrt_ratio_b * fa[1 : npts - 1] * fb[1 : npts - 1]
        d3 = sign * gb * fb[1 : npts - 1] ** 2
        block_diag_low1 = np.concatenate([d1, d2, d3])
        d1 = gab * sqrt_ratio_b * fa[1 : npts - 1] * fb[1 : npts - 1]
        d2 = sign * gab * sqrt_ratio_b * fa[1 : npts - 1] * fb[1 : npts - 1]
        block_diag_low2 = np.concatenate([d1, d2])
        block_diag_low3 = sign * (
            sqrt_ratio_b * gab * fa[1 : npts - 1] * fb[1 : npts - 1]
        )
        all_diags = [
            block_diag_low3,
            block_diag_low2,
            block_diag_low1,
            diag_low,
            diag_main,
            diag_upp,
            block_diag_upp1,
            block_diag_upp2,
            block_diag_upp3,
        ]
        offsets = [
            -3 * (npts - 2),
            -2 * (npts - 2),
            -(npts - 2),
            -1,
            0,
            1,
            (npts - 2),
            2 * (npts - 2),
            3 * (npts - 2),
        ]
        return sp.diags(
            all_diags,
            offsets,
            dtype=np.float64,
            shape=(4 * (npts - 2), 4 * (npts - 2)),
        )

    def sparse_matrix(self, m, mu_a, mu_b, ga, gb, gab):
        dtht = self.dtht
        nabl = self.nabla_factor
        sin_th = np.sin(self.tht)
        cos_th = np.cos(self.tht)
        sa, sb = self.vort_a, self.vort_b
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
