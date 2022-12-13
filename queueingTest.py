import unittest

from queuing import *


class ErrorsTest(unittest.TestCase):
    def test_lmda_ge_mu(self):
        self.assertIsNone(lmda_ge_mu(1, 2))
        self.assertRaises(ValueError, lmda_ge_mu, *[2, 1])
        self.assertRaises(ValueError, lmda_ge_mu, *[2, 2])

    def test_lmda_mu_le_0(self):
        self.assertIsNone(lmda_mu_le_0(1, 1))
        self.assertRaises(ValueError, lmda_mu_le_0, *[0, 0])
        self.assertRaises(ValueError, lmda_mu_le_0, *[1, 0])
        self.assertRaises(ValueError, lmda_mu_le_0, *[0, 1])

    def test_lmda_mu_lt_0(self):
        self.assertIsNone(lmda_mu_lt_0(0, 0))
        self.assertIsNone(lmda_mu_lt_0(0, 1))
        self.assertIsNone(lmda_mu_lt_0(1, 0))
        self.assertIsNone(lmda_mu_lt_0(1, 1))
        self.assertRaises(ValueError, lmda_mu_lt_0, *[-1, 1])

    def test_sigma2_lt_0(self):
        self.assertIsNone(sigma2_lt_0(1))
        self.assertIsNone(sigma2_lt_0(0))
        self.assertRaises(ValueError, sigma2_lt_0, *[-1])

    def test_lmda_ge_cmu(self):
        self.assertIsNone(lmda_ge_cmu(0, 2, 1))
        self.assertRaises(ValueError, lmda_ge_cmu, *[1, 1, 1])
        self.assertRaises(ValueError, lmda_ge_cmu, *[2, 1, 1])

    def test_lmda_et_cmu(self):
        self.assertIsNone(lmda_et_cmu(2, 1, 1))
        self.assertIsNone(lmda_et_cmu(0, 1, 1))
        self.assertRaises(ValueError, lmda_et_cmu, *[1, 1, 1])

    def test_c_lt_1(self):
        self.assertIsNone(c_lt_1(2))
        self.assertIsNone(c_lt_1(1))
        self.assertRaises(ValueError, c_lt_1, *[0])

    def test_n_lt_c(self):
        self.assertIsNone(n_lt_c(2, 1))
        self.assertIsNone(n_lt_c(1, 1))
        self.assertRaises(ValueError, n_lt_c, *[0, 1])

    def test_k_lt_c(self):
        self.assertIsNone(k_lt_c(2, 1))
        self.assertIsNone(k_lt_c(1, 1))
        self.assertRaises(ValueError, k_lt_c, *[0, 1])


class QueueingTest(unittest.TestCase):
    def test_eval_MG1(self):
        rho, l, w, wq, lq, p0 = eval_MG1(lmda=1.125, mu=2.35, sigma2=0.2)
        self.assertAlmostEqual(rho, 0.479, delta=0.001)
        self.assertAlmostEqual(l, 0.941, delta=0.001)
        self.assertAlmostEqual(w, 0.837, delta=0.001)
        self.assertAlmostEqual(wq, 0.411, delta=0.001)
        self.assertAlmostEqual(lq, 0.463, delta=0.001)
        self.assertAlmostEqual(p0, 0.521, delta=0.001)

    def test_eval_MMc(self):
        rho, l, w, wq, lq, p0 = eval_MMc(lmda=3.6, mu=2.15, c=3)
        self.assertAlmostEqual(rho, 0.558, delta=0.001)
        self.assertAlmostEqual(l, 2.057, delta=0.001)
        self.assertAlmostEqual(w, 0.571, delta=0.001)
        self.assertAlmostEqual(wq, 0.106, delta=0.001)
        self.assertAlmostEqual(lq, 0.383, delta=0.001)
        self.assertAlmostEqual(p0, 0.171, delta=0.001)

    def test_eval_MGc(self):
        rho, l, w, wq, lq = eval_MGc(lmda=5.42, mu=2.18, c=3, sigma2=0.56)
        self.assertAlmostEqual(rho, 0.829, delta=0.001)
        self.assertAlmostEqual(l, 8.640, delta=0.001)
        self.assertAlmostEqual(w, 1.594, delta=0.001)
        self.assertAlmostEqual(wq, 1.135, delta=0.001)
        self.assertAlmostEqual(lq, 6.153, delta=0.001)

    def test_eval_MMcN(self):
        rho, l, w, wq, lq, p0, pN, lmda_effective = eval_MMcN(lmda=12.98, mu=3.47, c=4, n=15)
        self.assertAlmostEqual(rho, 0.895, delta=0.001)
        self.assertAlmostEqual(l, 7.217, delta=0.001)
        self.assertAlmostEqual(w, 0.581, delta=0.001)
        self.assertAlmostEqual(wq, 0.293, delta=0.001)
        self.assertAlmostEqual(lq, 3.639, delta=0.001)
        self.assertAlmostEqual(p0, 0.011, delta=0.001)
        self.assertAlmostEqual(pN, 0.043, delta=0.001)
        self.assertAlmostEqual(lmda_effective, 12.417, delta=0.001)

    def test_eval_MMcKK(self):
        rho, l, w, wq, lq, p0, lmda_effective = eval_MMcK(lmda=2.65, mu=1.2, c=5, k=6)
        self.assertAlmostEqual(rho, 0.809, delta=0.001)
        self.assertAlmostEqual(l, 4.169, delta=0.001)
        self.assertAlmostEqual(w, 0.859, delta=0.001)
        self.assertAlmostEqual(wq, 0.026, delta=0.001)
        self.assertAlmostEqual(lq, 0.125, delta=0.001)
        self.assertAlmostEqual(p0, 0.001, delta=0.001)
        self.assertAlmostEqual(lmda_effective, 4.853, delta=0.001)


if __name__ == '__main__':
    unittest.main()
