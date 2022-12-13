import unittest
from queuing import *


class QueueingTest(unittest.TestCase):
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


if __name__ == '__main__':
    unittest.main()
