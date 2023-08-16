from .fable import fable
from qiskit_aer import AerSimulator
import numpy as np


def test_fable_real():
    n = 3
    a = np.random.randn(2**n, 2**n)

    simulator = AerSimulator(method="unitary")

    circ, alpha = fable(a)
    circ.save_state()
    u_be = simulator.run(circ).result().get_unitary().data
    np.testing.assert_array_almost_equal(
        a/alpha/2**n,  np.real(u_be[:2**n, :2**n])
    )


def test_fable_complex():
    n = 3
    a = np.random.randn(2**n, 2**n) + 1j * np.random.randn(2**n, 2**n)

    simulator = AerSimulator(method="unitary")

    circ, alpha = fable(a)
    circ.save_state()
    u_be = simulator.run(circ).result().get_unitary().data
    np.testing.assert_array_almost_equal(
        a/alpha/2**n,  u_be[:2**n, :2**n]
    )
