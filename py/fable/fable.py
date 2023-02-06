#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from qiskit import QuantumCircuit
from ._util import compressed_uniform_rotation, sfwht, gray_permutation


def fable(a, epsilon):
    '''FABLE - Fast Approximate BLock Encodings.

    Args:
        a: array
            matrix to be block encoded.
        epsilon: float >= 0
            compression threshold.
    Returns:
        circuit: qiskit circuit
            circuit that block encodes A
        alpha: float
            subnormalization factor
    '''

    epsm = np.finfo(a.dtype).eps
    alpha = np.linalg.norm(np.ravel(a), np.inf)
    if alpha > 1:
        alpha = alpha + np.sqrt(epsm)
        a = a/alpha
    else:
        alpha = 1.0

    n, m = a.shape
    if n != m:
        k = max(n, m)
        a = np.pad(a, ((0, k - n), (0, k - m)))
        n = k
    logn = int(np.ceil(np.log2(n)))
    if n < 2**logn:
        a = np.pad(a, ((0, 2**logn - n), (0, 2**logn - n)))
        n = 2**logn

    a = np.ravel(a)

    if all(np.abs(np.imag(a)) < epsm):  # real data
        a = gray_permutation(
                sfwht(
                    2.0 * np.arccos(np.real(a))
                )
            )
        # threshold the vector
        a[abs(a) <= epsilon] = 0
        # compute circuit
        OA = compressed_uniform_rotation(a)
    else:  # complex data
        # magnitude
        a_m = gray_permutation(
                sfwht(
                    2.0 * np.arccos(np.abs(a))
                )
            )
        a_m[abs(a_m) <= epsilon] = 0

        # phase
        a_p = gray_permutation(
                sfwht(
                    -2.0 * np.angle(a)
                )
            )
        a_p[abs(a_p) <= epsilon] = 0

        # compute circuit
        OA = compressed_uniform_rotation(a_m) + \
            compressed_uniform_rotation(a_p, ry=False)

    circ = QuantumCircuit(2*logn + 1)

    # diffusion on row indices
    for i in range(logn):
        circ.h(i+1)

    # matrix oracle
    circ = circ.compose(OA)

    # swap register
    for i in range(logn):
        circ.swap(i+1,  i+logn+1)

    # diffusion on row indices
    for i in range(logn):
        circ.h(i+1)

    # reverse bits because of little-endiannes
    circ = circ.reverse_bits()

    return circ, alpha
