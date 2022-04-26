# Fast Approximate BLock Encodings (FABLE)

FABLE can synthesize quantum circuits for approximate block-encodings of matrices. A block-encoding is the embedding of a matrix in the leading block of of a larger unitary matrix.

FABLE is a quantum representation for dense, unstructured matrices. The gate complexity of FABLE circuits scales linear in the number of matrix elements, which is optimal for the unstructured case. FABLE includes a circuit compression algorithm that can significantly reduce the gate complexity and works particularly well if there is certain structure available in the matrix to be block encoded.

We provide two reference implementations of the FABLE algorithm: 
* a Python implementation built on top of [Qiskit](https://qiskit.org/)
* a MATLAB implementation built on top of [QCLAB](https://github.com/QuantumComputingLab/qclab)


## Qiskit - Python Implementation

FABLE can be installed from PyPI as follows:

```
pip install fable-circuits
```

After installation, it can be loaded and used as follows:

```py
from fable import fable
import numpy as np
from qiskit import Aer
simulator = Aer.get_backend("unitary_simulator")


# generate a random matrix and block encode it
n = 3
N = 2**n
A = np.random.randn(N, N)
circ, alpha = fable(A, 0)
result = simulator.run(circ).result()
unitary = result.get_unitary(circ)
np.linalg.norm(alpha * N * unitary.data[0:N, 0:N] - A)/np.linalg.norm(A)
```

### QCLAB - MATLAB Implementation ###

In order to run the MATLAB implementation of FABLE:

1. Install [QCLAB](https://github.com/QuantumComputingLab/qclab)
2. Clone FABLE and add `fable-qclab` directory to your MATLAB path.

After installation, FABLE can be run for a target matrix `A` as either:

```
logging = true ;
[circuit, OA, alpha, info] = fable( A, 'cutoff', 1e-4, logging ) ;
[circuit, OA, alpha, info] = fable( A, 'percentage', 80, logging ) ;
```        
The first option (`'cutoff'`) ignores coefficients smaller than `1e-4` in absolute value, the second option
(`'percentage'`) applies an 80% compression and only retains the 20% largest coefficients. The `'percentage'` and `logging` options are only available in the MATLAB version of FABLE.

## Reference

Cite the following reference for FABLE:

*FABLE: Fast Approximate Quantum Circuits for Block-Encodings*, Daan Camps, Roel Van Beeumen, 2022.

## Developers - Lawrence Berkeley National Laboratory
- [Daan Camps](http://campsd.github.io/) - dcamps@lbl.gov
- [Roel Van Beeumen](http://www.roelvanbeeumen.be/) - rvanbeeumen@lbl.gov