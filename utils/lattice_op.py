#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"


"""
Basic lattice operation utils.
"""


def random_orthogonal_matrix(dim: int = 3) -> np.ndarray:
    """
    Generate a random orthogonal matrix using Gram-Schmidt process.

    Parameters:
        dim (int): Dimension of the orthogonal matrix.

    Returns:
        np.ndarray: A random orthogonal matrix of shape (dim, dim).
    """
    import numpy as np

    # Generate a random matrix
    A = np.random.rand(dim, dim)

    # Apply Gram-Schmidt process
    Q, R = np.linalg.qr(A)
    return Q


if __name__ == "__main__":
    main()
