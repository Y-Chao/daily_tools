#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"

import ase.calculators.vasp as vasp_calculator
from ase.io import read

# kpts parameter
kpt = (2, 2, 1)

# Generic VASP calculator common for all cases
calc = vasp_calculator.Vasp(
    # Performance
    istart=0,  # Start from WAVECAR
    kpar=1,
    ncore=6,  # Parallelization
    icharg=2,
    lreal="Auto",
    # kpar=1,    #Parallelization over nodes
    # ispin=1,
    # npar=4,
    isym=0,
    # DFT parameters
    encut=500,  # PW cutoff
    gamma=False,
    kpts=kpt,  # kpts
    xc="rpbe",  # xc for pseudopotential
    # gga='RP',         #actual xc
    ismear=1,  # Gaussian smearing
    sigma=0.2,  # Fermi temperature
    # Structure optimization
    ibrion=1,  # optimization method 2=cg
    # potim=0.2,
    # CONVERGENCE
    algo=False,
    ediff=1e-5,
    ediffg=2e-2,
    prec="Accurate",
    nelm=350,
    nsw=200,
    # OUTPUT
    lvhar=False,  # write hartree potential
    lwave=False,  # write WAVCAR
    lcharg=False,  # write CHARGECAR
    # laechg=True,
    # Dipole correction
    # ldipol=True,
    idipol=3,  # Apply field along z axis
    dipol=(0.5, 0.5, 0.5),
    # VASPSOL
    lsol=False,  # Vaspsol
    # lambda_d_k=3.0,    #Debye screening length
    # tau=0,             #surface tension (apparently creates problem if nonzero)
    # nelect= nelecs0-2*charge,      #Total number of electrons
    # OTHER
    ivdw=11,
    lplane=True,
    # lasph=True,
    # lorbit=11,
)

# Workflows for the coverage research
atom = read("init.vasp")
atom.calc = calc
energy = atom.get_potential_energy()

with open("out.energy", "w") as fd:
    fd.write(str(energy))
