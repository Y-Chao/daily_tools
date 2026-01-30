#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = ["Chao Yang", "Jinbo Zhu"]
__version__ = "1.1"


"""
Load the pre-trained model and use it to predict the energy and forces.
"""


import argparse
import os
import socket
import sys
from copy import deepcopy
from shutil import copy2

import numpy as np
from ase import units
from ase.filters import FrechetCellFilter
from ase.io import Trajectory, read, write
from ase.md.logger import MDLogger
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import BFGS


def parse_args():
    parser = argparse.ArgumentParser(
        description="Predict energy and forces using a pre-trained model."
    )
    parser.add_argument(
        "--umlp",
        type=str,
        required=True,
        help="Path of the UMLP model (MACE or Fairchem)",
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="Input file (e.g., .xyz, .cif)"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Output file name (no suffix)"
    )
    parser.add_argument(
        "-m",
        "--model",
        type=str,
        default="oc20",
        choices=["oc20", "omat", "omol", "odac", "omc"],
        help="Type of the model",
    )
    parser.add_argument(
        "-b", "--batch_size", type=int, default=32, help="Batch size for prediction"
    )
    parser.add_argument(
        "-t",
        "--task",
        type=str,
        default="sp",
        choices=["sp", "geo_opt", "cell_opt", "md"],
        help="Task to perform",
    )
    parser.add_argument(
        "-f",
        "--force",
        type=float,
        default=0.05,
        help="Force tolerance for optimization",
    )
    parser.add_argument(
        "-n", "--num_steps", type=int, default=1000, help="Number of optimization steps"
    )
    parser.add_argument(
        "-s",
        "--structure",
        type=str,
        default="bulk",
        choices=["bulk", "surface"],
        help="Type of structure for cell optimization",
    )
    parser.add_argument(
        "--temperature", type=float, default=300, help="Input MD temperature (K)"
    )
    parser.add_argument("--timestep", type=float, default=1, help="Time steps in fs")
    parser.add_argument(
        "--pressure",
        type=float,
        default=1.01325,
        help="The desired pressure, in bar (1 bar = 1e5 Pa). Default is 1.01325",
    )
    parser.add_argument(
        "--compressibility",
        type=float,
        default=4.57e-5,
        help="The compressibility of the material, in bar-1, Default is 4.57e-5",
    )
    parser.add_argument(
        "--ensemble",
        type=str,
        default="nve",
        choices=[
            "nve",
            "nvt_langevin",
            "nvt_andersen",
            "nvt_berendsen",
            "npt_berendsen",
            "npt_melchionna",
            "npt",
        ],
        help="Choose ensemble from: nve, nvt_langevin, nvt_andersen, nvt_berendsen, npt_berendsen",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cpu",
        choices=["cpu", "cuda"],
        help="Use GPU for prediction",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Load the pre-trained model
    if "mace" in args.umlp.lower():
        umlp_type = "MACE"
    elif "fairchem" in args.umlp.lower():
        umlp_type = "Fairchem"
    elif "vasp" in args.umlp.lower():
        umlp_type = "VASP"
    else:
        raise ValueError("Only 'MACE', 'Fairchem' and 'VASP' could be support")

    if umlp_type == "MACE":
        from mace.calculators import mace_mp

        calculator = mace_mp(model=args.umlp, device=args.device)

    elif umlp_type == "Fairchem":
        from fairchem.core import FAIRChemCalculator
        from fairchem.core.units.mlip_unit import load_predict_unit
        from fairchem.core.units.mlip_unit.api.inference import InferenceSettings

        # Ref: https://fair-chem.github.io/core/common_tasks/ase_calculator.html
        settings = InferenceSettings(
            tf32=True,
            activation_checkpointing=False,
            merge_mole=False,
            compile=False,
            wigner_cuda=True,
            external_graph_gen=False,
            internal_graph_gen_version=2,
        )
        model = load_predict_unit(
            args.umlp, device=args.device, inference_settings=settings
        )
        calculator = FAIRChemCalculator(model, task_name=args.model)

    elif umlp_type == "VASP":
        if not os.path.exists(args.umlp) and not args.umlp.endswith(".py"):
            hostname = socket.gethostname()
            if hostname.startswith("asp2a"):
                src_file = "/home/project/31010029/opt/scripts"
            elif hostname.startswith("stdct"):
                src_file = "/scratch/projects/CFP01/CFP01-SF-002/opt/fairchem/"
            else:
                raise ValueError("No such a server")
            copy2(os.path.join(src_file, "run_vasp_sp.py"), "./")
            sys.exit("Check the parameter for vasp calculatr")
        else:
            copy2(args.umlp, "./run_vasp_sp.py")

        try:
            from run_vasp_sp import vasp_calc
        except ImportError:
            raise ImportError("Cannot import vasp_calc from run_vasp_sp.py")

        calculator = deepcopy(vasp_calc)
    else:
        raise ValueError("Only 'MACE', 'Fairchem' and 'VASP' could be support")

    # Read the input structure
    atoms = read(args.input)
    atoms.calc = calculator

    # Set the calculator
    if args.task == "sp":
        # Single point energy calculation
        print("Single point energy calculation started.")
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        write(f"{args.output}.traj", atoms)

        print(f"Energy: {energy} eV")
        # print(f"Forces: {forces} eV/Å")
        print(f"Max forces: {np.max(np.linalg.norm(forces, axis=1))} eV/Å")
        print("Single point energy calculation completed.")

    elif args.task == "geo_opt":
        # Geometry optimization
        print("Geometry optimization started.")
        opt = BFGS(atoms)
        traj = Trajectory(f"{args.output}.traj", "w", atoms)
        opt.attach(traj)
        opt.run(fmax=args.force, steps=args.num_steps)
        print(f"Energy: {atoms.get_potential_energy()} eV")
        # print(f"Forces: {atoms.get_forces()} eV/Å")
        print(f"Max forces: {np.max(np.linalg.norm(atoms.get_forces(), axis=1))} eV/Å")
        print("Geometry optimization completed.")

    elif args.task == "cell_opt":
        # Cell optimization
        print("Cell optimization started.")
        if args.structure == "bulk":
            mask_c = [1, 1, 1, 1, 1, 1]
        elif args.structure == "surface":
            mask_c = [1, 1, 0, 0, 0, 1]
        else:
            raise ValueError("Invalid structure type. Use 'bulk' or 'surface'.")
        if args.pressure != 1.01325:
            ecf = FrechetCellFilter(
                atoms, mask=mask_c, scalar_pressure=args.pressure * units.GPa
            )
        else:
            ecf = FrechetCellFilter(atoms, mask=mask_c)
        opt = BFGS(ecf)
        traj = Trajectory(f"{args.output}.traj", "w", atoms)
        opt.attach(traj)
        opt.run(fmax=args.force, steps=args.num_steps)
        print(f"Energy: {atoms.get_potential_energy()} eV")
        # print(f"Forces: {atoms.get_forces()} eV/Å")
        print(f"Max forces: {np.max(np.linalg.norm(atoms.get_forces(), axis=1))} eV/Å")
        print("Cell optimization completed.")

    elif args.task == "md":
        # Molecular dynamic
        print("Molecular dynamic started.")
        MaxwellBoltzmannDistribution(atoms, temperature_K=args.temperature)

        if args.ensemble == "nve":
            from ase.md.verlet import VelocityVerlet

            dyn = VelocityVerlet(atoms, args.timestep * units.fs)
        elif args.ensemble == "nvt_langevin":
            print("Using Langevin thermostat for NVT ensemble.")
            print("   Suitable for mixed system, while harder to reproduce results.")
            print(
                "   Recommend friction=0.1 / fs with strong friction, fast converagence but lost physical properties; 0.005 / fs similar to NVE."
            )
            from ase.md.langevin import Langevin

            dyn = Langevin(
                atoms,
                timestep=args.timestep * units.fs,
                temperature_K=args.temperature,
                friction=0.02 / units.fs,
            )
        elif args.ensemble == "nose_hoover":
            print("Using Nose-Hoover chain thermostat for NVT ensemble.")
            print("   Suitable for precise NVT ensemble sampling and universal.")
            print(
                "   Recommend tdamp=100~200 fs for mixed system, 50~100 fs for solid system."
            )
            from ase.md.nose_hoover_chain import NoseHooverChainNVT

            dyn = NoseHooverChainNVT(
                atoms,
                timestep=args.timestep * units.fs,
                temperature_K=args.temperature,
                tdamp=200 * units.fs,
                tchain=3,
                tloop=1,
            )
        elif args.ensemble == "nvt_andersen":
            # modify "Andersen probability" based on your preference
            print("Using Andersen thermostat for NVT ensemble is not recommended.")

            from ase.md.andersen import Andersen

            dyn = Andersen(
                atoms, args.timestep * units.fs, args.temperature, andersen_prob=0.002
            )
        elif args.ensemble == "nvt_berendsen":
            # modify "taut: Time constant for Berendsen temperature coupling" based on your preference
            print("Using Berendsen thermostat for NVT ensemble.")
            print(
                "   Recommend taut=0.5~1 ps for mixed system, 0.1~0.5 ps for solid system. Smaller value means quicker temperature adjustment."
            )
            from ase.md.nvtberendsen import NVTBerendsen

            dyn = NVTBerendsen(
                atoms,
                args.timestep * units.fs,
                args.temperature,
                taut=0.5 * 1000 * units.fs,
            )
        elif args.ensemble == "npt_berendsen":
            # modify "taut: Time constant for Berendsen temperature coupling"
            #        "taup: Time constant for Berendsen pressure coupling" based on your preference
            print("Using Berendsen thermostat for NPT ensemble.")
            print(
                "   Recommend taut=0.5~1 ps for mixed system, 0.1~0.5 ps for solid system. Smaller value means quicker temperature adjustment."
            )
            print(
                "   Recommend taup=0.5~2 ps for solid system, 2-5 ps for mixed system. Smaller value means quicker volume adjustment."
            )
            print(
                "   Recommend compressibility_au=1e-5~1e-4 for solid system, 1e-5~5e-5 for mixed system."
            )
            from ase.md.nptberendsen import NPTBerendsen

            dyn = NPTBerendsen(
                atoms,
                timestep=args.timestep * units.fs,
                temperature_K=args.temperature,
                taut=0.5 * 1000 * units.fs,
                pressure_au=args.pressure * units.bar,
                taup=5 * 1000 * units.fs,
                compressibility_au=args.compressibility / units.bar,
            )
        elif args.ensemble == "npt_melchionna" or args.ensemble == "npt":
            # modify "tau: Time constant for the barostat" based on your preference
            print("Using Martyna-Tobias-Klein method for NPT ensemble.")
            print(
                "   Recommend ttime=0.01~2 ps for solid system, 0.5-2 ps for mixed system. Smaller value means quicker volume adjustment."
            )
            from ase.md.npt import NPT

            dyn = NPT(
                atoms,
                timestep=args.timestep * units.fs,
                temperature_K=args.temperature,
                ttime=1 * 1000 * units.fs,
                pfactor=None,
                externalstress=args.pressure * units.bar,
            )
        else:
            raise ValueError("Invalid ensemble setting.")

        dyn.run(steps=300)
        md_logger = MDLogger(
            dyn,
            atoms,
            logfile=f"{args.output}.log",
            header=True,
            stress=False,
            mode="w",
        )
        dyn.attach(md_logger, interval=100)
        trajectory = Trajectory(f"{args.output}.traj", "w", atoms)
        dyn.attach(trajectory.write, interval=1)
        dyn.run(steps=args.num_steps)
        print("Molecular dynamic completed.")
    else:
        raise ValueError("Invalid task specified.")


if __name__ == "__main__":
    main()
