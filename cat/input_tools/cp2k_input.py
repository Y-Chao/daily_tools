#!/usr/bin/env python
# -*- encoding: utf-8 -*-

"""
Refer to:
https://github.com/chenggroup/ChecMatE/checmate/pretask/sets.py
"""
__author__ = "Chao Yang"
__version__ = "1.0"

"""
Generate the input files from a template file and updated by user defined config file.
Classes:
    - Cp2kInput: Generate the input files for different cp2k calculations.
Functions:
    - load_config: Load the template config file
    - update_dict: Update the dictionary
    - iterdict: Iteratively parse the dictionary to a list
    - find_key_in_nested_dict: Find the key in the nested dictionary
    - reformat_list: Reformat the index list, such as for 1, 2, 3, 4, 5 -> "1..5"
To do:
    - All keys and subsections should be upper case, which is the lowercase in the value.
    - The basis set and potential set should be more flexible.
"""

import copy
import json
import sys
from pathlib import Path

import numpy as np
import yaml
from ase.atoms import Atoms
from ase.data import chemical_symbols
from ase.io import read

MODULE_DIR = Path(__file__).resolve().parent


def load_config(fpath: str):
    """
    Load template config file
    Parameters:
    ------------
    fpath: str
        The path of the template config file
    """
    fmt = Path(fpath).suffix

    with open(fpath, "r", encoding="utf8") as fd:
        if fmt == ".yaml":
            config = yaml.load(fd, Loader=yaml.FullLoader)
        elif fmt == ".json":
            config = json.load(fd)
        else:
            raise TypeError(
                "An unsupported file format. (Only json and yaml are supported."
            )

    return config


def update_dict(old_dict: dict, new_dict: dict):
    """
    Update the dictionary accroding to user defined directory
    old_dict: dict
        The old dictionary
    new_dict: dict
        The updated dictionary, should be set from the root_directory

    Return:
    -------
    None
    """
    # Here, we update the old_dict iteratively
    for k, v in new_dict.items():
        if k in old_dict and isinstance(old_dict[k], dict) and isinstance(v, dict):
            update_dict(old_dict[k], new_dict[k])
        else:
            # If key not in old_dict,
            #    value in both dict is not a dictory,
            # we enter this branch

            # Ensure the key is not case sensitive
            old_dict_keys = {}
            for i in old_dict:
                old_dict_keys[i.upper()] = i

            if k.upper() in old_dict_keys:
                old_dict[old_dict_keys[k.upper()]] = new_dict[k]
            else:
                old_dict[k] = new_dict[k]


def iterdict(input_dict: dict, outlist: list = ["\n"], loop_level: int = 0):
    """
    Iteratively parse the dictionary to a list.
    Parameters:
    ------------
    input_dict: dict
        The input dictionary
    outlist: list
        The output list
    loop_level: int
        The loop level. Generally, we can think it as the depth of the dictionary, which
        can be regarded as the number of the index between & and &END.
    Returns:
    ---------
    outlist
    """
    if len(outlist) == 0:
        outlist = ["\n"]

    start_index = (
        len(outlist) - loop_level - 2
    )  # len(outlist) - loop_level represents the index of the new &END
    for key, value in input_dict.items():
        key = key.upper()
        if isinstance(value, dict):
            outlist.insert(-1 - loop_level, f"{'  '*loop_level}&" + key)
            outlist.insert(-1 - loop_level, f"{'  '*loop_level}&END " + key)
            iterdict(value, outlist, loop_level + 1)
        elif isinstance(value, list):
            for v in value:
                outlist.insert(-1 - loop_level, f"{'  '*loop_level}&" + key)
                outlist.insert(-1 - loop_level, f"{'  '*loop_level}&END " + key)
                iterdict(v, outlist, loop_level + 1)
        else:
            if value not in chemical_symbols:
                if isinstance(value, str):
                    if len(value.split(".")) > 1:
                        value = str(value).lower()
                    elif value.lower() in ["revpbe"]:
                        value = "revPBE"
                    else:
                        value = str(value).upper()
                else:
                    value = str(value)
            if key == "-":
                outlist[start_index] = outlist[start_index] + " " + value
            else:
                outlist.insert(-1 - loop_level, "  " * loop_level + key + " " + value)

    return outlist


def find_key_in_nested_dict(d, key):
    if key in d:
        return True
    for v in d.values():
        if isinstance(v, dict):  # If value is dict, search in the value
            if find_key_in_nested_dict(v, key):
                return True
    return False


def reformat_list(index):
    """
    Reformat the index list, such as for 1, 2, 3, 4, 5 -> "1..5"
    parameters:
    -----------
    index: list
        The index list
    return:
    -------
    index_list: str
        The reformated index list
    """
    ranges = []
    start = index[0]
    for i in range(1, len(index)):
        if index[i] != index[i - 1] + 1:
            if start == index[i - 1]:
                ranges.append(f"{start}")
            else:
                ranges.append(f"{start}..{index[i-1]}")
            start = index[i]

    # the previous loop does not include the last element
    # we need to add it manually
    if start == index[-1]:
        ranges.append(f"{start}")
    else:
        ranges.append(f"{start}..{index[-1]}")
    return " ".join(ranges)


class Cp2kInput:
    """
    It is used to generate the input files for different cp2k calculations.
    Firstly,
    """

    def __init__(
        self,
        structure: Atoms,
        task: str = "geo_opt",
        xc: str = "PBE",
        vdw: str = None,
        precision: str = "normal",
        user_config: dict = {},
        template_path: str = None,
        **kwargs,
    ):
        self.structure = structure
        self.task = task
        self.precision = precision
        self.user_config = user_config
        self.xc = xc if xc is not None else "PBE"
        self.vdw = vdw
        self.template_config = load_config(template_path)
        self.kwargs = kwargs

    @property
    def cp2kinp(self):
        """Generate the input file"""
        cp2kinp = self.initialize_cp2kinp()

        # update the atomoic information
        cell_config = self.set_Cell()

        # Update the basis and potential information, update kind information
        dft_param_config = self.set_Basis_Potential()

        # Update the kind information
        kind_config = self.set_Kind()

        if len(self.structure.constraints) > 0:
            constraints = self.structure.constraints[0].index
            constraints_config = self.set_Constraints(constraints)
            update_dict(cp2kinp, constraints_config)

        if self.task == "sgcpmd":
            if not find_key_in_nested_dict(self.user_config, "define_region"):
                raise ValueError(
                    "The define_region should be defined in the user_config."
                )

        if self.user_config:
            update_dict(cp2kinp, self.user_config)

        update_dict(cp2kinp, cell_config)
        update_dict(cp2kinp, dft_param_config)
        update_dict(cp2kinp, kind_config)
        if self.xc.lower() != "pbe":
            if self.xc.lower() == "scan":
                xc_config = self.template_config["xc_scan_section"]
            elif self.xc.lower() == "revpbe":
                xc_config = self.template_config["xc_revpbe_section"]
            else:
                raise ValueError("The xc functional is not supported.")
            update_dict(cp2kinp, xc_config)
        if self.vdw is not None:
            vdw_config = self.set_vdw()
            update_dict(cp2kinp, vdw_config)
        return cp2kinp

    def initialize_cp2kinp(self):
        """
        Based on the task type, we initialize the cp2k input file.
        Currently, we support the following tasks:
            - geo_opt
            - cell_opt
            - bomd
            - sgcpmd
            - restart
            - energy

        return:
        --------
        init_cp2kinp: dict
            The initialized cp2k input file
        """

        init_cp2kinp = copy.deepcopy(self.template_config["default_section"])
        if self.task == "geo_opt":
            update_dict(init_cp2kinp, self.template_config["geoopt_section"])
        elif self.task == "cell_opt":
            update_dict(init_cp2kinp, self.template_config["cellopt_section"])
        elif self.task == "bomd":
            update_dict(init_cp2kinp, self.template_config["bomd_section"])
        elif self.task == "sgcpmd":
            update_dict(init_cp2kinp, self.template_config["sgcpmd_section"])
        elif self.task == "restart":
            update_dict(init_cp2kinp, self.template_config["restart_section"])
        elif self.task == "energy":
            pass
        else:
            raise ValueError("The task is not supported.")

        if self.kwargs.get("ot", False):
            print("Using the orbital transformation method.")
            update_dict(init_cp2kinp, self.template_config["ot_section"])
        else:
            print("Using the default diagonalization method.")
            update_dict(init_cp2kinp, self.template_config["diagonalization_section"])

        return init_cp2kinp

    def set_Cell(self):

        cell = self.structure.get_cell()
        cell_a = np.array2string(
            cell[0], formatter={"float_kind": lambda x: "%.5f" % x}
        )[1:-1]
        cell_b = np.array2string(
            cell[1], formatter={"float_kind": lambda x: "%.5f" % x}
        )[1:-1]
        cell_c = np.array2string(
            cell[2], formatter={"float_kind": lambda x: "%.5f" % x}
        )[1:-1]

        cell_config = {
            "force_eval": {
                "subsys": {
                    "cell": {"A": cell_a, "B": cell_b, "C": cell_c},
                    "topology": {
                        "coord_file_name": "coord.xyz",
                        "coord_file_format": "XYZ",
                    },
                }
            }
        }
        return cell_config

    def set_Constraints(self, constraints):
        """
        Update the constraints information
        """
        constraints = reformat_list(constraints)
        constraints_config = {
            "motion": {"constraint": {"fixed_atoms": {"list": constraints}}}
        }
        return constraints_config

    def set_Basis_Potential(self):
        """
        Update the basis and potential information.
        Here, we just use the basis_molopt as the basis set file name.
        It should be updated in the future.
        """

        basis_set_file_name = [
            "BASIS_MOLOPT",
            "BASIS_ADMM",
            "BASIS_ADMM_MOLOPT",
            "BASIS_MOLOPT_HSE06",
        ]
        potential_file_name = "GTH_POTENTIALS"

        basis_potential_config = {
            "force_eval": {
                "dft": {
                    "basis_set_file_name": basis_set_file_name[0],
                    "potential_file_name": potential_file_name,
                }
            }
        }
        return basis_potential_config

    def set_Kind(self):
        """
        Update the kind information
        """
        elements = list(set(self.structure.get_chemical_symbols()))
        kind_config = {"force_eval": {"subsys": {"kind": []}}}

        default_molopt_dict = {
            "Cu": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q11", "scan": "GTH-SCAN-q11"},
            },
            "Au": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q11", "scan": "GTH-SCAN-q11"},
            },
            "Ag": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q11", "scan": "GTH-SCAN-q11"},
            },
            "Pd": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q10", "scan": "GTH-SCAN-q10"},
            },
            "Pt": {
                "basis_set": "DZVP-A5-Q10-323-MOL-T1-DERIVED_SET-1",
                "potential": {"pbe": "GTH-PBE-q10", "scan": "GTH-SCAN-q10"},
            },
            "Ru": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE", "scan": "GTH-SCAN"},
            },
            "Co": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE", "scan": "GTH-SCAN"},
            },
            "Ni": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE", "scan": "GTH-SCAN"},
            },
            "Fe": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE", "scan": "GTH-SCAN"},
            },
            "O": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q6", "scan": "GTH-SCAN-q6"},
            },
            "H": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q1", "scan": "GTH-SCAN-q1"},
            },
            "C": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q4", "scan": "GTH-SCAN-q4"},
            },
            "N": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q5", "scan": "GTH-SCAN-q5"},
            },
            "S": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q6", "scan": "GTH-SCAN-q6"},
            },
            "Li": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q3", "scan": "GTH-SCAN-q3"},
            },
            "Na": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q9", "scan": "GTH-SCAN-q9"},
            },
            "K": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q9", "scan": "GTH-SCAN-q9"},
            },
            "Rb": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q9", "scan": "GTH-SCAN-q9"},
            },
            "Cs": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q9", "scan": "GTH-SCAN-q9"},
            },
            "Cl": {
                "basis_set": "DZVP-MOLOPT-SR-GTH",
                "potential": {"pbe": "GTH-PBE-q7", "scan": "GTH-SCAN-q7"},
            },
        }

        for ele in elements:
            if self.xc.lower() in ["pbe", "revpbe"]:
                xc = "pbe"
            elif self.xc.lower() == "scan":
                xc = "scan"
            else:
                raise ValueError(f"The xc functional {self.xc} is not supported.")
            kind_config["force_eval"]["subsys"]["kind"].append(
                {
                    "-": ele,
                    "basis_set": default_molopt_dict[ele]["basis_set"],
                    "potential": default_molopt_dict[ele]["potential"][xc],
                }
            )
        return kind_config

    def set_vdw(self):
        """
        Update the vdw information
        """
        if self.vdw.lower() == "dftd3":
            vdw_config = self.template_config[f"vdw_{self.xc.lower()}_d3_section"]
        else:
            raise ValueError(f"The vdw method {self.vdw} is not supported.")
        return vdw_config

    def write_cp2kinp(self, fpath: str):
        """
        Write the cp2k input file
        """
        cp2kinp = self.cp2kinp
        cp2kinp_list = iterdict(cp2kinp)

        with open(fpath, "w", encoding="utf8") as fd:
            for line in cp2kinp_list:
                fd.write(line + "\n")


def main():
    if len(sys.argv) > 4 or len(sys.argv) < 1:
        print("Usage: python inputs.py atoms_file system_config_file")
        sys.exit(1)

    atoms_file = sys.argv[1]
    atoms = read(atoms_file)

    if len(sys.argv) == 3:
        user_config = load_config(sys.argv[2])
    else:
        user_config = {}

    Cp2kInput(
        atoms,
        task="bomd",
        precision="normal",
        xc="revpbe",
        vdw="dftd3",
        user_config=user_config,
        template_path=MODULE_DIR / "template" / "cp2k_input.json",
        ot=True,
    ).write_cp2kinp("./cp2k_md.inp")


if __name__ == "__main__":
    main()
