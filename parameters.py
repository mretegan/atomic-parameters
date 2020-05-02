# coding: utf-8
# !/usr/bin/env python

import re
import os
import sys
import glob

import logging
import argparse
import subprocess
import collections

import xraydb

xdb = xraydb.XrayDB()


class odict(collections.OrderedDict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value


class Element(object):
    SUBSHELLS = {
        "3d": {"atomic_numbers_range": (21, 30 + 1), "core_electrons": 18},
        "4d": {"atomic_numbers_range": (39, 48 + 1), "core_electrons": 36},
        "4f": {"atomic_numbers_range": (57, 71 + 1), "core_electrons": 54},
        "5d": {"atomic_numbers_range": (72, 80 + 1), "core_electrons": 68},
        "5f": {"atomic_numbers_range": (89, 103 + 1), "core_electrons": 86},
    }

    def __init__(self, symbol, charge=None):
        self.symbol = symbol
        self.charge = charge

    @property
    def atomic_number(self):
        return xdb.atomic_number(self.symbol)

    @property
    def valence_subshell(self):
        """Name of the valence subshell"""
        atomic_number = self.atomic_number
        for subshell, prop in self.SUBSHELLS.items():
            if atomic_number in range(*prop["atomic_numbers_range"]):
                return subshell
        return None

    @property
    def valence_occupancy(self):
        """Occupancy of the valence subshell"""
        assert self.charge is not None, "The charge must be set."

        # Reverse the string holding the charge before changing it to
        # an integer.
        charge = int(self.charge[::-1])

        # Calculate the number of electrons of the ion.
        ion_electrons = xdb.atomic_number(self.symbol) - charge

        core_electorns = self.SUBSHELLS[self.valence_subshell]["core_electrons"]
        occupancy = ion_electrons - core_electorns
        return occupancy

    def __repr__(self):
        if self.charge is None:
            return "{:s}".format(self.symbol)
        return "{:s}{:s}".format(self.symbol, self.charge)


class Configuration(object):  # noqa
    OCCUPANCIES = {"s": 2, "p": 6, "d": 10, "f": 14}

    def __init__(self, name):
        self.name = name
        self.energy = None
        self.atomic_parameters = None

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):  # noqa
        PATTERNS = (r"^(\d)(\w)(\d+),(\d)(\w)(\d+)$", r"^(\d)(\w)(\d+)$")

        # Test the configuration string.
        tokens = [token for pattern in PATTERNS for token in re.findall(pattern, value)]
        if not tokens:
            raise ValueError("Invalid configuration string.")
        [tokens] = tokens

        if len(tokens) == 3:
            core = None
            valence = tokens
        elif len(tokens) == 6:
            core = tokens[:3]
            valence = tokens[-3:]

        valence_level, valence_shell, valence_occupancy = valence
        valence_level = int(valence_level)
        valence_occupancy = int(valence_occupancy)
        if valence_occupancy > self.OCCUPANCIES[valence_shell]:
            raise ValueError("Wrong number of electrons in the valence shell.")

        if core:
            core_level, core_shell, core_occupancy = core
            core_level = int(core_level)
            core_occupancy = int(core_occupancy)
            if core_occupancy > self.OCCUPANCIES[core_shell]:
                raise ValueError("Wrong number of electrons in the core shell.")

            self.levels = [core_level, valence_level]
            self.shells = [core_shell, valence_shell]
            self.occupancies = [core_occupancy, valence_occupancy]
        else:
            self.levels = [valence_level, ]
            self.shells = [valence_shell, ]
            self.occupancies = [valence_occupancy, ]

        self.subshells = [str(level) + shell for level, shell in zip(self.levels, self.shells)]

        self._name = value

    @classmethod
    def from_subshells_and_occupancies(cls, subshells, occupancies):
        name = ",".join(
            "{:s}{:d}".format(subshell, occupancy)
            for subshell, occupancy in zip(subshells, occupancies)
        )
        return cls(name)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __lt__(self, other):
        return self.name < other.name

    def __repr__(self):
        return self.name


class Cowan(object):
    """Calculate the parameters of an electronic configuration using Cowan's programs."""

    RCN_INPUT = "22 -9    2   10  1.0    5.E-06    1.E-09-2   130   1.0  0.65 11.0 0.50 0.0  0.70\n"
    RCN2_INPUT = """G5INP     000                 00        00000000  9999999999 .00       1229
        -1
    """
    RCN = "scripts/runrcn.sh"
    RCN2 = "scripts/runrcn2.sh"
    RCG = "scripts/runrcg.sh"

    RYDBER_TO_EV = 13.605693122994

    def __init__(self, element, configuration, basename="input"):
        self.element = element
        self.configuration = configuration
        self.basename = basename
        self.dirname = os.path.dirname(__file__)

    @staticmethod
    def normalize_configuration_name(configuration):
        """Configuration name expected by Cowan's programs."""
        occupancies = configuration.occupancies
        subshells = configuration.subshells

        name = str()
        for subshell, occupancy in zip(subshells, occupancies):
            # For 5d elements, the 4f subshell must be included explicitly.
            if "5d" in subshell and "4f" not in subshells:
                subshell = "4f14 5d"
            elif "4f" in subshell and "3d" not in configuration.name:
                subshell = "3d10 4f"
            name += "{0:s}{1:02d} ".format(subshell.upper(), occupancy)
        return name.rstrip()

    def run(self, command):
        """Run the "command"; discard stdout and stderr, but check the exit status."""
        try:
            subprocess.run(
                (command, self.basename),
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True
            )
        except subprocess.CalledProcessError:
            logging.critical("The command %s did not finish successfully.", command)
            sys.exit()

    def rcn(self):
        """Create the input and run the RCN program."""
        rcn_input = self.RCN_INPUT
        for configuration in (self.configuration,):
            line = "{:5d}           {:8s}         {:8s}\n".format(
                self.element.atomic_number,
                configuration.name,
                self.normalize_configuration_name(configuration),
            )
            rcn_input += line
        rcn_input += "{:5d}\n".format(-1)

        filename = "{:s}.rcn".format(self.basename)
        with open(filename, "w") as fp:
            fp.write(rcn_input)

        self.run(os.path.join(self.dirname, self.RCN))

    def rcn2(self):
        """Create the input and run the RCN2 program."""
        filename = "{:s}.rcn2".format(self.basename)
        with open(filename, "w") as fp:
            fp.write(self.RCN2_INPUT)

        self.run(os.path.join(self.dirname, self.RCN2))

    def rcg(self):
        """Create the input and run the RCG program."""
        filename = "{:s}.rcg".format(self.basename)
        # The input file has ".orig" appended to the end.
        os.rename(filename + ".orig", filename)
        with open(filename, "r") as fp:
            lines = fp.readlines()
        with open(filename, "w") as fp:
            for line in lines:
                fp.write(re.sub(r"80998080", r"99999999", line))
        self.run(os.path.join(self.dirname, self.RCG))

    def remove_calculation_files(self):
        filenames = sorted(glob.glob(self.basename + "*"))
        filenames.append("FTN02")
        for filename in filenames:
            try:
                os.remove(filename)
            except FileNotFoundError:
                pass

    def get_parameters(self, debug=False):  # noqa
        self.rcn()
        self.rcn2()
        self.rcg()

        subshells = self.configuration.subshells

        # Parse the output of the RCG program to get the names of the parameters.
        tmp = list()
        filename = "{:s}.rcg_out".format(self.basename)
        with open(filename) as fp:
            for line in fp:
                if "PARAMETER VALUES IN" in line:
                    # Skip two lines
                    for _ in range(2):
                        line = next(fp)
                    while line.split():
                        tokens = re.split(r"\s{2,}", line.strip())
                        tmp.extend(tokens)
                        line = next(fp)

        # Process the names.
        names = list()
        for name in tmp:
            if name.startswith("F") or name.startswith("G"):
                start = name[:2]
                idx1 = int(name[3]) - 1
                idx2 = int(name[4]) - 1
                subshell1 = subshells[idx1]
                subshell2 = subshells[idx2]
                name = "{}({},{})".format(start, subshell1, subshell2)
            elif name.startswith("ZETA"):
                idx = int(name.split()[-1]) - 1
                subshell = subshells[idx]
                name = "ζ({})".format(subshell)
            else:
                continue
            names.append(name)

        # Parse the output of the RCN program to get the values of the parameters.
        values = list()
        filename = "{:s}.rcn_out".format(self.basename)
        with open(filename) as fp:
            for line in fp:
                if "ETOT=" in line:
                    energy = float(line.split()[-1]) * self.RYDBER_TO_EV
                    line = next(fp)

                    # Skip a few empty lines.
                    while not line.split():
                        line = next(fp)

                    # Parse the atomic parameters.
                    tokens = map(float, line.split()[4::2])
                    values.extend(tokens)

                    # In some cases the parameters also span the next two lines.
                    for _ in range(2):
                        line = next(fp)
                        tokens = line.split()
                        if tokens:
                            values.extend(map(float, tokens[::2]))

        parameters = odict()
        for name, value in zip(names, values):
            parameters[name] = value

        # Remove files generated during the calculations.
        if not debug:
            self.remove_calculation_files()

        return energy, parameters


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--element", default="Fe")
    parser.add_argument("-c", "--configuration", default="3d5")
    parser.add_argument("-l", "--loglevel", default="debug")
    parser.add_argument("-d", "--debug", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(format="%(message)s", level=args.loglevel.upper())

    element = Element(args.element)
    conf = Configuration(args.configuration)

    cowan = Cowan(element, conf)
    conf.energy, conf.atomic_parameters = cowan.get_parameters(args.debug)

    logging.info("%2s %-8s", element.symbol, conf)
    logging.info("E = %-.4f eV", conf.energy)
    for parameter, value in conf.atomic_parameters.items():
        logging.debug("%-s = %-.4f eV", parameter, value)
    logging.debug("")


if __name__ == "__main__":
    main()
