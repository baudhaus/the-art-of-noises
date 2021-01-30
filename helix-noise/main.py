#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author ulrichtsblr

DESCRIPTION:
main.py provides driver code for the helix.py script.

FUNCTIONS:
main

USAGE:
$ python main.py

"""
import subprocess


def main():
    """
    Main function.

    ARGUMENTS:
        None
    RETURNS:
        None
    """
    genomes = [
        "Homo_sapiens.fasta",
        "Mus_musculus.fasta",
        "Drosophila_melanogaster.fasta",
        "Caenorhabditis_elegans.fasta",
        "Arabidopsis_thaliana.fasta",
        "Saccharomyces_cerevisiae.fasta",
        "Escherichia_coli.fasta",
        "Bacteriophage_phiX174.fasta"
    ]
    channels = ["A", "C", "G", "T"]
    for seq in genomes:
        for nt in channels:
            cmd = "python helix.py {} {}".format(seq, nt)
            subprocess.check_call(cmd, shell=True)
    return None


if __name__ == "__main__":
    main()
