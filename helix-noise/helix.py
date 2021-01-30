#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author ulrichtsblr

DESCRIPTION:
helix.py transforms a DNA sequence into audio.

FUNCTIONS:
import_data
read_fasta
decode
write_wav
main

USAGE:
$ python helix.py <fasta_file_name> <channel>

"""
import sys
import numpy as np
from scipy.io import wavfile


def import_data(file_name):
    """
    Imports a plain text file into the namespace.

    ARGUMENTS:
        file_name -- name of the file to be imported [string]
    RETURNS:
        data -- list of lines in the file, excluding line feed [list]
    """
    data = []

    try:
        wrapper = open(file_name, "r")
        for raw_line in wrapper:
            line = raw_line[:-1]
            data.append(line)
        wrapper.close()
    except IOError:
        sys.exit("invalid file name")

    return data


def read_fasta(file_name):
    """
    Imports and parses a FASTA file into a DNA sequence.

    ARGUMENTS:
        file_name -- name of the FASTA file to be imported [string]
    RETURNS:
        seq -- DNA sequence [string]
    """
    # import fasta file
    data = import_data(file_name)

    # parse fasta file
    seq = ""

    for line in data:
        if not line.strip():            # if empty line, then skip
            continue
        elif line.startswith(">"):      # if header, then skip
            continue
        else:                           # if sequence, then append
            subseq = "".join(line.strip())
            seq += subseq

    seq = seq.upper()

    return seq


def decode(seq, nt, window_length):
    """
    Decodes a DNA sequence into a channel.

    ARGUMENTS:
        seq -- DNA sequence [string]
        nt -- target channel nucleotide (A, C, G or T) [character]
        window_length -- length of the sliding window [integer]
    RETURNS:
        channel -- target channel (A, C, G or T) [numpy.ndarray]
    """
    nwindows = len(seq) - window_length + 1
    channel = np.zeros(nwindows)

    for i in range(nwindows):
        window = seq[i : i + window_length]     # sliding window
        amp = []                                # waveform amplitude
        for j in window:
            if j == nt:
                amp.append("1")
            else:
                amp.append("0")
        amp = "".join(amp)
        amp = int(amp, base=2)                  # cast base 2 as base 10
        channel[i] = amp

    return channel


def write_wav(file_name, sample_rate, dBFS, sample_data):
    """
    Writes a channel to a WAV file.

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.wavfile.write.html

    ARGUMENTS:
        file_name -- name of the output WAV file [string]
        sample_rate -- audio sample rate [integer]
        dBFS -- decibels relative to full scale [float]
        sample_data -- channel [numpy.ndarray]
    RETURNS:
        None
    """
    # internal parameters
    bit_depth = 16
    r = (10 ** (dBFS / 20))
    c = np.max(sample_data)

    # processing
    sample_data = (((sample_data * (2 ** bit_depth)) / c) - (2 ** (bit_depth - 1))) * r
    sample_data = np.int16(sample_data)

    # assertions
    print("\tassertions:")
    print("\t(bit_depth, r, c) = ({}, {}, {})".format(bit_depth, r, c))
    print("\tobject_type = {}".format(type(sample_data)))
    print("\tdtype = {}".format(sample_data.dtype))
    print("\thead4 = {}".format(sample_data[0:4]))

    # output
    wavfile.write(file_name, sample_rate, sample_data)

    return None


def main():
    """
    Main function.

    ARGUMENTS:
        argv -- command line arguments [list]
    RETURNS:
        None
    """
    # command line arguments
    args = sys.argv[1:]
    file_name = args[0]
    channel = args[1].upper()

    # parameters
    window_length = 16      # genomic look-ahead
    sample_rate = 44100     # audio sample rate [Hz]
    dBFS = -12.00           # decibels relative to full scale

    # parse genome
    print("parsing: {} ...".format(file_name))
    seq = read_fasta(file_name)

    # decode genome
    print("decoding: {} ...".format(file_name))
    sample_data = decode(seq, channel, window_length)

    # write to .wav file
    print("waving: {} ...".format(file_name))
    outfile_name = file_name.partition(".")[0] + "_channel_{}.wav".format(channel)
    write_wav(outfile_name, sample_rate, dBFS, sample_data)

    # done
    print("done: {}".format(outfile_name))

    return None


if __name__ == "__main__":
    main()
