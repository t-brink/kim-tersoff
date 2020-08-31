#!/usr/bin/env python3
#
# Copyright (c) 2015,2020 Tobias Brink
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""Converts parameters for Tersoff-type potentials to the form used by LAMMPS.

Supports the following "styles":

* original Tersoff       [Tersoff, PRB 37, 6991 (1988)]
* multicomponent Tersoff [Tersoff, PRB 39, 5566 (1989)]
* Albe                   [Albe et al., PRB65, 195124 (2002)]

Usage: tersoff2lammps.py [input] [output].

Input is an ini-file with the following sections:

[settings]
style = ...   # one of tersoff1988, tersoff1989, albe2002
comment = ... # a comment that appears at the top of the output file.

[X-Y] # One for each interaction, Si-Si, C-C, Si-C. tersoff1989 uses mixing
      # rules to calculate the cross-terms, only chi must be provided.
<params>


List of params:

tersoff1988           tersoff1989               albe2002
A (eV)                A (eV)                    gamma
B (eV)                B (eV)                    S
lam1 (A^-1)           lambda (A^-1)             beta (A^-1)
lam2 (A^-1)           mu (A^-1)                 D0 (eV)
beta                  beta                      r0 (A)
n                     n                         c
c                     c                         d
d                     d                         h
h                     h                         mu
lam3 (A^-1)           R (A)                     Rcut (A)
R (A)                 S (A)                     Dcut (A)
D (A)                 chi (cross-terms only)


Variants:

tersoff1989full: instead of giving chi for the cross terms, provide
                 the full parameter set.


ZBL:

If the parameters Z, ZBLcut, and ZBLexpscale are given (see LAMMPS
documentation for their meaning), the file will be compatible with
LAMMPS' tersoff/zbl variant.

"""

import configparser
import sys
import math

cp = configparser.ConfigParser()
cp.optionxform = lambda option: option # case-sensitive parser!
cp.read(sys.argv[1])


# Set style.
if cp["settings"]["style"] == "tersoff1988":
    elem = None
    sec = None
    for section in cp:
        if section in ("settings", "DEFAULT"): continue
        if elem is not None:
            print("tersoff1988 only supports a single element!")
            sys.exit(2)
        try:
            elem1, elem2 = section.split("-")
        except ValueError:
            print("section headings must follow the format element-element")
            sys.exit(3)
        if elem1 != elem2:
            print("tersoff1988 does not support cross-interactions ({})!"
                  "".format(section))
            sys.exit(4)
        elem = elem1
        sec = cp[section]
    m = 3
    gamma = 1
    lambda3 = sec["lam3"]
    c = sec["c"]
    d = sec["d"]
    costheta0 = sec["h"]
    n = sec["n"]
    beta = sec["beta"]
    lambda2 = sec["lam2"]
    B = sec["B"]
    R = sec["R"]
    D = sec["D"]
    lambda1 = sec["lam1"]
    A = sec["A"]
    with open(sys.argv[2], "w") as f:
        if "comment" in cp["settings"]:
            f.write("# {}\n\n".format(cp["settings"]["comment"]))
        f.write(" ".join(str(i)
                         for i in [elem, elem, elem,
                                   m, gamma, lambda3, c, d, costheta0,
                                   n, beta, lambda2, B, R, D, lambda1, A]))
        # Do we have ZBL parameters?
        if "Z" in sec:
            Zi = sec["Z"]
            Zj = sec["Z"]
            ZBLcut = sec["ZBLcut"]
            ZBLexpscale = sec["ZBLexpscale"]
            f.write(" ")
            f.write(" ".join(str(i) for i in [Zi, Zj,
                                              ZBLcut, ZBLexpscale]))
        f.write("\n")
elif cp["settings"]["style"] == "tersoff1989":
    elements = set()
    terms = set()
    params = {}
    # Find all elements.
    for section in cp:
        if section in ("settings", "DEFAULT"): continue
        try:
            elem1, elem2 = section.split("-")
        except ValueError:
            print("section headings must follow the format element-element")
            sys.exit(3)
        pair = frozenset([elem1, elem2])
        if pair in terms:
            print("Section is double: {}".format(section))
            sys.exit(5)
        elements.add(elem1)
        elements.add(elem2)
        terms.add(pair)
        #
        params[pair] = cp[section]
    # Check for completeness.
    for elem1 in elements:
        for elem2 in elements:
            p = frozenset([elem1, elem2])
            if p not in params:
                print("Parameters missing for {}-{}.".format(elem1, elem2))
                sys.exit(6)
    # Write.
    with open(sys.argv[2], "w") as f:
        if "comment" in cp["settings"]:
            f.write("# {}\n\n".format(cp["settings"]["comment"]))
        m = 3
        gamma = 1
        lambda3 = 0
        for elem1 in sorted(elements):
            for elem2 in sorted(elements):
                for elem3 in sorted(elements):
                    sec1 = params[frozenset([elem1, elem1])]
                    sec2 = params[frozenset([elem2, elem2])]
                    sec3 = params[frozenset([elem3, elem3])]
                    # Cutoff.
                    R1 = float(sec1["R"])
                    R3 = float(sec3["R"])
                    S1 = float(sec1["S"])
                    S3 = float(sec3["S"])
                    R_ = math.sqrt(R1 * R3)
                    S_ = math.sqrt(S1 * S3)
                    R = (R_ + S_) / 2
                    D = (S_ - R_) / 2
                    # Other params.
                    if elem1 == elem2:
                        sec = sec1
                        if elem2 != elem3:
                            n = beta = lambda2 = B = lambda1 = A = 0
                        else:
                            n = sec["n"]
                            beta = sec["beta"]
                            lambda2 = sec["mu"]
                            B = sec["B"]
                            lambda1 = sec["lambda"]
                            A = sec["A"]
                        f.write(
                            " ".join(
                                str(i) for i in [
                                    elem1, elem2, elem3, m, gamma, lambda3,
                                    sec["c"], sec["d"], sec["h"],
                                    n, beta, lambda2, B,
                                    R, D, lambda1, A
                                ]))
                    else:
                        chi = float(params[frozenset([elem1, elem2])]["chi"])
                        if elem2 != elem3:
                            n = beta = lambda2 = B = lambda1 = A = 0
                        else:
                            mu1 = float(sec1["mu"])
                            mu2 = float(sec2["mu"])
                            B1 = float(sec1["B"])
                            B2 = float(sec2["B"])
                            l1 = float(sec1["lambda"])
                            l2 = float(sec2["lambda"])
                            A1 = float(sec1["A"])
                            A2 = float(sec2["A"])
                            #
                            n = sec1["n"]
                            beta = sec1["beta"]
                            lambda2 = (mu1 + mu2) / 2
                            B = chi * math.sqrt(B1 * B2)
                            lambda1 = (l1 + l2) / 2
                            A = math.sqrt(A1 * A2)
                        f.write(
                            " ".join(
                                str(i) for i in [
                                    elem1, elem2, elem3, m, gamma, lambda3,
                                    sec1["c"], sec1["d"], sec1["h"],
                                    n, beta, lambda2, B,
                                    R, D, lambda1, A
                                ]))
                    # Do we have ZBL parameters?
                    sec12 = params[frozenset([elem1, elem2])]
                    if "Z" in sec12:
                        if elem2 != elem3:
                            Zi = Zj = ZBLcut = ZBLexpscale = 0
                        else:
                            Zi = sec1["Z"]
                            Zj = sec2["Z"]
                            ZBLcut = sec12["ZBLcut"]
                            ZBLexpscale = sec12["ZBLexpscale"]
                        f.write(" ")
                        f.write(" ".join(str(i) for i in [Zi, Zj,
                                                          ZBLcut, ZBLexpscale]))
                    f.write("\n")
elif cp["settings"]["style"] == "tersoff1989full":
    elements = set()
    terms = set()
    params = {}
    # Find all elements.
    for section in cp:
        if section in ("settings", "DEFAULT"): continue
        try:
            elem1, elem2 = section.split("-")
        except ValueError:
            print("section headings must follow the format element-element")
            sys.exit(3)
        pair = frozenset([elem1, elem2])
        if pair in terms:
            print("Section is double: {}".format(section))
            sys.exit(5)
        elements.add(elem1)
        elements.add(elem2)
        terms.add(pair)
        #
        params[pair] = cp[section]
    # Check for completeness.
    for elem1 in elements:
        for elem2 in elements:
            p = frozenset([elem1, elem2])
            if p not in params:
                print("Parameters missing for {}-{}.".format(elem1, elem2))
                sys.exit(6)
    # Write.
    with open(sys.argv[2], "w") as f:
        if "comment" in cp["settings"]:
            f.write("# {}\n\n".format(cp["settings"]["comment"]))
        m = 3
        gamma = 1
        lambda3 = 0
        for elem1 in sorted(elements):
            for elem2 in sorted(elements):
                for elem3 in sorted(elements):
                    sec = params[frozenset([elem1, elem2])]
                    sec11 = params[frozenset([elem1, elem1])]
                    sec22 = params[frozenset([elem2, elem2])]
                    sec13 = params[frozenset([elem1, elem3])]
                    # Cutoff.
                    R_ = float(sec13["R"])
                    S_ = float(sec13["S"])
                    R = (R_ + S_) / 2
                    D = (S_ - R_) / 2
                    # Other params.
                    if elem2 != elem3:
                        n = beta = lambda2 = B = lambda1 = A = 0
                    else:
                        n = sec11["n"]
                        beta = sec11["beta"]
                        lambda2 = sec["mu"]
                        B = sec["B"]
                        lambda1 = sec["lambda"]
                        A = sec["A"]
                    f.write(
                        " ".join(
                            str(i) for i in [
                                elem1, elem2, elem3, m, gamma, lambda3,
                                sec11["c"], sec11["d"], sec11["h"],
                                n, beta, lambda2, B,
                                R, D, lambda1, A
                            ]))
                    # Do we have ZBL parameters?
                    if "Z" in sec:
                        if elem2 != elem3:
                            Zi = Zj = ZBLcut = ZBLexpscale = 0
                        else:
                            Zi = sec11["Z"]
                            Zj = sec22["Z"]
                            ZBLcut = sec["ZBLcut"]
                            ZBLexpscale = sec["ZBLexpscale"]
                        f.write(" ")
                        f.write(" ".join(str(i) for i in [Zi, Zj,
                                                          ZBLcut, ZBLexpscale]))
                    f.write("\n")
elif cp["settings"]["style"] == "albe2002":
    elements = set()
    terms = set()
    params = {}
    # Find all elements.
    for section in cp:
        if section in ("settings", "DEFAULT"): continue
        try:
            elem1, elem2 = section.split("-")
        except ValueError:
            print("section headings must follow the format element-element")
            sys.exit(3)
        pair = frozenset([elem1, elem2])
        if pair in terms:
            print("Section is double: {}".format(section))
            sys.exit(5)
        elements.add(elem1)
        elements.add(elem2)
        terms.add(pair)
        #
        params[pair] = cp[section]
    # Check for completeness.
    for elem1 in elements:
        for elem2 in elements:
            p = frozenset([elem1, elem2])
            if p not in params:
                print("Parameters missing for {}-{}.".format(elem1, elem2))
                sys.exit(6)
    # Write.
    with open(sys.argv[2], "w") as f:
        if "comment" in cp["settings"]:
            f.write("# {}\n\n".format(cp["settings"]["comment"]))
        m = 1
        for elem1 in sorted(elements):
            for elem2 in sorted(elements):
                for elem3 in sorted(elements):
                    gamma = params[frozenset([elem1, elem3])]["gamma"]
                    lambda3 = 2 * float(params[frozenset([elem1, elem3])]["mu"])
                    c = params[frozenset([elem1, elem3])]["c"]
                    d = params[frozenset([elem1, elem3])]["d"]
                    costheta0 = -float(params[frozenset([elem1, elem3])]["h"])
                    if elem2 == elem3:
                        n = 1
                        beta = 1
                        _beta = float(params[frozenset([elem1, elem2])]["beta"])
                        _S = float(params[frozenset([elem1, elem2])]["S"])
                        _D0 = float(params[frozenset([elem1, elem2])]["D0"])
                        _r0 = float(params[frozenset([elem1, elem2])]["r0"])
                        lambda2 = _beta * math.sqrt(2 / _S)
                        B = _S * _D0 / (_S - 1) * math.exp(lambda2 * _r0)
                        lambda1 = _beta * math.sqrt(2 * _S)
                        A = _D0 / (_S - 1) * math.exp(lambda1 * _r0)
                    else:
                        n = beta = lambda2 = B = lambda1 = A = 0
                    R = params[frozenset([elem1, elem3])]["Rcut"]
                    D = params[frozenset([elem1, elem3])]["Dcut"]
                    f.write(
                        " ".join(
                            str(i) for i in [
                                    elem1, elem2, elem3, m, gamma, lambda3,
                                    c, d, costheta0, n, beta, lambda2, B,
                                    R, D, lambda1, A
                            ]))
                    # Do we have ZBL parameters?
                    sec12 = params[frozenset([elem1, elem2])]
                    sec11 = params[frozenset([elem1, elem1])]
                    sec22 = params[frozenset([elem2, elem2])]
                    if "Z" in sec12:
                        if elem2 != elem3:
                            Zi = Zj = ZBLcut = ZBLexpscale = 0
                        else:
                            Zi = sec11["Z"]
                            Zj = sec22["Z"]
                            ZBLcut = sec12["ZBLcut"]
                            ZBLexpscale = sec12["ZBLexpscale"]
                        f.write(" ")
                        f.write(" ".join(str(i) for i in [Zi, Zj,
                                                          ZBLcut, ZBLexpscale]))
                    f.write("\n")
else:
    print("Unsupported style '{}'".format(cp["settings"]["style"]))
    sys.exit(1)
