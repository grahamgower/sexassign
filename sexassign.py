#!/usr/bin/env python
# Copyright (c) 2016-2018 Graham Gower <graham.gower@gmail.com>
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

from __future__ import print_function
import os
import sys
import operator
import csv
from distutils.version import LooseVersion

import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.linalg as linalg
from scipy.stats import beta, binom, chi2, norm

def parse_idxstats(filename, min_length, exclude_contigs, include_contigs):
    N = {}
    L = {}
    with open(filename) as f:
        for line in f:
            fields = line.split()
            cid = fields[0]
            clen = int(fields[1])
            nreads = int(fields[2])
            if clen < min_length or cid in exclude_contigs:
                continue
            if include_contigs and cid not in include_contigs:
                continue
            N[cid] = nreads
            L[cid] = clen
    return N, L

def get_sex(sample, Nx, Na, Lx, La):
    Rx = float(Nx)/(Nx+Na)

    # Beta CI with non-informative prior, aka Jefferey's interval.
    # See Brown, Cai, and DasGupta (2001). doi:10.1214/ss/1009213286
    Rx_CI = beta.interval(0.99, Nx+0.5, Na+0.5)

    # expected ratios from the chromosome lengths
    Elx_X0 = float(Lx)/(Lx+2*La)
    Elx_XX = float(Lx)/(Lx+La)

    #ll_x0 = beta.logpdf(Elx_X0, Nx+0.5, Na+0.5)
    #ll_xx = beta.logpdf(Elx_XX, Nx+0.5, Na+0.5)
    ll_x0 = binom.logpmf(Nx, Nx+Na, Elx_X0)
    ll_xx = binom.logpmf(Nx, Nx+Na, Elx_XX)

    # likelihood ratio test
    alpha = 0.001
    if chi2.sf(2*(ll_x0-ll_xx), 1) < alpha:
        sex = 'M'
    elif chi2.sf(2*(ll_xx-ll_x0), 1) < alpha:
        sex = 'F'
    else:
        # indeterminate
        sex = 'U'

    if ll_x0 > ll_xx:
        Elx = 2*Elx_X0
    else:
        Elx = Elx_XX

    Mx = Rx/Elx
    Mx_CI = [Rx_CI[0]/Elx, Rx_CI[1]/Elx]

    if Mx < 0.4 or Mx > 1.2:
        #print("Warning: {} has unexpected Mx={:g}".format(sample, Mx), file=sys.stderr)
        pass

    if Mx > 0.6 and Mx < 0.8:
        # suspicious sample, may be contaminated
        sex = 'U'

    return Elx, Mx, Mx_CI, sex

def parse_multi(fn_list, chrX, min_length, min_reads, exclude_contigs, include_contigs):
    xset = set(chrX)
    cset = None
    data = []
    M = []

    for fn in fn_list:
        sample = os.path.basename(fn)
        if sample.endswith(".idxstats"):
            sample = sample[:-len(".idxstats")]

        N, L = parse_idxstats(fn, min_length, exclude_contigs, include_contigs)

        if cset is None:
            cset = set(N.keys())
            xset &= cset
            if len(xset) == 0:
                raise Exception("{}: cannot find X chromosome".format(fn))
        else:
            if cset & set(N.keys()) != cset:
                raise Exception("{}: reference mismatch".format(fn))

        Nx = Lx = 0
        for ch in xset:
            Nx += N[ch]
            Lx += L[ch]
        Nt = sum(N.values())
        Na = Nt-Nx
        Lt = sum(L.values())
        La = Lt-Lx


        if Nx+Na < min_reads:
            #print("{} has {} reads, ignoring".format(sample, Nx+Na), file=sys.stderr)
            continue

        Elx, Mx, Mx_CI, sex = get_sex(sample, Nx, Na, Lx, La)
        M.append([(float(n)/Nt)/(float(L[s])/Lt) for s,n in N.iteritems()])

        datum = (sample, Nx, Na, Lx, La, Mx, Mx_CI, sex, Elx)
        data.append(datum)

    return data, np.array(M)

def parse_pecnerova_csv(fn):
    """
    Pecnerova et al. 2017, mammoth data.
    """
    data = []
    with open(fn) as f:
        for row in csv.DictReader(f):
            sample = row["ID"]
            if not sample:
                continue
            #age = row["Age (calBP)"]
            #try:
            #    age = float(int(age))
            #except ValueError:
            #    age = None
            #tissue = row["Material"]
            #location = row["Locality"]
            Nx = int(row["chrX"])
            Na = int(row["Total"]) - Nx

            # length of chrX, length of autosome, from LoxAfr4
            Lx, La = 120050768, 3050898245

            Elx, Mx, Mx_CI, sex = get_sex(sample, Nx, Na, Lx, La)
            datum = (sample, Nx, Na, Lx, La, Mx, Mx_CI, sex, Elx)
            data.append(datum)
    return data

def PCA_eig(data, proj_data=None):

    X = data - np.mean(data, axis=0)
    C = np.cov(X, rowvar=0, ddof=1)
    evals, evecs = linalg.eigh(C)

    # sort eigenvalues in decreasing order
    idx = np.argsort(np.abs(evals))[::-1]
    evecs = evecs[:,idx]
    # sort eigenvectors according to same index
    evals = evals[idx]

    # transform data using eigenvectors
    if proj_data is None:
        proj_data = data
    trans_data = np.dot(proj_data, evecs)

    return trans_data, evecs, evals

def PCA_svd(data, proj_data=None):

    X = data - np.mean(data, axis=0)

    U, s, Vh = linalg.svd(X)
    evals = s # WTF? shouldn't this be s*s/(data.shape[0]-1)?
    evecs = Vh.T

    # transform data using eigenvectors
    if proj_data is None:
        proj_data = data
    trans_data = np.dot(proj_data, evecs)

    return trans_data, evecs, evals

def plot_PCA(ax, data, M):
    # calculate eigenv{al,ec}s
    #trans_data, evecs, evals = PCA_svd(M)
    trans_data, evecs, evals = PCA_eig(M)
    # principal components
    pc1 = trans_data[:,0]
    pc2 = trans_data[:,1]
    pc3 = trans_data[:,2]
    pc4 = trans_data[:,3]
    colours = {"M":"red", "F":"blue", "U":"gray"}
    markers = {"M":"o", "F":"^", "U":"x"}
    for i, (s,_,_,_,_,_,_,sex,_) in enumerate(data):
        col = colours[sex]
        m = markers[sex]
        ax.scatter(pc1[i], pc2[i], marker=m, facecolor=col, alpha=0.7, lw=0, s=60, label=s)
        #ax.scatter(pc3[i], pc4[i], marker=m, facecolor=col, alpha=0.7, lw=0, s=60, label=s)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")

    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    #sub_width = 0.4
    #ax2 = fig1.add_axes([1-sub_width, 0, sub_width, 0.3])
    ax2 = inset_axes(ax, width="30%", height="20%")#, loc=4)

    nscree = min(20, len(evals))
    bar_width = 0.8
    scree_x = np.arange(nscree)
    scree_y = 100.0 * evals[:nscree] / sum(evals)
    ax2.bar(scree_x+0.5 -bar_width/2.0, scree_y, width=bar_width)
    #ax2.set_xticks(scree_x+0.5)
    #ax2.set_xticklabels([x+1 for x in scree_x])
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])

def plot_samples(ax, data, colour="black"):

    samples = map(operator.itemgetter(0), data)
    Nx = np.array(map(operator.itemgetter(1), data), dtype=int)
    Na = np.array(map(operator.itemgetter(2), data), dtype=int)
    Mx = np.array(map(operator.itemgetter(5), data), dtype=float)
    Mx_CI = map(operator.itemgetter(6), data)
    sex = map(operator.itemgetter(7), data)
    Elx = np.array(map(operator.itemgetter(8), data), dtype=float)

    Rx_m = [x for x,sx in zip(Mx,sex) if sx=='M']
    Rx_f = [x for x,sx in zip(Mx,sex) if sx=='F']

    ax.vlines(0.5, -2, len(samples), linestyle=':')
    ax.vlines(1.0, -2, len(samples), linestyle=':')

    y_pos = np.arange(len(samples))
    ax.set_yticks(y_pos)
    ax.set_yticklabels(["{} ({})".format(s,x+a) for s,x,a in zip(samples,Nx,Na)])


    if len(Rx_m) > 1:
        ax.vlines(np.mean(Rx_m), -2, len(samples), linestyle='-', color='red')
        #ax.vlines(2*np.mean(Rx_m), -2, len(samples), linestyle='-.', color='blue')
        m_ci = norm.interval(0.99, np.mean(Rx_m), np.std(Rx_m))
        ax.fill_between(m_ci, -2, len(samples), alpha=0.2, color="red", edgecolor="none")
    else:
        ax.fill_between([0.4, 0.6], -2, len(samples), alpha=0.2, color="red", edgecolor="none")

    if len(Rx_f) > 1:
        ax.vlines(np.mean(Rx_f), -2, len(samples), linestyle='-', color='blue')
        f_ci = norm.interval(0.99, np.mean(Rx_f), np.std(Rx_f))
        ax.fill_between(f_ci, -2, len(samples), alpha=0.2, color="blue", edgecolor="none")
    else:
        ax.fill_between([0.8, 1.2], -2, len(samples), alpha=0.2, color="blue", edgecolor="none")

    sex_colour = {'M': "red", 'F': "blue", 'U': 'black'}
    ecol = [sex_colour[sx] for sx in sex]

    ax.scatter(Mx, y_pos, facecolor=ecol, edgecolors=colour, lw=0.5, s=60)

    err_low = Mx - np.array(map(operator.itemgetter(0), Mx_CI))
    err_high = np.array(map(operator.itemgetter(1), Mx_CI)) - Mx
    ax.errorbar(Mx, y_pos, xerr=[err_low, err_high], ecolor=colour, marker="none", fmt="none", capsize=0)

    ax.set_ylim([-0.5, len(samples)-0.5])
    ax.set_xlim([0, 1.5])

    ax.set_xlabel('Read dosage (X)', size=16)
    ax.set_ylabel('Sample (number of sequences)', size=16)


def parse_list(fn):
    l = []
    with open(fn) as f:
        for line in f:
            line = line.rstrip()
            if not line or line[0] == "#":
                continue
            l.append(line)
    return l

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="sex determination from read dosages")
    parser.add_argument("-w", "--wide", action="store_true", default=False, help="plot widescreen ratio (16x9) [%(default)s]")
    parser.add_argument("-o", "--opdf", type=str, default="out.pdf", help="output filename [%(default)s]")
    parser.add_argument("-n", "--min-reads", type=int, default=5000, help="exclude samples with fewer reads than this [%(default)s]")
    parser.add_argument("-l", "--min-length", type=int, default=1000*1000, help="exclude contigs shorter than this [%(default)s]")
    parser.add_argument("-X", "--chrX-list", default=None, help="""file containing list of X linked contigs,
                                                                if not specified, defaults to X, chrX, ChrX""")
    parser.add_argument("-e", "--exclude-list", default=None, help="""file containing list of contigs to exclude,
                                                                    if not specified, defaults to Y, chrY, ChrY, M, MT, Mt""")
    parser.add_argument("-i", "--include-list", default=None, help="""file containing list of contigs to include,
                                                                    if not specified, defaults to all contigs""")
    parser.add_argument("-p", "--pca", action="store_true", default=False, help="Plot PCA of read dosages")
    parser.add_argument("infiles", nargs="+", help="input file(s)")
    args = parser.parse_args()

    if args.min_length < 1:
        args.min_length = 1

    if args.chrX_list is not None:
        args.chrX = parse_list(args.chrX_list)
    else:
        args.chrX = ["X","chrX","ChrX"]

    if args.exclude_list is not None:
        args.exclude_contigs = set(parse_list(args.exclude_list))
    else:
        args.exclude_contigs = set(["Y","chrY","ChrY","M","MT","Mt"])

    if args.include_list:
        args.include_contigs = set(args.include_contigs.split(","))
    else:
        args.include_contigs = None

    return args

def get_id(row):
    sfields = row[0].split("_")
    s = sfields[0]
    if s.startswith("A") or s.startswith("M"):
        s = s[1:]
    if s.startswith("GilbertM"):
        s = s[len("GilbertM"):]
    return s

if __name__ == "__main__":
    args = parse_args()

    if len(args.infiles) == 1 and args.infiles[0].endswith(".csv"):
        data = parse_pecnerova_csv(args.infiles[0])
        M = None
    else:
        data, M = parse_multi(args.infiles, args.chrX, args.min_length, args.min_reads, args.exclude_contigs, args.include_contigs)

    data.sort(key=lambda x: LooseVersion(get_id(x)), reverse=True)

    print("sample", "Mx", "sex", "Nx", "Na", "Lx", "La", sep="\t")
    for row in sorted(data, key=get_id):
        (sample, Nx, Na, Lx, La, Mx, Mx_CI, sex, Elx) = row
        print(sample, Mx, sex, Nx, Na, Lx, La, sep="\t")

    if len(data) > 25:
        height = len(data)/15
    else:
        height = 1
    width = 2

    pdf = PdfPages(args.opdf)
    if args.wide:
        fig_w, fig_h = plt.figaspect(9.0/16.0)
    else:
        fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(width*fig_w, height*fig_h))
    gs1 = gridspec.GridSpec(1+height, 1)

    ax1 = fig1.add_subplot(gs1[0:])
    plot_samples(ax1, data)

    plt.tight_layout()
    pdf.savefig(figure=fig1)

    if M is None and args.pca:
        print("Cannot plot PCA, no M matrix", file=sys.stderr)
    elif args.pca:
        fig2 = plt.figure(figsize=(fig_w, fig_h))
        gs2 = gridspec.GridSpec(1, 1)
        ax2 = fig2.add_subplot(gs2[0])

        plot_PCA(ax2, data, M)

        plt.tight_layout()
        pdf.savefig(figure=fig2)

    pdf.close()

