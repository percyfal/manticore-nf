#!/usr/bin/env python
r"""Partition BEDFILE into equisized bins.

Partition BEDFILE regions into bins (subsets) of regions. Partitioning
is either done sequentially or greedily with the aim of generating
equisized partitions.

Sequential partitioning bins regions based on their index order.

Greedy partitioning first sorts regions by size and then iterates
through the remaining regions, adding to the smallest bin. Greedy
partitioning may be preferred when sequences are of very different
sizes.

"""
import argparse
import logging
import os
import sys
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


class Interval:
    """0-based bed interval"""

    def __init__(self, chrom, begin, end):
        self.chrom = chrom
        self.begin = begin
        self.end = end

    def __len__(self):
        return self.end - self.begin

    def __repr__(self):
        return f"{self.chrom}_{self.begin}-{self.end}"

    def __str__(self):
        return f"{self.chrom}\t{self.begin}\t{self.end}"


def bin_length(x):
    """Calculate total length of bin"""
    return sum(len(r) for r in x)


def cumsum(x):
    """Calculate total length of list of bins"""
    return sum(bin_length(b) for b in x)


def greedy_partition(intervals, npartitions):
    # Sort intervals by length
    ix = sorted(range(len(intervals)), key=lambda k: len(intervals[k]), reverse=True)
    # Seed output intervals with npartitions longest intervals
    out = [[intervals[i]] for i in ix[0:npartitions]]
    # Keep track of output lengths
    outlen = [len(r[0]) for r in out]
    for j in ix[npartitions : len(ix)]:  # noqa: E203
        # Get index of output set with shortest total region length
        # and add current region
        imin = outlen.index(min(outlen))
        out[imin].append(intervals[j])
        outlen[imin] += len(intervals[j])
    return out


def sequential_partition(intervals, npartitions):
    # Get total length
    length = bin_length(intervals)
    binsize = length / npartitions
    out = []
    rbin = []
    size = 0
    nremaining = len(intervals)
    for r in intervals:
        if (size >= binsize) or (nremaining < (npartitions - len(out))):
            out.append(rbin)
            # Adjust binsize
            if npartitions > len(out):
                binsize = (length - cumsum(out)) / (npartitions - len(out))
            rbin = []
        rbin.append(r)
        size = bin_length(rbin)
        nremaining = nremaining - 1

    out.append(rbin)
    return out


def parse_args(argv):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Partition BED_FILE into equisized bins.",
        epilog="Example: python partition_bed.py intervals.bed",
    )
    parser.add_argument(
        "bed_file",
        metavar="BED_FILE",
        type=Path,
        help="Input interval file in bed format",
    )
    parser.add_argument(
        "--output-prefix",
        metavar="OUTPUT_PREFIX",
        type=Path,
        default=Path(os.curdir),
        help="Output prefix",
    )
    parser.add_argument(
        "--npartitions",
        "-n",
        default=20,
        type=int,
        help="number of requested partitions",
    )
    parser.add_argument(
        "--greedy", "-g", action="store_true", help="do greedy partition"
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.bed_file.is_file():
        logger.error(f"The given input bed file {args.bed_file} was not found!")
        sys.exit(2)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    data = pd.read_table(args.bed_file, header=None)
    intervals = [Interval(x[0], x[1], x[2]) for x in data.to_dict("records")]

    try:
        assert len(intervals) >= args.npartitions
    except AssertionError:
        logger.warning(
            "Number of intervals smaller than number of partitions: "
            f"'{len(intervals)} < {args.npartitions}': "
            "lower the number of partitions with --npartitions"
        )
        raise

    func = sequential_partition
    if args.greedy:
        func = greedy_partition

    out = func(intervals, args.npartitions)

    assert (
        len(out) == args.npartitions
    ), "total number of output partitions less than requested"
    assert cumsum(out) == bin_length(intervals), (
        "total bin length not equal to input bed length"
        f"{cumsum(out)} != {bin_length(intervals)}"
    )
    for partition in out:
        fn = args.output_prefix / f"{repr(partition[0])}.bed"
        fn.write_text("".join([f"{x}\n" for x in partition]))


if __name__ == "__main__":
    sys.exit(main())
