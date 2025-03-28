from __future__ import absolute_import, print_function, division

from .genome import Genome
import pysam
import logging
import re
from collections import OrderedDict

_logger = logging.getLogger("cnvpytor.fasta")


class Fasta:

    def __init__(self, filename):
        """
        Opens FASTA file, reads chromosome names/lengths and detects reference genome
        """
        self.reference_genome = None
        self.filename = filename
        self.file = None
        try:
            self.file = pysam.FastaFile(filename)
        except IOError:
            _logger.error(f"Problem opening file '{filename}'")
            exit(0)
        except ValueError:
            _logger.error(f"Index for filename '{filename}' is missing!")
            exit(0)

        self.len = {}
        if self.file:
            _logger.info(f"File: {filename} successfully open")
            self.reference_genome = Genome.detect_genome(self.file.references, self.file.lengths)
            if self.reference_genome:
                _logger.info(f"Detected reference genome: {self.reference_genome}")
            for c, l in zip(self.file.references, self.file.lengths):
                self.len[c] = l

    def get_chr_len(self):
        """Get chromosome names and lengths."""
        return self.file.references, self.file.lengths

    def _get_sequence(self, chr_name):
        if chr_name not in self.len:
            _logger.warning(f"Can not find chromosome '{chr_name}' in fasta file '{self.filename}'.")
            return None
        _logger.debug(f"Reading chromosome: {chr_name}")
        return self.file.fetch(chr_name).upper()

    def read_chromosome_gc(self, chr_name):
        """Reads chromosome GC/AT content in 100bp bins."""
        seq = self._get_sequence(chr_name)
        if seq is None:
            return None, None

        gc = [seq.count("G", i, i + 100) + seq.count("C", i, i + 100)
              for i in range(0, len(seq), 100)]
        at = [seq.count("A", i, i + 100) + seq.count("T", i, i + 100)
              for i in range(0, len(seq), 100)]

        n = self.len[chr_name] // 100 + 1
        if len(gc) < n:
            gc.append(0)
            at.append(0)

        tot = len(seq)
        sgc = sum(gc)
        sat = sum(at)
        snn = tot - sgc - sat
        _logger.info(f"GC/AT/N content: {100. * sgc / tot:.1f}% / {100. * sat / tot:.1f}% / {100. * snn / tot:.1f}%")

        return gc, at

    def read_chromosome_mask_p_regions(self, chr_name):
        """Reads chromosome strict mask P regions."""
        seq = self._get_sequence(chr_name)
        if seq is None:
            return None
        p = re.compile("([P]+)")
        return [i.span() for i in p.finditer(seq)]

    def print_reference_genome_template(self):
        """Prints to stdout template for configuration file for reference genome."""
        prefix = (
            'import_reference_genomes = {\n'
            '    "NAME": {\n'
            '        "name": "FULL NAME",\n'
            '        "species": "SPECIES NAME",\n'
            '        "chromosomes": OrderedDict(['
        )
        suffix = (
            '        ]),\n'
            '        "gc_file": "/..PATH../GC_FILE.pytor",\n'
            '        "mask_file": "/..PATH../MASK_FILE.pytor"\n'
            '    }\n'
            '}'
        )
        print(prefix)
        s = ""
        for ix, (c, l) in enumerate(zip(self.file.references, self.file.lengths), 1):
            if ix % 4 == 1:
                s += '            '
            s += f'("{c}", ({l},"A")), '
            if ix % 4 == 0:
                s += "\n"
        print(s[:-2], end="")
        print(suffix)
