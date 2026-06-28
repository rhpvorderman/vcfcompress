# Copyright (C) 2023 Leiden University Medical Center
# This file is part of vcfcompress
#
# vcfcompress is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# vcfcompress is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with vcfcompress.  If not, see <https://www.gnu.org/licenses/
import io
import sys
from typing import TextIO

import xopen

DEFAULT_BLOCK_SIZE = 64 * 1024


class VCFFormatError(RuntimeError):
    pass


def vcf_read_header(vcf_file: TextIO):
    first_line = vcf_file.readline()
    if not first_line.startswith("##fileformat=VCFv4"):
        raise VCFFormatError("Only VCFv4 formatted files are supported.")
    header = [first_line]
    while True:
        line = vcf_file.readline()
        if line.startswith("#CHROM"):
            header.append(line)
            return "".join(header)


def chunk_vcf(vcf_file: TextIO, block_size = DEFAULT_BLOCK_SIZE):
    remainder = ""
    while True:
        block = vcf_file.read(block_size)
        if not (block or remainder):
            return
        last_newline_pos = block.rfind("\n")
        if last_newline_pos == -1:
            raise OverflowError(
                f"VCF record larger than block size ({block_size}): "
                f"{(remainder + block)[:100]}[...]")
        chunk = remainder + block[:last_newline_pos + 1]
        remainder = block[last_newline_pos + 1:]
        yield chunk


def transpose_chunk(vcf_chunk: str) -> str:
    result = io.StringIO()
    records = vcf_chunk.splitlines(keepends=False)
    result.write(f"##Block;records={len(records)}")
    transposed = zip(*(line.split("\t") for line in records))
    print(list(transposed))


if __name__ == "__main__":
    with xopen.xopen(sys.argv[1], "rt", encoding="utf-8") as vcf_file:
        header = vcf_read_header(vcf_file)
        print(header)
        for chunk in chunk_vcf(vcf_file):
            transpose_chunk(chunk)
            break

