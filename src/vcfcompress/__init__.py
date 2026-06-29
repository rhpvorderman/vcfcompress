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
import typing
from typing import Iterator, TextIO, Tuple

import xopen

DEFAULT_BLOCK_SIZE = 64 * 1024


class VCFFormatError(RuntimeError):
    pass

class VCFHeader(typing.NamedTuple):
    info_count: int
    format_count: int
    sample_count: int
    text: str

    @classmethod
    def from_file(cls, vcf_file: TextIO):
        first_line = vcf_file.readline()
        if not first_line.startswith("##fileformat=VCFv4"):
            raise VCFFormatError("Only VCFv4 formatted files are supported.")
        header = [first_line]
        format_count = 0
        info_count = 0
        sample_count = 0
        while True:
            line = vcf_file.readline()
            if line.startswith("##FORMAT"):
                format_count += 1
                continue
            if line.startswith("##INFO"):
                info_count += 1
                continue
            if line.startswith("#CHROM"):
                sample_count = len(line.split("\t")) - 8
                header.append(line)
                break
        else: # No break
            raise VCFFormatError("No sample line found")
        text = "".join(header)
        return cls(info_count, format_count, sample_count, text)


def split_vcf_line(line: str, info_count: int, format_count: int):
    chrom, pos, id, ref, alt, qual, filter, info, format, *samples = line.split("\t")
    info_fields = info.split(";")
    if len(info_fields) < info_count:
        info_fields.extend("" for _ in range(info_count - len(info_fields)))
    sample_columns = []
    for sample in samples:
        sample_fields = sample.split(":")
        if len(sample_fields) < format_count:
            sample_fields.extend("" for _ in range(format_count - len(sample_fields)))
        sample_columns.extend(sample_fields)
    ans=  (chrom, pos, id, ref, alt, qual, filter, *info_fields, format, *sample_columns)
    return ans

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


def read_vcf_block(vcf_file: TextIO, info_count: int, format_count: int,
                   block_size: int = 10_000) -> Iterator[Tuple[str, ...]]:
    for i in range(block_size):
        line = vcf_file.readline().strip()
        if line == "":
            sys.exit(0)
        yield split_vcf_line(line, info_count, format_count)

if __name__ == "__main__":
    with xopen.xopen(sys.argv[1], "rt", encoding="utf-8") as vcf_file:
        header = VCFHeader.from_file(vcf_file)
        while True:
            blocks = zip(*read_vcf_block(vcf_file, header.info_count, header.format_count))
            for block in blocks:
                print("\t".join(block))
