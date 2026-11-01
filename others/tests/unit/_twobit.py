import struct

_CODE = {"T": 0, "C": 1, "A": 2, "G": 3}


def write_twobit(path, sequences):
    """Write a minimal UCSC 2bit file (A/C/G/T only, no N blocks, no soft-mask blocks)."""
    names = list(sequences)
    header_size = 16
    index_size = sum(1 + len(name.encode()) + 4 for name in names)
    offset = header_size + index_size
    index = b""
    body = b""
    for name in names:
        seq = sequences[name].upper()
        if set(seq) - set(_CODE):
            raise ValueError("write_twobit supports A/C/G/T only")
        packed = bytearray((len(seq) + 3) // 4)
        for i, base in enumerate(seq):
            packed[i // 4] |= _CODE[base] << (6 - 2 * (i % 4))
        # dnaSize, nBlockCount, maskBlockCount, reserved (zero blocks, so no start/size arrays)
        record = struct.pack("<IIII", len(seq), 0, 0, 0) + bytes(packed)
        index += struct.pack("B", len(name.encode())) + name.encode() + struct.pack("<I", offset)
        offset += len(record)
        body += record
    header = struct.pack("<IIII", 0x1A412743, 0, len(names), 0)
    with open(path, "wb") as handle:
        handle.write(header + index + body)
    return path
