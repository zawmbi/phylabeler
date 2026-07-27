"""
archive_io.py - Transparent compressed and archived input handling.

Lets the rest of PhyLabeler accept files in whatever form they arrive in:

- plain text     names.dmp, tree.nwk
- gzip/bz2/xz    names.dmp.gz, tree.nwk.gz
- tar archives   taxdump.tar.gz, taxdump.tgz, taxdump.tar
- zip archives   taxdmp.zip (NCBI's other official packaging)

Formats are detected from magic bytes rather than file extensions, so an
NCBI taxdump.tar.gz still works after being renamed to .gz or saved with no
suffix at all. Archive members are read as streams; nothing is unpacked to
disk unless a caller explicitly asks for it via extract_matching().

Pure Python stdlib - no external dependencies.
"""

import bz2
import gzip
import io
import lzma
import os
import shutil
import tarfile
import zipfile

GZIP_MAGIC = b"\x1f\x8b"
BZ2_MAGIC = b"BZh"
XZ_MAGIC = b"\xfd7zXZ\x00"
ZIP_MAGIC = b"PK\x03\x04"
TAR_MAGIC = b"ustar"
TAR_MAGIC_OFFSET = 257

# Enough bytes to reach the tar magic in the header block.
_PROBE_SIZE = TAR_MAGIC_OFFSET + 8

_OPENERS = {
    "gzip": gzip.open,
    "bz2": bz2.open,
    "xz": lzma.open,
}

# Suffixes that only indicate compression, not content type.
COMPRESSION_SUFFIXES = (".gz", ".bz2", ".xz")


class ArchiveError(Exception):
    """Raised when an archive is unreadable or lacks an expected member."""


# ---- Format detection ----

def _compression_from_magic(head):
    if head.startswith(GZIP_MAGIC):
        return "gzip"
    if head.startswith(BZ2_MAGIC):
        return "bz2"
    if head.startswith(XZ_MAGIC):
        return "xz"
    return None


def _looks_like_tar(block):
    return block[TAR_MAGIC_OFFSET:TAR_MAGIC_OFFSET + 5] == TAR_MAGIC


def probe(path):
    """Inspect a file's magic bytes.

    Returns (compression, container) where compression is one of
    None/"gzip"/"bz2"/"xz" and container is one of None/"tar"/"zip".
    A taxdump.tar.gz reports ("gzip", "tar"); a names.dmp.gz reports
    ("gzip", None).
    """
    with open(path, "rb") as fh:
        head = fh.read(_PROBE_SIZE)

    if head.startswith(ZIP_MAGIC):
        return None, "zip"

    compression = _compression_from_magic(head)
    if compression is None:
        return None, "tar" if _looks_like_tar(head) else None

    # A compressed stream may still be a tar underneath.
    try:
        with _OPENERS[compression](path, "rb") as fh:
            inner = fh.read(_PROBE_SIZE)
    except (OSError, EOFError):
        return compression, None
    return compression, "tar" if _looks_like_tar(inner) else None


def is_archive(path):
    """True if the path is a tar or zip archive, compressed or not."""
    try:
        return probe(path)[1] is not None
    except OSError:
        return False


def strip_compression_suffix(name):
    """Drop a trailing .gz/.bz2/.xz so the real extension can be examined."""
    lower = name.lower()
    for suffix in COMPRESSION_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


# ---- Single files ----

def open_text(path, encoding="utf-8", errors="replace"):
    """Open a text file, transparently decompressing gzip/bz2/xz.

    Raises ArchiveError for tar/zip archives, which hold many files and must
    be read through stream_members() instead.
    """
    compression, container = probe(path)
    if container:
        raise ArchiveError(
            f"{path} is a {container} archive, not a single text file"
        )
    if compression:
        return _OPENERS[compression](path, "rt", encoding=encoding,
                                     errors=errors)
    return open(path, "r", encoding=encoding, errors=errors)


# ---- Archive members ----

def stream_raw_members(path, wanted=None):
    """Yield (basename, binary_stream) for each file member of an archive.

    Members arrive in archive order and each stream is only valid until the
    next iteration, so read it before advancing. `wanted` optionally limits
    the walk to a collection of basenames.
    """
    _, container = probe(path)
    wanted = set(wanted) if wanted is not None else None

    if container == "tar":
        with tarfile.open(path, "r:*") as tar:
            for member in tar:
                if not member.isfile():
                    continue
                # basename only: archive paths never steer where we read or
                # write, so a hostile "../../etc/passwd" member is inert.
                base = os.path.basename(member.name)
                if wanted is not None and base not in wanted:
                    continue
                stream = tar.extractfile(member)
                if stream is None:
                    continue
                with stream:
                    yield base, stream

    elif container == "zip":
        with zipfile.ZipFile(path) as zf:
            for info in zf.infolist():
                if info.is_dir():
                    continue
                base = os.path.basename(info.filename)
                if wanted is not None and base not in wanted:
                    continue
                with zf.open(info) as stream:
                    yield base, stream

    else:
        raise ArchiveError(f"{path} is not a tar or zip archive")


def stream_members(path, wanted=None, encoding="utf-8", errors="replace"):
    """Yield (basename, text_stream) for each file member of an archive.

    Same streaming contract as stream_raw_members(): read the stream before
    advancing the iterator.
    """
    for base, raw in stream_raw_members(path, wanted):
        with io.TextIOWrapper(raw, encoding=encoding, errors=errors) as text:
            yield base, text


def list_members(path):
    """Return the file member names inside a tar or zip archive."""
    _, container = probe(path)
    if container == "tar":
        with tarfile.open(path, "r:*") as tar:
            return [m.name for m in tar if m.isfile()]
    if container == "zip":
        with zipfile.ZipFile(path) as zf:
            return [i.filename for i in zf.infolist() if not i.is_dir()]
    raise ArchiveError(f"{path} is not a tar or zip archive")


def extract_matching(path, dest_dir, match):
    """Extract archive members whose basename satisfies match(basename).

    Members are written flat into dest_dir under their basename, so nested or
    hostile archive paths cannot escape it. Returns the extracted file paths.
    """
    os.makedirs(dest_dir, exist_ok=True)
    written = []
    for base, raw in stream_raw_members(path):
        if not match(base):
            continue
        out_path = os.path.join(dest_dir, base)
        with open(out_path, "wb") as out:
            shutil.copyfileobj(raw, out)
        written.append(out_path)
    return written
