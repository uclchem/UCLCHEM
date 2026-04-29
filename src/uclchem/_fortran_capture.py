"""Capture Fortran stdout/stderr output and route through logging/files.

F2PY-compiled Fortran code writes directly to file descriptors 1 (stdout)
and 2 (stderr), bypassing Python's sys.stdout/sys.stderr. This module
provides a context manager that redirects those file descriptors to a pipe,
reads the output in a background thread, and streams it to the terminal
and optionally to a per-model log file on disk in real time.
"""

import os
import sys
import threading
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path


def _reader_thread(
    read_fd: int, saved_stdout_fd: int, prefix: str, log_file: str | Path | None
):
    """Read from pipe and stream each line to terminal and log.

    Uses raw os.read() to bypass all Python IO buffering.
    The log file is lazily created on first non-empty line.

    Parameters
    ----------
    read_fd : int
        file descriptor of where fortran puts output
    saved_stdout_fd : int
        file descriptor of where to log output
    prefix : str
        prefix of what to add in front of the log
    log_file : str | Path | None
        path to write logs to.
        If None, do not write to file, only to stdout.
    """
    log_handle = None
    buf = b""
    prefix_bytes = prefix.encode()
    try:
        while True:
            chunk = os.read(read_fd, 4096)
            if not chunk:
                break
            buf += chunk
            while b"\n" in buf:
                line_bytes, buf = buf.split(b"\n", 1)
                line = line_bytes.decode("utf-8", errors="replace")
                if not line:
                    continue
                os.write(
                    saved_stdout_fd,
                    prefix_bytes + line_bytes + b"\n",
                )
                if log_file is not None:
                    if log_handle is None:
                        log_path = Path(log_file)
                        log_path.parent.mkdir(
                            parents=True,
                            exist_ok=True,
                        )
                        log_handle = log_path.open("w")
                    log_handle.write(f"{line}\n")
                    log_handle.flush()
        # Flush any trailing data without a final newline
        if buf:
            line = buf.decode("utf-8", errors="replace")
            if line:
                os.write(
                    saved_stdout_fd,
                    prefix_bytes + buf + b"\n",
                )
                if log_handle is not None:
                    log_handle.write(f"{line}\n")
    except Exception:  # noqa: S110
        pass
    finally:
        os.close(read_fd)
        if log_handle is not None:
            log_handle.close()


@contextmanager
def capture_fortran_output(
    label: str = "", log_file: str | Path | None = None
) -> Iterator[None]:
    """Capture Fortran output, stream to terminal and file.

    Redirects OS-level file descriptors 1 and 2 to a pipe.
    A background reader thread prints each line to the terminal
    in real time and optionally writes to a per-model log file.

    Parameters
    ----------
    label : str
        Identifier prepended to terminal lines. Default = "".
    log_file : str | Path | None
        Per-model log file path.
        Only created if there is at least one line of output. Default = None.

    Yields
    ------
    None
        nothing

    Examples
    --------
    >>> from uclchem.model import Cloud
    >>>
    >>> with capture_fortran_output(
    ...     label="model_3",
    ...     log_file="logs/model_3.log",
    ... ):
    ...    result = Cloud({})
    """
    # Flush Python buffers before redirecting
    sys.stdout.flush()
    sys.stderr.flush()

    # Save original file descriptors
    saved_stdout_fd = os.dup(1)
    saved_stderr_fd = os.dup(2)

    # Create pipe: Fortran writes to write_fd, reader thread reads from read_fd
    read_fd, write_fd = os.pipe()

    prefix = f"[{label}] " if label else ""

    try:
        # Redirect stdout and stderr to the pipe's write end
        os.dup2(write_fd, 1)
        os.dup2(write_fd, 2)
        os.close(write_fd)

        # Start reader thread — streams lines as they arrive
        reader = threading.Thread(
            target=_reader_thread,
            args=(read_fd, saved_stdout_fd, prefix, log_file),
            daemon=True,
        )
        reader.start()

        yield

    finally:
        # Flush any C-level buffers on the redirected fds
        sys.stdout.flush()
        sys.stderr.flush()

        # Restore original file descriptors
        os.dup2(saved_stdout_fd, 1)
        os.dup2(saved_stderr_fd, 2)
        os.close(saved_stderr_fd)

        # write_fd was closed above; fd 1 & 2 now point to originals,
        # so the pipe write end has no open references and reader sees EOF.
        reader.join(timeout=5.0)

        # Close saved stdout fd after reader thread is done with it
        os.close(saved_stdout_fd)
