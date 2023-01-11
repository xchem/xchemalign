import subprocess
import sys

from xchemalign import __version__


def test_cli_version():
    cmd = [sys.executable, "-m", "xchemalign", "--version"]
    assert subprocess.check_output(cmd).decode().strip() == __version__
