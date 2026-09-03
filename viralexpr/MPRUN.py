"""
Command execution utilities with logging and error handling.
"""

import subprocess
import sys
import logging

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('viralexpr')


def run_cmd(cmd, log=True, fail_exit=True):
    """Run a shell command with logging and error handling.

    Args:
        cmd: Command string or list to execute.
        log: Whether to log the command.
        fail_exit: Whether to exit on failure.

    Returns:
        subprocess.CompletedProcess instance.
    """
    if isinstance(cmd, list):
        cmd_str = ' '.join(cmd)
    else:
        cmd_str = cmd
        cmd = cmd_str.split()

    if log:
        logger.info('Running: {}'.format(cmd_str))

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error('Command failed: {}'.format(cmd_str))
        if result.stderr:
            logger.error('stderr: {}'.format(result.stderr.strip()))
        if fail_exit:
            sys.exit(result.returncode)

    return result
