#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__original__ = "ExPyRe"
__version__ = "1.0"

import os
import re
import shlex
import subprocess
import time
import warnings
from pathlib import Path


def _optionally_remote_args(
        args: list[str],
        shell: str,
        host:  str|None,
        remsh_cmd: str | list[str],
        in_dir: str = "_HOME_",
    ) -> list[str]:
    """
    Prepare command arguments for local or remote execution.
    Parameters
    ----------
    args : list[str]
        The list of arguments/commands to execute.
    shell : str
        The shell to use for executing the commands.
    host : str | None
        The hostname of the remote machine. If None, execute locally.
    remsh_cmd : str | list[str]
        The remote shell command.
    in_dir : str 
        The directory to change into before executing the command. Default, _HOME_.

    Returns
    -------
    list[str]
        The prepared list of command arguments.
    """

    if isinstance(remsh_cmd, str):
        remsh_cmd = shlex.split(remsh_cmd)

    if in_dir == "_PWD_":
        assert host is None, "in_dir='_PWD_' is only supported for local execution."
    
    # Parse args into a single command string
    # blackslash escape single quotes in args
    safe_cmd = ""
    if in_dir == "_HOME_":
        safe_cmd = 'cd $HOME'
    elif in_dir != "_PWD_":
        safe_cmd = f'cd {shlex.quote(in_dir)}'
    
    if len(args) > 0:
        safe_cmd += " && "
        safe_cmd += ' '.join([re.sub(r'([\'" \(\)])', r'\\\1', arg) for arg in args])
    
    args = shell.split() + [safe_cmd]
    if host is not None:
        args = remsh_cmd + [host] + args[0:-1] + ["'" + args[-1] + "'"]
    return args


def subprocess_run(
    host: str|None,
    args: list[str],
    script: str | None = None,
    shell: str = "bash -c",
    remsh_cmd: str | None = None,
    retry: tuple[int, int] | None = None,
    in_dir: str | None = None,
    dry_run: bool = False,
    verbose: bool = False,
) -> tuple[str, str]:
    """
    Run a subprocess, optionally via ssh on a remote machine.

    Parameters
    ----------
    host : str
        The hostname of the remote machine.
    args : list[str]
        The list of arguments/commands to execute.
    script : str | None
        The script content to execute. If provided, it will be written to a temporary file and
    shell: str
        The shell to use for executing the commands, default is 'bash -c'.
    remsh_cmd : str | None
        The remote shell command. If None, defaults to environment variable RE_RSH or 'ssh'.
    retry : tuple[int, int] | None
        The number of times to retry the command on failure. If None, no retries are performed.
    in_dir : str | None
        The directory to change into before executing the command. If None, uses the current directory.
    dry_run : bool
        If True, only print the command without executing it.
    verbose : bool
        If True, print verbose output.

    Returns
    -------
    stdout: bytes.decode()
    stderr: bytes.decode()
        The standard output and standard error of the executed command.
    """

    if remsh_cmd is None:
        remsh_cmd = os.environ.get("RE_RSH", "ssh")
    if retry is None:
        retry = tuple([int(ii) for ii in os.environ["RE_RETRY"].strip().split("")])
    else:
        retry = (3, 5)

    # ensure at least 1 retry with 0 delay
    retry = tuple(max(retry[0], 1), max(retry[1], 0))

    args = _optionally_remote_args(args, shell, host, remsh_cmd, in_dir)

    if verbose:
        if dry_run:
            print("[DRY-RUN]: ")
        else:
            print("[RUN]: ")
        print(f"{" ".join([shlex.quote(arg) for arg in args])}")

        if script is not None:
            print("[SCRIPT]: ")
            print(script.rstrip())

    if script is not None:
        script = script.encode()

    if dry_run:
        return args, script
    
    for i_try in range(retry[0]):
        try:
            p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE,
                                 close_fds=False, env=os.environ)
            stdout, stderr = p.communicate(script)
            if p.returncode != 0:
                raise RuntimeError(f'Failed to run command: {" ".join(args)} with error {stderr.decode()}\n')
            
            if i_try > 0:
                warnings.warn(f'Succeeded to run "{" ".join(args)}" on attempt {i_try} after failure(s), trying again')
            break
        except Exception:
            if i_try == retry[0] - 1:
                warnings.warn(f'Failed to run "{" ".join(args)}" on attempt {i_try}.\n STDERR: {stderr.decode()}\n')
                raise RuntimeError(f'Failed to run command: {" ".join(args)} after {retry[0]} attempts\n')
            warnings.warn(f'Failed to run "{" ".join(args)}" on attempt {i_try}, trying again.\n STDERR: {stderr.decode()}\n')
            time.sleep(retry[1])

    if verbose:
        print("[STDOUT]: ")
        print(stdout.decode().rstrip())
        print("[STDERR]: ")
        print(stderr.decode().rstrip())

    return stdout.decode(), stderr.decode()


def subprocess_copy(
        from_path: str,
        to_path: str,
        from_host: str = "_LOCAL_",
        to_host:str = "_LOCAL_",
        rcp_cmd: str = "rsync",
        rcp_args: str="-a",
        remsh_cmd: str | None = None,
        remsh_flags: str = "-e",
        retry: tuple[int, int] = False,
        delete: bool= False,
        dry_run: bool = False,
        verbose: bool = False
        ):
    """
    Run a remote copy ("e.g. rsync or scp") in a subprocess.

    Parameters
    ----------
    from_path : str
        The source path to copy from.
    to_path : str
        The destination path to copy to.
    from_host : str
        The hostname of the source machine. Use "_LOCAL_" for local machine.
    to_host : str
        The hostname of the destination machine. Use "_LOCAL_" for local machine.
    rcp_cmd : str
        The remote copy command, e.g., 'rsync' or 'scp'.
    rcp_args : str
        The arguments for the remote copy command.
    rmesh_cmd : str | None
        The remote shell command. If None, defaults to environment variable RE_RSH or 'ssh'.
    rmesh_flags: str | None
        The flags to pass to the remote shell command. 'rsync -e "ssh -i keyfile"'
    retry : tuple[int, int]
        The number of times to retry the command on failure.
    delete : bool
        Whether to delete source file after copy.
    dry_run : bool
        If True, only print the command without executing it.
    verbose : bool
        If True, print verbose output.
    """

    # avoide copy on remote to remote
    if from_host != "_LOCAL_" and to_host != "_LOCAL_":
        raise RuntimeError("Remote to remote copy is not supported.")
    
    if remsh_cmd is None:
        remsh_cmd = os.environ.get("RE_RSH", "ssh")

    rcp_args = remsh_flags + " " + remsh_cmd + " " + rcp_args

    if delete:
        rcp_args += " --delete"

    if isinstance(from_path, str) or isinstance(from_path, Path):
        from_path = [from_path]
    
    if from_host is None:
        abs_from_files = []
        for f in from_path:
            if not Path(f).is_absolute():
                abs_from_files.append(str(Path.cwd() / f))
            abs_from_files.append(str(f))
    else:
        abs_from_files = from_path

    if to_host is None and not Path(to_path).is_absolute():
        abs_to_files = Path.home() / to_path
    else:
        abs_to_files = to_path

    if from_host is None or from_host == "_LOCAL_":
        from_host = ""

    if to_host is None or to_host == "_LOCAL_":
        to_host = ""

    if len(from_host) > 0:
        from_host += ":"
    if len(to_host) > 0:
        to_host += ":"

    # Prepend host parts to from_path and to_path
    abs_from_files = [from_host + str(f) for f in abs_from_files]
    abs_to_files = to_host + str(abs_to_files)

    # do copy
    retval = subprocess_run(
        host = None,
        args = [rcp_cmd] + rcp_args.split() + abs_from_files + [abs_to_files],
        retry=retry,
        in_dir='_PWD_',
        dry_run=dry_run,
        verbose=verbose
    )

    if dry_run:
        return retval