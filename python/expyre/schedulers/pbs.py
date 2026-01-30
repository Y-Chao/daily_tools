#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__original__ = "ExPyRe"
__version__ = "1.0"

from .base import Scheduler


class PBS(Scheduler):
    """
    PBS scheduler class.

    Parameters
    ----------
    host : str
        The hostname of the remote scheduler.
    remsh_cmd : str | None
        The remote shell command. If None, defaults to environment variable RE_RSH or 'ssh'.

    """

    def __init__(self, host: str, remsh_cmd: str | None = None):
        super().__init__(host, remsh_cmd)
        self.hold_command = "qhold"
        self.release_command = "qrls"
        self.cancel_command = "qdel"

    def submit(
        self,
        ids: str,
        remote_dir: str,
        partition: str,
        commands: list[str],
        max_time: int,
        node_dict: dict,
        script_exec: str = "/bin/bash",
        pre_submit_cmds: list[str] = [],
        verbose: bool = False,
    ) -> str:
        """
        Submit a job on a remote machine using PBS scheduler.
        Parameters
        ----------
        ids: str
            The local job IDs to submit.
        remote_dir: str
            The remote directory where the job scripts are located.
        partition: str
            The partition to submit the job to.
        commands: list[str]
            The list of commands to execute in the job script.
        max_time: int
            The maximum time for the job in minutes.
        node_dict: dict
            The dictionary containing node information.
        script_exec: str
            The script executor, default is '/bin/bash'.
        pre_submit_cmds: list[str]
            The list of commands to execute before submission.
        verbose: bool
            Whether to print verbose output.

        Returns:
        --------
        remote_ids: str
            The remote job IDs after submission.
        """
        ...
