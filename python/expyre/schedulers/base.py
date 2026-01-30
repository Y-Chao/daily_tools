#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__original__ = "ExPyRe"
__version__ = "1.0"

import os

from ..subprocess import subprocess_run


class Scheduler:
    """
    Base class for all schedulers, which define what methods a scheduler should implement.
    """

    def __init__(self, host: str, remsh_cmd: str | None = None):
        self.host = host
        self.remsh_cmd = self.initialize_remsh(remsh_cmd)
        self.hold_command = None
        self.release_command = None
        self.cancel_command = None

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
    ):
        raise NotImplementedError("submit method not implemented.")

    def status(self, ids: str, verbose: bool = False):
        raise NotImplementedError("status method not implemented.")

    def hold(self, remote_ids: str, verbose: bool = False):
        if isinstance(remote_ids, str):
            remote_ids = [remote_ids]

        subprocess_run(
            self.host,
            args=self.hold_command + remote_ids,
            remsh_cmd=self.remsh_cmd,
            verbose=verbose,
        )

    def release(self, remote_ids: str, verbose: bool = False):
        if isinstance(remote_ids, str):
            remote_ids = [remote_ids]

        subprocess_run(
            self.host,
            args=self.release_command + remote_ids,
            remsh_cmd=self.remsh_cmd,
            verbose=verbose,
        )

    def cancel(self, remote_ids: str, verbose: bool = False):
        if isinstance(remote_ids, str):
            remote_ids = [remote_ids]

        subprocess_run(
            self.host,
            args=self.cancel_command + remote_ids,
            remsh_cmd=self.remsh_cmd,
            verbose=verbose,
        )

    def initialize_remsh(self, remsh_cmd: str | None) -> str:
        """
        Initialize the remote shell command.

        Args:
            remsh_cmd (str|None): The remote shell command provided by the user.

        Returns:
            str: The initialized remote shell command.
        """
        if remsh_cmd is None:
            return os.environ.get("RE_RSH", "ssh")
        else:
            return remsh_cmd
