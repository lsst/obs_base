# This file is part of obs_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Accessor for sphinx documentation generator that accesses the butler
subcommand plugins that are provided by this package.
"""

from __future__ import annotations

from typing import cast

import click

from lsst.utils import doImportType

from .. import cmd


class ButlerCmdDocGen(click.Group):
    """Provide access of butler subcommand plugins to Sphinx."""

    def list_commands(self, ctx: click.Context) -> list[str]:
        """List the click commands provided by this package.

        Parameters
        ----------
        ctx : click.Context
            The current Click context.

        Returns
        -------
        commands : `list` [`str`]
            The names of the commands that can be called by the butler command.
        """
        return cmd.__all__

    def get_command(self, ctx: click.Context, cmd_name: str) -> click.Command | None:
        """Get a click command provided by this package.

        Parameters
        ----------
        ctx : `click.Context`
            The current Click context.
        cmd_name : `str`
            The name of the command to return.

        Returns
        -------
        command : `click.Command`
            A Command that wraps a callable command function.
        """
        return cast(click.Command | None, doImportType("lsst.obs.base.cli.cmd." + cmd_name))


@click.command(cls=ButlerCmdDocGen)
def cli() -> None:
    """Run the command."""
    pass
