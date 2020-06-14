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


import click

from lsst.daf.butler.cli.utils import addArgumentHelp, ParameterType


class instrument_parameter:  # noqa: N801

    defaultHelp = "The name or fully-qualified class name of an instrument."

    def __init__(self, parameterType=ParameterType.OPTION, required=False, help=defaultHelp, multiple=False):
        self.parameterType = parameterType
        self.required = required
        self.help = help
        self.multiple = multiple

    def __call__(self, f):
        if self.parameterType == ParameterType.OPTION:
            return click.option("-i", "--instrument",
                                required=self.required,
                                help=self.help,
                                multiple=self.multiple)(f)
        if self.parameterType == ParameterType.ARGUMENT:
            f.__doc__ = addArgumentHelp(f.__doc__, self.help)
            return click.argument("instrument",
                                  required=self.required,
                                  nargs=-1 if self.multiple else 1,
                                  metavar="INSTRUMENT ..." if self.multiple else "INSTRUMENT")(f)
