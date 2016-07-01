#
# LSST Data Management System
# Copyright 2016 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
__all__ = ["ExposureIdInfo"]


class ExposureIdInfo(object):
    """!Exposure ID and number of bits used

    Attributes include:
    - expId  exposure ID as a long int
    - expBits  maximum number of bits allowed for exposure IDs
    - maxBits  maximum number of bits available for values that combine exposure ID
        with other information, such as source ID
    - unusedBits  maximum number of bits available for non-exposure info (maxBits - expBits)

    One common use is creating an ID factory for making a source table.
    For example, given a data butler `butler` and a data ID `dataId`:

        from lsst.afw.table import IdFactory, SourceTable
        exposureIdInfo = butler.get("expIdInfo", dataId)
        sourceIdFactory = IdFactory.makeSource(exposureIdInfo.expId, exposureIdInfo.unusedBits)
        schema = SourceTable.makeMinimalSchema()
        #...add fields to schema as desired, then...
        sourceTable = SourceTable.make(self.schema, sourceIdFactory)

    At least one bit must be reserved, even if there is no exposure ID, for reasons
    that are not entirely clear (this is DM-6664).
    """
    def __init__(self, expId=0L, expBits=1, maxBits=64):
        """!Construct an ExposureIdInfo

        See the class doc string for an explanation of the arguments.
        """
        expId = long(expId)
        expBits = int(expBits)
        maxBits = int(maxBits)

        if expId.bit_length() > expBits:
            raise RuntimeError("expId=%s uses %s bits > expBits=%s" % (expId, expId.bit_length(), expBits))
        if maxBits < expBits:
            raise RuntimeError("expBits=%s > maxBits=%s" % (expBits, maxBits))

        self.expId = expId
        self.expBits = expBits
        self.maxBits = maxBits

    @property
    def unusedBits(self):
        return self.maxBits - self.expBits
