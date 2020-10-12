# This file is part of obs_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import lsst.daf.persistence as dafPersist
from lsst.obs.base import CameraMapper

__all__ = ["CompositeMapper"]


class CompositeMapper(CameraMapper):
    packageName = "obs_base"

    def __init__(self, root, policy=None, **kwargs):
        if policy is None:
            policy = dafPersist.Policy()
        super(CompositeMapper, self).__init__(policy, repositoryDir=root, root=root, **kwargs)

    def _makeCamera(self, policy, repositoryDir):
        """Normally this makes a camera. For composite testing, we don't need a
        camera.
        """
        return None

    def std_stdTestType(self, item, dataId):
        item.standardized = True
        return item

    def bypass_bypassTestType(self, datasetType, pythonType, location, dataId):
        return set(dataId.keys())
