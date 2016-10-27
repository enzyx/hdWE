#
# This file is part of hdWE. 
# Copyright (C) 2016 Manuel Luitz <manuel.luitz@tum.de>
# Copyright (C) 2016 Rainer Bomblies <r.bomblies@tum.de>
# Copyright (C) 2016 Fabian Zeller
#
# hdWE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# hdWE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with hdWE. If not, see <http://www.gnu.org/licenses/>.
# 
"""
Natural/Physical Constants
"""
# kT in in kcal/mol*K at 298 K
kT = 0.0019872041 * 298 
# kT for langevin dynamics
# kT = 1.
# lower boundary for numerical correctness of divisions by zero
num_boundary = 1e-99

def getLogDirPath(WORKDIR, JOBNAME):
    """
    Global definition of the logdir name
    """
    return "{wd}/{jn}-log/".format(wd=WORKDIR, jn=JOBNAME)

def getRunDirPath(WORKDIR, JOBNAME):
    """
    Global definition of the logdir name
    """
    return "{wd}/{jn}-run/".format(wd=WORKDIR, jn=JOBNAME)

