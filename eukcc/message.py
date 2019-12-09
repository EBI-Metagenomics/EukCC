#!/usr/bin/env python3
#
# This file is part of the EukCC (https://github.com/openpaul/eukcc).
# Copyright (c) 2019 Paul Saary
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# provides all file operation functions
# used inthis package
import datetime


def log(m, quiet=False):
    ts = datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S")
    if not quiet:
        print("{}: {}".format(ts, m))
