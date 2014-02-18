#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Ryan Dwyer on 2013-10-11.
Copyright (c) 2013 Cornell University. All rights reserved.
"""

from autoassign import autoassign


class PropClass(object):
    @autoassign
    def __init__(self, x=2, y=4):
        pass

    @property
    def xy(self):
        return self.x * self.y

    @property
    def xdy(self):
        return self.x / self.y


def test():
    pc = PropClass(x=19294, y=3945755)
    pc.xdy ** 10 + pc.xy ** 10


class UpdateClass(object):
    @autoassign
    def __init__(self, x=2, y=4):
        self._update_params()

    def _update_params(self):
        self.xy = self.x*self.y
        self.xdy = self.x / self.y


def test2():
    uc = UpdateClass(x=19294, y=3945755)
    uc.xdy ** 10 + uc.xy ** 10


# if __name__ == '__main__':
#     from timeit import timeit
#     print(timeit("test()", setup="from __main__ import PropClass, test",
#                  number=1000000))
#     print(timeit("test2()", setup="from __main__ import UpdateClass, test2",
#                  number=1000000))

"""This shows that implementing propreties in python does not result
in a large overhead time cost."""
