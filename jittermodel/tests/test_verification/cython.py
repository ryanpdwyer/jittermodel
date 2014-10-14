#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import timeit


def format_time(t):
    if t > 1:
        formatted_time = "{0:.3f} s".format(t)
    elif t > 1e-3:
        formatted_time = "{0:.2f} ms".format(t * 1e3)
    elif t > 1e-6:
        formatted_time = "{0:.1f} us".format(t * 1e6)
    else:
        formatted_time = "{0:.0f} ns".format(t * 1e9)

    return formatted_time

def time_thetaI(call, n):
    t = timeit.timeit(
        """{call}(
        1, 0.1, 0.0001 + 0.0165j,3500, 2000+15j,3-0.001j, 3 - 100j)""".format(
        call=call),
        setup="import jittermodel._sim; import jittermodel.simulation",
        number=n)
    print("{call}: {t}".format(call=call, t=format_time(t/n)))


def time_thetaII(call, n):
    t = timeit.timeit(
        """{call}(
        1, 0.1, 3 - 0.001j, 3 - 0.001j, 3 - 100j, 0.00012 + 0.0165j)""".format(
        call=call),
        setup="import jittermodel._sim; import jittermodel.simulation",
        number=n)
    print("{call}: {t}".format(call=call, t=format_time(t/n)))


def time_im_dielectric_2(call, n):
    t = timeit.timeit(
        """{call}(1, 0.1, 0.001, 3 - 0.001j, 4.65, 3e-4, 300, 1e5, 300, 0.0138,
        0.1602, 0.0088542, 2)""".format(
        call=call),
        setup="import jittermodel._sim; import jittermodel.simulation",
        number=n)
    print("{call} model 2: {t}".format(call=call, t=format_time(t/n)))

def time_im_dielectric_1(call, n):
    t = timeit.timeit(
        """{call}(1, 0.1, 0.001, 3 - 0.001j, 4.65, 3e-4, 300, 1e5, 300, 0.0138,
        0.1602, 0.0088542, 1)""".format(
        call=call),
        setup="import jittermodel._sim; import jittermodel.simulation",
        number=n)
    print("{call} model 1: {t}".format(call=call, t=format_time(t/n)))


def main():
    time_thetaI("jittermodel.simulation._thetaI", 1000)
    time_thetaI("jittermodel._sim._thetaI_np", 10000)
    time_thetaI("jittermodel._sim._thetaI_math", 10000)
    time_thetaI("jittermodel._sim._thetaI_c", 1000000)

    time_thetaII("jittermodel.simulation._thetaII", 10000)
    time_thetaII("jittermodel._sim._thetaII_math", 100000)
    time_thetaII("jittermodel._sim._thetaII_c", 10000000)

    time_im_dielectric_2("jittermodel.simulation._im_dielectric", 1000)
    time_im_dielectric_2("jittermodel._sim._im_dielectric", 10000)
    time_im_dielectric_2("jittermodel._sim._im_dielectric_c", 10000)


    time_im_dielectric_1("jittermodel.simulation._im_dielectric", 1000)
    time_im_dielectric_1("jittermodel._sim._im_dielectric", 10000)
    time_im_dielectric_1("jittermodel._sim._im_dielectric_c", 10000)


if __name__ == '__main__':
    main()
