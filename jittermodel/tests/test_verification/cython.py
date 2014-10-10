# -*- coding: utf-8 -*-
from __future__ import division
import timeit


def _format_time(t):
    if t > 1:
        formatted_time = "{0:.3f} s".format(t)
    elif t > 1e-3:
        formatted_time = "{0:.2f} ms".format(t * 1e3)
    elif t > 1e-6:
        formatted_time = "{0:.1f} us".format(t * 1e6)
    else:
        formatted_time = "{0:.0f} ns".format(t * 1e9)

    return formatted_time


def time_theta(call, n):
    t = timeit.timeit(
        """{call}(
        1, 0.1, 0.0001 + 0.0165j,3500, 2000+15j,3-0.001j, 3 - 100j)""".format(
        call=call),
        setup="import jittermodel._sim; import jittermodel.simulation",
        number=n)
    print("{call}: {t}".format(call=call, t=_format_time(t/n)))


def main():
    time_theta("jittermodel.simulation._thetaI", 1000)
    time_theta("jittermodel._sim._thetaI_np", 10000)
    time_theta("jittermodel._sim._thetaI_math", 10000)


if __name__ == '__main__':
    main()
