# Copyright (c) 2021 Valerii Sukhorukov. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the
#    names of its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER ''AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------

"""Functions extracting empitical length and curvature distributions.

Source publication:
Zhen Zhang Yukako Nishimura, and Pakorn Kanchanawonga 'Extracting
microtubule networks from superresolution single-molecule localization
microscopy data' MBoC, 2017.
"""

from collections import namedtuple
from typing import Sequence

Extent = namedtuple('Range', ['begin', 'end'])


def curvature(density: bool):
    """Curvature data extracted from publication plots.
    """

    pix = Extent(begin=[177, 900], end=[1195, 81])
    val = Extent(begin=[0, 0],     end=[1.5, 1])
    # Blue:
    contro_pix = [0,   696, 251, 81,  109, 232, 458, 591, 685, 738, 797, 822,
                  832, 852, 865, 871, 874, 887, 880, 891, 895, 894, 897, 895,
                  895, 895, 895, 895, 899, 900, 899]
    # Red:
    ca_ras_pix = [0,   702, 228, 81,  151, 164, 330, 527, 616, 711, 770, 808,
                  817, 799, 820, 830, 867, 881, 876, 886, 888, 899, 897, 886,
                  899, 895, 895, 899, 899, 900, 899]

    return _process(pix, val, [contro_pix, ca_ras_pix],
                    start=1, density=density)


def length(density: bool) -> tuple[list, list[list]]:
    """Length data extracted from publication plots.
    """

    pix = Extent(begin=[121, 873], end=[1000, 75])
    val = Extent(begin=[0, 0],     end=[20, 1])
    # Green:
    contro_pix = [0,   0,   441, 75,  201, 339, 425, 477, 493, 567, 615, 604,
                  652, 656, 676, 702, 714, 720, 740, 759, 742, 788, 800, 788,
                  801, 809, 809, 815, 829, 813, 833, 822, 848, 832, 846, 833,
                  858, 859, 854, 854, 858, 854, 858, 858, 858, 857, 864, 856,
                  855, 858, 861, 860, 863, 863, 863, 873, 864, 868, 868, 870]
    # Red:
    ca_ras_pix = [0,   0,   453, 75,  256, 372, 388, 453, 460, 548, 589, 632,
                  620, 651, 612, 724, 718, 701, 738, 772, 749, 755, 758, 804,
                  786, 830, 786, 799, 814, 827, 811, 800, 833, 837, 816, 820,
                  833, 833, 854, 852, 829, 843, 847, 847, 847, 857, 855, 865,
                  868, 872, 872, 851, 868, 863, 862, 868, 872, 857, 873, 865]

    return _process(pix, val, [contro_pix, ca_ras_pix],
                    start=2, density=density)


def _process(
        pix: Extent,
        val: Extent,
        data: Sequence[Sequence],
        start: int = 0,
        density: bool = False,
) -> tuple[list, list[list]]:

    range_pix = [pix.end[0] - pix.begin[0],
                 pix.begin[1] - pix.end[1]]
    range_val = [val.end[i] - val.begin[i] for i in range(2)]
    scale = [range_val[i] / range_pix[i] for i in range(2)]

    n = len(data[0])
    bin_pix = range_pix[0] / n
    bin_val = bin_pix * scale[0]

    x = [i * bin_val for i in range(start, n)]
    res = [[(pix.begin[1] - v) * scale[1] for v in d[start:]] for d in data]
    if density:
        s = (sum(r) for r in res)
        res = ([r/ss/bin_val for r in rr] for rr, ss in zip(res, s))
    return x, res


def avg(
        xx: Sequence,
        yy: Sequence
) -> float:
    """Averages from frequency distribution.
    """

    dx = xx[1] - xx[0]
    a = sum(yy) * dx
    yy = [y / a for y in yy]

    return sum([x * y for x, y in zip(xx, yy)]) * dx


if __name__ == '__main__':

    as_density = True
    x_curv, (control_curv, ca_ras_curv) = curvature(as_density)
    control_curv_avg = avg(x_curv, control_curv)
    x_leng, (control_leng, ca_ras_leng) = length(as_density)
    control_leng_avg = avg(x_leng, control_leng)
