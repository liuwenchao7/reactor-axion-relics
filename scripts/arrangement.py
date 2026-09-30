from signal import raise_signal
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle, Rectangle

plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 100
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['font.size'] = 15

# https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARycenter/etd/R8520-406_TPMH1342E.pdf

class Arrangement():
    def __init__(self, **kargs):
        raise NotImplementedError

    def save(self, file, minN=0):
        data = np.vstack([np.arange(minN, self.N + minN), 
                            self.x, 
                            self.y, 
                            self.z,
                            self.rotY,
                            self.rotZ]).T

        # No. of Array, x shift, y shift, z shift
        # rotation angle around y axis, rotation angle around z axis
        fmt = ['%d', '%.4f', '%.4f', '%.4f', '%d', '%.3f']

        np.savetxt(file, 
                data, 
                fmt=fmt, 
                delimiter=' ', 
                comments='', 
                header=f'{self.N}')


class PMTArrangement(Arrangement):
    """
    PMT position manager
    """
    def __init__(self, detr, xeh, top=True) -> None:
        WindowThickness = 1.50  # in mm
        CasingHeight = 27.00  # in mm
        BaseThickness = 1.50  # in mm
        PMTToPMTBase = 3.00  # in mm
        CasingWidth = 25.40  # in mm
        self.pmtw = CasingWidth / 10  # in cm
        self.pmts = 2.78  # in cm
        self.pmts_large = 3.2
        self.offset = self.pmts / 2
        assert self.pmts > self.pmtw
        self.pmth = (WindowThickness + CasingHeight + BaseThickness + PMTToPMTBase) / 10  # in cm
        self.detr = detr
        self.top = top
        self.shift_z = (xeh + 1/2 * self.pmth) # * (1 if self.top else -1)

    def circular(self):
        x = []
        y = []
        angle = []

        x.append([-2 * self.pmts, -self.pmts, 0, self.pmts])
        y.append([-self.offset, -self.offset, -self.offset, -self.offset])
        angle.append([0, 0, 0, 0])

        xtemp = [-self.offset * 3, -self.offset, self.offset]
        ytemp = np.array([self.offset, self.offset * 1.75, self.offset])
        x.append(xtemp)
        y.append(ytemp)
        angle.append([0, 0, 0])
        x.append(xtemp)
        y.append(-ytemp - self.pmts)
        angle.append([0, 0, 0])

        for n in [5, 4, 3, 2]:
            r = ((n * self.pmts) ** 2 + self.offset ** 2) ** 0.5
            rad = 2 * np.arcsin(self.offset / r)
            n = 2 * np.pi // rad
            n = n // 2 * 2
            rt = 2 * np.pi / n
            a = np.arange(n) * rt / np.pi * 180
            angle.append(a)
            x.append(np.cos(a / 180 * np.pi - rad / 2) * r)
            y.append(np.sin(a / 180 * np.pi - rad / 2) * r)

        self.x = np.hstack(x)
        self.y = np.hstack(y)
        self.angle = np.hstack(angle)
        self.get_lrtb()

        self.x = self.xcenter[self.pmtindex]
        self.y = self.ycenter[self.pmtindex]
        self.rotZ = self.angle[self.pmtindex]
        if self.top:
            self.z = np.ones(self.N) * self.shift_z
            self.rotY = np.zeros(self.N)
        else:
            self.z = np.ones(self.N) * self.shift_z
            self.rotY = np.ones(self.N) * 180

    def veto(self, rotY=False):
        x = []
        y = []
        angle = []

        # each PMT separate by 30 degrees
        d = 15
        if rotY:
            r = self.detr - self.pmts_large
            addtional = d / 2
        else:
            r = self.detr - self.pmts
            addtional = 0
        rad = 2 * np.arcsin(self.offset / r)
        n = 360 // d
        rt = 2 * np.pi / n
        a = np.arange(n) * rt / np.pi * 180 + addtional
        angle.append(a)
        x.append(np.cos(a / 180 * np.pi - rad / 2) * r)
        y.append(np.sin(a / 180 * np.pi - rad / 2) * r)

        # (x,y) is the left-bottom corner of a PMT
        self.x = np.hstack(x)
        self.y = np.hstack(y)
        self.angle = np.hstack(angle)
        self.get_lrtb()

        self.x = self.xcenter[self.pmtindex]
        self.y = self.ycenter[self.pmtindex]
        self.rotZ = self.angle[self.pmtindex]
        if self.top:
            self.z = np.ones(self.N) * self.shift_z
            if rotY:
                self.rotY = np.ones(self.N) * 90
            else:
                self.rotY = np.zeros(self.N)
        else:
            self.z = np.ones(self.N) * self.shift_z
            if rotY:
                self.rotY = np.ones(self.N) * 90
            else:
                self.rotY = np.ones(self.N) * 180

    def get_lrtb(self):
        self.xreal = self.x + (self.pmts - self.pmtw) / 2 * np.sqrt(2) * np.cos(np.pi / 4 + self.angle / 180 * np.pi)
        self.yreal = self.y + (self.pmts - self.pmtw) / 2 * np.sqrt(2) * np.sin(np.pi / 4 + self.angle / 180 * np.pi)

        self.xcenter = self.x + self.pmts / 2 * np.sqrt(2) * np.cos(np.pi / 4 + self.angle / 180 * np.pi)
        self.ycenter = self.y + self.pmts / 2 * np.sqrt(2) * np.sin(np.pi / 4 + self.angle / 180 * np.pi)

        lbot = np.array([self.xreal, self.yreal])
        rbot = np.array([self.xreal + self.pmtw * np.cos(self.angle / 180 * np.pi), self.yreal + self.pmtw * np.sin(self.angle / 180 * np.pi)])
        rtop = np.array([self.xreal + self.pmtw * np.cos((self.angle + 45) / 180 * np.pi), self.yreal + self.pmtw * np.sin((self.angle + 45) / 180 * np.pi)])
        ltop = np.array([self.xreal + self.pmtw * np.cos((self.angle + 90) / 180 * np.pi), self.yreal + self.pmtw * np.sin((self.angle + 90) / 180 * np.pi)])
        self.lrtb = np.array([lbot, rbot, rtop, ltop])
        self.lrtb = np.transpose(self.lrtb, axes=[2, 0, 1])
        self.get_index()

    def rectangular(self):
        center = np.arange(-6, 7) * self.pmts
        xcenter, ycenter = np.meshgrid(center + self.offset, center + self.offset)
        self.xcenter = xcenter.flatten()
        self.ycenter = ycenter.flatten()
        self.angle = np.zeros_like(self.xcenter)
        self.xreal = self.xcenter - self.pmtw / 2
        self.yreal = self.ycenter - self.pmtw / 2

        lbot = np.array([self.xcenter - self.offset + (self.pmts - self.pmtw) / 2, self.ycenter - self.offset + (self.pmts - self.pmtw) / 2])
        rbot = np.array([self.xcenter + self.offset - (self.pmts - self.pmtw) / 2, self.ycenter - self.offset + (self.pmts - self.pmtw) / 2])
        rtop = np.array([self.xcenter + self.offset - (self.pmts - self.pmtw) / 2, self.ycenter + self.offset - (self.pmts - self.pmtw) / 2])
        ltop = np.array([self.xcenter - self.offset + (self.pmts - self.pmtw) / 2, self.ycenter + self.offset - (self.pmts - self.pmtw) / 2])
        lrtb = np.array([lbot, rbot, rtop, ltop])
        self.lbot = lbot
        self.lrtb = np.transpose(lrtb, axes=[2, 0, 1])
        self.get_index()

        self.x = self.xcenter[self.pmtindex]
        self.y = self.ycenter[self.pmtindex]
        self.rotZ = self.angle[self.pmtindex]
        if self.top:
            self.z = np.ones(self.N) * self.shift_z
            self.rotY = np.zeros(self.N)
        else:
            self.z = np.ones(self.N) * self.shift_z
            self.rotY = np.ones(self.N) * 180

    def get_index(self):
        self.overlap_check()
        self.detcenter = np.array([0, 0])
        self.pmtindex = []
        for i in range(len(self.lrtb)):
            if np.all(np.power(self.lrtb[i] - self.detcenter, 2).sum(axis=1) < self.detr ** 2):
                self.pmtindex.append(i)
        self.N = len(self.pmtindex)

    @np.errstate(invalid='ignore')
    def overlap_check(self):
        for ii in range(self.lrtb.shape[0]):
            for jj in range(self.lrtb.shape[1]):
                sub = np.delete(self.lrtb, ii, axis=0) - self.lrtb[ii][jj]
                accum = np.zeros(sub.shape[0])

                for i, j in zip([0, 1, 2, 3], [1, 2, 3, 0]):
                    l0 = sub[:, i, :]
                    l1 = sub[:, j, :]
                    cos = (l0 * l1).sum(axis=1) / np.sqrt((l0 ** 2).sum(axis=1) * (l1 ** 2).sum(axis=1))
                    accum += np.arccos(cos)

                # assert np.sum(np.isnan(accum)) <= 1
                assert np.sum(np.abs(accum - 2 * np.pi)[~np.isnan(accum)] < 1e-2) <= 0

    def draw(self, file=None):
        fig = plt.figure(figsize=(6, 6))
        # fig.tight_layout()
        gs = gridspec.GridSpec(1, 1, figure=fig, left=0.15, right=0.90, top=0.90, bottom=0.15, wspace=0.3, hspace=0.3)

        ax = fig.add_subplot(gs[0, 0])
        ax.scatter(*self.detcenter, label='Center', color='r')
        ax.add_patch(Circle(xy=[0, 0], radius=self.detr, alpha=1.0, facecolor='none', edgecolor='r'))
        for i in self.pmtindex:
            ax.add_patch(Rectangle(xy=[self.xreal[i], self.yreal[i]], 
                                width=self.pmtw, height=self.pmtw, angle=self.angle[i], 
                                alpha=1.0, facecolor='none', edgecolor='b'))
        ax.set_xlabel('x [cm]')
        ax.set_ylabel('y [cm]')
        ax.legend(title=f'N={self.N}')
        ax.grid()
        if not file is None:
            fig.savefig(file, transparent=True, dpi=400)


class TeflonArrangement(Arrangement):
    """
    Teflon position manager
    """
    def __init__(self, sides, radius, copper_w, teflon_w, xeh) -> None:
        self.N = sides
        self.radius = radius
        self.copper_w = copper_w
        self.teflon_w = teflon_w
        self.xeh = xeh

    def circular(self):
        rad = 2 * np.pi / self.N
        angle = (np.arange(self.N) + 0.5) * rad
        self.x = np.cos(angle) * (self.radius + self.copper_w / 2 + self.teflon_w)
        self.y = np.sin(angle) * (self.radius + self.copper_w / 2 + self.teflon_w)
        self.z = np.ones(self.N) * self.xeh
        self.rotY = np.zeros(self.N)
        self.rotZ = angle * 180 / np.pi + 90


class TankArrangement(Arrangement):
    """
    Tank outline manager
    """
    def __init__(self, radius, height, thickness=[0]) -> None:
        self.thickness = thickness
        self.radius = radius
        self.height = height

    def getTankOutline(self, radius_ratio=1/8, angle=np.pi*2/5, spline_points=11):
        thickness = np.array(self.thickness)[:, None]
        self.radius_ratio = radius_ratio
        x = []
        y = []
        self.small_radius = self.radius_ratio * self.radius
        self.small_angle = np.linspace(0, angle, spline_points)
        x.append(np.cos(self.small_angle) * (self.small_radius + thickness) + (1 - self.radius_ratio) * self.radius)
        y.append(np.sin(self.small_angle) * (self.small_radius + thickness))
        large_radius = (1 - self.radius_ratio) * self.radius / np.cos(angle) + self.radius_ratio * self.radius
        large_angle = np.linspace(angle, np.pi / 2, spline_points)
        x.append(np.cos(large_angle) * (large_radius + thickness))
        y.append(np.sin(large_angle) * (large_radius + thickness) - (1 - self.radius_ratio) * self.radius * np.tan(angle))
        x = np.hstack(x)
        y = np.hstack(y)
        # top to bottom
        self.x = np.hstack([x[:, ::-1], x])
        self.y = np.hstack([-y[:, ::-1] - self.height / 2, y + self.height / 2])
        self.N = self.x.shape[1]

    def draw(self, file=None):
        fig = plt.figure(figsize=(6, 6))
        # fig.tight_layout()
        gs = gridspec.GridSpec(1, 1, figure=fig, left=0.15, right=0.95, top=0.95, bottom=0.10, wspace=0.15, hspace=0.2)
        ax = fig.add_subplot(gs[0, 0])
        for i in range(len(self.thickness)):
            ax.plot(self.x[i], self.y[i])
        ax.plot([0, np.cos(self.small_angle[-1]) * self.small_radius + (1 - self.radius_ratio) * self.radius],
                [0, np.sin(self.small_angle[-1]) * self.small_radius + self.height / 2], color='k', linestyle='dashed')
        ax.plot([(1 - self.radius_ratio) * self.radius, self.radius],
                [self.height / 2, self.height / 2], color='k', linestyle='dashed')
        ax.scatter([0, (1 - self.radius_ratio) * self.radius], [0, self.height / 2], color='r', s=4.0)
        # ax.set_xscale('log')
        ax.set_xlabel('x [cm]')
        ax.set_ylabel('y [cm]')
        ax.grid()
        ax.set_aspect('equal', adjustable='box')
        # ax.legend()
        if not file is None:
            fig.savefig(file, transparent=True, dpi=400)

    def save(self, file, i, shift_z=0,
            hight_cut=None, upward=True,
            flange_height=None,
            flange_width=None,
            flange_to_welding=None):
        x = self.x[i].copy()
        y = self.y[i].copy()
        if not flange_height is None:
            x = np.hstack([x[:len(x)//2], [x.max(),
                        x.max() + flange_width,
                        x.max() + flange_width,
                        x.max()], x[len(x)//2:]])
            y = np.hstack([y[:len(y)//2], [self.height / 2 - flange_height - flange_to_welding,
                        self.height / 2 - flange_height - flange_to_welding,
                        self.height / 2 - flange_to_welding,
                        self.height / 2 - flange_to_welding], y[len(y)//2:]])
        y += shift_z
        if not hight_cut is None:
            assert hight_cut - shift_z < self.height / 2
            if upward:
                x = np.hstack([[0, self.thickness[i] + self.radius], x[y > hight_cut]])
                y = np.hstack([[hight_cut, hight_cut], y[y > hight_cut]])
            else:
                x = np.hstack([x[y < hight_cut], [self.thickness[i] + self.radius, 0]])
                y = np.hstack([y[y < hight_cut], [hight_cut, hight_cut]])
        # i-th thickness
        data = np.vstack([np.arange(len(x)), x, y]).T

        # No. of Array, r, z
        # rotation angle around y axis, rotation angle around z axis
        fmt = ['%d', '%.4f', '%.4f']

        np.savetxt(file, 
                data, 
                fmt=fmt, 
                delimiter=' ', 
                comments='', 
                header=f'{len(x)}')
