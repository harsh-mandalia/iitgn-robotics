import math
import os
import sys
from enum import IntEnum

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.relpath("../../Mapping/grid_map_lib/"))
try:
    from grid_map_lib import GridMap
except ImportError:
    raise

do_animation = True

class SweepSearcher:
    class SweepDirection(IntEnum):
        UP = 1
        DOWN = -1

    class MovingDirection(IntEnum):
        RIGHT = 1
        LEFT = -1

    def __init__(self, mdirection, sdirection, xinds_goaly, goaly):
        self.moving_direction = mdirection
        self.sweep_direction = sdirection
        self.turing_window = []
        self.update_turning_window()
        self.xinds_goaly = xinds_goaly
        self.goaly = goaly

    def move_target_grid(self, cxind, cyind, gmap):
        nxind = self.moving_direction + cxind
        nyind = cyind

        # found safe grid
        if not gmap.check_occupied_from_xy_index(nxind, nyind, occupied_val=0.5):
            return nxind, nyind
        else:  # occupided
            ncxind, ncyind = self.find_safe_turning_grid(cxind, cyind, gmap)
            if (ncxind is None) and (ncyind is None):
                # moving backward
                ncxind = -self.moving_direction + cxind
                ncyind = cyind
                if gmap.check_occupied_from_xy_index(ncxind, ncyind):
                    # moved backward, but the grid is occupied by obstacle
                    return None, None
            else:
                # keep moving until end
                while not gmap.check_occupied_from_xy_index(ncxind + self.moving_direction, ncyind, occupied_val=0.5):
                    ncxind += self.moving_direction
                self.swap_moving_direction()
            return ncxind, ncyind

    def find_safe_turning_grid(self, cxind, cyind, gmap):

        for (dxind, dyind) in self.turing_window:

            nxind = dxind + cxind
            nyind = dyind + cyind

            # found safe grid
            if not gmap.check_occupied_from_xy_index(nxind, nyind, occupied_val=0.5):
                return nxind, nyind

        return None, None

    def is_search_done(self, gmap):
        for ix in self.xinds_goaly:
            if not gmap.check_occupied_from_xy_index(ix, self.goaly, occupied_val=0.5):
                return False

        # all lower grid is occupied
        return True

    def update_turning_window(self):
        self.turing_window = [
            (self.moving_direction, 0.0),
            (self.moving_direction, self.sweep_direction),
            (0, self.sweep_direction),
            (-self.moving_direction, self.sweep_direction),
        ]

    def swap_moving_direction(self):
        self.moving_direction *= -1
        self.update_turning_window()

    def search_start_grid(self, grid_map):
        xinds = []
        y_ind = 0
        if self.sweep_direction == self.SweepDirection.DOWN:
            xinds, y_ind = search_free_grid_index_at_edge_y(grid_map, from_upper=True)
        elif self.sweep_direction == self.SweepDirection.UP:
            xinds, y_ind = search_free_grid_index_at_edge_y(grid_map, from_upper=False)

        if self.moving_direction == self.MovingDirection.RIGHT:
            return min(xinds), y_ind
        elif self.moving_direction == self.MovingDirection.LEFT:
            return max(xinds), y_ind

        raise ValueError("self.moving direction is invalid ")

def find_sweep_direction_and_start_posi(ox, oy):
    # find sweep_direction
    max_dist = 0.0
    vec = [0.0, 0.0]
    sweep_start_pos = [0.0, 0.0]
    for i in range(len(ox) - 1):
        dx = ox[i + 1] - ox[i]
        dy = oy[i + 1] - oy[i]
        d = np.sqrt(dx ** 2 + dy ** 2)

        if d > max_dist:
            max_dist = d
            vec = [dx, dy]
            sweep_start_pos = [ox[i], oy[i]]

    return vec, sweep_start_pos


def convert_grid_coordinate(ox, oy, sweep_vec, sweep_start_posi):
    tx = [ix - sweep_start_posi[0] for ix in ox]
    ty = [iy - sweep_start_posi[1] for iy in oy]

    th = math.atan2(sweep_vec[1], sweep_vec[0])

    c = np.cos(-th)
    s = np.sin(-th)

    rx = [ix * c - iy * s for (ix, iy) in zip(tx, ty)]
    ry = [ix * s + iy * c for (ix, iy) in zip(tx, ty)]

    return rx, ry


def convert_global_coordinate(x, y, sweep_vec, sweep_start_posi):
    th = math.atan2(sweep_vec[1], sweep_vec[0])
    c = np.cos(th)
    s = np.sin(th)

    tx = [ix * c - iy * s for (ix, iy) in zip(x, y)]
    ty = [ix * s + iy * c for (ix, iy) in zip(x, y)]

    rx = [ix + sweep_start_posi[0] for ix in tx]
    ry = [iy + sweep_start_posi[1] for iy in ty]

    return rx, ry


def search_free_grid_index_at_edge_y(grid_map, from_upper=False):
    yind = None
    xinds = []

    if from_upper:
        xrange = range(grid_map.height)[::-1]
        yrange = range(grid_map.width)[::-1]
    else:
        xrange = range(grid_map.height)
        yrange = range(grid_map.width)

    for iy in xrange:
        for ix in yrange:
            if not grid_map.check_occupied_from_xy_index(ix, iy):
                yind = iy
                xinds.append(ix)
        if yind:
            break

    return xinds, yind


def setup_grid_map(ox, oy, reso, sweep_direction, offset_grid=10):
    width = math.ceil((max(ox) - min(ox)) / reso) + offset_grid
    height = math.ceil((max(oy) - min(oy)) / reso) + offset_grid
    center_x = np.mean(ox)
    center_y = np.mean(oy)

    grid_map = GridMap(width, height, reso, center_x, center_y)

    grid_map.set_value_from_polygon(ox, oy, 1.0, inside=False)

    grid_map.expand_grid()

    xinds_goaly = []
    goaly = 0
    if sweep_direction == SweepSearcher.SweepDirection.UP:
        xinds_goaly, goaly = search_free_grid_index_at_edge_y(grid_map, from_upper=True)
    elif sweep_direction == SweepSearcher.SweepDirection.DOWN:
        xinds_goaly, goaly = search_free_grid_index_at_edge_y(grid_map, from_upper=False)

    return grid_map, xinds_goaly, goaly


def sweep_path_search(sweep_searcher, gmap, grid_search_animation=False):
    # search start grid
    cxind, cyind = sweep_searcher.search_start_grid(gmap)
    if not gmap.set_value_from_xy_index(cxind, cyind, 0.5):
        print("Cannot find start grid")
        return [], []

    x, y = gmap.calc_grid_central_xy_position_from_xy_index(cxind, cyind)
    px, py = [x], [y]

    if grid_search_animation:
        fig, ax = plt.subplots()

    while True:
        cxind, cyind = sweep_searcher.move_target_grid(cxind, cyind, gmap)

        if sweep_searcher.is_search_done(gmap) or (cxind is None or cyind is None):
            print("Done")
            break

        x, y = gmap.calc_grid_central_xy_position_from_xy_index(
            cxind, cyind)

        px.append(x)
        py.append(y)

        gmap.set_value_from_xy_index(cxind, cyind, 0.5)

        if grid_search_animation:
            gmap.plot_grid_map(ax=ax)
            plt.pause(1.0)

    gmap.plot_grid_map()

    return px, py


def planning(ox, oy, reso,
             moving_direction=SweepSearcher.MovingDirection.RIGHT,
             sweeping_direction=SweepSearcher.SweepDirection.UP,
             ):
    sweep_vec, sweep_start_posi = find_sweep_direction_and_start_posi(ox, oy)

    rox, roy = convert_grid_coordinate(ox, oy, sweep_vec, sweep_start_posi)

    gmap, xinds_goaly, goaly = setup_grid_map(rox, roy, reso, sweeping_direction)

    sweep_searcher = SweepSearcher(moving_direction, sweeping_direction, xinds_goaly, goaly)

    px, py = sweep_path_search(sweep_searcher, gmap)

    rx, ry = convert_global_coordinate(px, py, sweep_vec, sweep_start_posi)

    print("Path length:", len(rx))

    return rx, ry

def planning_animation(ox, oy, reso):  # pragma: no cover
    px, py = planning(ox, oy, reso)

    # animation
    if do_animation:
        for ipx, ipy in zip(px, py):
            plt.cla()
            plt.plot(ox, oy, "-xb")
            plt.plot(px, py, "-r")
            plt.plot(ipx, ipy, "or")
            plt.axis("equal")
            plt.grid(True)
            if plt.waitforbuttonpress(0.01):
                break

def main():  # pragma: no cover
    print("start!!")

    # ox = [0.0, 20.0, 50.0, 100.0, 130.0, 40.0, 0.0]
    # oy = [0.0, -20.0, 0.0, 30.0, 60.0, 80.0, 0.0]
    # reso = 5.0
    # planning_animation(ox, oy, reso)

    # ox = [0.0, 50.0, 50.0, 0.0, 0.0]
    # oy = [0.0, 0.0, 30.0, 30.0, 0.0]
    # reso = 1.3
    # planning_animation(ox, oy, reso)

    # ox = [0.0, 20.0, 50.0, 200.0, 130.0, 40.0, 0.0]
    # oy = [0.0, -80.0, 0.0, 30.0, 60.0, 80.0, 0.0]
    # reso = 5*2

    kml1="80.29821599780011,13.13286087620951,0 80.29796720705747,13.13182773034574,0 80.2991076628975,13.13169660518675,0 80.29944038024139,13.13165299234652,0 80.29937762620091,13.13109111234745,0 80.29793726328148,13.13128457469073,0 80.29761056937502,13.13006729251901,0 80.29806252902425,13.12999465932925,0 80.29835086427846,13.13024456717143,0 80.29941213185222,13.12998434135773,0 80.30023518680595,13.12965420556122,0 80.30044586533064,13.12906413112875,0 80.30050201356336,13.12866328487464,0 80.30036109015327,13.12804180049814,0 80.30025662226946,13.12776937285265,0 80.3000525602131,13.12752807411891,0 80.29968951597387,13.12732678106736,0 80.29850101748954,13.12753770530468,0 80.29819061015509,13.12698463906859,0 80.29826246128076,13.12663847836676,0 80.29863060875304,13.12616498023377,0 80.29866173571808,13.12568598023374,0 80.29872099323003,13.12514727753639,0 80.29866665209987,13.1245795774647,0 80.29799448684999,13.12443297917521,0 80.29745577682087,13.1246465082633,0 80.29777843797754,13.12406362090931,0 80.29901796325314,13.12452616431568,0 80.2994215614677,13.12418173688713,0 80.29890539382183,13.1237513033455,0 80.29843763117461,13.12349401672662,0 80.29966445016528,13.12313628366699,0 80.30020583464541,13.12344717000246,0 80.30023942717858,13.12346802009502,0 80.30066050587938,13.12375513411228,0 80.30091433792411,13.12427053640675,0 80.30122652457746,13.12483053454407,0 80.30162645319366,13.12545180733082,0 80.30201672601071,13.12617996201633,0 80.30261722067974,13.12696180815264,0 80.30289554473505,13.12760244671527,0 80.30330431624785,13.12807993910696,0 80.30333822855845,13.12813887406227,0 80.30371150284772,13.12883379102835,0 80.30402529417192,13.12914611634577,0 80.30388023777046,13.12927645834875,0 80.30334712621725,13.12996559499788,0 80.30301378651669,13.12975218051434,0 80.30107548260035,13.13287222952495,0 80.30069077457982,13.13320044898597,0 80.2995692553596,13.1334836836321,0 80.2994531642517,13.13318768783351,0 80.30037261724844,13.13287222535551,0 80.30060620803647,13.13244118435543,0 80.30019561062818,13.1323564180065,0 80.29969301828299,13.13235642075422,0 80.29913484053139,13.13250202156531,0 80.29843887307364,13.13277263561584,0 80.29821599780011,13.13286087620951,0"
    kml2="80.29896188177538,13.13457871394119,0 80.30179705144192,13.13393332239194,0 80.30251262471795,13.13258722535745,0 80.3037296832117,13.13018483182888,0 80.30439533014118,13.12995472237192,0 80.30455207341703,13.13097186008446,0 80.30469291859053,13.13205528046478,0 80.3074130726969,13.13881532179428,0 80.30969697138669,13.14515284587348,0 80.31109872275637,13.14839401014196,0 80.31550715366353,13.15925779755686,0 80.31269281527685,13.16050495443757,0 80.31218459464893,13.15938846639837,0 80.3096718932964,13.15980990169297,0 80.30836266443212,13.15689452707756,0 80.3072121392284,13.15451115489642,0 80.30830948366555,13.15412619754567,0 80.30801269745258,13.15313937038605,0 80.3077707179326,13.1526328941917,0 80.3064751044341,13.15288880296986,0 80.30517696832142,13.14955735572043,0 80.30657288759892,13.1490290278613,0 80.30605448482379,13.14784294170603,0 80.30462792442583,13.14825030518204,0 80.30289714996391,13.14510334409799,0 80.30425576758503,13.14446955608774,0 80.30385630492628,13.14330550518714,0 80.3025426440617,13.14327683884104,0 80.30063763636079,13.14084518677255,0 80.30074407799306,13.14083845736187,0 80.30250238236535,13.14030651352447,0 80.3022178732283,13.13946864712285,0 80.30060454392019,13.14006144559007,0 80.29896188177538,13.13457871394119,0"
    array=kml2.split(" ")
    ox=[]
    oy=[]
    for i in range(len(array)):
        lst=array[i].split(",")
        x=(float(lst[0])-80.29821599780011)*100000
        y=(float(lst[1])-13.13286087620951)*100000
        ox.append(x)
        oy.append(y)
    # reso=20
    reso=25
    # plt.plot(ox, oy, "-xb")
    planning_animation(ox, oy, reso)

    plt.show()

    print("done!!")

main()
