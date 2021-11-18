from numpy import sin,cos,tan,arctan2,pi,arctan,Inf
import shapely.geometry as gp
from math import sqrt
import matplotlib.pyplot as plt
import math
import os
import sys
from enum import IntEnum
import numpy as np
sys.path.append(os.path.relpath("../../Mapping/grid_map_lib/"))
try:
    from grid_map_lib import GridMap
except ImportError:
    raise

from ortools.constraint_solver import pywrapcp
from ortools.constraint_solver import routing_enums_pb2

def longestside(points):
  edgelengths = [sqrt((points[idx + 1][0] - points[idx][0])**2 + (points[idx + 1][1] - points[idx][1])**2) for idx in range(len(points) - 1)]
  #print(edgelengths)
  lastedge=sqrt((points[-1][0]-points[0][0])**2+(points[-1][1] - points[0][1])**2)
  edgelengths.append(lastedge)
  #print(edgelengths)
  
  if(max(edgelengths)>lastedge):
    maxpos = edgelengths.index(max(edgelengths))
    return points[maxpos],points[maxpos+1]
  else:
    return points[-1],points[0]
#longestside(PL)
# longestside(PointList)

# Operation: Shift origin and rotate axes to align with longest side

#Shift the origin
def originShift(points,newOrigin):
  newPoints=[(points[idx][0]-newOrigin[0],points[idx][1]-newOrigin[1]) for idx in range(len(points))]
  return newPoints
#originShift(PL,(100,101))

#Rotate the axis about origin along edge passing through rotation point
def AxisRotate(points,RotPoint):
  theta=arctan2(RotPoint[1],RotPoint[0])
  newPoints=[(points[idx][0]*cos(theta)+points[idx][1]*sin(theta),-points[idx][0]*sin(theta)+points[idx][1]*cos(theta)) for idx in range(len(points))]
  return newPoints,theta
#print(AxisRotate(PL,PL[1]))

def operatePolygon(points):
  #Extract longest side
  Longside=longestside(points)
  P1=Longside[0]
  P2=Longside[1]
  idxP1=points.index(P1)
  idxP2=points.index(P2)

  
  #Shift origin to P1
  pointsO1=originShift(points,P1)
  #Rotate axes to make P1-P2 as x-axis
  pointsO2,theta=AxisRotate(pointsO1,pointsO1[idxP2])


  return pointsO2,idxP1,idxP2,theta
# print(operatePolygon(PL))

#Decompose cell of operated polygon i.e pointsO2
def CellDecomposition(points):

  #xmin
  P_xmin=min(points)
  # print(P_xmin)
  P_xmin=gp.Point(P_xmin)

  #xmax
  P_xmax=max(points)
  # print(P_xmax)
  P_xmax=gp.Point(P_xmax)

  P_ymin=min(points,key = lambda i: i[1])
  # print(P_ymin)
  P_ymin=gp.Point(P_ymin)

  P_ymax=max(points,key = lambda i: i[1])
  # print(P_ymax)
  P_ymax=gp.Point(P_ymax)

  line_max = gp.LineString([gp.Point(P_xmax.x,P_ymax.y),gp.Point(P_xmax.x,P_ymin.y)])
  line_min = gp.LineString([gp.Point(P_xmin.x,P_ymax.y),gp.Point(P_xmin.x,P_ymin.y)])
  poly_points=gp.Polygon(points)
  sweep0=[]
  for i in points:
    sweep0.append(((i[0],P_ymax.y),(i[0],P_ymin.y)))
    plt.plot([i[0],i[0]], [P_ymax.y,P_ymin.y])
  Sweep=gp.MultiLineString(sweep0)
  intersection_list=[]
  for i in Sweep:
    abc=i.intersection(poly_points)  # abscissa
    x=0
    lines=[]
    points=[]
    if type(abc)==gp.collection.GeometryCollection:
      lines.append((abc[0].boundary[0].y, abc[0].boundary[1].y))
      points.append(abc[1].y)
      x=abc[0].boundary[0].x
    elif type(abc)==gp.multilinestring.MultiLineString:
      lines.append((abc[0].boundary[0].y, abc[0].boundary[1].y))
      lines.append((abc[1].boundary[0].y, abc[1].boundary[1].y))
      x=abc[0].boundary[0].x
    elif type(abc)==gp.linestring.LineString:
      lines.append((abc.boundary[0].y, abc.boundary[1].y))
      x=abc.boundary[0].x
    elif type(abc)==gp.point.Point:
      points.append(abc.y)
      x=abc.x
    lst=[x,points,lines]
    intersection_list.append(lst)
  intersection_list = sorted(intersection_list, key=lambda x: x[0])
  # print()
  # for i in intersection_list:
  #   print(i)
  # print()

  final_poly=[]
  first=[(intersection_list[0][0], intersection_list[0][1][0]), (intersection_list[1][0], intersection_list[1][2][0][0]), (intersection_list[1][0], intersection_list[1][2][0][1]), (intersection_list[0][0], intersection_list[0][1][0])]
  final_poly.append(first)
  for i in range(1, len(intersection_list)-2):
    x1,p1,l1=intersection_list[i]
    x2,p2,l2=intersection_list[i+1]
    if (len(l2)>len(l1) and len(l2)==2 and  l2[0][1]==l2[1][0]):
      l2=[(l2[0][0],l2[1][1])]
    elif (len(l1)>len(l2) and len(l1)==2 and l1[0][1]==l1[1][0]):
      l1=[(l1[0][0],l1[1][1])]
    elif (len(l1)==2 and len(l2)==2 and l2[0][1]==l2[1][0] and l1[0][1]==l1[1][0]):
      l1=[(l1[0][0],l1[1][1])]
      l2=[(l2[0][0],l2[1][1])]

    mid1=[]
    if(p1!=[]):
      mid1.append((p1[0],0))
    if(l1!=[]):
      for i in l1:
        mp=(i[0]+i[1])/2
        mid1.append((mp,1))

    mid2=[]
    if(p2!=[]):
      mid2.append((p2[0],0))
    if(l2!=[]):
      for i in l2:
        mp=(i[0]+i[1])/2
        mid2.append((mp,1))
    # print(mid1,mid2)

    for i in mid1:
      for j in mid2:
        points=[]
        midpoint_joint = gp.LineString([gp.Point(x1,i[0]),gp.Point(x2,j[0])])
        inter=midpoint_joint.intersection(poly_points)
        
        # print(type(inter)==gp.linestring.LineString)
        
        if type(inter)==gp.linestring.LineString:

          if i[1]==0 and j[1]==1:
            points.append((x1,i[0]))
            for y1,y2 in l2:
              if(j[0]==(y1+y2)/2):
                points.append((x2,y1))
                points.append((x2,y2))
            points.append((x1,i[0]))

          elif i[1]==1 and j[1]==0:
            points.append((x2,j[0]))
            for y1,y2 in l1:
              if(i[0]==(y1+y2)/2):
                points.append((x1,y1))
                points.append((x1,y2))
            points.append((x2,j[0]))
            
          elif i[1]==1 and j[1]==1:
            first=()
            for y1,y2 in l1:
              if(i[0]==(y1+y2)/2):
                first=(x1,y1)
                points.append((x1,y1))
                points.append((x1,y2))
            for y1,y2 in l2:
              if(j[0]==(y1+y2)/2):
                points.append((x2,y2))
                points.append((x2,y1))
            points.append(first)


        if len(points)>0:
          final_poly.append(points)

        
  n=len(intersection_list)
  last=[(intersection_list[n-1][0], intersection_list[n-1][1][0]), (intersection_list[n-2][0], intersection_list[n-2][2][0][0]), (intersection_list[n-2][0], intersection_list[n-2][2][0][1]), (intersection_list[n-1][0], intersection_list[n-1][1][0])]
  final_poly.append(last)

  return final_poly
  

# grid based sweep

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
    # sweep_start_posi=[400,-600]
    # print(sweep_vec, sweep_start_posi)

    rox, roy = convert_grid_coordinate(ox, oy, sweep_vec, sweep_start_posi)

    gmap, xinds_goaly, goaly = setup_grid_map(rox, roy, reso, sweeping_direction)

    sweep_searcher = SweepSearcher(moving_direction, sweeping_direction, xinds_goaly, goaly)

    px, py = sweep_path_search(sweep_searcher, gmap)

    rx, ry = convert_global_coordinate(px, py, sweep_vec, sweep_start_posi)

    print("Path length:", len(rx))

    return rx, ry
  
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
    # sweep_start_pos=[0, 0]
    vec=[0,1]
    return vec, sweep_start_pos

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

def smooth(kml_str):
  array=kml_str.split(" ")
  ox=[]
  oy=[]
  for i in range(len(array)):
      lst=array[i].split(",")
      x=(float(lst[0])-80.29821599780011)*100000
      y=(float(lst[1])-13.13286087620951)*100000
      ox.append(x)
      oy.append(y)
  try_x=[]
  try_y=[]

  try_x.append(ox[0])
  try_y.append(oy[0])

  for i in range(len(array)-2):
      temp1 = arctan((oy[i+1]-oy[i])/(ox[i+1]-ox[i])) *180/pi
      temp2 = arctan((oy[i+2]-oy[i+1])/(ox[i+2]-ox[i+1])) *180/pi
      temp = temp2 - temp1
      if abs(temp)>25:
          try_x.append(ox[i+1])
          try_y.append(oy[i+1])

  try_x.append(ox[i+1])
  try_y.append(oy[i+1])

  try_x.append(ox[i+2])
  try_y.append(oy[i+2])

  points=[]
  for i in range(len(try_x)):
      points.append([try_x[i],try_y[i]])
  return points

def distance(poly,cell1,cell2):
  cell1=cell1.copy()
  cell2=cell2.copy()
  cell1.pop(-1)
  cell2.pop(-1)
  dist=0

  cx1,cy1=0,0
  n1=len(cell1)
  for x,y in cell1:
    cx1+=x
    cy1+=y
  cx1=cx1/n1
  cy1=cy1/n1

  cx2,cy2=0,0
  n2=len(cell2)
  for x,y in cell2:
    cx2+=x
    cy2+=y
  cx2=cx2/n2
  cy2=cy2/n2

  poly_points=gp.Polygon(poly)
  midpoint_joint = gp.LineString([gp.Point(cx1,cy1),gp.Point(cx2,cy2)])
  inter=midpoint_joint.intersection(poly_points)
  
  flag = type(inter)==gp.linestring.LineString
  if (flag==False):
    dist=Inf
  else:
    dist=sqrt((cx2-cx1)**2 + (cy2-cy1)**2)

  return dist

# Gives the distance between two nodes
def distance_callback(from_index, to_index):
    
    # Convert routing variable Index to distance matrix NodeIndex.
    from_node = manager.IndexToNode(from_index)
    to_node = manager.IndexToNode(to_index)
    return data['distance_matrix'][from_node][to_node]

# Displays the solution to the output
def print_solution(manager, routing, solution):
    lst=[0]
    print('Objective: {} miles'.format(solution.ObjectiveValue()))
    index = routing.Start(0)
    
    plan_output = 'Best route found for vehicle 0:\n'
    route_distance = 0
    
    while not routing.IsEnd(index):
        plan_output += ' {} ->'.format(manager.IndexToNode(index))
        previous_index = index
        index = solution.Value(routing.NextVar(index))
        route_distance += routing.GetArcCostForVehicle(previous_index, index, 0)
        lst.append(manager.IndexToNode(index))
    plan_output += ' {}\n'.format(manager.IndexToNode(index))
    # lst.append(manager.IndexToNode(index))
    
    print(plan_output)
    plan_output += 'Total distance for route: {}miles\n'.format(route_distance)
    return lst



kml1="80.29821599780011,13.13286087620951,0 80.29796720705747,13.13182773034574,0 80.2991076628975,13.13169660518675,0 80.29944038024139,13.13165299234652,0 80.29937762620091,13.13109111234745,0 80.29793726328148,13.13128457469073,0 80.29761056937502,13.13006729251901,0 80.29806252902425,13.12999465932925,0 80.29835086427846,13.13024456717143,0 80.29941213185222,13.12998434135773,0 80.30023518680595,13.12965420556122,0 80.30044586533064,13.12906413112875,0 80.30050201356336,13.12866328487464,0 80.30036109015327,13.12804180049814,0 80.30025662226946,13.12776937285265,0 80.3000525602131,13.12752807411891,0 80.29968951597387,13.12732678106736,0 80.29850101748954,13.12753770530468,0 80.29819061015509,13.12698463906859,0 80.29826246128076,13.12663847836676,0 80.29863060875304,13.12616498023377,0 80.29866173571808,13.12568598023374,0 80.29872099323003,13.12514727753639,0 80.29866665209987,13.1245795774647,0 80.29799448684999,13.12443297917521,0 80.29745577682087,13.1246465082633,0 80.29777843797754,13.12406362090931,0 80.29901796325314,13.12452616431568,0 80.2994215614677,13.12418173688713,0 80.29890539382183,13.1237513033455,0 80.29843763117461,13.12349401672662,0 80.29966445016528,13.12313628366699,0 80.30020583464541,13.12344717000246,0 80.30023942717858,13.12346802009502,0 80.30066050587938,13.12375513411228,0 80.30091433792411,13.12427053640675,0 80.30122652457746,13.12483053454407,0 80.30162645319366,13.12545180733082,0 80.30201672601071,13.12617996201633,0 80.30261722067974,13.12696180815264,0 80.30289554473505,13.12760244671527,0 80.30330431624785,13.12807993910696,0 80.30333822855845,13.12813887406227,0 80.30371150284772,13.12883379102835,0 80.30402529417192,13.12914611634577,0 80.30388023777046,13.12927645834875,0 80.30334712621725,13.12996559499788,0 80.30301378651669,13.12975218051434,0 80.30107548260035,13.13287222952495,0 80.30069077457982,13.13320044898597,0 80.2995692553596,13.1334836836321,0 80.2994531642517,13.13318768783351,0 80.30037261724844,13.13287222535551,0 80.30060620803647,13.13244118435543,0 80.30019561062818,13.1323564180065,0 80.29969301828299,13.13235642075422,0 80.29913484053139,13.13250202156531,0 80.29843887307364,13.13277263561584,0 80.29821599780011,13.13286087620951,0"
kml2="80.29896188177538,13.13457871394119,0 80.30179705144192,13.13393332239194,0 80.30251262471795,13.13258722535745,0 80.3037296832117,13.13018483182888,0 80.30439533014118,13.12995472237192,0 80.30455207341703,13.13097186008446,0 80.30469291859053,13.13205528046478,0 80.3074130726969,13.13881532179428,0 80.30969697138669,13.14515284587348,0 80.31109872275637,13.14839401014196,0 80.31550715366353,13.15925779755686,0 80.31269281527685,13.16050495443757,0 80.31218459464893,13.15938846639837,0 80.3096718932964,13.15980990169297,0 80.30836266443212,13.15689452707756,0 80.3072121392284,13.15451115489642,0 80.30830948366555,13.15412619754567,0 80.30801269745258,13.15313937038605,0 80.3077707179326,13.1526328941917,0 80.3064751044341,13.15288880296986,0 80.30517696832142,13.14955735572043,0 80.30657288759892,13.1490290278613,0 80.30605448482379,13.14784294170603,0 80.30462792442583,13.14825030518204,0 80.30289714996391,13.14510334409799,0 80.30425576758503,13.14446955608774,0 80.30385630492628,13.14330550518714,0 80.3025426440617,13.14327683884104,0 80.30063763636079,13.14084518677255,0 80.30074407799306,13.14083845736187,0 80.30250238236535,13.14030651352447,0 80.3022178732283,13.13946864712285,0 80.30060454392019,13.14006144559007,0 80.29896188177538,13.13457871394119,0"
smoothened=smooth(kml2)
PointList=[tuple(l) for l in smoothened]

#AUV specs
AUV_length=0.2
AUV_breadth=0.1
AUV_Loc=[0,0.3]

poly= gp.Polygon(PointList)

a,b,c,d=operatePolygon(PointList)

ox=[]
oy=[]
for i in a:
    ox.append(i[0])
    oy.append(i[1])

plt.plot(ox,oy)

final_poly=CellDecomposition(a)
plt.show()

# for i in final_poly:
#   print(i)

n=len(final_poly)

distance_matrix=[]
for i in range(n+1):
  lst=[]
  for j in range(n+1):
    if i==n and j==n:
      dis=0
    elif i==n:
      if j==n-1:
        dis=0
      else:
        dis=Inf
    elif j==n:
      if i==n-1:
        dis=0
      else:
        dis=Inf
    else:
      dis=distance(a, final_poly[i], final_poly[j])
    lst.append(dis)
  distance_matrix.append(lst)

# print()
# for i in distance_matrix:
#   print(i)

# solving Travelling Salesman Problem

data={}
data['distance_matrix'] = distance_matrix
data['no_of_vehicles'] = 1
data['depot'] = 0

manager = pywrapcp.RoutingIndexManager(
    len(data['distance_matrix']), 
    data['no_of_vehicles'],
    data['depot'])

routing = pywrapcp.RoutingModel(manager)

transit_callback_index = routing.RegisterTransitCallback(distance_callback)

routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

search_parameters = pywrapcp.DefaultRoutingSearchParameters()
search_parameters.first_solution_strategy = (
    routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)

# Call the SolveWithParameters() and print the solution
solution = routing.SolveWithParameters(search_parameters)
# print()
# print(solution)
lst=[]
if solution:
    lst=print_solution(manager, routing, solution)

lst=lst[1:]
print(lst)
pxs,pys=[],[]
oxs,oys=[],[]
print(len(final_poly))
for i in lst:
  if i<n:
    ox=[]
    oy=[]
    for x,y in final_poly[i]:
      ox.append(x)
      oy.append(y)
    reso=(max(ox)-min(ox))/10
    if reso>2:
      reso=2
    if reso>0.01:
      px, py = planning(ox, oy, reso)
      pxs.append(px)
      pys.append(py)
      oxs.append(ox)
      oys.append(oy)
    plt.close("all")
for i in range(len(pxs)):
  plt.plot(oxs[i], oys[i], "-xb")
  plt.plot(pxs[i], pys[i],"-r")
  plt.pause(0.1)

plt.show()