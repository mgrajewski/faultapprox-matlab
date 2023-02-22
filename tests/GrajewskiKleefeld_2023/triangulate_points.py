# -*- coding: utf-8 -*-
"""
@author: Matthias Grajewski
This file is part of faultapprox-matlab
(https://github.com/mgrajewski/faultapprox-matlab)

"""
import numpy as np
from numpy import genfromtxt
import open3d as o3d
import csv

name = "Peanut"

csv_file = "./raw_results/points" + name + ".txt"
csv_file_normals = "./raw_results/normals" + name + ".txt"

# read data fromn the csv-files
points = genfromtxt(csv_file, delimiter=',')
normals = genfromtxt(csv_file_normals, delimiter=',')

pcloud = o3d.geometry.PointCloud()

pcloud.points = o3d.utility.Vector3dVector(points)
pcloud.normals = o3d.utility.Vector3dVector(normals)

distances = pcloud.compute_nearest_neighbor_distance()

distances = np.array(distances)

# heuristic
distances = 3.0*distances

mesh_ball = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcloud, \
                                                                            o3d.utility.DoubleVector(distances))

o3d.io.write_triangle_mesh('raw_results/' + name + '.ply', mesh_ball)
