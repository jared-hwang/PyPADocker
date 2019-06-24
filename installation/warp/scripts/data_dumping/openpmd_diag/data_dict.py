"""
This file defines useful correspondance dictionaries
which are used in the openPMD writer
"""
import numpy as np

# Correspondance between quantity and corresponding dimensions
# As specified in the openPMD standard, the arrays represent the
# 7 basis dimensions L, M, T, I, theta, N, J
unit_dimension_dict = {
    "rho" : np.array([-3., 0., 1., 1., 0., 0., 0.]),
    "J" : np.array([-2., 0., 0., 1., 0., 0., 0.]),
    "E" : np.array([ 1., 1.,-3.,-1., 0., 0., 0.]),
    "B" : np.array([ 0., 1.,-2.,-1., 0., 0., 0.]),
    "charge" : np.array([0., 0., 1., 1., 0., 0., 0.]),
    "mass" : np.array([1., 0., 0., 0., 0., 0., 0.]),
    "weighting" : np.array([0., 0., 0., 0., 0., 0., 0.]),
    "id"        : np.array([0., 0., 0., 0., 0., 0., 0.]),
    "position" : np.array([1., 0., 0., 0., 0., 0., 0.]),
    "positionOffset" : np.array([1., 0., 0., 0., 0., 0., 0.]),
    "momentum" : np.array([1., 1.,-1., 0., 0., 0., 0.]),
    "E" : np.array([ 1., 1.,-3.,-1., 0., 0., 0.]),
    "B" : np.array([ 0., 1.,-2.,-1., 0., 0., 0.]),
    "t" : np.array([ 0., 0., 1., 0., 0., 0., 0.])}

# Spatial offset of the different fields
x_offset_dict = {
    "rho":0., "Et":0., "Ey":0., "Ez":0., "Jt":0., "Jy":0., "Jz":0.,
    "Br":0., "Bx":0., "Er":0.5, "Ex":0.5, "Jr":0.5, "Jx":0.5, "Bt":0.5,
    "By":0.5, "Bz":0.5 }
y_offset_dict = {
    "rho":0., "Ex":0., "Ez":0., "Jx":0., "Jz":0., "By":0.,
    "Ey":0.5, "Jy":0.5, "Bx":0.5, "Bz":0.5  }
z_offset_dict = {
    "rho":0., "Er":0., "Ex":0., "Et":0., "Ey":0.,
    "Jr":0., "Jx":0., "Jt":0., "Jy":0., "Bz":0.,
    "Ez":0.5, "Jz":0.5, "Br":0.5, "Bx":0.5, "Bt":0.5, "By":0.5 }

# Typical weighting of different particle properties
macro_weighted_dict = {
    "charge": np.uint32(0),
    "mass": np.uint32(0),
    "weighting": np.uint32(1),
    "id": np.uint32(0),
    "position": np.uint32(0),
    "positionOffset": np.uint32(0),
    "momentum": np.uint32(0),
    "E": np.uint32(0),
    "B": np.uint32(0),
    "t": np.uint32(0) }

weighting_power_dict = {
    "charge": 1.,
    "mass": 1.,
    "weighting": 1.,
    "id": 0.,
    "position": 0.,
    "positionOffset": 0.,
    "momentum": 1.,
    "E": 0.,
    "B": 0.,
    "t": 0.}

# Correspondance between the names in OpenPMD and the names in Warp
circ_dict_quantity = { 'rho':'Rho', 'Er':'Ex', 'Et':'Ey', 'Ez':'Ez',
                        'Br':'Bx', 'Bt':'By', 'Bz':'Bz',
                        'Jr':'Jx', 'Jt':'Jy', 'Jz':'Jz' }
cart_dict_quantity = { 'rho':'Rho', 'Ex':'Ex', 'Ey':'Ey', 'Ez':'Ez',
                        'Bx':'Bx', 'By':'By', 'Bz':'Bz',
                        'Jx':'Jx', 'Jy':'Jy', 'Jz':'Jz' }

# Correspondance between openPMD path and short-hand names
particle_quantity_dict = { 'E' : 'e', 'B' : 'b', \
                           'position' : '', 'momentum' : 'u' }

# Correspondance between the boundary conditions in Warp,
# and the corresponding representative integer
field_boundary_dict = {
    0: np.string_("reflecting"),
    1: np.string_("reflecting"),
    2: np.string_("periodic"),
    3: np.string_("openbc") }
particle_boundary_dict = {
    0: np.string_("absorbing"),
    1: np.string_("reflecting"),
    2: np.string_("periodic") }
# Correspondance between the field solver in Warp,
# and the corresponding representative integer
field_solver_dict = {
    0: np.string_("Yee"),
    1: np.string_("CK"),
    3: np.string_("Lehe") }
