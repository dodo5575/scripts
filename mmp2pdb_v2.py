#!/usr/bin/env python3

# author: Chen-Yu Li cli56@illinois.edu

import sys, optparse, json, os, glob, re, random
from subprocess import Popen, PIPE
from math import *
import numpy as np

###############################################
## subroutines 
###############################################


def alignBaseNormal(pdb_dict, seq, tName_array, unitDestVec_array, direction):

    tail_center = measure_center_by_name(pdb_dict, tName_array)
    BaseNormal_array = unitVec(getBaseNormal(pdb_dict, seq))
    
    #rotate_matrix = rotateMatrix_2(unitDestVec_array, BaseNormal_array)
    
    if direction == 'x':
        rotate_matrix = rotateMatrix_x(unitDestVec_array, BaseNormal_array)
    elif direction == 'y':
        rotate_matrix = rotateMatrix_y(unitDestVec_array, BaseNormal_array)
    elif direction == 'z':
        rotate_matrix = rotateMatrix_z(unitDestVec_array, BaseNormal_array)

    return rotateby_nt(pdb_dict, rotate_matrix, tail_center)


def alignNt(pdb_dict, hName_array, tName_array, unitDestVec_array, direction):

    head_center = measure_center_by_name(pdb_dict, hName_array)
    tail_center = measure_center_by_name(pdb_dict, tName_array)
    unitOriVec_array = unitVec(head_center - tail_center)

    #print(unitDestVec_array, unitOriVec_array)
    if direction == 'x':
        rotate_matrix = rotateMatrix_x(unitDestVec_array, unitOriVec_array)
    elif direction == 'y':
        rotate_matrix = rotateMatrix_y(unitDestVec_array, unitOriVec_array)
    elif direction == 'z':
        rotate_matrix = rotateMatrix_z(unitDestVec_array, unitOriVec_array)
    elif direction == '2':
        rotate_matrix = rotateMatrix_2(unitDestVec_array, unitOriVec_array)
    else:
        raise Exception("alignNt(pdb_dict, hName_array, tName_array, destVec_array, direction): Unrecognized direction %s" % direction)

    return rotateby_nt(pdb_dict, rotate_matrix, tail_center)


def rotation_bp(Seg, k, strand_1, ind_1, seq_1, strand_2, ind_2, seq_2, bp_m):

    global data, anchor_atoms_bp

    # align base-pair
    center = measure_center_bp(bp_m)
    head = data[strand_1][ind_1]['xyz']
    dest_vec = unitVec(head - center)

    group = np.array([])
    for i in anchor_atoms_bp:
        group = np.append(group, bp_m[seq_1][i]) 
    group = np.resize(group, (len(group)/3, 3))    
    anchor = measure_center(group)
    unitOriVec_array = unitVec(anchor - center)

    rotate_axis = (dest_vec + unitOriVec_array) / 2

    rotate_matrix = rotateMatrix_by_axis(rotate_axis, pi) 
    bp_r = rotateby_bp(bp_m, rotate_matrix, center) 

    ## align base-pair normal
    if k + 1 in data[Seg]:
        head = data[Seg][k]['xyz']
        tail = data[Seg][k + 1]['xyz']
    elif k - 1 in data[Seg]:   
        head = data[Seg][k - 1]['xyz']
        tail = data[Seg][k]['xyz']

    axis_vec = unitVec(head - tail)

    center = measure_center_bp(bp_r)
    v1 = bp_r[seq_1][" C1'"] - center
    v2 = bp_r[seq_2][" C1'"] - center
    # angle of v1, v2 > 180
    bp_normal = unitVec(-1 * np.cross(v1, v2)) 

    theta = acos(np.dot(axis_vec, bp_normal))

    group = np.array([])
    for i in anchor_atoms_bp:
        group = np.append(group, bp_r[seq_1][i]) 
    group = np.resize(group, (len(group)/3, 3))    
    anchor = measure_center(group)
    rotate_axis = center - anchor 

    rotate_matrix = rotateMatrix_by_axis(rotate_axis, theta) 
    bp_r = rotateby_bp(bp_r, rotate_matrix, center) 

    # If the direction is wrong
    center = measure_center_bp(bp_r)
    v1 = bp_r[seq_1][" C1'"] - center
    v2 = bp_r[seq_2][" C1'"] - center
    bp_normal = unitVec(-1 * np.cross(v1, v2)) 
    if acos(np.dot(axis_vec, bp_normal)) > (15 * pi / 180):
        rotate_matrix = rotateMatrix_by_axis(rotate_axis, -1 * theta) 
        bp_r = rotateby_bp(bp_r, rotate_matrix, center) 
        bp_r = rotateby_bp(bp_r, rotate_matrix, center) 

    return bp_r

def rotation_nt(group, rid, pdb):

    global data, anchor_atoms_nt 
    
    direc = data[group][rid]['direc']

    seq = data[group][rid]['seq']

    central = data[group][rid]['bond1'][0]
    for key in data.keys():
        if central in data[key]:
            head = data[key][central]['xyz']
    
    dest_vec = unitVec(head - data[group][rid]['xyz'])

    
    ## alignNt
    head_center = measure_center_by_name(pdb, base_atoms(seq))
    tail_center = measure_center_by_name(pdb, anchor_atoms_nt)
    unitOriVec_array = unitVec(head_center - tail_center)

    rotate_axis = (dest_vec + unitOriVec_array) / 2

    rotate_matrix = rotateMatrix_by_axis(rotate_axis, pi) 
    pdb_r = rotateby_nt(pdb, rotate_matrix, tail_center) 

    #pdb_r = alignNt(pdb_r, base_atoms(seq), anchor_atoms_nt, dest_vec, 'z') 

    #while np.linalg.norm(dest_vec - unitOriVec_array) > 0.1: 
    #    pdb_r = alignNt(pdb_r, base_atoms(seq), anchor_atoms_nt, dest_vec, 'x') 
    #    head_center = measure_center_by_name(pdb_r, base_atoms(seq))
    #    tail_center = measure_center_by_name(pdb_r, anchor_atoms_nt)
    #    unitOriVec_array = unitVec(head_center - tail_center)


    ## alignBaseNormal
    for key in data.keys():
        central = data[group][rid]['bond1'][0]
        if (central in data[key]) and (central + 1 in data[key]):
            head = data[key][central]['xyz']
            tail = data[key][central + 1]['xyz']
        elif central in data[key] and central - 1 in data[key]:   
            head = data[key][central - 1]['xyz']
            tail = data[key][central]['xyz']
#        else:
#            raise Exception("rotation(group, rid, moved_pdb): Cannot define axis direction:\ngroup = %s\nrid = %s\ncentral = %s" % (group, rid, central))

    axis_vec = unitVec(head - tail)

    BaseNormal_array = direc * unitVec(getBaseNormal(pdb_r, seq))

    theta = acos(np.dot(axis_vec, BaseNormal_array))
    
    head_center = measure_center_by_name(pdb_r, base_atoms(seq))
    tail_center = measure_center_by_name(pdb_r, anchor_atoms_nt)
    nt_axis = head_center - tail_center

    rotate_matrix = rotateMatrix_by_axis(nt_axis, theta) 
    pdb_r = rotateby_nt(pdb_r, rotate_matrix, tail_center) 

    BaseNormal_array = direc * unitVec(getBaseNormal(pdb_r, seq))
    # If the direction is wrong
    if acos(np.dot(axis_vec, BaseNormal_array)) > (15 * pi / 180):
        rotate_matrix = rotateMatrix_by_axis(nt_axis, -1 * theta) 
        pdb_r = rotateby_nt(pdb_r, rotate_matrix, tail_center) 
        pdb_r = rotateby_nt(pdb_r, rotate_matrix, tail_center) 


    return pdb_r


def translation_bp(Segment, rid, bp_dict):

    global data, bp

    center = measure_center_bp(bp_dict)

    move_vec = data[Segment][rid]['xyz'] - center 
    bp_m = moveby_bp(bp_dict, move_vec)

    return bp_m


def translation_nt(group, rid):

    global data, anchor_atoms_nt

    anchor_atoms_nt_pdb = np.array([])
    seq = data[group][rid]['seq']
    for i in anchor_atoms_nt:
        anchor_atoms_nt_pdb = np.append(anchor_atoms_nt_pdb, pdb[seq][i], axis=0)
    anchor_center = measure_center(anchor_atoms_nt_pdb)
    move_vec = data[group][rid]['xyz'] - anchor_center 
    moved_pdb = moveby_nt(pdb[seq], move_vec)

    return moved_pdb


def getBaseNormal(coor_dict, nt):

    if nt == 'A' or nt == 'G':
    # N1 N3 x N1 C4
        BaseNormal = np.cross((coor_dict[" N3 "] - coor_dict[" N1 "]), (coor_dict[" C4 "] - coor_dict[" N1 "]))
        return BaseNormal

    elif nt == 'T' or nt == 'C':
    # N1 C4 x N1 N3
        BaseNormal = np.cross((coor_dict[" C4 "] - coor_dict[" N1 "]), (coor_dict[" N3 "] - coor_dict[" N1 "]))
        return BaseNormal

    else:
        raise Exception("getBaseNormal(coor_dict, nt): Unrecognized nucleotide %s" % nt)



def rotateMatrix_x(b, a):

    R = np.array([[1,0,0],
                  [0, np.dot(a, b), -1 * np.linalg.norm(np.cross(a,b))],
                  [0, np.linalg.norm(np.cross(a,b)), np.dot(a, b)]])

    return R


def rotateMatrix_y(b, a):

    R = np.array([[np.dot(a, b), 0, np.linalg.norm(np.cross(a,b))],
                  [0,1,0],
                  [-1 * np.linalg.norm(np.cross(a,b)), 0, np.dot(a, b)]])

    return R


def rotateMatrix_z(b, a):

    R = np.array([[np.dot(a, b), -1 * np.linalg.norm(np.cross(a,b)), 0],
                  [np.linalg.norm(np.cross(a,b)), np.dot(a, b), 0],
                  [0,0,1]])

    return R


def rotateMatrix_1(b, a):

    R = np.array([unitVec(a),
                  unitVec(np.cross(np.cross(a,b), a)),
                  unitVec(np.cross(a,b))])

    return R


def rotateMatrix_2(b, a):

    v = unitVec(np.cross(a,b))

    vx = np.array([[0, -1 * v[2], v[1]],
                   [v[2], 0, -1 * v[0]],
                   [-1 * v[1], v[0], 0]])
    
    #angle = vecAngle(a, b)
    #angle = acos(np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b)))

    #s = np.linalg.norm(v) * sin(angle)
    #c = np.dot(a,b) * cos(angle)
    #R = np.identity(3) + vx + ((1 - c) / s**2) * np.dot(vx, vx)
    sinTheta = np.cross(a,b) / (np.linalg.norm(a) * np.linalg.norm(b)) 
    cosTheta = np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b)) 

    R = np.identity(3) + sinTheta * vx + (1 - cosTheta) * np.dot(vx, vx)

    return R

def rotateMatrix_by_axis(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    Numerically Accurate!!
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = unitVec(axis)
    a = cos(theta / 2)
    b, c, d = -axis * sin(theta / 2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def vecAngle(l1, l2):

    return acos(vecdot(l1, l2) / (veclength(l1) * veclength(l2)))  


def nt1to3(nt):

    if nt == 'A' or nt == 'a':

        return 'ADE'

    elif nt == 'G' or nt == 'g':

        return 'GUA'

    elif nt == 'T' or nt == 't':

        return 'THY'

    elif nt == 'C' or nt == 'c':

        return 'CYT'

    else:
        raise Exception("nt1to3(nt): Unrecognized nucleotide %s" % nt)


def base_atoms(nt):

    if nt == 'A' or nt == 'a':

        return np.array([" N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 ", " N7 ", " C8 ", " N9 "])
        #return np.array([" C4 "])

    elif nt == 'G' or nt == 'g':

        return np.array([" N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 ", " N7 ", " C8 ", " N9 "])

    elif nt == 'T' or nt == 't':

        return np.array([" N1 "])
        #return np.array([" N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 "])

    elif nt == 'C' or nt == 'c':

        return np.array([" N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 "])

    else:
        raise Exception("base_atoms(nt): Unrecognized nucleotide %s" % nt)



def rotateby_bp(ori_dict, rotate_matrix, anchor_array):

    dest_dict = {}

    for key in ori_dict.keys():

        dest_dict[key] = {}
    
        for k, v in ori_dict[key].items():

            dest_dict[key][k] = np.dot(rotate_matrix, (v - anchor_array)) + anchor_array 

    return dest_dict


def rotateby_nt(ori_dict, rotate_matrix, anchor_array):

    dest_dict = {}

    for k, v in ori_dict.items():

        dest_dict[k] = np.dot(rotate_matrix, (v - anchor_array)) + anchor_array 

    return dest_dict


def moveby_bp(ori_dict, vec_array):

    dest_dict = {}

    for key in ori_dict.keys():
        dest_dict[key] = {}

        for k, v in ori_dict[key].items():

            dest_dict[key][k] = v + vec_array

    return dest_dict


def moveby_nt(ori_dict, vec_array):

    dest_dict = {}

    for k, v in ori_dict.items():

        dest_dict[k] = v + vec_array

    return dest_dict


def unitVec(l):

    vec = np.array([])
    vlength = veclength(l)

    for i in l:
        vec = np.append(vec, i / vlength)

    return vec


def vecadd(l1, l2):

    return [x + y for x, y in zip(l1,l2)]


def vecsub(l1, l2):

    return [x - y for x, y in zip(l1,l2)]


def vecdot(l1, l2):

    return sum([x * y for x, y in zip(l1,l2)])


def veccross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return c


def veclength(l):

    length = 0
    for i in l:
        length += i**2

    length = sqrt(length)

    return length


def measure_center_bp(bp_dict):
 
    group = np.array([])
    for nt in bp_dict.keys():
        for atom in bp_dict[nt].keys():
            group = np.append(group, bp_dict[nt][atom]) 
    
    group = np.resize(group, (len(group)/3, 3))
    group_center = measure_center(group)

    return group_center


def measure_center_by_name(pdb_dict, name_array):
 
    group = np.array([])
    for i in name_array:
        group = np.append(group, pdb_dict[i]) 
    
    group = np.resize(group, (len(group)/3, 3))
    group_center = measure_center(group)

    return group_center


def measure_center(xyz_array):

    (x, y, z) = (0, 0, 0)

    n = len(xyz_array)

    if xyz_array.ndim == 2:

        for xyz in xyz_array:
            x += xyz[0]
            y += xyz[1]
            z += xyz[2]
    
        return np.array([x / n, y / n, z / n])

    elif xyz_array.ndim == 1:

        return xyz_array

    else:
        raise Exception("measure_center(xyz_array): Unsupported dimension %s" % xyz_array.ndim)



def init_bp():
    bp = {
    "AT": {
    "A": {
            " P  ": np.array([   1.159,  9.516, -1.474 ]),
	    " O1P": np.array([   1.187, 10.852, -2.109 ]),
	    " O2P": np.array([   2.020,  9.346, -0.281 ]),
	    " O5'": np.array([  -0.348,  9.126, -1.108 ]),
	    " C5'": np.array([  -1.169,  8.499, -2.112 ]),
	    " H5'": np.array([  -0.578,  8.081, -2.802 ]),
	    "H5''": np.array([  -1.861,  9.152, -2.421 ]),
	    " C4'": np.array([  -1.947,  7.349, -1.501 ]),
	    " H4'": np.array([  -2.791,  7.245, -2.028 ]),
	    " O4'": np.array([  -1.206,  6.096, -1.540 ]),
	    " C1'": np.array([  -1.269,  5.500, -0.253 ]),
	    " H1'": np.array([  -2.135,  5.004, -0.194 ]),
	    " N9 ": np.array([  -0.105,  4.633, -0.098 ]),
	    " C5 ": np.array([   1.216,  2.883,  0.078 ]),
	    " N7 ": np.array([   2.034,  4.001,  0.187 ]),
	    " C8 ": np.array([   1.205,  5.010,  0.076 ]),
	    " H8 ": np.array([   1.721,  5.864,  0.144 ]),
	    " N1 ": np.array([   0.448,  0.657, -0.023 ]),
	    " C2 ": np.array([  -0.777,  1.168, -0.186 ]),
	    " H2 ": np.array([  -1.580,  0.581, -0.293 ]),
	    " N3 ": np.array([  -1.161,  2.442, -0.237 ]),
	    " C4 ": np.array([  -0.104,  3.258, -0.097 ]),
	    " C6 ": np.array([   1.490,  1.506,  0.115 ]),
	    " N6 ": np.array([   2.713,  0.996,  0.278 ]),
	    " H61": np.array([   2.827, -0.007,  0.294 ]),
	    " H62": np.array([   3.496,  1.625,  0.382 ]),
	    " C2'": np.array([  -1.298,  6.667,  0.728 ]),
	    "H2''": np.array([  -0.423,  7.151,  0.736 ]),
	    " H2'": np.array([  -1.688,  6.389,  1.606 ]),
	    " C3'": np.array([  -2.326,  7.509, -0.030 ]),
	    " H3'": np.array([  -2.195,  8.463,  0.238 ]),
	    " O3'": np.array([  -3.615,  6.944,  0.179 ])},
    "T": {  
	    " P  ": np.array([   1.127, -9.197,  1.503 ]),
	    " O1P": np.array([   1.154,-10.535,  2.137 ]),
	    " O2P": np.array([   1.987, -9.026,  0.311 ]),
	    " O5'": np.array([  -0.380, -8.807,  1.138 ]),
	    " C5'": np.array([  -1.202, -8.180,  2.142 ]),
	    " H5'": np.array([  -0.611, -7.763,  2.832 ]),
	    "H5''": np.array([  -1.894, -8.833,  2.450 ]),
	    " C4'": np.array([  -1.980, -7.030,  1.534 ]),
	    " H4'": np.array([  -2.824, -6.925,  2.060 ]),
	    " O4'": np.array([  -1.238, -5.777,  1.573 ]),
	    " C1'": np.array([  -1.301, -5.180,  0.287 ]),
	    " H1'": np.array([  -2.167, -4.683,  0.229 ]),
	    " N1 ": np.array([  -0.129, -4.306,  0.132 ]),
	    " C6 ": np.array([   1.129, -4.843, -0.036 ]),
	    " H6 ": np.array([   1.376, -5.812, -0.070 ]),
	    " C2 ": np.array([  -0.327, -2.945,  0.159 ]),
	    " O2 ": np.array([  -1.426, -2.436,  0.305 ]),
	    " N3 ": np.array([   0.815, -2.198,  0.007 ]),
	    " H3 ": np.array([   0.684, -1.116,  0.024 ]),
	    " C4 ": np.array([   2.103, -2.666, -0.164 ]),
	    " O4 ": np.array([   3.032, -1.873, -0.287 ]),
	    " C5 ": np.array([   2.238, -4.105, -0.182 ]),
	    " C5M": np.array([   3.596, -4.706, -0.364 ]),
	    " H51": np.array([   3.515, -5.703, -0.355 ]),
	    " H52": np.array([   4.130, -4.416,  0.430 ]),
	    " H53": np.array([   4.102, -4.482, -1.197 ]),
	    " C2'": np.array([  -1.331, -6.345, -0.695 ]),
	    "H2''": np.array([  -0.456, -6.829, -0.705 ]),
	    " H2'": np.array([  -1.722, -6.066, -1.572 ]),
	    " C3'": np.array([  -2.358, -7.189,  0.062 ]),
	    " H3'": np.array([  -2.225, -8.143, -0.207 ]),
	    " O3'": np.array([  -3.647, -6.623, -0.145 ])}
    },
    "TA": {
    "T": {  
	    " P  ": np.array([   1.127,  9.197, -1.503 ]),
	    " O1P": np.array([   1.154, 10.535, -2.137 ]),
	    " O2P": np.array([   1.987,  9.026, -0.311 ]),
	    " O5'": np.array([  -0.380,  8.807, -1.138 ]),
	    " C5'": np.array([  -1.202,  8.180, -2.142 ]),
	    " H5'": np.array([  -0.611,  7.763, -2.832 ]),
	    "H5''": np.array([  -1.894,  8.833, -2.450 ]),
	    " C4'": np.array([  -1.980,  7.030, -1.534 ]),
	    " H4'": np.array([  -2.824,  6.925, -2.060 ]),
	    " O4'": np.array([  -1.238,  5.777, -1.573 ]),
	    " C1'": np.array([  -1.301,  5.180, -0.287 ]),
	    " H1'": np.array([  -2.167,  4.683, -0.229 ]),
	    " N1 ": np.array([  -0.129,  4.306, -0.132 ]),
	    " C6 ": np.array([   1.129,  4.843,  0.036 ]),
	    " H6 ": np.array([   1.376,  5.812,  0.070 ]),
	    " C2 ": np.array([  -0.327,  2.945, -0.159 ]),
	    " O2 ": np.array([  -1.426,  2.436, -0.305 ]),
	    " N3 ": np.array([   0.815,  2.198, -0.007 ]),
	    " H3 ": np.array([   0.684,  1.116, -0.024 ]),
	    " C4 ": np.array([   2.103,  2.666,  0.164 ]),
	    " O4 ": np.array([   3.032,  1.873,  0.287 ]),
	    " C5 ": np.array([   2.238,  4.105,  0.182 ]),
	    " C5M": np.array([   3.596,  4.706,  0.364 ]),
	    " H51": np.array([   3.515,  5.703,  0.355 ]),
	    " H52": np.array([   4.130,  4.416, -0.430 ]),
	    " H53": np.array([   4.102,  4.482,  1.197 ]),
	    " C2'": np.array([  -1.331,  6.345,  0.695 ]),
	    "H2''": np.array([  -0.456,  6.829,  0.705 ]),
	    " H2'": np.array([  -1.722,  6.066,  1.572 ]),
	    " C3'": np.array([  -2.358,  7.189, -0.062 ]),
	    " H3'": np.array([  -2.225,  8.143,  0.207 ]),
	    " O3'": np.array([  -3.647,  6.623,  0.145 ])},
    "A": {
            " P  ": np.array([   1.159, -9.516,  1.474 ]),
	    " O1P": np.array([   1.187,-10.852,  2.109 ]),
	    " O2P": np.array([   2.020, -9.346,  0.281 ]),
	    " O5'": np.array([  -0.348, -9.126,  1.108 ]),
	    " C5'": np.array([  -1.169, -8.499,  2.112 ]),
	    " H5'": np.array([  -0.578, -8.081,  2.802 ]),
	    "H5''": np.array([  -1.861, -9.152,  2.421 ]),
	    " C4'": np.array([  -1.947, -7.349,  1.501 ]),
	    " H4'": np.array([  -2.791, -7.245,  2.028 ]),
	    " O4'": np.array([  -1.206, -6.096,  1.540 ]),
	    " C1'": np.array([  -1.269, -5.500,  0.253 ]),
	    " H1'": np.array([  -2.135, -5.004,  0.194 ]),
	    " N9 ": np.array([  -0.105, -4.633,  0.098 ]),
	    " C5 ": np.array([   1.216, -2.883, -0.078 ]),
	    " N7 ": np.array([   2.034, -4.001, -0.187 ]),
	    " C8 ": np.array([   1.205, -5.010, -0.076 ]),
	    " H8 ": np.array([   1.721, -5.864, -0.144 ]),
	    " N1 ": np.array([   0.448, -0.657,  0.023 ]),
	    " C2 ": np.array([  -0.777, -1.168,  0.186 ]),
	    " H2 ": np.array([  -1.580, -0.581,  0.293 ]),
	    " N3 ": np.array([  -1.161, -2.442,  0.237 ]),
	    " C4 ": np.array([  -0.104, -3.258,  0.097 ]),
	    " C6 ": np.array([   1.490, -1.506, -0.115 ]),
	    " N6 ": np.array([   2.713, -0.996, -0.278 ]),
	    " H61": np.array([   2.827,  0.007, -0.294 ]),
	    " H62": np.array([   3.496, -1.625, -0.382 ]),
	    " C2'": np.array([  -1.298, -6.667, -0.728 ]),
	    "H2''": np.array([  -0.423, -7.151, -0.736 ]),
	    " H2'": np.array([  -1.688, -6.389, -1.606 ]),
	    " C3'": np.array([  -2.326, -7.509,  0.030 ]),
	    " H3'": np.array([  -2.195, -8.463, -0.238 ]),
	    " O3'": np.array([  -3.615, -6.944, -0.179 ])}
    },
    "GC": {
    "G": {
	    " P  ": np.array([   1.500,  9.367, -1.357 ]),
	    " O1P": np.array([   1.571, 10.709, -1.977 ]),
	    " O2P": np.array([   2.327,  9.165, -0.146 ]),
	    " O5'": np.array([  -0.025,  9.008, -1.031 ]),
	    " C5'": np.array([  -0.837,  8.410, -2.060 ]),
	    " H5'": np.array([  -0.240,  7.986, -2.741 ]),
	    "H5''": np.array([  -1.507,  9.081, -2.378 ]),
	    " C4'": np.array([  -1.655,  7.272, -1.480 ]),
	    " H4'": np.array([  -2.488,  7.192, -2.027 ]),
	    " O4'": np.array([  -0.941,  6.003, -1.515 ]),
	    " C1'": np.array([  -1.047,  5.396, -0.235 ]),
	    " H1'": np.array([  -1.925,  4.918, -0.202 ]),
	    " N9 ": np.array([   0.092,  4.500, -0.063 ]),
	    " C4 ": np.array([   0.062,  3.126, -0.076 ]),
	    " N2 ": np.array([  -1.719,  0.159, -0.375 ]),
	    " H21": np.array([  -1.620, -0.846, -0.367 ]),
	    " H22": np.array([  -2.652,  0.519, -0.518 ]),
	    " N3 ": np.array([  -1.033,  2.358, -0.252 ]),
	    " C2 ": np.array([  -0.744,  1.068, -0.216 ]),
	    " N1 ": np.array([   0.521,  0.570, -0.020 ]),
	    " H1 ": np.array([   0.631, -0.454, -0.008 ]),
	    " C6 ": np.array([   1.662,  1.342,  0.163 ]),
	    " O6 ": np.array([   2.758,  0.792,  0.331 ]),
	    " C5 ": np.array([   1.367,  2.729,  0.126 ]),
	    " N7 ": np.array([   2.200,  3.832,  0.264 ]),
	    " C8 ": np.array([   1.403,  4.858,  0.145 ]),
	    " H8 ": np.array([   1.938,  5.698,  0.234 ]),
	    " C2'": np.array([  -1.074,  6.552,  0.756 ]),
	    "H2''": np.array([  -0.189,  7.016,  0.790 ]),
	    " H2'": np.array([  -1.491,  6.273,  1.621 ]),
	    " C3'": np.array([  -2.064,  7.426, -0.015 ]),
	    " H3'": np.array([  -1.916,  8.374,  0.266 ]),
	    " O3'": np.array([  -3.370,  6.889,  0.156 ])},
    "C": {
	    " P  ": np.array([   1.069, -9.374,  1.462 ]),
	    " O1P": np.array([   1.043,-10.714,  2.089 ]),
	    " O2P": np.array([   1.961, -9.220,  0.291 ]),
	    " O5'": np.array([  -0.418, -8.939,  1.064 ]),
	    " C5'": np.array([  -1.246, -8.295,  2.053 ]),
	    " H5'": np.array([  -0.660, -7.899,  2.760 ]),
	    "H5''": np.array([  -1.963, -8.930,  2.340 ]),
	    " C4'": np.array([  -1.975, -7.121,  1.433 ]),
	    " H4'": np.array([  -2.827, -6.996,  1.941 ]),
	    " O4'": np.array([  -1.201, -5.889,  1.496 ]),
	    " C1'": np.array([  -1.218, -5.284,  0.213 ]),
	    " H1'": np.array([  -2.068, -4.763,  0.138 ]),
	    " N1 ": np.array([  -0.019, -4.444,  0.089 ]),
	    " C6 ": np.array([   1.226, -4.989, -0.052 ]),
	    " H6 ": np.array([   1.505, -5.948, -0.092 ]),
	    " C5 ": np.array([   2.313, -4.213, -0.163 ]),
	    " H5 ": np.array([   3.289, -4.404, -0.271 ]),
	    " C2 ": np.array([  -0.182, -3.057,  0.120 ]),
	    " O2 ": np.array([  -1.323, -2.589,  0.248 ]),
	    " N3 ": np.array([   0.907, -2.260,  0.009 ]),
	    " C4 ": np.array([   2.118, -2.803, -0.129 ]),
	    " N4 ": np.array([   3.164, -1.979, -0.237 ]),
	    " H41": np.array([   3.005, -0.982, -0.211 ]),
	    " H42": np.array([   4.100, -2.344, -0.343 ]),
	    " C2'": np.array([  -1.258, -6.443, -0.776 ]),
	    "H2''": np.array([  -0.397, -6.951, -0.769 ]),
	    " H2'": np.array([  -1.621, -6.148, -1.660 ]),
	    " C3'": np.array([  -2.324, -7.261, -0.049 ]),
	    " H3'": np.array([  -2.212, -8.217, -0.321 ]),
	    " O3'": np.array([  -3.593, -6.658, -0.282 ])}
    },
    "CG": {
    "C": {
	    " P  ": np.array([   1.069,  9.374, -1.462 ]),
	    " O1P": np.array([   1.043, 10.714, -2.089 ]),
	    " O2P": np.array([   1.961,  9.220, -0.291 ]),
	    " O5'": np.array([  -0.418,  8.939, -1.064 ]),
	    " C5'": np.array([  -1.246,  8.295, -2.053 ]),
	    " H5'": np.array([  -0.660,  7.899, -2.760 ]),
	    "H5''": np.array([  -1.963,  8.930, -2.340 ]),
	    " C4'": np.array([  -1.975,  7.121, -1.433 ]),
	    " H4'": np.array([  -2.827,  6.996, -1.941 ]),
	    " O4'": np.array([  -1.201,  5.889, -1.496 ]),
	    " C1'": np.array([  -1.218,  5.284, -0.213 ]),
	    " H1'": np.array([  -2.068,  4.763, -0.138 ]),
	    " N1 ": np.array([  -0.019,  4.444, -0.089 ]),
	    " C6 ": np.array([   1.226,  4.989,  0.052 ]),
	    " H6 ": np.array([   1.505,  5.948,  0.092 ]),
	    " C5 ": np.array([   2.313,  4.213,  0.163 ]),
	    " H5 ": np.array([   3.289,  4.404,  0.271 ]),
	    " C2 ": np.array([  -0.182,  3.057, -0.120 ]),
	    " O2 ": np.array([  -1.323,  2.589, -0.248 ]),
	    " N3 ": np.array([   0.907,  2.260, -0.009 ]),
	    " C4 ": np.array([   2.118,  2.803,  0.129 ]),
	    " N4 ": np.array([   3.164,  1.979,  0.237 ]),
	    " H41": np.array([   3.005,  0.982,  0.211 ]),
	    " H42": np.array([   4.100,  2.344,  0.343 ]),
	    " C2'": np.array([  -1.258,  6.443,  0.776 ]),
	    "H2''": np.array([  -0.397,  6.951,  0.769 ]),
	    " H2'": np.array([  -1.621,  6.148,  1.660 ]),
	    " C3'": np.array([  -2.324,  7.261,  0.049 ]),
	    " H3'": np.array([  -2.212,  8.217,  0.321 ]),
	    " O3'": np.array([  -3.593,  6.658,  0.282 ])},
    "G": {
	    " P  ": np.array([   1.500, -9.367,  1.357 ]),
	    " O1P": np.array([   1.571,-10.709,  1.977 ]),
	    " O2P": np.array([   2.327, -9.165,  0.146 ]),
	    " O5'": np.array([  -0.025, -9.008,  1.031 ]),
	    " C5'": np.array([  -0.837, -8.410,  2.060 ]),
	    " H5'": np.array([  -0.240, -7.986,  2.741 ]),
	    "H5''": np.array([  -1.507, -9.081,  2.378 ]),
	    " C4'": np.array([  -1.655, -7.272,  1.480 ]),
	    " H4'": np.array([  -2.488, -7.192,  2.027 ]),
	    " O4'": np.array([  -0.941, -6.003,  1.515 ]),
	    " C1'": np.array([  -1.047, -5.396,  0.235 ]),
	    " H1'": np.array([  -1.925, -4.918,  0.202 ]),
	    " N9 ": np.array([   0.092, -4.500,  0.063 ]),
	    " C4 ": np.array([   0.062, -3.126,  0.076 ]),
	    " N2 ": np.array([  -1.719, -0.159,  0.375 ]),
	    " H21": np.array([  -1.620,  0.846,  0.367 ]),
	    " H22": np.array([  -2.652, -0.519,  0.518 ]),
	    " N3 ": np.array([  -1.033, -2.358,  0.252 ]),
	    " C2 ": np.array([  -0.744, -1.068,  0.216 ]),
	    " N1 ": np.array([   0.521, -0.570,  0.020 ]),
	    " H1 ": np.array([   0.631,  0.454,  0.008 ]),
	    " C6 ": np.array([   1.662, -1.342, -0.163 ]),
	    " O6 ": np.array([   2.758, -0.792, -0.331 ]),
	    " C5 ": np.array([   1.367, -2.729, -0.126 ]),
	    " N7 ": np.array([   2.200, -3.832, -0.264 ]),
	    " C8 ": np.array([   1.403, -4.858, -0.145 ]),
	    " H8 ": np.array([   1.938, -5.698, -0.234 ]),
	    " C2'": np.array([  -1.074, -6.552, -0.756 ]),
	    "H2''": np.array([  -0.189, -7.016, -0.790 ]),
	    " H2'": np.array([  -1.491, -6.273, -1.621 ]),
	    " C3'": np.array([  -2.064, -7.426,  0.015 ]),
	    " H3'": np.array([  -1.916, -8.374, -0.266 ]),
	    " O3'": np.array([  -3.370, -6.889, -0.156 ])}
    }
    }

    return bp

def init_pdb_charmm():
    pdb = {
    "A": {
            " P  ": np.array([   0.288, -9.220,  -1.848 ]),
	    " C4'": np.array([   3.212, -6.864,  -1.355 ]),
	    " H4'": np.array([   4.069, -6.713,  -1.848 ]),
	    " O4'": np.array([   2.387, -5.664,  -1.352 ]),
	    " C1'": np.array([   2.281, -5.198,  -0.016 ]),
	    " H1'": np.array([   3.030, -4.578,   0.218 ]),
	    " C2'": np.array([   2.304, -6.454,   0.850 ]),
	    " H2'": np.array([   1.451, -6.973,   0.796 ]),
	    "H2''": np.array([   2.514, -6.256,   1.807 ]),
	    " O1P": np.array([   0.421,-10.485,  -2.605 ]),
	    " O2P": np.array([  -0.692, -9.226,  -0.740 ]),
	    " O5'": np.array([   1.721, -8.770,  -1.295 ]),
	    " C5'": np.array([   2.585, -7.995,  -2.146 ]),
	    " H5'": np.array([   2.038, -7.632,  -2.900 ]),
	    "H5''": np.array([   3.293, -8.601,  -2.509 ]),
	    " N9 ": np.array([   1.021, -4.412,   0.101 ]),
	    " C5 ": np.array([  -0.409, -2.735,   0.190 ]),
	    " N7 ": np.array([  -1.161, -3.897,   0.354 ]),
	    " C8 ": np.array([  -0.268, -4.856,   0.294 ]),
	    " H8 ": np.array([  -0.503, -5.824,   0.381 ]),
	    " N1 ": np.array([   0.269, -0.516,  -0.025 ]),
	    " C2 ": np.array([   1.584, -0.914,  -0.172 ]),
	    " H2 ": np.array([   2.282, -0.210,  -0.302 ]),
	    " N3 ": np.array([   1.969, -2.189,  -0.149 ]),
	    " C4 ": np.array([   0.923, -3.037,   0.035 ]),
	    " C6 ": np.array([  -0.820, -1.378,   0.165 ]),
	    " N6 ": np.array([  -1.951, -0.913,   0.284 ]),
	    " H61": np.array([  -2.731, -1.524,   0.420 ]),
	    " H62": np.array([  -2.091,  0.076,   0.245 ]),
	    " C3'": np.array([   3.458, -7.146,   0.127 ]),
	    " H3'": np.array([   3.568, -8.115,   0.350 ]),
	    " O3'": np.array([   4.677, -6.520,   0.518 ])},
    "T": {  
	    " P  ": np.array([   0.249, -9.221,  -1.866 ]),
	    " C4'": np.array([   3.176, -6.866,  -1.373 ]),
	    " H4'": np.array([   4.034, -6.717,  -1.864 ]),
	    " O4'": np.array([   2.348, -5.666,  -1.370 ]),
	    " C1'": np.array([   2.243, -5.199,  -0.034 ]),
	    " H1'": np.array([   2.994, -4.581,   0.200 ]),
	    " C2'": np.array([   2.265, -6.453,   0.832 ]),
	    " H2'": np.array([   1.411, -6.971,   0.778 ]),
	    "H2''": np.array([   2.474, -6.254,   1.789 ]),
	    " O1P": np.array([   0.383,-10.486,  -2.623 ]),
	    " O2P": np.array([  -0.730, -9.227,  -0.758 ]),
	    " O5'": np.array([   1.683, -8.771,  -1.313 ]),
	    " C5'": np.array([   2.547, -7.996,  -2.164 ]),
	    " H5'": np.array([   2.000, -7.632,  -2.918 ]),
	    "H5''": np.array([   3.254, -8.602,  -2.528 ]),
	    " N1 ": np.array([   0.982, -4.412,   0.083 ]),
	    " C6 ": np.array([  -0.215, -5.030,   0.274 ]),
	    " H6 ": np.array([  -0.250, -6.027,   0.336 ]),
	    " C2 ": np.array([   1.072, -3.026,  -0.008 ]),
	    " O2 ": np.array([   2.183, -2.509,  -0.181 ]),
	    " N3 ": np.array([  -0.065, -2.290,   0.097 ]),
	    " H3 ": np.array([  -0.018, -1.293,   0.034 ]),
	    " C4 ": np.array([  -1.249, -2.888,   0.284 ]),
	    " O4 ": np.array([  -2.329, -2.129,   0.380 ]),
	    " C5 ": np.array([  -1.361, -4.314,   0.381 ]),
	    " C5M": np.array([  -2.245, -4.760,   0.521 ]),
	    " H51": np.array([  -2.964, -4.067,   0.574 ]),
	    " H52": np.array([  -2.432, -5.376,  -0.244 ]),
	    " H53": np.array([  -2.223, -5.283,   1.373 ]),
	    " C3'": np.array([   3.421, -7.147,   0.109 ]),
	    " H3'": np.array([   3.532, -8.115,   0.332 ]),
	    " O3'": np.array([   4.638, -6.523,   0.500 ])},
    "G": {
	    " P  ": np.array([   0.288, -9.220,  -1.848 ]),
	    " C4'": np.array([   3.212, -6.864,  -1.355 ]),
	    " H4'": np.array([   4.069, -6.713,  -1.848 ]),
	    " O4'": np.array([   2.387, -5.664,  -1.352 ]),
	    " C1'": np.array([   2.281, -5.198,  -0.016 ]),
	    " H1'": np.array([   3.030, -4.578,   0.218 ]),
	    " C2'": np.array([   2.304, -6.454,   0.850 ]),
	    " H2'": np.array([   1.451, -6.973,   0.796 ]),
	    "H2''": np.array([   2.514, -6.256,   1.807 ]),
	    " O1P": np.array([   0.421,-10.485,  -2.605 ]),
	    " O2P": np.array([  -0.692, -9.226,  -0.740 ]),
	    " O5'": np.array([   1.721, -8.770,  -1.295 ]),
	    " C5'": np.array([   2.585, -7.995,  -2.146 ]),
	    " H5'": np.array([   2.038, -7.632,  -2.900 ]),
	    "H5''": np.array([   3.293, -8.601,  -2.509 ]),
	    " N9 ": np.array([   1.021, -4.412,   0.101 ]),
	    " C4 ": np.array([   0.923, -3.037,   0.035 ]),
	    " N2 ": np.array([   2.478,  0.061,  -0.344 ]),
	    " H21": np.array([   2.180,  1.015,  -0.361 ]),
	    " H22": np.array([   3.446, -0.163,  -0.456 ]),
	    " N3 ": np.array([   1.969, -2.189,  -0.149 ]),
	    " C2 ": np.array([   1.584, -0.914,  -0.172 ]),
	    " N1 ": np.array([   0.269, -0.516,  -0.025 ]),
	    " H1 ": np.array([   0.079,  0.465,  -0.057 ]),
	    " C6 ": np.array([  -0.820, -1.378,   0.165 ]),
	    " O6 ": np.array([  -1.951, -0.913,   0.284 ]),
	    " C5 ": np.array([  -0.409, -2.735,   0.190 ]),
	    " N7 ": np.array([  -1.161, -3.897,   0.354 ]),
	    " C8 ": np.array([  -0.268, -4.856,   0.294 ]),
	    " H8 ": np.array([  -0.503, -5.824,   0.381 ]),
	    " C3'": np.array([   3.458, -7.146,   0.127 ]),
	    " H3'": np.array([   3.568, -8.115,   0.350 ]),
	    " O3'": np.array([   4.677, -6.520,   0.518 ])},
    "C": {
	    " P  ": np.array([   0.249, -9.221,  -1.866 ]),
	    " C4'": np.array([   3.176, -6.866,  -1.373 ]),
	    " H4'": np.array([   4.034, -6.717,  -1.864 ]),
	    " O4'": np.array([   2.348, -5.666,  -1.370 ]),
	    " C1'": np.array([   2.243, -5.199,  -0.034 ]),
	    " H1'": np.array([   2.994, -4.581,   0.200 ]),
	    " C2'": np.array([   2.265, -6.453,   0.832 ]),
	    " H2'": np.array([   1.411, -6.971,   0.778 ]),
	    "H2''": np.array([   2.474, -6.254,   1.789 ]),
	    " O1P": np.array([   0.383,-10.486,  -2.623 ]),
	    " O2P": np.array([  -0.730, -9.227,  -0.758 ]),
	    " O5'": np.array([   1.683, -8.771,  -1.313 ]),
	    " C5'": np.array([   2.547, -7.996,  -2.164 ]),
	    " H5'": np.array([   2.000, -7.632,  -2.918 ]),
	    "H5''": np.array([   3.254, -8.602,  -2.528 ]),
	    " N1 ": np.array([   0.982, -4.412,   0.083 ]),
	    " C6 ": np.array([  -0.215, -5.030,   0.274 ]),
	    " H6 ": np.array([  -0.250, -6.027,   0.336 ]),
	    " C5 ": np.array([  -1.361, -4.314,   0.381 ]),
	    " H5 ": np.array([  -2.245, -4.760,   0.521 ]),
	    " C2 ": np.array([   1.072, -3.026,  -0.008 ]),
	    " O2 ": np.array([   2.183, -2.509,  -0.181 ]),
	    " N3 ": np.array([  -0.065, -2.290,   0.097 ]),
	    " C4 ": np.array([  -1.249, -2.888,   0.284 ]),
	    " N4 ": np.array([  -2.329, -2.129,   0.380 ]),
	    " H41": np.array([  -3.225, -2.550,   0.520 ]),
	    " H42": np.array([  -2.249, -1.135,   0.312 ]),
	    " C3'": np.array([   3.421, -7.147,   0.109 ]),
	    " H3'": np.array([   3.532, -8.115,   0.332 ]),
	    " O3'": np.array([   4.638, -6.523,   0.500 ])},
    "A3": {
            " P  ": np.array([   0.288, -9.220,  -1.848 ]),
	    " C4'": np.array([   3.212, -6.864,  -1.355 ]),
	    " H4'": np.array([   4.069, -6.713,  -1.848 ]),
	    " O4'": np.array([   2.387, -5.664,  -1.352 ]),
	    " C1'": np.array([   2.281, -5.198,  -0.016 ]),
	    " H1'": np.array([   3.030, -4.578,   0.218 ]),
	    " C2'": np.array([   2.304, -6.454,   0.850 ]),
	    " H2'": np.array([   1.451, -6.973,   0.796 ]),
	    "H2''": np.array([   2.514, -6.256,   1.807 ]),
	    " C3'": np.array([   3.458, -7.146,   0.127 ]),
	    " H3'": np.array([   3.568, -8.115,   0.350 ]),
	    " O3'": np.array([   4.677, -6.520,   0.518 ]),
	    " O1P": np.array([   0.421,-10.485,  -2.605 ]),
	    " O2P": np.array([  -0.692, -9.226,  -0.740 ]),
	    " O5'": np.array([   1.721, -8.770,  -1.295 ]),
	    " C5'": np.array([   2.585, -7.995,  -2.146 ]),
	    " H5'": np.array([   2.038, -7.632,  -2.900 ]),
	    "H5''": np.array([   3.293, -8.601,  -2.509 ]),
	    " N9 ": np.array([   1.021, -4.412,   0.101 ]),
	    " C5 ": np.array([  -0.409, -2.735,   0.190 ]),
	    " N7 ": np.array([  -1.161, -3.897,   0.354 ]),
	    " C8 ": np.array([  -0.268, -4.856,   0.294 ]),
	    " H8 ": np.array([  -0.503, -5.824,   0.381 ]),
	    " N1 ": np.array([   0.269, -0.516,  -0.025 ]),
	    " C2 ": np.array([   1.584, -0.914,  -0.172 ]),
	    " H2 ": np.array([   2.282, -0.210,  -0.302 ]),
	    " N3 ": np.array([   1.969, -2.189,  -0.149 ]),
	    " C4 ": np.array([   0.923, -3.037,   0.035 ]),
	    " C6 ": np.array([  -0.820, -1.378,   0.165 ]),
	    " N6 ": np.array([  -1.951, -0.913,   0.284 ]),
	    " H61": np.array([  -2.731, -1.524,   0.420 ]),
	    " H62": np.array([  -2.091,  0.076,   0.245 ])},
    "T3": { 
	    " P  ": np.array([   0.249, -9.221,  -1.866 ]),
	    " C4'": np.array([   3.176, -6.866,  -1.373 ]),
	    " H4'": np.array([   4.034, -6.717,  -1.864 ]),
	    " O4'": np.array([   2.348, -5.666,  -1.370 ]),
	    " C1'": np.array([   2.243, -5.199,  -0.034 ]),
	    " H1'": np.array([   2.994, -4.581,   0.200 ]),
	    " C2'": np.array([   2.265, -6.453,   0.832 ]),
	    " H2'": np.array([   1.411, -6.971,   0.778 ]),
	    "H2''": np.array([   2.474, -6.254,   1.789 ]),
	    " C3'": np.array([   3.421, -7.147,   0.109 ]),
	    " H3'": np.array([   3.532, -8.115,   0.332 ]),
	    " O3'": np.array([   4.638, -6.523,   0.500 ]),
	    " O1P": np.array([   0.383,-10.486,  -2.623 ]),
	    " O2P": np.array([  -0.730, -9.227,  -0.758 ]),
	    " O5'": np.array([   1.683, -8.771,  -1.313 ]),
	    " C5'": np.array([   2.547, -7.996,  -2.164 ]),
	    " H5'": np.array([   2.000, -7.632,  -2.918 ]),
	    "H5''": np.array([   3.254, -8.602,  -2.528 ]),
	    " N1 ": np.array([   0.982, -4.412,   0.083 ]),
	    " C6 ": np.array([  -0.215, -5.030,   0.274 ]),
	    " H6 ": np.array([  -0.250, -6.027,   0.336 ]),
	    " C2 ": np.array([   1.072, -3.026,  -0.008 ]),
	    " O2 ": np.array([   2.183, -2.509,  -0.181 ]),
	    " N3 ": np.array([  -0.065, -2.290,   0.097 ]),
	    " H3 ": np.array([  -0.018, -1.293,   0.034 ]),
	    " C4 ": np.array([  -1.249, -2.888,   0.284 ]),
	    " O4 ": np.array([  -2.329, -2.129,   0.380 ]),
	    " C5 ": np.array([  -1.361, -4.314,   0.381 ]),
	    " C5M": np.array([  -2.245, -4.760,   0.521 ]),
	    " H51": np.array([  -2.964, -4.067,   0.574 ]),
	    " H52": np.array([  -2.432, -5.376,  -0.244 ]),
	    " H53": np.array([  -2.223, -5.283,   1.373 ])},
    "G3": {
	    " P  ": np.array([   0.288, -9.220,  -1.848 ]),
	    " C4'": np.array([   3.212, -6.864,  -1.355 ]),
	    " H4'": np.array([   4.069, -6.713,  -1.848 ]),
	    " O4'": np.array([   2.387, -5.664,  -1.352 ]),
	    " C1'": np.array([   2.281, -5.198,  -0.016 ]),
	    " H1'": np.array([   3.030, -4.578,   0.218 ]),
	    " C2'": np.array([   2.304, -6.454,   0.850 ]),
	    " H2'": np.array([   1.451, -6.973,   0.796 ]),
	    "H2''": np.array([   2.514, -6.256,   1.807 ]),
	    " C3'": np.array([   3.458, -7.146,   0.127 ]),
	    " H3'": np.array([   3.568, -8.115,   0.350 ]),
	    " O3'": np.array([   4.677, -6.520,   0.518 ]),
	    " O1P": np.array([   0.421,-10.485,  -2.605 ]),
	    " O2P": np.array([  -0.692, -9.226,  -0.740 ]),
	    " O5'": np.array([   1.721, -8.770,  -1.295 ]),
	    " C5'": np.array([   2.585, -7.995,  -2.146 ]),
	    " H5'": np.array([   2.038, -7.632,  -2.900 ]),
	    "H5''": np.array([   3.293, -8.601,  -2.509 ]),
	    " N9 ": np.array([   1.021, -4.412,   0.101 ]),
	    " C4 ": np.array([   0.923, -3.037,   0.035 ]),
	    " N2 ": np.array([   2.478,  0.061,  -0.344 ]),
	    " H21": np.array([   2.180,  1.015,  -0.361 ]),
	    " H22": np.array([   3.446, -0.163,  -0.456 ]),
	    " N3 ": np.array([   1.969, -2.189,  -0.149 ]),
	    " C2 ": np.array([   1.584, -0.914,  -0.172 ]),
	    " N1 ": np.array([   0.269, -0.516,  -0.025 ]),
	    " H1 ": np.array([   0.079,  0.465,  -0.057 ]),
	    " C6 ": np.array([  -0.820, -1.378,   0.165 ]),
	    " O6 ": np.array([  -1.951, -0.913,   0.284 ]),
	    " C5 ": np.array([  -0.409, -2.735,   0.190 ]),
	    " N7 ": np.array([  -1.161, -3.897,   0.354 ]),
	    " C8 ": np.array([  -0.268, -4.856,   0.294 ]),
	    " H8 ": np.array([  -0.503, -5.824,   0.381 ])},
    "C3": {
	    " P  ": np.array([   0.249, -9.221,  -1.866 ]),
	    " C4'": np.array([   3.176, -6.866,  -1.373 ]),
	    " H4'": np.array([   4.034, -6.717,  -1.864 ]),
	    " O4'": np.array([   2.348, -5.666,  -1.370 ]),
	    " C1'": np.array([   2.243, -5.199,  -0.034 ]),
	    " H1'": np.array([   2.994, -4.581,   0.200 ]),
	    " C2'": np.array([   2.265, -6.453,   0.832 ]),
	    " H2'": np.array([   1.411, -6.971,   0.778 ]),
	    "H2''": np.array([   2.474, -6.254,   1.789 ]),
	    " C3'": np.array([   3.421, -7.147,   0.109 ]),
	    " H3'": np.array([   3.532, -8.115,   0.332 ]),
	    " O3'": np.array([   4.638, -6.523,   0.500 ]),
	    " O1P": np.array([   0.383,-10.486,  -2.623 ]),
	    " O2P": np.array([  -0.730, -9.227,  -0.758 ]),
	    " O5'": np.array([   1.683, -8.771,  -1.313 ]),
	    " C5'": np.array([   2.547, -7.996,  -2.164 ]),
	    " H5'": np.array([   2.000, -7.632,  -2.918 ]),
	    "H5''": np.array([   3.254, -8.602,  -2.528 ]),
	    " N1 ": np.array([   0.982, -4.412,   0.083 ]),
	    " C6 ": np.array([  -0.215, -5.030,   0.274 ]),
	    " H6 ": np.array([  -0.250, -6.027,   0.336 ]),
	    " C5 ": np.array([  -1.361, -4.314,   0.381 ]),
	    " H5 ": np.array([  -2.245, -4.760,   0.521 ]),
	    " C2 ": np.array([   1.072, -3.026,  -0.008 ]),
	    " O2 ": np.array([   2.183, -2.509,  -0.181 ]),
	    " N3 ": np.array([  -0.065, -2.290,   0.097 ]),
	    " C4 ": np.array([  -1.249, -2.888,   0.284 ]),
	    " N4 ": np.array([  -2.329, -2.129,   0.380 ]),
	    " H41": np.array([  -3.225, -2.550,   0.520 ]),
	    " H42": np.array([  -2.249, -1.135,   0.312 ])}
    }

    return pdb


#####################################################
## Main function ####################################
#####################################################

pdb = init_pdb_charmm()
bp = init_bp()
anchor_atoms_bp = np.array([" H1'"])
anchor_atoms_nt = np.array([" H3'"])


#######################
## Parsing the mmp file
#######################

token_1 = 'group'
token_2 = 'atom'
chunks = []
current_chunk_1 = []
current_chunk_2 = ''

with open(sys.argv[1], 'r') as inFile:

    for line in inFile:

        if line.startswith(token_1):
            current_chunk_1.append(current_chunk_2) 
            current_chunk_2 = '' 
            chunks.append(current_chunk_1)
            current_chunk_1 = []

        if line.startswith(token_2) and current_chunk_2: 
        # if line starts with token and the current chunk is not empty
            current_chunk_1.append(current_chunk_2) #  add not empty chunk to chunks
            current_chunk_2 = '' #  make current chunk blank
        # just append a line to the current chunk on each iteration
        current_chunk_2 += line
    
    current_chunk_1.append(current_chunk_2) #  add not empty chunk to chunks
    chunks.append(current_chunk_1)  #  append the last chunk outside the loop    



data = {}
# find group name
pattern1 = re.compile('group \((.*)\) (DnaStrand|DnaSegment)')
# find atom index and coordinate
pattern2 = re.compile('atom ([0-9]*) \([0-9][0-9]+\) \((-*[0-9]*), (-*[0-9]*), (-*[0-9]*)\) def')
# find bonds
pattern3 = re.compile('bond1 ([0-9]*) *([0-9]*)')
# find bond direction
pattern4 = re.compile('bond_direction ([0-9]*) ([0-9]*)')
# find sequence
pattern5 = re.compile('info atom dnaBaseName = ([ATCG])')

for i in chunks:

    # initialize each group
    if (pattern1.search(i[0])):
        m1 = pattern1.search(i[0])
        group = m1.group(1)
        data[group] = {}

    for j in i:
        # If there is a valid atom, get atom index and coordinate
        if (pattern2.search(j)):
            m2 = pattern2.search(j)
            aid = int(m2.group(1))
            data[group][aid] = {}
            data[group][aid]['xyz'] = np.array([float(m2.group(2)) / 1000.0, float(m2.group(3)) / 1000.0, float(m2.group(4)) / 1000.0])
            # Then, get bonds
            if (pattern3.search(j)):
                m3 = pattern3.search(j)
#                print(m3.group(0), m3.groups())
                if not m3.group(2) == '':
                    data[group][aid]['bond1'] = np.array([int(m3.group(1)), int(m3.group(2))])
                else:
                    data[group][aid]['bond1'] = np.array([int(m3.group(1))])

            # Then, get bond direction
            if (pattern4.search(j)):
                m4 = pattern4.search(j)
                data[group][aid]['bond_direction'] = np.array([int(m4.group(1)), int(m4.group(2))])

            # Then, get sequence
            if (pattern5.search(j)):
                m5 = pattern5.search(j)
                data[group][aid]['seq'] = m5.group(1)

#print(data)


###################
## Process the data
###################
###################################################
# data format:
#
# (from mmp file)
# ----- 1 layer -----
# Group: DnaSegment[0-9]* -> DNA segment, central axis
#        DnaStrand[0-9]*  -> actual DNA strand with sequence
# ----- 2 layer -----
# unique indices of nucleodites
# ----- 3 layer -----
# bond1: bonding information
# bond_direction: 5' to 3' direction
# xyz: coordinate
# seq: sequence
#
# (custom defined)
# direc on DnaStrand: n to n+1 == 5' to 3' -> direc = 1
# ntBond on DnaSegment: connected nucleotides indices on each elements of DnaSegment; first index has direc = 1
###################################################

## get DnaSegment and DnaStrand group names
key_DnaSegment = np.array([])
key_DnaStrand  = np.array([])
for key in data.keys():
    if re.search('DnaSegment', key):
        key_DnaSegment = np.append(key_DnaSegment, key)
    if re.search('DnaStrand', key):
        key_DnaStrand = np.append(key_DnaStrand, key)


## determine direction: 
## n to n+1 == 5' to 3' -> direc = 1
for key in key_DnaStrand:
    for k, v in data[key].items():
       if 'bond_direction' in data[key][k]:
           bond_direction = data[key][k]['bond_direction']
           if bond_direction[0] == k:
               data[key][k]['direc'] = -1 
           else:
               data[key][k]['direc'] = 1
       else:
           bond_direction = data[key][k + 1]['bond_direction']
           if bond_direction[0] == k:
               data[key][k]['direc'] = 1 
           else:
               data[key][k]['direc'] = -1

## get connected nucleotides on each elements of DnaSegment
for key in key_DnaSegment:
    for k, v in data[key].items():
        data[key][k]['ntBond'] = np.array([], dtype = np.uint16)

        for Strand in key_DnaStrand:
            for index in data[Strand]:
                if k in data[Strand][index]['bond1']: 
                    data[key][k]['ntBond'] = np.append(data[key][k]['ntBond'], index)

        # make sure the ['direc'] in ['ntBond'][0] is 1
        for Strand in key_DnaStrand:
            if data[key][k]['ntBond'][0] in data[Strand]:
                if data[Strand][ data[key][k]['ntBond'][0] ]['direc'] == -1:
                    data[key][k]['ntBond'] = data[key][k]['ntBond'][::-1]


## fill in the coordinates
## start with double-stranded region
for Seg in key_DnaSegment:
    for k, v in data[Seg].items():

        # if double-strand
        if len(data[Seg][k]['ntBond']) == 2:

            # find coorsponding seq
            ind_1 = data[Seg][k]['ntBond'][0]
            ind_2 = data[Seg][k]['ntBond'][1]
            for Strand in key_DnaStrand:
                if ind_1 in data[Strand]:
                    seq_1 = data[Strand][ind_1]['seq']
                    strand_1 = Strand
                if ind_2 in data[Strand]:
                    seq_2 = data[Strand][ind_2]['seq']
                    strand_2 = Strand

            # find coorsponding base-pair name
            for bp_key in bp.keys():
                if re.match(seq_1, bp_key):
                    seq_12 = bp_key

            # translation
            bp_m = translation_bp(Seg, k, bp[seq_12])            

            # rotation
            bp_m_r = rotation_bp(Seg, k, strand_1, ind_1, seq_1, strand_2, ind_2, seq_2, bp_m)            

            # assign coordinates to each residue in DnaStrand
            data[strand_1][ind_1]['xyz_aa'] = bp_m_r[seq_1]
            data[strand_2][ind_2]['xyz_aa'] = bp_m_r[seq_2]
            

## then, single-stranded region
for Strand in key_DnaStrand:
    for k, v in data[Strand].items():

        if not 'xyz_aa' in data[Strand][k]:

            # Translation
            pdb_m = translation_nt(Strand, k) 
    
            # Rotation
            data[Strand][k]['xyz_aa'] = rotation_nt(Strand, k, pdb_m) 



#print(data)


#############
## output pdb
#############

index = 0
with open(sys.argv[2], 'w') as out:

    for key in data.keys():

        resid = 0
        if key in key_DnaStrand:
    
            m = re.search('DnaStrand([0-9]*)', key)
            segn = int(m.group(1))

            k_array = np.sort(np.array(list(data[key].keys())))

            direc = data[key][k_array[0]]['direc']
            k_array = k_array[::direc]

            for k in k_array:
                v = data[key][k]

                resid += 1
                print(k)
                seq = data[key][k]['seq']

                #resid = k_array[::direc][k_array == k] 
                if 'xyz_aa' in v:
                    for aname, coor in v['xyz_aa'].items():
                        index += 1
                        out.write("%-6s%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" % ("ATOM", index, aname, nt1to3(seq), 'A', resid, coor[0], coor[1], coor[2], 1,0, 'D%s' % segn))

## Test
#p1   = np.array([3,3,-5])
#ori  = np.array([3,4,-5])
#dest = np.array([-4,-3,5])
#print(ori-p1)
##print(dest)
#anchor = np.array([0,0,0])
#tmp = {'Y': np.array([3,4,-5])}
#r_m = rotateMatrix_2(dest, tmp['Y'])
#
#tmp = rotateby_nt(tmp, r_m, anchor)
#print(tmp, np.linalg.norm(dest - tmp['Y']))
#p1 = rotateby_nt({'p1':p1}, r_m, anchor)
#print(p1, np.linalg.norm(tmp['Y'] - p1['p1']))

#center = measure_center_bp(bp['AT'])
#v1 = bp['AT']['A'][" H1'"] - center
#v2 = bp['AT']['T'][" H1'"] - center
#
#print(vecAngle(v1, v2) * 180 / pi)

