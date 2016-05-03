#!/usr/bin/env python3

# author: Chen-Yu Li cli56@illinois.edu

import sys, optparse, json, os, glob, re, random, numpy
from subprocess import Popen, PIPE
from math import *
import numpy as np
import pandas as pd


###################################################
# data format:
#
# (from mmp file)
# ----- 1 layer -----
# Group: DnaSegment[0-9]* -> DNA segment, central axis
#        DnaStrand[0-9]*  -> actual DNA strand with sequence
# There may be more than one group layer!!
# ----- 2 layer -----
# unique indices of nucleodites
# ----- 3 layer -----
# bond1: bonding information
# bond_direction: 5' to 3' direction
# xyz: coordinate
# seq: sequence
#
#
# (flatten to pandas DataFrame 'data_df')
# index: point index
#
# columns: 
# In addition to infomation from mmp file:
# (DNA strand)
# segname: str, T*
# direc: int, n to n+1 == 5' to 3' -> direc = 1
# n_nt, p_nt, b_nt: int, next, previous, base-paired nucleotide index
# E_seg: str, connected DNA segment index
# xyz_aa: dict, all-atom coordinate dictionary
# residue: int, residue id in the order of output pdb, 1-based
# CW: str, "%s%s%s" % (sname, rname, aname) 
# (DNA segment)
# segname: E*
# ntBond: np.array, connected nucleotides indices; first index has direc = 1
# n_seg, p_seg: int, next, previous segment index 
# c_seg: np.array, cross-over segment index 
# E_ind: int, index order in its DNA segment
#
#
# Else
# seg_map: key=segname; value=group path
# seg_head: key=segname; value=5' nucleotide index
# 
# 
# E_df
# index: sorted DNA segment E*
# 
# columns:
# indice: np.array, sorted point index of each DNA segment E*
# vec: np.array, unit vector of indice[1] -> indice[0]
# para_group: int, group number, distinguished by vec
# nn: nested np.array, [nearest neighbor DNA segment E*, index of itself, index of nearest neighbor]
#
###################################################


###############################################
## subroutines 
###############################################

def wc(seq):
    if re.search("[atgcATGC]",seq):
        intab  = "ATGC"
        outtab = "TACG"
        trantab = str.maketrans(intab, outtab)
        wc = seq.upper().translate(trantab)

    else:
        raise Exception("wc(): wrong nucleotide '%s'." % seq)
    return wc


def find_central(rid): 

    global data_df

    ## find two continuous nt with 'E_seg' in upstream
    if data_df.ix[rid]['p_nt'] != 'NA':
        k = data_df.ix[rid]['p_nt']
        interval = 1
        while data_df.ix[k]['p_nt'] != 'NA':
            if data_df.ix[k]['E_seg'] != 'NA' and data_df.ix[ data_df.ix[k]['p_nt'] ]['E_seg'] != 'NA':
                c0 = data_df.ix[ data_df.ix[k]['p_nt'] ]['E_seg']
                c1 = data_df.ix[k]['E_seg']
                coor0 = data_df.ix[c0]['xyz']
                coor1 = data_df.ix[c1]['xyz']
                vec = coor1 - coor0
                head = vec * interval + coor1
                if c0 > c1:
                    axis_vec = unitVec(coor1 - coor0)            
                else:
                    axis_vec = unitVec(coor0 - coor1)

                break
            else:
                k = data_df.ix[k]['p_nt']
                interval += 1

    ## if not found, try downstream 
    if 'head' not in locals():
        
        k = data_df.ix[rid]['n_nt']
        interval = 1
        while data_df.ix[k]['n_nt'] != 'NA':
            if data_df.ix[k]['E_seg'] != 'NA' and data_df.ix[ data_df.ix[k]['n_nt'] ]['E_seg'] != 'NA':
                c0 = data_df.ix[ data_df.ix[k]['n_nt'] ]['E_seg']
                c1 = data_df.ix[k]['E_seg']
                coor0 = data_df.ix[c0]['xyz']
                coor1 = data_df.ix[c1]['xyz']
                vec = coor1 - coor0
                head = vec * interval + coor1
                if c0 > c1:
                    axis_vec = unitVec(coor1 - coor0)            
                else:
                    axis_vec = unitVec(coor0 - coor1)

                break
            else:
                k = data_df.ix[k]['n_nt']
                interval += 1

    return (head, axis_vec)

def close_atom(array, index):

    global data_df

    close = 0
    minimum = 99999
    for i in array:
        dist = np.linalg.norm(data_df.ix[i]['xyz'] - data_df.ix[index]['xyz']) 
        if dist < minimum:
            minimum = dist
            close = i

    return close

def trace_seg(array1):

    global data_df

    ## find one of the two ends
    for i in array1:
        bond = 0
        if str(data_df.ix[i]['bond1']) != 'NA':
            for j in data_df.ix[i]['bond1']:
                if j in array1:
                    bond += 1
        
        for j in array1:
            if i != j and str(data_df.ix[j]['bond1']) != 'NA' and i in data_df.ix[j]['bond1']:
                bond += 1

        # if the one of the terminal is found
        if bond == 1:
            start = i
            break


    ## traverse from the found end
    array2 = np.array([start], dtype=np.int32)
    
    for count in range(len(array1)-1):
        record=0
        now = array2[-1]
        # start from 'bond1' of itself
        if str(data_df.ix[now]['bond1']) != 'NA':
            for i in data_df.ix[now]['bond1']:
                if i in array1 and i not in array2:
                    array2 = np.append(array2, i)
                    record=1

        # if nothing found, find others which bind to itself        
        if record != 1:
            for i in array1:
                if i not in array2:
                    if str(data_df.ix[i]['bond1']) != 'NA' and now in data_df.ix[i]['bond1']:
                        array2 = np.append(array2, i)

    return array2


def arrayRearrange(array1, element, destRow):

    t,interval = array1.shape
    array2 = array1.copy()

    if element in array2:
        row, col = np.where(array2 == element)
        if row != destRow:
            tmp = array2[row]
            array2 = np.delete(array2, row, axis=0)
            array2 = np.insert(array2, destRow, tmp, axis=0)
        
        #if col != 0:
        #    array2[0] = array2[0][::-1]

    return array2


def parseString(s):

    # s should be a string start with 'atom'
    atom_dict = {}

    # find group name
    #pattern1 = re.compile('group \((.*)\)')
    # find atom index and coordinate
    pattern2 = re.compile('atom ([0-9]*) \([0-9][0-9]+\) \((-*[0-9]*\.*[0-9]*), (-*[0-9]*\.*[0-9]*), (-*[0-9]*\.*[0-9]*)\) def')
    # find bonds
    pattern3 = re.compile('bond1 ([0-9]*) *([0-9]*) *([0-9]*) *([0-9]*)')
    # find bond direction
    pattern4 = re.compile('bond_direction ([0-9]*) ([0-9]*)')
    # find sequence
    pattern5 = re.compile('info atom dnaBaseName = ([ATCG])')

    # If there is a valid atom, get atom index and coordinate
    if (pattern2.search(s)):
        m2 = pattern2.search(s)
        #aid = int(m2.group(1))
        atom_dict['xyz'] = np.array([float(m2.group(2)) / 1000.0, float(m2.group(3)) / 1000.0, float(m2.group(4)) / 1000.0])
        # Then, get bonds
        if (pattern3.search(s)):
            m3 = pattern3.search(s)
            bonds = [x for x in m3.groups() if x != '']
            #print(m3.groups() ,bonds)
            atom_dict['bond1'] = np.array(bonds, dtype=np.uint16)

        # Then, get bond direction
        if (pattern4.search(s)):
            m4 = pattern4.findall(s)
            atom_dict['bond_direction'] = np.array([item for sublist in m4 for item in sublist], dtype = np.int32)

        # Then, get sequence
        if (pattern5.search(s)):
            m5 = pattern5.search(s)
            atom_dict['seq'] = m5.group(1)

    return atom_dict


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


def rotation_bp(k, ind_1, seq_1, ind_2, seq_2, bp_m):

    global data_df, anchor_atoms_bp, seg_cutOff

    # align base-pair
    center = measure_center_bp(bp_m)
    head = data_df.ix[ind_1]['xyz']
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
    #if 'n_seg' in data[k]:
    #    head = data[k]['xyz']
    #    tail = data[ data[k]['n_seg'] ]['xyz']
    #elif 'p_seg' in data[k]:
    #    head = data[ data[k]['p_seg'] ]['xyz']
    #    tail = data[k]['xyz']

    if k + 1 in data_df.index.values and np.linalg.norm(data_df.ix[k]['xyz'] - data_df.ix[k+1]['xyz']) < seg_cutOff:
        head = data_df.ix[k]['xyz']
        tail = data_df.ix[k + 1]['xyz']
    elif k - 1 in data_df.index.values:   
        head = data_df.ix[k - 1]['xyz']
        tail = data_df.ix[k]['xyz']
    elif data_df.ix[k]['p_seg'] != 'NA':
        head = data_df.ix[ data_df.ix[k]['p_seg'] ]['xyz']
        tail = data_df.ix[k]['xyz']
    elif data_df.ix[k]['n_seg'] != 'NA':
        head = data_df.ix[k]['xyz']
        tail = data_df.ix[ data_df.ix[k]['n_seg'] ]['xyz']

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

def rotation_nt(rid, pdb):

    global data_df, anchor_atoms_nt 
    
    direc = data_df.ix[rid]['direc']

    seq = data_df.ix[rid]['seq']

    central = data_df.ix[rid]['E_seg']
    if central == 'NA':
        head, axis_vec = find_central(rid)
    else:
        head = data_df.ix[central]['xyz']
    
    dest_vec = unitVec(head - data_df.ix[rid]['xyz'])

    
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
    if 'axis_vec' not in locals():
        if central + 1 in data_df.index.values:
            head = data_df.ix[central]['xyz']
            tail = data_df.ix[central + 1]['xyz']
            axis_vec = unitVec(head - tail)
        elif central - 1 in data_df.index.values:   
            head = data_df.ix[central - 1]['xyz']
            tail = data_df.ix[central]['xyz']
            axis_vec = unitVec(head - tail)

#        else:
#            raise Exception("rotation(group, rid, moved_pdb): Cannot define axis direction:\ngroup = %s\nrid = %s\ncentral = %s" % (group, rid, central))

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


def translation_bp(rid, bp_dict):

    global data_df, bp

    center = measure_center_bp(bp_dict)

    move_vec = data_df.ix[rid]['xyz'] - center 
    bp_m = moveby_bp(bp_dict, move_vec)

    return bp_m


def translation_nt(rid):

    global data_df, anchor_atoms_nt, pdb

    anchor_atoms_nt_pdb = np.array([])
    seq = data_df.ix[rid]['seq']
    for i in anchor_atoms_nt:
        anchor_atoms_nt_pdb = np.append(anchor_atoms_nt_pdb, pdb[seq][i], axis=0)
    anchor_center = measure_center(anchor_atoms_nt_pdb)
    move_vec = data_df.ix[rid]['xyz'] - anchor_center 
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


def traverse(o, tree_types=(list, tuple)):
    if isinstance(o, tree_types):
        for value in o:
            for subvalue in traverse(value):
                yield subvalue
    else:
        yield o


def depth_list(l):
    if isinstance(l, list):
        return 1 + max(depth_list(item) for item in l)
    else:
        return 0


def depth_dict(d, depth=0):

    if not isinstance(d, dict) or not d:
        return depth
    return max(depth_dict(v, depth+1) for k, v in d.items())


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
seg_cutOff = 5


#######################
## Parsing the mmp file
#######################

seg_map = {}
Segment_n = 0
Strand_n  = 0
data = {}
current_group = np.array([])
current_line = ''
pattern1 = re.compile('^group (.*)')
pattern2 = re.compile('DnaSegment')
pattern3 = re.compile('DnaStrand')
pattern4 = re.compile('egroup \((.*)\)')
pattern5 = re.compile('atom ([0-9]*) \([0-9][0-9]+\) \((-*[0-9]*\.*[0-9]*), (-*[0-9]*\.*[0-9]*), (-*[0-9]*\.*[0-9]*)\) def')
# find bonds
pattern6 = re.compile('bond1 ([0-9]*) *([0-9]*)')
# find bond direction
pattern7 = re.compile('bond_direction ([0-9]*) ([0-9]*)')
# find sequence
pattern8 = re.compile('info atom dnaBaseName = ([ATCG])')
with open(sys.argv[1], 'r') as inFile:
    for num, line in enumerate(inFile, 1):

        # if it starts with 'atom'
        if re.search('^atom', line):

            # if current line start with 'atom', parse it
            if re.search('^atom', current_line) and pattern5.search(current_line):
                m5 = pattern5.search(current_line)
                aid = int(m5.group(1))       
                #print(num, repr(current_line))
                data[aid] = parseString(current_line)
                current_line = ''

                last_group = current_group[-1]                
                if pattern2.search(last_group):
                    segname = 'E' + str(Segment_n)
                    data[aid]['segname'] = segname
                elif pattern3.search(last_group):
                    segname = 'T' + str(Strand_n)
                    data[aid]['segname'] = segname

            # else, remove current line
            else:
                current_line = ''

        # if it reaches 'group'
        if (pattern1.search(line)):
            m1 = pattern1.search(line)
            if re.search('\(Clipboard\)', m1.group(1)):
            #if m1.group(1) == 'Clipboard':
                break

            current_group = np.append(current_group, m1.group(1))
        #print(num, current_group)

        # append the line
        current_line += line
        #print(num, repr(current_line))

        # if it reaches 'egroup'
        if (pattern4.search(line)):
            
            # parse current line
            if re.search('^atom', current_line) and pattern5.search(current_line):
                m5 = pattern5.search(current_line)
                aid = int(m5.group(1))       
                data[aid] = parseString(current_line)
                current_line = ''

                last_group = current_group[-1]                
                if pattern2.search(last_group):
                    segname = 'E' + str(Segment_n)
                    data[aid]['segname'] = segname
                elif pattern3.search(last_group):
                    segname = 'T' + str(Strand_n)
                    data[aid]['segname'] = segname
            else:
                current_line = ''

            # add current group to seg_map, increment and leave the last group
            m4 = pattern4.search(line)
            if re.search(m4.group(1), current_group[-1]): 

                if 'segname' in locals():
                    group_path = ''
                    for i in current_group:
                        group_path = group_path + '/' + i
                    #print(num, segname, group_path)
                    seg_map[segname] = group_path
                    del segname

                last_group = current_group[-1]                
                if pattern2.search(last_group):
                    Segment_n += 1
                elif pattern3.search(last_group):
                    Strand_n += 1

                current_group = np.delete(current_group, len(current_group)-1)


print('Finished reading %s' % sys.argv[1])
#for k in np.sort(list(seg_map.keys())):
#    print('%s: %s' % (k, seg_map[k]))

#print(data[12446])

## convert to pandas dataframe
data_df = pd.DataFrame(data)
data_df = data_df.transpose()
data_df = data_df.fillna('NA')

#print(data_df.columns.values)
#print(data_df.index.values)
#print(data_df)
#sys.exit(0)


###################
## Process the data
###################

## determine direction: 
## n to n+1 == 5' to 3' -> direc = 1
## get p_nt, n_nt
data_df['direc'] = 'NA'
data_df['p_nt']  = 'NA'
data_df['n_nt']  = 'NA'
for i in data_df[ data_df['segname'].str.contains('T[0-9]*') ].index.values:

    # if it is not the first nucleotide
    if str(data_df.ix[i]['bond_direction']) != 'NA':

        #print(i, data_df.ix[i]['bond_direction'])
        bond_direction_all = data_df.ix[i]['bond_direction']

        # if there are more than one pair, use the continuous pair
        if len(bond_direction_all) > 2:
            bond_direction_all = np.reshape(bond_direction_all, (int(len(bond_direction_all) / 2), 2))
            for sub in bond_direction_all:
                if np.absolute(int(sub[0] - sub[1])) == 1:
                    bond_direction = sub
        else:
                bond_direction = bond_direction_all
        
        if bond_direction[0] == i:
            data_df.ix[i, 'direc'] = -1 
        else:
            data_df.ix[i, 'direc'] = 1

    # if it is the first nucleotide
    else:
        bond_direction_all = data_df.ix[i + 1]['bond_direction']

        if len(bond_direction_all) > 2:
            bond_direction_all = np.reshape(bond_direction_all, (int(len(bond_direction_all) / 2), 2))

            for sub in bond_direction_all:
                if np.absolute(sub[0] - sub[1]) == 1:
                    bond_direction = sub
        else:
            bond_direction = bond_direction_all

        if bond_direction[0] == i:
            data_df.ix[i, 'direc'] = 1 
        else:
            data_df.ix[i, 'direc'] = -1

    # get p_nt, n_nt
    if str(data_df.ix[i]['bond_direction']) != 'NA':
        bond_direction = data_df.ix[i]['bond_direction']
        indices, = np.where(bond_direction == i)
        for k in indices:
            if k % 2 == 0:
                data_df.ix[i, 'n_nt'] = bond_direction[k + 1]
                data_df.ix[ bond_direction[k + 1], 'p_nt'] = i
            else:
                data_df.ix[i, 'p_nt'] = bond_direction[k - 1]
                data_df.ix[ bond_direction[k - 1], 'n_nt'] = i


## get connected nucleotides on each elements of DnaSegment
## get b_nt, E_seg
data_df['ntBond'] = 'NA'
data_df['b_nt']  = 'NA'
data_df['E_seg'] = 'NA'
# parse 'bond1' of the DNA strands
for i in data_df[ data_df['segname'].str.contains('T[0-9]*') ].index.values:
    if str(data_df.ix[i]['bond1']) != 'NA':
        for j in data_df.ix[i]['bond1']:
            if re.search('E[0-9]*', data_df.ix[j]['segname']):
                if data_df.ix[j]['ntBond'] == 'NA':
                    data_df.ix[j, 'ntBond'] = np.array([i], dtype = np.uint16) 

                else:
                    data_df.set_value(j, 'ntBond', np.append(data_df.ix[j]['ntBond'], i))

                data_df.ix[i, 'E_seg'] = j


# parse 'bond1' of the DNA segments
for i in data_df[ data_df['segname'].str.contains('E[0-9]*') ].index.values:

    if str(data_df.ix[i]['bond1']) != 'NA':
        for j in data_df.ix[i]['bond1']:
            if re.search('T[0-9]*', data_df.ix[j]['segname']):
                if str(data_df.ix[i]['ntBond']) == 'NA':
                    data_df.ix[i, 'ntBond'] = np.array([j], dtype = np.uint16)
                else:
                    data_df.set_value(i, 'ntBond', np.append(data_df.ix[i]['ntBond'], j))
                data_df.ix[j, 'E_seg'] = i


    # make sure the ['direc'] in ['ntBond'][0] is 1
    if data_df.ix[i]['ntBond'] != 'NA' and data_df.ix[ data_df.ix[i]['ntBond'][0] ]['direc'] == -1:
        data_df.set_value(i, 'ntBond', data_df.ix[i]['ntBond'][::-1])


    # get b_nt
    if len(data_df.ix[i]['ntBond']) == 2 and isinstance(data_df.ix[i]['ntBond'], np.ndarray):
        data_df.ix[ data_df.ix[i]['ntBond'][0], 'b_nt'] = data_df.ix[i]['ntBond'][1]
        data_df.ix[ data_df.ix[i]['ntBond'][1], 'b_nt'] = data_df.ix[i]['ntBond'][0]


## get c_seg
data_df['c_seg'] = 'NA'
for i in data_df[ data_df['segname'].str.contains('E[0-9]*') ].index.values:
    if str(data_df.ix[i]['ntBond']) != 'NA':
        for j in data_df.ix[i]['ntBond']:


            if data_df.ix[j]['p_nt'] != 'NA' and \
               data_df.ix[ data_df.ix[j]['p_nt'] ]['E_seg'] != 'NA' and\
               data_df.ix[ data_df.ix[j]['E_seg'] ]['segname'] != data_df.ix[ data_df.ix[ data_df.ix[j]['p_nt'] ]['E_seg'] ]['segname']:   

                #print(i,j, data_df.ix[j]['p_nt'], data_df.ix[ data_df.ix[j]['p_nt'] ]['E_seg'])
                if str(data_df.ix[i]['c_seg']) == 'NA':
                    data_df.ix[i, 'c_seg'] = np.array([ data_df.ix[ data_df.ix[j]['p_nt'] ]['E_seg'] ])
                else:
                    data_df.set_value(i, 'c_seg', np.append(data_df.ix[i]['c_seg'], data_df.ix[ data_df.ix[j]['p_nt'] ]['E_seg']))


            if data_df.ix[j]['n_nt'] != 'NA' and \
               data_df.ix[ data_df.ix[j]['n_nt'] ]['E_seg'] != 'NA' and\
               data_df.ix[ data_df.ix[j]['E_seg'] ]['segname'] != data_df.ix[ data_df.ix[ data_df.ix[j]['n_nt'] ]['E_seg'] ]['segname']:    

                #print(i,j, data_df.ix[j]['n_nt'], data_df.ix[ data_df.ix[j]['n_nt'] ]['E_seg'])
                if str(data_df.ix[i]['c_seg']) == 'NA':
                    data_df.ix[i, 'c_seg'] = np.array([ data_df.ix[ data_df.ix[j]['n_nt'] ]['E_seg'] ])
                else:
                    data_df.set_value(i, 'c_seg', np.append(data_df.ix[i]['c_seg'], data_df.ix[ data_df.ix[j]['n_nt'] ]['E_seg']))


count = 0
for i in data_df[ data_df['c_seg'] != 'NA' ].index.values:
    count += len(data_df.ix[i]['c_seg'])

print('Total number of single crossover: %d' % int(count/2))


#for i in data_df[ data_df['c_seg'] != 'NA' ].index.values:
#    if i not in data_df.ix[ data_df.ix[i]['c_seg'] ]['c_seg'].values[0]:
#        print(i, data_df.ix[i]['segname'], data_df.ix[i]['c_seg'], data_df.ix[ data_df.ix[i]['c_seg'] ]['segname'].values, data_df.ix[ data_df.ix[i]['c_seg'] ]['c_seg'].values) 

#print( data_df.ix[ data_df.ix[5200]['c_seg'] ]['c_seg'])
#print( data_df.ix[ data_df.ix[5200]['c_seg'] ]['c_seg'].values)
#print( data_df.ix[ data_df.ix[5200]['c_seg'] ]['c_seg'].values[0])
#sys.exit(0)


## Get the correct arrangement of each E*, p_seg and n_seg
data_df['n_seg'] = 'NA'
data_df['p_seg'] = 'NA'
E_index = {}
for i in data_df[ data_df['segname'].str.contains('E[0-9]*') ].index.values:
    seg = data_df.ix[i]['segname']
    if seg not in E_index:
        E_index[seg] = np.array([], dtype=np.int32)
    E_index[seg] = np.append(E_index[seg], i)


for k1, v1 in E_index.items():
    E_index[k1] = trace_seg(v1)


## construct E_df
E_df = pd.DataFrame(index=np.sort(list(E_index.keys())), columns=['indice'])
for k, v in E_index.items():
    E_df.ix[k, 'indice'] = v

# calculate 'vec' of each E*
E_df['vec'] = 'NA'
for i in E_df.index.values:
    vec = unitVec(data_df.ix[ E_df.ix[i]['indice'][0] ]['xyz'] - data_df.ix[ E_df.ix[i]['indice'][1] ]['xyz'])
    E_df.ix[i, 'vec'] = vec

#print(E_df)

# group E* by 'vec'
E_df['para_group'] = 'NA'
group = 0
while len(E_df[ E_df['para_group'] == 'NA' ].index.values) != 0:
    unGroup = E_df[ E_df['para_group'] == 'NA' ].index.values
    if len(unGroup) == 1:
        E_df.ix[ unGroup[0], 'para_group'] = group
        break
        
    i = unGroup[0]
    
    tmp_group = np.array([i])
    for j in unGroup[1:]:
    
        # if they are parallel
        if np.absolute(np.absolute(np.dot(E_df.ix[i]['vec'], E_df.ix[j]['vec'])) - 1) <= 0.01:
            tmp_group = np.append(tmp_group, j)
    
    for k in tmp_group:
        E_df.ix[k, 'para_group'] = group
    
    group += 1


#print(np.unique(E_df['para_group'].values))
#sys.exit(0)

# For each 'para_group', arbitrarily choose one end of the first DNA segment to be the 'start'
for i in np.unique(E_df['para_group'].values):

    E_model = np.sort(E_df[ E_df['para_group'] == i ].index.values)[0]
    start = E_df.ix[E_model]['indice'][0]
    for j in E_df[ E_df['para_group'] == i ].index.values:
        # set the end close to 'start' to be the head
        if E_df.ix[j]['indice'][0] != close_atom(E_df.ix[j]['indice'], start):
            E_df.ix[j, 'indice'] = E_df.ix[j]['indice'][::-1]


# update 'vec' of each E*
E_df['vec'] = 'NA'
for i in E_df.index.values:
    vec = unitVec(data_df.ix[ E_df.ix[i]['indice'][0] ]['xyz'] - data_df.ix[ E_df.ix[i]['indice'][1] ]['xyz'])
    E_df.ix[i, 'vec'] = vec


# get 'E_ind'
data_df['E_ind'] = 'NA'

for k1 in E_df.index.values:
    v1 = E_df.ix[k1]['indice']

    for i in range(len(v1)):
        if i == 0:
            data_df.ix[ v1[i], 'n_seg'] = v1[i + 1] 
        elif i == len(v1)-1:
            data_df.ix[ v1[i], 'p_seg'] = v1[i - 1] 
        else:
            data_df.ix[ v1[i], 'n_seg'] = v1[i + 1] 
            data_df.ix[ v1[i], 'p_seg'] = v1[i - 1] 

        data_df.ix[ v1[i], 'E_ind'] = i

#print(data_df)
#print(data_df[ data_df['segname'] == 'T22' ])
#print(data_df.ix[802])
#print(data_df.ix[5545])
#print(data_df.ix[5595])
#print(data_df['E_ind']) 
#sys.exit(0)


# find nearest neighbor
E_df['nn'] = 'NA'
for i in E_df.index.values:
    for j in E_df.ix[i]['indice']:

        if str(data_df.ix[j]['c_seg']) != 'NA':
            for k in data_df.ix[j]['c_seg']:

                # if they are parallel
                if E_df.ix[i]['para_group'] == E_df.ix[ data_df.ix[k]['segname'] ]['para_group']:
                    if isinstance(E_df.ix[i]['nn'], str):
                        E_df.ix[i, 'nn'] = np.array([data_df.ix[k]['segname'], j, data_df.ix[j]['E_ind'], k, data_df.ix[k]['E_ind']])
                    else:
                        E_df.set_value(i, 'nn', np.append(E_df.ix[i]['nn'], [data_df.ix[k]['segname'], j, data_df.ix[j]['E_ind'], k, data_df.ix[k]['E_ind']]))


for i in E_df.index.values:
    if isinstance(E_df.ix[i]['nn'], np.ndarray):
        tmp = np.reshape(E_df.ix[i]['nn'], (int(len(E_df.ix[i]['nn']) / 5), 5))
        E_df.ix[i, 'nn'] = tmp[ np.argsort(tmp[:,0]) ]

#print(E_df) 
#sys.exit(0)


## check if there is any nucleotide without 'seq'
count = 0
for i in data_df[ data_df['segname'].str.contains('T[0-9]*') ].index.values:
    if data_df.ix[i]['seq'] == 'NA':
        #print(i)
        if data_df.ix[i]['b_nt'] != 'NA' and data_df.ix[ data_df.ix[i]['b_nt'] ]['seq'] != 'NA':
            data_df.ix[i, 'seq'] = wc(data_df.ix[ data_df.ix[i]['b_nt'] ]['seq'])

        else:
            print('Sequence not found: ', end="")
            print('%d in segname %s, '% (i, data_df.ix[i]['segname']), end="")
            print('replace with T')
            data_df.ix[i, 'seq'] = 'T'
            count += 1

print('%d nucleoditides are guessed.' % count)


## fill in the coordinates
## start with double-stranded region
data_df['xyz_aa'] = 'NA'
for i in data_df[ data_df['segname'].str.contains('E[0-9]*') ].index.values:

    # if double-strand
    if len(data_df.ix[i]['ntBond']) == 2 and isinstance(data_df.ix[i]['ntBond'], np.ndarray):

        # find coorsponding seq
        ind_1 = data_df.ix[i]['ntBond'][0]
        ind_2 = data_df.ix[i]['ntBond'][1]
        #print(ind_1, ind_2)
        seq_1 = data_df.ix[ind_1]['seq']
        seq_2 = data_df.ix[ind_2]['seq']

        # find coorsponding base-pair name
        for bp_key in bp.keys():
            if re.match(seq_1, bp_key):
                seq_12 = bp_key

        # translation
        bp_m = translation_bp(i, bp[seq_12])            

        # rotation
        bp_m_r = rotation_bp(i, ind_1, seq_1, ind_2, seq_2, bp_m)            

        # assign coordinates to each residue in DnaStrand
        data_df.set_value(ind_1, 'xyz_aa', bp_m_r[seq_1])
        data_df.set_value(ind_2, 'xyz_aa', bp_m_r[seq_2])
            

## then, single-stranded region
for i in data_df[ data_df['segname'].str.contains('T[0-9]*') ].index.values:

    if data_df.ix[i]['xyz_aa'] == 'NA':

        # Translation
        pdb_m = translation_nt(i) 
    
        # Rotation
        #print(i)
        data_df.set_value(i, 'xyz_aa', rotation_nt(i, pdb_m)) 

#print(data_df)
#sys.exit(0)


#############
## output pdb
#############

## find the 5' of each DNA strand
seg_head = {}
for i in data_df[ data_df['p_nt'] == 'NA' ].index.values:

    if re.search('T[0-9]*', data_df.ix[i]['segname']):
        seg_head[ data_df.ix[i]['segname'] ] = i

seg_head_key = np.sort(list(seg_head.keys()))
#print(seg_head)
#print(seg_map)
#sys.exit(0)

chain = 64
index = 0
residue = 0
data_df['residue'] = 'NA'

with open(sys.argv[2], 'w') as out:

    for seg in seg_head_key:
        if chain == 90:
            chain = 65
        else:
            chain += 1

        resid = 0
        rid = seg_head[seg]
        while rid in data_df.index.values:

            resid += 1
            seq = data_df.ix[rid]['seq']
 
            if data_df.ix[rid]['xyz_aa'] != 'NA':

                residue += 1
                data_df.ix[rid, 'residue'] = residue

                for aname, coor in data_df.ix[rid]['xyz_aa'].items():
                    index += 1
                    if index <= 99999:
                        out.write("%-6s%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" % ("ATOM", index, aname, nt1to3(seq), chr(chain), resid, coor[0], coor[1], coor[2], 1,0, seg))
                    else:
                        out.write("%-6s%5s %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" % ("ATOM", '*****', aname, nt1to3(seq), chr(chain), resid, coor[0], coor[1], coor[2], 1,0, seg))

            # if there is next nucleotide, continue, else, break
            if 'n_nt' in data_df.ix[rid]:
                rid = data_df.ix[rid]['n_nt']        
            else:
                break

print('Finished writing %s' % sys.argv[2])

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


####################################################
# inter-helical push 
# usage:
# > cat push.helix.for.make_ndx  | make_ndx  -f hextube.pdb -o tmp.ndx
# > awk '/^[0-9]/ && NF == 2 {printf "bond %8d %8d %8d %8d\n", $1-1, $2-1, 1, 30}' tmp.ndx > push.extrabonds

# Number of base-pairs within HJ which would be excluded from extrabonds
cutOff = 7 

push_done = []
#print(E_df) 

with open('push.helix.for.make_ndx', 'w') as FH:
    FH.write("del 0 - 30\n")

    #print(E_df.ix['E6']['nn'])
    for seg1 in E_df.index.values:
    #seg1 = 'E32'

        if isinstance(E_df.ix[seg1]['nn'], np.ndarray):
            segs = np.unique(E_df.ix[seg1]['nn'][:,0])
            #print(segs)

            for seg2 in segs:

                # if this pair has not been considered
                if (seg1, seg2) not in push_done and (seg2, seg1) not in push_done:
            
                    row, col = np.where(E_df.ix[seg1]['nn'] == seg2)
                    #print(len(row))
                    
                    tmp = E_df.ix[seg1]['nn'][row]
                    crossOver = tmp[ np.argsort( tmp[:,2].astype(int) ) ]
                    #print(crossOver)
                    
                    push_array = np.array([], dtype=np.uint16)
            
            
                    # for residues before crossovers
                    for a, b in zip(range(int(crossOver[0, 2]) - cutOff, -1, -1), range(int(crossOver[0, 4]) - cutOff, -1, -1)):
                    
                        push_array = np.append(push_array, [E_df.ix[seg1]['indice'][a], a, E_df.ix[ seg2 ]['indice'][b], b])
            
            
                    # for residues between crossovers
                    for i in range(len(crossOver) - 1):
                    
                        for a, b in zip(range(int(crossOver[i, 2]) + cutOff, int(crossOver[i+1, 2]) - cutOff), range(int(crossOver[i, 4]) + cutOff, int(crossOver[i+1, 4]) - cutOff)):
                            push_array = np.append(push_array, [E_df.ix[seg1]['indice'][a], a, E_df.ix[ seg2 ]['indice'][b], b])
                    
                            
                    # for residues after crossovers
                    for a, b in zip(range(int(crossOver[-1, 2]) + cutOff, len(E_df.ix[seg1]['indice'])), range(int(crossOver[-1, 4]) + cutOff, len(E_df.ix[ seg2 ]['indice']))):
                    
                        push_array = np.append(push_array, [E_df.ix[seg1]['indice'][a], a, E_df.ix[ seg2 ]['indice'][b], b])
                    
                    push_array = np.reshape(push_array, (int(len(push_array)/4), 4))
            
                    #print(push_array)
            
            
                    # print to output
                    for i in push_array:
                        if isinstance(data_df.ix[ i[0] ]['ntBond'], np.ndarray) and\
                           len(data_df.ix[ i[0] ]['ntBond']) == 2 and\
                           isinstance(data_df.ix[ i[2] ]['ntBond'], np.ndarray) and\
                           len(data_df.ix[ i[2] ]['ntBond']) == 2:
                            for j, k in zip(data_df.ix[ i[0] ]['ntBond'], data_df.ix[ i[2] ]['ntBond'][::-1]):
            
                                if isinstance(data_df.ix[j]['residue'], int) and isinstance(data_df.ix[k]['residue'], int):
                                    FH.write('ri %d %d & a P\n' % (data_df.ix[j]['residue'], data_df.ix[k]['residue']))


                    # after it is done, add this pair to 'push_done'
                    if len(push_done) == 0:
                        push_done = [(seg1, seg2), (seg2, seg1)]
                    else:
                        push_done.append((seg1, seg2))
                        push_done.append((seg2, seg1))

    FH.write('q\n')

print('Finished writing push.helix.for.make_ndx')
#print('Nearest neighbor DNA segments considered for writing "push.helix.for.make_ndx":')
#print(push_done)


####################################################
# chickenwire (CW)

data_df['CW'] = 'NA'

with open('chickenwire.for.make_ndx', 'w') as FH:

    FH.write("del 0 - 30\n")
    FH.write("!a C1' C2'  C3' C4' O4' C5'  P  O1P  O2P  O3'  O5' H*\n")
    FH.write("name 0 BASE\n")

    CWatom_num   = 0
    CWatoms      = []
    CWres        = []
    CWseg        = []
    CWatoms_hash = {}
    CWbonds      = []

    # iterate strand by strand
    for seg in seg_head_key:
        rid = seg_head[seg]

        while rid in data_df.index.values: 
        
            ri1 = rid
            if data_df.ix[ri1]['b_nt'] != 'NA':
                ri2 = data_df.ix[ri1]['b_nt']
            else:
                ri2 = ''

            # if the coorsponding CW atom has not been made
            if str(data_df.ix[ri1]['CW']) == 'NA':

                # if it has an associated DNA segment E*
                if str(data_df.ix[ri1]['E_seg']) != 'NA':
                    E_seg_index = data_df.ix[ri1]['E_seg']
                    E_seg_name = data_df.ix[E_seg_index]['segname']
                    E_seg_E_ind = data_df.ix[E_seg_index]['E_ind']

                    sname = "%s" % E_seg_name
                    rname = "B%s" % E_seg_E_ind
                    aname = "L0"

                    CWatom_num += 1
                    CWatoms.append(aname)
                    CWres.append(rname)
                    CWseg.append(sname)
                    key = "%s%s%s" % (sname, rname, aname)

                    if key in CWatoms_hash:
                        print('Warning: CWatoms_hash %s: %s is going to be replaced by %d' % (key, CWatoms_hash[key], CWatom_num))

                    CWatoms_hash[key] = CWatom_num

                    data_df.ix[ri1, 'CW'] = key
                    if ri2 != '' and data_df.ix[ri2]['CW'] == 'NA':
                        data_df.ix[ri2, 'CW'] = key
                        FH.write("ri %s %s & 0\n" % (data_df.ix[ri1]['residue'], data_df.ix[ri2]['residue']) )
                        FH.write("name %d %s\n" % (CWatom_num, key))

                    if ri2 == '':
                        FH.write("ri %s %s & 0\n" % (data_df.ix[ri1]['residue'], '') )
                        FH.write("name %d %s\n" % (CWatom_num, key))

                else:
                    if 'E_seg_index' in locals():
                        del E_seg_index

                    if data_df.ix[ri1]['p_nt'] != 'NA':
                        ri_temp = data_df.ix[ri1]['p_nt']
                        count = 1
                        while ri_temp in data_df.index.values:
                            if data_df.ix[ri_temp]['E_seg'] != 'NA':
                                E_seg_index = data_df.ix[ri_temp]['E_seg']
                                E_seg_name = data_df.ix[E_seg_index]['segname']
                                E_seg_E_ind = data_df.ix[E_seg_index]['E_ind']
                                break

                            if data_df.ix[ri_temp]['p_nt'] != 'NA':
                                ri_temp = data_df.ix[ri_temp]['p_nt']
                                count += 1
                            else:
                                break

                    if 'E_seg_index' not in locals() and data_df.ix[ri1]['n_nt'] != 'NA':
                        ri_temp = data_df.ix[ri1]['n_nt']
                        count = -1
                        while ri_temp in data_df.index.values:
                            if data_df.ix[ri_temp]['E_seg'] != 'NA':
                                E_seg_index = data_df.ix[ri_temp]['E_seg']
                                E_seg_name = data_df.ix[E_seg_index]['segname']
                                E_seg_E_ind = data_df.ix[E_seg_index]['E_ind'] 
                                break

                            if data_df.ix[ri_temp]['n_nt'] != 'NA':
                                ri_temp = data_df.ix[ri_temp]['n_nt']
                                count -= 1
                            else:
                                break


                    sname = "%s" % E_seg_name
                    rname = "B%s" % E_seg_E_ind
                    aname = "L%d" % count

                    CWatom_num += 1
                    CWatoms.append(aname)
                    CWres.append(rname)
                    CWseg.append(sname)
                    key = "%s%s%s" % (sname, rname, aname)

                    if key in CWatoms_hash:
                        print('Warning: CWatoms_hash %s: %s is going to be replaced by %d' % (key, CWatoms_hash[key], CWatom_num))

                    CWatoms_hash[key] = CWatom_num

                    # if E_seg does not exist, it must be ssDNA
                    data_df.ix[ri1, 'CW'] = key
                    FH.write("ri %s %s & 0\n" % (data_df.ix[ri1]['residue'], '') )
                    FH.write("name %d %s\n" % (CWatom_num, key))



            # if there is next nucleotide, continue, else, break
            if 'n_nt' in data_df.ix[rid]:
                rid = data_df.ix[rid]['n_nt']        
            else:
                break

    FH.write("quit")

print('Finished writing chickenwire.for.make_ndx')
#for i in data_df[ data_df['segname'].str.contains('T[0-9]*') ].index.values:
#    if data_df.ix[i]['CW'] == 'NA':
#        print(data_df.ix[i])


#################
# ib-to-ib bonds.
# only with next one to avoid double counting.
for seg in seg_head_key:
    rid = seg_head[seg]

    while rid in data_df.index.values: 
    
        ri1 = rid
        if data_df.ix[ri1]['n_nt'] != 'NA':
            txt1 = data_df.ix[ri1]['CW']
            txt2 = data_df.ix[ data_df.ix[ri1]['n_nt'] ]['CW'] 
            CWbonds.append([txt1, txt2])

        # if there is next nucleotide, continue, else, break
        if 'n_nt' in data_df.ix[rid]:
            rid = data_df.ix[rid]['n_nt']        
        else:
            break


# for debug
prefix = sys.argv[2][:-4]

CWbonds_F = prefix + '_CWbonds' 
with open(CWbonds_F, 'w') as FH:
    FH.write('%s' % str(CWbonds))
print('Finished writing %s' % CWbonds_F)

CWatoms_hash_F = prefix + '_CWatoms_hash' 
with open(CWatoms_hash_F, 'w') as FH:
    FH.write('%s' % str(CWatoms_hash))
print('Finished writing %s' % CWatoms_hash_F)


## write psf
with open('chickenwire.psf', 'w') as FPSF:
    FPSF.write("PSF EXT\n")
    FPSF.write("%10d !NTITLE\n" % 1)
    FPSF.write("* Chickenwire representation made by %s EXT\n\n" % sys.argv[0])

    FPSF.write("%10d !NATOM\n" % (len(CWatoms)) )

    index=0
    #resid=1
    for a in range(len(CWatoms)):
        index += 1
        FPSF.write("%10d " % index)
        FPSF.write("%-8s " % CWseg[a])
        FPSF.write("%-8d " % index)
        FPSF.write("%-8s " % CWres[a])
        FPSF.write("%-6s " % CWatoms[a])
        FPSF.write("%6d "  % 6)
        FPSF.write("  0.000000       1.00800")
        FPSF.write("           0   0.00000     -0.301140E-02\n")

    ## replace CWbonds using hash
    CWbonds_num = []
    for i in range(len(CWbonds)):

        if (CWbonds[i][0] not in CWatoms_hash):
            sys.stderr.write("Undefined %s\n" % CWbonds[i][0])
            continue

        if (CWbonds[i][1] not in CWatoms_hash):
            sys.stderr.write("Undefined %s\n" % CWbonds[i][1])
            continue

        bi = CWatoms_hash[ CWbonds[i][0] ]
        bj = CWatoms_hash[ CWbonds[i][1] ]
        if (bi < 1 or bi > len(CWatoms)+1):
            sys.stderr.write("Bond(%s,%s) out of range\n" % (bi,bj))
        if (bj < 1 or bj > len(CWatoms)+1):
            sys.stderr.write("Bond(%s,%s) out of range\n" % (bi,bj))
        CWbonds_num.append([bi,bj])

    FPSF.write("\n%10d !NBOND: bonds\n" % len(CWbonds_num))
    count = 0
    for i in CWbonds_num:
        for j in i:
            FPSF.write("%10d" % j)
            count+=1
            if count % 8 == 0:
                FPSF.write("\n")

    FPSF.write("\n")

    FPSF.write("\n%10d !NTHETA: angles\n" % 0)
    FPSF.write("\n%10d !NPHI: dihedrals\n" % 0)
    FPSF.write("\n%10d !NIMPHI: impropers\n" % 0)
    FPSF.write("\n%10d !NDON: donors\n" % 0)
    FPSF.write("\n%10d !NACC: acceptors\n" % 0)

print('Finished writing chickenwire.psf')


####################################################
# export temporary data

prefix = sys.argv[2][:-4]

data_df_F = prefix + '_data_df.csv'
data_df.to_csv(data_df_F)
print('Export data_df to %s' % data_df_F)

E_df_F = prefix + '_E_df.csv'
E_df.to_csv(E_df_F)
print('Export E_df to %s' % E_df_F)

seg_map_F = prefix + '_seg_map.csv'
with open(seg_map_F, 'w') as FH:
    for k in np.sort(list(seg_map.keys())):
        FH.write('%s,%s\n' % (k, seg_map[k]))

print('Export seg_map to %s' % seg_map_F)




