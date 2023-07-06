import numpy as np
import os
from Geometry3D import *
from sympy import Point3D, Line, Segment3D, Plane
import math
import random
from numba import *
from scipy.linalg import eigh
from scipy import linalg
from scipy.optimize import minimize_scalar

# Dimensions of the RVE
L = 100.0
W = 100.0
B = 100.0

# Distribution of length
# generate log-normal distributed random variable
b = 1
w = 5
volume_fraction_total = 0.002
volume_fraction = 0
plate_count = 0
length_lognormal = []

while volume_fraction < volume_fraction_total:
    fresh_lognormal_length = np.random.lognormal(mean=math.log(5, math.exp(1)), sigma=math.log(2, math.exp(1)), size=None)
    if fresh_lognormal_length > 9 or fresh_lognormal_length < 1:
        continue
    volume_fraction += (fresh_lognormal_length*b*w) / (L*B*W)
    plate_count += 1
    length_lognormal.append(fresh_lognormal_length)

    if volume_fraction > volume_fraction_total:
        volume_fraction -= (fresh_lognormal_length*b*w) / (L*B*W)
        plate_count -= 1
        length_lognormal.pop()
        break
print(len(length_lognormal))
np.savetxt('plate_lengths.csv', length_lognormal, delimiter=',')

# Final Position & Orientation
max_incl = len(length_lognormal)
num_incl = 0
alpha = []
beta = []
gamma = []
x_coordinate = []
y_coordinate = []
z_coordinate = []

# Rotation matrix
def matrix_11(alpha, beta, gamma):
    return math.cos(alpha) * math.cos(beta)
def matrix_12(alpha, beta, gamma):
    return math.cos(alpha) * math.sin(beta) * math.sin(gamma) - math.sin(alpha) * math.cos(gamma)
def matrix_13(alpha, beta, gamma):
    return math.cos(alpha) * math.sin(beta) * math.cos(gamma) + math.sin(alpha) * math.sin(gamma)
def matrix_21(alpha, beta, gamma):
    return math.sin(alpha) * math.cos(beta)
def matrix_22(alpha, beta, gamma):
    return math.sin(alpha) * math.sin(beta) * math.sin(gamma) + math.cos(alpha) * math.cos(gamma)
def matrix_23(alpha, beta, gamma):
    return math.sin(alpha) * math.sin(beta) * math.cos(gamma) - math.cos(alpha) * math.sin(gamma)
def matrix_31(alpha, beta, gamma):
    return -math.sin(beta)
def matrix_32(alpha, beta, gamma):
    return math.cos(beta) * math.sin(gamma)
def matrix_33(alpha, beta, gamma):
    return math.cos(beta) * math.cos(gamma)

def rotation_matrix(alpha, beta, gamma):
    return [[matrix_11(alpha, beta, gamma), matrix_12(alpha, beta, gamma), matrix_13(alpha, beta, gamma)],
            [matrix_21(alpha, beta, gamma), matrix_22(alpha, beta, gamma), matrix_23(alpha, beta, gamma)],
            [matrix_31(alpha, beta, gamma), matrix_32(alpha, beta, gamma), matrix_33(alpha, beta, gamma)]]

# Creating 8 Vertices
# confidence_offset = 0.5
# def vertex_1(w,l,b):
#     return [(w+confidence_offset)/2, (l+confidence_offset)/2, (b+confidence_offset)/2]
# def vertex_2(w,l,b):
#     return [-(w+confidence_offset)/2, (l+confidence_offset)/2, (b+confidence_offset)/2]
# def vertex_3(w,l,b):
#     return [-(w+confidence_offset)/2, (l+confidence_offset)/2, -(b+confidence_offset)/2]
# def vertex_4(w,l,b):
#     return [(w+confidence_offset)/2, (l+confidence_offset)/2, -(b+confidence_offset)/2]
# def vertex_5(w,l,b):
#     return [(w+confidence_offset)/2, -(l+confidence_offset)/2, (b+confidence_offset)/2]
# def vertex_6(w,l,b):
#     return [-(w+confidence_offset)/2, -(l+confidence_offset)/2, (b+confidence_offset)/2]
# def vertex_7(w,l,b):
#     return [-(w+confidence_offset)/2, -(l+confidence_offset)/2, -(b+confidence_offset)/2]
# def vertex_8(w,l,b):
#     return [(w+confidence_offset)/2, -(l+confidence_offset)/2, -(b+confidence_offset)/2]


# intersect_flag = False
mu_A_list_all = []
u1_list_all = []
u2_list_all = []
u3_list_all = []


def P(u1, u2, u3):
    return np.array([[u1[0], u2[0], u3[0]],
                     [u1[1], u2[1], u3[1]],
                     [u1[2], u2[2], u3[2]]])


def D(l, b, w):
    return np.array([[2 / (l ** 2), 0, 0],
                    [0, 2 / (b ** 2), 0],
                    [0, 0, 2 / (w ** 2)]])


def A(P, D):
    return np.dot(np.dot(P, D), P.T)


def A_inv(A):
    return np.linalg.inv(A)


def intersection(mu_A, mu_B, Sigma_A, Sigma_B):
    v = mu_B-mu_A
    def k(x):
        return 1 - np.dot(np.dot(v, np.linalg.inv((1 / (1 - x)) * Sigma_B + (1 / x) * Sigma_A)), v.T)
    res = minimize_scalar(k, bounds=(0.001, 1), method='bounded')
    return k(res.x) > 0

while (num_incl < max_incl):
    random_alpha = random.uniform(0, 180)
    random_alpha_radian = random_alpha * math.pi/180
    random_beta = random.uniform(0, 180)
    random_beta_radian = random_beta * math.pi / 180
    random_gamma = random.uniform(0, 180)
    random_gamma_radian = random_gamma * math.pi / 180
    random_x = random.uniform(L/10, 9*L/10)
    random_y = random.uniform(W/10, 9*W/10)
    random_z = random.uniform(B/10, 9*B/10)
    translation = [random_x, random_y, random_z]
    l = length_lognormal[len(x_coordinate)]

    # point_list = []
    # point_object_list = []
    #
    # point_1 = list(np.array(np.dot(rotation_matrix(random_alpha_radian, random_beta_radian, random_gamma_radian),
    #                  vertex_1(w, length_lognormal[len(x_coordinate)], b))) + np.array(translation))
    # point_list.append(point_1)
    # # print(point_1)
    # point_2 = list(np.array(np.dot(rotation_matrix(random_alpha_radian, random_beta_radian, random_gamma_radian),
    #                  vertex_2(w, length_lognormal[len(x_coordinate)], b))) + np.array(translation))
    # point_list.append(point_2)
    # point_3 = list(np.array(np.dot(rotation_matrix(random_alpha_radian, random_beta_radian, random_gamma_radian),
    #                  vertex_3(w, length_lognormal[len(x_coordinate)], b))) + np.array(translation))
    # point_list.append(point_3)
    # point_4 = list(np.array(np.dot(rotation_matrix(random_alpha_radian, random_beta_radian, random_gamma_radian),
    #                  vertex_4(w, length_lognormal[len(x_coordinate)], b))) + np.array(translation))
    # point_list.append(point_4)
    # point_5 = list(np.array(np.dot(rotation_matrix(random_alpha_radian, random_beta_radian, random_gamma_radian),
    #                  vertex_5(w, length_lognormal[len(x_coordinate)], b))) + np.array(translation))
    # point_list.append(point_5)
    # point_6 = list(np.array(np.dot(rotation_matrix(random_alpha_radian, random_beta_radian, random_gamma_radian),
    #                  vertex_6(w, length_lognormal[len(x_coordinate)], b))) + np.array(translation))
    # point_list.append(point_6)
    # point_7 = list(np.array(np.dot(rotation_matrix(random_alpha_radian, random_beta_radian, random_gamma_radian),
    #                  vertex_7(w, length_lognormal[len(x_coordinate)], b))) + np.array(translation))
    # point_list.append(point_7)
    # point_8 = list(np.array(np.dot(rotation_matrix(random_alpha_radian, random_beta_radian, random_gamma_radian),
    #                  vertex_8(w, length_lognormal[len(x_coordinate)], b))) + np.array(translation))
    # point_list.append(point_8)
    #
    # for j in range(len(point_list)):
    #     p_geometry3D = Point(point_list[j])
    #     point_object_list.append(p_geometry3D)


    # if num_incl == 0:
    #     plate_list_all.append(plate)
    #     x_coordinate.append(random_x)
    #     y_coordinate.append(random_y)
    #     z_coordinate.append(random_z)
    #     alpha.append(random_alpha)
    #     beta.append(random_beta)
    #     gamma.append(random_gamma)
    #     num_incl = num_incl + 1
    #     print(num_incl)

    # intersect_flag = False
    # for i in range(len(plate_list_all)):
    #     if intersection(plate, plate_list_all[i]) != None:
    #         intersect_flag = True
    #         break


    # if intersect_flag == False:
    #     plate_list_all.append(plate)
    #     x_coordinate.append(random_x)
    #     y_coordinate.append(random_y)
    #     z_coordinate.append(random_z)
    #     alpha.append(random_alpha)
    #     beta.append(random_beta)
    #     gamma.append(random_gamma)
    #     num_incl = num_incl + 1
    #     print(num_incl)

# # Creating CSV data for Abaqus
# params_list = [x_coordinate, y_coordinate, z_coordinate, alpha, beta, gamma]
# params_array = np.array(params_list)
# params_array_T = np.transpose(params_array)
# # params_data = open("plate_params.dat", "w")
# # for i in range(len(x_coordinate)):
# #     params_data.write('%6.6f %6.6f \n' %(params_array[0][i], params_array[1][i], params_array[2][i], params_array[3][i],
# #                                          params_array[4][i], params_array[5][i]))
# # params_data.close()
# np.savetxt('plate_params_data.csv', params_array_T, delimiter=',')

    mu_A = np.array(translation)
    u1 = list(np.dot(rotation_matrix(random_alpha_radian, random_beta_radian, random_gamma_radian), [1, 0, 0])
              + np.array(translation))
    u2 = list(np.dot(rotation_matrix(random_alpha_radian, random_beta_radian, random_gamma_radian), [0, 1, 0])
              + np.array(translation))
    u3 = list(np.dot(rotation_matrix(random_alpha_radian, random_beta_radian, random_gamma_radian), [0, 0, 1])
              + np.array(translation))
    P_new = P(u1, u2, u3)
    D_new = D(l, b, w)
    A_new = A(P_new, D_new)
    A_inv_new = A_inv(A_new)

    if num_incl == 0:
        mu_A_list_all.append(mu_A)
        u1_list_all.append(u1)
        u2_list_all.append(u2)
        u3_list_all.append(u3)
        x_coordinate.append(random_x)
        y_coordinate.append(random_y)
        z_coordinate.append(random_z)
        alpha.append(random_alpha)
        beta.append(random_beta)
        gamma.append(random_gamma)
        num_incl = num_incl + 1
        print(num_incl)

    intersect_flag = False
    for i in range(len(mu_A_list_all)):
        P_initial = P(u1_list_all[i], u2_list_all[i], u3_list_all[i])
        D_initial = D(l, b, w)
        A_initial = A(P_initial, D_initial)
        A_inv_initial = A_inv(A_initial)
        if intersection(mu_A_list_all[i], mu_A, A_inv_initial, A_inv_new) == True:
            intersect_flag = True
            break

    if intersect_flag == False:
        mu_A_list_all.append(mu_A)
        u1_list_all.append(u1)
        u2_list_all.append(u2)
        u3_list_all.append(u3)
        x_coordinate.append(random_x)
        y_coordinate.append(random_y)
        z_coordinate.append(random_z)
        alpha.append(random_alpha)
        beta.append(random_beta)
        gamma.append(random_gamma)
        num_incl = num_incl + 1
        print(num_incl)


# Creating CSV data for Abaqus
params_list = [x_coordinate, y_coordinate, z_coordinate, alpha, beta, gamma]
params_array = np.array(params_list)
params_array_T = np.transpose(params_array)
# params_data = open("plate_params.dat", "w")
# for i in range(len(x_coordinate)):
#     params_data.write('%6.6f %6.6f \n' %(params_array[0][i], params_array[1][i], params_array[2][i], params_array[3][i],
#                                          params_array[4][i], params_array[5][i]))
# params_data.close()
np.savetxt('plate_params_data.csv', params_array_T, delimiter=',')


