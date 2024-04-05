#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 20:33:36 2022

@author: luoxini
"""

from gurobipy import *
import gurobipy as gp
import time
import math
import random

from data_generator import Data

'''
# auxiliary parameters
max_col = 10
max_box = 100
max_bin = 50

total_bin_type = 5
total_box_type = 10
total_bay_type = 3

h_space = 1.5
v_space = 1.5
# bin parameter
#bin_types = tuplelist(range(total_bin_type))

# bay_parameter
#bay_types = tuplelist(range(total_bay_type))
bay_width=120
bay_height=260

# boxes parameters
box_types, box_width, box_height, box_weight, box_quantities = multidict(
    {'big_pizza':[34.1,5.2,15,100],'small_pizza':[18.6,2.2,10,100],
    'tube':[10.9,4.3,10,100],'tray':[27.6,3.1,15,100]})  
    # weight,quantities?
box_types, box_width, box_height, box_weight, box_quantities = multidict(
    {0:[34.1,5.2,30,200],1:[18.6,2.2,20,200],2:[10.9,4.3,20,200],3:[27.6,3.1,30,200]})
# beam parameters
beam_weight_limit = 500
beam_height = 4

# model
model = gp.Model('nonlinear warehouse model')
# add decision variables
y = model.addVars(total_bin_type, total_box_type, vtype=GRB.BINARY, name='y')
n = model.addVars(total_bin_type, total_box_type, lb=0, vtype=GRB.INTEGER, name='n')
m = model.addVars(total_bin_type, total_box_type, lb=0, vtype=GRB.INTEGER, name='m')

x = model.addVars(total_bin_type, lb=0, vtype=GRB.INTEGER, name='x')
c = model.addVars(total_bin_type, lb=0, vtype=GRB.INTEGER, name='c')
b = model.addVars(total_bin_type, lb=0, vtype=GRB.INTEGER, name='b')

z = model.addVars(total_bin_type, total_bay_type, vtype=GRB.BINARY, name='z')
d = model.addVars(total_bin_type, total_bay_type, lb=0, vtype=GRB.INTEGER, name='d')
o = model.addVars(total_bay_type, lb=0, vtype=GRB.INTEGER, name='o')

# set objective
model.setObjective(o.sum(), GRB.MINIMIZE)


# add constraints
model.addConstrs((n[i,j] <= max_col * y[i,j] for i in range(total_bin_type) for j in range(total_box_type)), name='1')
model.addConstrs((n.sum(i, '*') == c[i] for i in range(total_bin_type)), name='2')
model.addConstrs((m[i,j] <= max_box * y[i,j] for i in range(total_bin_type) for j in range(total_box_type)), name='3')
model.addConstrs((m.sum(i,'*') <= max_box for i in range(total_bin_type)), name='4')
model.addConstrs((d[i,t] <= max_bin * z[i,t] for i in range(total_bin_type) for t in range(total_bay_type)),name='5')
model.addConstrs((d.sum('*',t) <= max_bin for t in range(total_bay_type)),name='6')
model.addConstrs((gp.quicksum(n[i,j]*box_width[j] for j in range(total_box_type))+h_space*(c[i]+1) <= bay_width for i in range(total_bin_type)),name='7')
model.addConstrs((gp.quicksum(m[i,j]*box_weight[j] for j in range(total_box_type)) <= beam_weight_limit for i in range(total_bin_type)), name='8')
model.addConstrs(((x[i]-v_space)*n[i,j]>=m[i,j]*box_height[j] for i in range(total_bin_type) for j in range(total_box_type)),name='9')
model.addConstrs((gp.quicksum(d[i,t]*(x[i]+beam_height) for i in range(total_bin_type)) == bay_height for t in range(total_bay_type)),name='10')
model.addConstrs((gp.quicksum(m[i,j]*b[i] for i in range(total_bin_type))==box_quantities[j] for j in range(total_box_type)),name='11')
model.addConstrs((gp.quicksum(d[i,t]*o[t] for t in range(total_bay_type))==b[i] for i in range(total_bin_type)),name='12')

# 7.1 
model.addConstrs((c[i] <= max_col for i in range(total_bin_type)),name='col limit')
model.addConstrs((n[i,j]>=y[i,j] for i in range(total_bin_type) for j in range(total_box_type)),name='n>=y')
model.addConstrs((m[i,j]>=y[i,j] for i in range(total_bin_type) for j in range(total_box_type)),name='m>=y')
model.addConstrs((d[i,t]>=z[i,t] for i in range(total_bin_type) for t in range(total_bay_type)),name='d>=z')
model.addConstrs((b[i]>=gp.quicksum(y[i,j] for j in range(total_box_type)) for i in range(total_bin_type)),name='b>=sum_y')
model.addConstrs((gp.quicksum(y[i,j] for j in range(total_box_type))>=1 for i in range(total_bin_type)),name='sum_y>=1')
model.addConstrs((o[t]>=gp.quicksum(z[i,t] for i in range(total_bin_type)) for t in range(total_bay_type)),name='o>=sum_z')
model.addConstrs((gp.quicksum(z[i,t] for i in range(total_bin_type))>=1 for t in range(total_bay_type)),name='sum_z>=1')

model.update()
model.write('nonlinear.lp')
#model.Params.TimeLimit = 60

# Solve
try:
    model.optimize()
except gp.GurobiError:
    print("Optimize failed due to non-convexity")

model.Params.NonConvex = 2
model.optimize()
'''


class NIP:
    def __init__(self, I, T, data):
        self.I = I  # total bin type
        # max_col_num = math.floor((data.W-data.a)/(min(data.w_j)+data.a))

        self.T = T  # total bay type
        # max_bin_num = math.floor(data.H/(min(data.h_j)+data.a))
        self.data = data
        self.model = gp.Model('nonlinear')
        self.o_t = {}

    def construct_NIP(self):
        # model
        # model = gp.Model('nonlinear')

        # decision varaibles
        x_i = self.model.addVars(self.I, lb=0, vtype=GRB.INTEGER, name='x')  # bin height
        c_i = self.model.addVars(self.I, lb=0, vtype=GRB.INTEGER, name='c')  # column num
        b_i = self.model.addVars(self.I, lb=0, vtype=GRB.INTEGER, name='b')  # used number of bin i

        m_ij = self.model.addVars(self.I, self.data.J, lb=0, vtype=GRB.INTEGER, name='m')  # stored column of type j
        n_ij = self.model.addVars(self.I, self.data.J, lb=0, vtype=GRB.INTEGER,
                                  name='n')  # stored quant in every column of type j
        # n_ij = [[0 for j in range(J)] for i in range(I)]
        # for i in range(I):
        # for j in range(J):
        # n_ij[i][j] = int(f_j[j]/e_j[j])+1
        p_ij = self.model.addVars(self.I, self.data.J, lb=0, vtype=GRB.INTEGER, name='p')  # p = m*n

        d_it = self.model.addVars(self.I, self.T, lb=0, vtype=GRB.INTEGER,
                                  name='d')  # quant of bin type i in bay type t

        for t in range(self.T):
            self.o_t[t + 1] = self.model.addVar(lb=0, vtype=GRB.INTEGER, name='o')  # used number of bay t

        # objective
        self.model.setObjective(quicksum(self.o_t[t+1] for t in range(self.T)), GRB.MINIMIZE)

        self.model.update()
        # consraints
        self.model.addConstrs((gp.quicksum(m_ij[i, j] for j in range(self.data.J)) == c_i[i] for i in range(self.I)),
                              name='column_number_equal_cons')
        self.model.addConstrs((quicksum((self.data.w_j[j] * m_ij[i, j]) for j in range(self.data.J)) + self.data.a * (
                c_i[i] + 1) <= self.data.W for i in range(self.I)), name='width_cons')
        self.model.addConstrs(
            (self.data.h_j[j] * n_ij[i, j] + self.data.b <= x_i[i] for i in range(self.I) for j in range(self.data.J)),
            name='quant_of_one_col_and_bin_height_cons')
        self.model.addConstrs(((n_ij[i, j] - 1) * self.data.e_j[j] <= self.data.f_j[j] for i in range(self.I) for j in
                               range(self.data.J)), name='max_allow_stack_cons')
        self.model.addConstrs(
            (quicksum(d_it[i, t] * (x_i[i] + self.data.B) for i in range(self.I)) <= self.data.H for t in
             range(self.T)), name='bay_height_cons')
        self.model.addConstrs(
            (quicksum(p_ij[i, j] * b_i[i] for i in range(self.I)) >= self.data.q_j[j] for j in range(self.data.J)),
            name='SUK_total_quant_cons')
        self.model.addConstrs(
            (quicksum(d_it[i, t] * self.o_t[t + 1] for t in range(self.T)) == b_i[i] for i in range(self.I)),
            name='used_bin_cons')

        self.model.addConstrs(p_ij[i, j] == m_ij[i, j] * n_ij[i, j] for i in range(self.I) for j in range(self.data.J))

        self.model.update()
        self.model.write('nonlinear model.lp')
        return

    def solve_NIP(self):
        start = time.time()
        self.model.Params.NonConvex = 2
        self.model.Params.TimeLimit = 1800
        self.model.optimize()
        end = time.time()
        lb = self.model.getAttr('ObjBoundC')
        ub = self.model.getAttr('ObjVal')
        return lb, ub, end-start
        # if self.model.Status == 2:
        #     return True, self.model.ObjVal, end-start
        # else:
        #     return False, None, 1800

    def check_feasible(self, a):
        # self.model.remove(self.model.getConstrByName('check_feasible'))
        print(self.model.getAttr('NumConstrs'))
        self.model.addConstr(quicksum(self.o_t[t+1] for t in range(self.T)) <= a, name='check_feasible')
        print(self.model.getAttr('NumConstrs'))
        #self.model.getConstrByName('check_feasible').setAttr('rhs', a)
        self.model.setObjective(0)

        self.model.update()
        self.model.write('nonlinear model.lp')
        self.model.Params.NonConvex = 2
        self.model.Params.SolutionLimit = 1
        self.model.optimize()

        if self.model.Status == 3: # infeasible
            return False  # infeasible
        else:
            return True  # feasible

    def check_feasible_range(self,a_1,a_o):
        left = a_1
        right = a_o
        while left < right:
            mid = math.floor(left + (right - left)/2)
            if self.check_feasible(mid): # feasible
                right = mid
            else:
                left = mid + 1
        return left



