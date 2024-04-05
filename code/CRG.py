from gurobipy import *
import gurobipy as gp
import math
import time
import numpy as np

# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 16:37:39 2022

@author: luoxini
"""


class CRG:
    def __init__(self, data):
        self.data = data
        self.SRMP = None
        self.y_PSP = None
        self.row_PSP = None
        self.k_i = []
        self.J_constrs = []
        self.I_bar_constrs = []
        self.P_bar = 0
        self.Q_bar = 0
        self.I_bar = 0
        self.x = {}
        self.y = {}
        self.d = {}
        self.m_pj = {}
        self.pattern_d = {}
        self.y_linexpr = None
        self.y_obj_linexpr = None
        self.bay_pattern_constr = None
        self.RANDOM = True

    # 计算v_i=\sum_{j}{u_j m_j n_j}
    def CalV_i(self, u, k):  # m_pattern一直固定，不列入参数
        # u 为一个len=J的list，存储当前所要计算用到的u_j值
        # k为一个list
        v_i = {}  # dict
        for kk in k:
            best_v_under_k = [-1] * 2
            for m_pattern in self.data.m_pattern_set:
                current_v = sum(u[j] * m_pattern[j] * self.data.lambda_j_under_k[(kk, j)] for j in range(self.data.J))
                if current_v >= best_v_under_k[1]:
                    best_v_under_k[0] = m_pattern
                    best_v_under_k[1] = current_v
            v_i[(kk, tuple(best_v_under_k[0]))] = best_v_under_k[1]

        if len(v_i) == 1:
            return v_i
        # print('v_i='+str(v_i))
        '''
        v_i =   每个k下v最大的bin pattern及对应v  SRMP的dual不同会变化
        {(28, (9, 0, 0, 0)): 0.3333333333333333, 
         (23, (9, 0, 0, 0)): 0.3333333333333333, 
         (19, (9, 0, 0, 0)): 0.25, 
         (18, (9, 0, 0, 0)): 0.25, 
         (15, (9, 0, 0, 0)): 0.16666666666666666, 
         (14, (9, 0, 0, 0)): 0.16666666666666666, 
         (13, (9, 0, 0, 0)): 0.16666666666666666, 
         (12, (9, 0, 0, 0)): 0.16666666666666666, 
         (11, (9, 0, 0, 0)): 0.08333333333333333, 
         (9, (9, 0, 0, 0)): 0.08333333333333333, 
         (8, (9, 0, 0, 0)): 0.08333333333333333, 
         (7, (9, 0, 0, 0)): 0.08333333333333333, 
         (6, (0, 4, 0, 1)): 0.019772376543209874, 
         (5, (0, 0, 0, 3)): 0.013020833333333334, 
         (4, (1, 5, 0, 0)): 0.009645061728395061}
        '''

        # elimination bad bin pattern  希望k更小 v更大
        alist = []
        for k, v in v_i.items():
            alist.append([k, v])
        for i in range(len(alist)):
            alist[i].insert(0, i)
        # print(alist)
        '''
        alist = 
        [[0, (28, (9, 0, 0, 0)), 0.3333333333333333],
         [1, (23, (9, 0, 0, 0)), 0.3333333333333333],
         [2, (19, (9, 0, 0, 0)), 0.25],
         [3, (18, (9, 0, 0, 0)), 0.25],
         [4, (15, (9, 0, 0, 0)), 0.16666666666666666],
         [5, (14, (9, 0, 0, 0)), 0.16666666666666666],
         [6, (13, (9, 0, 0, 0)), 0.16666666666666666],
         [7, (12, (9, 0, 0, 0)), 0.16666666666666666],
         [8, (11, (9, 0, 0, 0)), 0.08333333333333333],
         [9, (9, (9, 0, 0, 0)), 0.08333333333333333],
         [10, (8, (9, 0, 0, 0)), 0.08333333333333333],
         [11, (7, (9, 0, 0, 0)), 0.08333333333333333],
         [12, (6, (0, 4, 0, 1)), 0.019772376543209874],
         [13, (5, (0, 0, 0, 3)), 0.013020833333333334],
         [14, (4, (1, 5, 0, 0)), 0.009645061728395061]]
        '''

        del_list = []
        for i in range(len(alist) - 1):
            if alist[i + 1][2] >= alist[i][2]:
                del_list.append(alist[i][1:])
        # print('del_list='+str(del_list))
        '''
        del_list = 
        [[(28, (9, 0, 0, 0)), 0.3333333333333333], 
         [(19, (9, 0, 0, 0)), 0.25], 
         [(15, (9, 0, 0, 0)), 0.16666666666666666], 
         [(14, (9, 0, 0, 0)), 0.16666666666666666], 
         [(13, (9, 0, 0, 0)), 0.16666666666666666], 
         [(11, (9, 0, 0, 0)), 0.08333333333333333], 
         [(9, (9, 0, 0, 0)), 0.08333333333333333], 
         [(8, (9, 0, 0, 0)), 0.08333333333333333]]
        '''
        i = 0
        if len(del_list) > 0:
            while del_list[i][0] in v_i.keys():
                del v_i[del_list[i][0]]
                i += 1
                if i == len(del_list):
                    break
            # print('v_i='+str(v_i))
        '''
        基于目前的SRMP dual结果，得到的在row generating中有必要考虑的column
        v_i = 
        {(23, (9, 0, 0, 0)): 0.3333333333333333,
         (18, (9, 0, 0, 0)): 0.25,
         (12, (9, 0, 0, 0)): 0.16666666666666666,
         (7, (9, 0, 0, 0)): 0.08333333333333333,
         (6, (0, 4, 0, 1)): 0.019772376543209874,
         (5, (0, 0, 0, 3)): 0.013020833333333334,
         (4, (1, 5, 0, 0)): 0.009645061728395061}
        '''

        return v_i

    def knapsack_dp(items, C):
        WEIGHT, VALUE = range(2)

        # order by max value per item weight
        items = sorted(items, key=lambda item: item[VALUE] / float(item[WEIGHT]), reverse=True)
        if not isinstance(C, int):
            items_copy = []
            for item in items:
                items_copy.append([int(item[0] * 10), item[1]])
            C_copy = int(C * 10)
        else:
            items_copy = items
            C_copy = C

        tab = [[0, [0 for _ in items]] for _ in range(0, C_copy + 1)]  # value, [item tab]
        for i, item in enumerate(items_copy):
            weight, value = item
            for c in range(weight, C_copy + 1):
                tabbefore = tab[c - weight]  # previous max tab to try adding this item to
                new_value = tabbefore[0] + value
                used = tabbefore[1][i]
                if tab[c][0] < new_value:
                    # old max tab with this added item is better
                    tab[c] = (new_value, tabbefore[1][:])
                    tab[c][1][i] += 1  # use one more

        value, bagged = tab[C_copy]
        numbagged = sum(bagged)
        weight = sum(items[i][0] * n for i, n in enumerate(bagged))
        # convert to (iten, count) pairs) in name order
        bagged = sorted([items_copy[i][WEIGHT], n] for i, n in enumerate(bagged) if n)
        if not isinstance(C, int):
            for bag in bagged:
                bag[0] /= 10
        return value, weight, bagged

    #############################求解x_PSP，只要找对应k下\sum(u*m*n)最大的bin pattern即可
    # 直接找到每个x_PSP的解，加入SRMP
    def add_x_column(self, which_i, pattern):  # pattern为m_pj,list len=J
        global P_bar

        x_newcolumn = Column(-1, self.SRMP.getConstrByName('I_constr_bin_height=' + str(self.k_i[which_i - 1])))
        for j in range(self.data.J):
            mj = pattern[j]
            nj = self.data.lambda_j_under_k[(self.k_i[which_i - 1], j)]
            x_newcolumn.addTerms(mj * nj, self.J_constrs[j])

        self.x[self.P_bar + 1] = self.SRMP.addVar(obj=0, lb=0, vtype=GRB.CONTINUOUS, column=x_newcolumn,
                                                  name='x' + str(self.P_bar + 1))
        self.P_bar += 1
        self.m_pj[self.P_bar] = list(pattern)

        self.SRMP.write('SRMP.lp')
        self.SRMP.update()
        self.SRMP.optimize()
        return self.SRMP.getVars(), self.SRMP.getAttr('ObjVal')

    # 更新x_PSP，判断是否还要继续True
    def update_x_PSP(self, which_i):
        items = []
        for j in range(self.data.J):
            value = self.J_constrs[j].pi * self.data.lambda_j_under_k[(self.k_i[which_i - 1], j)]
            weight = self.data.w_j[j] + self.data.a
            items.append([weight, value])
        capacity = self.data.W - self.data.a
        result = CRG.knapsack_dp(items, float(capacity))
        obj = self.I_bar_constrs[which_i - 1].pi - result[0]

        pattern_to_j = {}
        for j in range(self.data.J):
            pattern_to_j[self.data.w_j[j] + self.data.a] = j

        pattern = [0 for _ in range(self.data.J)]
        for item, num in result[2]:
            related_j = pattern_to_j[item]
            pattern[related_j] = num
        if obj < 0 and pattern not in self.m_pj.values():
            return (True, pattern)
        else:
            return (False, None)

        # u = []
        # for j in range(self.data.J):
        #     u.append(self.J_constrs[j].pi)
        # pattern = list(list(CRG.CalV_i(self, u, [self.k_i[which_i - 1]]).keys())[0][1])  # tuple
        # if list(pattern) not in self.m_pj.values():
        #     return (True, pattern)  # not in m_pj, add
        # else:
        #     return (False, None)  # in m_pj, don't add

    '''
    # 每个i有相应x_PSP
    #list 长度为I
    x_PSP_list = []
    m_list = []
    n_list = []
    h_constr_list = []
    W_constr_list = []
    E_constr_list = []
    big_M_constr_list = []


    for i in range(I_bar):
        x_PSP_list.append(gp.Model('x_PSP'+' for bin height='+str(k_i[i])))

        m_list.append(x_PSP_list[i].addVars(J, lb=0, vtype=GRB.INTEGER, name='m for bin height='+str(k_i[i])))
        n_list.append(x_PSP_list[i].addVars(J, lb=0, vtype=GRB.INTEGER, name='n for bin height='+str(k_i[i])))

        h_constr_list.append(x_PSP_list[i].addConstrs(n_list[i][j] * h_j[j] + b <= k_i[i] for j in range(J)))
        E_constr_list.append(x_PSP_list[i].addConstrs((n_list[i][j]-1) * e_j[j] <= f_j[j] for j in range(J)))
        x_PSP_list[i].setObjective(I_bar_constrs[i].pi - quicksum(m_list[i][j]*n_list[i][j]*J_constrs[j].pi for j in range(J)), GRB.MINIMIZE)
        # 上面是涉及i的，下面是不涉及i的
        W_constr_list.append(x_PSP_list[i].addLConstr(quicksum(m_list[i][j] * w_j[j] for j in range(J))+a*(quicksum(m_list[i][j] for j in range(J))+1)<=W))
        big_M_constr_list.append(x_PSP_list[i].addConstrs(n_list[i][j]<= (int(f_j[j]/e_j[j])+1) * m_list[i][j] for j in range(J)))
        # x_PSP[i].Params.NonConvex = 2
        x_PSP_list[i].update()

        for j in range(J):
            x_PSP_list[i].setAttr('ConstrName',h_constr_list[i][j],'h_constr for bin_height='+str(k_i[i])+' and j='+str(j))
            x_PSP_list[i].setAttr('ConstrName',E_constr_list[i][j],'E_constr for bin_height='+str(k_i[i])+' and j='+str(j))
            x_PSP_list[i].setAttr('ConstrName',big_M_constr_list[i][j],'big_M_constr for bin height='+str(k_i[i])+' and j='+str(j))
        x_PSP_list[i].setAttr('ConstrName',W_constr_list[i],'W_constr for bin height='+str(k_i[i]))

        x_PSP_list[i].update()
        x_PSP_list[i].optimize()
    # x_PSP_list[i]为每个i对应的x_PSP
    # m_list[i],n_list[i]为每个x_PSP的变量
    for i in range(I_bar):
        print(x_PSP_list[i].getVars())
        print(x_PSP_list[i].getAttr('ObjVal'))



    # 判断x_PSP.objval<0之后，到包括更新x_PSP参数,并解x_PSP,返回x_PSP的目标函数值
    def add_x_column(which_i):  # which_i = 1,2,3...
        global P_bar    
        x_newcolumn = Column(-1, SRMP.getConstrByName('I_constr_bin_height='+str(k_i[which_i-1])))
        for j in range(J):
            mj = m_list[which_i-1][j].getAttr('X')
            nj = n_list[which_i-1][j].getAttr('X')
            x_newcolumn.addTerms(mj*nj ,J_constrs[j])

        x[P_bar+1] = SRMP.addVar(obj=0, lb=0, vtype=GRB.CONTINUOUS, column=x_newcolumn, name='x'+str(P_bar+1))
        P_bar += 1
        SRMP.write('SRMP.lp')
        SRMP.update()
        SRMP.optimize()

        print(SRMP.getVars())
        return SRMP.getAttr('ObjVal')


        for i in range(I_bar):  # 全部的x_PSP目标函数都需要更新
            x_PSP_list[i].setObjective(I_bar_constrs[i].pi - quicksum(m_list[i][j]*n_list[i][j]*J_constrs[j].pi for j in range(J)),GRB.MINIMIZE)

        x_PSP_list[which_i-1].optimize()

        return x_PSP_list[which_i-1].getAttr('ObjVal')


    count_x = 0
    for i in range(I_bar):
        while x_PSP_list[i].getAttr('ObjVal') < 0:
            if count_x == 0:
                count_x += 1
                previous_objval = x_PSP_list[i].getAttr('ObjVal')
                add_x_column(i+1, count_x)
            else:
                if x_PSP_list[i].getAttr('ObjVal') != previous_objval:
                    count_x += 1
                    previous_objval = x_PSP_list[i].getAttr('ObjVal')
                    add_x_column(i+1, count_x)
    '''

    ### 添加新的y_column到SRMP
    def add_y_column(self, y_newcolumnCoeff):  # 参数从optimize函数return得来
        y_newcolumn = Column(y_newcolumnCoeff, self.I_bar_constrs.values())
        self.y[self.Q_bar + 1] = self.SRMP.addVar(obj=1, lb=0, vtype=GRB.CONTINUOUS, column=y_newcolumn,
                                                  name='y' + str(self.Q_bar + 1))
        self.Q_bar += 1
        self.pattern_d[self.Q_bar] = y_newcolumnCoeff

        self.SRMP.update()
        self.SRMP.write('SRMP.lp')
        self.SRMP.optimize()
        return self.SRMP.getVars(), self.SRMP.getAttr('ObjVal')

    ### 更新y_PSP变量和系数,包括求解y_PSP
    def update_y_PSP(self, method):
        if self.I_bar > len(self.d):  # I_bar增加的话，增加vars
            self.y_PSP.remove(self.y_PSP.getConstrs())
            ld = len(self.d)
            for i in range(self.I_bar - ld):
                self.d[i + ld + 1] = self.y_PSP.addVar(lb=0, vtype=GRB.INTEGER, name='d' + str(i + ld + 1))
                self.y_PSP.update()
                self.y_linexpr.addTerms((self.k_i[i + ld] + self.data.B), self.d[i + ld + 1])
                self.y_obj_linexpr.addTerms(self.I_bar_constrs[i + ld].pi, self.d[i + ld + 1])
            self.bay_pattern_constr = self.y_PSP.addLConstr(self.y_linexpr <= self.data.H)
            self.y_PSP.setObjective(1 - self.y_obj_linexpr, GRB.MINIMIZE)
            self.y_PSP.update()

        for i in range(self.I_bar):
            self.d[i + 1].obj = -self.I_bar_constrs[i].pi
        self.y_PSP.update()

        # dp求解y_PSP
        # items
        items = []
        for var in self.y_PSP.getVars():
            items.append([int(self.y_PSP.getCoeff(self.bay_pattern_constr, var)), -var.getAttr('Obj')])
        '''
        print(items)
        '''
        if method == 'dp':
            objval_y_PSP = 1 - CRG.knapsack_dp(items, self.data.H)[0]
            result = CRG.knapsack_dp(items, self.data.H)[2]

        y_coeff = [0] * self.y_PSP.getAttr('NumVars')
        for i in result:
            index = self.k_i.index(i[0] - self.data.B)
            y_coeff[index] = i[1]

        if objval_y_PSP < -0.0001:
            return (True, y_coeff, objval_y_PSP)  # add
        else:
            return (False, 'objval_y_PSP >= 0')  # stop
        '''
        ### 直接求解y_PSP
        y_PSP.optimize()
        if y_PSP.getAttr('ObjVal') < 0:
            return True
        else:
            return False
        '''

    # 如果一个k存在对应多个bin pattern，他们的v_i是相等的
    def add_column_and_row(self, add_k, remain_v, row_y_coeff):
        for k in add_k:
            # add row
            self.I_bar_constrs[self.I_bar] = self.SRMP.addLConstr(0, GRB.GREATER_EQUAL, 0,
                                                                  name='I_constr_bin_height=' + str(k))
            self.I_bar += 1
            self.k_i.append(int(k))
            self.SRMP.update()
        # 加入new x column
        for kk in remain_v:
            if kk[0] in add_k:
                new_row_x_column = Column(-1, self.SRMP.getConstrByName('I_constr_bin_height=' + str(kk[0])))
                for j in range(self.data.J):
                    mj = kk[1][j]
                    nj = self.data.lambda_j_under_k[(kk[0], j)]
                    new_row_x_column.addTerms(mj * nj, self.SRMP.getConstrByName('J_constr_j=' + str(j)))
                self.x[self.P_bar + 1] = self.SRMP.addVar(obj=0, lb=0, column=new_row_x_column, vtype=GRB.CONTINUOUS,
                                                          name='x' + str(self.P_bar + 1))
                self.P_bar += 1
                self.m_pj[self.P_bar] = list(kk[1])

        # y column coeffs
        new_row_y_column = Column()
        for k in row_y_coeff:
            new_row_y_column.addTerms(k[1], self.SRMP.getConstrByName('I_constr_bin_height=' + str(k[0] - self.data.B)))
        # 加入new y column(一个)
        self.y[self.Q_bar - 1] = self.SRMP.addVar(obj=1, lb=0, column=new_row_y_column, vtype=GRB.CONTINUOUS,
                                                  name='y' + str(self.Q_bar))
        self.Q_bar += 1

        self.pattern_d[self.Q_bar] = [0 for ib in range(self.I_bar)]
        for coeff in row_y_coeff:
            if coeff[0] - self.data.B in self.k_i:
                self.pattern_d[self.Q_bar][self.k_i.index(coeff[0] - self.data.B)] = coeff[1]

        self.SRMP.update()
        self.SRMP.write('SRMP.lp')
        self.SRMP.optimize()
        return self.SRMP.getVars(), self.SRMP.getAttr('ObjVal')

    def update_row_PSP(self, method):  # phase 1+2 重新计算要考虑的列，并求解row_PSP
        if self.SRMP.getAttr('NumConstrs') == self.data.J + len(self.data.possible_k):
            return (False, 'No new row to add', None, None)
        ###### phase1 重新计算考虑的列
        alist = []  ### row PSP第二阶段需要考虑的列
        for k in self.data.possible_k:  ### 在每个bin height下找到最好的bin pattern
            items = []
            for j in range(self.data.J):
                value = self.J_constrs[j].pi * self.data.lambda_j_under_k[(k, j)]
                weight = self.data.w_j[j] + self.data.a
                items.append([weight, value])
            capacity = self.data.W - self.data.a
            result = CRG.knapsack_dp(items, float(capacity))  ### 得到每个bin height下的bin patern情况 result

            pattern_to_j = {}
            for j in range(self.data.J):
                pattern_to_j[self.data.w_j[j] + self.data.a] = j

            pattern = [0 for _ in range(self.data.J)]
            for item, num in result[2]:
                related_j = pattern_to_j[item]
                pattern[related_j] = num

            alist.append([k, pattern, result[0]])

        # eliminate bad bin pattern
        for i in range(len(alist)):
            alist[i].insert(0, i)

        del_list = []
        for i in range(len(alist) - 1):
            if alist[i + 1][3] >= alist[i][3]:
                del_list.append(alist[i][1:])

        if len(del_list) > 0:
            for aa in alist:
                if aa[1:] in del_list:
                    del aa
        remain_v = {}
        for aaa in alist:
            key = [aaa[1]]
            key.append(tuple(aaa[2]))
            key = tuple(key)
            remain_v[key] = aaa[3]
        for k in list(remain_v.keys()):
            if k[0] in self.k_i:
                remain_v.pop(k)

        remain_v_only_k = {}
        for k in remain_v.keys():
            remain_v_only_k[k[0]] = remain_v[k]

        # uu = []
        # for j in range(self.data.J):
        #     uu.append(self.J_constrs[j].pi)
        # v = CRG.CalV_i(self, uu, self.data.possible_k)
        #
        # # 除去I_bar中的k_i 已经存在的
        # remain_v = {}
        # for k in v.keys():
        #     if k[0] not in self.k_i:
        #         remain_v[k] = v[k]
        # # print('remain_v ='+str(remain_v))
        # '''
        # remain_v =
        # {(18, (9, 0, 0, 0)): 0.25,
        #  (12, (9, 0, 0, 0)): 0.16666666666666666,
        #  (7, (9, 0, 0, 0)): 0.08333333333333333,
        #  (6, (0, 4, 0, 1)): 0.019772376543209874,
        #  (5, (0, 0, 0, 3)): 0.013020833333333334,
        #  (4, (1, 5, 0, 0)): 0.009645061728395061}
        # '''
        # remain_v_only_k = {}
        # for k in remain_v.keys():
        #     remain_v_only_k[k[0]] = remain_v[k]
        # # print('remain_v_only_k='+str(remain_v_only_k))
        # '''
        # remain_v_only_k =
        # {18: 0.25,
        #  12: 0.16666666666666666,
        #  7: 0.08333333333333333,
        #  6: 0.019772376543209874,
        #  5: 0.013020833333333334,
        #  4: 0.009645061728395061}
        # '''
        ##### phase 1 end

        ##### phase 2 start
        self.row_PSP.remove(self.row_PSP.getVars())
        self.row_PSP.remove(self.row_PSP.getConstrs())
        row_d = {}
        row_d_set = set(self.k_i + list(remain_v_only_k.keys()))
        for k in row_d_set:
            row_d[k] = self.row_PSP.addVar(lb=0, vtype=GRB.INTEGER, name='d_k=' + str(k))
        self.row_PSP.update()

        coeff_for_d = {}  # v值
        for k in row_d.keys():
            if k in self.k_i:
                coeff_for_d[k] = self.I_bar_constrs[self.k_i.index(k)].pi
            else:
                coeff_for_d[k] = remain_v_only_k[k]
        # print(coeff_for_d)
        '''
        coeff_for_d=
        {23: 0.1111111111111111, 
         15: 0.06944444444444445, 
         28: 0.125, 
         14: 0.07142857142857142, 
         18: 0.25, 
         12: 0.16666666666666666, 
         7: 0.08333333333333333, 
         6: 0.019772376543209874, 
         5: 0.013020833333333334, 
         4: 0.009645061728395061}
        '''
        row_bay_pattern_constr = self.row_PSP.addLConstr(
            quicksum(row_d[k] * (k + self.data.B) for k in row_d.keys()) <= self.data.H)
        self.row_PSP.setObjective(
            1 - quicksum(row_d[k] * coeff_for_d[k] for k in row_d.keys()), GRB.MINIMIZE)
        self.row_PSP.update()
        '''
        row_PSP.optimize()
        '''
        ### dp求解row_PSP
        items = []
        for var in self.row_PSP.getVars():
            items.append([int(self.row_PSP.getCoeff(row_bay_pattern_constr, var)), -var.getAttr('Obj')])
        '''
        print(items)
        '''
        if method == 'dp':
            objval_row_PSP = 1 - CRG.knapsack_dp(items, self.data.H)[0]
            result = CRG.knapsack_dp(items, self.data.H)[2]

        add_k = []
        for i in result:
            if (i[0] - self.data.B) not in self.k_i and (i[0] - self.data.B) not in add_k:
                add_k.append(i[0] - self.data.B)

        if objval_row_PSP < 0:
            if len(add_k) > 0:
                return (True, add_k, remain_v, result)
            else:
                return (False, 'add_k is empty', None, None)
        else:
            return (False, 'objval row_PSP >= 0', None, None)

        '''
        # 找到value>0的row_d,并整理全部new k
        add_k = []   #所要加入的k，考虑一个k有多个bin pattern的情况
        for var in row_PSP.getVars():
            if var.getAttr('X') > 0:   
                if (int(var.getAttr('VarName')[4:]) not in k_i) and (int(var.getAttr('VarName')[4:]) not in add_k):
                    add_k.append(int(var.getAttr('VarName')[4:]))

        if (row_PSP.getAttr('ObjVal') < 0):
            if len(add_k)>0:
                return (True, add_k, remain_v)  # add
            else:
                return (False, 'add_k is empty')
        else:
            return (False, 'objval >= 0')  # stop
        '''

    ## 仅有新的k加入时才有row_PSP，否则在y_PSP和x_PSP加入

    def solve_CRG(self):
        # SRMP---------------------
        # Initialize

        # p^bar(4)  q^bar(4)
        start = time.time()

        self.P_bar = self.data.J  # pattern p set
        self.Q_bar = self.data.J  # pattern q set
        self.I_bar = self.data.J

        # pattern q(a bay) contains how many times of bin type i

        # Initialize: one pattern p only have one type i (p=i)  d_pi
        # column number
        for j in range(self.data.J):
            self.m_pj[j + 1] = [0] * self.data.J
            self.m_pj[j + 1][j] = int((self.data.W - self.data.a) / (self.data.a + self.data.w_j[j]))

        n_pj = {}  # SKU quantities in one column
        for j in range(self.data.J):
            n_pj[j + 1] = [0] * self.data.J
            if self.RANDOM == True:
                n_pj[j + 1][j] = int(self.data.f_j[j] / self.data.e_j[j]) + 1
            else:
                n_pj[j + 1][j] = int((self.data.possible_k[0] - self.data.b) / self.data.h_j[j])

        self.k_i = []  # bin height
        for i in range(self.I_bar):
            j = i
            if self.RANDOM == True:
                self.k_i.append(math.ceil(n_pj[i + 1][i] * self.data.h_j[j] + self.data.b))
            else:
                self.k_i.append(self.data.possible_k[0])
        self.k_i = list(set(self.k_i))
        self.I_bar = len(self.k_i)
        self.Q_bar = self.I_bar

        s_ip = tupledict()  # 后面不会再用到了
        for p in range(self.P_bar):
            for i in range(self.I_bar):
                if self.RANDOM == True:
                    if math.ceil(n_pj[p + 1][p] * self.data.h_j[p] + self.data.b) == self.k_i[i]:
                        s_ip[(i, p)] = 1
                    else:
                        s_ip[(i, p)] = 0
                else:
                    s_ip[(i, p)] = 1

        d_qi = tupledict()
        for i in range(self.I_bar):
            for q in range(self.Q_bar):
                if i == q:
                    d_qi[(q, i)] = int(self.data.H / (self.k_i[i] + self.data.B))
                else:
                    d_qi[(q, i)] = 0

        # 记录已经加入的pattern y
        for q in range(self.Q_bar):
            new = []
            for i in range(self.I_bar):
                new.append(d_qi[(q, i)])

            self.pattern_d[q + 1] = new
        # ——————————————————————————————————————————————————————————————
        # decision variables
        self.SRMP = gp.Model('SRMP')
        self.SRMP.setParam('OutputFlag', 0)
        for p in range(self.P_bar):
            self.x[p + 1] = self.SRMP.addVar(lb=0, obj=0, vtype=GRB.CONTINUOUS, name='x' + str(p + 1))
        for q in range(self.Q_bar):
            self.y[q + 1] = self.SRMP.addVar(lb=0, obj=1, vtype=GRB.CONTINUOUS, name='y' + str(q + 1))

        # add constraints
        self.I_bar_constrs = self.SRMP.addConstrs(
            quicksum(d_qi[q, i] * self.y[q + 1] for q in range(self.Q_bar)) >= quicksum(
                s_ip[i, p] * self.x[p + 1] for p in range(self.P_bar)) for i in range(self.I_bar))
        self.J_constrs = self.SRMP.addConstrs(
            quicksum(self.m_pj[p + 1][j] * n_pj[p + 1][j] * self.x[p + 1] for p in range(self.P_bar)) >= self.data.q_j[
                j] for j in range(self.data.J))

        for i in range(self.I_bar):
            self.SRMP.setAttr('ConstrName', self.I_bar_constrs[i], 'I_constr_bin_height=' + str(self.k_i[i]))
        for j in range(self.data.J):
            self.SRMP.setAttr('ConstrName', self.J_constrs[j], 'J_constr_j=' + str(j))

        self.SRMP.setObjective(quicksum(self.y[q + 1] for q in range(self.Q_bar)), GRB.MINIMIZE)

        self.SRMP.update()
        self.SRMP.write('SRMP.lp')
        self.SRMP.optimize()

        # print('SRMP.getVars()------------------------')
        # print(str(self.SRMP.getVars()))
        # print('SRMP.ObjVal------------------------')
        # print(str(self.SRMP.getAttr('ObjVal')))

        # generating---------------------------------------------

        # construct_y_PSP
        self.y_PSP = gp.Model('y_PSP')
        self.y_PSP.setParam('OutputFlag', 0)
        # decision variable
        for i in range(self.I_bar):
            self.d[i + 1] = self.y_PSP.addVar(lb=0, vtype=GRB.INTEGER, name='d' + str(i + 1))

        self.y_linexpr = LinExpr(quicksum(self.d[i + 1] * (self.k_i[i] + self.data.B) for i in range(self.I_bar)))
        self.bay_pattern_constr = self.y_PSP.addLConstr(self.y_linexpr <= self.data.H)

        self.y_obj_linexpr = quicksum(self.I_bar_constrs[i].pi * self.d[i + 1] for i in range(self.I_bar))
        self.y_PSP.setObjective(1 - self.y_obj_linexpr, GRB.MINIMIZE)

        self.y_PSP.update()
        self.y_PSP.write('y_PSP.lp')

        # construct_row_PSP
        self.row_PSP = gp.Model('row_PSP')
        self.row_PSP.setParam('OutputFlag', 0)

        y_iteration = 1
        x_iteration = 1
        row_iteration = 1

        # 是否还需要继续生成flag
        flag_y = True
        flag_x = [False for i in range(self.I_bar)]
        flag_row = False

        while (flag_y or max(flag_x) or flag_row):
            flag_y, y_coeff = CRG.update_y_PSP(self, 'dp')[0], CRG.update_y_PSP(self, 'dp')[1]
            while flag_y:
                if (y_coeff in self.pattern_d.values()):
                    flag_y = False
                else:
                    CRG.add_y_column(self, y_coeff)
                    # print('------------------')
                    # print('y iteration:', y_iteration)
                    # print('------------------')
                    # print('SRMP.Vars:')
                    # print(self.SRMP.getVars())
                    # print('SRMP.objVal:')
                    # print(self.SRMP.getAttr('ObjVal'))
                    flag_y, y_coeff = CRG.update_y_PSP(self, 'dp', )[0], CRG.update_y_PSP(self, 'dp')[1]

                y_iteration += 1

            # print('complete y generating!!!')
            flag_x = []
            pattern = []
            for i0 in range(self.I_bar):
                f, p = CRG.update_x_PSP(self, i0 + 1)[0], CRG.update_x_PSP(self, i0 + 1)[1]
                flag_x.append(f)
                pattern.append(p)

            # print('start x generating!!!')
            while max(flag_x):
                for i in range(self.I_bar):
                    # print('start i=', i + 1)
                    while flag_x[i] == True:
                        this_pattern = pattern[i]
                        # print('pattern get!!!', this_pattern)
                        CRG.add_x_column(self, i + 1, this_pattern)
                        # print('------------------')
                        # print('x iteration:', x_iteration)
                        # print('------------------')
                        # print('SRMP.Vars:')
                        # print(self.SRMP.getVars())
                        # print('SRMP.objVal:')
                        # print(self.SRMP.getAttr('ObjVal'))

                        x_iteration += 1
                        for ii in range(self.I_bar):
                            flag_x[ii], pattern[ii] = CRG.update_x_PSP(self, ii + 1)[0], CRG.update_x_PSP(self, ii + 1)[
                                1]

            # print('complete x generating!!!')
            flag_row, add_k, remain_v, row_y_coeff = CRG.update_row_PSP(self, 'dp')[0], CRG.update_row_PSP(self, 'dp')[
                1], CRG.update_row_PSP(self, 'dp')[2], CRG.update_row_PSP(self, 'dp')[3]
            # print('start row generating!!!')
            while flag_row:
                CRG.add_column_and_row(self, add_k, remain_v, row_y_coeff)
                # print('------------------')
                # print('row iteration:', row_iteration)
                # print('------------------')
                # print('add k=', add_k)
                # print('SRMP.Vars:')
                # print(self.SRMP.getVars())
                # print('SRMP.objVal:')
                # print(self.SRMP.getAttr('ObjVal'))

                row_iteration += 1
                flag_row, add_k, remain_v, row_y_coeff = CRG.update_row_PSP(self, 'dp')[0], \
                                                         CRG.update_row_PSP(self, 'dp')[1], \
                                                         CRG.update_row_PSP(self, 'dp')[2], \
                                                         CRG.update_row_PSP(self, 'dp')[3]
            flag_y, y_coeff = CRG.update_y_PSP(self, 'dp')[0], CRG.update_y_PSP(self, 'dp')[1]

        end1 = time.time()
        # print('CRG runtime: %f s' % (end1 - start))
        # print('total y generation times:', y_iteration - 1)
        # print('total x generation times:', x_iteration - 1)
        # print('total row generation times:', row_iteration - 1)
        # print('---------------------')

        sol = self.SRMP.getAttr('ObjVal')
        print("start IP solution ")
        ## SRMP松弛的解有小数，说明不是MP最优解，否则解都应该为整数    提供上界UB
        ### 设置为integer求解

        for var in self.SRMP.getVars():
            var.setAttr('VType', GRB.INTEGER)
        self.SRMP.update()
        self.SRMP.setParam('OutputFlag', 0)
        self.SRMP.optimize()

        end2 = time.time()

        # print('SRMP.getVars()------------------------')
        # print(str(self.SRMP.getVars()))

        bin_pattern_num, bay_pattern_num = 0, 0
        for var in self.SRMP.getVars():
            if var.X > 0:
                # print(var.VarName, ' : ')
                # print(self.SRMP.getCol(var))
                # print(var.x)
                if 'x' in var.VarName:
                    bin_pattern_num += 1
                else:
                    bay_pattern_num += 1

        # print('SRMP.ObjVal------------------------')
        # print(str(self.SRMP.getAttr('ObjVal')))

        int_sol = self.SRMP.getAttr('ObjVal')

        CRG_time = end1 - start
        IP_time = end2 - end1

        return y_iteration - 1, x_iteration - 1, row_iteration - 1, CRG_time, IP_time, sol, int_sol, bin_pattern_num, bay_pattern_num
