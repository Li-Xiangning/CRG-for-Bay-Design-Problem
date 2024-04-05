import copy
import numpy as np
from gurobipy import *
import gurobipy as gp
from data_generator import Data as data
import time


class Node:
    def __init__(self):
        self.local_LB = 0
        self.local_UB = np.inf
        self.sol = {}
        self.int_sol = {}
        self.branch_var_list = []
        self.RMP = None
        # self.PSP = None
        self.cnt = None
        self.is_integer = False
        self.is_feasible = False
        self.is_pruned = False

    def deepcopy_node(node):
        new_node = Node()
        new_node.local_LB = 0
        new_node.local_UB = np.inf
        new_node.sol = copy.deepcopy(node.sol)
        new_node.int_sol = copy.deepcopy(node.int_sol)
        new_node.branch_var_list = []
        new_node.RMP = node.RMP.copy()
        #new_node.PSP = node.PSP.deepcopy()
        new_node.cnt = node.cnt
        new_node.is_integer = node.is_integer
        new_node.is_feasible = node.is_feasible
        new_node.is_pruned = node.is_pruned
        new_node.local_LB = -np.inf
        new_node.local_UB = np.inf
        return new_node

class BP:
    def __init__(self, data):
        self.data = data

    def init_rmp(self):
        # ----------construct RMP----------
        init_RMP = gp.Model('RMP_0')

        # decision variables
        T_bar = self.data.J
        y = {}
        for t in range(T_bar):
            y[t + 1] = init_RMP.addVar(lb=0, vtype=GRB.CONTINUOUS, obj=1, name='y' + str(t + 1))
        # 记录出现过的列
        A_t = {}
        for t in range(T_bar):
            A_t[t + 1] = [0 for _ in range(self.data.J)]
            term1 = math.floor((self.data.W - self.data.a) / (self.data.w_j[t] + self.data.a))
            term2 = math.floor(self.data.f_j[t] / self.data.e_j[t]) + 1
            term3 = math.floor(self.data.H / math.ceil(self.data.h_j[t] * term2 + self.data.b))
            A_t[t + 1][t] = term1 * term2 * term3

        # add constraints
        J_constrs = init_RMP.addConstrs(
            quicksum(A_t[t + 1][j] * y[t + 1] for t in range(T_bar)) >= self.data.q_j[j] for j in
            range(self.data.J))

        for j in range(self.data.J):
            init_RMP.setAttr('ConstrName', J_constrs[j], 'J_constr_j=' + str(j))

        init_RMP.setObjective(quicksum(y[t + 1] for t in range(T_bar)), GRB.MINIMIZE)

        init_RMP.update()
        init_RMP.write('RMP.lp')
        init_RMP.setParam('OutputFlag', 0)
        init_RMP.optimize()
        return init_RMP

    def model_solve_psp(self, rmp):

        return


    def dp_solve_psp(self, rmp):
        # 用新的possible bin height计算 1 ~ max_bin_height
        possible_k = []
        lambda_j_under_k = {}
        possible_k, lambda_j_under_k = self.data.CalAllLambda_j_under_k()

        psp_possible_bin = {}
        psp_best_u_under_k = {}
        psp_best_v_under_k = {}

        for z_i in possible_k:
            # cur_best = 0
            # self.PSP_possible_bin[z_i + self.data.B] = 0
            # self.PSP_best_u_under_k[z_i] = []
            # self.PSP_best_v_under_k[z_i] = []

            items = []
            for j in range(self.data.J):
                const = rmp.getConstrByName('J_constr_j=' + str(j))
                # print(const)
                value = const.Pi * lambda_j_under_k[(z_i, j)]
                weight = self.data.w_j[j] + self.data.a
                items.append([weight, value])
            capacity = self.data.W - self.data.a
            result = BP.knapsack_dp(items, float(capacity))

            pattern_to_j = {}
            for j in range(self.data.J):
                pattern_to_j[self.data.w_j[j] + self.data.a] = j

            pattern = [0 for _ in range(self.data.J)]
            for item, num in result[2]:
                related_j = pattern_to_j[item]
                pattern[related_j] = num


            psp_possible_bin[z_i + self.data.B] = result[0]
            psp_best_u_under_k[z_i] = pattern
            psp_best_v_under_k[z_i] = [lambda_j_under_k[z_i, j] for j in range(self.data.J)]

            # for u_j in self.data.m_pattern_set:
            #     cur = sum(u_j[j]*self.data.lambda_j_under_k[z_i, j]*self.J_constrs[j].pi for j in range(self.data.J))
            #     if cur > cur_best:
            #         self.PSP_possible_bin[z_i + self.data.B] = cur
            #         cur_best = cur
            #         self.PSP_best_u_under_k[z_i] = u_j
            #         self.PSP_best_v_under_k[z_i] = [self.data.lambda_j_under_k[z_i,j] for j in range(self.data.J)]

        value, weight, bagged = BP.knapsack_dp(list(psp_possible_bin.items()), self.data.H)
        A_j = [0 for _ in range(self.data.J)]
        for j in range(self.data.J):
            A_j[j] = sum(
                psp_best_u_under_k[item[0] - self.data.B][j] * psp_best_v_under_k[item[0] - self.data.B][j] * item[1]
                for
                item in bagged)
        PSP_obj = 1 - value
        return PSP_obj, A_j

    def knapsack_dp(items, C):
        WEIGHT, VALUE = range(2)

        # order by max value per item weight
        items = sorted(items, key=lambda item: item[VALUE] / float(item[WEIGHT]), reverse=True)
        if not isinstance(C, int):
            items_copy = []
            for item in items:
                items_copy.append([int(item[0]*10), item[1]])
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

    def solve(self):
        start = time.time()
        # initial rmp0
        rmp0 = self.init_rmp()
        eps = -0.001
        # node list
        Queue = []

        global_ub = np.inf
        global_lb = 0
        incumbent_node = None

        # root node
        root_node = Node()
        root_node.RMP = rmp0
        root_node.cnt = 0
        Queue.append(root_node)

        while len(Queue) > 0 and global_ub - global_lb > -eps:

            # depth first strategy
            current_node = Queue.pop()
            # column generation solve current_node
            current_node.RMP.optimize()
            while True:
                if current_node.RMP.Status != 2:
                    break
                psp_obj, A_j = self.dp_solve_psp(current_node.RMP)
                if psp_obj > eps:
                    break
                # add column to RMP
                cur_col = A_j
                J_constrs = []  # 加分支后rmp中有多于J个约束，提取出J_constrs
                for j in range(self.data.J):
                    J_constrs.append(current_node.RMP.getConstrByName('J_constr_j=' + str(j)))
                y_newcolumn = Column(cur_col, J_constrs)

                # start (判断新的列是否已经存在)
                cur_cols = []
                for v in current_node.RMP.getVars():
                    cc = []
                    for j in range(self.data.J):
                        cons = current_node.RMP.getConstrByName('J_constr_j=' + str(j))
                        cc.append(current_node.RMP.getCoeff(cons, v))
                    cur_cols.append(cc)
                if cur_col in cur_cols:
                    break
                # end
                cur_var_num = current_node.RMP.NumVars
                current_node.RMP.addVar(obj=1, lb=0, vtype=GRB.CONTINUOUS, column=y_newcolumn,
                           name='y' + str(cur_var_num + 1))
                current_node.RMP.update()
                current_node.RMP.write('RMP.lp')
                current_node.RMP.optimize()
                # print(current_node.RMP.getVars())
                # print(current_node.RMP.getAttr('ObjVal'))

                # update PSP
                psp_obj, A_j = self.dp_solve_psp(current_node.RMP)
                # print(psp_obj)

            # column generation solve end

            # 检查不同情况
            is_integer = True
            is_feasible = True
            is_pruned = False
            if current_node.RMP.Status == 2: # feasible
                is_feasible = True
                for var in current_node.RMP.getVars():
                    current_node.sol[var.VarName] = var.x
                    # 向下取整
                    current_node.int_sol[var.VarName] = (int)(var.x)
                    # check is integer
                    if abs((int)(var.x) - var.x) >= -eps:
                        is_integer = False
                        current_node.branch_var_list.append(var.VarName)

                if is_integer == True: # integer feasible solution
                    # update local LB and UB
                    current_node.local_UB = current_node.RMP.getAttr('ObjVal')
                    current_node.local_LB = current_node.RMP.getAttr('ObjVal')
                    # update global LB and UB
                    if current_node.RMP.getAttr('ObjVal') < global_ub:
                        global_ub = current_node.local_UB
                        incumbent_node = Node.deepcopy_node(current_node)
                        incumbent_node.local_UB = current_node.local_UB
                        incumbent_node.local_LB = current_node.local_LB

                    is_pruned = True  # 最优性剪枝
                else: # 小数解
                    current_node.local_LB = current_node.RMP.getAttr('ObjVal')
                    if current_node.local_LB > global_ub:
                        is_pruned = True  # 界限剪枝
                        is_integer = False
                    else:
                        # ------------- branching ----------------#
                        # find the branching var
                        selected_var = None
                        cur_dist = 1
                        for varname in current_node.branch_var_list:
                            if abs((current_node.sol[varname] - current_node.int_sol[varname]) - 0.5) < cur_dist:
                                cur_dist = abs((current_node.sol[varname] - current_node.int_sol[varname]) - 0.5)
                                selected_var = varname
                        left_bound = current_node.int_sol[selected_var]
                        right_bound = current_node.int_sol[selected_var] + 1
                        # create two child nodes
                        left_node = Node.deepcopy_node(current_node)
                        right_node = Node.deepcopy_node(current_node)
                        # left node
                        branch_var = left_node.RMP.getVarByName(selected_var)
                        left_node.RMP.addConstr(branch_var <= left_bound, name='branch_left_' + str(left_node.cnt))
                        left_node.RMP.update()
                        left_node.cnt += 1
                        # right node
                        branch_var = right_node.RMP.getVarByName(selected_var)
                        right_node.RMP.addConstr(branch_var >= right_bound, name='branch_right_' + str(right_node.cnt))
                        right_node.RMP.update()
                        right_node.cnt += 2

                        Queue.append(left_node)
                        Queue.append(right_node)

                    # # update bound
                    # current_node.local_LB = current_node.RMP.getAttr('ObjVal')
                    # if current_node.RMP.getAttr('ObjVal') > global_lb:
                    #     global_lb = current_node.RMP.getAttr('ObjVal')
            else: # infeasible
                is_feasible = False
                is_pruned = True # 非可行性剪枝

            print('lb: ', global_lb, 'ub: ', global_ub)
            # if global_ub - global_lb < -eps:
            #     break

            # update lower bound
            temp_global_lb = np.inf
            for node in Queue:
                node.RMP.optimize()
                if (node.RMP.status == 2):
                    psp_obj, A = self.dp_solve_psp(node.RMP)
                    if (node.RMP.ObjVal + psp_obj < temp_global_lb):
                        temp_global_lb = node.RMP.ObjVal
            global_lb = temp_global_lb

            # 1800s terminate
            # if time.time() - start >= 1800:
            #     break

            # if abs(global_ub - global_lb) < 1:
            #     break

        end = time.time()
        solve_time = end - start
        return solve_time, global_lb, global_ub




