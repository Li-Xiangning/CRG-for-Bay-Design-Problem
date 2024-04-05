from gurobipy import *
import gurobipy as gp
from collections import defaultdict
import time
from data_generator import Data as data

class Item:
    def __init__(self, w, h, d):
        self.width = w
        self.height = h
        self.demand = d
        self.demand_num = 1
        self.sku_class = -1


class PlateType:
    def __init__(self, w, h, stage):
        self.width = w
        self.height = h
        self.stage = stage  # stage: 1, 2, 3(waste)

    def check_same(self, plate):
        if self.width == plate.width and self.height == plate.height and self.stage == plate.stage:
            return True # 相同
        else:
            return False  # 不相同

class Cut:
    def __init__(self, j, p):  # p is a PlateType, j is item type(number)
        self.item_type = j
        self.plate_type = p



def generate_cut_and_plates(items, plate0): # items 为item列表
    # initial cuts set (C) and plates set (R)
    C = [] # queue
    R = [plate0]
    A = {}
    # A[(-1, -1, 0)] = 1

    # initialize some cuts
    item_type = len(items)
    for j0 in range(item_type):
        C.append(Cut(j0, plate0))

    while len(C) > 0:
        cut = C.pop(0)
        # add two plates(top and right)
        ## define stage
        new_top_stage, new_right_stage = 0, 0
        if cut.plate_type.stage == 1:
            new_top_stage = 1
            new_right_stage = 2
        else:
            new_top_stage = 3
            new_right_stage = 2
        ## top plate
        top_plate = None

        if cut.plate_type.height > items[cut.item_type].height:
            if cut.plate_type.stage == 1 or cut.plate_type.stage == 3:
                top_plate = PlateType(cut.plate_type.width, cut.plate_type.height-items[cut.item_type].height, new_top_stage)
            else:
                top_plate = PlateType(items[cut.item_type].width, cut.plate_type.height - items[cut.item_type].height,
                                      new_top_stage)

            # check if is already in
            top_flag = False
            for pla in R:
                if top_plate.check_same(pla):
                    top_flag = True
                    break

            # is not in R, add
            if top_flag == False and top_plate.stage < 3:
                R.append(top_plate)

            if top_plate.stage < 3:
                if top_flag == False:
                    A[(cut.item_type, cut.plate_type, top_plate)] = 1
                else:
                    A[(cut.item_type, cut.plate_type, pla)] = 1

        ## right plate
        right_plate = None
        if cut.plate_type.width > items[cut.item_type].width:
            if cut.plate_type.stage == 1 or cut.plate_type.stage == 3:
                right_plate = PlateType(cut.plate_type.width-items[cut.item_type].width, items[cut.item_type].height, new_right_stage)
            else:
                right_plate = PlateType(cut.plate_type.width - items[cut.item_type].width, cut.plate_type.height,
                                        new_right_stage)

            # check if already in
            right_flag = False
            for pla in R:
                if right_plate.check_same(pla):
                    right_flag = True
                    break

            # is not in R, add
            if right_flag == False and right_plate.stage < 3:
                R.append(right_plate)

            if right_plate.stage < 3:
                if right_flag == False:
                    A[(cut.item_type, cut.plate_type, right_plate)] = 1
                else:
                    A[(cut.item_type, cut.plate_type, pla)] = 1


        # for all item type
        for j in range(item_type):
            item = items[j]
            # check can be cut from the top plate ?
            if (top_plate != None) and (top_flag == False) and (item.width <= top_plate.width) and (item.height <= top_plate.height) and top_plate.stage < 3:
                C.append(Cut(j, top_plate))
            # check can be cut from the right plate ?
            if (right_plate != None) and (right_flag == False)and (item.width <= right_plate.width) and (item.height <= right_plate.height) and right_plate.stage < 3:
                C.append(Cut(j, right_plate))

    for plate in R:
        flag = False
        for item in items:
            if plate.width >= item.width and plate.height >= item.height:
                flag = True
                break
        if flag == False:
            R.remove(plate)
    return R, A


def model(R, A, items, data):
    time1 = time.time()
    # item type
    I = len(items)

    # change items demand to sku class demand
    Class = data.J
    class_demand = []
    for j in range(data.J):
        class_demand.append(data.q_j[j])

    # plate type
    J = len(R)
    # set parameters a
    a = defaultdict(int)  # a[(i,j,k)] = 0/1
    # for i, pla_j, pla_k in A.keys():
    #     a[(i, R.index(pla_j), R.index(pla_k))] = 1
    for i in range(I):
        for j in range(J):
            for k in range(J):
                if (i, R[j], R[k]) in A.keys():
                    a[(i, j, k)] = 1

    # for k in a.keys():
    #     if a[k] == 1:
    #         print(k)
    # set params: if item i can be cut from plate j
    b = {}

    model = gp.Model('one-cut model')
    # decision variables
    x = {}
    for i in range(I):
        for j in range(J):
            x[(i, j)] = model.addVar(lb=0, vtype=GRB.INTEGER, name='x_' + str(i) + '_' + str(j))
            # set param if item i can be cut from plate j
            if (R[j].width >= items[i].width) and (R[j].height >= items[i].height):
                b[(i, j)] = 1
            else:
                b[(i, j)] = 0


    # constraints
    # model.addConstrs((quicksum(b[(i, j)] * x[(i, j)] for j in range(J)) >= items[i].demand for i in range(I)), name='item demand')

    # -> class demand constraint
    model.addConstrs((quicksum(b[(i, j)] * items[i].demand_num * x[(i, j)] for j in range(J) for i in range(I) if items[i].sku_class == sku) >= class_demand[sku] for sku in range(Class)), name='demand')

    model.addConstrs((quicksum(a[(i, j, k)] * x[(i, j)] for j in range(J) for i in range(I)) >= quicksum(b[(i, k)] * x[(i, k)] for i in range(I)) for k in range(1, J)), name='plate_num')

    # objective
    model.setObjective(quicksum(x[(i, 0)] for i in range(I)), GRB.MINIMIZE)

    # solve
    model.Params.TimeLimit = 1800
    model.update()
    print('model construct completed')
    time2 = time.time()
    construct_time = time2 - time1
    # model.write('one-cut model.lp')
    model.optimize()
    time3 = time.time()
    solve_time = time3 - time2
    print('one cut model construct time: ', construct_time)
    print('one cut model solve time: ', solve_time)
    print('construct and solve time: ', time3 - time1)

    # for var in model.getVars():
    #     print(var.getAttr('VarName'), var.getAttr('X'))
    return


def SKU_to_items(data):
    items = []

    for j in range(data.J):
        max_stack_num = int(data.f_j[j] / data.e_j[j]) + 1
        for stack_num in range(1, max_stack_num + 1):
            w = data.w_j[j] + data.a
            h = math.ceil(data.h_j[j] * stack_num + data.b + data.B)
            new_item = Item(w, h, data.q_j[j])
            new_item.sku_class = j
            new_item.demand_num = stack_num
            items.append(new_item)

    return items


# items = [Item(4,3,5), Item(2,2,5)]
# plate0 = PlateType(6,6,1)
# R, A = generate_cut_and_plates(items, plate0)
#
# for k,v in A.items():
#     print(k[0],' ', k[1].width, k[1].height, k[1].stage, '  ',k[2].width, k[2].height, k[2].stage)
#
# model(R, A, items)


def check_feasible(R, A, items, data, a):
    time1 = time.time()
    # item type
    I = len(items)

    # change items demand to sku class demand
    Class = data.J
    class_demand = []
    for j in range(data.J):
        class_demand.append(data.q_j[j])

    # plate type
    J = len(R)
    # set parameters a
    a = defaultdict(int)  # a[(i,j,k)] = 0/1
    # for i, pla_j, pla_k in A.keys():
    #     a[(i, R.index(pla_j), R.index(pla_k))] = 1
    for i in range(I):
        for j in range(J):
            for k in range(J):
                if (i, R[j], R[k]) in A.keys():
                    a[(i, j, k)] = 1

    # for k in a.keys():
    #     if a[k] == 1:
    #         print(k)
    # set params: if item i can be cut from plate j
    b = {}

    model = gp.Model('one-cut model')
    # decision variables
    x = {}
    for i in range(I):
        for j in range(J):
            x[(i, j)] = model.addVar(lb=0, vtype=GRB.INTEGER, name='x_' + str(i) + '_' + str(j))
            # set param if item i can be cut from plate j
            if (R[j].width >= items[i].width) and (R[j].height >= items[i].height):
                b[(i, j)] = 1
            else:
                b[(i, j)] = 0


    # constraints
    # model.addConstrs((quicksum(b[(i, j)] * x[(i, j)] for j in range(J)) >= items[i].demand for i in range(I)), name='item demand')

    # -> class demand constraint
    model.addConstrs((quicksum(b[(i, j)] * items[i].demand_num * x[(i, j)] for j in range(J) for i in range(I) if items[i].sku_class == sku) >= class_demand[sku] for sku in range(Class)), name='demand')

    model.addConstrs((quicksum(a[(i, j, k)] * x[(i, j)] for j in range(J) for i in range(I)) >= quicksum(b[(i, k)] * x[(i, k)] for i in range(I)) for k in range(1, J)), name='plate_num')

    model.addConstr(quicksum(x[(i, 0)] for i in range(I)) <= a, name='check_feasible')
    # objective
    model.setObjective(0)

    # solve
    model.Params.TimeLimit = 1800
    model.update()
    model.write('one_cut_model.lp')
    model.optimize()
    return


