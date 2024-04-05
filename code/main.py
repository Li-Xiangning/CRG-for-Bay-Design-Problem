from data_generator import Data as Data
# from nonlinear import NIP as NIP
from CRG import CRG as CRG
from BP import BP as BP
from NIP import NIP
from one_cut_model import *
import math
import time
import numpy as np
import pandas as pd
import random

# try_data_list = {'S1': Data(180,120,1.5,1.5,4), 'S2':Data(240,120,1.5,1.5,4), 'S3':Data(300,100,1.5,1.5,4), 'S4':Data(252,120,1.5,1.5,4)}
try_data_list = {'S1': Data(180,120,1.5,1.5,4)}  # change
# try_J = [5,10,20,30,40,50]
try_J = [8]  # change
instance_num = 1  # change
RANDOM = True # change

for setting, data in try_data_list.items():
    for j in try_J:
        avg_calls_y, avg_calls_x, avg_calls_row = 0, 0, 0
        avg_time_CRG, avg_time_BP, avg_time_NIP = 0, 0, 0
        avg_time_IP = 0
        avg_int_sol = 0
        avg_bin_pat_num, avg_bay_pat_num = 0, 0

        opt_num_CRG = 0
        opt_num_BP = 0
        opt_num_NIP = 0

        avg_gap_CRG = 0
        avg_gap_BP = 0
        avg_gap_NIP = 0

        for inst in range(instance_num):
            if RANDOM == True:
                data.RandomData(j)
                data.CalPossible_k()
            else:
                data.J = j
                # 用于 case study------------------------------
                # data.q_j = [31,144,150,31016,31016,35668,26364,31016,54586,93919,110355,30524]
                # data.w_j = [8.0,13.3,34.1,18.6,43.0,37.4,33.9,26.3,27.6,10.9,13.2,22.4 ]   # 单位mm
                # data.h_j = [4.4,6.8,4.5,2.2,3.0,7.1,3.6,5.1,3.7,2.2,4.8,14.0 ]
                # data.e_j = [417.7,547.0,1784.8,1141.4,1141.4,1236.6,1196.0,744.9,1236.6,1196.0,1196.0,1196.0 ]
                # data.f_j = [100000 for _ in range(data.J)]
                # data.CalPossible_k_for_case(2, 48)

                # 用于 SKU size sensitive analysis------------------------
                # data.q_j = [15178, 10709, 17887, 14485, 13955, 18387, 18387, 11224, 13579, 13579, 14594, 13294, 19307, 16451, 15856, 10832, 15346, 15083, 19211, 11456, 11900, 19161, 18805, 13285, 15978, 16213, 17855, 17880, 18815, 10569]
                # data.w_j = [round(random.uniform(20, 60), 1) for _ in range(j)]
                # data.h_j = [round(random.uniform(4, 12), 1) for _ in range(j)]
                # data.e_j = [39.9, 59.1, 41.3, 56.0, 48.2, 31.1, 31.1, 72.0, 78.0, 78.0, 72.0, 69.7, 70.3, 21.7, 64.6, 53.5, 76.4, 68.2, 31.3, 62.1, 62.5, 72.2, 61.1, 20.3, 67.1, 78.9, 53.1, 31.8, 39.8, 73.5]
                # data.f_j = [56.7, 79.4, 78.4, 87.7, 54.9, 71.3, 71.3, 57.0, 142.6, 142.6, 110.1, 51.1, 124.7, 84.6, 125.4, 70.7, 139.9, 148.5, 96.2, 74.3, 97.2, 101.4, 112.0, 122.4, 75.3, 86.9, 62.3, 51.2, 70.1, 78.8]
                # data.CalPossible_k()
                # --------------------------------------------

                # 用于 reserved space sensitive analysis------------------------

                data.q_j = [17188, 13993, 12283, 13391, 19090]
                data.w_j = [13, 22, 29, 29, 42]
                data.h_j = [2.8, 9.9, 6.8, 6.8, 3.5]
                data.e_j = [51.3, 74.8, 42.8, 42.8, 46.3]
                data.f_j = [107.7, 52.6, 89.0, 89.0, 51.3]

                data.CalPossible_k()
                # --------------------------------------------
            data.CalLambda_j_under_k()
            # data.Findm_j()

            #print('Start Column-and-row generation !!!')
            crg = CRG(data)
            if RANDOM == False:
                crg.RANDOM = False
            y_iteration, x_iteration, row_iteration, CRGtime, IPtime, sol, int_sol, bin_pat_num, bay_pat_num = crg.solve_CRG()

            avg_calls_y += y_iteration
            avg_calls_x += x_iteration
            avg_calls_row += row_iteration
            avg_time_CRG += CRGtime
            avg_time_IP += IPtime
            avg_int_sol += int_sol
            avg_bin_pat_num += bin_pat_num
            avg_bay_pat_num += bay_pat_num
            if abs(sol-int_sol) <= 1:
                opt_num_CRG += 1
            else:
                avg_gap_CRG += math.floor(abs(sol-int_sol))
            print('instance', inst+1)
            print('----------------------')
            print('Column-and-row generation runtime: ', CRGtime)
            print('CRG y iteration time: ', y_iteration)
            print('CRG x iteration time: ', x_iteration)
            print('CRG row iteration time: ', row_iteration)
            print('CRG sol: ', sol)
            print('CRG int sol: ', int_sol)
            print('IP time: ', IPtime)
            print('----------------------')


            print('start one cut model !!!')
            items = SKU_to_items(data)
            onecut_start = time.time()
            R, a = generate_cut_and_plates(items, PlateType(data.W - data.a, data.H, 1))
            print('generation completed, start to solve ')
            thistime = time.time()
            onecut_generate_time = thistime - onecut_start
            print('--------------------one cut model generate time: ', onecut_generate_time)
            model(R, a, items, data)  # 里面print了construct和solve的时间
            # find optimal IP by cut model
            


            # # find optimal solution (check feasible)
            # a = []
            # for aa in range(math.ceil(sol), int(int_sol)):
            #     a.append(aa)
            # for aa in a:
            #     nip = NIP(50, 20, data)
            #     nip.construct_NIP()
            #     flag = nip.check_feasible(aa)
            #     print("check feasible: ", aa, flag)
"""
            bp = BP(data)
            time, lb, ub = bp.solve()

            avg_time_BP += time
            # if abs(ub - lb) < 1:
            #     opt_num_BP += 1
            # else:
            #     avg_gap_BP += math.floor(abs(ub - lb))
            print('BP runtime: ', time)
            print('BP lower bound: ', lb)
            print('BP upper bound: ', ub)
            print('----------------------')

            # # print('Start nonlinear model !!!')
            # nip = NIP(5, 3, data)   # total bin type, total bay type
            # nip.construct_NIP()
            # lb, ub, time = nip.solve_NIP()
            #
            # avg_time_NIP += time
            # if abs(ub - lb) < 1:
            #     opt_num_NIP += 1
            # else:
            #     avg_gap_NIP += abs(ub - lb)
            #
            # print('NIP runtime: ', time)
            # print('NIP lower bound: ', lb)
            # print('NIP upper bound: ', ub)
            # print('----------------------')


        avg_calls_y = avg_calls_y/instance_num
        avg_calls_x = avg_calls_x/instance_num
        avg_calls_row = avg_calls_row/instance_num
        avg_time_CRG = avg_time_CRG/instance_num
        avg_time_IP = avg_time_IP/instance_num
        avg_int_sol = avg_int_sol/instance_num
        avg_bin_pat_num = avg_bin_pat_num/instance_num
        avg_bay_pat_num = avg_bay_pat_num/instance_num
        if opt_num_CRG != instance_num:
            avg_gap_CRG = avg_gap_CRG/(instance_num-opt_num_CRG)
        else:
            avg_gap_CRG  = None

        avg_time_BP = avg_time_BP/instance_num
        # if opt_num_BP != instance_num:
        #     avg_gap_BP = avg_gap_BP/instance_num
        # else:
        #     avg_gap_BP = None

        avg_time_NIP = avg_time_NIP/instance_num
        if opt_num_NIP != instance_num:
            avg_gap_NIP = avg_gap_NIP/instance_num
        else:
            avg_gap_NIP = None
        avg_time_NIP = avg_time_NIP/instance_num

        print(setting, '&', j)
        print('------------CRG summary-------------')
        print('average calls of y: ', avg_calls_y)
        print('average calls of x: ', avg_calls_x)
        print('average calls of row: ', avg_calls_row)
        print('average runtime of CRG: ', avg_time_CRG)
        print('average runtime of IP: ', avg_time_IP)
        print('optimal instance number of CRG: ', opt_num_CRG)
        print('average gap for non optimal instances: ', avg_gap_CRG)
        print('average integer solution: ', avg_int_sol)
        print('average bin pattern num used: ', avg_bin_pat_num)
        print('average bay pattern num used: ', avg_bay_pat_num)

        print('------------BP summary---------------')
        print('average runtime of BP: ', avg_time_BP)
        print('optimal instance number of CRG: ', opt_num_BP)
        # print('average gap of BP : ', avg_gap_BP)

        print('------------NIP summary---------------')
        print('average runtime of NIP: ', avg_time_NIP)
        print('optimal instance number of NIP: ', opt_num_NIP)
        print('average gap of NIP : ', avg_gap_NIP)
"""