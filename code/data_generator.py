import math
import random
####目前按照k_i integer

'''
H = 260  # bay height
W = 120  # bay width
a = 1.5  # horizontal space
b = 1.5  # vertical space
B = 4  # beam height

J = 50

q_j = [20000, 10000, 15000, 20000, 15000]  # product number
w_j = [10.9, 18.6, 27.6, 34.1, 45]  #product width
h_j = [4.3, 2.2, 3.1, 5.2, 8]  # product height
e_j = [20, 20, 30, 30, 70]  #product weight
f_j = [100, 100, 100, 100, 300]   #maximum allowable weight
'''

class Data:
    def __init__(self,H,W,a,b,B):
        self.H = H
        self.W = W
        self.a = a
        self.b = b
        self.B = B
        self.q_j = []
        self.w_j = []
        self.h_j = []
        self.e_j = []
        self.f_j = []
        self.J = 0
        self.possible_k = []
        self.m_pattern_set = []
        self.lambda_j_under_k = {}

    def RandomData(self, J):
        self.J = J
        self.q_j = [random.randint(10000, 20000) for _ in range(J)]
        self.w_j = [round(random.uniform(10, 50), 1) for _ in range(J)]
        self.h_j = [round(random.uniform(2, 10), 1) for _ in range(J)]
        self.e_j = [round(random.uniform(20, 80), 1) for _ in range(J)]
        self.f_j = [round(random.uniform(50, 150), 1) for _ in range(J)]

        # 按宽度从小到大
        w_dict = {}
        for j in range(J):
            w_dict[self.w_j[j]] = j

        self.w_j.sort()
        q_j1, h_j1, e_j1, f_j1 = [],[],[],[]
        for j in range(J):
            q_j1.append(self.q_j[w_dict[self.w_j[j]]])
            h_j1.append(self.h_j[w_dict[self.w_j[j]]])
            e_j1.append(self.e_j[w_dict[self.w_j[j]]])
            f_j1.append(self.f_j[w_dict[self.w_j[j]]])

        self.q_j, self.h_j, self.e_j, self.f_j = q_j1, h_j1, e_j1, f_j1

        print('q_j------------')
        print(self.q_j)
        print('w_j------------')
        print(self.w_j)
        print('h_j------------')
        print(self.h_j)
        print('e_j------------')
        print(self.e_j)
        print('f_j------------')
        print(self.f_j)

        return

    def CalMaxBinHeight(self):
        max_bin_h = 0
        for j in range(self.J):
            max_stack_num = int(self.f_j[j] / self.e_j[j]) + 1
            cur_bin_h = math.ceil(self.h_j[j] * max_stack_num + self.b)
            max_bin_h = max(max_bin_h, cur_bin_h)
        return max_bin_h

    def CalAllLambda_j_under_k(self):
        #n_j / lambdd_j--------------------------
        all_lambda_j_under_k={}
        possible_k = []
        max_bin_height = self.CalMaxBinHeight()
        for aa in range(1,max_bin_height+1):
            possible_k.append(aa)

        for k in possible_k:
            for j in range(self.J):
                #maximum quantity in a column, n_j
                all_lambda_j_under_k[(k,j)] = min(int((k-self.b)/self.h_j[j]), int(self.f_j[j]/self.e_j[j])+1)
        #print('lambda_j_under_k-----------')
        #print(lambda_j_under_k)
        return possible_k, all_lambda_j_under_k


    def CalPossible_k(self):
        # possible bin height values k_i
        possible_k=[]
        for j in range(self.J):
            for n in range(1,int(self.f_j[j]/self.e_j[j])+1+1):
                if math.ceil(n * self.h_j[j]+self.b) not in possible_k:
                    possible_k.append(math.ceil(n * self.h_j[j]+self.b))
        possible_k.sort()
        possible_k.reverse()

        self.possible_k = possible_k
        print('possible_k------------')
        print(possible_k)

        return
        '''
        possible_k =     全部/15   不变
        [28, 23, 19, 18, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 4]   
        '''

    def CalPossible_k_for_case(self, step, max_height):  # step, max_height
        # possible_k=[]
        # for j in range(self.J):
        #     for n in range(1,int(self.f_j[j]/self.e_j[j])+1+1):
        #         if math.ceil(n * self.h_j[j]+self.b) not in possible_k:
        #             possible_k.append(math.ceil(n * self.h_j[j]+self.b))
        # possible_k.sort()
        # possible_k.reverse()

        step_k = []
        for h in range(step, max_height+step, step):
            step_k.append(h)
        step_k.reverse()

        self.possible_k = step_k
        return


    def CalLambda_j_under_k(self):

        #n_j / lambdd_j--------------------------
        lambda_j_under_k={}
        for k in self.possible_k:
            for j in range(self.J):
                #maximum quantity in a column, n_j
                lambda_j_under_k[(k,j)] = min(int((k-self.b)/self.h_j[j]), int(self.f_j[j]/self.e_j[j])+1)
        print('lambda_j_under_k-----------')
        print(lambda_j_under_k)
        self.lambda_j_under_k = lambda_j_under_k
        return self.possible_k
        '''
        lambda_j_under_k = (k,j):lambda_j    不变
        {(4, 0): 0, (4, 1): 1, (4, 2): 0, (4, 3): 0, 
         (5, 0): 0, (5, 1): 1, (5, 2): 0, (5, 3): 1, 
         (6, 0): 0, (6, 1): 2, (6, 2): 1, (6, 3): 1, 
         (7, 0): 1, (7, 1): 2, (7, 2): 1, (7, 3): 1, 
         (8, 0): 1, (8, 1): 2, (8, 2): 1, (8, 3): 2, 
         (9, 0): 1, (9, 1): 3, (9, 2): 1, (9, 3): 2, 
         (11, 0): 1, (11, 1): 4, (11, 2): 2, (11, 3): 3, 
         (12, 0): 2, (12, 1): 4, (12, 2): 2, (12, 3): 3, 
         (13, 0): 2, (13, 1): 5, (13, 2): 2, (13, 3): 3, 
         (14, 0): 2, (14, 1): 5, (14, 2): 2, (14, 3): 4, 
         (15, 0): 2, (15, 1): 6, (15, 2): 3, (15, 3): 4, 
         (18, 0): 3, (18, 1): 6, (18, 2): 3, (18, 3): 4, 
         (19, 0): 3, (19, 1): 6, (19, 2): 4, (19, 3): 4, 
         (23, 0): 4, (23, 1): 6, (23, 2): 5, (23, 3): 4, 
         (28, 0): 4, (28, 1): 6, (28, 2): 6, (28, 3): 4}
        '''

    def Findm_j(self):
        #————————————————————————————————————————
        #寻找横向组合m_j

        w_j_plus_a = []
        for i in range(len(self.w_j)):
            w_j_plus_a.append(self.w_j[i])
            w_j_plus_a[i] += self.a

        def find(start, sum, candidates):
            if sum < min(candidates) and sum >= 0:
                result.append(path[:])

            elif sum <= 0:
                return
            else:
                for j in range(start, len(candidates)):
                    path.append(j)
                    find(j, sum - candidates[j], candidates)
                    path.pop()
                    if candidates[j] > sum:
                        break
            return result
        path=[]
        result=[]
        m_pattern_set = find(0, self.W-self.a, w_j_plus_a)
        #print('m_pattern_set-------------')
        #print(m_pattern_set)
        print('len(m_pattern_set)------')
        print(len(m_pattern_set))
        '''  
        m_pattern_set =       33
        [[0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1], 
         [0, 0, 0, 0, 0, 0, 0, 2], [0, 0, 0, 0, 0, 0, 1, 1], 
         [0, 0, 0, 0, 0, 0, 3], [0, 0, 0, 0, 0, 1, 2], 
         [0, 0, 0, 0, 0, 1, 3], [0, 0, 0, 0, 1, 1, 1], 
         [0, 0, 0, 0, 2, 2], [0, 0, 0, 0, 2, 3], 
         [0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 1, 1, 2], 
         [0, 0, 0, 1, 1, 3], [0, 0, 0, 1, 2, 2], 
         [0, 0, 0, 3, 3], [0, 0, 1, 1, 1, 2], 
         [0, 0, 1, 2, 3], [0, 0, 1, 3, 3], 
         [0, 0, 2, 2, 2], [0, 1, 1, 1, 1, 1], 
         [0, 1, 1, 1, 3], [0, 1, 1, 2, 2], 
         [0, 1, 1, 2, 3], [0, 2, 2, 3], 
         [0, 2, 3, 3], [1, 1, 1, 1, 2], 
         [1, 1, 1, 1, 3], [1, 1, 1, 2, 2], 
         [1, 1, 3, 3], [1, 2, 2, 2], 
         [1, 2, 2, 3], [2, 2, 2, 2], [3, 3, 3]]
        '''
        # 整理m_patttern
        m_pattern_set_count = []
        for i in range(len(m_pattern_set)):
            count_j_column = [0]*self.J #计数 每个m_pattern中每种SKU的列数 即m_j
            for m in m_pattern_set[i]:
                for j in range(self.J):
                    if m == j:
                        count_j_column[j] += 1
            m_pattern_set_count.append(count_j_column)
        #print('m_pattern_set_count--------------')
        #print(m_pattern_set_count)
        self.m_pattern_set = m_pattern_set_count
        return
        '''
        m_pattern_set_count = [m_j]      33    不变
        [[9, 0, 0, 0], [7, 1, 0, 0], [7, 0, 1, 0], [6, 2, 0, 0], [6, 0, 0, 1], 
        [5, 1, 1, 0], [5, 1, 0, 1], [4, 3, 0, 0], [4, 0, 2, 0], [4, 0, 1, 1], 
        [3, 4, 0, 0], [3, 2, 1, 0], [3, 2, 0, 1], [3, 1, 2, 0], [3, 0, 0, 2], 
        [2, 3, 1, 0], [2, 1, 1, 1], [2, 1, 0, 2], [2, 0, 3, 0], [1, 5, 0, 0], 
        [1, 3, 0, 1], [1, 2, 2, 0], [1, 2, 1, 1], [1, 0, 2, 1], [1, 0, 1, 2], 
        [0, 4, 1, 0], [0, 4, 0, 1], [0, 3, 2, 0], [0, 2, 0, 2], [0, 1, 3, 0], 
        [0, 1, 2, 1], [0, 0, 4, 0], [0, 0, 0, 3]]
        '''



