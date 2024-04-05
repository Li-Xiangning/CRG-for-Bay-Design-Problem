\ Model RMP_0_copy_copy_copy_copy_copy_copy_copy_copy_copy_copy_copy_copy_copy
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10
Subject To
 J_constr_j=0: 160 y1 + 12 y7 + 36 y8 + 27 y9 >= 13749
 J_constr_j=1: 72 y2 + 20 y6 + 12 y7 + 2 y8 >= 16223
 J_constr_j=2: 30 y3 + 32 y6 + 18 y7 >= 14047
 J_constr_j=3: 48 y4 + 18 y7 + 36 y8 + 40 y9 + 40 y10 >= 16447
 J_constr_j=4: 40 y5 >= 17517
 branch_right_0: y2 >= 94
 branch_right_2: y7 >= 781
 branch_right_4: y8 >= 42
 branch_right_6: y1 >= 15
 branch_right_8: y9 >= 18
 branch_left_10: y4 <= 3
 branch_right_11: y9 >= 19
 branch_left_13: y4 <= 2
 branch_right_14: y9 >= 20
 branch_left_16: y4 <= 1
 branch_right_17: y9 >= 21
 branch_left_19: y4 <= 0
 branch_left_20: y9 <= 21
Bounds
End
