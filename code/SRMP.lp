\ Model SRMP
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10 + y11 + y12 + y13 + y14
   + y15 + y16 + y17 + y18 + y19 + y20 + y21 + y22 + y23 + y24 + y25 + y26
   + y27 + y28 + y29 + y30 + y31 + y32
Subject To
 I_constr_bin_height=37: - x20 + 4 y1 + 4 y21 - x21 - x22 - x46 >= 0
 I_constr_bin_height=8: - x19 + 15 y2 + 2 y25 - x23 - x24 + y31 + y32 >= 0
 I_constr_bin_height=9: - x14 + 13 y3 + 5 y23 + y24 + 12 y25 >= 0
 I_constr_bin_height=11: - x3 - x13 + 12 y4 >= 0
 I_constr_bin_height=12: - x12 + 11 y5 + y21 + 4 y26 + 5 y28 >= 0
 I_constr_bin_height=14: - x7 - x8 + 10 y6 + y20 + y27 - x25 - x26 - x27
   >= 0
 I_constr_bin_height=15: - x2 + 9 y7 + y18 + y23 + 8 y30 - x28 - x29 >= 0
 I_constr_bin_height=16: - x5 - x18 + 9 y8 - x30 - x31 - x32 - x33 - x34
   - x43 >= 0
 I_constr_bin_height=17: - x1 + 8 y9 + y19 + 8 y31 >= 0
 I_constr_bin_height=19: - x6 + 7 y10 + y17 + 7 y18 + y24 - x35 - x36 - x44
   >= 0
 I_constr_bin_height=20: - x17 + 7 y11 + y19 + 4 y22 + 4 y23 + 6 y24 + y29
   - x37 - x38 - x39 - x45 + 7 y32 >= 0
 I_constr_bin_height=21: - x15 + 7 y12 + 3 y27 + 4 y28 - x40 >= 0
 I_constr_bin_height=22: - x11 + 6 y13 + 6 y17 + 6 y29 >= 0
 I_constr_bin_height=23: - x16 + 6 y14 + 5 y19 + 6 y20 >= 0
 I_constr_bin_height=24: - x4 + 6 y15 + 3 y22 + y30 >= 0
 I_constr_bin_height=25: - x9 - x10 + 6 y16 + 4 y26 + 3 y27 - x41 - x42
   >= 0
 J_constr_j=0: 32 x1 + 4 x22 + x24 + 3 x34 + 4 x35 + 4 x39 >= 10786
 J_constr_j=1: 28 x2 + 4 x21 + 2 x23 + 8 x29 >= 11358
 J_constr_j=2: 7 x3 + x25 + x46 >= 18713
 J_constr_j=3: 18 x4 + 4 x31 + 2 x33 + 6 x34 + 2 x36 + 2 x37 + 2 x43
   >= 18289
 J_constr_j=4: 10 x5 + 2 x30 + 8 x31 + 2 x32 + 2 x43 + 2 x45 >= 16895
 J_constr_j=5: 8 x6 + 2 x35 + 2 x36 + 4 x39 + 2 x40 + 2 x44 >= 18170
 J_constr_j=6: 8 x7 + 2 x28 >= 15918
 J_constr_j=7: 6 x8 + 6 x25 + 2 x26 + 4 x27 + 2 x43 >= 17552
 J_constr_j=8: 9 x9 + 2 x36 + 2 x38 + 6 x41 + 3 x42 >= 11073
 J_constr_j=9: 12 x10 + 6 x35 + 3 x36 + 3 x38 + 4 x42 + 3 x44 >= 16193
 J_constr_j=10: 6 x11 >= 14384
 J_constr_j=11: 6 x12 + 2 x26 >= 11741
 J_constr_j=12: 6 x13 >= 18335
 J_constr_j=13: 9 x14 >= 18405
 J_constr_j=14: 6 x15 + 4 x28 + 4 x29 + 6 x40 >= 15590
 J_constr_j=15: 6 x16 + 4 x32 + 4 x33 + 2 x34 + 3 x41 + 3 x42 + 2 x43
   >= 16757
 J_constr_j=16: 6 x17 + 2 x26 + 2 x27 + 6 x37 + 3 x38 + 3 x39 + 6 x45
   >= 14387
 J_constr_j=17: 4 x18 + 4 x30 >= 18005
 J_constr_j=18: 4 x19 + 4 x23 + 4 x24 >= 11612
 J_constr_j=19: 8 x20 + 8 x21 + 8 x22 + 2 x44 + 8 x46 >= 12299
Bounds
End
