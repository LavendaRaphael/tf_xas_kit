#!/bin/env python
import xas_module
import os

str_exp=os.environ['goto_pto_exp']
os.chdir(str_exp)
#--------------------------------------------------[extract]
list2d_alpha = []
list2d_alpha.append( [20, '20210924.pto110_a20_postscaling.csv',[1,2]] )
list2d_alpha.append( [41, '20210924.pto110_a41_postscaling.csv',[1,2]] )
#list2d_alpha.append( [ 25, '20211113.Angel-Pt110-OXAS.csv', [0,1] ])
#list2d_alpha.append( [ 30, '20211113.Angel-Pt110-OXAS.csv', [0,2] ])
#list2d_alpha.append( [ 35, '20211113.Angel-Pt110-OXAS.csv', [0,3] ])
#list2d_alpha.append( [ 40, '20211113.Angel-Pt110-OXAS.csv', [0,4] ])
#list2d_alpha.append( [ 45, '20211113.Angel-Pt110-OXAS.csv', [0,5] ])
#list2d_alpha.append( [ 50, '20211113.Angel-Pt110-OXAS.csv', [0,6] ])


xas_module.def_exp_xyzfit( 
    list2d_alpha = list2d_alpha,
    str_outfile = 'xyzfit.csv'
)
