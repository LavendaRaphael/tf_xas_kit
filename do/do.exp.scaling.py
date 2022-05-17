#!/bin/env python
import os
from tf_xas_kit import exp_class_paras
from tf_xas_kit import exp_class_structure

list1d_key = exp_class_structure.def_list1d_key()
dict_structure = exp_class_structure.def_dict_structure()
class_paras = exp_class_paras.class_paras()

for str_key in list1d_key:
    class_structure = dict_structure[ str_key ]
    str_chdir = class_structure.str_chdir
    os.chdir(str_chdir)
    print(os.getcwd())

    xas_module.def_exp_scaling( 
        class_structure = class_structure,
        str_outfile = 'xas.exp_scaling.csv',
        class_paras = class_paras,
        )
'''
#--------------------------------------------------[Pt.110]
xas_module.def_exp_scaling( 
    str_datfile = '20210924.Pt110-XAS.CSV',
    int_xcolumn = 6,
    int_ycolumn = 8,
    str_outfile = '20210924.pto110_a20_postscaling.csv'
    )
xas_module.def_exp_scaling( 
    str_datfile = '20210924.Pt110-XAS.CSV',
    int_xcolumn = 10,
    int_ycolumn = 12,
    str_outfile = '20210924.pto110_a41_postscaling.csv'
    )
'''
