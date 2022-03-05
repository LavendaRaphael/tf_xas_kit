import xas_module
import os
import local_module

dict_structure = local_module.def_dict_structure()

str_workdir = 'Pt.110.x2y12_O22_vac/'
str_chdir = dict_structure[ str_workdir ].str_chdir
str_chdir += 'vasp_sch/atom_11/'
os.chdir( str_chdir )

#list1d_angle = [00, 90]
list1d_angle = [90, 45]
xas_module.def_atom_findpeak( 
    list1d_angle = list1d_angle,
    str_workdir = 'Pt.110.x2y12_O22_vac/',
    str_jsonfile = 'xas.a20_b90.json',
    )