import local_module
import os
import xas_module

list1d_workdir=[]
#----------------------------------
list1d_workdir.append('110.x2y1.a1b2_O3_a1b2')

#----------------------------------[loop]
dict_structure = local_module.def_dict_structure()

for str_workdir in list1d_workdir:
    class_structure = dict_structure[ str_workdir ]
    str_chdir = class_structure.str_chdir
    os.chdir(str_chdir)
    print(os.getcwd())

    xas_module.def_vasp_jobinit( class_structure )
