import os
from XasPtO import local_dict_structure
from XasPtO import local_class_structure
from XasPtO import local_class_paras_pp
from tf_xas_kit import pp

#----------------------------------
list1d_key = local_dict_structure.def_list1d_key()
dict_structure = local_dict_structure.def_dict_structure()

instance_paras = local_class_paras_pp.def_instance_paras_alphabeta()

for str_key in list1d_key:
    dict_input = dict_structure[ str_key ]
    instance_structure = local_class_structure.def_class_structure_alphabeta( dict_input )

    str_chdir = instance_structure.str_chdir
    os.chdir(str_chdir)
    print(os.getcwd())

    pp.def_alphabeta_workflow( 
        instance_structure, 
        instance_paras,
        ) 
    
