import os
import def_class_paras
import def_dict_structure
from tf_xas_kit import pp


#----------------------------------
dict_structure = def_dict_structure.def_dict_structure()
list1d_key = def_dict_structure.def_list1d_key()
instance_paras = def_class_paras.def_class_paras_alphabeta()

for str_key in list1d_key:
    dict_input = dict_structure[ str_key ]
    instance_structure = def_dict_structure.def_class_structure_scaling_json( dict_input )

    str_chdir = instance_structure.str_chdir
    os.chdir(str_chdir)
    print(os.getcwd())

    pp.def_alphabeta_workflow( 
        class_structure = instance_structure, 
        class_paras = instance_paras,
        ) 
    
