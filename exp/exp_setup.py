def def_list1d_key():
    list1d_key=[]

    #---------------------------------- 
    #list1d_key.append('exp.20210926.pto111')
    #list1d_key.append('exp.20210924.pto110_a20')
    #list1d_key.append('exp.20210924.pto110_a41')

    return list1d_key

def def_dict_structure():

    dict_structure = {}
    #===================================================================================
    str_key='exp.20210926.pto111'
    dict_structure[ str_key ] = def_pto_class(
        marker = ['exp','111'],
        str_datfile = '20210926.Pt111-XAS.CSV',
        list1d_column = [0,2]
        )
    str_key='exp.20210924.pto110_a20'
    dict_structure[ str_key ] = def_pto_class(
        marker = ['exp','110'],
        str_datfile = '20210924.Pt110-XAS.CSV',
        list1d_column = [6,8]
        )
    str_key='exp.20210924.pto110_a41'
    dict_structure[ str_key ] = def_pto_class(
        marker = ['exp','110'],
        str_datfile = '20210924.Pt110-XAS.CSV',
        list1d_column = [10,12]
        )

    
