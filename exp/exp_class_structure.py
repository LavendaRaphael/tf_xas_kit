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


def def_pto_class(
        marker = None,
        str_workdir='',
        str_datfile = None,
        list1d_column = None,
    ):
    str_exp=os.environ['goto_pto_exp']
    goto_pto_work_110=os.environ['goto_pto_work_110']
    goto_pto_work_111=os.environ['goto_pto_work_111']
    class_structure = class_structure()

    if ('110' in marker):
        tuple_mainxrange = (527.0, 540.0)
        tuple_postxrange = (539, 544)
    elif ( '111' in marker):
        tuple_mainxrange = (527.0, 540.0)
        tuple_postxrange = (534, 544)
    else:
        raise

    if ( 'exp' in marker ):
        goto_pto_work = str_exp
        class_structure.list1d_column = list1d_column
    else:
        raise
    class_structure.str_chdir = goto_pto_work + str_workdir
    class_structure.tuple_mainxrange = tuple_mainxrange
    class_structure.tuple_postxrange = tuple_postxrange
    class_structure.str_datfile = str_datfile

    return class_structure

class class_structure():

    @property
    def list1d_column(self):
        return self._list1d_column
    @list1d_column.setter
    def list1d_column(self, list1d_temp):
        self._list1d_column = list1d_temp

    @property
    def str_datfile(self):
        return self._str_datfile
    @str_datfile.setter
    def str_datfile(self, str_temp):
        self._str_datfile = str_temp

    @property
    def tuple_mainxrange(self):
        return self._tuple_mainxrange
    @tuple_mainxrange.setter
    def tuple_mainxrange(self, tuple_temp):
        self._tuple_mainxrange = tuple_temp

    @property
    def tuple_postxrange(self):
        return self._tuple_postxrange
    @tuple_postxrange.setter
    def tuple_postxrange(self, tuple_temp):
        self._tuple_postxrange = tuple_temp

    @property
    def str_chdir(self):
        return self._str_chdir
    @str_chdir.setter
    def str_chdir(self, str_temp):
        self._str_chdir = str_temp

