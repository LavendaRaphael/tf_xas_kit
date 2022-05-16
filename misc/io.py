def def_extract( str_datfile, list1d_column, log_head=True, dtype=float ):
#------------------------------[]
#------------------------------[]
    def_startfunc( locals() )

    with open( str_datfile, 'r', newline='' ) as obj_datfile:
        obj_datreader = csv.reader( filter( lambda row: (row.strip() and (row.strip()[0]!='#')), obj_datfile ), delimiter= ' ', skipinitialspace=True )
        if log_head:
            list1d_line = next(obj_datreader)
        list1d_line = next(obj_datreader)
        if (',' in list1d_line[0]):
            delimiter=','
        else:
            delimiter=' '
    def_print_paras( locals(),['delimiter'] )

    list1d_header = []
    list2d_data = []
    with open( str_datfile, 'r', newline='' ) as obj_datfile:
        obj_datreader = csv.reader( filter( lambda row: (row.strip() and (row.strip()[0]!='#')), obj_datfile ), delimiter=delimiter, skipinitialspace=True )
        if log_head:
            list1d_line = next(obj_datreader)
            list1d_header = [ list1d_line[i] for i in list1d_column ]
        else:
            list1d_header = [ '' for i in list1d_column ]
        def_print_paras( locals(),['list1d_header'] )

        for list1d_line in obj_datreader:
            if ( not list1d_line[ list1d_column[0] ] ): continue
            list_temp = []
            for int_i in list1d_column:
                list_temp.append( list1d_line[int_i] )
            list2d_data.append( list_temp )
    array2d_data = numpy.array( list2d_data, dtype=dtype )
    tup_shapedata = numpy.shape(array2d_data)
    def_print_paras( locals(),['tup_shapedata'] )

    def_endfunc()
    return list1d_header, array2d_data

def def_print_paras( dict_localpara, list_paraname ):
    dict_paraprint = {}
    for str_paraname in list_paraname:
        dict_paraprint[ str_paraname ] = dict_localpara[ str_paraname ]
    print(json.dumps( dict_paraprint, indent=4, cls=NumpyEncoder ))

class NumpyEncoder(json.JSONEncoder):
#----------------------------------------------[]
# https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
#----------------------------------------------[]
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        elif isinstance(obj, type):
            return str(obj)
        return json.JSONEncoder.default(self, obj)

def def_startfunc( dict_args={}, list1d_del=[] ):
    this_function_name = inspect.currentframe().f_back.f_code.co_name
    print("#"+'-'*20+"["+this_function_name+"]\n")
    for str_del in list1d_del:
        del dict_args[ str_del ]
    print(json.dumps( obj=dict_args, indent=4, cls=NumpyEncoder ))

def def_endfunc():
    print("#"+'-'*20+"<<\n")

