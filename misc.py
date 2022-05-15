def def_writedata( list2d_header, list3d_data, str_outfile): 
    def_startfunc( locals(), ['list3d_data'] ) 
 
    with open( str_outfile, 'w', newline='' ) as obj_outfile: 
        obj_outwriter = csv.writer( obj_outfile, delimiter=',' ) 
        # header 
        list1d_header = [] 
        for list1d_temp in list2d_header:  
            list1d_header.extend( list1d_temp ) 
        obj_outwriter.writerow( list1d_header ) 
         
        int_lenline = numpy.shape( list3d_data[0] )[0] 
        def_print_paras( locals(),['int_lenline'] ) 
        for int_i in range(len(list3d_data)): 
            if ( list3d_data[int_i].ndim == 2 ): continue 
            list3d_data[int_i] = list3d_data[int_i].reshape(int_lenline,1) 
        for int_i in range(int_lenline): 
            list1d_data = [] 
            for array2d_temp in list3d_data: 
                list1d_data.extend( array2d_temp[int_i]  )  
            obj_outwriter.writerow( list1d_data ) 
 
    def_endfunc() 
    return 

def def_interp(list2d_data):
#------------------------------[]
# list2d_data = []
# list_data.append( [ array1d_xdata, array2d_ydata ] )
#------------------------------[]
    def_startfunc()

    for list1d_data in list2d_data:
        list1d_data[0] = list1d_data[0].reshape( -1 )

    float_xl = float('-inf')
    float_xr = float('+inf')
    for list1d_data in list2d_data:
        array1d_xdata = list1d_data[0]
        float_xl = max( float_xl, numpy.amin( array1d_xdata ))
        float_xr = min( float_xr, numpy.amax( array1d_xdata ))

    dict_json = {}
    dict_json[ 'float_xl' ] = float_xl
    dict_json[ 'float_xr' ] = float_xr
    print( json.dumps( dict_json, indent=4 ) )

    int_shape1dxdata = numpy.shape(list2d_data[0][0])[0]
    dict_json = {}
    dict_json[ 'int_shape1dxdata' ] = int_shape1dxdata
    print( json.dumps( dict_json, indent=4 ) )

    array1d_xdata_interp = numpy.linspace( start=float_xl, stop=float_xr, num=int_shape1dxdata )

    list1d_ydata_interp = []
    for list1d_data in list2d_data:
        array1d_xdata = list1d_data[0]
        array2d_ydata = list1d_data[1]
        int_shape2dydata_1 = numpy.shape( array2d_ydata )[1]
        array2d_temp = numpy.empty( shape=(int_shape1dxdata, int_shape2dydata_1) )
        for int_i in range(int_shape2dydata_1):
            array2d_temp[:,int_i] = numpy.interp( x=array1d_xdata_interp, xp=array1d_xdata, fp=array2d_ydata[:,int_i] )

        list1d_ydata_interp.append( array2d_temp )

    def_endfunc()
    return array1d_xdata_interp, list1d_ydata_interp

def def_findpeaks(
        array1d_xdata,
        array1d_ydata,
        float_relheight = 0.4,
        float_relprominence = 0.02
        ):
#------------------------------[]
#------------------------------[]
    def_startfunc( locals(), ['array1d_xdata', 'array1d_ydata'] )

    array1d_xdata = numpy.reshape( a=array1d_xdata, newshape=-1 )
    array1d_ydata = numpy.reshape( a=array1d_ydata, newshape=-1 )

    float_y_max = max(array1d_ydata)

    array1d_peak_indices, dict_properties = scipy.signal.find_peaks(
        array1d_ydata,
        height = float_relheight * float_y_max,
        prominence = float_relprominence * float_y_max
        )

    dict_peaks = {}
    dict_peaks[ 'E(eV)' ] = array1d_xdata[ array1d_peak_indices ]
    dict_peaks[ 'relheight' ] = dict_properties[ 'peak_heights' ] / float_y_max
    dict_peaks[ 'relprominence' ] = dict_properties[ 'prominences' ] / float_y_max
    def_print_paras( locals(),['dict_peaks'])

    def_endfunc()
    return dict_peaks

def def_findarea( array1d_xdata, array1d_ydata, tuple_xrange): 
    def_startfunc( locals(), ['array1d_xdata', 'array1d_ydata'] ) 
 
    array1d_xdata = numpy.reshape( a=array1d_xdata, newshape=-1 ) 
    array1d_ydata = numpy.reshape( a=array1d_ydata, newshape=-1 ) 
  
    array_x_new=[] 
    array_y_new=[] 
 
    for int_i in range(len(array1d_xdata)): 
        float_x = array1d_xdata[int_i] 
        float_y = array1d_ydata[int_i] 
        if ( float_x > tuple_xrange[0] and float_x < tuple_xrange[1] ): 
            array_x_new.append( float_x ) 
            array_y_new.append( float_y ) 
 
    float_area = numpy.trapz( y=array_y_new, x=array_x_new ) 
    print(json.dumps( obj={'float_area':float_area}, indent=4)) 
 
    def_endfunc() 
    return float_area

def def_findscaling( 
        array1d_xdata, 
        array1d_ydata, 
        tuple_xrange 
        ): 
 
    if ( numpy.amax( array1d_xdata ) < tuple_xrange[1] ): 
        raise ValueError('tuple_xrange is too large!') 
    float_area = def_findarea( 
        array1d_xdata = array1d_xdata, 
        array1d_ydata = array1d_ydata, 
        tuple_xrange = tuple_xrange 
        ) 
    float_scaling = ( tuple_xrange[1] - tuple_xrange[0] ) / float_area 

    return float_scaling

def def_findscaling_dict(
        array1d_xdata,
        array1d_ydata,
        class_structure
        ):
    float_mainscaling = def_findscaling(
        array1d_xdata = array1d_xdata,
        array1d_ydata = array1d_ydata,
        tuple_xrange = class_structure.tuple_mainxrange
        )

    float_postscaling = def_findscaling(
        array1d_xdata = array1d_xdata,
        array1d_ydata = array1d_ydata,
        tuple_xrange = class_structure.tuple_postxrange
        )
    dict_scaling = {
        'float_mainscaling' : float_mainscaling,
        'float_postscaling': float_postscaling
        }

    def_print_paras( locals(),['dict_scaling'] )
    return dict_scaling

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

