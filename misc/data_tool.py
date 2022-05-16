def def_broad(
        array1d_xdata,
        array2d_ydata,
        float_hwhm,
        str_method,
        int_nbin,
        ):
    float_xl = numpy.amin( array1d_xdata )
    float_xr = numpy.amax( array1d_xdata )
    array1d_xdata_new = numpy.linspace(
        start = float_xl,
        stop = float_xr,
        num = int_nbin
        )

    array2d_delta =  array1d_xdata_new[:,None] - array1d_xdata[None,:]
    array2d_delta = def_lineshape(
        str_method = str_method,
        float_x = array2d_delta,
        float_hwhm = float_hwhm,
        )

    int_xdata_new_shape = array1d_xdata_new.shape[0]
    int_ydata_shape1 = array2d_ydata.shape[1]
    array2d_ydata_new = numpy.empty( shape=( int_xdata_new_shape, int_ydata_shape1 ) )
    for int_z in range( int_ydata_shape1 ):
        array1d_temp = array2d_delta * array2d_ydata[ :,int_z ][None,]
        array2d_ydata_new[ :,int_z ] = numpy.sum( array1d_temp, axis=1 )

    return array1d_xdata_new, array2d_ydata_new

def def_lineshape(
        str_method,
        float_x,
        float_hwhm,
        ):
    if ( str_method == 'gaussian' ):
        float_y = def_gaussian(
            float_x = float_x,
            float_hwhm = float_hwhm,
            )
    elif ( str_method == 'lorentzian' ):
        float_y = def_lorentzian(
            float_x = float_x,
            float_hwhm = float_hwhm,
            )
    else:
        raise ValueError()

    return float_y

def def_gaussian(
        float_x,
        float_hwhm,
        ):
    return numpy.sqrt(numpy.log(2) / numpy.pi) / float_hwhm\
                             * numpy.exp(-(float_x / float_hwhm)**2 * numpy.log(2))
def def_lorentzian(
        float_x,
        float_hwhm,
        ):
    return float_hwhm / numpy.pi / (float_x**2 + float_hwhm**2)

def def_mix(list2d_data):
#------------------------------[]
# list2d_data = []
# list2d_data.append( [ array1d_xdata, array2d_ydata, [0,2], 0.7 ] )
#------------------------------[]
    dict_args = copy.deepcopy(locals())
    for list_temp in dict_args['list2d_data']:
        del list_temp[0:2]
    def_startfunc( dict_args )

    list2d_xydata = []
    for list1d_data in list2d_data:
        list2d_xydata.append( [ list1d_data[0], list1d_data[1] ] )
    array1d_xdata_interp, list1d_ydata_interp = def_interp( list2d_xydata )

    int_len1dxdata = len( array1d_xdata_interp )
    int_lenycolumn = len( list2d_data[0][2] )
    array2d_ydata_mix = numpy.zeros( shape=(int_len1dxdata, int_lenycolumn) )

    int_len2ddata = len(list2d_data)
    for int_i in range(int_len2ddata):
        list1d_ycolumn = list2d_data[int_i][2]
        array2d_ydata = list1d_ydata_interp[ int_i][:, list1d_ycolumn ]
        float_scaling = list2d_data[int_i][3]
        array2d_ydata_mix += array2d_ydata * float_scaling

    def_endfunc()
    return array1d_xdata_interp, array2d_ydata_mix

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
        tuple_mainxrange,
        tuple_postxrange,
        ):
    float_mainscaling = def_findscaling(
        array1d_xdata = array1d_xdata,
        array1d_ydata = array1d_ydata,
        tuple_xrange = tuple_mainxrange
        )

    float_postscaling = def_findscaling(
        array1d_xdata = array1d_xdata,
        array1d_ydata = array1d_ydata,
        tuple_xrange = tuple_postxrange
        )
    dict_scaling = {
        'float_mainscaling' : float_mainscaling,
        'float_postscaling': float_postscaling
        }

    def_print_paras( locals(),['dict_scaling'] )
    return dict_scaling

