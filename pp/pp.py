import numpy
from tf_xas_kit import misc
from tf_xas_kit import vasp
from tf_xas_kit import feff

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

def def_code2xas(
        str_code = 'vasp',
        str_outfile = 'xas.csv',
        log_tm2xas = False,
        str_broadmethod = None,
        float_hwhm = None,
        int_broadnbin = None,
        ):
    if ( str_code == 'vasp' ):
        list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata = vasp.def_vasp2xas(
            log_tm2xas = log_tm2xas,
            str_broadmethod = str_broadmethod,
            float_hwhm = float_hwhm,
            int_broadnbin = int_broadnbin,
            )
    elif ( str_code == 'feff' ):
        list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata = feff.def_feff2xas()
    else:
        raise

    misc.def_writedata(
        list2d_header = [list1d_xheader, list1d_yheader],
        list3d_data = [array2d_xdata, array2d_ydata],
        str_outfile = str_outfile)

    return list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata

