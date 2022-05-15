import pandas

def def_feff2xas():

    str_xheader = 'E(eV)'
    list1d_yheader = [ 'x','y','z' ]

    df_xas_x = pandas.read_csv(
        'ellip_z/xmu.dat',
        sep=' ',
        skipinitialspace=True,
        header=None,
        names= [str_xheader, list1d_yheader[0]] ,
        usecols=[0,3],
        comment='#',
        )
    print(df_xas_x)

    df_xas_y = df_xas_x.rename( columns={ list1d_yheader[0]: list1d_yheader[1] } )

    df_xas_z = pandas.read_csv(
        'polar_z/xmu.dat',
        sep=' ',
        skipinitialspace=True,
        header=None,
        names=[str_xheader, list1d_yheader[2]],
        usecols=[0,3],
        comment='#',
        )
    print(df_xas_z)

    df_xas = pandas.merge( df_xas_x, df_xas_y, on=str_xheader )
    df_xas = pandas.merge( df_xas, df_xas_z, on=str_xheader )
    print(df_xas)

    # df_xas.to_csv( str_outfile )

    list1d_xheader = [str_xheader]
    array2d_xdata = df_xas[list1d_xheader].to_numpy()
    array2d_ydata = df_xas[list1d_yheader].to_numpy()

    return list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata

