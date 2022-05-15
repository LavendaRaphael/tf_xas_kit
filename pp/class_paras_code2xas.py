def def_class_paras():

    class_paras = class_paras()

    # generate xas from MYCARXAS
    #class_paras.log_tm2xas = False
    class_paras.log_tm2xas = True
    class_paras.str_broadmethod = 'gaussian'
    #class_paras.str_broadmethod = 'lorentzian'
    #class_paras.float_hwhm = 0.45*math.sqrt( math.log(2) )
    #class_paras.float_hwhm = 0.225
    class_paras.float_hwhm = 0.4*math.sqrt( 2*math.log(2) )
    class_paras.int_broadnbin = 3000

    str_temp = ''
    class_paras.str_xasfile = 'xas'+str_temp+'.csv'

    return class_paras

class class_paras():

    @property
    def log_tm2xas(self):
        return self._log_tm2xas
    @log_tm2xas.setter
    def log_tm2xas(self, obj_temp):
        self._log_tm2xas = obj_temp

    @property
    def str_broadmethod(self):
        return self._str_broadmethod
    @str_broadmethod.setter
    def str_broadmethod(self, obj_temp):
        self._str_broadmethod = obj_temp

    @property
    def float_hwhm(self):
        return self._float_hwhm
    @float_hwhm.setter
    def float_hwhm(self, obj_temp):
        self._float_hwhm = obj_temp

    @property
    def str_xasfile(self):
        return self._str_xasfile
    @str_xasfile.setter
    def str_xasfile(self, obj_temp):
        self._str_xasfile = obj_temp

    @property
    def int_broadnbin(self):
        return self._int_broadnbin
    @int_broadnbin.setter
    def int_broadnbin(self, obj_temp):
        self._int_broadnbin = obj_temp

