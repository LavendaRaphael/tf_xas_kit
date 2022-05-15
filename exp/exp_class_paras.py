def def_class_paras():

    class_paras = class_paras()

    # scaling method
    #class_paras.str_scalingmethod = 'float_mainscaling'
    class_paras.str_scalingmethod = 'float_postscaling'

    return class_paras

class class_paras(object):

    @property
    def str_scalingmethod(self):
        return self._str_scalingmethod
    @str_scalingmethod.setter
    def str_scalingmethod(self, str_temp):
        self._str_scalingmethod = str_temp

