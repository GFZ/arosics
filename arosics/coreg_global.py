# -*- coding: utf-8 -*-
__author__='Daniel Scheffler'

from .coreg_baseclass import CoReg

import inspect

def inheritdocstring(name, bases, attrs):
    if not '__doc__' in attrs:
        # create a temporary 'parent' to (greatly) simplify the MRO search
        temp = type('temporaryclass', bases, {})
        for cls in inspect.getmro(temp):
            if cls.__doc__ is not None:
                attrs['__doc__'] = cls.__doc__
                break

    return type(name, bases, attrs)



class COREG_GLOBAL(CoReg):
    #__metaclass__ = inheritdocstring

    def __init__(self, im_ref, *args, win_pos=(None,None), **kwargs):
        """
        :param win_pos(tuple):          custom matching window position as map values in the same projection like the
                                        reference image (default: central position of image overlap)
        """
        super(COREG_GLOBAL, self).__init__(im_ref, *args, **kwargs)
        self.__doc__ += CoReg.__doc__
        self.win_pos_XY = win_pos  # updated by self.get_opt_winpos_winsize()



