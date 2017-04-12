from pytriqs.plot.protocol import clip_array
import numpy

def plot_base(self, opt_dict, xlabel, ylabel, X, allow_spectral_mode=False):
    r"""
    Plot protocol for Green's function objects.

    Parameters
    ----------
    opt_dict: dictionary
              Can contain:
              - mode: string, default None
                      Mode to plot the Green's function in:
                      -- 'R': real part only
                      -- 'I': imaginary part only
                      -- 'S': spectral function
              - x_window: tuple, default None
                          (xmin,xmax)
              - name: string, default = ''
                      If not '', it remplaces the name of the function just for this plot.
    xlabel: string
            Label to apply to the x axis.
    ylabel: string
            Label to apply to the y axis.
    X: list
       The x values the object can take, i.e. the mesh.
    allow_spectral_mode: boolean, default False
                         Can the spectral function be measured for this type of Green's function?

    Returns
    -------
    plot_data: list of dict
               Object passed to oplot to plot.
    """
    name = opt_dict.pop('name', self.name)  # consume it
    name_prefix = opt_dict.pop('name_prefix', '')  # consume it
    if name and name_prefix:
        raise ValueError, 'name and name_prefix cannot be used at the same time'
    if name_prefix:
        name_save, self.name = self.name, name or name_prefix

    rx = opt_dict.pop('x_window', None)  # consume it
    X = numpy.array(X).real
    sl = clip_array(X, *rx) if rx else slice(len(X)) # the slice due to clip option x_window

    def mdic(prefix, f):
        from itertools import product
        ind_range = product(*map(range,reversed(self.target_shape)))
        make_label = lambda ind: "%s%s %s" % (prefix,name,"_".join(map(str, reversed(ind))))
        make_data_sl = lambda ind: (sl,) + tuple(reversed(ind))
        return [{'xlabel': xlabel,
                 'ylabel': ylabel(self.name or name),
                 'xdata': X[sl],
                 'label': make_label(ind),
                 'ydata': f(self.data[make_data_sl(ind)])} for ind in ind_range]

    # backward compat.
    if 'RI' in opt_dict:
        assert 'mode' not in opt_dict, "Can not have both flags RI and mode."
        import warnings
        warnings.warn("oplot: 'RI' flag is deprecated, use 'mode' instead")
        opt_dict['mode'] = opt_dict.pop('RI', '')

    # if data is real, overrule
    mode = opt_dict.pop('mode', '')
    if self.data.dtype == numpy.float64 : 
        res = mdic('', lambda x: x)    
    elif mode == '':
        res = mdic('Re ', lambda x: x.real) + mdic('Im ', lambda x: x.imag)
    elif mode == 'R':
        res = mdic('Re ', lambda x: x.real)
    elif mode == 'I':
        res = mdic('Im ', lambda x: x.imag)
    elif mode == 'S':
        if allow_spectral_mode:
            res = mdic('', lambda x: -1 / numpy.pi * x.imag)
        else:
            raise ValueError, "Cannot measure the spectral function for this type of Green's function."
    else:
        raise ValueError, "The 'mode' flag is meaningless. Expected 'R', 'I', or 'S' and I got %s." % mode

    res[0].update(opt_dict) # Add all other unused parameters to the dict
    if name_prefix:
        self.name = name_save
    return res

#------------------
