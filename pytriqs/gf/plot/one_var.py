import plot_base

#---------------------------------------------------------------
# A list of plot functions for 
#---------------------------------------------------------------

def imfreq(self, opt_dict):
    r"""
    Plot protocol for GfImFreq objects.

    Parameters
    ----------
    opt_dict: dictionary
              Can contain:
              - mode: string, default None
                      Mode to plot the Green's function in:
                      -- 'R': real part only
                      -- 'I': imaginary part only
              - x_window: tuple, default None 
                          (xmin,xmax)
              - name: string, default = ''
                      If not '', it remplaces the name of the function just for this plot.
    """
    return plot_base.plot_base( self, opt_dict,  r'$\omega_n$',
            lambda name : r'%s$(i\omega_n)$'%name, [x.imag for x in self.mesh] )

#----------------------------------------------------------------

def imtime(self, opt_dict):
    r"""
    Plot protocol for GfImTime objects.

    Parameters
    ----------
    opt_dict: dictionary
              Can contain:
              - mode: string, default None
                      Mode to plot the Green's function in:
                      -- 'R': real part only
                      -- 'I': imaginary part only
              - x_window: tuple, default None 
                          (xmin,xmax)
              - name: string, default = ''
                      If not '', it remplaces the name of the function just for this plot.
    """
    return plot_base.plot_base( self, opt_dict,  r'$\tau$', lambda name : r'%s$(\tau)$'%name,  list(self.mesh) )

#----------------------------------------------------------------

def legendre(self, opt_dict):
    r"""
    Plot protocol for GfLegendre objects.

    Parameters
    ----------
    opt_dict: dictionary
              Can contain:
              - mode: string, default None
                      Mode to plot the Green's function in:
                      -- 'R': real part only
                      -- 'I': imaginary part only
              - x_window: tuple, default None 
                          (xmin,xmax)
              - name: string, default = ''
                      If not '', it remplaces the name of the function just for this plot.
    """
    return plot_base.plot_base( self, opt_dict,  r'$l_n$', lambda name : r'%s$(l_n)$'%name, list(self.mesh) )

#----------------------------------------------------------------

def refreq(self, opt_dict):
    r"""
    Plot protocol for GfReFreq objects.

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
    """
    return plot_base.plot_base(self, opt_dict,  r'$\omega$', lambda name : r'%s$(\omega)$'%name, list(self.mesh), allow_spectral_mode = True)

#----------------------------------------------------------------

def retime (self, opt_dict):
    r"""
    Plot protocol for GfReTime objects.

    Parameters
    ----------
    opt_dict: dictionary
              Can contain:
              - mode: string, default None
                      Mode to plot the Green's function in:
                      -- 'R': real part only
                      -- 'I': imaginary part only
              - x_window: tuple, default None 
                          (xmin,xmax)
              - name: string, default = ''
                      If not '', it remplaces the name of the function just for this plot.
    """
    return plot_base.plot_base(self, opt_dict,  r'$\t$', lambda name : r'%s$(\t)$'%name, list(self.mesh))

