from typing import List, Tuple
from astropy.io import fits
import numpy as np

DEFAULT_SCALE_RADIUS, DEFAULT_SCALE_WL = 4787, 18


class Spectrum( dict ):
    __z = float( )
    __gmag = float( )
    __namestring = str( )

    def __init__( self, **kwargs ):
        """
        Spectrum constructor.  Possible kwargs:
            z: redshift value
        
            gmag: Fiber magnitude in g
        
            namestring or ns: Spectrum ID namestring
        
            dict: A premade { wl: (flux, err) } dictionary
        
        All other values will raise a KeyError
        
        Contructor CAN be called with no values passed.  Values are assigned types, but not initialized.
        Getters and setters are available for z, gmag and redshift.  Use the dictionary access form for the
        spectrographic data.
        
        :param kwargs:
        :type kwargs: dict
        :raises: KeyError
        """
        super( Spectrum, self ).__init__( )
        if len( kwargs ) != 0:
            for arg, val in kwargs.items( ):
                if arg in 'z':
                    self.__z = val
                elif arg in 'gmag':
                    self.__gmag = val
                elif arg in "namestring" or arg in "ns":
                    self.__namestring = val
                elif arg in "dict":
                    self.update( val )
                else:
                    raise KeyError(
                        "Unknown key value in Spectrum constructor( **kwargs )\nKey: %s\nValue: %s" % (arg, val) )
        return

    def __repr__( self ):
        return '%s  z: %s   gmag: %s\n%s    %s' % (
            self.getNS( ), self.getRS( ), self.getGmag( ), self.getWavelengths( )[ 0 ], self.getWavelengths( )[ -1 ])

    def abErr( self, wl_range: Tuple[ float, float ] = (None, None) ) -> float:
        """
        Determines the AB magnitude at every point within the wl_range individually.
        Returns the standard deviation of that value.

        :param wl_range: ( low, high ), defaults to DEFAULT_SCALE_WL +/- DEFAULT_SCALE_RADIUS
        :type wl_range: tuple
        :return: AB Magnitude error
        :rtype: float
        """
        from numpy import log10, nanstd
        minwl = wl_range[ 0 ] or DEFAULT_SCALE_WL - DEFAULT_SCALE_RADIUS
        maxwl = wl_range[ 1 ] or DEFAULT_SCALE_WL + DEFAULT_SCALE_RADIUS
        err_v = list( )
        for wl in self.getWavelengths( ):
            if minwl <= wl <= maxwl:
                f_v = 3.34E4 * pow( wl, 2 ) * 1E-17 * self[ wl ][ 0 ]
                if f_v < 0:
                    continue
                err_v.append( -2.5 * log10( f_v ) + 8.9 )
        return float( nanstd( err_v ) )

    @DeprecationWarning
    def align( self, wlList: List ):
        wls = self.getWavelengths()
        for wl in wls:
            if wl not in wlList:
                del self[ wl ]

    @DeprecationWarning
    def alignToSpec(self, spec ):
        self_wls = set( self.getWavelengths() )
        spec_wls = set( spec.getWavelengths() )

        for wl in ( self_wls - spec_wls ):
            del self[ wl ]
        for wl in ( spec_wls - self_wls ):
            del spec[ wl ]

    def aveFlux( self, central_wl=DEFAULT_SCALE_WL, radius=DEFAULT_SCALE_RADIUS ) -> float:
        """
        Determines the average flux density within a radius of a central wavelength.
        
        :param central_wl:
        :type central_wl: float
        :param radius: 
        :type radius: float
        :return: 
        :rtype: float
        """
        central_wl = central_wl or DEFAULT_SCALE_WL
        radius = radius or DEFAULT_SCALE_RADIUS
        s = 0
        n = 0
        for wl in self.getWavelengths():
            if(central_wl - radius <= wl <= central_wl + radius):
                s += self.getFlux( wl )
                n += 1
        try:
            return s / n
        except ZeroDivisionError as e:
            print(
                f"Spectrum.aveFlux: ZeroDivisionError - unable to determine average flux for spectrum {self.getNS() }.  Is the region of interest loaded?" )
            print( f"central_wl: {central_wl}     radius: {radius}" )
            print( self, flush=True )
            from sys import exit
            exit( 1 )

    def bin(self, size: float = None, init_wl: float = None, inplace=False):
        """
        Bin the spectrum wavelengths and conserving the flux and error.
        :param size: The size of the bins delta lumbda.
        :param init_wl: Where to start the binning
        :param inplace:
        :return:
        """

        z = self.getRS()
        size = size or (z + 1) * 1  # 1 is step size at rest

        wavelengths_source = np.array(self.getWavelengths())
        flux_source = np.array(self.getFluxlist())
        error_source = np.array(self.getErrList())

        init_wl = init_wl or int(wavelengths_source[0]/size) * size
        binned_wls = np.arange(init_wl, int(wavelengths_source[-1]), size)

        binned_flx = []
        binned_err = []

        binflux_append = binned_flx.append
        binerr_append =binned_err.append

        for wl_bin in binned_wls:

            lst = [i for i in range(len(wavelengths_source))  # list of indices of wavelengths in a given bin
                   if wl_bin <= wavelengths_source[i] < wl_bin + size]

            if len(lst) == 0:
                index_of_prev_wl = [i for i in range(len(wavelengths_source))
                                    if wavelengths_source[i] <= wl_bin][-1]
                binflux_append(flux_source[index_of_prev_wl])
                binerr_append(error_source[index_of_prev_wl])
                continue

            last_wl_in_bin = wavelengths_source[lst[-1]]

            # special cases for flux
            first_flux_frac = flux_source[lst[0] - 1] * (wavelengths_source[lst[0]] - wl_bin) / size
            last_flux_frac = flux_source[lst[-1]] * (wl_bin + size - last_wl_in_bin) / size
            total_flux = last_flux_frac + first_flux_frac

            # error calculations
            first_flux_error = error_source[lst[0] - 1] * (wavelengths_source[lst[0]] - wl_bin) / size
            last_flux_error = error_source[lst[-1]] * (wl_bin + size - last_wl_in_bin) / size
            total_error = first_flux_error + last_flux_error

            # if len(lst) == 1:
            #     binned_flx.append(total_flux)
            #     continue

            for i in lst[:-1]:
                total_flux += flux_source[i] * (wavelengths_source[i + 1] - wavelengths_source[i]) / size
                total_error += error_source[i]
            if wl_bin == binned_wls[0]:
                indices = [i for i in range(3) if wavelengths_source[i] < init_wl]
                total_flux += flux_source[indices].sum() / size
                total_error += error_source[indices].sum()

            binflux_append(total_flux)
            binerr_append(total_error)
        if not inplace:
            spec2 = self.cpy_info()
            spec2.setDict(binned_wls, binned_flx, binned_err)
            return spec2
        else:
            self.setDict(wavelengthList=binned_wls,fluxList=binned_flx,errList=binned_err)
            return None

    def chi(self, main, wls=None, scaled=True, central_wl=4000, radius=18, normalized=False):
        """
        Chi analysis of the flux of self against main.
        :param main:
        :param wls: Analysis on some subset of wavelengths.
        :param scaled: Scale, by multiplying by the right coefficient,  the two spectra in some region
        :param central_wl:
        :param radius:
        :param normalized: With respect to the number of data points used.
        :return:
        """
        if scaled:
            spec = self.scale2(to=main, central_wl=central_wl, radius=radius,inplace=False)
        else:
            spec = self
        if wls is None:
            wls = main.getWavelengths()
        specwls = spec.getWavelengths()
        wl0, wlf = max(specwls[0], wls[0]), min(wls[-1], specwls[-1])
        chi = sum([(spec[wl][0] - main[wl][0]) ** 2 / (main[wl][1] ** 2 + spec[wl][1] ** 2)
                   for wl in np.arange(wl0, wlf + 1) if not np.isnan(main[wl][1]) and not np.isnan(spec[wl][1])])
        # chi = sum([(spec[wl][0] - main[wl][0])**2/(spec[wl][1]**2) for wl in np.arange(wl0,wlf+1) if spec[wl][1]!=0])
        if normalized:
            chi = chi / len(wls)
        return chi

    def cpy( self ):
        """
        Returns a deep copy of this spectrum

        :rtype: Spectrum
        """
        from copy import deepcopy
        return deepcopy( self )

    def cpy_info( self ):
        """
        Resturns a spectrum devoid of wl/flux data, but with the same namestring, redshift and gmag as this one.

        :rtype: Spectrum
        """
        spec = Spectrum( ns=self.getNS( ), z=self.getRS( ), gmag=self.getGmag( ) )
        return spec

    def dim_to_ab( self, to_mag_ab: float, scale_wl: float = DEFAULT_SCALE_WL ) -> None:
        """
        Determines the desired flux that would be exhibited at a given AB Magnitude and wavelength,
        then passes that value to the Spectrum.scale() method, scaling the spectrum to that flux density
        
        :param to_mag_ab: AB Magnitude desired to be dimmed to
        :type to_mag_ab: float
        :param scale_wl: Wavelength to scale around.  Defaults to common.constants.DEFAULT_SCALE_WL
        :type scale_wl: float
        :return: None
        :rtype: None
        """
        exponent = (8.9 - to_mag_ab) / 2.5
        f_v = pow( 10, exponent )
        f_lambda = f_v / (3.34E4 * 1E-17 * pow( scale_wl, 2 ))
        self.scale( scaleflx=f_lambda )

    def getFlux( self, wavelength: float ) -> float:
        """
        :param wavelength:
        :type wavelength: float
        :return: The flux density at wavelength.  Equivalent to Spectrum[ wavelength ][ 0 ]
        :rtype: float
        """
        return self[ wavelength ][ 0 ]

    def getErr( self, wavelength: float ) -> float:
        """
        :param wavelength:
        :type wavelength: float
        :return: The flux density error at wavelength.  Equivalent to Spectrum[ wavelength ][ 1 ]
        :rtype: float
        """
        return self[ wavelength ][ 1 ]

    def getFluxlist( self ) -> List[ float ]:
        """
        Returns the flux densities in a list, ordered by wavelength
        
        :rtype: list 
        """
        return [ self.getFlux( wl ) for wl in self.getWavelengths( ) ]

    def getErrList( self ) -> List[ float ]:
        """
        Returns the flux densities in a list, ordered by wavelength

        :rtype: list 
        """
        return [ self.getErr( wl ) for wl in self.getWavelengths( ) ]

    def getGmag( self ) -> float:
        """
        Returns the magnitude in G filter of the spectrum

        :rtype: float
        """
        return self.__gmag

    def getNS( self ) -> str:
        """
        Returns the namestring of the spectrum object
        :rtype: str
        """
        return self.__namestring

    def getRS( self ) -> float:
        """
        Returns the stored redshift of the spectrum

        :rtype: float
        """
        return self.__z

    def getWavelengths( self ) -> List[ float ]:
        """
        :return: A list of the wavelengths in this Spectrum, sorted by increasing value
        :rtype: list
        """
        return sorted( self.keys( ) )

    def lineDict( self, wavelength: float ) -> dict:
        """
        Returns a simple dictionary of { wavelength : float, flux density : float, error : float }
        to allow writing with a CSV dict writer
        
        :param wavelength:
        :type wavelength: float
        :rtype: dict 
        """
        return { 'wavelength': wavelength, 'flux density': self.getFlux( wavelength ),
                 'error': self.getErr( wavelength ) }

    def lineDictList( self ) -> List[ dict ]:
        """
        Returns a list of Spectrum.lineDict values
        
        :rtype: list
        """
        return [ self.lineDict( wl ) for wl in self.getWavelengths( ) ]

    def magAB( self, wl_range: Tuple[ float, float ] = (None, None) ) -> float:
        """
        Determines the average AB Magnitude over the given band of interest.

        If wl_range is not specified, defaults to common.constants DEFAULT_SCALE_WL +/- DEFAULT_SCALE_RADIUS

        :param wl_range: Band range over which to determine AB magntiude
        :type wl_range: tuple
        :return: AB Magnitude
        :rtype: float
        """
        from numpy import mean, log10

        minwl = wl_range[ 0 ] or DEFAULT_SCALE_WL - DEFAULT_SCALE_RADIUS
        maxwl = wl_range[ 1 ] or DEFAULT_SCALE_WL + DEFAULT_SCALE_RADIUS

        f_vlist = list( )
        f_v = None
        for wl in self.getWavelengths( ):
            if minwl <= wl <= maxwl:
                f_v = 3.34E4 * pow( wl, 2 ) * 1E-17 * self[ wl ][ 0 ]
                f_vlist.append( f_v )
        try:
            f_v = mean( f_vlist )
        except RuntimeWarning as e:
            print(
                f"Spectrum.magAB(): {self.getNS()} got a RuntimeWarning when trying to form flux mean: {f_v} \n fluxlist: {f_vlist}" )
            raise e
        return -2.5 * log10( f_v ) + 8.9

    def nearest( self, wavelength: float ) -> float:
        """
        Simple wrapper for tools.find_nearest_wavelength.  Passes this Spectrum's sorted wavelength list (via
        Spectrum.getWavelengths()) and wavelength to find_nearest_wavelength, and returns that value.

        :param wavelength: Wavelength of interest
        :type wavelength: float
        :return: Value of the nearest wavelength
        :rtype: float
        """
        from spectrum.utils import find_nearest_wavelength
        return find_nearest_wavelength( self.getWavelengths( ), wavelength )

    def scale2(self, to, central_wl=4767, radius=18, inplace = False):
        """
        Scales the flux of this spectrum to a certain value
        :param to: Either a spectrum or a float
        :param central_wl:
        :param radius:
        :return:
        """

        if not isinstance(to, float):
            if isinstance(to,Spectrum):
                aveflux_to = to.aveFlux(central_wl, radius)
            else:
                print('Cannot scale. to must be either a Spectrum object or a float.')
                exit(1)
        else:
            aveflux_to = to
        r = aveflux_to / self.aveFlux(central_wl, radius)

        wls = self.getWavelengths()
        fluxlist = [r * f for f in self.getFluxlist()]
        errlist = [r * f for f in self.getErrList()]
        if inplace:
            self.setDict(wls, fluxlist, errlist)
            return None
        else:
            spec2 = self.cpy()
            spec2.setDict(wls, fluxlist, errlist)
            return spec2

    def setDict( self, wavelengthList, fluxList, errList ):
        """
        Replace the current wavelength dictionary with the data passed to method.

        :param wavelengthList: wavelength values
        :param fluxList: flux density values
        :param errList: flux density error values
        :type wavelengthList: list
        :type fluxList: list
        :type errList: list
        :return: None
        :rtype: None
        """
        self.clear( )
        for i in range( len( wavelengthList ) ):
            self[ wavelengthList[ i ] ] = (fluxList[ i ], errList[ i ])

    def setRS(self, redshift ):
        """
        Manually set the redshift of the spectrum

        :type redshift: float
        :return: None
        """
        assert type( redshift ) == float
        self.__z = redshift

    def setNS(self, namestring : str ):
        self.__namestring = namestring


    def shiftToRest( self, z: float = None, inplace = False ) -> None:
        """
        
        Shifts this spectrum to rest frame using z.  If z is None, uses the stored value of z.
        
        :param z: Redshift to use to shift to rest frame
        :type z: float
        :rtype: None
        """
        if z is None:
            z = self.__z

        wls = [ wl / (1 + z) for wl in self.getWavelengths( ) ]
        fluxlist = self.getFluxlist( )
        errlist = self.getErrList( )
        if inplace:
            self.clear( )
            for i in range( len( wls ) ):
                self[ wls[ i ] ] = (fluxlist[ i ], errlist[ i ])
            return
        else:
            return Spectrum(dict={wl : (f,e) for wl,f,e in zip(wls,fluxlist,errlist)},z=z, ns=self.__namestring)

    def scale( self, scale_spec=None, scaleflux: float = None, scaleWL: float = DEFAULT_SCALE_WL,
               radius: float = DEFAULT_SCALE_RADIUS ):
        """
        Simple scaling process.  At minimum, pass either scale_spec or scaleflux.  If scale_spec is passed, the
        scaling flux density will be determined from it via scale_spec.aveFlux().
        
        :param scale_spec: Spectrum to scale to.  If not used, pass scaleflux.
        :type scale_spec: Spectrum
        :param scaleflux: Flux density to scale to.  If not used, pass scale_spec
        :type scaleflux: float
        :param scaleWL: Central wavelength to scale around.  Defaults to common.constants.DEFAULT_SCALE_WL
        :type scaleWL: float
        :param radius: Radius around central wavelength.  Defaults to common.constants.DEFAULT_SCALE_RADIUS
        :type radius: float
        :rtype: Spectrum
        :raises: AssertionError
        """
        assert scale_spec is not None or scaleflux is not None

        if scale_spec is not None:
            scaleflux = scale_spec.aveFlux( scaleWL, radius )

        scalar = scaleflux / self.aveFlux( scaleWL, radius )
        if scalar == 1.0: return self
        for wl in self:
            flux, err = self[ wl ]
            self[ wl ] = (flux * scalar, err * scalar)
        return self

    def trim( self, wlLow: float = None, wlHigh: float = None ) -> None:
        """
        Deletes any values exclusive of the wlLow or wlHigh range.
        
        :param wlLow: Minimum wavelength to keep
        :type wlLow: float
        :param wlHigh: Maximum wavelength to keep
        :type wlHigh: float
        :return: None
        """
        for wl in self.getWavelengths():
            if wlLow is not None and wl < wlLow:
                del self[ wl ]
            elif wlHigh is not None and wlHigh < wl:
                del self[ wl ]


    def plot_i(self, label=None,lw=0.3, **kwargs):
        from matplotlib import pyplot
        y = self.getFluxlist()
        x = self.getWavelengths()
        pyplot.plot(x, y, lw=lw, label=label or str(self.getNS()), **kwargs)
        pyplot.xlabel('Wavelength ($Å$)')
        pyplot.ylabel('Flux Density ($erg \ cm^{-2}s^{-1} Å^{-1}$)')
        pyplot.legend()

    def residual(self,func,**kargs):
        """
        Uses func to find the background and residual of the spectrum.
        :param func: Callable operates on a flux array.
        :param func:
        :param kargs:
        :return:
        """
        bcgrnd = Spectrum(dict={wl: (flux, 0) for wl, flux in
                                zip(self.getWavelengths(), func(self.getFluxlist(), **kargs))},
                          ns='Background')
        res = {wl: (flux1 - flux2, 0) for wl, flux1, flux2 in
               zip(self.getWavelengths(), self.getFluxlist(), bcgrnd.getFluxlist())}
        res = Spectrum(dict=res, ns='Residual')
        return res, bcgrnd

    def step(self, label=None, where='pre',**kwargs):
        from matplotlib import pyplot
        y = self.getFluxlist()
        x = self.getWavelengths()
        pyplot.step(x, y, where=where, lw=.3, label=label or str(self.getNS()),**kwargs)
        pyplot.xlabel('Wavelength ($Å$)')
        pyplot.ylabel('Flux Density ($erg \ cm^{-2}s^{-1} Å^{-1}$)')
        pyplot.legend()


def from_fits(file, z=None, ns=None):
    file = fits.open(file)
    try:
        loglam = file[1].data['loglam']
        wls = [10 ** l for l in loglam]
    except:
        coef0 = file[0].header['COEFF0']
        coef1 = file[0].header['COEFF0']
        naxis = file[1].header['NAXIS2']
        wls = [10 ** (coef0 + i * coef1) for i in range(1, naxis)]
    if z is None:
        try:
            z = file[0].header['Z']
        except KeyError:
            lst = file[3].data['LINEZ']
            z = lst[lst != 0][0]

    flux = file[1].data['flux']
    err = file[1].data['ivar']
    dicti = {w: (a, b) for w, a, b in zip(wls, flux, err)}

    return Spectrum(dict=dicti, z=z, ns=ns or 'Spectrum')


def snip(x, niter, smoothWindow=None):
    """Implementation of the SNIP algorithm for estimating background
    in a noisy spectrum with peaks."""
    ssize = len(x)
    workspace = np.zeros(2 * ssize, dtype=float)
    workspace[0:ssize] = x
    workspace[ssize:] = x
    if smoothWindow is not None:
        if smoothWindow not in [3, 5, 7, 9, 11, 13, 15]:
            raise ValueError('smoothWindow={:d} must be in 3,5,7,..,15'.format(smoothWindow))
        bw = (smoothWindow - 1) // 2

    i = niter
    while i > 0:
        for j in range(i, ssize - i):
            if smoothWindow is None:
                a = workspace[ssize + j]
                b = 0.5 * (workspace[ssize + j - i] + workspace[ssize + j + i])
                if b < a:
                    a = b
                workspace[j] = a
            else:
                a = workspace[ssize + j]
                av, men = 0., 0.
                for w in range(j - bw, j + bw + 1):
                    if w >= 0 and w < ssize:
                        av += workspace[ssize + w]
                        men += 1
                av /= men

                b, men = 0., 0.
                for w in range(j - i - bw, j - i + bw + 1):
                    if w >= 0 and w < ssize:
                        b += workspace[ssize + w]
                        men += 1
                b /= men

                c, men = 0., 0.
                for w in range(j + i - bw, j + i + bw + 1):
                    if w >= 0 and w < ssize:
                        c += workspace[ssize + w]
                        men += 1
                c /= men

                b = (b + c) / 2
                if b < a:
                    av = b
                workspace[j] = av

        for j in range(i, ssize - i):
            workspace[ssize + j] = workspace[j]

        i -= 1
    return workspace[ssize:]


def average(speclist : List[Spectrum], ns='composite', full_range=False) -> Spectrum:
    from pandas import DataFrame, Series
    """
    Makes a composite spectrum from the given speclist by getting the geometric mean at each point. The composite wavelength
    range results from an 'outer' join of all wavelength ranges of the spectra. 
    :param speclist: List of spectra to be used.
    :param ns: Name for the output Spectrum
    :return:
    """
    df = DataFrame()
    if not full_range:
        df.dropna(inplace=True)

    for spec in speclist:
        df[spec.getNS()] = Series(index=spec.getWavelengths(), data=spec.getFluxlist())
    flux = df.mean(axis=1).values
    error = df.std(axis=1).values
    wls = df.index
    d = {wl: (fl, err) for wl,fl,err in zip(wls,flux,error)}
    return Spectrum(dict=d, ns=ns)
