import numpy as np
from astropy.io import fits
from astropy import units as u
import os

class LSSInputs:
    def __init__(self):
        self.slit_length = 2.*u.deg
        self.slit_width = 2.*u.arcsec
        self.x_0 = 3.5*u.deg
        self.y_0 = 0.*u.deg
        self.pixel_scale = 0.80*u.arcsec
        self.plate_scale = 80.*u.arcsec/u.mm
        self.gap_size = 100 # in pixels
        self.num_pixels = 4096 # in spatial direction
    
    @staticmethod
    def make_spectral_efficiency(infile="inputs/1150_3550_1000_4p3_1420.txt", outfile="UVIM_LSS_spectral_efficiency.fits"):
        # first convert the spectral efficiency file to fits file
        spec_eff = np.loadtxt(infile)
        spec_eff_dict = {"wavelength": spec_eff[:, 0] * u.nm, "efficiency": spec_eff[:, 1]}
        # convert from nm to microns
        spec_eff_dict["wavelength"] = spec_eff_dict["wavelength"].to(u.um).value

        # only one trace
        # required fits structure is located in spectral_efficiency in scopesim
        hdu0 = fits.PrimaryHDU()
        hdu0.header["ECAT"] = 1
        hdu0.header["EDATA"] = 2
        hdu0.header["DATE"] = np.datetime64('today', 'D').astype(str)
        hdu1 = fits.BinTableHDU.from_columns(
            [fits.Column(name="description", format="20A", array=["UVIM_LSS_spectral_efficiency"]),
            fits.Column(name="extension_id", format="I", array=[2])]
        )
        hdu2 = fits.BinTableHDU.from_columns(
            [fits.Column(name="wavelength", format="E", array=spec_eff_dict["wavelength"]),
            fits.Column(name="efficiency", format="E", array=spec_eff_dict["efficiency"])]
        )
        hdu2.header["EXTNAME"] = "UVIM_LSS_spectral_efficiency"
        hdul = fits.HDUList([hdu0, hdu1, hdu2])
        hdul.writeto(outfile, overwrite=True)

    def make_slit_geometry(self, outfile="UVIM_LSS_slit_geometry.dat"):
        slit_length = (self.slit_length).to(u.arcsec) # (1 deg/72 mm ?)
        slit_width = (self.slit_width).to(u.arcsec) # 2 pixels or 20 microns
        # relative to the field, located at 3.5 deg in x direction, and centered in y direction +/- 0.5 deg
        # need four coords to define rectangular aperture
        # y is the spatial direction, x is the spectral
        x_0 = (self.x_0).to(u.arcsec)
        y_0 = (self.y_0).to(u.arcsec)
        slit_coords = np.array([[x_0.value - slit_width.value/2, y_0.value - slit_length.value/2],
                                [x_0.value - slit_width.value/2, y_0.value + slit_length.value/2],
                                [x_0.value + slit_width.value/2, y_0.value + slit_length.value/2],
                                [x_0.value + slit_width.value/2, y_0.value - slit_length.value/2]])
        # write to dat file (but don't overwrite if it already exists)
        if not os.path.exists(outfile):
            np.savetxt(outfile, slit_coords, fmt="%f", delimiter="    ", header="x    y")
        else:
            # just write to a different text file
            np.savetxt(outfile.replace(".dat", "_new.dat"), slit_coords, fmt="%f", delimiter="    ", header="x    y")

    def make_spectral_trace(self, slit_geometry="UVIM_LSS_slit_geometry.dat", 
                            infile="inputs/UVEXS_Spectral_Resolution_R2000.txt", 
                            outfile="UVIM_LSS_spectral_trace.fits"):
        # 2 4k by 4k detectors with long side along the spectral direction, so 4096 pixels in spatial direction 
        # and 8192 pixels in spectral direction
        # spectral file contains wavelength to position mapping
        # set y as dispersion direction
        # eventually should check that the set_dispersion function in the SpectralTraceList is consistent with the true dispersion values 
        data = np.loadtxt(infile, skiprows=2, unpack=True)
        wavelength = data[0] * u.nm
        y_pos = data[1] * u.mm
        dispersion = data[2] * u.nm # per pixel
        wavelength = wavelength.to(u.um) # convert to microns

        # get slit geometry in spatial direction for centering the trace
        slit_coords = np.loadtxt(slit_geometry, skiprows=1)
        slit_s_min = np.min(slit_coords[:,1]) # in arcsec
        slit_s_max = np.max(slit_coords[:,1]) # in arcsec
        slit_s_center = (slit_s_min + slit_s_max) / 2 
        # assume the slit is centered on detector, so 2048 pixels in each direction
        s_min = -self.num_pixels/2 * self.pixel_scale + slit_s_center # in arcsec
        s_max = self.num_pixels/2 * self.pixel_scale + slit_s_center # in arcsec
        x_det_min = s_min / self.plate_scale # in mm
        x_det_max = s_max / self.plate_scale # in mm
        print(wavelength)
        print(y_pos)
        # assume linear mapping between s and x
        s_grid = np.linspace(s_min, s_max, len(wavelength)) # in arcsec
        x_grid = np.linspace(x_det_min, x_det_max, len(wavelength)) # in mm
        print(s_grid)
        print(x_grid)
        # write to fits file
        hdu0 = fits.PrimaryHDU()
        hdu0.header["ECAT"] = 1
        hdu0.header["EDATA"] = 2
        hdu1 = fits.BinTableHDU.from_columns(
            [fits.Column(name="description", format="20A", array=["UVIM_LSS_trace"]),
            fits.Column(name="extension_id", format="I", array=[2]),
            fits.Column(name="aperture_id", format="I", array=[0]),
            fits.Column(name="image_plane_id", format="I", array=[0])]
        )
        hdu2 = fits.BinTableHDU.from_columns(
            [fits.Column(name="wavelength", format="E", array=wavelength.value),
            fits.Column(name="s", format="E", array=s_grid.value),
            fits.Column(name="x", format="E", array=x_grid.value),
            fits.Column(name="y", format="E", array=y_pos.value)]
        )
        hdu2.header["EXTNAME"] = "UVIM_LSS_trace"
        hdu2.header["DISPDIR"] = "y"
        hdu2.header["x_unit"] = "mm"
        hdu2.header["y_unit"] = "mm"
        hdul = fits.HDUList([hdu0, hdu1, hdu2])
        hdul.writeto(outfile, overwrite=True)

    @staticmethod
    def make_filter_response(infile="inputs/graded_overcoat_00nm.csv", outfile="UVIM_LSS_filter_response.dat"):
        # filter response file contains wavelength to transmission mapping
        data = np.loadtxt(infile, skiprows=1, unpack=True, delimiter=",")
        wavelength = data[0] * u.nm
        transmission = data[1]
        transmission = np.array(transmission) / 100.0 # convert from percentage to fraction

        if not os.path.exists(outfile):
            with open(outfile, 'w') as f:
                f.write("wavelength    transmission\n")
                for wl, trans in zip(wavelength, transmission):
                    f.write(f"{wl.value}    {trans}\n")
        else:
            with open(outfile.replace(".dat", "_new.dat"), 'w') as f:
                f.write("wavelength    transmission\n")
                for wl, trans in zip(wavelength, transmission):
                    f.write(f"{wl.value}    {trans}\n")

if __name__ == "__main__":
    # run  python3 make_LSS_inputs.py from command line
    lss_inputs = LSSInputs()
    lss_inputs.make_spectral_efficiency()
    lss_inputs.make_slit_geometry()
    lss_inputs.make_spectral_trace()
    lss_inputs.make_filter_response()