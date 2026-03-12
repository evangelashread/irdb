import numpy as np
from astropy.io import fits
from astropy import units as u
import os

class LSSInputs:
    def __init__(self):
        self.slit_length = 1.*u.deg
        self.slit_width = 2.*u.arcsec
        self.y_0 = 0.*u.deg
        self.x_0 = 0.*u.deg 
        self.pixel_scale = 0.80*u.arcsec
        self.plate_scale = 80.*u.arcsec/u.mm
        self.gap_size = 100 # in pixels, haven't done anything with this yet
        self.num_pixels = 4096 # in spatial direction
        self.inputs_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "inputs/"))
        self.outputs_dir = os.path.abspath(os.path.dirname(__file__))
    
    def make_spectral_efficiency(self, infile="1150_3550_1000_4p3_1420.txt", outfile="UVIM_LSS_spectral_efficiency.fits"):
        # first convert the spectral efficiency file to fits file
        spec_eff = np.loadtxt(os.path.join(self.inputs_dir, infile))
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
            [fits.Column(name="description", format="20A", array=["UVIM_LSS_trace"]),
            fits.Column(name="extension_id", format="I", array=[2])]
        )
        hdu2 = fits.BinTableHDU.from_columns(
            [fits.Column(name="wavelength", format="E", array=spec_eff_dict["wavelength"]),
            fits.Column(name="efficiency", format="E", array=spec_eff_dict["efficiency"])]
        )
        hdu2.header["EXTNAME"] = "UVIM_LSS_trace"
        hdul = fits.HDUList([hdu0, hdu1, hdu2])
        if not os.path.exists(os.path.join(self.outputs_dir, outfile)):
            hdul.writeto(os.path.join(self.outputs_dir, outfile))
        else:
            # write to different file
            outfile_new = outfile.replace(".fits", "_new.fits")
            hdul.writeto(os.path.join(self.outputs_dir, outfile_new))

    def make_slit_geometry(self, outfile="UVIM_LSS_slit_geometry.dat"):
        slit_length = (self.slit_length).to(u.arcsec) # (1 deg/72 mm ?)
        slit_width = (self.slit_width).to(u.arcsec) # 2 pixels or 20 microns
        # relative to the field, located at 3.5 deg in y direction, and centered in x direction +/- 0.5 deg
        # need four coords to define rectangular aperture
        # x is the spatial direction, y is the spectral (to be consistent with ScopeSim)
        x_0 = (self.x_0).to(u.arcsec)
        y_0 = (self.y_0).to(u.arcsec)
        slit_coords = np.array([[x_0.value - slit_length.value/2, y_0.value - slit_width.value/2],
                                [x_0.value + slit_length.value/2, y_0.value - slit_width.value/2],
                                [x_0.value + slit_length.value/2, y_0.value + slit_width.value/2],
                                [x_0.value - slit_length.value/2, y_0.value + slit_width.value/2]])
        # write to dat file (but don't overwrite if it already exists)
        if not os.path.exists(os.path.join(self.outputs_dir, outfile)):
            with open(os.path.join(self.outputs_dir, outfile), 'w') as f:
                f.write("x    y\n")
                for x, y in zip(slit_coords[:,0], slit_coords[:,1]):
                    f.write(f"{x}    {y}\n")
        else:
            # just write to a different dat file
            outfile_new = outfile.replace(".dat", "_new.dat")
            with open(os.path.join(self.outputs_dir, outfile_new), 'w') as f:
                f.write("x    y\n")
                for x, y in zip(slit_coords[:,0], slit_coords[:,1]):
                    f.write(f"{x}    {y}\n")

    def make_spectral_trace(self, slit_geometry="UVIM_LSS_slit_geometry.dat", 
                            infile="UVEXS_Spectral_Resolution_R2000.txt", 
                            outfile="UVIM_LSS_spectral_trace.fits",
                            n_slit_positions=400):
        
        data = np.loadtxt(os.path.join(self.inputs_dir, infile), skiprows=2, unpack=True)
        wavelength = data[0] * u.nm
        y_pos = (data[1] * u.mm).value
        wavelength = wavelength.to(u.um).value # convert to microns

        # get slit geometry in spatial direction for centering the trace
        slit_coords = np.loadtxt(os.path.join(self.outputs_dir, slit_geometry), skiprows=1)
        # Determine which direction is spatial (longer dimension)
        x_extent = abs(slit_coords[:,0].max() - slit_coords[:,0].min())
        y_extent = abs(slit_coords[:,1].max() - slit_coords[:,1].min())
        spatial_col = 0 if x_extent > y_extent else 1  # 0=horizontal, 1=vertical
        slit_s_min = np.min(slit_coords[:,spatial_col]) # in arcsec
        slit_s_max = np.max(slit_coords[:,spatial_col]) # in arcsec
        slit_s_center = (slit_s_min + slit_s_max) / 2 
        
        # assume the slit is centered on detector, so 2048 pixels in each direction
        s_min = -self.num_pixels/2 * self.pixel_scale.value + slit_s_center # in arcsec
        s_max = self.num_pixels/2 * self.pixel_scale.value + slit_s_center # in arcsec
        x_det_min = s_min / self.plate_scale.value # in mm
        x_det_max = s_max / self.plate_scale.value # in mm

        # for a long-slit spectrograph, each position in the slit creates a vertical trace
        # this means we effectively have a grid of traces
        s_positions = np.linspace(s_min, s_max, n_slit_positions) # in arcsec
        x_positions = np.linspace(x_det_min, x_det_max, n_slit_positions) # in mm
        
        # grid w/ N_slit_positions * N_wavelengths rows
        # y varies with wavelength, but s and x do not
        wavelength_grid = np.tile(wavelength, n_slit_positions)
        y_grid = np.tile(y_pos, n_slit_positions)
        s_grid = np.repeat(s_positions, len(wavelength)) # in arcsec
        x_grid = np.repeat(x_positions, len(wavelength)) # in mm
        # write to fits file in the format SpectralTraceList expects
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
            [fits.Column(name="wavelength", format="E", array=wavelength_grid),
            fits.Column(name="s", format="E", array=s_grid),
            fits.Column(name="x", format="E", array=x_grid),
            fits.Column(name="y", format="E", array=y_grid)]
        )
        hdu2.header["EXTNAME"] = "UVIM_LSS_trace"
        hdu2.header["DISPDIR"] = "y"
        hdu2.header["TUNIT1"] = "um"
        hdu2.header["TUNIT2"] = "arcsec"
        hdu2.header["TUNIT3"] = "mm"
        hdu2.header["TUNIT4"] = "mm"
        hdu2.header["WAVECOLN"] = "wavelength"
        hdu2.header["SLITPOSN"] = "s"
        hdul = fits.HDUList([hdu0, hdu1, hdu2])
        if not os.path.exists(os.path.join(self.outputs_dir, outfile)):
            hdul.writeto(os.path.join(self.outputs_dir, outfile))
        else:
            # write to different file
            outfile_new = outfile.replace(".fits", "_new.fits")
            hdul.writeto(os.path.join(self.outputs_dir, outfile_new))

    def make_filter_response(self, infile="graded_overcoat_00nm.csv", outfile="UVIM_LSS_filter_response.dat"):
        # filter response file contains wavelength to transmission mapping
        data = np.loadtxt(os.path.join(self.inputs_dir, infile), skiprows=1, unpack=True, delimiter=",")
        wavelength = data[0] * u.nm
        transmission = data[1]
        transmission = np.array(transmission) / 100.0 # convert from percentage to fraction

        if not os.path.exists(os.path.join(self.outputs_dir, outfile)):
            with open(os.path.join(self.outputs_dir, outfile), 'w') as f:
                f.write("wavelength    transmission\n")
                for wl, trans in zip(wavelength, transmission):
                    f.write(f"{wl.value}    {trans}\n")
        else:
            outfile_new = outfile.replace(".dat", "_new.dat")
            with open(os.path.join(self.outputs_dir, outfile_new), 'w') as f:
                f.write("wavelength    transmission\n")
                for wl, trans in zip(wavelength, transmission):
                    f.write(f"{wl.value}    {trans}\n")
                    
    def make_dispersion_file(self, infile="UVEXS_Spectral_Resolution_R2000.txt", outfile="UVIM_LSS_dispersion.dat"):
        data = np.loadtxt(os.path.join(self.inputs_dir, infile), skiprows=2, unpack=True)
        wavelength = data[0] * u.nm
        dispersion = data[2] * u.nm # per pixel
        wavelength = wavelength.to(u.um) # convert to microns
        dispersion = dispersion.to(u.um) # convert to microns per pixel

        # write to dat file (but don't overwrite if it already exists)
        if not os.path.exists(os.path.join(self.outputs_dir, outfile)):
            with open(os.path.join(self.outputs_dir, outfile), 'w') as f:
                f.write("# wavelength_unit: um\n")
                f.write("# dispersion_unit: um\n")
                f.write("wavelength    dispersion\n")
                for wl, d in zip(wavelength, dispersion):
                    f.write(f"{wl.value}    {d.value}\n")
        else:
            outfile_new = outfile.replace(".dat", "_new.dat")
            with open(os.path.join(self.outputs_dir, outfile_new), 'w') as f:
                f.write("# wavelength_unit: um\n")
                f.write("# dispersion_unit: um\n")
                f.write("wavelength    dispersion\n")
                for wl, d in zip(wavelength, dispersion):
                    f.write(f"{wl.value}    {d.value}\n")
                    
class NUVInputs:
    def __init__(self):
        self.pixel_scale = 1.03 * u.arcsec
        self.plate_scale = 103.0 * u.arcsec / u.mm
        self.inputs_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "inputs/"))
        self.outputs_dir = os.path.abspath(os.path.dirname(__file__))
        
    def make_qe_curve(self, infile="nuv_qe_Hf02.csv", outfile="UVIM_NUV_QE.dat"):
        data = np.loadtxt(os.path.join(self.inputs_dir, infile), delimiter=',', skiprows=4, unpack=True)
        wavelength = data[0] * u.nm
        qe = data[3] # already a fraction
        wavelength = wavelength.to(u.um) # convert to microns

        if not os.path.exists(os.path.join(self.outputs_dir, outfile)):
            with open(os.path.join(self.outputs_dir, outfile), 'w') as f:
                f.write("# wavelength_unit: um\n")
                f.write("wavelength    transmission\n")
                for wl, q in zip(wavelength, qe):
                    f.write(f"{wl.value}    {q}\n")
        else:
            outfile_new = outfile.replace(".dat", "_new.dat")
            with open(os.path.join(self.outputs_dir, outfile_new), 'w') as f:
                f.write("# wavelength_unit: um\n")
                f.write("wavelength    transmission\n")
                for wl, q in zip(wavelength, qe):
                    f.write(f"{wl.value}    {q}\n")
    
    def make_dichroic_response(self, infile="dichroic_bandpass.csv", outfile="UVIM_dichroic_response.dat"):
        # Note: this same file should be used for the FUV surfaces list, too
        data = np.loadtxt(os.path.join(self.inputs_dir, infile), delimiter=',', skiprows=4, unpack=True)
        wavelength = data[0] * u.nm
        wavelength = wavelength.to(u.um) # convert to microns
        reflection = data[1]
        transmission = data[2] # already a fraction
        
        if not os.path.exists(os.path.join(self.outputs_dir, outfile)):
            with open(os.path.join(self.outputs_dir, outfile), 'w') as f:
                f.write("# wavelength_unit: um\n")
                f.write("wavelength    reflection    transmission\n")
                for wl, re, tr in zip(wavelength, reflection, transmission):
                    f.write(f"{wl.value}    {re}    {tr}\n")
        else:
            outfile_new = outfile.replace(".dat", "_new.dat")
            with open(os.path.join(self.outputs_dir, outfile_new), 'w') as f:
                f.write("# wavelength_unit: um\n")
                f.write("wavelength    reflection    transmission\n")
                for wl, re, tr in zip(wavelength, reflection, transmission):
                    f.write(f"{wl.value}    {re}    {tr}\n")
        
if __name__ == "__main__":
    # run python3 make_LSS_inputs.py from command line
    
    lss_inputs = LSSInputs()
    lss_inputs.make_spectral_efficiency()
    lss_inputs.make_slit_geometry()
    lss_inputs.make_spectral_trace()
    lss_inputs.make_filter_response()
    lss_inputs.make_dispersion_file()
    
    nuv_inputs = NUVInputs()
    nuv_inputs.make_qe_curve()
    nuv_inputs.make_dichroic_response()