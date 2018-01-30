http://www.eso.org/~rsiebenm/agn_models/



SED library
A complete set of 3600 SEDs can be downloaded as zip files in plain ASCII format or as an IDL structure. In ASCII format the data have two columns. The first is the wavelength in microns, the second column is the flux in Jy. SEDs are computed for AGNs at a distance of 50 Mpc and a luminosity of 10^11 Lo.

Model grid
We consider five basic model paramters:
The inner radius of the dusty torus: R = 300, 514, 772, 1000, 1545 in units: [10^15 cm]
The cloud volume filling factor: Vc = 1.5, 7.7, 38.5, 77.7 (%).
The optical depth (in V) of the individual clouds: Ac = 0, 4.5, 13.5, 45.
The optical depth (in V) of the disk midplane: Ad = 0, 30, 100, 300, 1000.
The viewing angle: th = 1 ,.., 9 (corresponding to bins at 86, 80, 73, 67, 60, 52, 43, 33, and 19 degree measured from the pole (z-axis)

File notation:
Rxxxx_Vcxxx_Acxxxx_Adxxxx.thx 
R1545_Vc777_Ac4050_Ad1000.th9






Ralf.Siebenmorgen@eso.org











By extensive testing
I found that
(1) low cloud filling factor + large viewing angle = high cloud filling factor + large viewing angle
(2) viewing angle 3-8 are very similar
(3) parameter 2 and 3 are very similar
...

So I decide to take 4 models:

    get_one_sed 0 0 1 1 0  -- low filling factor, zero viewing angle (similar to large viewing angle)
    get_one_sed 4 3 3 3 0  -- high filling factor, zero viewing angle (similar to large viewing angle)
    get_one_sed 4 2 1 3 2  -- mid filling factor, intermidate viewing angle
    get_one_sed 4 2 1 3 8  -- mid filling factor, large viewing angle




