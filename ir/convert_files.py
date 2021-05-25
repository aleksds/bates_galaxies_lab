import numpy as np
import os
import glob
from astropy.io import ascii
import shutil

projpath = os.getcwd() + '/Brown2014/An_Atlas_of_Galaxy_SEDs/An_Atlas_of_Galaxy_SEDs/'
files = glob.glob(projpath + '*.dat')
new_projpath = projpath + "Converted/"

def remove_files(file_extension, path = None):
    if path == None:
        r_proj_path = os.getcwd() + "/"
        r_files = glob.glob(r_proj_path + file_extension)
    else:
        r_proj_path = os.getcwd() + path + "/"
        r_files = glob.glob(r_proj_path + file_extension)
    for r in r_files:
        os.remove(r)


#remove_files('*.pdf' )

for i in range(0, len(files)):

    table = ascii.read(files[i])

    wave = np.array(table['col1'])
    micro_wave = wave * 1e-4 #Converts to Microns from Angstroms
    flam = np.array(table['col2'])
    flux = wave * flam
    fnu = flam * wave ** 2 / 3e18  # * 1e46
    con_data = [wave, fnu, table[2]]
    broken_path = files[i].split("/")
    remove_filetype = broken_path[len(broken_path) - 1].split(".")
    new_name = remove_filetype[0] + "_conv.dat"
    print(new_name)
    ascii.write([micro_wave,fnu], new_name, names=['Column 1: Rest Wavelength (microns)','Column 2: Flux (ergs/s/cm^2/Hz)'])
    if os.path.exists(new_projpath):
        file_path = os.getcwd() + "/" + new_name
        shutil.move(file_path, new_projpath)