# -*- coding: utf-8 -*-
# =============================================================================
# Created By  : Benjamin Wriedt
# Copyright   : Copyright 2021, Institute of Chemical Engineering, Prof. Dr. Dirk Ziegenbalg, Ulm University
# License     : GNU LGPL
# ============================================================================
"""This program is part of the evaluation procedure of actinometry developed 
at the Institute of Chemical Engineering, University Ulm.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
#=======================
#Imports
import numpy as np  # import of numpy
from scipy.integrate import odeint
from scipy.optimize import minimize as scymin
from scipy import constants  # physical constants
import math  # import for rounding
import datetime
from lmfit import Parameter, Parameters, report_fit, fit_report
from lmfit import minimize
from lmfit import Minimizer as Lmfitminer


# Program to solve the ODE
def actin_feox(ls_file, quant_file, P_rad, V_inj, X_all, t_all, **kwargs): # input of experimental data

    y0 = kwargs.get('ode_y0', [1e-5])  # initial guess for odeint, variation helps if numerics have problems
    initial_params = kwargs.get('mini_ini', np.array([0.00]))
    P_all = kwargs.get('P_all', len(X_all)*[1])
    bounds = kwargs.get('mini_bound', [(0, 1)]) # defines parameter liimts for minimize optimization
    cropped_filename = quant_file[3:]
    kappa_file = 'kappa' + str(cropped_filename)
    file_out = kwargs.get('file_out', False)
    c_act = kwargs.get('c_act', 0.08)  # mol/L concentration of actinometer
    l = kwargs.get('l', 0.15875)  # cm optical path length
    comment = kwargs.get('comment', '-')  # comment of single measurement
    shell_out = kwargs.get('shell_out', True)

    # preparation of output file
    if file_out is True and shell_out is True: 
        date = datetime.datetime.now()
        def write(x):
            with open('acFeOxout' + date.strftime('%Y%m%d%H%H%S') + comment + '.txt', 'a') as out:
                out.write(x + '\n')
            print(x)
    if file_out is True and shell_out is False:
        date = datetime.datetime.now()
        def write(x):
            with open('acFeOxout' + date.strftime('%Y%m%d%H%H%S') + comment + '.txt', 'a') as out:
                out.write(x + '\n')
        write('#Output File for the Actinometric Calculations with Ferrioxalate\nComment: ' + comment +
              '\nLamp File Used: ' + ls_file + '\nRadiant Power: ' + str(P_rad) + '\nDate: ' + str(date))
    if file_out is False and shell_out is True:
        def write(x):
            print(x)
    else:
        def write(x):
            pass


    # csv import from ls_file
    ls_data = np.genfromtxt(ls_file, delimiter=",")
    wavelength_ls_raw = ls_data[:, 0]
    arb_intensity_raw = ls_data[:, 1]
    trans_FEP_raw = ls_data[:, 2]

    # csv import from quant_file
    quant_data = np.genfromtxt(quant_file, delimiter=";")
    wavelength_quant_raw = quant_data[:, 0]
    quantum_yield_raw = quant_data[:, 1]

    # csv import from kappa_file
    kappa_data = np.genfromtxt(kappa_file, delimiter=";")
    wavelength_kappa_raw = kappa_data[:, 0]
    absorption_coefficient_raw = kappa_data[:, 1]

    #convert wavelength to naturals and interpolation of corresponding values
    x1 = math.ceil(wavelength_ls_raw[0])
    x2 = math.floor(wavelength_ls_raw[-1]+1)
    wavelength_ls = np.arange(x1, x2, 1)  # set up wavelengths as evenly spaced
    arb_intensity = np.interp(wavelength_ls, wavelength_ls_raw, arb_intensity_raw)
    trans_FEP = np.interp(wavelength_ls, wavelength_ls_raw, trans_FEP_raw)

    x3 = math.ceil(wavelength_quant_raw[0])
    x4 = math.floor(wavelength_quant_raw[-1]+1)
    wavelength_quant = np.arange(x3, x4, 1)  # set up wavelengths as evenly spaced
    quantum_yield = np.interp(wavelength_quant, wavelength_quant_raw, quantum_yield_raw)

    x5 = math.ceil(wavelength_kappa_raw[0])
    x6 = math.floor(wavelength_kappa_raw[-1]+1)
    wavelength_kappa = np.arange(x5, x6, 1)  # set up wavelengths as evenly spaced
    absorption_coefficient = np.interp(wavelength_kappa, wavelength_kappa_raw, absorption_coefficient_raw)

    # spectral calculations
    energy_of_photons = constants.c * constants.h / wavelength_ls / 1e-9  # J
    fraction_of_radiant_flux = arb_intensity / sum(arb_intensity)  # 1
    spectral_radiant_flux = fraction_of_radiant_flux * P_rad  # W/nm
    spectral_photon_flux = spectral_radiant_flux / (energy_of_photons * constants.N_A)  # mol/(s nm)
    fraction_of_spectral_photon_flux = spectral_photon_flux / sum(spectral_photon_flux)  # 1
    total_photon_flux = sum(spectral_photon_flux)  # mol/s

    write("total_photon_flux according to the given technical data: "+str(total_photon_flux)+" mol/s")

    # get index of defined max and min values
    def min_max(data, mi, ma):
        pos = [0, 0]
        for counter, option in enumerate(data):
            if option >= mi and pos[0] == 0:
                pos[0] = counter
            elif counter == len(data):
                pos[0] = 0
            if option >= ma and pos[1] == 0 and counter != len(data):
                pos[1] = len(data)
            else:
                pos[1] = counter
        return pos

    # cut to common wavelength range (discards ls, kappa and quantum yield not mutually available)
    def common_member(wavelength_ls, wavelength_quant, wavelength_kappa, arb_intensity, trans_FEP, quantum_yield, absorption_coefficient, spectral_photon_flux, fraction_of_spectral_photon_flux):
        wavelength_ls = wavelength_ls.tolist()
        wavelength_quant = wavelength_quant.tolist()
        wavelength_kappa = wavelength_kappa.tolist()
        arb_intensity = arb_intensity.tolist()
        trans_FEP = trans_FEP.tolist()
        absorption_coefficient = absorption_coefficient.tolist()
        quantum_yield = quantum_yield.tolist()
        spectral_photon_flux = spectral_photon_flux.tolist()
        fraction_of_spectral_photon_flux = fraction_of_spectral_photon_flux.tolist()

        a_set = set(wavelength_ls)
        b_set = set(wavelength_quant)
        c_set = set(wavelength_kappa)
        same_elements = list(a_set & b_set & c_set)
        same_elements.sort()

        low_cut_a = wavelength_ls.index(same_elements[0])
        high_cut_a = wavelength_ls.index(same_elements[-1])+1 #+1 for list starting at position 0
        low_cut_b = wavelength_quant.index(same_elements[0])
        high_cut_b = wavelength_quant.index(same_elements[-1])+1
        low_cut_c = wavelength_kappa.index(same_elements[0])
        high_cut_c = wavelength_kappa.index(same_elements[-1])+1

        wavelength_ls_cut = wavelength_ls[low_cut_a:high_cut_a]
        arb_intensity_cut = arb_intensity[low_cut_a:high_cut_a]
        trans_FEP_cut = trans_FEP[low_cut_a:high_cut_a]
        spectral_photon_flux_cut = spectral_photon_flux [low_cut_a:high_cut_a]
        fraction_of_spectral_photon_flux_cut = fraction_of_spectral_photon_flux [low_cut_a:high_cut_a]
        absorption_coefficient_cut = absorption_coefficient[low_cut_c:high_cut_c]
        wavelength_kappa_cut = wavelength_kappa[low_cut_c:high_cut_c]
        wavelength_quant_cut = wavelength_quant[low_cut_b:high_cut_b]
        quantum_yield_cut = quantum_yield[low_cut_b:high_cut_b]

        cut_dict = {'wavelength_ls_cut': wavelength_ls_cut, 'wavelength_quant_cut': wavelength_quant_cut, 'wavelength_kappa_cut': wavelength_kappa_cut, 'arb_intensity_cut': arb_intensity_cut, 'trans_FEP_cut': trans_FEP_cut, 'quantum_yield_cut': quantum_yield_cut, 'absorption_coefficient_cut': absorption_coefficient_cut, 'spectral_photon_flux_cut': spectral_photon_flux_cut, 'fraction_of_spectral_photon_flux_cut': fraction_of_spectral_photon_flux_cut}
        return cut_dict

    # call of function, transformation and back-transformation to lists (function) and arrays (afterwards)
    cut_dict = common_member(wavelength_ls, wavelength_quant, wavelength_kappa, arb_intensity, trans_FEP, quantum_yield, absorption_coefficient, spectral_photon_flux, fraction_of_spectral_photon_flux)

    arb_intensity = cut_dict['arb_intensity_cut']
    absorption_coefficient = cut_dict['absorption_coefficient_cut']
    trans_FEP = cut_dict['trans_FEP_cut']
    quantum_yield = cut_dict['quantum_yield_cut']
    spectral_photon_flux = cut_dict['spectral_photon_flux_cut']
    fraction_of_spectral_photon_flux = cut_dict['fraction_of_spectral_photon_flux_cut']

    absorption_coefficient = np.array(absorption_coefficient)
    trans_FEP = np.array(trans_FEP)
    quantum_yield = np.array(quantum_yield)
    spectral_photon_flux = np.array(spectral_photon_flux)
    fraction_of_spectral_photon_flux = np.array(fraction_of_spectral_photon_flux)
    specific_photon_flux = sum(spectral_photon_flux)
    print("specific_photon_flux between "+str(wavelength[0])+" nm and "+str(wavelength[-1])+" nm according to the given technical data: "+str(specific_photon_flux)+" mol/s")

    # differential equation to determine irradiated photon flux
    def difeq(x, t, q, p):
        dxdt = (q / (c_act * V_inj)) * sum(trans_FEP * quantum_yield * fraction_of_spectral_photon_flux * (1 - np.exp(-absorption_coefficient * c_act * l * (1 - x))))
        return dxdt  # watch out! q is in mmol/s

    def fcn2min(params, P_all, t_all, data): #objective function, returns array to be minimized
        model = []
        #unpack parameters
        parvals = params.valuesdict()
        q = parvals['q']
        for m, o in zip(P_all, t_all):  # optimization to multiple points at all irradiation times, not only one point
            t_range = np.linspace(0, o, 100)  # define 100 points for the given time interval (required for numerics)
            fit_data = odeint(difeq, y0, t_range, args=(q, m))[:, 0][-1] # calculate the fitted data (100 data points)
            model.append(fit_data)  #list with results from odeint
        return np.asarray(model) - np.asarray(data) #return value to be minimized via least squares

    params = Parameters()
    params.add("q", value=0.00, min=-1.00, max=1.00) #photon flux q as parameter to optimize

    # define experimental conditions:
    general_params = {"X": X_all, "P": P_all, "t": t_all}
    
    #experimental data for X:
    data = list(general_params["X"])
    minner = Lmfitminer(fcn2min, params, fcn_args=(P_all, t_all, data)) #do fit with least squares model
    result = minner.minimize()

    # calculate final result
    final = data + result.residual

    # optional: write error report
    #report_fit(result)

    #unpack parameter
    qval = result.params.valuesdict()  #dict with q and its value
    q_calc = qval['q']

    # print all data returned by minimize function
    write("absorbed photon flux: "+str(q_calc/1000)+" mol/s")

    # calculate conversions predicted by the model by using the photon flux previously determined
    x_calc = []
    for i, j, k in zip(general_params["P"], general_params["X"], general_params["t"]):
        t_i = np.linspace(0, k, 100)
        x_calc.append(odeint(difeq, y0, t_i, args=(q_calc, i))[:, 0][-1])
    # print lists of calculated and experimental data
    X_all_round = ["%.3f" % x for x in X_all] # Round results to 3 digits
    x_calc_round = ["%.3f" % x for x in x_calc] # Round results to 3 digits
    std_exp_model = np.std(np.asarray(X_all_round,dtype=float)-np.asarray(x_calc_round,dtype=float))

    ext_phot_effi = float(q_calc/(total_photon_flux*1000))  #  WATCH OUT! q_calc ist in mmol/s
    write("Total external photonic efficiency: "+str(round(ext_phot_effi, 3)))
    specific_phot_effi = float(q_calc / (specific_photon_flux * 1000))  # WATCH OUT! q_calc ist in mmol/s
    write("External photonic efficiency between "+str(wavelength[0])+" nm and "+str(wavelength[-1])+" nm is "+ str(round(specific_phot_effi, 3)))

    info_dict={'ls_file': ls_file, 'kappa_file': kappa_file}

    # Write results in output file
    out = (comment, q_calc/1000, specific_phot_effi, '{:0.3e}'.format(float(std_exp_model)), total_photon_flux, wavelength[0], wavelength[-1], specific_photon_flux)
    out_lists = [X_all_round,x_calc_round]
    final_out = (out,out_lists,info_dict)
    return final_out