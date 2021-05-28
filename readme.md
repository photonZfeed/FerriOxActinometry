# =============================================================================
# Created By  : Benjamin Wriedt
# Copyright   : Copyright 2021, Institute of Chemical Engineering, Prof. Dr. Dirk Ziegenbalg, Ulm University
# License     : GNU LGPL
# ============================================================================

1. Origin and Scope

	This set of programming code and data is part of the evaluation procedure of actinometry developed 
	at the Institute of Chemical Engineering, University Ulm. It correlates with the recommendations for
	the experimental procedure described in the full paper "Common pitfalls in chemical actinometry", 
	Journal of Flow Chemistry (2020) 10:295–306, DOI:10.1007/s41981-019-00072-7.

2. Installation

	The original programming codes were aligned to Qtiplot 0.9.9.13-win32 and Python 2.7.14 used by QtiPlot.
	The required programs are available at the following addresses:
	-	QtiPlot full version 32-bit: https://www.qtiplot.com/download.html or the 32-bit demo: https://www.qtiplot.com/demo.html
	-   Python 2.7.14: https://www.python.org/downloads/release/python-2714/ (works independently of newer versions)
	
	To make sure that all programs run properly and to avoid registry problems, please follow the installation order:
	-	install Qtiplot
	-	install Python 2.7.14
	-	install missing packages for Python 2.7 (get-pip.py is recommended as it is already installed with Python)

3. Included Modules
	
	3.1 photon_flux_calculation.qti
	 
	 This main file is used as a GUI to enter experimental data (table 1) and experimental conditions (table 2).
	 All dependencies must be provided in the same folder as photon_flux_calculation.qti. 
	 Some examples are located in the data subfolder beforehand.
	 The program photon_flux_calculation.py must be opened in this template to get the correct input data.
	 After the excecution of photon_flux_calculation.py, the result output is automatically entered in table 1.
	
	3.2 photon_flux_calculation.py

	 The script must be opened and executed within QtiPlot, compatibility is given for 0.9.9.13 32-bit version.
	 Based on the experimental data input via the graphical interface of QtiPlot, the available photon flux (q_abs) in 
	 the evaluated volume is calculated by passing the variables to lmfit_single.py. Overall and partial emitted photon range 
	 flux in the evaluated wavelength are calculated while the latter is used to determine the external photonic efficiency.

	 Dependencies:
	   lmfit_single.py
	   ls_file (to be specified)
	   quant_file (to be specified)
	   Input from QtiPlot tables (see Input variables)

	 Program options (boolean):
	   make_graph: creates a scaled graph with every run
	   csv_export: creates a .csv in \data with photon flux and ext. phot. eff. (separator: \t)
	   clear_others: Erases measurements with a different comment (assumingly from another data set)

	 Input variables:
	   P_emit: electrical energy need of the light source
	   P_eff: efficiency of the light source
	   
	 Output variables:
	   P_rad: emitted radiant power of the light source

	3.3 lmfit_single.py

	 This python script is called while executing the evaluation script within QtiPlot while all required data is exchanged. 
	 It calculates the available photon flux by solving the appropriate ODE. Spectra are imported and the wavelength intervall 
	 is interpolated to 1 nm for all sets. Parts of the spectrum with no corresponding data are cut off. Lightsource
	 data is converted from electric to photonic radiant power to calculate the external photonic efficiency.
	 
	 Called from:
	   data_to_graph.py
	   ls_file (to be specified)
	   quant_file (to be specified)
	   Input from QtiPlot tables (see Input variables)

	 Program options (boolean):
	   make_graph: creates a scaled graph with every run
	   csv_export: creates a .csv in \data with photon flux and ext. phot. eff. (separator: \t)
	   clear_others: Erases measurements with a different comment (assumingly from another data set)

	 Input variables:
	   y0: Start value for variation of the available photon flux
	   ls_file: filename of the .csv for the light source emission and reactor material transmission spectrum
	   quant_file: filename of the chi_*.csv for the quantum yield spectrum, sets kappa_file as kappa_*.csv
	   P_rad: electric radiant power emitted from the light source(s)
	   c_act: starting concentration of the used actinometer
	   l: optical path length (constant)
	   [V_inj]: list of irradiated reactor volume
	   [X_all]: list of experimentally determined conversion
	   [t_all]: list of experimental irradiation time
	   
	 Args:
	   wavelength_ls, arb_intensity: spectral data of the arbitrary light source intensity from ls_file
	   trans_FEP: spectral data of the reactor material transmission (based on wavelength_ls) from ls_file
	   wavelength_quant, quantum_yield: spectral data of the quantum yield from quant_file
	   wavelength_kappa, absorption_coefficient: spectral data of the absorption coefficient from kappa_file
	   spectral_photon_flux: absolute emitted photon flux according to ls_file
	   fraction_of_spectral_photon_flux: relative emitted photon flux in the evaluated spectrum according to ls_file

	 Output variables:
	   final_out(out,out_lists,info_dict): results data table passed to QtiPlot
		   out [0,i]:
			   i=0: comment from table 1
			   i=1: available photon flux in the irradiated volume according to actinometry (q_abs)
			   i=2: external phononic efficiency in evaluated wavelenth range
			   i=3: standard deviation between experimental value and lmfit model
			   i=4: overall emitted photon flux of the light source
			   i=5: minimum wavelength of the evaluated interval
			   i=6: maximum wavelength of the evaluated interval
			   i=7: emitted photon flux of the light source in the evaluated interval
		   out_lists[X_all_round,x_calc_round]:
			   X_all_round: input conversion    
			   x_calc_round: fitted conversion
		   info_dict[ls_file, kappa_file]:
			   ls_file: filename of light source spectrum to be saved with the results (name specifies the light source)
			   kappa_file: filename of absorbance spectrum to be saved with the results (name specifies the actinometer)
