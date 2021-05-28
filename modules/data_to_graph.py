# -*- coding: utf-8 -*-
# =============================================================================
# Created By  : Benjamin Wriedt
# Copyright   : Copyright 2021, Institute of Chemical Engineering, Prof. Dr. Dirk Ziegenbalg, Ulm University'
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
from lmfit_single import actin_feox
import fnmatch
import os
import sys
import numpy as np
import os.path
import shutil

make_graph = False
csv_export =True
clear_others = True

## send standard strings to detailed_calculations.txt
save_stdout = sys.stdout # stdout is saved
fh = open("detailed_calculations.txt","w")
sys.stdout = fh

## advanced setting of variables and checks
qti.app.setLocale(QtCore.QLocale.c()) # make sure that the decimal separator is the dot character
y0 = 1e-5
comm=[]
ext_eff=[]
q_calc=[]
nit = 0
plotcom = []
plotphotef = []
ls_file = ""
quant_file = ""
comment = ""

# get input data from the tables, clear data with an irradiation time of 0
t1 = table("experimental_data")

index = 1
while index <= t1.numRows():
	irr_time_value = t1.text("irr. time/s", index)
	if irr_time_value == "0" or irr_time_value == "":
		del t1[index-1]
		print("removed col (no irr): ", index)
		continue
	index += 1


t2 = table("setup_properties")
P_emit = t2.cell(2, 3)
P_eff = t2.cell(2, 4)
P_rad = P_emit*P_eff
ls_file = t2.text("specific input", 5)
quant_file = t2.text("specific input", 6)

t3 = table("setup_properties-output")

# sort out different experimental data sets
if clear_others is True:
	i = 1
	first_annotation = t1.text("annotation", 1)
	while i <= t1.numRows():
		i_annotation = t1.text("annotation", i)
		if not first_annotation == i_annotation:
			del t1[i-1]
			print("removed col (others): ", i)
			continue
		i += 1
else:
	pass

# check if calculation is already done
if "photon flux/ mol/s" in t1.colNames():
	pass
else:
	q = t1.addColumn(2)
	q.setName("photon flux/ mol/s")
	t1.setColNumericFormat("photon flux/ mol/s", 2, 2)


if "ext. phot. eff./%" in t1.colNames():
	pass
else:
	q = t1.addColumn(7)
	q.setName("ext. phot. eff./%")
	t1.setColNumericFormat("ext. phot. eff./%", 1, 2)

# call of lmfit_single, pass of variables
print(t1.numRows())
for i in range (1, t1.numRows()+1):
		comment = t1.text("annotation", i)
		output = actin_feox(ls_file = ls_file, quant_file=quant_file, P_rad=P_rad, V_inj=t1.cell("irr. Vol/mL", i), X_all=[t1.cell("conversion/-", i)], t_all=np.array([t1.cell("irr. time/s", i)]), file_out=False, comment=comment, shell_out=False, plot_parity=False, c_act=t2.cell(2, 1), l=t2.cell(2, 2))
		t1.setCell("photon flux/ mol/s", i, output[0][1])	
		t1.setCell("ext. phot. eff./%", i, output[0][2])
		epe_percent = t1.cell("ext. phot. eff./%", i) * 100
		t1.setCell("ext. phot. eff./%", i, epe_percent)
		nit = nit+1
		print(output)

# print output in setup_properties-output
t3.setText("quantity", 1, "emitted photon flux/ mol/s")
t3.setCell("value", 1, output[0][4])
t3.setText("quantity", 2, "minimum wavelength / nm")
t3.setCell("value", 2, output[0][5])
t3.setText("quantity", 3, "maximum wavelength / nm")
t3.setCell("value", 3, output[0][6])
t3.setText("quantity", 4, "specific photon flux / mol/s")
t3.setCell("value", 4, output[0][7])
t3.notifyChanges()

## create a new graph without title
if make_graph is True:
	g =newGraph("photon flux over irradiation time")
	l = g.activeLayer()
	l.showGrid()
	l.setAntialiasing(True)
	l.enableAxis(Layer.Bottom, True)
	l.enableAxis(Layer.Left, True)
	l.enableAxis(Layer.Right, False)
	l.enableAxis(Layer.Top, False)
	l.setCanvasFrame(0)
	l.removeTitle()

	## x-axis setup 
	l.setAxisTitle(Layer.Bottom, "time/s")
	l.setAxisTitle(Layer.Right, "")
	l.setAxisTitleFont(Layer.Bottom, QtGui.QFont("Arial", 10))
	l.setAxisFont(Layer.Bottom, QtGui.QFont("Arial", 10))
	l.setAxisTicksLength(Layer.Bottom, 1, 1, 3, 6) # arguments are (axis, majTicksType, minTicksType, minLength, majLength)

	## y-axis setup 
	l.setAxisTitle(Layer.Left, "photon flux/mol/s")
	l.setAxisTitle(Layer.Top, "")
	l.setAxisTitleFont(Layer.Left, QtGui.QFont("Arial", 10))
	l.setAxisFont(Layer.Left, QtGui.QFont("Arial", 10))
	l.setAxisNumericFormat(Layer.Left, 3, 1) # arguments are (axis, format, precision, formula)
	l.setAxisTicksLength(Layer.Left, 1, 1, 3, 6) # arguments are (axis, majTicksType, minTicksType, minLength, majLength)

	## input of data from experimental_data
	c = g.activeLayer()
	c = l.insertCurve (table("experimental_data"), "experimental_data_irr. time/s","experimental_data_photon flux/ mol/s")
	s = c.symbol()
	s.setSize(QtCore.QSize(7, 7))# or s.setSize(7)
	s.setBrush(QtGui.QBrush(Qt.darkYellow))
	s.setPen(QtGui.QPen(Qt.darkYellow, 3))
	s.setStyle(PlotSymbol.XCross)
	l.replot() # redraw the plot layer object


	## general axis settings - must be at the end of the code to not being overwwritten
	l.drawAxesBackbones(True)
	l.setAxesLinewidth(2)

print("Emitted photon flux of the lightsource according to the available data: " + str(output[0][4]))

## end writing print command to detailed_calculations.txt
sys.stdout = save_stdout # return to normal:
fh.close()

## .csv export
if csv_export is True:
	csv_name = t1.text("annotation", 1)
	with open("data/" + str(csv_name) + ".csv", "w") as fout:
		content = []
		for i in range (1, t1.numRows()):
			time = t1.text("irr. time/s", i)
			flux = t1.text("photon flux/ mol/s", i)
			fout.write(str(time) + "\t" + str(flux) + "\n")