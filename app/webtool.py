#! /usr/bin/env python
#-*- coding: utf-8 -*-
#import os

#from __future__ import print_function
############################################## standard libs
import sys
import os
import time
from datetime import datetime
from copy import deepcopy
from math import degrees, radians, floor, ceil
##############################################  numpy & argparse
import numpy
import argparse
##############################################  pystrain
from pystrain.strain import *
from pystrain.geodesy.utm import *
from pystrain.iotools.iparser import *
import pystrain.grid
from pystrain.station import Station

############################################## ploting
from scipy.spatial import Delaunay
from math import sqrt, radians, sin, cos, atan2, pi, asin


from flask import Flask, flash, request, redirect, url_for, render_template
#from werkzeug import secure_filename

from math import sqrt, radians, sin, cos, atan2, pi, asin

app = Flask(__name__, template_folder="templates")

UPLOAD_FOLDER = '/uploads'
ALLOWED_EXTENSIONS = set(['txt', 'vel'])

#sta_lst_ell = []
#x_mean = 0
#y_mean = 0
#NoSta = 0
#input_filename = ""

@app.route('/StrainTool/website')
def website():
  return  render_template('website/index.html')

@app.route('/StrainWebTool')
def webtool():
  return render_template('webtool/index.html')

@app.route('/StrainWebTool/inputs')
def webtool_inputs():
  return render_template('webtool/tmpl_inputs.html')

@app.route('/StrainWebTool/parameters', methods=['GET', 'POST'])
def webtool_params():
  sta_lst_ell = []
  x_mean = 0
  y_mean = 0
  NoSta = 0
  input_filename = ""
  global NoSta
  global input_filename
  global x_mean
  global y_mean
  if request.method == 'POST':
    file = request.files['file']
    stations = []
    for line in file.readlines():
      stations.append(Station(line))
    #if len(stations):
      #return stations
    #else:
        #return None
  #sta_lst_ell = []
  for sta in stations:
     sta_lst_ell.append(sta)
  for idx, sta in enumerate(sta_lst_ell):
    sta_lst_ell[idx].lon = round(degrees(sta.lon), 3)
    sta_lst_ell[idx].lat = round(degrees(sta.lat), 3)
    sta_lst_ell[idx].vn = round(sta.vn*1.e3, 1)
    sta_lst_ell[idx].ve = round(sta.ve*1.e3, 1)
    sta_lst_ell[idx].sn = round(sta.sn*1.e3, 1)
    sta_lst_ell[idx].se = round(sta.se*1.e3, 1)
    sta_lst_ell[idx].rho = round(sta.rho*1.e3, 1)
    sta_lst_ell[idx].t = round(sta.t, 2)
  input_filename=file.filename
  NoSta = format(len(stations))
  x_mean, y_mean = barycenter(stations)
  #print('[DEBUG] Number of stations parsed: {}'.format(len(stations)))
  return render_template('webtool/tmpl_params.html', content = sta_lst_ell, input_file=file.filename, NoSta = NoSta, clon = x_mean, clat = y_mean)

@app.route('/StrainWebTool/results', methods=['GET', 'POST'])
def webtool_results():
  global NoSta
  global input_filename
  global x_mean
  global y_mean
  if request.method == 'POST':
    lonmin = request.form['lonmin']
    lonmax = request.form['lonmax']
    latmin = request.form['latmin']
    latmax = request.form['latmax']
    x_step = request.form['x_step']
    y_step = request.form['y_step']
    par_wt = request.form['par_wt']
    par_dmin = request.form['par_dmin']
    par_dmax = request.form['par_dmax']
    par_dstep = request.form['par_dstep']
  else:
    lonmin = 0
  return render_template('webtool/tmpl_results.html', input_file = input_filename, NoSta = NoSta,clon = x_mean, clat = y_mean, lonmin = lonmin, lonmax = lonmax, latmin = latmin, latmax = latmax, x_step = x_step, y_step = y_step, par_wt = par_wt, par_dmin = par_dmin, par_dmax = par_dmax, par_dstep = par_dstep)
