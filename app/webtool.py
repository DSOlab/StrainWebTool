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
  if request.method == 'POST':
    file = request.files['file']
    stations = []
    for line in file.readlines():
      stations.append(Station(line))
    #if len(stations):
      #return stations
    #else:
        #return None
  sta_lst_ell = []
  for sta in stations:
     sta_lst_ell.append(sta)
  for idx, sta in enumerate(sta_lst_ell):
    sta_lst_ell[idx].lon = degrees(sta.lon)
    sta_lst_ell[idx].lat = degrees(sta.lat)
    sta_lst_ell[idx].vn = sta.vn*1.e3
    sta_lst_ell[idx].ve = sta.ve*1.e3
    sta_lst_ell[idx].sn = sta.sn*1.e3
    sta_lst_ell[idx].se = sta.se*1.e3
    sta_lst_ell[idx].rho = sta.rho*1.e3
  NoSta = format(len(stations))
  x_mean, y_mean = barycenter(stations)
  #print('[DEBUG] Number of stations parsed: {}'.format(len(stations)))
  return render_template('webtool/tmpl_params.html', content = sta_lst_ell, input_file=file.filename, NoSta = NoSta, clon = x_mean, clat = y_mean)


@app.route('/StrainWebTool/results')
def webtool_results():
  return render_template('webtool/tmpl_results.html')
