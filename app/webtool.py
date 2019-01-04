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
from flask_restful import reqparse

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
  global NoSta
  global input_filename
  global x_mean
  global y_mean
  sta_lst_ell = []
  x_mean = 0
  y_mean = 0
  NoSta = 0
  input_filename = ""
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
    #x_step = request.form['x-grid-step']
    #y_step = request.form['y_step']
    #par_wt = request.form['par_wt']
    #par_dmin = request.form['par_dmin']
    #par_dmax = request.form['par_dmax']
    #par_dstep = request.form['par_dstep']

    parser = reqparse.RequestParser()
    
    if request.form.get('shen'):
      parser.add_argument('shen',
        location='form',
        default='shen',
        #metavar='METHOD',
        dest='method',
        choices=['shen', 'veis'],
        required=False,
        help='Choose a method for strain estimation. If \'shen\' is passed in, the estimation will follow the algorithm described in Shen et al, 2015, using a weighted least squares approach. If \'veis\' is passed in, then the region is going to be split into delaneuy triangles and a strain estimated in each barycenter.')
    
    if request.form.get('veis'):
      parser.add_argument('veis',
        location='form',
        default='veis',
        #metavar='METHOD',
        dest='method',
        choices=['shen', 'veis'],
        required=False,
        help='Choose a method for strain estimation. If \'shen\' is passed in, the estimation will follow the algorithm described in Shen et al, 2015, using a weighted least squares approach. If \'veis\' is passed in, then the region is going to be split into delaneuy triangles and a strain estimated in each barycenter.')
          
    parser.add_argument('x-grid-step',
      location='form',
      default=0.5,
      #metavar='X_GRID_STEP',
      dest='x_grid_step',
      type=float,
      required=False,
      help='The x-axis grid step size in degrees. This option is only relevant if the program computes more than one strain tensors.')
    
    parser.add_argument('y-grid-step',
      location='form',
      default=0.5,
      #metavar='X_GRID_STEP',
      dest='y_grid_step',
      type=float,
      required=False,
      help='The y-axis grid step size in degrees. This option is only relevant if the program computes more than one strain tensors.')
    
    parser.add_argument('Wt',
      location='form',
      default=24,
      #metavar='Wt',
      dest='Wt',
      type=int,
      required=False,
      help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. Let W=Σ_i*G_i, the total reweighting coefficients of the data, and let Wt be the threshold of W. For a given Wt, the smoothing constant D is determined by Wd=Wt . It should be noted that W is a function of the interpolation coordinate, therefore for the same Wt assigned, D varies spatially based on the in situ data strength; that is, the denser the local data array is, the smaller is D, and vice versa.')
    
    parser.add_argument('dmin',
      location='form',
      default=1,
      #metavar='D_MIN',
      dest='dmin',
      type=int,
      required=False,
      help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. This is the lower limit for searching for an optimal D-parameter value. Unit is km.')
    
    parser.add_argument('dmax',
      location='form',
      default=500,
      #metavar='D_MAX',
      dest='dmax',
      type=int,
      required=False,
      help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. This is the upper limit for searching for an optimal d-param value. Unit is km.')
    
    parser.add_argument('dstep',
      location='form',
      default=2,
      #metavar='D_STEP',
      dest='dstep',
      type=int,
      required=False,
      help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. This is the step size for searching for an optimal d-param value. Unit is km.')

    if request.form.get('cut-excess-stations'):
      parser.add_argument('cut-excess-stations',
        location='form',
        dest='cut_outoflim_sta',
        help='This option is only considered if the \'-r\' option is set. If this this option is enabled, then any station (from the input file) outside the region limit (passed in via the \'-r\' option) is not considered in the strain estimation.',
        action='store_true')

    if request.form.get('generate-statistics'):
      parser.add_argument('generate-statistics',
        location='form',
        dest='generate_stats',
        help='Only relevant when \'--mehod=shen\' and \'--barycenter\' is not set. This option will create an output file, named \'strain_stats.dat\', where estimation info and statistics will be written.',
        action='store_true')

    
  else:
    lonmin = 0
    
  ##  Parse command line arguments and stack them in a dictionary
  args  = parser.parse_args()
  print(args)
  #dargs = vars(args)
  #print(dargs)
  return render_template('webtool/tmpl_results.html', input_file = input_filename, NoSta = NoSta,clon = x_mean, clat = y_mean, method = args.method, lonmin = lonmin, lonmax = lonmax, latmin = latmin, latmax = latmax, x_step = args.x_grid_step, y_step = args.y_grid_step, par_wt = args.Wt, par_dmin = args.dmin, par_dmax = args.dmax, par_dstep = args.dstep)
