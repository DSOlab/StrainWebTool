#! /usr/bin/env python
#-*- coding: utf-8 -*-

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
#from math import sqrt, radians, sin, cos, atan2, pi, asin

############################################## Flask
from flask import Flask, flash, request, redirect, url_for, render_template, send_file
from flask_restful import reqparse

#from werkzeug import secure_filename

app = Flask(__name__, template_folder="templates")

UPLOAD_FOLDER = '/uploads'
ALLOWED_EXTENSIONS = set(['txt', 'vel'])

#sta_lst_ell = []
#x_mean = 0
#y_mean = 0
#NoSta = 0
#input_filename = ""
Version = 'StrainTensor.py Version: 1.0-rc4.1'


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
    #global x_mean
    #global y_mean
    global sta_list_ell
    sta_list_ell = []
    x_mean = 0
    y_mean = 0
    NoSta = 0
    input_filename = ""
    if request.method == 'POST':
        file = request.files['file']
        stations = []
        for line in file.readlines():
            stations.append(Station(line))
    for sta in stations:
        sta_list_ell.append(sta)

    sta_list_ell_tmpl = deepcopy(sta_list_ell)
    for idx, sta in enumerate(sta_list_ell_tmpl):
        sta_list_ell_tmpl[idx].lon = round(degrees(sta.lon), 3)
        sta_list_ell_tmpl[idx].lat = round(degrees(sta.lat), 3)
        sta_list_ell_tmpl[idx].vn = round(sta.vn*1.e3, 1)
        sta_list_ell_tmpl[idx].ve = round(sta.ve*1.e3, 1)
        sta_list_ell_tmpl[idx].sn = round(sta.sn*1.e3, 1)
        sta_list_ell_tmpl[idx].se = round(sta.se*1.e3, 1)
        sta_list_ell_tmpl[idx].rho = round(sta.rho*1.e3, 1)
        sta_list_ell_tmpl[idx].t = round(sta.t, 2)
    input_filename=file.filename
    NoSta = format(len(stations))
    #x_mean, y_mean = barycenter(stations)
    #x_mean = degrees(x_mean)
    #y_mean = degrees(y_mean)
    grd = pystrain.grid.generate_grid(sta_list_ell, 0.5 , 0.5, True)
    x_mean = (grd.x_min + grd.x_max)/2.
    y_mean = (grd.y_min + grd.y_max)/2.
    #print('[DEBUG] Number of stations parsed: {}'.format(len(stations)))
    return render_template('webtool/tmpl_params.html', content = sta_list_ell_tmpl, input_file=file.filename, NoSta = NoSta, clon = x_mean, clat = y_mean, grd = grd)

def cut_rectangle(xmin, xmax, ymin, ymax, sta_lst, sta_list_to_degrees=False):
    new_sta_lst = []
    for sta in sta_lst:
        if sta_list_to_degrees:
            slon = degrees(sta.lon)
            slat = degrees(sta.lat)
        else:
            slon = sta.lon
            slat = sta.lat
        if slon >= xmin and slon <= xmax and slat >= ymin and slat <= ymax:
            new_sta_lst.append(sta)
    return new_sta_lst

def write_station_info(sta_lst, filename='station_info.dat'):
    with open(filename, 'wb') as fout:
        #print('{:^10s} {:^10s} {:^10s} {:7s} {:7s} {:7s} {:7s}'.format('Station', 'Longtitude', 'Latitude', 'Ve', 'Vn', 'sVe', 'sVn'), file=fout)
        #print('{:^10s} {:^10s} {:^10s} {:7s} {:7s} {:7s} {:7s}'.format('', 'deg.', 'deg', 'mm/yr', 'mm/yr', 'mm/yr', 'mm/yr'), file=fout)
        fout.write('{:^10s} {:^10s} {:^10s} {:7s} {:7s} {:7s} {:7s} \n'.format('Station', 'Longtitude', 'Latitude', 'Ve', 'Vn', 'sVe', 'sVn'))
        fout.write('{:^10s} {:^10s} {:^10s} {:7s} {:7s} {:7s} {:7s} \n'.format('', 'deg.', 'deg', 'mm/yr', 'mm/yr', 'mm/yr', 'mm/yr'))
        for idx, sta in enumerate(sta_lst):
            #print('{:10s} {:+10.5f} {:10.5f} {:+7.2f} {:+7.2f} {:+7.3f} {:+7.3f}'.format(sta.name, degrees(sta.lon), degrees(sta.lat), sta.ve*1e03, sta.vn*1e03, sta.se*1e03, sta.sn*1e03), file=fout)
            #print('{:10s} {:+10.5f} {:10.5f} {:+7.2f} {:+7.2f} {:+7.3f} {:+7.3f}'.format(sta.name, degrees(sta.lon), degrees(sta.lat), sta.ve*1e03, sta.vn*1e03, sta.se*1e03, sta.sn*1e03))
            fout.write('{:10s} {:+10.5f} {:10.5f} {:+7.2f} {:+7.2f} {:+7.3f} {:+7.3f} \n'.format(sta.name, degrees(sta.lon), degrees(sta.lat), sta.ve*1e03, sta.vn*1e03, sta.se*1e03, sta.sn*1e03))
    return

def print_model_info(fout, cmd, clargs):
    fout.write('{:} \n'.format(Version))
    fout.write('Command used:\n\t{:}\n'.format(' '.join(cmd)))
    fout.write('Run at: {:}\n'.format(datetime.now().strftime('%c')))
    fout.write('Command line switches/options parsed:\n')
    for key in clargs:
        fout.write('\t{:20s} -> {:}\n'.format(key, clargs[key]))
    return

class get_strain_param:
    strain_param_names = ['lat', 'lon', 'vx', 'dvx', 'vy', 'dvy', 'w', 'dw', 'exx', 'dexx', 'exy', 'dexy', 'eyy', 'deyy', 'emax', 'demax', 'emin', 'demin', 'shr', 'dshr', 'azi', 'dazi', 'dilat', 'ddilat', 'secinv', 'dsecinv' ]

    def __init__(self, *args, **kargs):
        self.set_none()

        if len(args) is not 0:
            self.init_from_ascii_line(args[0])

        if len(kargs) is not 0:
            for key, val in kargs.items():
                if key in station_member_names:
                    setattr(self, key, val)
    
    
    
    def init_from_ascii_line(self, input_line):

        l = input_line.split()
        try:
            self.lat     = float(l[0])
            self.lon     = float(l[1])
            self.vx      = float(l[2])
            self.dvx     = float(l[3])
            self.vy      = float(l[4])
            self.dvy     = float(l[5])
            self.w       = float(l[6])
            self.dw      = float(l[7])
            self.exx     = float(l[8])
            self.dexx    = float(l[9])
            self.exy     = float(l[10])
            self.dexy    = float(l[11])
            self.eyy     = float(l[12])
            self.deyy    = float(l[13])
            self.emax    = float(l[14])
            self.demax   = float(l[15])
            self.emin    = float(l[16])
            self.demin   = float(l[17])
            self.shr     = float(l[18])
            self.dshr    = float(l[19])
            self.azi     = float(l[20])
            self.dazi    = float(l[21])
            self.dilat   = float(l[22])
            self.ddilat  = float(l[23])
            self.secinv  = float(l[24])
            self.dsecinv = float(l[25])
        except:
            print('[DEBUG] Invalid Station instance constrution.')
            print('[DEBUG] Input line \"{}\"'.format(input_line.strip()))
            #raise RuntimeError

    def set_none(self):
        self.lat     = None
        self.lon     = None
        self.vx      = None
        self.dvx     = None
        self.vy      = None
        self.dvy     = None
        self.w       = None
        self.dw      = None
        self.exx     = None
        self.dexx    = None
        self.exy     = None
        self.dexy    = None
        self.eyy     = None
        self.deyy    = None
        self.emax    = None
        self.demax   = None
        self.emin    = None
        self.demin   = None
        self.shr     = None
        self.dshr    = None
        self.azi     = None
        self.dazi    = None
        self.dilat   = None
        self.ddilat  = None
        self.secinv  = None
        self.dsecinv = None

@app.route('/StrainWebTool/results', methods=['GET', 'POST'])
def webtool_results():
    global NoSta
    global input_filename
    #global x_mean
    #global y_mean
    global sta_list_ell
    
    if request.method == 'POST':
        lonmin = 0#request.form['lonmin']
        lonmax = 0#request.form['lonmax']
        latmin = 0#request.form['latmin']
        latmax = 0#request.form['latmax']

        parser = reqparse.RequestParser()
        
        if request.form.get('shen'):
            parser.add_argument('shen',
                location='form',
                default='shen',
                dest='method',
                choices=['shen', 'veis'],
                required=False,
                help='Choose a method for strain estimation. If \'shen\' is passed in, the estimation will follow the algorithm described in Shen et al, 2015, using a weighted least squares approach. If \'veis\' is passed in, then the region is going to be split into delaneuy triangles and a strain estimated in each barycenter.')
        
        if request.form.get('veis'):
            parser.add_argument('veis',
                location='form',
                default='veis',
                dest='method',
                choices=['shen', 'veis'],
                required=False,
                help='Choose a method for strain estimation. If \'shen\' is passed in, the estimation will follow the algorithm described in Shen et al, 2015, using a weighted least squares approach. If \'veis\' is passed in, then the region is going to be split into delaneuy triangles and a strain estimated in each barycenter.')
        
        if request.form.get('barycenter'):
            parser.add_argument('barycenter',
                location='form',
                dest='one_tensor',
                action='store_true',
                help='Only estimate one strain tensor, at the region\'s barycentre.')
        
        parser.add_argument('region',
            location='form',
            default='18/23/32/43', #reqparse.SUPPRESS,
            dest='region',
            help='Specify a region; any station (in the input file) falling outside will be ommited. The region should be given as a rectangle, specifying min/max values in longtitude and latitude (using decimal degrees). E.g. \"[...] --region=21.0/23.5/36.0/38.5 [...]\"',
            required=False)
        
        parser.add_argument('x-grid-step',
            location='form',
            default=0.5,
            dest='x_grid_step',
            type=float,
            required=False,
            help='The x-axis grid step size in degrees. This option is only relevant if the program computes more than one strain tensors.')
        
        parser.add_argument('y-grid-step',
            location='form',
            default=0.5,
            dest='y_grid_step',
            type=float,
            required=False,
            help='The y-axis grid step size in degrees. This option is only relevant if the program computes more than one strain tensors.')
        
        parser.add_argument('Wt',
            location='form',
            default=24,
            dest='Wt',
            type=int,
            required=False,
            help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. Let W=Σ_i*G_i, the total reweighting coefficients of the data, and let Wt be the threshold of W. For a given Wt, the smoothing constant D is determined by Wd=Wt . It should be noted that W is a function of the interpolation coordinate, therefore for the same Wt assigned, D varies spatially based on the in situ data strength; that is, the denser the local data array is, the smaller is D, and vice versa.')
        
        parser.add_argument('dmin',
            location='form',
            default=1,
            dest='dmin',
            type=int,
            required=False,
            help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. This is the lower limit for searching for an optimal D-parameter value. Unit is km.')
        
        parser.add_argument('dmax',
            location='form',
            default=500,
            dest='dmax',
            type=int,
            required=False,
            help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. This is the upper limit for searching for an optimal d-param value. Unit is km.')
        
        parser.add_argument('dstep',
            location='form',
            default=2,
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

    args.verbose_mode = False
    if 'one_tensor' not in args:
        args.one_tensor = False
    if 'cut_outoflim_sta' not in args:
        args.cut_outoflim_sta = False
    if 'generate_stats' not in args:
        args.generate_stats = False
    if 'max_beta_angle' not in args:
        args.max_beta_angle = 180
    
    args.d_coef = None
    args.ltype = 'gaussian'
    dargs = vars(args)
    
    ##  Time the program (for opt/ing purpose only)
    #start_time = time.time()

    ## Verbose print (function only exists in verbose mode)
    #vprint = print if args.verbose_mode else lambda *a, **k: None
    
    ## If needed, open a file to write model info and statistics
    fstats = open('strain_stats.dat', 'w') if args.generate_stats else None
    if fstats: print_model_info(fstats, sys.argv, dargs)

    ##  If a region is passed in, resolve it.
    ##+ If cutting out-of-limits stations option is set, or method is veis, then 
    ##+ only keep the stations that fall within it.
    ##  The region coordinates (min/max pairs) should be given in decimal degrees.
    if 'region' in args:
        try:
            lonmin, lonmax, latmin, latmax = [ float(i) for i in args.region.split('/') ]
            if args.cut_outoflim_sta or args.method == 'veis':
                Napr = len(sta_list_ell)
                #  Note that we have to convert radians to degrees for station 
                #+ coordinates, hence 'sta_list_to_degrees=True'
                sta_list_ell = cut_rectangle(lonmin, lonmax, latmin, latmax, sta_list_ell, True)
                Npst = len(sta_list_ell)
                #vprint('[DEBUG] Stations filtered to fit input region: {:7.3f}/{:7.3f}/{:7.3f}/{:7.3f}'.format(lonmin, lonmax, latmin, latmax))
                #vprint('[DEBUG] {:4d} out of original {:4d} stations remain to be processed.'.format(Npst, Napr))
        except:
            ## TODO we should exit with error here
            print('[ERROR] Failed to parse region argument \"{:}\"'.format(args.region))
 
    ##  Filter out stations that are never going to be used. This is an opt!
    if 'region' in args and not args.method == 'veis' and not args.cut_outoflim_sta:
        #vprint('[DEBUG] Filtering stations based on their distance from region barycentre.')
        Napr = len(sta_list_ell)
        mean_lon, mean_lat = radians(lonmin+(lonmax-lonmin)/2e0), radians(latmin+(latmax-latmin)/2e0)
        # print('[DEBUG] Barycentre set to {:+10.5f}, {:10.5f}'.format(lonmin+(lonmax-lonmin)/2e0, latmin+(latmax-latmin)/2e0))
        bc =  Station(lon=mean_lon, lat=mean_lat)
        endpt = Station(lon=radians(lonmax), lat=radians(latmax))
        cutoffdis = abs(endpt.haversine_distance(bc)/1e3) # barycentre to endpoint (km)
        # print('[DEBUG] Distance from barycentre to endpoint is {:10.3f}km'.format(cutoffdis))
        d = 2e0*(args.d_coef if args.d_coef is not None else args.dmax)
        cutoffdis += d * (2.15e0 if args.ltype == 'gaussian' else 10e0) # in km
        #vprint('[DEBUG] Using cut-off distance {:10.3f}km'.format(cutoffdis))
        #for s in sta_list_ell:
        #    d = s.haversine_distance(bc)
        #    print('station {:} is {:7.1f}km away {:}'.format(s.name, d/1e3, 'Acc' if s.haversine_distance(bc)/1e3 <= cutoffdis else 'Rej'))
        sta_list_ell = [ s for s in sta_list_ell if s.haversine_distance(bc)/1e3 <= cutoffdis ]
        Npst = len(sta_list_ell)
        print('[DEBUG] {:4d} out of original {:4d} stations remain to be processed.'.format(Npst, Napr))
        Npst = len(sta_list_ell)
    
    ##  Make a new station list (copy of the original one), where all coordinates
    ##+ are in UTM. All points should belong to the same ZONE.
    ##  Note that station ellipsoidal coordinates are in radians while the cartesian
    ##+ coordinates are in meters.
    ##
    ##  TODO is this mean_lon the optimal?? or should it be the region's mean longtitude
    ##
    mean_lon = degrees(sum([ x.lon for x in sta_list_ell ]) / len(sta_list_ell))
    utm_zone = floor(mean_lon/6)+31
    utm_zone = utm_zone + int(utm_zone<=0)*60 - int(utm_zone>60)*60
    #vprint('[DEBUG] Mean longtitude is {} deg.; using Zone = {} for UTM'.format(mean_lon, utm_zone))
    sta_list_utm = deepcopy(sta_list_ell)
    for idx, sta in enumerate(sta_list_utm):
        N, E, Zone, lcm = ell2utm(sta.lat, sta.lon, Ellipsoid("wgs84"), utm_zone)
        sta_list_utm[idx].lon = E
        sta_list_utm[idx].lat = N
        assert Zone == utm_zone, "[ERROR] Invalid UTM Zone."
    #vprint('[DEBUG] Station list transformed to UTM.')
    
    ##  Open file to write Strain Tensor estimates; write the header
    fout = open('strain_info.dat', 'w')
    #vprint('[DEBUG] Strain info written in file: {}'.format('strain_info.dat'))
    fout.write('{:^9s} {:^9s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s}\n'.format('Latitude', 'Longtitude', 'vx+dvx', 'vy+dvy', 'w+dw', 'exx+dexx', 'exy+dexy', 'eyy+deyy', 'emax+demax', 'emin+demin', 'shr+dshr', 'azi+dazi', 'dilat+ddilat', 'sec. invariant'))
    fout.write('{:^9s} {:^9s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s}\n'.format('deg', 'deg', 'mm/yr', 'mm/yr', 'deg/Myr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'deg.', 'nstrain/yr', 'nstrain/yr'))
    
    ##  Compute only one Strain Tensor, at the region's barycenter; then exit.
    if args.one_tensor:
        print('[DEBUG] Estimating Strain Tensor at region\'s barycentre.')
        if args.method == 'shen':
            sstr = ShenStrain(0e0, 0e0, sta_list_utm, **dargs)
        else:
            sstr = ShenStrain(0e0, 0e0, sta_list_utm, weighting_function='equal_weights')
        sstr.set_to_barycenter()
        sstr.estimate()
        sstr.print_details(fout, utm_zone)
        fout.close()
        write_station_info(sta_list_ell)
        NoTensors = 1
        #print('[DEBUG] Total running time: {:10.2f} sec.'.format((time.time() - start_time)))      
        #sys.exit(0)
    
    # strain_list = [] Probably we do not need to keep the tensors ...
    if args.method == 'shen' and not args.one_tensor:  ## Going for Shen algorithm ...
        ##  Construct the grid, in ellipsoidal coordinates --degrees--. If a region
        ##+ is not passed in, the grid.generate_grid will transform lon/lat pairs 
        ##+ to degrees and produce a grid from extracting min/max crds from the
        ##+ station list.
        if 'region' in args:
            grd = pystrain.grid.Grid(lonmin, lonmax, args.x_grid_step, latmin, latmax, args.y_grid_step)
        else:
            grd = pystrain.grid.generate_grid(sta_list_ell, args.x_grid_step, args.y_grid_step, True)
        print('[DEBUG] Grid Information:')
        print('[DEBUG]\tLongtitude : from {} to {} with step {} (deg)'.format(grd.x_min, grd.x_max, grd.x_step))
        print('[DEBUG]\tLatitude   : from {} to {} with step {} (deg)'.format(grd.y_min, grd.y_max, grd.y_step))
        print('[DEBUG] Number of Strain Tensors to be estimated: {}'.format(grd.xpts*grd.ypts))
        if fstats: fstats.write('{:^10s} {:^10s} {:^10s} {:^12s} {:^12s} {:^12s}\n'.format('Longtitude','Latitude','# stations', 'D (optimal)','CutOff dis.', 'Sigma'))
        if fstats: fstats.write('{:^10s} {:^10s} {:^10s} {:^12s} {:^12s} {:^12s}\n'.format('deg.','deg.','#', 'Km','#', '/'))
        #vprint('[DEBUG] Estimating strain tensor for each cell center:')
        ##  Iterate through the grid (on each cell center). Grid returns cell-centre
        ##+ coordinates in lon/lat pairs, in degrees!
        node_nr, nodes_estim = 0, 0
        for x, y in grd:
            clat, clon =  radians(y), radians(x)
            N, E, ZN, _ = ell2utm(clat, clon, Ellipsoid("wgs84"), utm_zone)
            assert ZN == utm_zone
            #vprint('[DEBUG] Grid point at {:+8.4f}, {:8.4f} or E={:}, N={:}'.format(x, y, E, N))
            #print('[DEBUG] {:5d}/{:7d}'.format(node_nr+1, grd.xpts*grd.ypts), end="\r")
            ## Construct the Strain instance, with all args (from input)
            sstr = ShenStrain(E, N, sta_list_utm, **dargs)
            ## check azimouth coverage (aka max β angle)
            if degrees(max(sstr.beta_angles())) <= args.max_beta_angle:
                try:
                    sstr.estimate()
                    #vprint('[DEBUG] Computed tensor at {:+8.4f} {:+8.4f} for node {:3d}/{:3d}'.format(x, y, node_nr+1, grd.xpts*grd.ypts))
                    sstr.print_details(fout, utm_zone)
                    if fstats: fstats.write('{:+9.4f} {:+10.4f} {:6d} {:14.2f} {:10.2f} {:12.3f}\n'.format(x,y,len(sstr.__stalst__), sstr.__options__['d_coef'],sstr.__options__['cutoff_dis'], sstr.__sigma0__))
                    # strain_list.append(sstr)
                    nodes_estim += 1
                except RuntimeError:
                    print('[DEBUG] Too few observations to estimate strain at {:+8.4f}, {:8.4f}. Point skipped.'.format(x,y))
                except ArithmeticError:
                    print('[DEBUG] Failed to compute parameter VcV matrix for strain at {:+8.4f}, {:8.4f}. Point skipped'.format(x,y))
            else:
                print('[DEBUG] Skipping computation at {:+8.4f},{:8.4f} because of limited coverage (max_beta= {:6.2f}deg.)'.format(x, y, degrees(max(sstr.beta_angles()))))
            node_nr += 1
        print('[DEBUG] Estimated Strain Tensors for {} out of {} nodes'.format(nodes_estim, node_nr))
        NoTensors = nodes_estim
    elif args.method == 'veis' and not args.one_tensor:
        ## Open file to write delaunay triangles.
        print('[DEBUG] Estimating Strain Tensors at the barycentre of Delaunay triangles')
        dlnout = open('delaunay_info.dat', 'w')
        points = numpy.array([ [sta.lon, sta.lat] for sta in sta_list_utm ])
        tri = Delaunay(points)
        print('[DEBUG] Number of Delaunay triangles: {}'.format(len(tri.simplices)))
        NoTensors = len(tri.simplices)
        for idx, trng in enumerate(tri.simplices):
            #print('[DEBUG] {:5d}/{:7d}'.format(idx+1, len(tri.simplices)), end="\r")
            ## triangle barycentre
            cx = (sta_list_utm[trng[0]].lon + sta_list_utm[trng[1]].lon + sta_list_utm[trng[2]].lon)/3e0
            cy = (sta_list_utm[trng[0]].lat + sta_list_utm[trng[1]].lat + sta_list_utm[trng[2]].lat)/3e0
            ##  Construct a strain instance, at the triangle's barycentre, with only
            ##+ 3 points (in UTM) and equal_weights weighting scheme.
            sstr = ShenStrain(cx, cy, [sta_list_utm[trng[0]], sta_list_utm[trng[1]], sta_list_utm[trng[2]]], weighting_function='equal_weights')
            sstr.estimate()
            sstr.print_details(fout, utm_zone)
            ## Print the triangle in the corresponding file (ellipsoidal crd, degrees)
            dlnout.write('> {:}, {:}, {:}\n'.format(sta_list_utm[trng[0]].name, sta_list_utm[trng[1]].name, sta_list_utm[trng[2]].name))
            dlnout.write('{:+8.5f} {:8.5f}\n{:+8.5f} {:8.5f}\n{:+8.5f} {:8.5f}\n{:+8.5f} {:8.5f}\n'.format(*[ degrees(x) for x in [sta_list_ell[trng[0]].lon, sta_list_ell[trng[0]].lat, sta_list_ell[trng[1]].lon, sta_list_ell[trng[1]].lat, sta_list_ell[trng[2]].lon, sta_list_ell[trng[2]].lat, sta_list_ell[trng[0]].lon, sta_list_ell[trng[0]].lat]]))
            # strain_list.append(sstr)
        dlnout.close()

    fout.close()
    write_station_info(sta_list_ell)
    
    sta_list_ell_tmpl = deepcopy(sta_list_ell)
    for idx, sta in enumerate(sta_list_ell_tmpl):
        sta_list_ell_tmpl[idx].lon = round(degrees(sta.lon), 3)
        sta_list_ell_tmpl[idx].lat = round(degrees(sta.lat), 3)
        sta_list_ell_tmpl[idx].vn = round(sta.vn*1.e3, 1)
        sta_list_ell_tmpl[idx].ve = round(sta.ve*1.e3, 1)
        sta_list_ell_tmpl[idx].sn = round(sta.sn*1.e3, 1)
        sta_list_ell_tmpl[idx].se = round(sta.se*1.e3, 1)
        sta_list_ell_tmpl[idx].rho = round(sta.rho*1.e3, 1)
        sta_list_ell_tmpl[idx].t = round(sta.t, 2)
    #print('[DEBUG] Total running time: {:10.2f} sec.'.format((time.time() - start_time)))
    grd_tmpl = pystrain.grid.generate_grid(sta_list_ell, 0.5 , 0.5, True)
    x_mean = (grd_tmpl.x_min + grd_tmpl.x_max)/2.
    y_mean = (grd_tmpl.y_min + grd_tmpl.y_max)/2.
    
    
    file = open('strain_info.dat', 'r')
    strain = []
    for line in file.readlines()[2:]:
        #print(line)
        strain.append(get_strain_param(line))
    strain_info = []
    for sta in strain:
        strain_info.append(sta)
    #print(strain_info)
        
    return render_template('webtool/tmpl_results.html', input_file = input_filename, NoSta = Npst,clon = x_mean, clat = y_mean, args = args, lonmin = lonmin, lonmax = lonmax, latmin = latmin, latmax = latmax, x_step = args.x_grid_step, y_step = args.y_grid_step, NoTensors = NoTensors, content = sta_list_ell_tmpl, strinfo = sstr, grd = grd_tmpl, strain_info = strain_info )

@app.route('/StrainWebTool/outputs/<filename>', methods=['GET', 'POST'])
def dowloadfile(filename):
    try:
        #Boto3 downloading the file file.csv here
        return send_file(filename, attachment_filename=filename)
    except Exception as ermsg:
        print(ermsg)
        return render_template('#', ermsg=ermsg)
