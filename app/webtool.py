import os
from flask import Flask, flash, request, redirect, url_for, render_template
from werkzeug import secure_filename

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

def set_none(self):
        '''Set to None.

            Set all instance member values to None.
        '''
        self.name = None                                                        
        self.lon  = None
        self.lat  = None
        self.ve   = None
        self.vn   = None
        self.se   = None
        self.sn   = None
        self.rho  = None
        self.t    = None
        
@app.route('/StrainWebTool/parameters', methods=['GET', 'POST'])
def webtool_params():
  if request.method == 'POST':
    file = request.files['file']
    a = file.read()
    result = a
    self.set_none()

    #stations = []
    for line in file.readlines():
      l = input_line.split()
      try:
        self.name = l[0]
        self.lon  = radians(float(l[1]))
        self.lat  = radians(float(l[2]))
        self.ve   = float(l[3]) / 1e3
        self.vn   = float(l[4]) / 1e3
        self.se   = float(l[5]) / 1e3
        self.sn   = float(l[6]) / 1e3
        self.rho  = float(l[7]) / 1e3
        self.t    = float(l[8])
        print self.name
      except:
        print('[DEBUG] Invalid Station instance constrution.')
        print('[DEBUG] Input line \"{}\"'.format(input_line.strip()))
        raise RuntimeError
    #if len(stations):
      #return stations
    #else:
      #return None

  return render_template('webtool/tmpl_params.html', content = l, input_file=file.filename)


@app.route('/StrainWebTool/results')
def webtool_results():
  return render_template('webtool/tmpl_results.html')
