import os
from flask import Flask, flash, request, redirect, url_for, render_template
from werkzeug import secure_filename
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
    a = file.read()
    result = a
  return render_template('webtool/tmpl_params.html', content = result, input_file=file.filename)


@app.route('/StrainWebTool/results')
def webtool_results():
  return render_template('webtool/tmpl_results.html')
