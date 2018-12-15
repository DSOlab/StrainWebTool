from flask import Flask, render_template
app = Flask(__name__, template_folder="templates")

@app.route('/StrainTool')
def website():
  return  render_template('website/index.html')

@app.route('/WebTool')
def webtool():
  return render_template('webtool/index.html')
