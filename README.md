# StrainWebTool
Online Portal for StrainTool software 

**Simple Guide**

**Install Flask-virtual inviroment**
<pre id="block-samp"><samp>
pip install Flask
**Debian, Ubuntu**
sudo apt-get install python-virtualenv
**CentOS, Fedora**
sudo yum install python-virtualenv
</samp></pre>

**Run app**
<pre id="block-samp"><samp>
**create virtual enviroment**
virtualenv flask

**activate venv**
. flask/bin/activate

export FLASK_APP=run.py
flask run
* Running on http://127.0.0.1:5000/

**OR**
flask run --host=0.0.0.0
* will run on your public IP
</samp></pre>


**Files tree**
<pre id="block-samp"><samp>
app (main application)
 |--> templates (html files)
 |--> static (images, style files, css)
 |--> website.py : app to serve main website for straintool
 |--> webtool.py " app to serve template for webtool.ONLY template
</samp></pre>
