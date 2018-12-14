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
**activate venv**
. flask/bin/activate

export FLASK_APP=run.py
flask run
* Running on http://127.0.0.1:5000/

**OR**
flask run --host=0.0.0.0
* will run on your public IP
</samp></pre>
