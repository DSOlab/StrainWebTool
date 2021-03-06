import os
from flask import Flask, flash, request, redirect, url_for, render_template
from werkzeug import secure_filename

UPLOAD_FOLDER = '/uploads'
ALLOWED_EXTENSIONS = set(['txt', 'vel'])

app = Flask(__name__, template_folder="templates")
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

#@app.route('/upload_file', methods=['GET', 'POST'])
#def upload_file():
    #if request.method == 'POST':
        ## check if the post request has the file part
        #if 'file' not in request.files:
            #flash('No file part')
            #return redirect(request.url)
        #file = request.files['file']
        ## if user does not select file, browser also
        ## submit an empty part without filename
        #if file.filename == '':
            #flash('No selected file')
            #return redirect(request.url)
        #if file and allowed_file(file.filename):
            #filename = secure_filename(file.filename)
            #file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            #return redirect(url_for('uploaded_file',
                                    #filename=filename))
    #return render_template('upload/upload.html')

@app.route('/upload', methods=['GET', 'POST'])
def upload_file():
    return render_template("upload/upload.html")

@app.route('/result', methods=['GET', 'POST'])
def result():
    if request.method == 'POST':
        file = request.files['file']
        a = file.read()
        result = a
    return render_template("templ.html", content = result)


@app.route('/getfile')
def getfile():
    with open('tmp/CNRS_midas.vel') as f:
       content = f.read()
    return render_template("templ.html", content=content)


'''
@app.route('/upload')
def upload_file():
   return render_template('upload/upload.html')
	
@app.route('/uploader', methods = ['GET', 'POST'])
def upload_file():
   if request.method == 'POST':
      f = request.files['file']
      f.save(secure_filename(f.filename))
      return 'file uploaded successfully'
		
if __name__ == '__main__':
   app.run(debug = True)
'''
