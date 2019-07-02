# StrainWebTool
Online Portal for StrainTool software 

[![License MIT](http://img.shields.io/badge/license-MIT-brightgreen.svg)](https://github.com/DSOlab/StrainWebTool/blob/master/LICENSE)
[![](https://img.shields.io/github/release/DSOlab/StrainWebTool.svg)](https://github.com/DSOlab/StrainWebTool/releases/latest)
[![](https://img.shields.io/github/tag/DSOlab/StrainWebTool.svg)](https://github.com/DSOlab/StrainWebTool/tags) 
[![](https://img.shields.io/github/stars/DSOlab/StrainWebTool.svg)](https://github.com/DSOlab/StrainWebTool/stargazers)
[![](https://img.shields.io/github/forks/DSOlab/StrainWebTool.svg)](https://github.com/DSOlab/StrainWebTool/network)
[![](https://img.shields.io/github/issues/DSOlab/StrainWebTool.svg)](https://github.com/DSOlab/StrainWebTool/issues)


## General

StrainWebTool is a web application developed to estimate strain tensor parameters using StrainTool Software.  The development of the application was based on [Flask](http://flask.pocoo.org/)  microframework for Python. [Bootstrap](https://getbootstrap.com/) open source toolkit was used to enable a responsive web design and [Leaflet](https://leafletjs.com/) open-source JavaScript library for producing interactive maps.

## Structure and background

The application consists of three basic parts:
1. webtool.py: the main python source code file, that includes all necessary functions utilizing StrainTool software and enables the  building of HTML templates.
2. static files: a folder including static files such as headers, footers, images, javascript source code files.
3. templates: a folder including the main templates for the application.

Three different HTML templates have been formulated to implement the application:
- tmpl_inputs.html is the template where the user uploads their input files.
- tmpl_params.html is the template where the user chooses the parameters that StrainTool will later on use to estimate strain tensor parameters
- tmpl_results.html is the template where the results are presented. The user can download result files and see the results visualized on an interactive map.
Each template consists of three basic columns. The first column contains all input forms for the parameters needed to estimate strain tensors. The second column contains the plot tools and options to generate GMT maps;  these however are not active in the current beta version. The third column is where the interactive map is placed and the results are visualized

## Installation


**Install StrainTool - Flask - virtual enviroment**
<pre id="block-samp"><samp>
- Install **StrainTool** Software with all prerequisites needed

**Install Flask**
$> pip install Flask

**Install Virtual Enviroment**
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
* Running on http://127.0.0.1:5000/inputs

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


## User Guidlines



## Active Map visualize strain tensor results


## Contributing

1. Create an issue and describe your idea
2. [Fork it](https://github.com/DSOlab/StrainTool/network#fork-destination-box)
3. Create your feature branch (`git checkout -b my-new-idea`)
4. Commit your changes (`git commit -am 'Add some feature'`)
5. Publish the branch (`git push origin my-new-idea`)
6. Create a new Pull Request
7. Profit! :white_check_mark:


## License
The work is licensed under [MIT-license](LICENSE)


## Authors & Bug Reports
**Dimitrios G. Anastasiou**
> Dr. Rural & Surveying Engineer | Dionysos Satellite Observatory - NTUA | dganastasiou@gmail.com

**Xanthos Papanikolaou**
> Rural & Surveying Engineer | Dionysos Satellite Observatory - NTUA | [xanthos@mail.ntua.gr](mailto:xanthos@mail.ntua.gr)

**Dr. Athanassios Ganas**
> Research Director | Institute of Geodynamics | National Observatory of Athens | [aganas@gein.noa.gr](mailto:aganas@gein.noa.gr)

**Prof. Demitris Paradissis**
> Professor NTUA |  Dionysos Satellite Observatory - NTUA | [dempar@central.ntua.gr](dempar@central.ntua.gr)

## ChangeLog

The history of releases can be viewed at [ChangeLog](.github/CHANGELOG.md)

## Acknowlegments
**EPOS IP - EPOS Implementation Phase**

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement N° 676564

Disclaimer: the content of this website reflects only the author’s view and the Commission is not responsible for any use that may be made of the information it contains.

## References
* Anastasiou D., Ganas A., Legrand J., Bruyninx C., Papanikolaou X., Tsironi V. and Kapetanidis V. (2019). Tectonic strain distribution over Europe from EPN data. EGU General Assembly 2019, Geophysical Research Abstracts, Vol. 21, EGU2019-17744-1 [Abstract](https://meetingorganizer.copernicus.org/EGU2019/EGU2019-17744-1.pdf)

* Shen, Z.-K., M. Wang, Y. Zeng, and F. Wang, (2015), Strain determination using spatially discrete geodetic data, Bull. Seismol. Soc. Am., 105(4), 2117-2127, doi: 10.1785/0120140247

* Veis, G., Billiris, H., Nakos, B., and Paradissis, D. (1992), Tectonic strain in greece from geodetic measurements, C.R.Acad.Sci.Athens, 67:129--166

* Python Software Foundation. Python Language Reference, version 2.7. Available at http://www.python.org

* [The Generic Mapping Tools - GMT](http://gmt.soest.hawaii.edu/)









