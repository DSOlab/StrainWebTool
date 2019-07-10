# StrainWebTool
Online Portal for [StrainTool](https://github.com/DSOlab/StrainTool) software 

- estimate strain tensor parameters using different methods.

- Visualize results on interactive map

- Download output files for further processing

 Check out the [demo](http://83.212.103.160/StrainWebTool/inputs) !!

[![License MIT](http://img.shields.io/badge/license-MIT-brightgreen.svg)](https://github.com/DSOlab/StrainWebTool/blob/master/LICENSE)
[![](https://img.shields.io/github/release/DSOlab/StrainWebTool.svg)](https://github.com/DSOlab/StrainWebTool/releases/latest)
[![](https://img.shields.io/github/tag/DSOlab/StrainWebTool.svg)](https://github.com/DSOlab/StrainWebTool/tags) 
[![](https://img.shields.io/github/stars/DSOlab/StrainWebTool.svg)](https://github.com/DSOlab/StrainWebTool/stargazers)
[![](https://img.shields.io/github/forks/DSOlab/StrainWebTool.svg)](https://github.com/DSOlab/StrainWebTool/network)
[![](https://img.shields.io/github/issues/DSOlab/StrainWebTool.svg)](https://github.com/DSOlab/StrainWebTool/issues)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3308499.svg)](https://doi.org/10.5281/zenodo.3308499)

## General

StrainWebTool is a web application developed to estimate strain tensor parameters using StrainTool Software.  The development of the application was based on [Flask](http://flask.pocoo.org/)  microframework for Python. [Bootstrap](https://getbootstrap.com/) open source toolkit was used to enable a responsive web design and [Leaflet](https://leafletjs.com/) open-source JavaScript library for producing interactive maps.



## Structure and background

The application consists of three basic parts:
1. `webtool.py`: the main python source code file, that includes all necessary functions utilizing StrainTool software and enables the  building of HTML templates.
2. static files: a folder including static files such as headers, footers, images, javascript source code files.
3. templates: a folder including the main templates for the application.

Three different HTML templates have been formulated to implement the application:
- `tmpl_inputs.html` is the template where the user uploads their input files.
- `tmpl_params.html` is the template where the user chooses the parameters that StrainTool will later on use to estimate strain tensor parameters
- `tmpl_results.html` is the template where the results are presented. The user can download result files and see the results visualized on an interactive map.
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

### Input Files

To perform the computations, StrainWebTool needs an input file, that holds input data. Usually, this implies a list of GPS/GNSS stations with their ellipsoidal coordinates (aka longitude and latitude) and their respective tectonic velocities (usually estimated using position time-series) along with the corresponding standard deviation values. The format of these files, should follow the convention:


<pre id="block-samp" <samp="">        station-name  longtitude   latitude   Ve     Vn    SigmaVe  SigmaVn  Sne  time-span 
         string           deg.        deg.   mm/yr  mm/yr   mm/yr    mm/yr    /   dec. years</pre>

Station coordinates are provided in longitude/latitude pairs in decimal degrees. Velocities and velocity standard deviations are provided in mm per years (mm/yr). Sne is the correlation coefficient between East and North velocity components and time-span is the total time span of the station timeseries in decimal degrees. Note that at his point the last two columns (aka Sne and time-span) are not used, so they could have random values.
There are no strict formatting rules on how the individual elements should be printed (i.e. how many fields, decimal places, etc.). The only condition is that fields are separated by whitespace(s). 
Note that the input file format is identical to what is used in StrainTool Software (Anastasiou et al., 2019 ); users can browse its dedicated web page (https://dsolab.github.io/StrainTool/) for a more detailed description.


### Option and parameters

After the uploading of the input-file, all the options for the estimation of strain tensor are unlocked.

The first part is  the selection of the method for strain estimation. If 'shen' is passed in, the estimation will follow the algorithm described in Shen et al, 2015, using a weighted least squares approach. If 'veis' is passed in, then the region is going to be split into delaneuy triangles and a strain estimated in each barycenter. Default is 'shen'. If  ‘One Tensor’ checked, then only one strain tensor will be estimated, at the region’s barycentre.
In the second part, user specifies the region as a rectangle and x-axis/y-axis grid steps. Any station falling outside this region will be omitted.
In the third part, the user selects the interpolation model parameters for ‘shen’ method. The options are:

- Wt: Let W=Σi*Gi, the total reweighting coefficients of the data, and let Wt be the threshold of W. For a given Wt, the smoothing constant D is determined by Wd=Wt . It should be noted that W is a function of the interpolation coordinate, therefore for the same Wt assigned, D varies spatially based on the in situ data strength; that is, the denser the local data array is, the smaller is D, and vice versa. Default is Wt=24.
- D min: This is the lower limit for searching for an optimal d-param value. Unit is km. Default is dmin=1km.
- D max: This is the upper limit for searching for an optimal d-param value. Unit is km. Default is dmax=500km.
- D step: This is the step size for searching for an optimal d-param value. Unit is km. Default is dstep=2km.
- D parameter: This is the 'D' parameter for computing the spatial weights. If this option is used, then the parameters: dmin, dmax, dstep and Wt are not used.
Final there are two special argument as:
- cut excess stations: If this option is enabled, then any station (from the input file) outside the region limit (passed in via the 'region' option) is not considered in the strain estimation.
- generate statistics: This option will create an output file, named 'strain_stats.dat', where estimation info and statistics will be written.

Note thast users can browse StrainTool’s dedicated web page (https://dsolab.github.io/StrainTool/) for a more detailed description.

### Output files

Results of `StrainWebTool` are recorded in the following three files:

*   **strain_info.dat :** This file includes strain tensor parameters, principal axis, rotational rates, dilatation etc.  
    The columns of the file are structured as below:

<pre id="block-samp" <samp="">Latitude  Longtitude     vx+dvx          vy+dvy           w+dw          exx+dexx        exy+dexy        eyy+deyy       emax+demax      emin+demin       shr+dshr        azi+dazi      dilat+ddilat   sec.inv.+dsec.inv.
   deg       deg         mm/yr           mm/yr          deg/Myr       nstrain/yr      nstrain/yr      nstrain/yr      nstrain/yr      nstrain/yr      nstrain/yr         deg.         nstrain/yr      nstrain/yr   
	        </pre>

*   **station_info.dat :** Stations' data used for the calculation of strain tensor are written at htis file. Format is:

<pre id="block-samp" <samp="">
Code    Longtitude Latitude  Ve   Vn   dVe  dVn 
string      deg       deg          mm/yr
	      </pre>

*   **strain_stats.dat :** Output file for statistics:
<pre id="block-samp" <samp=""> --HEADER-- 
Parameters and arguments used for estimation of strain tensors.
--statistics--
Longtitude  Latitude  # stations D (optimal)  CutOff dis.     Sigma
 deg.       deg.        #           Km           #            / 
 </pre>

## Active Map visualize strain tensor results

For the results visualization, the application uses an active map developed using the Leaflet javascript library. In this map, user can choose to plot principal axes, shear strain, dilatation or second invariant results. In addition, user can add as separate layer the stations and their respective velocities used to estimate strain tensors.

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

*  [Flask](http://flask.pocoo.org/)  microframework for Python (v1.0.2). 

* [Bootstrap](https://getbootstrap.com/) open source toolkit (v4.2).

* [Leaflet](https://leafletjs.com/) open-source JavaScript library (v1.4.0).









