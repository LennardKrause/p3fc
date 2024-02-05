# p3fc

## 2024 update for pilatus3_frame_conversion

p3fc.py offers an easy-to-use GUI to convert PILATUS3 images to the Bruker .sfrm format and to generate X-ray aperture masks for the Bruker SAINT+ integration engine of the [APEX3 Software Package](https://www.bruker.com/products/x-ray-diffraction-and-elemental-analysis/single-crystal-x-ray-diffraction/sc-xrd-software/overview/sc-xrd-software/apex3.html).

It is currently designed to convert data collected at the following synchrotron facilities / beamlines:
  - [Advanced Photon Source](https://www.aps.anl.gov/) / 15ID-D
  - [SPring-8](http://www.spring8.or.jp/en/) / BL02B1
  - [Diamond Light Source](https://www.diamond.ac.uk/Home.html) / I19-1

The data can then be integrated using Brukers SAINT+ Integration Engine V8.35A or later
  - requirement that frames consist of multiples of 512 pixels has been lifted
  - any frame size is now allowed

## Installation
Install directly from the wheel (check releases)
```
python3 -m pip install p3fc-1.0.0-py3-none-any.whl
```
it is installed as a script, run ```p3fc``` from a terminal.

 ## Add circular region masks
   - black circle masks the area it covers
   - green circle unmasks the area it covers
   - artistic arrangements are welcome
   - follows hierarchal order (small > large)
 
 ## Good to know:
 - use ```Options -> Set Wavelength``` to overwrite the .inf info
 - the ```Mask``` menu offers some useful functions 
 - the initial rectangle & circle will always be on top
 - objects are allowed to be placed anywhere
 - to mask a corner
   - cover the corner with a black circle (as illustrated below)
   - move the green circle outside of the image area (not necessary)
   - the mask will have the proper dimensions
 - use the slider to adjust the contrast
 - statusbar shows the intensity at the current mouse position x, y
 
 ## Can learn new formats:
  - currently needs:
    - SAINT detector name tag (e.g. 'CMOS-PHOTONII')
    - detector dimensions (in pixels)
    - header size, offset in bytes (can be zero)
    - datatype (e.g. np.int32)
    - rotation necessary (e.g. is the omega circle horizontal)?
  - limited flexibility but nothing is impossible

## Important
   - This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
   - Bruker is not associated with this software and will not support this

## Requirements

#### [Python](https://www.python.org/) 3.5 or later

#### Libraries (tested with):
  - [numpy (1.26.3)](https://www.numpy.org/)
  - [scipy (1.11.4)](https://scipy.org/)
  - [pyqtgraph (0.13.3)](https://pyqtgraph.readthedocs.io/en/latest/)
  - [PyQt6 (6.6.1)](https://www.riverbankcomputing.com/software/pyqt/intro/)
   
## Program Layout
Filebrowser / Image Conversion - Tab | Draw Beamstop - Tab
------------------------------------ | -------------------
![img_gui_convert](https://user-images.githubusercontent.com/48315771/57973478-82a81c00-79a9-11e9-88e6-2addb86d70c7.png) | ![img_gui_draw](https://user-images.githubusercontent.com/48315771/57973484-9a7fa000-79a9-11e9-9144-379d21f10f01.png)

#### Filebrowser / Image Conversion
Use the filebrowser to navigate to the frame folder, it may take a while to read a folder containing a large number of frames! The output folder line (*Output Directory*) can be edited freely and non-existing folders will be created recursively. By default, the output directory is linked to the input directory and a suffix (*_sfrm*) is added automatically. If the *link?* box is unchecked the input and output fields (*Input* and *Output Directory*) can be selected manually to be controlled by the filebrowser, a green ring indicates the currently active field. The *ow* box toggles between overwrite/skip if the converted frame is already existing.

#### Draw Beamstop
Once a folder with valid frames is selected, the *Draw Beamstop* tab becomes available. The filebrowser is disabled during conversion, however, the drawing tab is not. It is recommended to start the frame conversion prior to drawing masks as it assures that the mask files are stored in the same folder as the converted frames. The image is shown in native resolution, use the scroll bars to navigate to the beamstop shadow. Drag and adjust the patches (rectangle, circle) to where they are needed. The dead areas of the PILATUS3 detector and bad pixels are masked automatically. If a patch is not needed, simply put it outside the image area. Add/Remove circles using the ```Mask```menu. A saved mask is indicated by a green dot in the lower right corner and a color change of the patches. Saving a mask stores the position and the shape of the patches.

The masks are saved to *Output Directory* and follow the naming convention used by SAINT so no further steps are needed in order to use the masks.
