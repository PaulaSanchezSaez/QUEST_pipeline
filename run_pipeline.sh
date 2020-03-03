#!/bin/bash

period="201203"


first="1st"
second="2nd"
stf="01"
finf="15"
sts="16"
fins="31"
mid="15"


python dark_comb_parallel.py -p $period
python flat_dark_sub_m.py -p $period
python flat_dark_sub_e.py -p $period
python flat_comb_parallel.py -p $period$first -s $period$stf -f $period$finf
#python flat_comb_parallel.py -p $period$second -s $period$sts -f $period$fins
python reduce_science_img.py -p $period -m $period$mid
python astrometry_parallel.py -p $period
python sel_sci_img.py -p $period
python calibration_parallel.py -p $period
python copy_reduced_images.py -p $period
