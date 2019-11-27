#!/bin/bash -f
###### a script to use dynamical analogue method to predict short-term climate ######
###### based on NCEP CFSv2 operational products by Li Zhenkun, Jul. 13, 2016   ######

###### set some environment varibles ################################################
export CP_INTEL_DIR=/opt/software/compiler/intel/composer_xe_2015.2.164
export LD_LIBRARY_PATH=${CP_INTEL_DIR}/compiler/lib/intel64:${CP_INTEL_DIR}/mkl/lib/intel64
export LIBRARY_PATH=${CP_INTEL_DIR}/compiler/lib/intel64:${CP_INTEL_DIR}/mkl/lib/intel64
#export GRIB_API_LIB=/opt/software/libs/grib_api/1.13.0/intel/lib
#export GRIB_API_INCLUDE=/opt/software/libs/grib_api/1.13.0/intel/include
#export LD_LIBRARY_PATH=${GRIB_API_LIB}:${LD_LIBRARY_PATH}
export ECCODE_LIB=/home/lizhenkun/eccodes/intel/lib
export ECCODE_INCLUDE=/home/lizhenkun/eccodes/intel/include
export LD_LIBRARY_PATH=${ECCODE_LIB}:${LD_LIBRARY_PATH}
export NCARG_ROOT=/opt/software/libs/ncl/6.1.0/dap
export PATH=$NCARG_ROOT/bin:/opt/software/libs/GraphicsMagick/1.3.20/bin:$PATH
#####################################################################################

###### set some useful parameters for convenience ######
dt=3
date=`date -d " $dt days ago " +%Y%m%d`
year=`echo $date | cut -c 1-4`
mon=`echo $date | cut -c 5-6`
day=`echo $date | cut -c 7-8`
echo 'the current date is '$year $mon $day
if [ $mon != 12 ]; then
  ave_days=`cal $((10#$mon+1)) $year | xargs | awk '{print $NF}'`
else
  ave_days=31
fi
lead_days=10
domain_west=0
domain_east=180
domain_south=0
domain_north=85
start_year=1999
end_year=2010
samp_interval=2
samp_sgl_len=45
dap_dir=/home/lizhenkun/cfs/dap
input_dir_curt=/home/lizhenkun/cfs/cfs_operation_output
input_dir_hist=/data/cfs_reforecast_output
input_dir_prec=/home/lizhenkun/cfs/dap

daysofmon=`cal $mon $year | xargs | awk '{print $NF}'`
echo 'the days of current month is '$daysofmon
if [ $((daysofmon-10#$day+1)) != $lead_days ]; then
  echo 'it is not time to do dynamical analogue prediction, exit!'
  exit
fi

if [ ! -d $dap_dir/$date ]; then
  mkdir $dap_dir/$date
fi
output_dir=$dap_dir/$date

raw_input_dir=/data/cfs
decode_output_dir=/home/lizhenkun/cfs/cfs_operation_output
decode_exe_dir=/home/lizhenkun/cfs/cfs_decode
decode_exe=cfs_decode.exe
########################################################

echo '=====> Check if current forecast data exists <====='
for var_name in z500 z200
do
  if [ ! -e $decode_output_dir/$var_name/$var_name-$date-cfs-forecast.dat ]; then
    if [ -e $decode_output_dir/namelist ]; then
      rm -f $decode_output_dir/namelist
    fi
    cat > $decode_output_dir/namelist << EOF
&params
start_date    = $date
end_date      = $date
var_name      = '$var_name'
input_dir     = '$raw_input_dir'
output_dir    = '$decode_output_dir/$var_name'
/
EOF
    echo '  decoding '$var_name' now ...'
    $decode_exe_dir/$decode_exe $decode_output_dir/namelist
  fi
done

cd $dap_dir
echo '=====> Begin to do dynamical analogue prediction <====='
if [ -e namelist ]; then
  rm -f namelist
fi
echo '  create namelist file'
cat > namelist << EOF
&params
fct_date        = $date
lead_days       = $lead_days
ave_days        = $ave_days
domain_west     = $domain_west
domain_east     = $domain_east
domain_south    = $domain_south
domain_north    = $domain_north
start_year      = $start_year
end_year        = $end_year
samp_interval   = $samp_interval
samp_sgl_len    = $samp_sgl_len
anal_fct_num    = 10
var_num         = 2
var_names       = 'z200', 'z500'
input_dir_curt  = '$input_dir_curt'
input_dir_hist  = '$input_dir_hist'
input_dir_prec  = '$input_dir_prec'
output_dir      = '$output_dir'
path_separator  = '/'
/
EOF
./dap.exe namelist

echo '  plot dynamical analogue prediction results'
lead_days_str=`echo $(( $lead_days+100 )) | cut -c 2-3`
ave_days_str=`echo $(( $ave_days+100 )) | cut -c 2-3`
ncl 'date="'$date'"' 'lead_days="'$lead_days_str'"' 'ave_days="'$ave_days_str'"'  plot_precip_forecast.ncl
cd $output_dir
for((i=1;i<=12;i++))
do
  istr=`echo $(( $i+100 )) | cut -c 2-3`
  gm convert -rotate -90 -units PixelsPerInch -density 100 $date-lead-$lead_days_str-for-$ave_days_str-forecast-$istr.ps $date-lead-$lead_days_str-for-$ave_days_str-forecast-$istr.gif
done
echo 'All jobs successfully finished'
