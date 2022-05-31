"""
for VAR in tas pr rsds
do
    for ripf in r1 r2 r4 r5 r6 r7
    do
    OUTPUT=/data/wiel/timeseries/resampleset_Aug2021/global/mon/EC-Earth3bis/
    INPUT=/data2/Common/ecearth3bis/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3/historical/${ripf}i1p5f1/Amon/${VAR}/gr/v20210701/
    cdo mergetime ${INPUT}*.nc ${OUTPUT}${VAR}_mon_EC-Earth3bis_historical_${ripf}i1p5f1.nc
    done
done

for VAR in tas pr rsds
do
    for ripf in r1 r2 r4 r5 r6 r7
    do
    OUTPUT=/data/wiel/timeseries/resampleset_Aug2021/global/mon/EC-Earth3bis/
    INPUT=/data2/Common/ecearth3bis/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3/ssp585/${ripf}i1p5f1/Amon/${VAR}/gr/v20210701/
    cdo mergetime ${INPUT}*.nc ${OUTPUT}${VAR}_mon_EC-Earth3bis_ssp585_${ripf}i1p5f1.nc
    done
done
"""
# CMIP6

#pr tas.  pr #tas #hus850 hus925 ts ua200 ua850 va200 va850 wap500 wap850 zg850

# for VAR in pr tas psl
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL
#     do
#         for SSP in 585 #126 245 585
#         do
#             INPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/${VAR}/
#             cdo remapbil,r360x180 ${INPUT}${MODEL}_Amon_historical-ssp585_185001-210012_${VAR}.nc ${OUTPUT}${MODEL}_Amon_historical-ssp585_185001-210012_${VAR}.nc
#         done
#     done
# done


for VAR in ua700 #hus700 hus1000 ua1000 ua925 ua700 va1000 va925 va700 #hfls hus850 hus925 ts ua200 ua850 va200 va850 wap500 wap850 zg850 
do
    for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL
    do
        for SSP in 585 #126 245 585
        do
            INPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
            OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/${VAR}/
            cdo remapbil,r360x180 ${INPUT}${MODEL}_Amon_historical-ssp585_195001-210012_${VAR}.nc ${OUTPUT}${MODEL}_Amon_historical-ssp585_195001-210012_${VAR}.nc
        done
    done
done