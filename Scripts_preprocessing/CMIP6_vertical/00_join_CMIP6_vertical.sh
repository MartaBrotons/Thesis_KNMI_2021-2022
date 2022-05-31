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

## CMIP6

# hfls

for VAR in hur #hus850 hus925 #pr tas ts ua200 ua850 va200 va850 wap500 wap850 zg
do
    for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
    do
        for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
        do
            INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/${VAR}/areamean/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
            OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6_vertical/00_join_CMIP6_vertical_data/${VAR}/
            cdo mergetime ${INPUT}average_${VAR}_Lon270_310_Lat5_20_Amon_${MODEL}_historical_*.nc ${INPUT}average_${VAR}_Lon270_310_Lat5_20_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_average_Lon270_310_Lat5_20_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
        done
    done
done

# #hfss

# for VAR in hfss #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/hfss/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done

# #hur850

# for VAR in hur850 #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/hur/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done

# #hur925

# for VAR in hur925 #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/hur/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done


#hus850

# for VAR in hus850 #hus925 #pr tas ts ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/hus/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done

# #hus925

# for VAR in hus925 #hus925 #pr tas ts ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/hus/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done

# #ta200

# for VAR in ta200 #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/ta/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done

# #ta850

# for VAR in ta850 #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/ta/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done

# #ts

# for VAR in ts #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/ts/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done

# #ua200

# for VAR in ua200 #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/ua/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done

# #ua850

# for VAR in ua850 #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/ua/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done

# #va200

# for VAR in va200 #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/va/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done

# #va850

# for VAR in va850 #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/va/ #INPUT_${VAR} #/data/brotons/resampling/Method-2014/data_symlinks/tas/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done

# #wap500

# for VAR in wap500 #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/wap/
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo  mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done


# #wap850

# for VAR in wap850 #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/wap/ 
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done


# #zg850

# for VAR in zg850 #ua200 ua850 va200 va850 wap500 wap850 zg
# do
#     for MODEL in ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-CSM2-MR CanESM5 CanESM5-CanOE CESM2 CESM2-WACCM CIESM CMCC-CM2-SR5 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3 EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G HadGEM3-GC31-LL INM-CM4-8  INM-CM5-0 IPSL-CM6A-LR KACE-1-0-G MIROC6 MIROC-ES2L MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorESM2-MM UKESM1-0-LL 
#     do
#         for SSP in 585 #585MODEL in CanESM5 CanESM5-CanOE
#         do
#             INPUT=~/shared_data/volume_2/rhaarsma/DATA_MARTA/CMIP6_data/zg/ 
#             OUTPUT=~/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/00_join_CMIP6_models/${VAR}/
#             cdo mergetime ${INPUT}${VAR}_Amon_${MODEL}_historical_*.nc ${INPUT}${VAR}_Amon_${MODEL}_ssp${SSP}_*.nc ${OUTPUT}${MODEL}_Amon_historical-ssp${SSP}_195001-210012_${VAR}.nc
#         done
#     done
# done