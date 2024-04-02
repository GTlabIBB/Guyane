#!/usr/bin/env bash
touch commands.tmp
while read line;
do  start_day=`echo $line | cut -f1 -d","`;
end_day=`echo $line | cut -f2 -d","`;
echo -e "/usr/bin/time -v /home/etd530/Documents/Hysplit_R_interface/Hysplit_script/Hysplit_wind_analysis.R --from 2013-10-${start_day}-06-00 --to 2013-10-${end_day}-06-00 --lat=5.745974 --lon=-53.934047 --altitude 500,1000,2000 --byhour 1 --duration -200 --out guyane_vcardui_2013_byhour_Oct_500-1000-2000m_00-23h_10KmRes_200h_backwards_two-days_${start_day}-${end_day}.pdf --resolution 10000 --verbose --timezone="Brazil/East" --no_raster --margin=-60,-12,10,45 --windrose_times=-200 --cores 10 --run_id ${start_day} > guyane_vcardui_${start_day}-${end_day}.stdout.log 2>guyane_vcardui_${start_day}-${end_day}.stderr.log" >> commands.tmp
done < ../day_blocks_two-day.txt;
parallel -j 5 {} :::: commands.tmp && rm commands.tmp && touch DONE || touch FAILED
