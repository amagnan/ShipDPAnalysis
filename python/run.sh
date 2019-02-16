python discard.py -p qcd    -d 190214      > ~/DarkPhotonAnalysis/data/190214/qcd_Ana_rate.dat 
echo "qcd rate done"
python discard.py -p qcd    -d 190214 -l 1 > ~/DarkPhotonAnalysis/data/190214/qcd_Ana_rate_L.dat 
echo "qcd rate_L done"
python discard.py -p pbrem  -d 190214      > ~/DarkPhotonAnalysis/data/190214/pbrem_Ana_rate.dat 
echo "pbrem rate done" 
python discard.py -p pbrem  -d 190214 -l 1 > ~/DarkPhotonAnalysis/data/190214/pbrem_Ana_rate_L.dat 
echo "pbrem rate_L done"
python discard.py -p meson  -d 190214      > ~/DarkPhotonAnalysis/data/190214/meson_Ana_rate.dat
echo "meson rate done"
python discard.py -p meson  -d 190214 -l 1 > ~/DarkPhotonAnalysis/data/190214/meson_Ana_rate_L.dat 
echo "meson rate_L done"
