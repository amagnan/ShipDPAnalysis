import ROOT as r
from rootpyPickler import Unpickler
import shipRoot_conf
from ShipGeoConfig import ConfigRegistry 
shipRoot_conf.configure()
fgeo = r.TFile.Open("$ALP/geofile_full.conical.ALPACA-TGeant4.root")
sGeo = r.gGeoManager
upkl    = Unpickler(fgeo)
ShipGeo = upkl.load('ShipGeo')
fairGeo   = fgeo.FAIRGeom
top = sGeo.GetTopVolume()
