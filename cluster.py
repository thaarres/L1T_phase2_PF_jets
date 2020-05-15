import util, algos
from tests import compareCollections
import argparse
import pandas as pd
import sys
import regions
import numpy as np

def getEvents(input,loadFromh5=False):
    evs, jets = pd.DataFrame(),pd.DataFrame()
    if loadFromh5:
        print("Loading jets and pf candidates from h5")
        assert input.endswith('.h5'), "Wrong file type, need .h5!"
        evs  = pd.read_hdf(input, 'events')
        jets = pd.read_hdf(input, 'jets')
    else:
        print("Loading jets and pf candidates from CMSSW EDM file")
        assert input.endswith('.root'), "Wrong file type, need .root!"
        evs,jets = util.loadFromFile(input,writeToh5=True)   
     
    return evs,jets
        

if __name__ == "__main__":
    
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest='input', type=str, default="/afs/cern.ch/user/t/thaarres/public/TT_PU200.root")
  parser.add_argument('--fromh5', action='store_true')
  parser.add_argument('--doRegions', action='store_true')
  args = parser.parse_args()
  
  #Get all events, PF candidates and reconstructed jets
  events, ak4jets = getEvents(args.input,args.fromh5)
  
  ## Construct jets without splitting detector into regions
  # jets_nonregionized = algos.seedConeJetsAllEvents(events,jet_threshold=5.)
  
  #Construct jets splitting detector into regions
  if args.doRegions:
    eta_edges_all  = [[-3,3.],[-3, 0., 3.],[-3, 0., 3.],[-3, -1.0, 1.0, 3],[-3, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3]]
    phi_ndivisions = [2,2,3,3,10]
  
    jetCollections = { "df%i"%i: pd.DataFrame() for i in range(len(phi_ndivisions))}
    legends = []
    for i,(eta_edges,phi) in enumerate(zip(eta_edges_all,phi_ndivisions)):
      phi_edges = np.linspace(-np.pi, np.pi, phi)
      PFRegions = np.array([[regions.Region(eta_edges[i], eta_edges[i+1], phi_edges[j], phi_edges[j+1]) for j in range(len(phi_edges) - 1)] for i in range(len(eta_edges) - 1)])
      nRegions = (len(PFRegions)*len(PFRegions[0]))
      legends.append(r'Seed, %i regions'%nRegions)
      print("Using %i regions for jet reconstruction:"%nRegions)
      print("Eta bins: ")
      print(eta_edges)
      print("Phi bins: ")
      print(phi_edges)
    
      for index, event in events.groupby(level=0):
        if event.empty: continue
        particles_regionized = regions.regionize(event, PFRegions)
        for eta_bins in particles_regionized:
          for phi_bins in eta_bins:
            particles = phi_bins
            if not particles.empty:
              particles = particles.copy(deep=True)
              particles.sort_values('pt', inplace=True, ascending=False)
              jets = algos.seedConeJets(particles,jet_threshold=5.)

              if not jets.empty:
                jetCollections["df%i"%i] = jetCollections["df%i"%i].append(jets)
              
      jetCollections["df%i"%i].sort_values(by = ['i','pt'], inplace=True, ascending = [True, False])
    
    newJets=[ v for v in jetCollections.values() ]
    compareCollections(ak4jets,newCollections=newJets,baseCollectionName='AK4',newCollectionNames=legends,outpath="/eos/home-t/thaarres/www/L1T_phase2_PF_jets/regionCompare/")

  else:
     jets_nonregionized = algos.seedConeJetsAllEvents(events,jet_threshold=5.)
     compareCollections(ak4jets,newCollections=[jets_nonregionized],baseCollectionName='AK4',newCollectionNames=['Seed (1 region)'])
     