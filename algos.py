import pandas as pd
import numpy as np

def deltaR(a, b):
    dEta2 = (a['eta'] - b['eta'])**2
    dPhi = abs(a['phi'] - b['phi'])
    wrdPhi = 2 * np.pi - dPhi # wrapped around
    dPhi = np.minimum(dPhi, wrdPhi)
    dPhi2 = dPhi**2
    return np.sqrt(dEta2 + dPhi2)

def seedConeJetsAllEvents(events, cluster_distance=0.4, seed_threshold=0., jet_threshold = 0., pt_weighted_pos=True):
    jets = seedConeJets(events.loc[0],
                        cluster_distance = cluster_distance,
                        seed_threshold = seed_threshold,
                        jet_threshold = jet_threshold,
                        pt_weighted_pos = pt_weighted_pos)

    for i in range(1, events.index[-1]+1):
        jets = pd.concat([jets, seedConeJets(events.loc[i],
                            cluster_distance = cluster_distance,
                            seed_threshold = seed_threshold,
                            jet_threshold = jet_threshold,
                            pt_weighted_pos = pt_weighted_pos)])
    return jets

def seedConeJets(particles, cluster_distance=0.4, seed_threshold=0., jet_threshold = 0., pt_weighted_pos=True):
    index = particles.index[0]
    particles = particles.copy(deep=True)
    particles.sort_values('pt', ascending=False, inplace=True)

    jets = pd.DataFrame({'pt':[], 'eta':[], 'phi':[]})
    ijet = 0
    done = False
    while(not done):
        seed = particles.iloc[0]
        jet = {'pt' : 0., 'eta' : seed['eta'], 'phi' : seed['phi'], 'i':index}
        # sum the pt of particles within the cone
        incone = deltaR(particles, seed) < cluster_distance
        jet['pt'] = sum(particles[incone]['pt'])
        if pt_weighted_pos:
            jet['eta'] = np.average(particles[incone]['eta'], weights=particles[incone]['pt'])
            jet['phi'] = np.average(particles[incone]['phi'], weights=particles[incone]['pt'])
        if jet['pt'] > jet_threshold:
            jets = jets.append(jet, ignore_index=True)
        # remove the particles which were clustered to the seed
        particles = particles[~incone]
        if len(particles) == 0 or particles.iloc[0]['pt'] < seed_threshold:
            done = True
    jets.set_index('i', inplace=True)
    jets.sort_values('pt', inplace=True, ascending=False)
    return jets

