from more_itertools import powerset
import math
import numpy as np

def single_site_statmech_model(KdRatio, conc=16, consensusKd=16, concO=0.6):
    """
    This method takes as input the affinity (ratio of Kd to Kd consensus) of a single transcription factor binding site, and outputs
    the mean occupancy based on a statistical mechanical model of transcription factor binding (Philips).
    

    Parameters
    ----------
    KdRatio : float
        Kd-site/Kd-consensus.
    conc : float, optional
        Concentration of transcription factor for the model. The default is 16, corresponding to the Kd of the consensus sequence
        (therefore, setting the occupancy of the consensus sequence to 0.5 -- similar to Maerkl 2007 and before)
    consensusKd : float, optional
        Kd of the consensus sequence, used to convert the KdRatio to actual Kd. The default is 16.
    concO : float, optional
        Reference concentration for the Philips'ian stat.mech. model. The default is 0.6.

    Returns
    -------
    meanOcc. The mean occupancy of the binding site, according to the supplied parameters.

    """
    delE = math.log(KdRatio*consensusKd/concO)
    
    relMult = (conc/concO)*math.exp(-delE)
    meanOcc = relMult/(1+relMult)
    
    return meanOcc
    


def general_cluster_statmech_model(regObj, lowAffinityOcc=None, conc=16, consensusKd=16, concO=0.6):
    
    """
    This method takes as input a regulatory element (module) object of binding sites 
    with their different properties (spacing, scores, etc.), as well as a concentration parameter,
    and returns the mean occupancy of transcription factor to the binding sites. It 
    works by assuming that non-overlapping sites bind transcription factor independently,
    and it deals with dependent (overlapping) binding sites by first identifying 
    all of the possible bound states of the system, followed by determining the 
    weighted multiplicities (relative to unbound state) based on Philips PBOC 2nd edition. 
    Then the mean occ. is calculated as a weighted average of the probability of occupancy of 
    each state (ratio of the state's rel. weighted multiplicity to the partition function (also
    rel.), and by the number of transcription factors that are bound in that state. 
    
    Based on PBOC - 2nd, p.270, equation 6.112:
    delE (change in energy with binding, relative to solution) in units Boltzmann constant*T (1/beta) = ln(Kd/concO)
    When calculating occupancy, I first convert Kds into this delE value in terms of Kb*T
    
    lowAffinityOcc: float
        If specified, this means we only care about the mean occupancy coming from low-affinity sites,
        which have Kd values above the threshold ratio stored in this variable. In order to implement
        this, for independent binding sites, we skip the occupancy contribution if the site's Kd is 
        below this threshold. For overlapping clusters, we multiply each state's weighted rel. mult.
        by the number of low-affinity binding sites instead of the total number of binding sites, before
        dividing by the partition function.
    
    Returns
    -------
    Mean occupancy of the input regulatory module.

    """
    
    aggOcc = 0 # aggregate occupancy, which will accumulate independently from 
               # the modules (single sites and overlapping clusters)
    
    # iterate over spacings list, separate into different cases if it's an 
    # overlapping site (dependent) or non-overlapping site (independent)
    
    i=0
    
    while i < len(regObj.startPos):
        
        # -------------- if it's a non-overlapping site, just assume independent binding and add its occupancy to aggOcc
        if i == len(regObj.spacing): # account for the last site (note: we'd only arrive here if it's not part of an ovlp cluster (but can't check spacing),
                                     # so it's fine not to check, and just assume non-overlapping)
            
            if lowAffinityOcc: # see parameter definition
                
                if regObj.siteAffinities[i]>=lowAffinityOcc:
                    
                    fracOcc = single_site_statmech_model(regObj.siteAffinities[i], conc, consensusKd, concO)
                    aggOcc = aggOcc + fracOcc
            
            else: # note: this gets skipped if lowAffinityOcc=True, but site affinity < lowAffinityOcc, and therefore higher affinity non-ovlp sites are not added
                
                fracOcc = single_site_statmech_model(regObj.siteAffinities[i], conc, consensusKd, concO)
                aggOcc = aggOcc + fracOcc
            
        elif regObj.spacing[i] >= 0: # this is only fine if never end here with index corresponding to the final member of an ovlp cluster (else it would trigger)
            
            if lowAffinityOcc:
                
                if regObj.siteAffinities[i]>=lowAffinityOcc:
                    
                    fracOcc = single_site_statmech_model(regObj.siteAffinities[i], conc, consensusKd, concO)
                    aggOcc = aggOcc + fracOcc
            
            else:
                    fracOcc = single_site_statmech_model(regObj.siteAffinities[i], conc, consensusKd, concO)
                    aggOcc = aggOcc + fracOcc
            
        else: # -------------- must be an overlapping site (negative spacing), note that an ovlp cluster can't start on the last site
            
            start = i
            
            allCombosList = [] # list of tuples of all combos of sites (powerset)
            eachSitesOverlaps  = [] # list of lists (one corresponding to each site), identifying which sites overlap
            
            while regObj.spacing[i]<0: #find the last site in the overlapping cluster
                i += 1
                
                if i == len(regObj.spacing): #if we're on the last member of both the ovlp cluster and also the entire reg element
                    break
            
            end = i
            #identify all binding sites that overlap with each binding site (create a separate list for each site), and append these to eachSitesOverlaps
            for j in range(start,end+1): #for a given site j in the ovlp cluster
                
                tempList = []
                
                for k in range(start,end+1): #do the other sites k in the cluster
                       
                    if k is not j: #ignoring that site j itself
                    
                        if (k<j and regObj.endPos[k] >= regObj.startPos[j]) or (k>j and regObj.startPos[k] <= regObj.endPos[j]): #overlap with site j?
                            
                            tempList.append(k) #for each site j, append all sites k that do overlap, into tempList
                    
                if tempList:
                    eachSitesOverlaps.append(tempList) #then append this list of sites overlapping with j to eachSitesOverlaps, and move to next site j+1
        
            allCombosList = list(powerset(np.arange(start,end+1))) #generate all possible combinations of sites (including states impossible due to binding exclusivity)
            possibleStates = restrict_to_possible_states(allCombosList, eachSitesOverlaps) #list of possible states (sites that con be concurrently occupied in the ovlp cluster)
            # iterate over possible states and generate a relative boltzmann-weighted multiplicity for each, by summing the delta binding energies (relative to Esol)
            aggWeightedMeanOcc = 0 # Accumulate the numerator in the: weighted average (Pocc) mean occupancy calculation  (note that unbound state is multiplied by 0 sites)
            aggPartitionFunction = 0 # Accumulate the denominator (note that the unbound state results in +1)
            
            for y in range(len(possibleStates)):
                
                aggDelE = 0
                if lowAffinityOcc:
                    numLowAffSites=0
                numSites = len(possibleStates[y])
                
                for z in range(numSites):
                    
                    siteIndex = possibleStates[y][z] # index (in regEle object) of current site in current binding state
                    aggDelE = aggDelE + math.log(regObj.siteAffinities[siteIndex]*consensusKd/concO) # iterate and sum the deltaEnergy terms for that state
                    
                    if lowAffinityOcc and (regObj.siteAffinities[siteIndex]>=lowAffinityOcc):
                        numLowAffSites += 1
                        
                relMultOvlp = ((conc/concO)**numSites)*math.exp(-aggDelE) # calculate the weighted relative multiplicity term for that state (see PBOC or notes for derivation)
                
                aggPartitionFunction = aggPartitionFunction + relMultOvlp
                
                if lowAffinityOcc:
                    aggWeightedMeanOcc = aggWeightedMeanOcc + numLowAffSites*relMultOvlp
                else:
                    aggWeightedMeanOcc = aggWeightedMeanOcc + numSites*relMultOvlp
                
            meanOcc = aggWeightedMeanOcc/aggPartitionFunction # do the mean occupancy calculation
            
            aggOcc = aggOcc + meanOcc # add the occupancy from the ovlp cluster to the occupancy for the overall regulatory element
                    
        # -------------- after we're done with the case (non-ovlp site or ovlp cluster), move on to the next site
        i += 1
        
    return aggOcc
        
        
def restrict_to_possible_states(allCombosList, eachSitesOverlaps):
    """
    

    Parameters
    ----------
    allCombosList : list of tuples
        Powerset (all possible combinations) of the sites in an overlapping cluster.
    eachSitesOverlaps : list of lists
        List of lists of which sites overlap with each other.

    Returns
    -------
    possibleStates : List of lists
        This method processes the allCombosList into a list of possible states (no overlapping sites in a given state),
        using the lists sites that overlap with each given site (eachSitesOverlaps)

    """
    
    possibleStates = [] #list of possible states (sites that con be concurrently occupied in the ovlp cluster)

    
    #This 4th order nested for loop looks complicated but it's actually relatively simple 
    
    for j in range(len(allCombosList)): # for each tuple j in allCombosList
        
        keep = True # start by assuming we will keep this state (i.e. all of the sites in the list can be occupied at once)
        
        for x in range(len(eachSitesOverlaps)): # for each list x of overlapping sites (same as numSites)
            
            applicable=False # if the site is present in the j state, to which list x actually corresponds
            overlapping=False # assuming the site corresponding to list x is present, is one of the sites that overlaps with it present?
            
            for k in range(len(allCombosList[j])): # for a given site k in state j
                
                for n in range(len(eachSitesOverlaps[x])): # iterate through the sites in list x
                    
                    if allCombosList[j][k]-allCombosList[1][0] == x and applicable == False: # is the site present in state j, to which the list X (of its overlapping sites) actually corresponds?
                        applicable = True                                                    # - allCombosList[1][0] is because we need to calibrate to that index
                    
                    if allCombosList[j][k] == eachSitesOverlaps[x][n] and overlapping == False: # if a site in state j is found in list x
                        overlapping=True
                    
                    if applicable == True and overlapping == True: # if both conditions are met at any point, the state is invalid, break out of loop
                        keep = False
                        break
                    
            if applicable == True and overlapping == True: # if the state is invalid at any point, break out of loop
                break
            
        if keep == True: #if keep is still true, that state is possible, append it to possibleStates. Regardless, move to the next state j+1
            possibleStates.append(allCombosList[j])
            
    return possibleStates
    





# ----- CODE TESTING -------
# def threeSiteWithOverlaps(c, de, c0=0.6):
#     de1, de2, de3 = de
#     z0 = 1
#     z1 = (c/c0)*np.exp(-de1)
#     z2 = (c/c0)*np.exp(-de3)
#     z3 = (c/c0)**2*np.exp(-de1-de3)
#     z4 = (c/c0)*np.exp(-de2)
#     ztot = z0 + z1 + z2 + z3 + z4
    
#     p0 = z0 / ztot
#     p1 = z1 / ztot
#     p2 = z2 / ztot
#     p3 = z3 / ztot
#     p4 = z4 / ztot
#     print((z1+z3)/ztot)
#     occ = (z1 + z2 + 2*z3 + z4) / ztot
#     res = {'p0': p0, 'p1': p1, 'p2': p2, 'p3': p3, 'p4': p4, 'occ':occ}
#     return(res)

# consensusKd=16

# c0=0.6
# de = [math.log(28.500423084880122*consensusKd/c0), math.log(12.339657959232085*consensusKd/c0), math.log(24.912849062155786*consensusKd/c0)]
# c=16
# res_trg = threeSiteWithOverlaps(c,de)
# occ = res_trg['occ']
# print(occ)
# p1 = res_trg['p1']

# plt.legend(loc='best', title='Site #2 Kd')
# plt.xscale('log')
# plt.ylim([0,2])
# plt.ylabel('<N>')
# plt.xlabel('TF concentration (M)')

# # plt.show()
# plt.savefig('threeSiteWithOverlaps')
# -------------------------------------------
    