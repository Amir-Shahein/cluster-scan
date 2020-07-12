from more_itertools import powerset
import math
import numpy as np



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
    
#     occ = (z1 + z2 + 2*z3 + z4) / ztot
#     res = {'p0': p0, 'p1': p1, 'p2': p2, 'p3': p3, 'p4': p4, 'occ':occ}
#     return(res)

# de = [-10, 0, -10]
# c = np.logspace(-7, -2, 100)
# c0 = 0.6
# res_ref = threeSiteWithOverlaps(c, de)
# occ_ref = res_ref['occ']
# plt.plot(c/c0, occ_ref, 'k', label=0)

# for i in [-8, -10, -12]:
#     de[1] = i
#     res = threeSiteWithOverlaps(c, de)
#     occ = res['occ']
#     plt.plot(c/c0, occ, label=i)

# plt.legend(loc='best', title='Site #2 Kd')
# plt.xscale('log')
# plt.ylim([0,2])
# plt.ylabel('<N>')
# plt.xlabel('TF concentration (M)')

# # plt.show()
# plt.savefig('threeSiteWithOverlaps')



def general_statmech_model(regObj, conc, concO=0.6):
    
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
        
        if i == len(regObj.spacing): # account for the last site (note: we'd only arrive here if it's not part of an ovlp cluster)
            
            delE = math.log(regObj.Score[i]/concO)
            
            relMult = (conc/concO)*math.exp(-delE)
            
            aggOcc = aggOcc + relMult/(1+relMult)
                                
        elif regObj.spacing[i] >= 0:    #if it's a non-overlapping site, just assume independent binding and add its occupancy to aggOcc
                                        # this is only fine if never end up on index corresponding to the final member of an ovlp cluster (else it would initiate here)
            
            delE = math.log(regObj.Score[i]/concO)
            
            relMult = (conc/concO)*math.exp(-delE)
            
            aggOcc = aggOcc + relMult/(1+relMult)
            
        else: #must be an overlapping site (negative spacing), note this is fine because an ovlp cluster can't start on the last site
            
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
                
                for k in range(start,end+1): #do the other sites k in the cluster
                    
                    tempList = []
                    
                    if k is not j: #ignoring that site j itself
                    
                        if (k<j and regObj.endPos[k] >= regObj.startPos[j]) or (k>j and regObj.startPos[k] <= regObj.endPos[j]): #overlap with site j?
                            
                            tempList.append(k) #for each site j, append all sites k that do overlap, into tempList
                    
                    eachSitesOverlaps.append(tempList) #then append this list of sites overlapping with j to eachSitesOverlaps, and move to next site j+1
            
            allCombosList = list(powerset(np.arange(start,end+1))) #generate all possible combinations of sites (including states impossible due to binding exclusivity)
                        
            possibleStates = restrict_to_possible_states(allCombosList, eachSitesOverlaps) #list of possible states (sites that con be concurrently occupied in the ovlp cluster)
        
        #after we're done with the overlapping cluster, move on
        i += 1

        
        
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
        using the lists of which sites overlap with which (eachSitesOverlaps)

    """
    
    possibleStates = [] #list of possible states (sites that con be concurrently occupied in the ovlp cluster)

    
    #This 4th order nested for loop looks complicated but it's actually relatively simple 
    
    for j in range(len(allCombosList)): # for each tuple j in allCombosList
        
        keep = True # start by assuming we will keep this state (i.e. all of the sites in the lits can be occupied at once)
        
        for x in range(len(eachSitesOverlaps)): # for each list x of overlapping sites
            
            count=0 # the number of sites from the state (j) that are in a given list of overlapping sites
            
            for k in range(len(allCombosList[j])): # for a given site k in state j
                
                for n in range(len(eachSitesOverlaps[x])): # iterate through the sites in list x
                    
                    if allCombosList[j][k] == eachSitesOverlaps[x][n]: # +1 for each site in state j found in list x
                        
                        count += 1
                        
            if count>1: # if more than one site in state j is found in any of the individual lists (x's) in eachSsitesOverlaps, game over for that state j
                keep = False
                break
            
        if keep == True: #if keep is still true, that state is possible, append it to possibleStates, otherwise move to the next state j+1
            possibleStates.append(allCombosList[j])
            
    return possibleStates
    
    
    