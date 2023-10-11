from ROOT import *
import numpy as np
gROOT.SetBatch(True)
import datetime
#####The fitting of the fon from the decay: Bc -> BKPI

FITTING = False

chain = TChain('mytree')
chain.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2016_part_1.root')
chain.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2016_part_2.root')
chain.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2017_part_1.root')
chain.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2017_part_2.root')
chain.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2018_part_1.root')
chain.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2018_part_2.root')
#chain.Add('/eos/home-d/dshmygol/SimpleFile_MC_Bc_in_BKPI_v7_trigger_matching.root')
nEvt = chain.GetEntries()
print("Number of events:" , nEvt)

_fileOUT = 'cuts_B_pt.root'
fileOUT = TFile (_fileOUT, "recreate")

mytree = TTree("mytree","mytree")

my_vars = [
            'SAMEEVENT',
            #Mass Bc
            'Bc_mass',
            'Bc_pt', 'Bc_eta',
            'Bc_pvdistsignif2', 'Bc_pvdistsignif3',
            'Bc_pvdist', 'Bc_vtxprob',
            'Bc_pvcos2', 'Bc_pvcos3',
            #
            'BandK_mass',
            'BandPI_mass',
            'PIandK_mass',
            #Mass B+
            'B_mass',  'B_masc0',
            'B_pt', 'B_eta','BCV_pt',
            'B_pvdistsignif2', 'B_pvdistsignif3',
            'B_pvdist', 'B_vtxprob',
            'B_pvcos2', 'B_pvcos3',
            'B_Bcdistsignif2', 'B_Bcdistsignif3',
            'B_Bcdist',
            'B_Bccos2', 'B_Bccos3',
            #J/psi
            'J_psi_mass', 'J_psi_prob',
            'J_psi_pt',
            #Kaon1
            'K1_pt', 'K1_ips', 'K1_chrg','K1CV_pt',
            #Kaon2
            'K2_pt', 'K2_ips', 'K2_chrg','K2CV_pt',
            #pion3
            'PI3_pt', 'PI3_ips', 'PI3_chrg','PI3CV_pt',
            #mu1
            'mu1_pt', 'mu1_CV_pt',
            #mu2
            'mu2_pt', 'mu2_CV_pt',
			#trigger_matching
			'KA1_drTRG', 'KA1_dptTRG',
        ]

for _var in my_vars:
    exec(_var+'=np.zeros(1, dtype=float)')
    exec('mytree.Branch("' + _var + '"' + ' '*(25-len(_var)) + ',' + _var + ' '*(25-len(_var)) + ', "'+ _var + '/D")')

BBB = 0
# looping all over the events
for evt in range(nEvt):
    if chain.GetEntry(evt) <= 0:
        break
    #cuts
    if chain.Bc_mass < 5.92: continue
    if chain.Bc_mass > 6.6: continue
    
    if chain.Bc_pt < 65     : continue #30
    if chain.Bc_pt > 80     : continue #30
    """
    if chain.Bc_pvcos2 < 0.99 : continue
    
    if chain.Bc_eta > 0.26 : continue #or Bc_eta < 0.134757: continue
    
    if chain.Bc_pvdistsignif2 < 3: continue
    
    if chain.B_Bcdistsignif2 < 3  : continue
    
    if chain.K1_pt < 4: continue
    
    if chain.K2_pt < 2.5: continue #1.5
    
    if chain.PI3_pt < 1: continue
    #if KA2_pt < PI3_pt : continue #test &
    
    if chain.B_pt < 14.5     : continue
    
    if chain.Bc_pt < 20     : continue #30
    
    if chain.B_mass > 5.29 or chain.B_mass < 5.27: continue
    """
    '''
    if chain.BandK_mass > 6.134: continue
    
    if chain.BandPI_mass > 5.79 : continue
    
    if chain.PIandK_mass > 1: continue
    '''
    #
    
    if chain.KA1_dptTRG > 0.1: continue

    KA1_drTRG[0] = chain.KA1_drTRG
    KA1_dptTRG[0] = chain.KA1_dptTRG
    # assigning Kaon1 info
    K1_pt[0] = chain.K1_pt
    K1CV_pt[0] = chain.K1CV_pt
    K1_ips[0] = chain.K1_ips
    K1_chrg[0] = chain.K1_chrg
    # assigning Kaon2 info
    K2_pt[0] = chain.K2_pt
    K2CV_pt[0] = chain.K2CV_pt
    K2_ips[0] = chain.K2_ips
    K2_chrg[0] = chain.K2_chrg
    # assigning Pion3 info
    PI3_pt[0] = chain.PI3_pt
    PI3CV_pt[0] = chain.PI3CV_pt 
    PI3_ips[0] = chain.PI3_ips
    PI3_chrg[0] = chain.PI3_chrg
    # assigning muon+ info
    mu1_pt[0] = chain.mu1_pt
    mu1_CV_pt[0] = chain.mu1_CV_pt
    # assigning muon- info
    mu2_pt[0] = chain.mu2_pt
    mu2_CV_pt[0] = chain.mu2_CV_pt
    #    # assigning J/psi info
    J_psi_mass[0] = chain.J_psi_mass
    J_psi_prob[0] = chain.J_psi_prob
    J_psi_pt[0] = chain.J_psi_pt
    # assigning info about B+
    B_pt[0]= chain.B_pt
    BCV_pt[0] = chain.BCV_pt
    B_eta[0]= chain.B_eta
    B_vtxprob[0] = chain.B_vtxprob
    #B_mass[0]= chain.B_mass
    B_mass[0]= chain.B_mass
    B_masc0[0]= chain.B_masc0
    B_pvdistsignif2[0] = chain.B_pvdistsignif2
    B_pvdistsignif3[0] = chain.B_pvdistsignif3
    B_pvdist[0] = chain.B_pvdist
    B_pvcos2[0] = chain.B_pvcos2
    B_pvcos3[0] = chain.B_pvcos3

    Bc_pt[0]= chain.Bc_pt
    Bc_eta[0]= chain.Bc_eta
    Bc_vtxprob[0] = chain.Bc_vtxprob
    Bc_mass[0]= chain.Bc_mass
    Bc_pvdistsignif2[0] = chain.Bc_pvdistsignif2
    Bc_pvdistsignif3[0] = chain.Bc_pvdistsignif3
    Bc_pvdist[0] = chain.Bc_pvdist
    Bc_pvcos2[0] = chain.Bc_pvcos2
    Bc_pvcos3[0] = chain.Bc_pvcos3

    #otlet B+ Bc

    B_Bcdistsignif2[0] = chain.B_Bcdistsignif2
    B_Bcdistsignif3[0] = chain.B_Bcdistsignif3
    B_Bcdist[0] = chain.B_Bcdist
    B_Bccos2[0] = chain.B_Bccos2
    B_Bccos3[0] = chain.B_Bccos2

    
    SAMEEVENT[0] = chain.SAMEEVENT
    
    # fill _var branches in the tree if passes preselection cuts
    mytree.Fill()
                #
    if (evt % 10000 == 0):
        print('Running... Now in event %i / %i ' %(evt, nEvt))

print('Total entries stored in MyTree: %s' %(mytree.GetEntries()))
fileOUT.Write()
