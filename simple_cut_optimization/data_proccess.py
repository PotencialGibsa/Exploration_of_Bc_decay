from ROOT import *
#gROOT.SetBatch(True)
import datetime
#####The fitting of the fon from the decay: Bc -> BKPI

FITTING = True

ch = TChain('mytree')
#ch.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2016_part_1.root')
#ch.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2016_part_2.root')
#ch.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2017_part_1.root')
#ch.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2017_part_2.root')
#ch.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2018_part_1.root')
#ch.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2018_part_2.root')
ch.Add('/eos/home-d/dshmygol/SimpleFile_MC_Bc_in_BKPI_v7_trigger_matching.root')

nEvt = ch.GetEntries()
print("Number of events:" , nEvt)

x = RooRealVar('x','x', 4.8, 6)
#x.setRange("SBL", 5.92, 6.179)
#x.setRange("SBR", 6.32, 6.6)
#stringforsidebands = 'Bc_mass<6.179||Bc_mass>6.32'

varset = RooArgSet(x)

dh = 0
dh = RooDataSet ('dh','Dataset', varset)
N_cand = 0
for evt in range(nEvt):
    if ch.GetEntry(evt) <= 0: break
    if evt % 100000 == 0:
        print ("Events proceeded:" , evt)
        
    Bc_mass = ch.Bc_mas1

    if Bc_mass < 5.92: continue
    if Bc_mass > 6.6: continue

    var_write = ch.Bc_masc0

    x.setVal(var_write)
    dh.add(varset)
    N_cand += 1

print("The number of candidates on the picture is : ", N_cand)

fi = TFile('Bc_masc0_dataset_MC.root', 'recreate')
dh.Write()
fi.Close()