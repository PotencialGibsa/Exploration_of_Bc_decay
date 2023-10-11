from ROOT import *
import datetime
#####The fitting of the B+ from the decay: B+ -> K+ J/psi
#####The best parameters of the fit are in the picture "Parameters_of_fit.jpg" in this directory

ch1 = TChain('mytree')
ch1.Add('/eos/home-d/dshmygol/SimpleFile_MC_Bc_in_BKPI_v7_trigger_matching.root')

ch2 = TChain('mytree')
ch2.Add('/eos/home-d/dshmygol/SimpleFile_MC_Bc_in_BKPI_WS_v7_trigger_matching.root')

nEvt1 = ch1.GetEntries()
print("Number of events ch1:" , nEvt1)

nEvt2 = ch2.GetEntries()
print("Number of events ch2:" , nEvt2)

h1 = TH1F('Bc', 'Bc', 70, 5.9, 6.6) #5.9
h2 = TH1F('Bc_WS', 'Bc_WS', 70, 5.9, 6.6)

ch1.Draw("Bc_mas1 >> Bc", "gen_pion1_delta < 0.01 && gen_kaon1_delta < 0.01 && gen_kaon1_delta < 0.01")#&& (Bc_mas1 < 6.179 || Bc_mas1 > 6.32)")
ch2.Draw("Bc_mas1 >> Bc_WS", "gen_pion1_delta < 0.01 && gen_kaon2_delta < 0.01 && gen_kaon1_delta < 0.01")#&& (Bc_mas1 < 6.179 || Bc_mas1 > 6.32)")#, "Bc_mas1 > 5.9 && Bc_mas1 < 6.5 && (Bc_mas1 < 6.179 || Bc_mas1 > 6.283) && (Bc_mas1 < 6.232 || Bc_mas1 > 6.32)")

c=TCanvas()
h1.SetLineColor(2); h1.SetMarkerColor(2); #red right sign
h2.SetLineColor(4); h2.SetMarkerColor(4); #blue    wrong sign
h1.Draw('e0')
h2.Draw('e0same')
c.SaveAs('mass_comp.png')
