from ROOT import *
import math
gROOT.SetBatch(True)
import datetime

chain = TChain("mytree")

chain.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2016_part_1.root')
chain.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2016_part_2.root')
chain.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2017_part_1.root')
chain.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2017_part_2.root')
chain.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2018_part_1.root')
chain.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2018_part_2.root')

#cuts = "(K1_pt > 2.5) && (Bc_pt > 8) && (K2_ips > 0.1) && (Bc_pvdistsignif2 > 2) && (Bc_vtxprob > 0.01) && (B_Bcdistsignif2 > 0) && (K1_ips > 1.5) && (B_vtxprob > 0.04) && (PI3_pt > 0.25) && (PI3_ips > 0.1) && (B_Bccos2 > 0.99) && (B_pt > 5) && (K2_pt > 0.25) && (Bc_pvcos2 > 0.997)"
#cuts = "( K1_pt > 2.5 ) && ( Bc_pt > 8 ) && ( K2_ips > 0.01 ) && ( Bc_pvdistsignif2 > 2 ) && ( Bc_vtxprob > 0.001 ) && ( B_Bcdistsignif2 > 0 ) && ( K1_ips > 1.5 ) && ( B_vtxprob > 0.046 ) && ( PI3_ips > 0.01 )"
#cuts = "(B_mass < 5.33) && (B_mass > 5.23) && ( K1_pt > 2.5 ) && ( Bc_pt > 8 ) && ( K2_ips > 0.01 ) && ( Bc_pvdistsignif2 > 2 ) && ( Bc_vtxprob > 0.001 ) && ( B_Bcdistsignif2 > 0 ) && ( K1_ips > 1.5 ) && ( B_vtxprob > 0.046 ) && ( PI3_ips > 0.01 ) && ( PI3_pt > 0.05 ) && ( B_Bccos2 > 0.99 ) && (                B_pt > 7.5 ) && ( K2_pt > 0.05 ) && ( Bc_pvcos2 > 0.997 )"
#cuts = "(B_mass < 5.33) && (B_mass > 5.23) && ( K1_pt > 2 ) && ( Bc_pt > 8 ) && ( K2_ips > 0.01 ) && ( Bc_pvdistsignif2 > 2 ) && ( Bc_vtxprob > 0.001 ) && ( B_Bcdistsignif2 > 0 ) && ( K1_ips > 1.5 ) && ( B_vtxprob > 0.01 ) && ( PI3_ips > 0.01 ) && ( PI3_pt > 0.05 ) && ( B_Bccos2 > 0.99 ) && (                B_pt > 5 ) && ( K2_pt > 0.05 ) && ( Bc_pvcos2 > 0.996 )"   
#cuts = "(B_mass < 5.33) && (B_mass > 5.23) && (KA1_dptTRG < 0.1) && (KA1_drTRG < 0.1) && ( PI3_ips > 1.05 ) && ( Bc_pvdistsignif2 > 5.0 ) && ( Bc_vtxprob > 0.1 ) && ( B_pt > 10.0 ) && ( B_vtxprob > 0.02 ) && ( Bc_pvcos2 > 0.999 ) && ( K2_pt > 1.3 ) && ( K1_pt > 1.2 ) && ( B_Bcdistsignif2 > 2.0 ) && ( Bc_pt > 5.0 ) && ( K1_ips > 2.25 ) && ( B_Bccos2 > 0.99 ) && ( K2_ips > 0.04 ) && ( PI3_pt > 0.8 )"
#cuts = "(B_mass < 5.33) && (B_mass > 5.23) && ( Bc_vtxprob > 0.1 ) && ( B_Bccos2 > 0.998 ) && ( B_pt > 10.0 ) && ( K1_pt > 1.2 ) && ( K2_ips > 0.05 ) && ( Bc_pvdistsignif2 > 5.0 ) && ( PI3_pt > 0.75 ) && ( B_vtxprob > 0.02 ) && ( K2_pt > 1.25 ) && ( B_Bcdistsignif2 > 0.8 ) && ( K1_ips > 2.5 ) && ( PI3_ips > 1.1 ) && ( Bc_pvcos2 > 0.999 ) && ( Bc_pt > 5.0 )"
cuts = ""
variable = "Bc_mass"

canv = TCanvas("canvas", "canvas", 500,500, 2500,1500)

#pad3.SetGridx()
#pad3.SetGridy()
#pad3.GetFrame().SetFillColor( 18 )
#h1.SetMarkerStyle( 21 )
#h1.Draw( 'e1p' )
#label3 = TPaveLabel( 2, 600, 3.5, 650, 'option e1p' )
#label3.SetFillColor( 42 )
#label3.Draw()
num_bins = 68
min_range = 5.92
max_range = 6.6


signalN = chain.Draw(variable + ">>hist("+ str(num_bins) + "," + str(min_range) + "," + str(max_range) + ")", cuts)
hist = gDirectory.Get("hist")
hist.SetMarkerStyle( 21 )
hist.Draw( 'e1p' )
hist.SetTitle("")
y_axe = hist.GetYaxis()
y_axe.SetTitle("Candidates/ " + str((max_range - min_range)/num_bins *1000) + "MeV")
x_axe = hist.GetXaxis()
x_axe.SetTitle("\ M(B_c)  GeV") 

# Turn off the statistics display
gStyle.SetOptStat(0)

# Get the stats box object and modify its label
#stats = hist.GetListOfFunctions().FindObject("stats")
#stats.SetY1NDC(0.75)
#stats.SetY2NDC(0.95)
#stats.SetTextColor(kBlack)
#stats.SetOptStat(0)
#stats.SetOptFit(0)
#stats.SetLabel("Custom Label")

print("Evt = " , signalN)
canv.SaveAs("Bc_mass_plot.png")

###BK mass
"""
canv = TCanvas("canvas", "canvas", 500,500, 2500,1500)
h1 = TH1F('BandK_mass', 'BandK_mass', 136, 5.7, 6.5) #5.9
chain.Draw("BandK_mass >> BandK_mass", cuts )
canv.SaveAs("BandK_mass_plot.png")
"""
