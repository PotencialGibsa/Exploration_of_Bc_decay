from ROOT import *
import datetime
import glob
gROOT.SetBatch(True)
gStyle.SetOptStat(0000)
#####The fitting of the B+ from the decay: B+ -> K+ J/psi
#####The best parameters of the fit are in the picture "Parameters_of_fit.jpg" in this directory

ch1 = TChain('mytree')
#ch1_names = glob.glob("/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/*")
ch1_names = glob.glob("/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_*")
#print(ch1_names)
for ch1_name in ch1_names:
	ch1.Add(ch1_name)

ch2 = TChain('mytree')
#ch2_names = glob.glob("/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI_wrong_charge/SimpleFile_Bc_in_BKPI_wrong_charge_trigger_matching_*")
ch2_names = glob.glob("/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI_wrong_charge/*")
ch2_names = glob.glob("/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI_WS/*")


#print(ch2_names)
for ch2_name in ch2_names:
	ch2.Add(ch2_name)

nEvt1 = ch1.GetEntries()
print("Number of events ch1:" , nEvt1)

nEvt2 = ch2.GetEntries()
print("Number of events ch2:" , nEvt2)

c=TCanvas("canvas", "canvas", 1200, 1000)
c.cd()

pad1 = TPad("pad1","pad1", 0.0,0.33,1.0,1.0)
#pad1.SetLeftMargin(0.0)
#pad1.SetRightMargin(0.0)
#pad1.SetTopMargin(0.0)
pad1.SetBottomMargin(0.0)
pad1.Draw()
pad2 = TPad("pad2","pad2", 0.0,0.0,1.0,0.33)
#pad1.SetLeftMargin(0.0)
#pad1.SetRightMargin(0.0)
#pad1.SetTopMargin(0.0)
pad1.SetBottomMargin(0.011)
pad2.Draw()

pad1.cd()

h1 = TH1F('Bc', 'Bc', 68, 5.915, 6.595) #5.9
h2 = TH1F('Bc_WS', 'Bc_WS', 68, 5.915, 6.595)

ch1.Draw("Bc_mass >> Bc","(Bc_mass < 6.179 || Bc_mass > 6.32)")
#ch2.Draw("Bch1.Draw("Bc_mass >> Bc_c")c_mass >> Bc_WS")#, "(Bc_mass < 6.179 || Bc_mass > 6.32)")#, "Bc_mas1 > 5.9 && Bc_mas1 < 6.5 && (Bc_mas1 < 6.179 || Bc_mas1 > 6.283) && (Bc_mas1 < 6.232 || Bc_mas1 > 6.32)")
ch2.Draw("Bc_mass >> Bc_WS", "(Bc_mass < 6.179 || Bc_mass > 6.32)")#, "Bc_mas1 > 5.9 && Bc_mas1 < 6.5 && (Bc_mas1 < 6.179 || Bc_mas1 > 6.283) && (Bc_mas1 < 6.232 || Bc_mas1 > 6.32)")


h2.Scale(h1.Integral()/h2.Integral())

h1.SetTitle("")
label = TPaveLabel( 2, 600, 3.5, 650, 'red - B+K-Pi+\nblue - B+K+Pi-' )
h1.SetLineColor(2); h1.SetMarkerColor(2); #red
h2.SetLineColor(4); h2.SetMarkerColor(4); #blue


h1.Draw('e1p')
h2.Draw('e1psame')
#c.SaveAs('BKPI_fons_compare.png')

pad2.cd()
h3 = TH1F('Bc_c', 'Bc_c', 68, 5.915, 6.595)
ch1.Draw("Bc_mass >> Bc_c")
h3.Divide(h2)

h3.SetTitle("")
h3.SetXTitle("\ M(B_c) [GeV]")
h3.SetTitleOffset(0.0,"x")
h3.SetMarkerColor(kBlack)

h3_y_ax = h3.GetYaxis()
h3_y_ax.SetLabelSize(0.1)
h3_y_ax.SetTitleSize(0.1)
h3_y_ax.SetNdivisions(10)

h3_x_ax = h3.GetXaxis()
h3_x_ax.SetLabelSize(0.1)
h3_x_ax.SetTitleSize(0.045)

h3_x_ax.SetTitle("\ M(B_c) GeV")

h3.Draw()
#c.SaveAs('ratio.png')
c.SaveAs("test.png")
"""
x = RooRealVar('Bc','Bc',5.9, 6.6)
varset = RooArgSet(x)
dh = 0
dh = RooDataSet ('dh','Dataset', varset)
N_cand = 0
for evt in range(nEvt):
	if ch.GetEntry(evt) <= 0: break
	if evt % 100000 == 0:
	       print ("Events proceeded:" , evt)

	Bc_mass = ch.Bc_mas1

	if Bc_mass < 5.9: continue
	if Bc_mass > 6.5: continue

	if (Bc_mass > 6.231 - 0.052) and (Bc_mass < 6.231 + 0.052):
		continue
	if (Bc_mass > 6.276 - 0.044) and (Bc_mass < 6.276 + 0.044):
		continue

	x.setVal(Bc_mass)
	dh.add(varset)
	N_cand += 1

binN = 100
xframe = x.frame(binN)

xframe.SetTitle("Background_for_MC_BKPI")
dh.plotOn(xframe, RooFit.MarkerSize(0.6), RooFit.Name('data'))
"""
"""
print("The number of candidates on the picture is : ", N_cand)

canvas = TCanvas("canvas")
canvas.cd()
xframe.Draw()
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.05)
#latex.DrawLatex(0.55,0.85, "MC_BKPI_Background")

#canvas.SaveAs('Bc_mass_fit_decatype_5.png')
canvas.SaveAs('MC_BKPI_Background.png')
"""
