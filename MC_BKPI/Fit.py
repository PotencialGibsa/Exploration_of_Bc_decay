from ROOT import *
import datetime
#####The fitting of the B+ from the decay: B+ -> K+ J/psi
#####The best parameters of the fit are in the picture "Parameters_of_fit.jpg" in this directory

PDG_B_MASS = 5.279
PDG_BC_MASS = 6.2749

DECATYPE = 6
time = datetime.datetime.now().timetuple()
time_st = str(time[3])+'_'+str(time[4])+'_'+str(time[5])+'_'

ch = TChain('mytree')
ch.Add('/eos/home-d/dshmygol/SimpleFile_MC_Bc_in_BKPI_v0v1_test2.root')

nEvt = ch.GetEntries()
print("Number of events:" , nEvt)

x = RooRealVar('Bc_mass','Bc_mass',0, 15)
varset = RooArgSet(x)
dh = 0
dh = RooDataSet ('dh','Dataset', varset)
N_cand = 0
for evt in range(nEvt):
	if ch.GetEntry(evt) <= 0: break
	if evt % 100000 == 0:
	       print ("Events proceeded:" , evt)
	Bc_mass = ch.Bc_mass
	gen_muonP_delta = ch.gen_muonP_delta
	#cuts

	decatype = ch.decatype

	if decatype != DECATYPE: continue

	x.setVal(gen_muonP_delta)
	dh.add(varset)
	N_cand += 1
#############################################
'''
S = RooRealVar('S','Signal',130203,0,9999999)
S_mean  = RooRealVar ( "S_mean" , "mean "   , PDG_BC_MASS, 6.18  , 6.32       )
S1_sigma = RooRealVar ("S1_sigma", "sigma"   , 0.009, 0.001, 0.012   )
S2_sigma = RooRealVar ("S2_sigma", "sigma"   , 0.016, 0.009, 0.07    )
S_f     = RooRealVar ("S_f"     , "f"       , 0.76  , 0.000 , 1.0       )


pdfS1 = RooGaussian( "pdfS1"  , "gaus"    , x   , S_mean, S1_sigma)
pdfS2 = RooGaussian( "pdfS2"  , "gaus"    , x   , S_mean, S2_sigma)
pdfSig =RooAddPdf  ("pdfSig", "pdfSig", RooArgList(pdfS1, pdfS2), RooArgList(S_f))
#pdfSig = pdfS1

B       = RooRealVar ( "B"      , "B"       , 1 , 0     , 999999)
B_al    = RooRealVar ( "B_al"   , "B_{#alpha}", -0.326, -10, 100.0   )
k    = RooRealVar ( "k"   , "k", 1.0, -1.0, 2.0   )
b    = RooRealVar ( "b"   , "b", 0.0, -100, 100   )
pdfG    = RooGenericPdf ("pdfG"    , "@1 * @0 + @2 ", RooArgList (x, k,b) )

#nul = RooRealVar ( "nul"      , "nul"       , 0 , 0     , 0)

alist1  = RooArgList (pdfSig, pdfG)
alist2 = RooArgList (S, B)
pdfSum  = RooAddPdf  ("model", "model", alist1, alist2)

#nul.setConstant(True)
#S.setConstant(True)
k.setConstant(True)
b.setConstant(True)
S_mean.setConstant(True)
S_f.setConstant(True)
S1_sigma.setConstant(True)
S2_sigma.setConstant(True)
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
S_mean.setConstant(False)
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
S_f.setConstant(False)
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
S1_sigma.setConstant(False)
k.setConstant(False)
b.setConstant(False)
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
S2_sigma.setConstant(False)
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.Save() )
rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.Save() )
rrr.Print()

#latex.SetTextSize(0.04)
'''
'''
latex.DrawText(0.5,0.75,'B = ' + str(pdfSum.getParameters(dh)[0].getValV()))
latex.DrawText(0.5,0.75,'S = ' + str(pdfSum.getParameters(dh)[2].getValV()))
latex.DrawText(0.5,0.70,'S_mean = ' + str(pdfSum.getParameters(dh)[6].getValV()))
latex.DrawText(0.5,0.65,'S1_sigma = ' + str(pdfSum.getParameters(dh)[3].getValV()))
latex.DrawText(0.5,0.60,'S2_sigma = ' + str(pdfSum.getParameters(dh)[4].getValV()))
'''

#latex.DrawText(0.5,0.60,'Background = ' + str(pdfSum.getParameters(dh)[7].getValV()))
####################################
print("The number of candidates on the picture is : ", N_cand)

binN = 100
xframe = x.frame(binN)


print('chiSquare: ', xframe.chiSquare())

xframe.SetTitle("Decatype " + str(DECATYPE));
dh.plotOn(xframe, RooFit.MarkerSize(0.6), RooFit.Name('data'))

canvas = TCanvas("canvas")
canvas.cd()
xframe.Draw()
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.05)
latex.DrawLatex(0.15,0.85, " ")
latex.SetTextSize(0.04)


canvas.SaveAs('Delta_decatype_' + str(DECATYPE) + '.png')
