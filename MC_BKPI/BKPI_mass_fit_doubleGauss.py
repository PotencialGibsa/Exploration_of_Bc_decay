from ROOT import *
import datetime
#####The fitting of the B+ from the decay: B+ -> K+ J/psi
#####The best parameters of the fit are in the picture "Parameters_of_fit.jpg" in this directory

FITTING = False
DECATYPE = 1

decatypes_str = ["B_{c} -> anti(K^{0}) B^{+}", "B_{c} -> anti(K^{*0}) B^{+} " , "B_{c} -> B^{+} #pi^{+} K^{-}", "B_{c} -> B^{*+} #pi^{+} K^{-}",
 				 "B_{c} -> B^{*0}_{s2} #pi^{+} ; B^{*0}_{s2} -> B^{*+} K^{-}", "B_{c} -> B^{*0}_{s2} #pi^{+} ; B^{*0}_{s2} -> B^{+} K^{-}",
				 "B_{c} -> B^{0}_{s1} #pi^{+}"]

PDG_B_MASS = 5.279
PDG_BC_MASS = 6.2749

time = datetime.datetime.now().timetuple()
time_st = str(time[3])+'_'+str(time[4])+'_'+str(time[5])+'_'

ch = TChain('mytree')
ch.Add('/eos/home-d/dshmygol/SimpleFile_MC_Bc_in_BKPI_v0v1.root')

nEvt = ch.GetEntries()
print("Number of events:" , nEvt)

x = RooRealVar('Bc','Bc',5.9, 6.60)
varset = RooArgSet(x)
dh = 0
dh = RooDataSet ('dh','Dataset', varset)
N_cand = 0
for evt in range(nEvt):
	if ch.GetEntry(evt) <= 0: break
	if evt % 100000 == 0:
	       print ("Events proceeded:" , evt)

	gen_muonP_delta = ch.gen_muonP_delta
	gen_muonN_delta = ch.gen_muonN_delta
	gen_kaon1_delta = ch.gen_kaon1_delta
    	gen_kaon2_delta = ch.gen_kaon2_delta
    	gen_pion1_delta = ch.gen_pion1_delta
    	Bc_mass = ch.Bc_mass

	if Bc_mass < 5.9: continue #6.2
	if Bc_mass > 6.6: continue #6.35

	#cuts
	if gen_muonP_delta > 0.1 : continue
	if gen_muonN_delta > 0.1 : continue
    	if gen_kaon1_delta > 0.1 : continue
    	if gen_kaon2_delta > 0.1 : continue
    	if gen_pion1_delta > 0.1 : continue


	decatype = ch.decatype

	if decatype != DECATYPE: continue

	x.setVal(Bc_mass)
	dh.add(varset)
	N_cand += 1

#############################################
if FITTING :
	S = RooRealVar('S','Signal',130203,0,9999999)
	S_mean  = RooRealVar ( "S_mean" , "mean "   , PDG_BC_MASS, 6.18  , 6.32       )
	S1_sigma = RooRealVar ("S1_sigma", "sigma"   , 0.01, 0.0001, 10   )
	S2_sigma = RooRealVar ("S2_sigma", "gamma"   , 0.06, 0.001, 10    )
	S_f     = RooRealVar ("S_f"     , "S_f"       , 0.54  , 0 , 1       )


	pdfS1 = RooGaussian ("pdfS1", "gaus", x, S_mean, S1_sigma)
	pdfS2 = RooGaussian ("pdfS2", "gaus", x, S_mean, S2_sigma)

	pdfSig =RooAddPdf  ("pdfSig", "pdfSig", RooArgList(pdfS1, pdfS2), RooArgList(S_f))

	B       = RooRealVar ( "B"      , "B"       , 30000 , 0     , 999999)
	#B_al    = RooRealVar ( "B_al"   , "B_{#alpha}", -0.326, -10, 100.0   )
	k    = RooRealVar ( "k"   , "k", -0.326, -10.0, 100.0   )
	#b    = RooRealVar ( "b"   , "b", 0.0, -100, 100   )
	#pdfG    = RooGenericPdf ("pdfG"    , "exp(@0* @1) ", RooArgList (x, k) )
	pdfG    = RooGenericPdf ("pdfG"    , "exp(@0* @1)", RooArgList (x, k) )

	alist1  = RooArgList (pdfSig, pdfG)
	alist2 = RooArgList (S, B)
	pdfSum  = RooAddPdf  ("model", "model", alist1, alist2)

	S_mean.setConstant(True)
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	S_mean.setConstant(False)
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	S1_sigma.setConstant(False)
    	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
    	S2_sigma.setConstant(False)
    	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	S_f.setConstant(False)
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
    	k.setConstant(False)
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.Save() )
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.Save() )
	rrr.Print()

####################################

binN = 100
xframe = x.frame(binN)

xframe.SetTitle("Decatype " + str(DECATYPE));
dh.plotOn(xframe, RooFit.MarkerSize(0.6), RooFit.Name('data'))

if FITTING:
	pdfSum.plotOn(xframe, RooFit.Name('fitt'), RooFit.Range(x.getMin(), x.getMax()),RooFit.LineColor(kGreen) )

	H_FITRES = TH1F('H_FITRES', 'c2df-ndof-chi2-prob', 4, 1, 5)
	FIT_c2df    = xframe.chiSquare(rrr.floatParsFinal().getSize() )
	FIT_ndof    = (binN - rrr.floatParsFinal().getSize())
	FIT_chi2    = FIT_c2df * FIT_ndof
	H_FITRES[3] = FIT_chi2
	H_FITRES[1] = FIT_c2df; H_FITRES[2] = FIT_ndof
	H_FITRES[4] = TMath.Prob(FIT_chi2, FIT_ndof)
	####

	pdfSum.plotOn(xframe,RooFit.Components('pdfG'), RooFit.LineColor(kRed), RooFit.LineWidth(3))
	pdfSum.plotOn(xframe,RooFit.Components('pdfSig'), RooFit.LineColor(kBlue), RooFit.LineWidth(3))
	dh.plotOn(xframe, RooFit.MarkerSize(0.6), RooFit.Name('data'))

	print ('chisq %f, ndf %i, binN %i, chisqndf %f, prob %f' % (H_FITRES[3], H_FITRES[2], binN, H_FITRES[1], H_FITRES[4]))

print("The number of candidates on the picture is : ", N_cand)

canvas = TCanvas("canvas")
canvas.cd()
xframe.Draw()
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.05)
latex.DrawLatex(0.55,0.85, decatypes_str[DECATYPE])

if FITTING :
	latex.SetTextSize(0.04)
	latex.DrawText(0.55,0.8,'B = ' + str(pdfSum.getParameters(dh)["B"].getValV()))
	latex.DrawText(0.55,0.75,'S = ' + str(pdfSum.getParameters(dh)["S"].getValV()))
	latex.DrawLatex(0.55,0.70,'S_mean = ' + str(pdfSum.getParameters(dh)["S_mean"].getValV()))
	latex.DrawLatex(0.55,0.65,'S1_sigma = ' + str(pdfSum.getParameters(dh)["S1_sigma"].getValV()))
	latex.DrawLatex(0.55,0.60,'S2_sigma = ' + str(pdfSum.getParameters(dh)["S2_sigma"].getValV()))
	latex.DrawText(0.55,0.55,'k = ' + str(pdfSum.getParameters(dh)["k"].getValV()))
	latex.DrawText(0.55,0.45,'chi/ndf = ' + str(H_FITRES[3]) + "/" + str(H_FITRES[2]))

#canvas.SaveAs('Bc_mass_fit_decatype_5.png')
canvas.SaveAs('Bc_mass_fit_decatype_' + str(DECATYPE) + '.png')
