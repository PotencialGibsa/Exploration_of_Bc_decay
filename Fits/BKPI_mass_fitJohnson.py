from ROOT import *
import datetime
gROOT.SetBatch(True)
#####The fitting of the B+ from the decay: B+ -> K+ J/psi
#####The best parameters of the fit are in the picture "Parameters_of_fit.jpg" in this directory

FITTING = False
DECATYPE = 6

decatypes_str = ["B_{c} -> anti(K^{0}) B^{+}", "B_{c} -> anti(K^{*0}) B^{+} " , "B_{c} -> B^{+} #pi^{+} K^{-}", "B_{c} -> B^{*+} #pi^{+} K^{-}",
 				 "B_{c} -> B^{*0}_{s2} #pi^{+} ; B^{*0}_{s2} -> B^{*+} K^{-}", "B_{c} -> B^{*0}_{s2} #pi^{+} ; B^{*0}_{s2} -> B^{+} K^{-}",
				 "B_{c} -> B^{0}_{s1} #pi^{+}"]

PDG_B_MASS = 5.279
PDG_BC_MASS = 6.2749

time = datetime.datetime.now().timetuple()
time_st = str(time[3])+'_'+str(time[4])+'_'+str(time[5])+'_'

ch = TChain('mytree')
ch.Add('/eos/home-d/dshmygol/SimpleFile_MC_Bc_in_BKPI_WS_v7_trigger_matching.root')

nEvt = ch.GetEntries()
print("Number of events:" , nEvt)

x = RooRealVar('Bc_mass','Bc_mass',5.92, 6.6)
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
	Bc_mass = ch.Bc_mas1

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
	mu  = RooRealVar ( "mu" , "mean "   , 6.1416, 6.0  , 6.32       ) #PDG_BC_MASS
	lam = RooRealVar ("lam", "sigma"   , 0.154, 0.00001, 1   )
	gamma = RooRealVar ("gamma", "gamma"   , -0.04, -2, 1    )
	delta     = RooRealVar ("delta"     , "delta"       , 1.657  , 0.000001 , 10       )


	pdfSig = RooJohnson ("Johnson", "Johnson", x, mu, lam, gamma, delta)
	#pdfSig = pdfS1

	B       = RooRealVar ( "B"      , "B"       , 30000 , 0     , 999999)
	#B_al    = RooRealVar ( "B_al"   , "B_{#alpha}", -0.326, -10, 100.0   )
	k    = RooRealVar ( "k"   , "k", 17.0, 0.0, 20.0   )
	#b    = RooRealVar ( "b"   , "b", 0.0, -100, 100   )
	pdfG    = RooGenericPdf ("pdfG"    , "exp(@0* @1) ", RooArgList (x, k) )
	#pdfG    = RooGenericPdf ("pdfG"    , "@0* @1 + @2 ", RooArgList (x, k, b) )

	#nul = RooRealVar ( "nul"      , "nul"       , 0 , 0     , 0)

	alist1  = RooArgList (pdfSig, pdfG)
	alist2 = RooArgList (S, B)
	pdfSum  = RooAddPdf  ("model", "model", alist1, alist2)

	#nul.setConstant(True)
	#S.setConstant(True)
	k.setConstant(True)
	#b.setConstant(True)
	gamma.setConstant(True)
	delta.setConstant(True)
	mu.setConstant(True)
	lam.setConstant(True)
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	mu.setConstant(False)
	lam.setConstant(False)
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	gamma.setConstant(False)
	delta.setConstant(False)
	#lam.setConstant(False)
	k.setConstant(False)
	#b.setConstant(False)
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())

	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.Save() )
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.Save() )
	rrr.Print()

####################################

binN = 100
xframe = x.frame(binN)

#xframe.SetTitle("Decatype " + str(DECATYPE));
xframe.SetTitle("")
xframe.GetXaxis().SetTitle("M(Bc) [GeV]")
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
sdvig_x = 0.15
latex.SetTextSize(0.05)
latex.DrawLatex(sdvig_x,0.85, decatypes_str[DECATYPE])

if FITTING :
	latex.SetTextSize(0.04)
	latex.DrawText(sdvig_x,0.8,'background = ' + str(pdfSum.getParameters(dh)["B"].getValV()))
	latex.DrawText(sdvig_x,0.75,'signal = ' + str(pdfSum.getParameters(dh)["S"].getValV()))
	#latex.DrawLatex(sdvig_x,0.70,'#mu = ' + str(pdfSum.getParameters(dh)["mu"].getValV()))
	#latex.DrawLatex(sdvig_x,0.65,'#lambda = ' + str(pdfSum.getParameters(dh)["lam"].getValV()))
	#latex.DrawLatex(sdvig_x,0.60,'#gamma = ' + str(pdfSum.getParameters(dh)["gamma"].getValV()))
	#xlatex.DrawText(sdvig_x,0.55,'k = ' + str(pdfSum.getParameters(dh)["k"].getValV()))
	latex.DrawText(sdvig_x,0.70,'chi/ndf = ' + str(H_FITRES[3]) + "/" + str(H_FITRES[2]))

#canvas.SaveAs('Bc_mass_fit_decatype_5.png')
canvas.SaveAs('Bc_mass_fit_Johnson_decatype_' + str(DECATYPE) + '.png')
