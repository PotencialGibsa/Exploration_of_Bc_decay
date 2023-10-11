from ROOT import *
import datetime
#####The fitting of the B+ from the decay: B+ -> K+ J/psi
#####The best parameters of the fit are in the picture "Parameters_of_fit.jpg" in this directory

FITTING = False
DECATYPE = 2

PDG_B_MASS = 5.279
PDG_BC_MASS = 6.2749

time = datetime.datetime.now().timetuple()
time_st = str(time[3])+'_'+str(time[4])+'_'+str(time[5])+'_'

ch = TChain('mytree')
ch.Add('/eos/home-d/dshmygol/SimpleFile_MC_Bc_in_JpsiPi_v5.root')

nEvt = ch.GetEntries()
print("Number of events:" , nEvt)

x = RooRealVar('Bc_mass','Bc_mass', 6.0, 6.6)
varset = RooArgSet(x)
dh = 0
dh = RooDataSet ('dh','Dataset', varset)
N_cand = 0
for evt in range(nEvt):
	if ch.GetEntry(evt) <= 0: break
	if evt % 100000 == 0:
	       print ("Events proceeded:" , evt)
	decatype = ch.decatype

	Bc_mass = ch.B_mass
	B_pt = ch.B_pt
	B_pvdistsignif2 = ch.B_pvdistsignif2
	B_pvdistsignif3 = ch.B_pvdistsignif3
	B_pvdist = ch.B_pvdist
	B_pvcos2 = ch.B_pvcos2
	B_pvcos3 = ch.B_pvcos3
	J_psi_pt = ch.J_psi_pt
	PI1_pt = ch.PI1_pt
	PI1_ips = ch.PI1_ips
	B_vtxprob = ch.B_vtxprob
	PI1_ch = ch.PI1_chrg

	delta_PI1 = ch.delta_PI1
	delta_mu1 = ch.delta_mu1
	delta_mu2 = ch.delta_mu2

	if Bc_mass < 6.0: continue
	if Bc_mass > 6.6: continue

	#cuts from gen_mathcing
	if delta_PI1 > 0.02 : continue #0.002
	if delta_mu1 > 0.02 : continue #0.0015
 	if delta_mu2 > 0.02 : continue #0.0015

	#cuts
	if B_pt < 15 : continue
	if PI1_ips > 15 : continue
	if B_pvdistsignif2 < 5 : continue
	if B_pvcos2 < 0.995 : continue
	if PI1_pt < 3.5: continue
	if B_vtxprob < 0.1: continue

	if decatype != DECATYPE: continue

	x.setVal(Bc_mass)
	dh.add(varset)
	N_cand += 1
#############################################
if FITTING :
	S = RooRealVar('S','Signal', 858443,0,9999999)
	S_mean  = RooRealVar ( "S_mean" , "mean "   , PDG_BC_MASS, 6.18  , 6.32       )
	S1_sigma = RooRealVar ("S1_sigma", "sigma"   , 0.017959, 0.001, 0.07   )
	S2_sigma = RooRealVar ("S2_sigma", "sigma"   , 0.037969, 0.009, 0.07    )
	S_f     = RooRealVar ("S_f"     , "f"       , 0.32604  , 0.000 , 1.0       )


	pdfS1 = RooGaussian( "pdfS1"  , "gaus"    , x   , S_mean, S1_sigma)
	pdfS2 = RooGaussian( "pdfS2"  , "gaus"    , x   , S_mean, S2_sigma)
	pdfSig =RooAddPdf  ("pdfSig", "pdfSig", RooArgList(pdfS1, pdfS2), RooArgList(S_f))
	#pdfSig = pdfS1

	B       = RooRealVar ( "B"      , "B"       , 85800 , 0     , 999999)
	B_al    = RooRealVar ( "B_al"   , "B_{#alpha}", -0.326, -100, 10.0   )
	#k    = RooRealVar ( "k"   , "k", -2.0, -10.0, 10.0   )
	#b    = RooRealVar ( "b"   , "b", 12.6, -100, 100   )
	pdfG    = RooGenericPdf ("pdfG"    , "exp(@1*@0) ", RooArgList (x, B_al) )

	#nul = RooRealVar ( "nul"      , "nul"       , 0 , 0     , 0)

	alist1  = RooArgList (pdfSig, pdfG)
	alist2 = RooArgList (S, B)
	pdfSum  = RooAddPdf  ("model", "model", alist1, alist2)

	#nul.setConstant(True)
	#S.setConstant(True)
	#k.setConstant(True)
	#b.setConstant(True)
	S_mean.setConstant(True)
	S_f.setConstant(True)
	S1_sigma.setConstant(True)
	S2_sigma.setConstant(True)
	#k.setConstant(True)
	#b.setConstant(True)
	#B_al.setConstant(True)
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	S_mean.setConstant(False)
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	S1_sigma.setConstant(False)
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	S2_sigma.setConstant(False)
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	S_f.setConstant(False)
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.PrintLevel(-1), RooFit.Save())
	#k.setConstant(False)
	#b.setConstant(False)
	#rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.Save() )
	rrr = pdfSum.fitTo( dh, RooFit.NumCPU(7), RooFit.Save() )
	rrr.Print()

binN = 45
xframe = x.frame(binN)

xframe.SetTitle("Decatype " + str(DECATYPE));
dh.plotOn(xframe, RooFit.MarkerSize(0.6), RooFit.Name('data'))


if FITTING :
	pdfSum.plotOn(xframe, RooFit.Name('fitt'), RooFit.Range(x.getMin(), x.getMax()),RooFit.LineColor(kGreen) )

	H_FITRES = TH1F('H_FITRES', 'c2df-ndof-chi2-prob', 4, 1, 5)
	FIT_c2df    = xframe.chiSquare(rrr.floatParsFinal().getSize() )
	FIT_ndof    = (binN - rrr.floatParsFinal().getSize())
	FIT_chi2    = FIT_c2df * FIT_ndof; H_FITRES[3] = FIT_chi2;
	H_FITRES[1] = FIT_c2df; H_FITRES[2] = FIT_ndof
	H_FITRES[4] = TMath.Prob(FIT_chi2, FIT_ndof)

	pdfSum.plotOn(xframe,RooFit.Components('pdfG'), RooFit.LineColor(kRed), RooFit.LineWidth(3))
	pdfSum.plotOn(xframe,RooFit.Components('pdfSig'), RooFit.LineColor(kBlue), RooFit.LineWidth(3))
	pdfSum.plotOn(xframe,RooFit.Components('pdfS1'), RooFit.LineColor(kYellow), RooFit.LineWidth(3))
	pdfSum.plotOn(xframe,RooFit.Components('pdfS2'), RooFit.LineColor(kOrange), RooFit.LineWidth(3))
	dh.plotOn(xframe, RooFit.MarkerSize(0.6), RooFit.Name('data'))

	print ('chisq %f, ndf %i, binN %i, chisqndf %f, prob %f' % (H_FITRES[3], H_FITRES[2], binN, H_FITRES[1], H_FITRES[4]))

print("The number of candidates on the picture is : ", N_cand)

canvas = TCanvas("canvas")
canvas.cd()
xframe.Draw()
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.05)
latex.DrawLatex(0.1,0.85, "B_{c} -> J//#psi #rho^{+}")

if FITTING :
	latex.SetTextSize(0.04)
	latex.DrawText(0.1,0.8,'B = ' + str(pdfSum.getParameters(dh)["B"].getValV()))
	latex.DrawText(0.1,0.75,'S = ' + str(pdfSum.getParameters(dh)["S"].getValV()))
	latex.DrawText(0.1,0.70,'S_mean = ' + str(pdfSum.getParameters(dh)["S_mean"].getValV()))
	latex.DrawText(0.1,0.65,'S1_sigma = ' + str(pdfSum.getParameters(dh)["S1_sigma"].getValV()))
	latex.DrawText(0.1,0.60,'S2_sigma = ' + str(pdfSum.getParameters(dh)["S2_sigma"].getValV()))
	#latex.DrawText(0.15,0.55,'k = ' + str(pdfSum.getParameters(dh)["k"].getValV()))
	#latex.DrawText(0.15,0.50,'b = ' + str(pdfSum.getParameters(dh)["b"].getValV()))
	latex.DrawText(0.1,0.55,'B_al = ' + str(pdfSum.getParameters(dh)["B_al"].getValV()))
	latex.DrawText(0.1,0.45,'chi/ndf = ' + str(H_FITRES[3]) + "/" + str(H_FITRES[2]))


canvas.SaveAs('Bc_mass_fit_decatype_' + str(DECATYPE) + '.png')
