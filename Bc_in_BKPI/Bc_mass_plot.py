from ROOT import *
gROOT.SetBatch(True)
import datetime
#####The fitting of the fon from the decay: Bc -> BKPI

FITTING = False # is responsible for making fit

# Adding processed files after Xb_frame and MySelector
ch = TChain('mytree')
ch.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2016_part_1.root')
ch.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2016_part_2.root')
ch.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2017_part_1.root')
ch.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2017_part_2.root')
ch.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2018_part_1.root')
ch.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2018_part_2.root')
nEvt = ch.GetEntries()
print("Number of events:" , nEvt)

x = RooRealVar('Bc_mass','Bc_mass', 5.92, 6.6)

varset = RooArgSet(x)
# Reading the data from root files and making cuts in needed
dh = 0
dh = RooDataSet ('dh','Dataset', varset)
N_cand = 0
for evt in range(nEvt):
    if ch.GetEntry(evt) <= 0: break
    if evt % 100000 == 0:
        print ("Events proceeded:" , evt)
        
    Bc_mass = ch.Bc_mass

    if Bc_mass < 5.92: continue
    if Bc_mass > 6.6: continue

    Bc_vtxprob = ch.Bc_vtxprob
    if Bc_vtxprob < 0.01 : continue

    B_vtxprob = ch.B_vtxprob
    if B_vtxprob < 0.06 : continue

    Bc_pvdistsignif2 = ch.Bc_pvdistsignif2
    if Bc_pvdistsignif2 < 2 : continue

    B_Bcdistsignif2 = ch.B_Bcdistsignif2
    if B_Bcdistsignif2 < 1 : continue
    PI3_pt = ch.PI3_pt
    if PI3_pt < 0.25 : continue
    KA2_pt = ch.K2_pt
    if KA2_pt < 0.25 : continue
    KA1_pt = ch.K1_pt
    if KA1_pt < 2.25: continue
   
    
    
    Bc_pvcos2 = ch.Bc_pvcos2
    B_Bccos2 = ch.B_Bccos2
    if B_Bccos2 < 0.998 : continue
    
    B_m = ch.B_mass
    if B_m < 5.23 :continue
    if B_m > 5.33 :continue

    x.setVal(Bc_mass)
    dh.add(varset)
    N_cand += 1
print("The number of candidates on the picture is : ", N_cand)

# fitting
if FITTING :
    PDG_B = 5.27934
    PDG_K = 0.493677
    PDG_PI = 0.13957061
    B       = RooRealVar ( "B"      , "B"       ,1.37983*10**7, 0, 14000000) 
    B_al    = RooRealVar ( "B_al"   , "B_{#alpha}", 1.27, 0.001, 5.0   )
    B_be    = RooRealVar ( "B_be"   , "B_{#beta}", 0.9753, 0.001, 5.0   )
    H_al = RooRealVar ( "H_al"   , "H_{#alpha}", 199.64, 0, 500)
    Bc_m = RooRealVar ( "Bc_m"   , "Bc_m",  PDG_B + PDG_K + PDG_PI, PDG_B + PDG_K + PDG_PI-0.9, PDG_B + PDG_K + PDG_PI+0.5 )
    C   = RooRealVar ( "C"   , "C", -0.8, -5.0, 3.0   )
    D   = RooRealVar ( "D"   , "D", 0.1, -50.0, 300.0   )
    f   = RooRealVar ( "f"   , "f", 0.1, .0, 1.0   )
    E_m = RooRealVar ( "E_m" , "E_m", 6.0081,5.9 , 6.1   )
    E_s = RooRealVar ( "E_s" , "E_s", 0.058365, 0.04, 0.06   )
    E_m_2 = RooRealVar ( "E_m_2" , "E_m_2", 5.9344, 4, 6   )
    E_s_2 = RooRealVar ( "E_s_2" , "E_s_2", 0.0061488, 0.0001, 0.01   )
    Bc_m1 = RooRealVar ( "Bc_m1"   , "Bc_m1", 6.0, 5.0 , 7.0 )
    x_zn = RooRealVar ( "x_zn"   , "x_zn", 5.3287, 4.5 , 5.6 )
    
    pdfA    = RooFormulaVar ("pdfA"    , "@0-@1 > 0 ? (@0-@1)^(@2) : 0", RooArgList (x,Bc_m1, B_be ))
    pdfB    = RooFormulaVar ("pdfB"    , "1/(1+exp((@1-@0)/@2))", RooArgList (x,E_m, E_s ))
    pdfC    = RooFormulaVar ("pdfC"    , "1/(1+exp((@1-@0)/@2))", RooArgList (x,E_m_2, E_s_2 ))
    
    pdfG    = RooGenericPdf ("pdfG"    , " (@0-@1)^(@2) * (1+@3*(@0-6.2)) *@4 * @5", RooArgList (x,Bc_m, B_al,C,pdfB,pdfC))

    Bc_m.setConstant(True)
    B_al.setConstant(True)
    C.setConstant(True)
    f.setConstant(True)
    E_s.setConstant(True)
    E_s_2.setConstant(True)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )
    C.setConstant(False)
    f.setConstant(False)
    B_al.setConstant(False)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )
    E_s.setConstant(False)
    E_s_2.setConstant(False)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7), RooFit.Range('SBL,SBR'), RooFit.Save() )
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7), RooFit.Range('SBL,SBR'), RooFit.Save() )

    rrr.Print()

# making plots

binN = 60
xframe = x.frame(binN)

xframe.SetTitle("fon BKPI");
dh.plotOn(xframe, RooFit.MarkerSize(0.6), RooFit.Name('data'))

if FITTING :
	pdfG.plotOn(xframe, RooFit.Name('pdfG'),RooFit.LineColor(kRed),
                RooFit.Normalization(dh.reduce(stringforsidebands).numEntries(), 2), RooFit.Range('SBL,SBR') )

	H_FITRES = TH1F('H_FITRES', 'c2df-ndof-chi2-prob', 4, 1, 5)
	FIT_c2df    = xframe.chiSquare(rrr.floatParsFinal().getSize() )
	FIT_ndof    = (binN - rrr.floatParsFinal().getSize())
	FIT_chi2    = FIT_c2df * FIT_ndof; H_FITRES[3] = FIT_chi2;
	H_FITRES[1] = FIT_c2df; H_FITRES[2] = FIT_ndof
	H_FITRES[4] = TMath.Prob(FIT_chi2, FIT_ndof)

	dh.plotOn(xframe, RooFit.MarkerSize(0.6), RooFit.Name('data'), RooFit.Range('SBL,SBR'))

	print ('chisq %f, ndf %i, binN %i, chisqndf %f, prob %f' % (H_FITRES[3], H_FITRES[2], binN, H_FITRES[1], H_FITRES[4]))



canvas = TCanvas("canvas")
canvas.cd()
xframe.Draw()
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.05)
latex.DrawLatex(0.1,0.85, "Bc -> BKPI")

if FITTING :
    latex.SetTextSize(0.04)
    latex.DrawText(0.1,0.55,'C = ' + str(pdfG.getParameters(dh)["C"].getValV()))
    latex.DrawText(0.1,0.35,'Bc_m = ' + str(pdfG.getParameters(dh)["Bc_m"].getValV()))
    latex.DrawText(0.1,0.50,'B_al = ' + str(pdfG.getParameters(dh)["B_al"].getValV()))
    latex.DrawText(0.1,0.45,'chi/ndf = ' + str(H_FITRES[3]) + "/" + str(H_FITRES[2]))


canvas.SaveAs('Bc_mass.png')

