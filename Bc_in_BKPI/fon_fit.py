from ROOT import *
gROOT.SetBatch(True)
import datetime
#####The fitting of the fon from the decay: Bc -> BKPI

FITTING = True

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
x.setRange("SBL", 5.92, 6.179)
x.setRange("SBR", 6.32, 6.6)
stringforsidebands = 'Bc_mass<6.179||Bc_mass>6.32'

varset = RooArgSet(x)
'''
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

    PI3_pt = ch.PI3_pt
    KA2_pt = ch.K2_pt
    B_pt = ch.B_pt
    #KA1_dptTRG = ch.KA1_dptTRG
    #if KA1_dptTRG > 10: continue

    #if B_pt < 1 :continue
    #if PI3_pt < 1 :continue
    #if KA2_pt <1 :continue
    #Cut the signal
    #if (Bc_mass > 6.179) and (Bc_mass < 6.32) : continue

    """
    KA1_dptTRG = ch.KA1_dptTRG
    if KA1_dptTRG > 10: continue
    PI3_chrg = ch.PI3_chrg
    if PI3_chrg < 0 : continue
    PI3_pt = ch.PI3_pt
    if PI3_pt < 1 : continue
    B_pvcos2 = ch.B_pvcos2
    if B_pvcos2 < 0.99 : continue
    PIandK_mass = ch.PIandK_mass
    #if PIandK_mass < 1 :continue
    BandK_mass = ch.BandK_mass
    if BandK_mass > 5.9 : continue
    """

    x.setVal(Bc_mass)
    dh.add(varset)
    N_cand += 1
print("The number of candidates on the picture is : ", N_cand)

fi = TFile('mydataset_all.root', 'recreate')
dh.Write()
fi.Close()
'''
fr = TFile('mydataset_all.root', 'read')
dh = fr.dh.reduce('1>0')
fr.Close()
#############################################
if FITTING :
    PDG_B = 5.27934
    PDG_K = 0.493677
    PDG_PI = 0.13957061
    B       = RooRealVar ( "B"      , "B"       ,1.37983*10**7, 0, 14000000) #all 14418848
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
    #b    = RooRealVar ( "b"   , "b", 12.6, -100, 100   )
    # pdfG    = RooGenericPdf ("pdfG"    , "@4*@0/@7*(@0-@1)^(@2) * (1+@3*(@0-@6))^@5", RooArgList (x,Bc_m,  B_al,C,  B , H_al, Bc_m1, x_zn) )
    pdfA    = RooFormulaVar ("pdfA"    , "@0-@1 > 0 ? (@0-@1)^(@2) : 0", RooArgList (x,Bc_m1, B_be ))
    pdfB    = RooFormulaVar ("pdfB"    , "1/(1+exp((@1-@0)/@2))", RooArgList (x,E_m, E_s ))
    pdfC    = RooFormulaVar ("pdfC"    , "1/(1+exp((@1-@0)/@2))", RooArgList (x,E_m_2, E_s_2 ))
    #pdfG    = RooGenericPdf ("pdfG"    , "(1.0-@4) * (@0-@1)^(@2) * (1+@3*(@0-6.2)) + @4*@5", RooArgList (x,Bc_m, B_al,C,f,pdfA))
    #pdfG    = RooGenericPdf ("pdfG"    , " (@0-@1)^(@2) * (1+@3*(@0-6.2)) * @4", RooArgList (x,Bc_m, B_al,C,pdfB))
    pdfG    = RooGenericPdf ("pdfG"    , " (@0-@1)^(@2) * (1+@3*(@0-6.2)) *@4 * @5", RooArgList (x,Bc_m, B_al,C,pdfB,pdfC))

    Bc_m.setConstant(True)
    B_al.setConstant(True)
    C.setConstant(True)
    #D.setConstant(True)
    f.setConstant(True)
    E_s.setConstant(True)
    E_s_2.setConstant(True)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )
    C.setConstant(False)
    f.setConstant(False)
    B_al.setConstant(False)
    # Bc_m.setConstant(False)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )
    E_s.setConstant(False)
    E_s_2.setConstant(False)
    #D.setConstant(False)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7), RooFit.Range('SBL,SBR'), RooFit.Save() )
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7), RooFit.Range('SBL,SBR'), RooFit.Save() )

    rrr.Print()

binN = 60
xframe = x.frame(binN)

xframe.SetTitle("fon BKPI");
dh.plotOn(xframe, RooFit.MarkerSize(0.6), RooFit.Name('data'), RooFit.Range('SBL,SBR'))

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
    # latex.DrawText(0.1,0.8,'B = ' + str(pdfG.getParameters(dh)["B"].getValV()))
    latex.DrawText(0.1,0.55,'C = ' + str(pdfG.getParameters(dh)["C"].getValV()))
    # latex.DrawText(0.1,0.60,'x_zn = ' + str(pdfG.getParameters(dh)["x_zn"].getValV()))
    # latex.DrawText(0.1,0.40,'H_al = ' + str(pdfG.getParameters(dh)["H_al"].getValV()))
    latex.DrawText(0.1,0.35,'Bc_m = ' + str(pdfG.getParameters(dh)["Bc_m"].getValV()))
    # latex.DrawText(0.1,0.30,'Bc_m1 = ' + str(pdfG.getParameters(dh)["Bc_m1"].getValV()))
    latex.DrawText(0.1,0.50,'B_al = ' + str(pdfG.getParameters(dh)["B_al"].getValV()))
    latex.DrawText(0.1,0.45,'chi/ndf = ' + str(H_FITRES[3]) + "/" + str(H_FITRES[2]))


canvas.SaveAs('Bc_mass_fit_fon1.png')

'''
    B.setConstant(True)
    Bc_m1.setConstant(True)
    B_al.setConstant(True)
    H_al.setConstant(True)
    Bc_m.setConstant(True)
    C.setConstant(True)
    x_zn.setConstant(True)
    #rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )
    B.setConstant(True)
    Bc_m1.setConstant(True)
    B_al.setConstant(True)
    H_al.setConstant(True)
    Bc_m.setConstant(True)
    C.setConstant(False)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )
    B.setConstant(True)
    Bc_m1.setConstant(True)
    B_al.setConstant(False)
    H_al.setConstant(True)
    Bc_m.setConstant(True)
    C.setConstant(True)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )
    B.setConstant(True)
    Bc_m1.setConstant(True)
    B_al.setConstant(True)
    H_al.setConstant(False)
    Bc_m.setConstant(True)
    C.setConstant(True)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )
    B.setConstant(True)
    Bc_m1.setConstant(False)
    B_al.setConstant(True)
    H_al.setConstant(True)
    Bc_m.setConstant(True)
    C.setConstant(True)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )
    B.setConstant(True)
    Bc_m1.setConstant(True)
    B_al.setConstant(True)
    H_al.setConstant(True)
    Bc_m.setConstant(False)
    C.setConstant(True)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )
    x_zn.setConstant(False)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )
    B.setConstant(False)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )
    B.setConstant(False)
    Bc_m1.setConstant(False)
    B_al.setConstant(False)
    H_al.setConstant(False)
    Bc_m.setConstant(False)
    C.setConstant(False)
    rrr = pdfG.fitTo( dh, RooFit.NumCPU(7),RooFit.PrintLevel(-1), RooFit.Range('SBL,SBR'), RooFit.Save() )


'''
