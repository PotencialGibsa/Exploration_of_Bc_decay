from ROOT import *
gROOT.SetBatch(True)
import datetime
import glob

name = "Bc_pvcos2"

fr = TFile(name + '_dataset.root', 'read')
dh_data = fr.dh.reduce('1>0')
fr.Close()

fr = TFile(name + '_dataset_MC.root', 'read')
dh_MC = fr.dh.reduce('1>0')
fr.Close()

x = RooRealVar ( "x"      , "x"       , 0.995,  1)
binN = 40
x_data_frame = x.frame(binN)
x_MC_frame = x.frame(binN)

dh_data.plotOn(x_data_frame, RooFit.MarkerSize(1), RooFit.Name('data'))
dh_MC.plotOn(x_MC_frame, RooFit.MarkerSize(1), RooFit.Name('data'))


canvas = TCanvas("canvas", name, 500,500, 2500,1500)
canvas.Divide(2)
canvas.cd(1)
x_data_frame.SetTitle("Data")
x_data_frame.Draw()
canvas.cd(2)
x_MC_frame.SetTitle("MC")
x_MC_frame.Draw()

canvas.SaveAs('var_analyse/' + name + 'step2.png')
