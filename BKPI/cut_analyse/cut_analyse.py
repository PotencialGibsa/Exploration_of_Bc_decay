#Coomon cuts from XBframe and MySelector
#muon_pt > 3
#muons eta < 2.4
#K1_pt > 1.2
#if (p4k1.DeltaR(p4mum_0c) < 0.01) continue;
#if (p4k1.DeltaR(p4mup_0c) < 0.01) continue;
#eta K1 < 3
#K2_pt > 1
#K2_eta < 3
#PI3_pt > 0.4
#PI3_eta < 3
#Bc_vtxprob > 0.05
#B_vtxprob > 0.02

from ROOT import *
gROOT.SetBatch(True)
import datetime
#import matplotlib.pyplot as plt
import numpy as np
import random
from array import array
import os


def cut_to_string(branch, cut): #branch :: string , cut :: array of arrays like [">", 0.25]
    cut_string = ""
    for cut_ in cut:
        cut_string += " && ( " + branch + " " + cut_[0] + " " + str(cut_[1]) + " )" 
    return cut_string

def del_cut (branch, cuts):
    if cuts.find(branch) != -1:
        cut_string_to_replace = cuts[cuts.find(branch)-5 : cuts.find(branch)-5 + cuts[cuts.find(branch)-5 : len(cuts)].find(")")+1]
        if cut_string_to_replace != "" :
            cuts = cuts.replace(cut_string_to_replace, "")
   
    return cuts

def cut_proccess(chain_signal, chain_background, cuts):

    h1 = TH1F('Bc', 'Bc', 68, 5.915, 6.595) #5.9
    h2 = TH1F('Bc_WS', 'Bc_WS', 68, 5.915, 6.595)

    gen_matching_cuts = " && gen_muonP_delta < 0.2 && gen_muonN_delta < 0.2 && gen_kaon1_delta < 0.2 && gen_kaon2_delta < 0.2 && gen_pion1_delta < 0.2"

    signalN = chain_signal.Draw("Bc_mas1 >> Bc", "(Bc_mas1 > 6.179 && Bc_mas1 < 6.32) " + gen_matching_cuts + cuts) #window on Bc mass
    #backgroundN = chain_background.Draw("Bc_mass >> Bc_WS", "(Bc_mass > 6.179 && Bc_mass < 6.32) " + cuts) #the same window on data
    backgroundN = chain_background.Draw("Bc_mass >> Bc_WS", "((Bc_mass < 6.179 && Bc_mass > (6.038)) || (Bc_mass < (6.461) && Bc_mass > 6.32)) " + cuts) # the halfsum of sidebands on data
    backgroundN = backgroundN * 0.5

    f = signalN/(463.0/13.0 + 4.0 * backgroundN**0.5 + 5.0 * (25.0 + 8.0 * backgroundN**0.5 + 4.0 * backgroundN)**0.5) 

    return signalN, backgroundN, f


def get_branch_range(branch, cuts):
    True


def plots_make(data, range_cuts, cut):
    True


def plot_draw(x_range, y_range, title):
    g = TGraphErrors(len(x_range), x_range, y_range )
    g.SetTitle(title)
    g.SetLineColor( 1 )
    g.SetLineWidth( 2 )
    g.SetMarkerColor( 1 )
    g.SetMarkerStyle( 21 )
    #g.Draw("AP")
    return g


if __name__ == "__main__":

    REMEMBER_CUTS = True
    chain_signal = TChain("mytree")
    chain_background = TChain("mytree")

    chain_signal.Add("/eos/home-d/dshmygol/SimpleFile_MC_Bc_in_BKPI_v7_trigger_matching.root") #gen matching 

    #chain_background.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI_wrong_charge/SimpleFile_Bc_in_BKPI_wrong_charge_trigger_matching_2016_part_2.root')
    #chain_background.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI_wrong_charge/SimpleFile_Bc_in_BKPI_wrong_charge_trigger_matching_2016_part_1.root')
    #chain_background.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI_wrong_charge/SimpleFile_Bc_in_BKPI_wrong_charge_trigger_matching_2017_part_2.root')
    #chain_background.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI_wrong_charge/SimpleFile_Bc_in_BKPI_wrong_charge_trigger_matching_2017_part_1.root')
    #chain_background.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI_wrong_charge/SimpleFile_Bc_in_BKPI_wrong_charge_trigger_matching_2018_part_2.root')
    #chain_background.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI_wrong_charge/SimpleFile_Bc_in_BKPI_wrong_charge_trigger_matching_2018_part_1.root')
    
    chain_background.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2016_part_1.root')
    chain_background.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2016_part_2.root')
    chain_background.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2017_part_1.root')
    chain_background.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2017_part_2.root')
    chain_background.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2018_part_1.root')
    chain_background.Add('/eos/home-d/dshmygol/SimpleFile_Bc_in_BKPI/SimpleFile_Bc_in_BKPI_v6_trigger_matching_2018_part_2.root')
    
    #structure : {branch_name::string : [ cut::array]}
    #cut :: array of arrays, for example : [[">", 0.25] , ["<", 0.5] ]
    #means variable > 0.25 && variable < 0.5
    # sign is always first, the number is second
    cuts = {
                "K1_pt" :   [
                            [[">", 1.2]],
                            [[">", 1.25]],
                            [[">", 1.3]],
                            [[">", 1.35]],
                            [[">", 1.5]],
                            [[">", 1.75]],
                            [[">", 2.0]],
                            [[">", 2.25]],
                            [[">", 2.5]],
                            [[">", 2.75]],
                            [[">", 3.0]],
                            [[">", 3.25]],
                            [[">", 3.5]],
                            [[">", 3.75]],
                            [[">", 4.0]]
                            ]
            ,      
            
                "K2_pt" :   [ #check XbFrame
                            [[">", 1.0]],
                            [[">", 1.1]],
                            [[">", 1.2]],
                            [[">", 1.25]],
                            [[">", 1.3]],
                            [[">", 1.4]],
                            [[">", 1.5]],
                            [[">", 1.75]],
                            [[">", 2.0]],
                            [[">", 2.25]],
                            [[">", 2.5]],
                            [[">", 2.75]],
                            [[">", 3.0]],
                            [[">", 3.25]],
                            [[">", 3.5]],
                            [[">", 3.75]],
                            [[">", 4.0]]
                            ]
            ,
                "PI3_pt" :  [ #smaller distribute 0.1
                            [[">", 0.4]],
                            [[">", 0.5]],
                            [[">", 0.6]],
                            [[">", 0.75]],
                            [[">", 0.8]],
                            [[">", 0.9]],
                            [[">", 1.0]],
                            [[">", 1.25]],
                            [[">", 1.5]],
                            [[">", 1.75]],
                            [[">", 2.0]],
                            [[">", 2.25]],
                            [[">", 2.5]],
                            [[">", 2.75]],
                            [[">", 3.0]],
                            [[">", 3.25]],
                            [[">", 3.5]],
                            [[">", 3.75]],
                            [[">", 4.0]]
                            ]
            ,
                "Bc_pt" :   [
                            [[">", 5.0]],
                            [[">", 6.0]],
                            [[">", 7.0]],
                            [[">", 8.0]],
                            [[">", 8.5]],
                            [[">", 9.0]],
                            [[">", 10.0]],
                            [[">", 12.0]],
                            [[">", 15.0]],
                            [[">", 18.0]],
                            [[">", 20.0]],
                            [[">", 23.0]],
                            [[">", 26.0]],
                            [[">", 30.0]]
                            ]
            ,   
               "B_pt" :   [
                            [[">", 5.0]],
                            [[">", 5.5]],
                            [[">", 6.0]],
                            [[">", 6.5]],
                            [[">", 7.0]],
                            [[">", 7.5]],
                            [[">", 8.0]],
                            [[">", 8.5]],
                            [[">", 9.0]],
                            [[">", 10.0]],
                            [[">", 10.5]],
                            [[">", 11.0]],
                            [[">", 12.0]],
                            [[">", 15.0]],
                            [[">", 18.0]],
                            [[">", 20.0]]
                            ]
            ,
                "Bc_pvdistsignif2" : [
                                       [[">", 0]],
                                       [[">", 0.2]],
                                       [[">", 0.4]],
                                       [[">", 0.6]],
                                       [[">", 0.8]],
                                       [[">", 1.0]],
                                       [[">", 2.0]],
                                       [[">", 2.2]],
                                       [[">", 2.4]],
                                       [[">", 2.6]],
                                       [[">", 2.8]],
                                       [[">", 3.0]],
                                       [[">", 4.0]],
                                       [[">", 4.5]],
                                       [[">", 5.0]],
                                       [[">", 5.5]],
                                       [[">", 6.0]],
                                       [[">", 7.0]],
                                       [[">", 8.0]],
                                       [[">", 9.0]],
                                       [[">", 10.0]], 
                                    ]
            ,
                "B_Bcdistsignif2"   : [
                                       [[">", 0]],
                                       [[">", 0.2]],
                                       [[">", 0.4]],
                                       [[">", 0.6]],
                                       [[">", 0.7]],
                                       [[">", 0.8]],
                                       [[">", 0.9]],
                                       [[">", 1.0]],
                                       [[">", 2.0]],
                                       [[">", 3.0]],
                                       [[">", 4.0]],
                                       [[">", 5.0]],
                                       [[">", 6.0]],
                                       [[">", 7.0]],
                                       [[">", 8.0]],
                                       [[">", 9.0]],
                                       [[">", 10.0]], 
                                    ]
            ,
                "B_Bccos2"  :   [
                                [[">", 0.99]],
                                [[">", 0.991]],
                                [[">", 0.992]],
                                [[">", 0.993]],
                                [[">", 0.994]],
                                [[">", 0.995]],
                                [[">", 0.996]],
                                [[">", 0.997]],
                                [[">", 0.9975]],
                                [[">", 0.998]],
                                [[">", 0.9985]],
                                [[">", 0.999]],
                                ]
            ,
                "Bc_pvcos2" :   [
                                [[">", 0.99]],
                                [[">", 0.991]],
                                [[">", 0.992]],
                                [[">", 0.993]],
                                [[">", 0.994]],
                                [[">", 0.995]],
                                [[">", 0.996]],
                                [[">", 0.997]],
                                [[">", 0.998]],
                                [[">", 0.999]],
                                ]  
            ,
                "B_vtxprob" :   [
                                [[">", 0.02]],
                                [[">", 0.03]],
                                [[">", 0.04]],
                                [[">", 0.042]],
                                [[">", 0.044]],
                                [[">", 0.046]],
                                [[">", 0.048]],
                                [[">", 0.05]],
                                [[">", 0.06]],
                                [[">", 0.07]],
                                [[">", 0.08]],
                                [[">", 0.09]],
                                [[">", 0.1]],
                                [[">", 0.11]],
                                [[">", 0.12]],
                                [[">", 0.13]],
                                [[">", 0.14]],
                                [[">", 0.15]],
                                ]
            ,
                "Bc_vtxprob" :  [   
                                    [[">", 0.05]],
                                    [[">", 0.06]],
                                    [[">", 0.07]],
                                    [[">", 0.08]],
                                    [[">", 0.09]],
                                    [[">", 0.1]],
                                ]
            ,
                 "K1_ips" :  [      
                                    [[">", 0.01]],
                                    [[">", 0.03]],
                                    [[">", 0.05]],
                                    [[">", 0.07]],
                                    [[">", 0.09]],
                                    [[">", 0.1]],
                                    [[">", 0.2]],
                                    [[">", 0.3]],
                                    [[">", 0.4]],
                                    [[">", 0.5]],
                                    [[">", 0.6]],
                                    [[">", 0.7]],
                                    [[">", 0.8]],
                                    [[">", 0.9]],
                                    [[">", 1.0]],
                                    [[">", 1.1]],
                                    [[">", 1.2]],
                                    [[">", 1.3]],
                                    [[">", 1.4]],
                                    [[">", 1.5]],
                                    [[">", 1.8]],
                                    [[">", 2]],
                                    [[">", 2.25]],
                                    [[">", 2.5]],
                                    [[">", 2.75]],
                                    [[">", 3]],
                                    [[">", 3.5]],
                                    [[">", 4]],
                                ]
            ,
                "K2_ips" :  [       
                                    [[">", 0.01]],
                                    [[">", 0.03]],
                                    [[">", 0.04]],
                                    [[">", 0.05]],
                                    [[">", 0.06]],
                                    [[">", 0.07]],
                                    [[">", 0.09]],
                                    [[">", 0.1]],
                                    [[">", 0.2]],
                                    [[">", 0.3]],
                                    [[">", 0.4]],
                                    [[">", 0.5]],
                                    [[">", 0.6]],
                                    [[">", 0.7]],
                                    [[">", 0.8]],
                                    [[">", 0.9]],
                                    [[">", 1.0]],
                                    [[">", 1.1]],
                                    [[">", 1.2]],
                                    [[">", 1.3]],
                                    [[">", 1.4]],
                                    [[">", 1.5]],
                                ]
            ,
                "PI3_ips" :  [      
                                    [[">", 0.01]],
                                    [[">", 0.03]],
                                    [[">", 0.05]],
                                    [[">", 0.07]],
                                    [[">", 0.09]],
                                    [[">", 0.1]],
                                    [[">", 0.2]],
                                    [[">", 0.3]],
                                    [[">", 0.4]],
                                    [[">", 0.5]],
                                    [[">", 0.6]],
                                    [[">", 0.7]],
                                    [[">", 0.8]],
                                    [[">", 0.9]],
                                    [[">", 1.0]],
                                    [[">", 1.05]],
                                    [[">", 1.1]],
                                    [[">", 1.15]],
                                    [[">", 1.2]],
                                    [[">", 1.3]],
                                    [[">", 1.4]],
                                    [[">", 1.5]],
                                ]
        #raznost massi B+ i tablicnoy
        }
    # Repeating the main cycle
    NUM_repeats = 10

    cuts_to_remember = ""
    for N_rep in range(NUM_repeats):

        branches_iterate = cuts.keys()
        random.shuffle(branches_iterate)
        for branch in branches_iterate:
            cuts_axes = array( 'd' )
            signalsN = array( 'd' )
            backgroundsN = array( 'd' )
            fs = array( 'd' )
            for cut in cuts[branch]:
                cut_string = cut_to_string(branch, cut)
                
                #Create the cut with the memory about the previous good cuts
                if REMEMBER_CUTS:
                    #Replace the already existed cuts
                    cuts_to_remember_changed = del_cut(branch , cuts_to_remember)
                    cut_string += cuts_to_remember_changed
                    

                #Add the default cuts like trigger cuts
                cut_string += " && (KA1_dptTRG < 0.1) && (KA1_drTRG < 0.1)"
                cut_string += " && (B_mass < 5.33) && (B_mass > 5.23)"
                #strange cut on BandK_mass to check the decays in B0s or without B0s
                #cut_string += " && (BandK_mass < 5.9)"
                signalN, backgroundN, f = cut_proccess( chain_signal , chain_background, cuts = cut_string)
                #print("N_rep = " , N_rep, " branch = ", branch, " cut = ", cut_to_string(branch, cut) , " signalN = ", signalN, " backgroundN = ", backgroundN, " f = ", f)
                #Attention work only with cut[0] of length 1
                cuts_axes.append(cut[0][1])
                #
                signalsN.append(signalN)
                backgroundsN.append(backgroundN)
                fs.append(f)

            #find max f and enter this to remember this cuta
            max_f = max(fs)
            index_max_f = fs.index(max_f)
            if REMEMBER_CUTS:
                cuts_to_remember = del_cut(branch, cuts_to_remember) + cut_to_string(branch, cuts[branch][index_max_f])
            #

            if not os.path.isdir("./branch_cuts/" + branch):
                os.mkdir("./branch_cuts/" + branch)
            
            canv = TCanvas("canvas", "canvas", 500,500, 2500,1500)
            canv_signal = plot_draw(cuts_axes, signalsN, title = branch + "_signal_from_MC")
            canv_signal.Draw("AP")
            canv.SaveAs("./branch_cuts/" + branch + "/" + branch + "_signal_from_MC.png")

            canv = TCanvas("canvas", "canvas", 500,500, 2500,1500)
            canv_background = plot_draw(cuts_axes, backgroundsN, title = branch + "_background_from_data")
            canv_background.Draw("AP")
            canv.SaveAs("./branch_cuts/" + branch + "/" + branch + "_background_from_data.png")

            canv = TCanvas("canvas", "canvas", 500,500, 2500,1500)
            canv_f = plot_draw(cuts_axes, fs, title = branch + "_f_value")
            canv_f.Draw("AP")
            canv.SaveAs("./branch_cuts/" + branch + "/" + branch + "_f_value.png")

            if REMEMBER_CUTS:
        
                canv = TCanvas("canvas", "canvas", 500,500, 2500,1500)
                canv_f = plot_draw(cuts_axes, fs, title = branch + "_f_value" + cuts_to_remember)
                canv_f.Draw("AP")
                canv.SaveAs("./branch_cuts/_f_value.png")
                print("Cuts = " + cuts_to_remember)

