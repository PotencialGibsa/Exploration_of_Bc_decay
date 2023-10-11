from ROOT import *
#from variables import *
import glob
import numpy as np
import math
from Constants import *

def sqrt(var):
      return math.sqrt(abs(var))

def DetachSignificance2(vtx, vtxE1, vtxE2):
    return sqrt( vtx.X()**2 / (vtxE1.X()**2 + vtxE2.X()**2) + vtx.Y()**2 / (vtxE1.Y()**2 + vtxE2.Y()**2))

def DetachSignificance3(vtx, vtxE1, vtxE2):
    return sqrt( vtx.X()**2 / (vtxE1.X()**2 + vtxE2.X()**2) + vtx.Y()**2 / (vtxE1.Y()**2 + vtxE2.Y()**2) + vtx.Z()**2 / (vtxE1.Z()**2 + vtxE2.Z()**2))

def DirectionCos2 (v1, v2):
    r1 = sqrt(v1.X()**2 + v1.Y()**2)
    r2 = sqrt(v2.X()**2 + v2.Y()**2)
    return ( v1.X() * v2.X() + v1.Y() * v2.Y() ) / (r1*r2 + 0.0000001)

def DirectionCos3 (v1, v2):
    r1 = sqrt(v1.X()**2 + v1.Y()**2 + v1.Z()**2)
    r2 = sqrt(v2.X()**2 + v2.Y()**2 + v2.Z()**2)
    return ( v1.X() * v2.X() + v1.Y() * v2.Y() + v1.Z() * v2.Z() ) / (r1*r2 + 0.0000001)

def DirectionChi22 (vtx0, vtx0E, vtx1, vtx1E, P, PE):
    dvtx    = vtx1 - vtx0
    dvtxE   = TVector3( sqrt(vtx0E.X()**2 + vtx1E.X()**2), sqrt(vtx0E.Y()**2 + vtx1E.Y()**2), sqrt(vtx0E.Z()**2 + vtx1E.Z()**2) )
    Pscaled = P * (dvtx.Mag() / P.Mag())
    PscaledE= PE * (dvtx.Mag() / P.Mag())
    return DetachSignificance2 (Pscaled - dvtx, PscaledE, dvtxE)

def DirectionChi23 (vtx0, vtx0E, vtx1, vtx1E, P, PE):
    dvtx    = vtx1 - vtx0 ## vertex difference
    dvtxE   = TVector3( sqrt(vtx0E.X()**2 + vtx1E.X()**2), sqrt(vtx0E.Y()**2 + vtx1E.Y()**2), sqrt(vtx0E.Z()**2 + vtx1E.Z()**2) ) ## its error
    Pscaled = P * (dvtx.Mag() / P.Mag()) ## scaled momentum to be the same length as vertex difference
    PscaledE= PE * (dvtx.Mag() / P.Mag()) ## its error
    return DetachSignificance3 (Pscaled - dvtx, PscaledE, dvtxE)


def TH1NormBinWidth(hh):
    _htemp = hh.Clone()
    _int = _htemp.Integral()
    for _i in range(1, hh.GetNbinsX()+1):
        _htemp[_i] = _htemp[_i] / _htemp.GetBinWidth(_i)
    #
    _nam = hh.GetName() + '_no'
    _htemp.SetName('_nam')
    _htemp.Scale(_int / _htemp.Integral())
    return _htemp

def stri(n):
    if n>=0 and n<10:
        return '0'+str(n)
    else:
        return str(n)
years = [6,7,8]
for year in years:
    MyFileNames = glob.glob('/eos/home-d/dshmygol/Bc_Bfinder/201'+str(year)+'_Bc_in_JpsiPI_v6_trigger_matching/*/Charmonium/*/*/000?/BFinder*.root')
    print('The number of files is %i for %s' % (len(MyFileNames), iter))
    print('\nAdding the files to the chain')
    chain = TChain('wztree')
    Name_count = 0
    for name in MyFileNames[0:1000]:
        if Name_count % 500 == 0:
            print("Added: ", Name_count)
        chain.Add(name)
        Name_count+=1
    _fileOUT = '/eos/home-d/dshmygol/SimpleFile_Bc_in_JpsiPi/SimpleFile_Bc_in_JpsiPI_v6_trigger_matching_201'+str(year)+'.root'
    fileOUT = TFile (_fileOUT, "recreate")
    mytree = TTree("mytree","mytree")

    nevt = chain.GetEntries()
    print('Number of events: %s' % (nevt))
    my_vars = [
        'SAMEEVENT',
        #Mass BK
        'B_mass',
        'B_pt', 'B_eta',
        'B_pvdistsignif2', 'B_pvdistsignif3',
        'B_pvdist', 'B_vtxprob',
        'B_pvcos2', 'B_pvcos3',
        #J/psi
        'J_psi_mass', 'J_psi_prob',
        'J_psi_pt',
        #Pion1
        'PI1_pt', 'PI1_ips', 'PI1_chrg','PI1CV_pt',
        #muon1
        'mu1_pt','mu1CV_pt',
        #muon2
        'mu2_pt', 'mu2CV_pt',

		'PI1_drTRG', 'PI1_dptTRG',
        ]

    for _var in my_vars:
        exec(_var+'=np.zeros(1, dtype=float)')
        exec('mytree.Branch("' + _var + '"' + ' '*(25-len(_var)) + ',' + _var + ' '*(25-len(_var)) + ', "'+ _var + '/D")')

    ### finished with deploying vars for SimpleFile }}}

    # P4 ROOT Lorentz vector
    PI1P4,PI1P4_CV, JPSIP4, BP4  = [TLorentzVector() for i in range(4)]

    BBB = 0
    # looping all over the events
    for evt in range(nevt):
        if chain.GetEntry(evt) <= 0:
            break
        # securing that _nCand is set to 0 for every event before assign the real chain value
        _nCand = 0
        _nCand = chain.nCand
        if _nCand < 1:
            continue
        for cand in range(_nCand):
            #it only works when we set variable with [0] by its side
            # assigning PV info
            PV = TVector3(chain.PV_becos_XX[cand], chain.PV_becos_YY[cand], chain.PV_becos_ZZ[cand])
            PVE = TVector3(sqrt(chain.PV_becos_EX[cand]), sqrt(chain.PV_becos_EY[cand]), sqrt(chain.PV_becos_EZ[cand]))

            # assigning Pion info
            PI1P4.SetXYZM (chain.PI1_px[cand], chain.PI1_py[cand], chain.PI1_pz[cand], PDG_PION_MASS)
            PI1P4_CV.SetXYZM (chain.PI1_px_CV[cand], chain.PI1_py_CV[cand], chain.PI1_pz_CV[cand], PDG_PION_MASS)
            PI1_pt[0] = PI1P4.Pt()
            PI1CV_pt[0] = PI1P4_CV.Pt()
            PI1_ips[0] = chain.PI1_ips[cand]
            PI1_chrg[0] = chain.PI1_ch[cand]

            # assigning muon1 info
            mu1P4.SetXYZM (chain.B_mu_px1[cand], chain.B_mu_py1[cand], chain.B_mu_pz1[cand], PDG_MUON_MASS)
            mu1P4_CV.SetXYZM (chain.B_mu_px1_cjp[cand], chain.B_mu_py1_cjp[cand], chain.B_mu_pz1_cjp[cand], PDG_MUON_MASS)
            mu1_pt[0] = mu1P4.Pt()
            mu1CV_pt[0] = mu1P4_CV.Pt()

            # assigning muon2 info
            mu2P4.SetXYZM (chain.B_mu_px2[cand], chain.B_mu_py2[cand], chain.B_mu_pz2[cand], PDG_MUON_MASS)
            mu2P4_CV.SetXYZM (chain.B_mu_px2_cjp[cand], chain.B_mu_py2_cjp[cand], chain.B_mu_pz2_cjp[cand], PDG_MUON_MASS)
            mu2_pt[0] = mu2P4.Pt()
            mu2CV_pt[0] = mu2P4_CV.Pt()


            #    # assigning J/psi info
            J_psi_mass[0] = chain.B_J_mass[cand]
            J_psi_prob[0] = chain.B_J_Prob[cand]
            JPSIP4.SetXYZM(chain.B_J_px[cand], chain.B_J_py[cand], chain.B_J_pz[cand], chain.B_J_mass[cand])
            J_psi_pt[0] = JPSIP4.Pt()
            # assigning info about B+
            BV = TVector3(chain.B_DecayVtxX[cand], chain.B_DecayVtxY[cand], chain.B_DecayVtxZ[cand])
            BVE = TVector3(sqrt(chain.B_DecayVtxXE[cand]), sqrt(chain.B_DecayVtxYE[cand]), sqrt(chain.B_DecayVtxZE[cand]))
            BP4.SetXYZM( chain.B_px[cand], chain.B_py[cand], chain.B_pz[cand], chain.B_mass[cand])
            BP3 = BP4.Vect()
            B_pt[0]= BP4.Pt()
            B_eta[0]= BP4.Eta()
            B_vtxprob[0] = chain.B_Prob[cand]
            B_mass[0]= chain.B_mass[cand]
            B_pvdistsignif2[0] = DetachSignificance2( BV - PV, BVE, PVE)
            B_pvdistsignif3[0] = DetachSignificance3( BV - PV, BVE, PVE)
            B_pvdist[0] = (BV - PV).Mag()
            B_pvcos2[0] = DirectionCos2 ( BV - PV, BP3)
            B_pvcos3[0] = DirectionCos3 ( BV - PV, BP3); PI1_drTRG[0] = chain.PI1_drTRG[cand];PI1_dptTRG[0] = chain.PI1_dptTRG[cand]


            SAMEEVENT[0] = 0;
            if (BBB > -1):
                SAMEEVENT[0] = 1
            BBB = 1; ## the next candidate is not the 1st in event

            # fill _var branches in the tree if passes preselection cuts
            mytree.Fill()
                    #
        BBB = -1 ## when loop over candidates in event is finished, set this to -1, so the next candidate has SAMEEVENT=0
        if (evt % 10000 == 0):
            print('Running... Now in event %i / %i ' %(evt, nevt))



    print('Total entries stored in MyTree: %s' %(mytree.GetEntries()))

    #print("N_cand_raw = ", N_cand_raw, " N_cand_Bvtp = ", N_cand_Bvtp, " N_cand_Bcvtp = ", N_cand_Bcvtp)

    fileOUT.Write()
