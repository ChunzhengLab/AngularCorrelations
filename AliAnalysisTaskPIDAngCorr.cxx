/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//--------------------------------------------------------------------------------
// CVE analysis
// Contributor: Chunzheng Wang, <chunzheng.wang@cern.ch>, Shanghai
//--------------------------------------------------------------------------------

#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <algorithm>
// ROOT classes
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
// Alice analysis base class
#include "AliAnalysisTaskSE.h"
// Alice analysis additional classes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
// Alice AOD classes
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
// Alice classes
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
// Alice MC classes
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"
// Alice "V" classes
#include "AliVParticle.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVVZERO.h"
// Alice PID classes
#include "AliAODPid.h"
#include "AliAODpidUtil.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliOADBContainer.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskPIDAngCorr.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskPIDAngCorr);

//---------------------------------------------------
AliAnalysisTaskPIDAngCorr::AliAnalysisTaskPIDAngCorr() :
  AliAnalysisTaskSE(),
  fDebug(0),
  fTrigger("kMB"),
  fPeriod("LHC10h"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fFilterBit(768),
  fNclsCut(70),
  fChi2Max(4.0),
  fChi2Min(0.1),
  fDcaCutz(3.2),
  fDcaCutxy(2.4),
  fPtMin(0.2),
  fPtMax(5.0),
  fProtonPtMax(3.0),
  fNSigmaTPCCut(4),
  fNSigmaTOFCut(4),
  fV0PtMin(0.5),
  fV0CPAMin(0.995),
  fV0RapidityMax(0.5),
  fV0DecayLengthMin(3.),
  fV0DecayLengthMax(100.),
  fV0DCAToPrimVtxMax(1.5),
  fV0DcaBetweenDaughtersMax(1.),
  fDaughtersPtMax(20.),
  fDaughtersEtaMax(0.8),
  fDaughtersTPCNclsMin(70),
  fDaughtersDCAToPrimVtxMin(0.02),
  fV0PosProtonTPCNsigma(4.),
  fV0NegPionTPCNsigma(4.),
  fV0NegProtonTPCNsigma(4.),
  fV0PosPionTPCNsigma(4.),
  IsV0DaughterUseTOF(false),
  fV0PosProtonTOFNsigma(4.),
  fV0NegPionTOFNsigma(4.),
  fV0NegProtonTOFNsigma(4.),
  fV0PosPionTOFNsigma(4.),
  fLambdaMassCut(0.005),
  fMassMean(1.115683),
  fEtaCut(0.8),
  fDedxCut(10.0),
  fAOD(nullptr),         //! aod Event
  fPIDResponse(nullptr), //! PID Handler
  fUtils(nullptr),       //! Event Selection Options
  runNumList(0),
  fRunNum(-999), // runnumber
  fOldRunNum(-999), // old runnumber
  fRunNumBin(-999), // runnumer bin, 10:139510..., 11:170387..., 15HIR:246994...
  fVzBin(-999), // vertex z bin
  fCentBin(-999), // centrality bin: 0-10
  fCent(-999), // value of centrality
  fSPDCutPU(nullptr),
  fV0CutPU(nullptr),
  fCenCutLowPU(nullptr),
  fCenCutHighPU(nullptr),
  fMultCutPU(nullptr),
  fOutputList(nullptr),
  fEvtCount(nullptr),
  fHistRunNumBin(nullptr),
  fHistPt(nullptr),
  fHistEta(nullptr),
  fHistNhits(nullptr),
  fHist2DPDedx(nullptr),
  fHistDcaXY(nullptr),
  fHistDcaZ(nullptr),
  fHistV0Pt(nullptr),
  fHistV0Eta(nullptr),
  fHistV0DcatoPrimVertex(nullptr),
  fHistV0CPA(nullptr),
  fHistV0DecayLength(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;

  for (int i = 0; i < 2; ++i) fHistCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fHistVz[i]   = nullptr;
  for (int i = 0; i < 8; ++i) fHist2DCentQA[i]     = nullptr;
  for (int i = 0; i < 2; ++i) fHist2DMultCentQA[i] = nullptr;
  for (int i = 0; i < 6; ++i) fHist2DMultMultQA[i] = nullptr;

  for (int i = 0; i < 2; ++i) fHistPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2DEtaPhi[i] = nullptr;

  for (int i = 0; i < 2; ++i) fHistLambdaPt[i]                  = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaEta[i]                 = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaPhi[i]                 = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaDcaToPrimVertex[i]     = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaCPA[i]                 = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaDecayLength[i]         = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaMass[i]                = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaPt[i]              = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaEta[i]             = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaPhi[i]             = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaDcaToPrimVertex[i] = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaCPA[i]             = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaDecayLength[i]     = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaMass[i]            = nullptr;
  for (int i = 0; i < 2; ++i) fProfileLambdaMassVsPt[i]         = nullptr;
  for (int i = 0; i < 2; ++i) fProfileAntiLambdaMassVsPt[i]     = nullptr;

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      fHistSig[i][j] = nullptr;
      fHistBkg[i][j] = nullptr;
    }
  }
}

//---------------------------------------------------
AliAnalysisTaskPIDAngCorr::AliAnalysisTaskPIDAngCorr(const char *name) :
  AliAnalysisTaskSE(name),
  fDebug(0),
  fTrigger("kMB"),
  fPeriod("LHC10h"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fFilterBit(768),
  fNclsCut(70),
  fChi2Max(4.0),
  fChi2Min(0.1),
  fDcaCutz(3.2),
  fDcaCutxy(2.4),
  fPtMin(0.2),
  fPtMax(5.0),
  fProtonPtMax(3.0),
  fNSigmaTPCCut(4),
  fNSigmaTOFCut(4),
  fV0PtMin(0.5),
  fV0CPAMin(0.995),
  fV0RapidityMax(0.5),
  fV0DecayLengthMin(3.),
  fV0DecayLengthMax(100.),
  fV0DCAToPrimVtxMax(1.5),
  fV0DcaBetweenDaughtersMax(1.),
  fDaughtersPtMax(20.),
  fDaughtersEtaMax(0.8),
  fDaughtersTPCNclsMin(70),
  fDaughtersDCAToPrimVtxMin(0.02),
  fV0PosProtonTPCNsigma(4.),
  fV0NegPionTPCNsigma(4.),
  fV0NegProtonTPCNsigma(4.),
  fV0PosPionTPCNsigma(4.),
  IsV0DaughterUseTOF(false),
  fV0PosProtonTOFNsigma(4.),
  fV0NegPionTOFNsigma(4.),
  fV0NegProtonTOFNsigma(4.),
  fV0PosPionTOFNsigma(4.),
  fLambdaMassCut(0.005),
  fMassMean(1.115683),
  fEtaCut(0.8),
  fDedxCut(10.0),
  fAOD(nullptr),         //! aod Event
  fPIDResponse(nullptr), //! PID Handler
  fUtils(nullptr),       //! Event Selection Options
  runNumList(0),
  fRunNum(-999), // runnumber
  fOldRunNum(-999), // old runnumber
  fRunNumBin(-999), // runnumer bin, 10:139510..., 11:170387..., 15HIR:246994...
  fVzBin(-999), // vertex z bin
  fCentBin(-999), // centrality bin: 0-10
  fCent(-999), // value of centrality
  fSPDCutPU(nullptr),
  fV0CutPU(nullptr),
  fCenCutLowPU(nullptr),
  fCenCutHighPU(nullptr),
  fMultCutPU(nullptr),
  fOutputList(nullptr),
  fEvtCount(nullptr),
  fHistRunNumBin(nullptr),
  fHistPt(nullptr),
  fHistEta(nullptr),
  fHistNhits(nullptr),
  fHist2DPDedx(nullptr),
  fHistDcaXY(nullptr),
  fHistDcaZ(nullptr),
  fHistV0Pt(nullptr),
  fHistV0Eta(nullptr),
  fHistV0DcatoPrimVertex(nullptr),
  fHistV0CPA(nullptr),
  fHistV0DecayLength(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;

  for (int i = 0; i < 2; ++i) fHistCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fHistVz[i]   = nullptr;
  for (int i = 0; i < 8; ++i) fHist2DCentQA[i]     = nullptr;
  for (int i = 0; i < 2; ++i) fHist2DMultCentQA[i] = nullptr;
  for (int i = 0; i < 6; ++i) fHist2DMultMultQA[i] = nullptr;

  for (int i = 0; i < 2; ++i) fHistPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2DEtaPhi[i] = nullptr;

  for (int i = 0; i < 2; ++i) fHistLambdaPt[i]                  = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaEta[i]                 = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaPhi[i]                 = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaDcaToPrimVertex[i]     = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaCPA[i]                 = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaDecayLength[i]         = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaMass[i]                = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaPt[i]              = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaEta[i]             = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaPhi[i]             = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaDcaToPrimVertex[i] = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaCPA[i]             = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaDecayLength[i]     = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaMass[i]            = nullptr;
  for (int i = 0; i < 2; ++i) fProfileLambdaMassVsPt[i]         = nullptr;
  for (int i = 0; i < 2; ++i) fProfileAntiLambdaMassVsPt[i]     = nullptr;

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      fHistSig[i][j] = nullptr;
      fHistBkg[i][j] = nullptr;
    }
  }

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//------------------------------------------------

AliAnalysisTaskPIDAngCorr::~AliAnalysisTaskPIDAngCorr()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  if (fOutputList) delete fOutputList;
}

//---------------------------------------------------

void AliAnalysisTaskPIDAngCorr::Terminate(Option_t *)
{
  // Terminate loop
  Printf("Terminate");
}

//---------------------------------------------------

void AliAnalysisTaskPIDAngCorr::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList -> SetName(GetName());
  fOutputList -> SetOwner(kTRUE);
  ////////////////////////
  // Run Number Info
  ////////////////////////
  TString runNumList10h[90] = {
    "139510","139507","139505","139503","139465","139438","139437","139360","139329","139328","139314","139310",
    "139309","139173","139107","139105","139038","139037","139036","139029","139028","138872","138871","138870",
    "138837","138732","138730","138666","138662","138653","138652","138638","138624","138621","138583","138582",
    "138579","138578","138534","138469","138442","138439","138438","138396","138364","138275","138225","138201",
    "138197","138192","138190","137848","137844","137752","137751","137724","137722","137718","137704","137693",
    "137692","137691","137686","137685","137639","137638","137608","137595","137549","137546","137544","137541",
    "137539","137531","137530","137443","137441","137440","137439","137434","137432","137431","137430","137243",
    "137236","137235","137232","137231","137162","137161"};
  TString runNumList15o[138] = {
    "246994","246991","246989","246984","246982","246948","246945","246928","246871","246870","246867","246865",
    "246864","246859","246858","246851","246847","246846","246845","246844","246810","246809","246808","246807",
    "246805","246804","246766","246765","246763","246760","246759","246758","246757","246751","246750","246434",
    "246431","246424","246392","246391","246276","246275","246272","246271","246225","246222","246217","246185",
    "246182","246181","246180","246178","246153","246152","246151","246148","246115","246113","246089","246087",
    "246053","246052","246049","246048","246042","246037","246036","246012","246003","246001","245963","245954",
    "245952","245949","245923","245833","245831","245829","245793","245785","245775","245766","245759","245752",
    "245731","245729","245705","245702","245692","245683","245554","245545","245544","245543","245542","245540",
    "245535","245507","245505","245504","245501","245497","245496","245454","245453","245450","245446","245441",
    "245411","245410","245409","245407","245401","245397","245396","245353","245349","245347","245346","245345",
    "245343","245259","245233","245232","245231","245152","245151","245146","245145","245068","245066","245064",
    "244983","244982","244980","244975","244918","244917"};
  TString runNumList18q[125] = {
    "296623","296622","296621","296619","296618","296616","296615","296594","296553","296552","296551","296550",
    "296548","296547","296516","296512","296511","296510","296509","296472","296433","296424","296423","296420",
    "296419","296415","296414","296383","296381","296380","296379","296378","296377","296376","296375","296312",
    "296309","296304","296303","296280","296279","296273","296270","296269","296247","296246","296244","296243",
    "296242","296241","296240","296198","296197","296196","296195","296194","296192","296191","296143","296142",
    "296135","296134","296133","296132","296123","296074","296066","296065","296063","296062","296060","296016",
    "295942","295941","295937","295936","295913","295910","295909","295861","295860","295859","295856","295855",
    "295854","295853","295831","295829","295826","295825","295822","295819","295818","295816","295791","295788",
    "295786","295763","295762","295759","295758","295755","295754","295725","295723","295721","295719","295718",
    "295717","295714","295712","295676","295675","295673","295668","295667","295666","295615","295612","295611",
    "295610","295589","295588","295586","295585"};
  TString runNumList18r[89] = {
    "297595","297590","297588","297558","297544","297542","297541","297540","297537","297512","297483","297479",
    "297452","297451","297450","297446","297442","297441","297415","297414","297413","297406","297405","297380",
    "297379","297372","297367","297366","297363","297336","297335","297333","297332","297317","297311","297310",
    "297278","297222","297221","297218","297196","297195","297193","297133","297132","297129","297128","297124",
    "297123","297119","297118","297117","297085","297035","297031","296966","296941","296938","296935","296934",
    "296932","296931","296930","296903","296900","296899","296894","296852","296851","296850","296848","296839",
    "296838","296836","296835","296799","296794","296793","296790","296787","296786","296785","296784","296781",
    "296752","296694","296693","296691","296690"};
  runNumList = new std::map<int,int>;
  if      (fPeriod.EqualTo("LHC10h")) for (int i = 0; i < 90; i++) runNumList->insert(std::pair<int,int>(runNumList10h[i].Atoi(),i));
  else if (fPeriod.EqualTo("LHC15o")) for (int i = 0; i <138; i++) runNumList->insert(std::pair<int,int>(runNumList15o[i].Atoi(),i));
  else if (fPeriod.EqualTo("LHC18q")) for (int i = 0; i <125; i++) runNumList->insert(std::pair<int,int>(runNumList18q[i].Atoi(),i));
  else if (fPeriod.EqualTo("LHC18r")) for (int i = 0; i < 89; i++) runNumList->insert(std::pair<int,int>(runNumList18r[i].Atoi(),i));
  else return;

  fHistRunNumBin = new TH1I("runNumBin","",(int)runNumList->size(),0,(int)runNumList->size());
  std::map<int,int>::iterator iter;
  for (auto runNum : *runNumList) fHistRunNumBin->GetXaxis()->SetBinLabel(runNum.second, Form("%i",runNum.first));
  fOutputList->Add(fHistRunNumBin);

  ////////////////////////
  // Pile up Function
  ////////////////////////
  // Dobrin 15o pass2 Pile-up function
  if (fPeriod.EqualTo("LHC15o")) {
    fSPDCutPU = new TF1("fSPDCutPU", "450. + 3.9*x", 0, 50000);

    Double_t parV0[8] = {33.4237, 0.953516, 0.0712137, 227.923, 8.9239, -0.00319679, 0.000306314, -7.6627e-07};
    fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
    fV0CutPU->SetParameters(parV0);

    Double_t parV0CL0[6] = {0.0193587, 0.975914, 0.675714, 0.0292263, -0.000549509, 5.86421e-06};
    fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutLowPU->SetParameters(parV0CL0);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU->SetParameters(parV0CL0);

    Double_t parFB32[9] = {-812.822, 6.41796, 5421.83, -0.382601, 0.0299686, -26.6249, 321.388, -0.82615, 0.0167828};
    fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 6.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCutPU->SetParameters(parFB32);
  }

  // Rihan 18q/r Pile-up function
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

    Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
    fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
    fV0CutPU->SetParameters(parV0);

    Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
    fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutLowPU->SetParameters(parV0CL0);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU->SetParameters(parV0CL0);

    Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
    fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
    fMultCutPU->SetParameters(parFB32);
  }


  //------------------
  // QA
  //------------------
  // event-wise
  fEvtCount = new TH1D("EvtCount", "Event Count", 25, 1, 26);
  fEvtCount->GetXaxis()->SetBinLabel(1,"All");
  fEvtCount->GetXaxis()->SetBinLabel(2,"Read in");
  fEvtCount->GetXaxis()->SetBinLabel(3,"Event");
  fEvtCount->GetXaxis()->SetBinLabel(4,"Run Number");
  fEvtCount->GetXaxis()->SetBinLabel(5,"Vertex");
  fEvtCount->GetXaxis()->SetBinLabel(6,"Centrality");
  fEvtCount->GetXaxis()->SetBinLabel(7,"Pile up");
  fEvtCount->GetXaxis()->SetBinLabel(8,"Get VZERO Plane");
  fEvtCount->GetXaxis()->SetBinLabel(9,"Get ZDC Plane");
  fEvtCount->GetXaxis()->SetBinLabel(10,"Reset Vector");
  fEvtCount->GetXaxis()->SetBinLabel(11,"Loop Track");
  fEvtCount->GetXaxis()->SetBinLabel(12,"Get TPC Plane");
  fEvtCount->GetXaxis()->SetBinLabel(13,"Resolution");
  fEvtCount->GetXaxis()->SetBinLabel(14,"Loop V0");
  fEvtCount->GetXaxis()->SetBinLabel(15,"Pair Lambda");
  fEvtCount->GetXaxis()->SetBinLabel(20,"Manager");
  fEvtCount->GetXaxis()->SetBinLabel(21,"Handler");
  fEvtCount->GetXaxis()->SetBinLabel(22,"fAOD");
  fEvtCount->GetXaxis()->SetBinLabel(23,"fPID");
  fEvtCount->GetXaxis()->SetBinLabel(24,"fUtils");
  fEvtCount->GetXaxis()->SetBinLabel(25,"fMultSel");
  fOutputList->Add(fEvtCount);

  ////////////
  //QA Plots//
  ////////////
  // Event-wise QA
  fHistCent[0] = new TH1D("fHistCentBfCut", " Dist. of Centrality Before Cut ", 100, 0., 100.);
  fHistCent[1] = new TH1D("fHistCentAfCut", " Dist. of Centrality After Cut ", 100, 0., 100.);
  fOutputList->Add(fHistCent[0]);
  fOutputList->Add(fHistCent[1]);

  fHistVz[0] = new TH1D("fHistVzBfCut", "Dist of Centrality Before Cut", 200, -50., 50.);
  fHistVz[1] = new TH1D("fHistVzAfCut", "Dist of Centrality After Cut", 200, -50., 50.);
  fOutputList->Add(fHistVz[0]);
  fOutputList->Add(fHistVz[1]);

  fHist2DCentQA[0] = new TH2D("fHist2DCentQA_V0M_SPD1_BfCut", ";centV0M;centSPD1", 100, 0, 100, 100, 0, 100);
  fHist2DCentQA[1] = new TH2D("fHist2DCentQA_V0M_SPD1_AfCut", ";centV0M;centSPD1", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(fHist2DCentQA[0]);
  fOutputList->Add(fHist2DCentQA[1]);

  fHist2DCentQA[2] = new TH2D("fHist2DCentQA_V0M_TRK_BfCut", ";centV0M;centTRK", 100, 0, 100, 100, 0, 100);
  fHist2DCentQA[3] = new TH2D("fHist2DCentQA_V0M_TRK_AfCut", ";centV0M;centTRK", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(fHist2DCentQA[2]);
  fOutputList->Add(fHist2DCentQA[3]);

  fHist2DCentQA[4] = new TH2D("fHist2DCentQA_V0M_SPD0_BfCut", ";centV0M;centSPD0", 100, 0, 100, 100, 0, 100);
  fHist2DCentQA[5] = new TH2D("fHist2DCentQA_V0M_SPD0_AfCut", ";centV0M;centSPD0", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(fHist2DCentQA[4]);
  fOutputList->Add(fHist2DCentQA[5]);

  fHist2DCentQA[6] = new TH2D("fHist2DCentQA_SPD1_SPD0_BfCut", ";centSPD1;centSPD0", 100, 0, 100, 100, 0, 100);
  fHist2DCentQA[7] = new TH2D("fHist2DCentQA_SPD1_SPD0_AfCut", ";centSPD1;centSPD0", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(fHist2DCentQA[6]);
  fOutputList->Add(fHist2DCentQA[7]);

  fHist2DMultCentQA[0] = new TH2D("fHist2DMultCentQA_BfCut", ";centV0M;multFB32", 100, 0, 100, 20, 0, 5000);
  fHist2DMultCentQA[1] = new TH2D("fHist2DMultCentQA_AfCut", ";centV0M;multFB32", 100, 0, 100, 20, 0, 5000);
  fOutputList->Add(fHist2DMultCentQA[0]);
  fOutputList->Add(fHist2DMultCentQA[1]);

  // if (fMultComp.EqualTo("pileupByEDSTPC128") || fMultComp.EqualTo("pileupByGlobalTPC1")) {
  //   hMultMultQA[0] = new TH2D("hMultMultQAmTPCmESDPileupBefCut", "befCut;multTPC;multESD", 50, 0, 5000, 160, 0, 16000);
  //   hMultMultQA[1] = new TH2D("hMultMultQAmClobalmTPCEFPileupBefCut", "befCut;multGlobal;multTPCFE", 50, 0, 5000, 50, 0, 5000);
  //   hMultMultQA[2] = new TH2D("hMultMultQAmFB32mTrkTOFBefCut", "befCut;multTrkTOF;nTrk", 201, 0, 20000, 201, 0, 20000);
  //   hMultMultQA[3] = new TH2D("hMultMultQAmTPCmESDPileupAftCut", "aftCut;multTPC;multESD", 50, 0, 5000, 160, 0, 16000);
  //   hMultMultQA[4] = new TH2D("hMultMultQAmClobalmTPCEFPileupAftCut", "aftCut;multGlobal;multTPCFE", 50, 0, 5000, 50, 0, 5000);
  //   hMultMultQA[5] = new TH2D("hMultMultQAmFB32mTrkTOFAftCut", "aftCut;multTrkTOF;nTrk", 201, 0, 20000, 201, 0, 20000);
  //   fOutputList->Add(hMultMultQA[0]);
  //   fOutputList->Add(hMultMultQA[1]);
  //   fOutputList->Add(hMultMultQA[2]);
  //   fOutputList->Add(hMultMultQA[3]);
  //   fOutputList->Add(hMultMultQA[4]);
  //   fOutputList->Add(hMultMultQA[5]);
  // }

  // track-wise QA
  fHistPt  = new TH1D("fHistPt", "", 200, 0., 20.);
  fHistEta = new TH1D("fHistEta", "", 200, -10., 10.);
  fHistNhits = new TH1D("fHistNhits", "", 200, 0., 200.);
  fHist2DPDedx = new TH2D("fHist2DPDedx", "", 400, -10., 10., 400, 0, 1000);
  fHistDcaXY = new TH1D("fHistDcaXY", "", 100, 0., 10.);
  fHistDcaZ  = new TH1D("fHistDcaZ", "", 100, 0., 10.);
  fHistPhi[0] = new TH1D("fHistPhi", "", 100, 0, TMath::TwoPi());
  fHistPhi[1] = new TH1D("fHistPhi_afterNUA", "", 100, 0, TMath::TwoPi());
  fHist2DEtaPhi[0] = new TH2D("fHistEtaPhi", "", 16,-0.8,0.8, 100, 0, TMath::TwoPi());
  fHist2DEtaPhi[1] = new TH2D("fHistEtaPhi_afterfNUA", "", 16,-0.8,0.8, 100, 0, TMath::TwoPi());
  fOutputList->Add(fHistPt);
  fOutputList->Add(fHistEta);
  fOutputList->Add(fHistNhits);
  fOutputList->Add(fHist2DPDedx);
  fOutputList->Add(fHistDcaXY);
  fOutputList->Add(fHistDcaZ);
  for (int i = 0; i < 2; i++) fOutputList->Add(fHistPhi[i]);
  for (int i = 0; i < 2; i++) fOutputList->Add(fHist2DEtaPhi[i]);


  //V0s QA
  fHistV0Pt = new TH1D("hV0Pt","", 200, 0., 20.);
  fHistV0Eta = new TH1D("hV0Eta","", 200, -10., 10.);
  fHistV0DcatoPrimVertex = new TH1D("hV0DcaToPrimVertex","",200, 0., 20.);
  fHistV0CPA = new TH1D("hV0CPA","", 1000, 0.9, 1.);
  fHistV0DecayLength = new TH1D("hV0DecayLength","",500,0,500.);
  fOutputList->Add(fHistV0Pt);
  fOutputList->Add(fHistV0Eta);
  fOutputList->Add(fHistV0DcatoPrimVertex);
  fOutputList->Add(fHistV0CPA);
  fOutputList->Add(fHistV0DecayLength);

  //Lambda QA
  char chCut[5];
  for (int i = 0; i < 2; i++) {
    ///// Case 0 = before cut, case 1 = afterCut.
    if (i==0) sprintf(chCut,"Bf");
    if (i==1) sprintf(chCut,"Af");
    /// Lambdas:
    fHistLambdaPt[i] = new TH1D(Form("hLambdaPt_%sMassCut",chCut),"", 200, 0., 20.);
    fHistLambdaEta[i] = new TH1D(Form("hLambdaEta_%sMassCut",chCut),"",200, -10., 10.);
    fHistLambdaPhi[i] = new TH1D(Form("hLambdaPhi_%sMassCut",chCut),"", 360, 0., TMath::TwoPi());
    fHistLambdaDcaToPrimVertex[i] = new TH1D(Form("hLambdaDcaToPrimVertex_%sMassCut",chCut),"",200, 0., 20.);
    fHistLambdaCPA[i] = new TH1D(Form("hLambdaCPA_%sMassCut",chCut),"",200, 0.9, 1.);
    fHistLambdaDecayLength[i] = new TH1D(Form("hLambdaDecayLength_%sMassCut",chCut),"", 250, 0., 500.);
    fHistLambdaMass[i] = new TH1D(Form("hLambdaMass_%sMassCut",chCut),"",1000,1.,1.25); //  Current bin size = 0.00025
    fProfileLambdaMassVsPt[i] = new TProfile(Form("pLambdaMassVsPt_%sMassCut",chCut),"",200,0,20);
    fOutputList->Add(fHistLambdaPt[i]);
    fOutputList->Add(fHistLambdaEta[i]);
    fOutputList->Add(fHistLambdaPhi[i]);
    fOutputList->Add(fHistLambdaDcaToPrimVertex[i]);
    fOutputList->Add(fHistLambdaCPA[i]);
    fOutputList->Add(fHistLambdaDecayLength[i]);
    fOutputList->Add(fHistLambdaMass[i]);
    fOutputList->Add(fProfileLambdaMassVsPt[i]);

    // AntiLambdas
    fHistAntiLambdaPt[i] = new TH1D(Form("hAntiLambdaPt_%sMassCut",chCut),"", 200, 0., 20.);
    fHistAntiLambdaEta[i] = new TH1D(Form("hAntiLambdaEta_%sMassCut",chCut),"",200, -10., 10.);
    fHistAntiLambdaPhi[i] = new TH1D(Form("hAntiLambdaPhi_%sMassCut",chCut),"", 360, 0, TMath::TwoPi());
    fHistAntiLambdaDcaToPrimVertex[i] = new TH1D(Form("hAntiLambdaDcaToPrimVertex_%sMassCut",chCut),"",200, 0., 20.);
    fHistAntiLambdaCPA[i] = new TH1D(Form("hAntiLambdaCPA_%sMassCut",chCut),"",200, 0.9, 1.);
    fHistAntiLambdaDecayLength[i] = new TH1D(Form("hAntiLambdaDecayLength_%sMassCut",chCut),"", 250, 0., 500.);
    fHistAntiLambdaMass[i] = new TH1D(Form("hAntiLambdaMass_%sMassCut",chCut),"",1000,1.,1.25); // Current bin size = 0.00025
    fProfileAntiLambdaMassVsPt[i] = new TProfile(Form("pAntiLambdaMassVsPt_%sMassCut",chCut),"",200,0,20);
    fOutputList->Add(fHistAntiLambdaPt[i]);
    fOutputList->Add(fHistAntiLambdaEta[i]);
    fOutputList->Add(fHistAntiLambdaPhi[i]);
    fOutputList->Add(fHistAntiLambdaDcaToPrimVertex[i]);
    fOutputList->Add(fHistAntiLambdaCPA[i]);
    fOutputList->Add(fHistAntiLambdaDecayLength[i]);
    fOutputList->Add(fHistAntiLambdaMass[i]);
    fOutputList->Add(fProfileAntiLambdaMassVsPt[i]);
  }



  ////////////////////////
  // Results
  ////////////////////////

  ///Lambda-X correlators

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      fHistSig[i][j] = new TH1F(Form("hSig_%d_%d",i,j),"",30,-0.5*TMath::Pi(),1.5*TMath::Pi());
      fHistBkg[i][j] = new TH1F(Form("hBkg_%d_%d",i,j),"",30,-0.5*TMath::Pi(),1.5*TMath::Pi());
      fOutputList->Add(fHistSig[i][j]);
      fOutputList->Add(fHistBkg[i][j]);
    }
  }
  


  PostData(1,fOutputList);
  if (fDebug) Printf("UserCreateOutputObjects() Post Data Success!");
}

//------------------------------------------------

void AliAnalysisTaskPIDAngCorr::UserExec(Option_t *)
{
  if (fDebug) Printf("===============================We are in UserExec!================================");
  fEvtCount->Fill(1);
  //----------------------------
  // Handle
  //----------------------------
  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  } else fEvtCount->Fill(20);

  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  } else fEvtCount->Fill(21);

  fAOD = dynamic_cast <AliAODEvent*> (InputEvent());
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
  } else fEvtCount->Fill(22);

  fPIDResponse = handler->GetPIDResponse();
  if (!fPIDResponse) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
  } else fEvtCount->Fill(23);

  fUtils = new AliAnalysisUtils();
  if (!fUtils) {
    AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  } else fEvtCount->Fill(24);

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    if (!fMultSel) {
      AliError(Form("%s: Could not get AliMultSelection", GetName()));
    } else fEvtCount->Fill(25);
    if (!manager || !handler || !fAOD || !fPIDResponse || !fUtils || !fMultSel) return;
  }

  if (!manager || !handler || !fAOD || !fPIDResponse || !fUtils) return;
  fEvtCount->Fill(2);
  if (fDebug) Printf("Handles done!");

  //----------------------------
  // Trigger
  //----------------------------
  UInt_t mask = handler->IsEventSelected();
  bool isTrigselected = false;
  if (fTrigger.EqualTo("kMB"))
  isTrigselected = mask & AliVEvent::kMB;
  else if (fTrigger.EqualTo("kINT7"))
  isTrigselected = mask & AliVEvent::kINT7;
  else if (fTrigger.EqualTo("kINT7+kCentral+kSemiCentral"))
  isTrigselected = mask & (AliVEvent::kINT7 + AliVEvent::kCentral + AliVEvent::kSemiCentral);
  if (isTrigselected == false) return;
  fEvtCount->Fill(3);
  if (fDebug) Printf("trigger done!");

  //----------------------------
  // Run Number
  //----------------------------
  fRunNum = fAOD->GetRunNumber();
  // if (fRunNum != fOldRunNum) {
  //    // Load the run dependent calibration hist
  //     if (!LoadCalibHistForThisRun()) return;
  //     fRunNumBin = runNumList->at(fRunNum);
  //     fOldRunNum = fRunNum;
  //     if (fRunNumBin < 0) return;
  // }
  fHistRunNumBin->Fill(fRunNumBin);
  fEvtCount->Fill(4);
  if (fDebug) Printf("run nummbr done!");

  //----------------------------
  // Vertex
  //----------------------------
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  fVtx -> GetXYZ(fVertex);
  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  double vx = fVertex[0];
  double vy = fVertex[1];
  double vz = fVertex[2];
  if (fabs(fVertex[0])<1e-6 || fabs(fVertex[1])<1e-6 || fabs(fVertex[2])<1e-6) return;
  double dz = vz - fAOD->GetPrimaryVertexSPD()->GetZ();
  if (fabs(vz) > fVzCut) return;
  if (!fVtx || fVtx->GetNContributors() < 2 || vtSPD->GetNContributors()<1) return;
  fHistVz[0]->Fill(vz);
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp
  // fEventCuts.SetCentralityEstimators("V0M","CL1");
  // if (!fEventCuts->AcceptEvent(fAOD) ) return;
  if (fPeriod.EqualTo("LHC10h")) if (fabs(dz)>0.5) return;
  if (fPeriod.EqualTo("LHC15o")) {
      double covTrc[6],covSPD[6];
      fVtx->GetCovarianceMatrix(covTrc);
      fAOD->GetPrimaryVertexSPD()->GetCovarianceMatrix(covSPD);
      double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
      double errTrc = TMath::Sqrt(covTrc[5]);
      double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
      if (fabs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return;
  }
  fHistVz[1]->Fill(vz);
  for (int i = 0; i < 20; ++i) {
      if (vz > -10+i*1 && vz < -10+(i+1)*1) {fVzBin = i; break;}
  }
  if (fVzBin<-990) return;
  fEvtCount->Fill(5);
  if (fDebug) Printf("vertex done!");

  //----------------------------
  // Centrality
  //----------------------------
  double centV0M = -1, centTRK = -1, centSPD0 = -1, centSPD1 = -1, centV0A = -1;
  if (fPeriod.EqualTo("LHC10h")) {
    centV0M  = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
    centTRK  = fAOD->GetCentrality()->GetCentralityPercentile("TRK");
    centSPD0 = fAOD->GetCentrality()->GetCentralityPercentile("CL0");
    centSPD1 = fAOD->GetCentrality()->GetCentralityPercentile("CL1");
    centV0A  = fAOD->GetCentrality()->GetCentralityPercentile("V0A");
  } else if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    centV0M  = fMultSel->GetMultiplicityPercentile("V0M");
    centTRK  = fMultSel->GetMultiplicityPercentile("TRK");
    centSPD0 = fMultSel->GetMultiplicityPercentile("CL0");
    centSPD1 = fMultSel->GetMultiplicityPercentile("CL1");
  }
  //we use centV0M as the default centrality
  fCent = centV0M;
  fHist2DCentQA[0]->Fill(centV0M,centSPD1);
  fHist2DCentQA[2]->Fill(centV0M,centTRK);
  fHist2DCentQA[4]->Fill(centV0M,centSPD0);
  fHist2DCentQA[6]->Fill(centSPD1,centSPD0);
  if (fabs(fCent-centSPD1)>fCentDiffCut) return;
  fHist2DCentQA[1]->Fill(centV0M,centSPD1);
  fHist2DCentQA[3]->Fill(centV0M,centTRK);
  fHist2DCentQA[5]->Fill(centV0M,centSPD0);
  fHist2DCentQA[7]->Fill(centSPD1,centSPD0);
  if (fCent < 40 || fCent >= 50) return;
  // cent bin
  fCentBin = (int)fCent/10;
  fHistCent[0]->Fill(fCent);
  fEvtCount->Fill(6);
  if (fDebug) Printf("centrality done!");

  //----------------------------
  // Pile up
  //----------------------------
  if (fPeriod.EqualTo("LHC10h")) if (!RemovalForRun1()) return;
  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    // hMultCentQA[0]->Fill(fCent, fAOD->GetNumberOfTracks()); // raw Trk Multi Vs Cent(V0M)
    // if (PileUpMultiVertex(fAOD)) return;
    // if (!RejectEvtMultComp(fAOD)) return;
    // hMultCentQA[1]->Fill(fCent, fAOD->GetNumberOfTracks()); // Mult_Cent QA
    // if (!AODPileupCheck (fAOD)) return;
    if (!RejectEvtTFFit()) return; // 15o_pass2
  }
  fHistCent[1]->Fill(fCent);
  fEvtCount->Fill(7);
  if (fDebug) Printf("pile-up done!");

  //----------------------------
  // Loop Tracks / Fill Vectors
  //----------------------------
  //Reset vectors
  ResetVectors();
  fEvtCount->Fill(10);
  if (!LoopTracks()) return;
  fEvtCount->Fill(11);
  if (fDebug) Printf("Loop Tracks done!");
  //----------------------------
  // Get Lambda Vector
  //----------------------------
  if (!LoopV0s()) return;
  fEvtCount->Fill(14);
  if (fDebug) Printf("Get Lambda Vector done!");
  //----------------------------
  // Pair
  //----------------------------
  if (!Pair()) return;
  fEvtCount->Fill(15);
  if (fDebug) Printf("Pair done!");
  //----------------------------
  // Mix Event
  //----------------------------
  for (int j = 0; j < 6; j++) {
      if (vecNumMixBuffer[j].size() == 10) {
      if (fDebug) cout<<"No."<<j<<" particle is ready to mix pair"<<endl;
      for (int i = 0; i < 6; i++) {
        if (vecPhi[i].size() == 0) continue;
        for (float phi_1 : vecPhi[i]) {
          for (float phi_2 : vecPhiMixBuffer[j]) {
            if(phi_1 == phi_2) cout<<"tell me why?"<<endl;
            fHistBkg[i][j] -> Fill(RangePhi(phi_1 - phi_2));
          }
        }
      }
    }
  }

  for (int i = 0; i < 6; i++) {
    if (vecPhi[i].size() != 0) {
      vecPhiMixBuffer[i].insert(vecPhiMixBuffer[i].end(), vecPhi[i].begin(), vecPhi[i].end());
      vecNumMixBuffer[i].push_back(vecPhi[i].size());
    }
  }
  for (int i = 0; i < 6; i++) {
    if (vecNumMixBuffer[i].size() == 11) {
      vecPhiMixBuffer[i].erase(vecPhiMixBuffer[i].begin(), vecPhiMixBuffer[i].begin()+vecNumMixBuffer[i][0]);
      vecNumMixBuffer[i].erase(vecNumMixBuffer[i].begin());
    }
    if (i == 5) 
    {for (double phi : vecPhiMixBuffer[i]) cout<<phi<<endl; for (int n : vecNumMixBuffer[i]) cout<<n<<endl;}
  }

  fEvtCount->Fill(15);
  if (fDebug) Printf("MixPair done!");
  //------------------
  // Post output data.
  //------------------
  if (fDebug) Printf("analysis done!");
  PostData(1,fOutputList);
}


//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::LoopTracks()
{
  int nTrks = fAOD->GetNumberOfTracks();
  if (nTrks < 4) return false;
  for (int iTrk = 0; iTrk < nTrks; ++iTrk) {
    AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTrk);
    if (!track) {
      AliError(Form("%s: Could not get Track", GetName()));
      continue;
    }
    if (!track->TestFilterBit(fFilterBit)) continue;
    if (!AcceptAODTrack(track)) continue;
    //------------------
    // NUE & NUA
    //------------------
    double phi = track->Phi();
    double  pt = track->Pt();
    double eta = track->Eta();
    int charge = track->Charge();
    int     id = track->GetID();

    fHistPhi[0]->Fill(phi);
    fHist2DEtaPhi[0]->Fill(eta,phi);

    bool isItPiontrk = CheckPIDofParticle(track,1); // 1=pion
    bool isItProttrk = CheckPIDofParticle(track,3); // 3=proton
    isItProttrk *= (pt < fProtonPtMax);

    if(isItPiontrk) isItProttrk = false;

    if (charge > 0 && isItPiontrk) {
      vecID[0].push_back(id);
      vecPt[0].push_back(pt);
      vecPhi[0].push_back(phi);
    }
    if (charge < 0 && isItPiontrk) {
      vecID[1].push_back(id);
      vecPt[1].push_back(pt);
      vecPhi[1].push_back(phi);
    }
    if (charge > 0 && isItProttrk) {
      vecID[2].push_back(id);
      vecPt[2].push_back(pt);
      vecPhi[2].push_back(phi);
    }
    if (charge < 0 && isItProttrk) {
      vecID[3].push_back(id);
      vecPt[3].push_back(pt);
      vecPhi[3].push_back(phi);
    }
  }
  return true;
}


//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::LoopV0s()
{
  int nV0s = fAOD->GetNumberOfV0s();
  for (int iV0 = 0; iV0 < nV0s; iV0++) {
    AliAODv0 *v0 = fAOD->GetV0(iV0);
    if (!v0) {
      AliError(Form("%s: Could not get v0s", GetName()));
      continue;
    }
    //Basic kinematic variable
    double pt      = v0->Pt();
    double eta     = v0->PseudoRapV0();
    double dcaToPV = v0->DcaV0ToPrimVertex();//DCA to Primary Vertex
    double CPA     = v0->CosPointingAngle(fVertex);//cosine pointing angle
    double dl      = v0->DecayLengthV0(fVertex);
    fHistV0Pt              -> Fill(pt);
    fHistV0Eta             -> Fill(eta);
    fHistV0DcatoPrimVertex -> Fill(dcaToPV);
    fHistV0CPA             -> Fill(CPA);
    fHistV0DecayLength     -> Fill(dl);
    //V0 cut
    if (!IsGoodV0(v0)) continue;
    //V0 daughters cut
    AliAODTrack *nTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));
    AliAODTrack *pTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
    if (!(IsGoodDaughterTrack(nTrack)) || !(IsGoodDaughterTrack(pTrack))) continue;
    float nDcaPV = v0->DcaNegToPrimVertex();
    float pDcaPV = v0->DcaPosToPrimVertex();
    if ( nDcaPV<fDaughtersDCAToPrimVtxMin || pDcaPV<fDaughtersDCAToPrimVtxMin) continue;
    int code = GetLambdaCode(pTrack,nTrack);
    if (TMath::Abs(code) != 3122) continue;
    TVector2 Vt(v0->MomV0X(), v0->MomV0Y());
    double phi = Vt.Phi() > 0 ? Vt.Phi() : Vt.Phi() + TMath::TwoPi();
    int id_posDaughter = v0->GetPosID();
    int id_negDaughter = v0->GetNegID();

    if (code == 3122) {
      double massLambda  = v0->MassLambda();
      fHistLambdaPt[0]              -> Fill(pt);
      fHistLambdaEta[0]             -> Fill(eta);
      fHistLambdaPhi[0]             -> Fill(phi);
      fHistLambdaDcaToPrimVertex[0] -> Fill(dcaToPV);
      fHistLambdaCPA[0]             -> Fill(CPA);
      fHistLambdaDecayLength[0]     -> Fill(dl);
      fHistLambdaMass[0]            -> Fill(massLambda);
      fProfileLambdaMassVsPt[0]     -> Fill(pt, massLambda);

      if (TMath::Abs(massLambda - fMassMean) < fLambdaMassCut) {
        // //if a particle has been used as daughter particle before(It happends), we have to refuse a new one.
        // if (find(vecDaughterPosID[4].begin(), vecDaughterPosID[0].end(), id_posDaughter) != vecDaughterPosID[0].end()) continue;
        // if (find(vecDaughterNegID[4].begin(), vecDaughterNegID[0].end(), id_negDaughter) != vecDaughterNegID[0].end()) continue;
        fHistLambdaPt[1]              -> Fill(pt);
        fHistLambdaEta[1]             -> Fill(eta);
        fHistLambdaPhi[1]             -> Fill(phi);
        fHistLambdaDcaToPrimVertex[1] -> Fill(dcaToPV);
        fHistLambdaCPA[1]             -> Fill(CPA);
        fHistLambdaDecayLength[1]     -> Fill(dl);
        fHistLambdaMass[1]            -> Fill(massLambda);
        fProfileLambdaMassVsPt[1]     -> Fill(pt, massLambda);

        vecDaughterPosID[4].push_back(id_posDaughter); // Pos Daughter ID
        vecDaughterNegID[4].push_back(id_negDaughter); // Neg Daughter ID
        vecPt[4].push_back(pt);
        vecPhi[4].push_back(phi);
      }
    }

    if (code == -3122) {
      double massAntiLambda  = v0->MassAntiLambda();
      fHistAntiLambdaPt[0]              -> Fill(pt);
      fHistAntiLambdaEta[0]             -> Fill(eta);
      fHistAntiLambdaPhi[0]             -> Fill(phi);
      fHistAntiLambdaDcaToPrimVertex[0] -> Fill(dcaToPV);
      fHistAntiLambdaCPA[0]             -> Fill(CPA);
      fHistAntiLambdaDecayLength[0]     -> Fill(dl);
      fHistAntiLambdaMass[0]            -> Fill(massAntiLambda);
      fProfileAntiLambdaMassVsPt[0]     -> Fill(pt, massAntiLambda);

      if (TMath::Abs(massAntiLambda - fMassMean) < fLambdaMassCut) {
        // if (find(vecDaughterPosID[1].begin(), vecDaughterPosID[1].end(), id_posDaughter) != vecDaughterPosID[1].end()) continue;
        // if (find(vecDaughterNegID[1].begin(), vecDaughterNegID[1].end(), id_negDaughter) != vecDaughterNegID[1].end()) continue;
        fHistAntiLambdaPt[1]              -> Fill(pt);
        fHistAntiLambdaEta[1]             -> Fill(eta);
        fHistAntiLambdaPhi[1]             -> Fill(phi);
        fHistAntiLambdaDcaToPrimVertex[1] -> Fill(dcaToPV);
        fHistAntiLambdaCPA[1]             -> Fill(CPA);
        fHistAntiLambdaDecayLength[1]     -> Fill(dl);
        fHistAntiLambdaMass[1]            -> Fill(massAntiLambda);
        fProfileAntiLambdaMassVsPt[1]     -> Fill(pt, massAntiLambda);

        vecDaughterPosID[5].push_back(id_posDaughter); // Pos Daughter ID
        vecDaughterNegID[5].push_back(id_negDaughter); // Neg Daughter ID
        vecPt[5].push_back(pt);
        vecPhi[5].push_back(phi);
      }
    }
  }//loop V0 end
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::Pair()
{
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      for (int iPhi = 0; iPhi < (int)vecPhi[i].size(); iPhi++) {
        for (int jPhi = 0; jPhi < (int)vecPhi[j].size(); jPhi++) {
          if (i == j && iPhi == jPhi) continue;
          if (i >= 4 && j < 4) {
            int id_posDaughter = vecDaughterPosID[i][iPhi];
            int id_negDaughter = vecDaughterNegID[i][iPhi];
            int id = vecID[j][jPhi];
            if (id == id_posDaughter || id == id_negDaughter) continue;
          }
          if (i < 4 && j >= 4) {
            int id_posDaughter = vecDaughterPosID[j][jPhi];
            int id_negDaughter = vecDaughterNegID[j][jPhi];
            int id = vecID[i][iPhi];
            if (id == id_posDaughter || id == id_negDaughter) continue;
          }

          if (i >= 4 && j >= 4) {
            int id_posDaughter_0 = vecDaughterPosID[i][iPhi];
            int id_negDaughter_0 = vecDaughterNegID[i][iPhi];
            int id_posDaughter_1 = vecDaughterPosID[j][jPhi];
            int id_negDaughter_1 = vecDaughterNegID[j][jPhi];
            if (id_posDaughter_0 == id_posDaughter_1 || id_negDaughter_0 == id_negDaughter_1 ||
                id_posDaughter_0 == id_negDaughter_1 || id_negDaughter_0 == id_posDaughter_1
            ) continue;
          }
          double phi1 = vecPhi[i][iPhi];
          double phi2 = vecPhi[j][jPhi];
          fHistSig[i][j] -> Fill(RangePhi(phi1 - phi2));
        }
      }
    }
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::MixPair()
{
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      for (int iPhi = 0; iPhi < (int)vecPhi[i].size(); iPhi++) {
        for (int jPhi = 0; jPhi < (int)vecPhi[j].size(); jPhi++) {
          double phi1 = vecPhi[i][iPhi];
          double phi2 = vecPhiMixBuffer[j][jPhi];
          fHistBkg[i][j] -> Fill(RangePhi(phi1 - phi2));
        }
      }
    }
  }
  return true;
}

//---------------------------------------------------

void AliAnalysisTaskPIDAngCorr::ResetVectors()
{
  for (int i = 0; i < 6; i++) {
    std::vector<int>().swap(vecID[i]);
    std::vector<float>().swap(vecPt[i]);
    std::vector<float>().swap(vecPhi[i]);
    std::vector<int>().swap(vecDaughterPosID[i]);
    std::vector<int>().swap(vecDaughterNegID[i]);
  }
}


bool AliAnalysisTaskPIDAngCorr::RemovalForRun1()
{
  // pileup
  fUtils->SetUseOutOfBunchPileUp(true);
  fUtils->SetUseMVPlpSelection(true);
  // fUtils->SetMinPlpContribMV(5);
  bool isPileup = fUtils->IsPileUpEvent(fAOD);
  // bool isPileup = fUtils->IsPileUpMV(fAOD); // pp, p-Pb
  if (isPileup) return false;
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::RejectEvtMultComp() // 15o_pass1, old pile-up
{
   // TPC cluster cut
    Int_t multEsd = ((AliAODHeader*)fAOD->GetHeader())->GetNumberOfESDTracks(); // multESD
    const Int_t nTrk = fAOD->GetNumberOfTracks();
    Int_t multTPC=0; // FB128 + Common Track Cuts
    Int_t multTPCFE=0; // FB1 + Common Track Cuts + chi2 > 0.2
    Int_t multGlobal=0; // FB16 + Common Track Cuts + Dca Cuts
    for (Int_t it1 = 0; it1 < nTrk; it1++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it1);
      if (!aodTrk) continue;
      if (aodTrk->TestFilterBit(128)) multTPC++;
      double Eta  = aodTrk->Eta();
      double Pt    = aodTrk->Pt();
      double Phi  = aodTrk->Phi();
      if (Pt<0.2 || Pt>5.0 || TMath::Abs(Eta)>0.8 || aodTrk->GetTPCNcls()<fNclsCut || aodTrk->GetTPCsignal()<10.0) continue;
      if (aodTrk->TestFilterBit(1) && aodTrk->Chi2perNDF()>0.2)  multTPCFE++;
      if (!aodTrk->TestFilterBit(16) || aodTrk->Chi2perNDF()<0.1)   continue;
      Double_t dca[2] = {-99., -99.};
      Double_t cov[3] = {-99., -99., -99.};
      Double_t magField = fAOD->GetMagneticField();
      if (magField!=0) {
        if (aodTrk->PropagateToDCA(fAOD->GetPrimaryVertex(), magField, 100., dca, cov) && TMath::Abs(dca[0]) < 0.3 && TMath::Abs(dca[1]) < 0.3) multGlobal++;
      }
    }

    fHist2DMultMultQA[0]->Fill(multTPC,multEsd);
    fHist2DMultMultQA[1]->Fill(multGlobal,multTPCFE);

    //TODO
    TString fMultComp = "pileupByGlobalTPC1";

    if (fMultComp.EqualTo("pileupByEDSTPC128")) { // Rihan
      if ((Double_t)(multEsd*(1/3.45) - 90.) < (Double_t)multTPC)
      {
        fHist2DMultMultQA[3]->Fill(multTPC,multEsd);
        return true;
      }
      else return false;
    }

    if (fMultComp.EqualTo("pileupByGlobalTPC1")) { // A.Dobrin
      if (multTPCFE-1.78*multGlobal<62.87 && multTPCFE-1.48*multGlobal>-36.73) {
        fHist2DMultMultQA[4]->Fill(multGlobal,multTPCFE);
        return true;
      }
      else return false;
    }

    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::RejectEvtTFFit()
{
  Float_t centV0M = -999;
  Float_t centCL1 = -999;
  Float_t centCL0 = -999;

  AliMultSelection* fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
  if (!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(1);
  }
  centV0M = (Float_t) fCent;
  centCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centCL0 = fMultSelection->GetMultiplicityPercentile("CL0");

  Int_t nITSClsLy0 = fAOD->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = fAOD->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();

  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multTrk = 0;
  for (Int_t it = 0; it < nTracks; it++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
      if (!aodTrk) {
          delete aodTrk;
          continue;
      }
      if (aodTrk->TestFilterBit(32)) {
        if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
      }
  }

  fHist2DMultCentQA[0]->Fill(centV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)

  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  Float_t  multV0a = aodV0->GetMTotV0A();
  Float_t  multV0c = aodV0->GetMTotV0C();
  Float_t  multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;
  
  // //
  // Int_t tpcClsTot = fAOD->GetNumberOfTPCClusters();
  // Float_t nclsDif = Float_t(tpcClsTot) - (53182.6 + 113.326*multV0Tot - 0.000831275*multV0Tot*multV0Tot);

  // pile-up cuts
  if (centCL0 < fCenCutLowPU->Eval(centV0M)) return false;
  if (centCL0 > fCenCutHighPU->Eval(centV0M)) return false;
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) return false;
  if (multV0On < fV0CutPU->Eval(multV0Tot)) return false;
  if (Float_t(multTrk) < fMultCutPU->Eval(centV0M)) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fAOD->IsIncompleteDAQ()) return false;

  fHist2DMultCentQA[1]->Fill(centV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::RejectEvtTPCITSfb32TOF ()
{
    //TOD+FB32 pile-up removal
    // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp
    Int_t multTrk=0;
    Int_t multTrkTOF=0;
    int nTrk = fAOD->GetNumberOfTracks();
    for (Int_t it2 = 0; it2 < nTrk; it2++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it2);
      if (!aodTrk) continue;
      if (aodTrk->TestFilterBit(32)) {
        multTrk++;
        if ( TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10 && aodTrk->GetTOFsignal() >= 12000 && aodTrk->GetTOFsignal() <= 25000) multTrkTOF++;
        else return false;
      }
    }
    fHist2DMultMultQA[2]->Fill(multTrkTOF, nTrk);
    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::AODPileupCheck()
{
  Int_t isPileup = fAOD->IsPileupFromSPD(3);
  if (isPileup !=0 && fPeriod.EqualTo("LHC16t")) return false; // LHC16t : pPb
  if (fAOD->IsIncompleteDAQ()) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fPeriod.EqualTo("LHC15o")) {
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    if (!fMultSel->GetThisEventIsNotPileup()) return false;
    if (!fMultSel->GetThisEventIsNotPileupMV()) return false;
    if (!fMultSel->GetThisEventIsNotPileupInMultBins()) return false;
    if (!fMultSel->GetThisEventHasNoInconsistentVertices()) return false;
    if (!fMultSel->GetThisEventPassesTrackletVsCluster()) return false;
    if (!fMultSel->GetThisEventIsNotIncompleteDAQ()) return false;
    if (!fMultSel->GetThisEventHasGoodVertex2016()) return false;
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::PileUpMultiVertex()
{
  // check for multi-vertexer pile-up
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2    = 5.0;
  const double kMinWDist      = 15;

  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;

  int nPlp = 0;

  if (!(nPlp=fAOD->GetNumberOfPileupVerticesTracks()))
  return false;

  vtPrm = fAOD->GetPrimaryVertex();
  if (vtPrm == fAOD->GetPrimaryVertexSPD())
  return true;  // there are pile-up vertices but no primary

  //int bcPrim = vtPrm->GetBC();

  for (int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)fAOD->GetPileupVertexTracks(ipl);
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF()    > kMaxPlpChi2)    continue;
    //int bcPlp = vtPlp->GetBC();
    //if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2)
    // return kTRUE; // pile-up from other BC

    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist)        continue;

    return true; // pile-up: well separated vertices
  }
  return false;
}

//---------------------------------------------------

double AliAnalysisTaskPIDAngCorr::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
    // calculate sqrt of weighted distance to other vertex
    if (!v0 || !v1) {
        printf("One of vertices is not valid\n");
        return 0;
    }
    static TMatrixDSym vVb(3);
    double dist = -1;
    double dx = v0->GetX()-v1->GetX();
    double dy = v0->GetY()-v1->GetY();
    double dz = v0->GetZ()-v1->GetZ();
    double cov0[6],cov1[6];
    v0->GetCovarianceMatrix(cov0);
    v1->GetCovarianceMatrix(cov1);
    vVb(0,0) = cov0[0]+cov1[0];
    vVb(1,1) = cov0[2]+cov1[2];
    vVb(2,2) = cov0[5]+cov1[5];
    vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
    vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
    vVb.InvertFast();
    if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
    dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
        +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
    return dist>0 ? TMath::Sqrt(dist) : -1;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::AcceptAODTrack(AliAODTrack *track)
{
    //------------------
    // track cut
    //------------------
    double pt     = track->Pt();
    double eta    = track->Eta();
    int    nhits  = track->GetTPCNcls();
    double dedx   = track->GetTPCsignal();
    double chi2   = track->Chi2perNDF();
    int    charge = track->Charge();

    if ( pt < fPtMin
      || pt > fPtMax
      || fabs(eta) > fEtaCut
      || fabs(nhits) < fNclsCut
      || chi2 < fChi2Min 
      || chi2 > fChi2Max
      || dedx < fDedxCut) 
    return false;

    if (fFilterBit != 768){
      if (fPeriod.EqualTo("LHC10h")) {
        //------------------
        // dca cut
        //------------------
        double mag = fAOD->GetMagneticField();
        double dcaxy  = 999.;
        double dcaz   = 999.;
        double r[3];
        double dca[2];
        double cov[3];
        double vx = fVertex[0];
        double vy = fVertex[1];
        double vz = fVertex[2];
        bool proptodca = track->PropagateToDCA(fAOD->GetPrimaryVertex(), mag, 100., dca, cov);
        if (track->GetXYZ(r)) {
          dcaxy = r[0];
          dcaz  = r[1];
        } else {
          double dcax = r[0] - vx;
          double dcay = r[1] - vy;
          dcaz  = r[2] - vz;
          dcaxy = sqrt(dcax*dcax + dcay*dcay);
        }
        if (fabs(dcaxy)>fDcaCutxy) return false;
        if (fabs(dcaz)>fDcaCutz) return false;
        fHistDcaXY->Fill(dcaxy);
        fHistDcaZ->Fill(dcaz);
      }
    }

    fHistPt->Fill(pt);
    fHistEta->Fill(eta);
    fHistNhits->Fill(nhits);
    fHist2DPDedx->Fill(track->P()*charge, dedx);

    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck)
{
  if (pidToCheck==0) return kTRUE;    //// Charge Particles do not need PID check
  bool bPIDokay = kFALSE;

  if (!fPIDResponse) {
    Printf("\n Could Not access PIDResponse Task, please Add the Task...\n return with kFALSE pid\n");
    return kFALSE;
  }

  /// Rihan todo: To set the low pT cuts for nSigmaTPC from AddTaskMacro!
  /// Although someone barely needs to change it given the purity..

  float nSigTPC = 0, nSigTOF = 0, nSigRMS = 0;
  bool trkPtPID  = ftrack->Pt();
  int  trkChargePID = ftrack->Charge();

  ///Pion =>
  if (pidToCheck==1) {
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kPion);//Some warning show here (***TDatabasePDG::AddParicle: particle with PDGcode = 3124 already defind),I don't understand what happended. --chunzheng
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kPion);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if (trkPtPID<=0.5 && TMath::Abs(nSigTPC)<=fNSigmaTPCCut)     bPIDokay = kTRUE;
    // Using TPCTOF RMS cut for higher pt:
    else if (trkPtPID>0.5 && TMath::Abs(nSigRMS)<=fNSigmaTOFCut) bPIDokay = kTRUE;
    return bPIDokay;
  }
  ///Kaon =>
  else if (pidToCheck==2) {
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kKaon);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kKaon);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if (trkPtPID<=0.45 && TMath::Abs(nSigTPC)<=fNSigmaTPCCut)     bPIDokay = kTRUE;
    else if (trkPtPID>0.45 && TMath::Abs(nSigRMS)<=fNSigmaTOFCut) bPIDokay = kTRUE;
    return bPIDokay;
  }
  ///proton =>
  else if (pidToCheck==3) {///
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kProton);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kProton);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if (trkPtPID<=0.6 && TMath::Abs(nSigTPC)<=fNSigmaTPCCut) {
      bPIDokay = kTRUE;
      if (trkChargePID>0 && trkPtPID<0.4) bPIDokay = kFALSE;
    }
    else if (trkPtPID>0.6 && TMath::Abs(nSigRMS)<=fNSigmaTOFCut) {
      bPIDokay = kTRUE;
    }
    return bPIDokay;
  }
  else{
    Printf("\n -Ve number not allowed! Choose among: 0,1,2,3 (Charge Pion, Kaon, Proton)\n return with kFALSE \n");
    return kFALSE;
  }

  return kFALSE;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::IsGoodV0(AliAODv0 *aodV0)
{
  // Offline reconstructed V0 only
  if (aodV0->GetOnFlyStatus()) return false;
  // Get daughters and check them
  AliAODTrack *myTrackNegTest = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(1));
  AliAODTrack *myTrackPosTest = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(0));
  if (!myTrackPosTest || !myTrackNegTest) {
    Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
    return false;
  }
  // Unlike signs of daughters
  if (myTrackNegTest->Charge() == myTrackPosTest->Charge()) return false;
  // Cosinus of pointing angle
  double dCPA = aodV0->CosPointingAngle(fVertex);
  // cut on Cosinus of pointing angle
  if (dCPA < fV0CPAMin) return false;
  // DCA of V0
  double dV0Dca = aodV0->DcaV0ToPrimVertex();
  if (TMath::Abs(dV0Dca) > fV0DCAToPrimVtxMax) return false;
  // V0 path length before decay
  double dDecayLength = aodV0->DecayLengthV0(fVertex);
  if (dDecayLength > fV0DecayLengthMax) return false;
  if (dDecayLength < fV0DecayLengthMin) return false;
  // DCA between daughters
  double dDCA = aodV0->DcaV0Daughters();
  if (dDCA > fV0DcaBetweenDaughtersMax) return false;
  double dPt = aodV0->Pt();
  if (dPt < fV0PtMin ) return false;
  double dRapidity = aodV0->RapLambda();
  if (TMath::Abs(dRapidity) > fV0RapidityMax) return false;
  return kTRUE;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDAngCorr::IsGoodDaughterTrack(const AliAODTrack *track)
{
  // TPC refit
  if (!track->IsOn(AliAODTrack::kTPCrefit)) return false;
  // No kinks
  if (int(track->GetProdVertex()->GetType()) == AliAODVertex::kKink) return false;
  // Maximum value of transverse momentum
  double dPt = track->Pt();
  if (dPt > fDaughtersPtMax) return false;
  // Maximum value of pseudorapidity
  double dEta = track->Eta();
  if (TMath::Abs(dEta) > fDaughtersEtaMax) return false;
  // Minimum number of clusters
  float nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < fDaughtersTPCNclsMin) return false;
  // Findable clusters > 0
  int findable = track->GetTPCNclsF();
  if (findable <= 0) return false;
  // [number of crossed rows]>0.8  [number of findable clusters].
  if (nCrossedRowsTPC/findable < 0.8) return false;
  return true;
}

//---------------------------------------------------

int AliAnalysisTaskPIDAngCorr::GetLambdaCode(const AliAODTrack *pTrack, const AliAODTrack *nTrack)
{
  bool IsLambda     = kFALSE;
  bool IsAntiLambda = kFALSE;
  int  code = 0;

  //Λ-->(p+)+(π-)
  float nSigTPCPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton));//TPC p+
  float nSigTPCNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion));//TPC π-
  //(Λ-)-->(p-)+(π+)
  float nSigTPCPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion));//TPC π+
  float nSigTPCNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton));//TPC p-

  IsLambda     = (nSigTPCPosProton < fV0PosProtonTPCNsigma) && (nSigTPCNegPion < fV0NegPionTPCNsigma);
  IsAntiLambda = (nSigTPCNegProton < fV0NegProtonTPCNsigma) && (nSigTPCPosPion < fV0PosPionTPCNsigma);

  if (IsV0DaughterUseTOF) {
    float nSigTOFPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack, AliPID::kProton));//TOF p+
    float nSigTOFNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack, AliPID::kPion));//TOF π-
    float nSigTOFPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack, AliPID::kPion));//TOF π+
    float nSigTOFNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack, AliPID::kProton));//TOF p-

    IsLambda     *= (nSigTOFPosProton < fV0PosProtonTOFNsigma && nSigTOFNegPion < fV0NegPionTOFNsigma);
    IsAntiLambda *= (nSigTOFNegProton < fV0NegProtonTOFNsigma && nSigTOFPosPion < fV0PosPionTOFNsigma);
  }

  if (IsLambda)     code =  3122;
  if (IsAntiLambda) code = -3122;
  if (IsLambda && IsAntiLambda) code = 0;
  return code;
}

//---------------------------------------------------

double AliAnalysisTaskPIDAngCorr::RangePhi(double phi)
{
  while (phi >= 1.5*TMath::Pi())
  phi -= TMath::TwoPi();
  while (phi < -0.5*TMath::Pi())
  phi += TMath::TwoPi();
  return phi;
}