#ifndef AliAnalysisTaskPIDAngCorr_cxx
#define AliAnalysisTaskPIDAngCorr_cxx
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliAODv0.h"
#include "AliAODZDC.h"
#include "TProfile3D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TComplex.h"
#include "TList.h"
#include "TFile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "AliEventCuts.h"

class AliAnalysisTaskPIDAngCorr : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskPIDAngCorr();
  AliAnalysisTaskPIDAngCorr(const char *name);
  virtual ~AliAnalysisTaskPIDAngCorr();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  //Global
  void SetDebug(int debug)                 {this->fDebug   =   debug;}
  void SetTrigger(TString trigger)         {this->fTrigger = trigger;}
  void SetPeriod(TString period)           {this->fPeriod  =  period;}
  //Event
  void SetVzCut(double vzCut)              {this->fVzCut       =       vzCut;}
  void SetCentCut(float centDiffCut)       {this->fCentDiffCut = centDiffCut;}
  //Track
  void SetFilterBit(int filterBit)         {this->fFilterBit   =   filterBit;}
  void SetNclsCut(int nclsCut)             {this->fNclsCut     =     nclsCut;}
  void SetChi2Max(float chi2Max)           {this->fChi2Max     =     chi2Max;}
  void SetChi2Min(float chi2Min)           {this->fChi2Min     =     chi2Min;}
  void SetDCAcutZ(float dcaCutz)           {this->fDcaCutz     =     dcaCutz;}
  void SetDCAcutXY(float dcaCutxy)         {this->fDcaCutxy    =    dcaCutxy;}
  void SetPtMin(float ptMin)               {this->fPtMin       =       ptMin;}
  void SetPtMax(float ptMax)               {this->fPtMax       =       ptMax;}
  void SetProtonPtMax(double protonPtMax)  {this->fProtonPtMax = protonPtMax;}
  //PID
  void SetNSigmaTPCCut(double nSigmaTPC) {this->fNSigmaTPCCut = nSigmaTPC;}
  void SetNSigmaTOFCut(double nSigmaTOF) {this->fNSigmaTOFCut = nSigmaTOF;}
  //V0
  void SetV0PtMin(double v0PtMin)                                   {this->fV0PtMin                  =                  v0PtMin;}
  void SetV0CPAMin(double v0CPAMin)                                 {this->fV0CPAMin                 =                 v0CPAMin;}
  void SetV0RapidityMax(double v0RapidityMax)                       {this->fV0RapidityMax            =            v0RapidityMax;}
  void SetV0DecayLengthMax(double v0DecayLengthMax)                 {this->fV0DecayLengthMax         =         v0DecayLengthMax;}
  void SetV0DecayLengthMin(double v0DecayLengthMin)                 {this->fV0DecayLengthMin         =         v0DecayLengthMin;}
  void SetV0DCAToPrimVtxMax(double v0DCAToPrimVtxMax)               {this->fV0DCAToPrimVtxMax        =        v0DCAToPrimVtxMax;}
  void SetV0DcaBetweenDaughtersMax(double v0DcaBetweenDaughtersMax) {this->fV0DcaBetweenDaughtersMax = v0DcaBetweenDaughtersMax;}
  //V0 Daughter Cut
  void SetDaughtersPtMax(double daughtersPtMax)                     {this->fDaughtersPtMax           =           daughtersPtMax;}
  void SetDaughtersEtaMax(double daughtersEtaMax)                   {this->fDaughtersEtaMax          =          daughtersEtaMax;}
  void SetDaughtersTPCNclsMin(double daughtersTPCNclsMin)           {this->fDaughtersTPCNclsMin      =      daughtersTPCNclsMin;}
  void SetDaughtersDCAToPrimVtxMin(double daughtersDCAToPrimVtxMin) {this->fDaughtersDCAToPrimVtxMin = daughtersDCAToPrimVtxMin;}
  void SetPosProtonTPCNsigma(float V0PosProtonTPCNsigma)            {this->fV0PosProtonTPCNsigma     =     V0PosProtonTPCNsigma;}
  void SetNegPionTPCNsigma(float   V0NegPionTPCNsigma)              {this->fV0NegPionTPCNsigma       =       V0NegPionTPCNsigma;}
  void SetNegProtonTPCNsigma(float V0NegProtonTPCNsigma)            {this->fV0NegProtonTPCNsigma     =     V0NegProtonTPCNsigma;}
  void SetPosPionTPCNsigma(float   V0PosPionTPCNsigma)              {this->fV0PosPionTPCNsigma       =       V0PosPionTPCNsigma;}
  void IfDaughtersPIDUseTOF(bool daughterPIDUseTOF)                 {this->IsV0DaughterUseTOF        =        daughterPIDUseTOF;}
  void SetPosProtonTOFNsigma(float V0PosProtonTOFNsigma)            {this->fV0PosProtonTOFNsigma     =     V0PosProtonTOFNsigma;}
  void SetNegPionTOFNsigma(float   V0NegPionTOFNsigma)              {this->fV0NegPionTOFNsigma       =       V0NegPionTOFNsigma;}
  void SetNegProtonTOFNsigma(float V0NegProtonTOFNsigma)            {this->fV0NegProtonTOFNsigma     =     V0NegProtonTOFNsigma;}
  void SetPosPionTOFNsigma(float   V0PosPionTOFNsigma)              {this->fV0PosPionTOFNsigma       =       V0PosPionTOFNsigma;}
  //Lambda Mass Cut
  void SetMassMean(double massMean)           {this->fMassMean      =      massMean;}
  void SetLambdaMassCut(double lambdaMassCut) {this->fLambdaMassCut = lambdaMassCut;}

private:
  ////////////////////////
  // Procedural function
  ////////////////////////
  bool    GetVZEROPlane();
  bool      GetZDCPlane();
  bool GetZDCPlaneLsFit();
  bool       LoopTracks();
  bool      GetTPCPlane();
  bool          LoopV0s();
  bool             Pair();
  bool          MixPair();
  void     ResetVectors();

  ////////////////////////
  //Functional function
  ////////////////////////
  // Read in
  //bool      LoadCalibHistForThisRun();
  // Pile-up
  bool           RejectEvtMultComp();
  bool              RejectEvtTFFit();
  bool      RejectEvtTPCITSfb32TOF();
  bool              AODPileupCheck();
  bool           PileUpMultiVertex();
  bool              RemovalForRun1();
  double                  GetWDist(const AliVVertex* v0, const AliVVertex* v1);
  // Track
  bool              AcceptAODTrack(AliAODTrack *track);
  bool          CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck);
  // V0
  bool                    IsGoodV0(AliAODv0 *aodV0);
  bool         IsGoodDaughterTrack(const AliAODTrack *track);
  int                GetLambdaCode(const AliAODTrack *pTrack, const AliAODTrack *ntrack);

  //
  double                  RangePhi(double phi);

  //////////////////////
  // Cuts and options //
  //////////////////////
  //Global
  int                       fDebug; // debug level controls amount of output statements
  TString                 fTrigger; //
  TString                  fPeriod; // period

  //Event
  double                    fVzCut; // vz cut
  float               fCentDiffCut; // centrality restriction for V0M and TRK

  //Calib
  bool              IsVZEROCalibOn; // switch for VZERO qn calib
  bool                IsZDCCalibOn;
  bool                IsTPCCalibOn; // switch for TPC qn calib
  bool                   IsQAVZERO; // flag for V0 qn QA
  bool                     IsQAZDC;
  bool                     IsQATPC;

  //Plane
  float                fPlanePtMin;
  float                fPlanePtMax;
  float                 fEtaGapPos; // value for the Eta Gap Pos
  float                 fEtaGapNeg;

  //Track
  int                   fFilterBit; // AOD filter bit selection
  int                     fNclsCut; // ncls cut for all tracks
  float                   fChi2Max; // upper limmit for chi2
  float                   fChi2Min; // lower limmit for chi2
  float                   fDcaCutz; // dcaz cut for all tracks
  float                  fDcaCutxy; // dcaxy cut for all tracks
  float                     fPtMin; // minimum pt for track
  float                     fPtMax; // maximum pt for track
  bool                     IsDoNUE; // switch for NUE
  bool                     IsDoNUA; // switch for NUA
  double              fProtonPtMax; // Max pt for proton
  //PID
  float             fNSigmaTPCCut;
  float             fNSigmaTOFCut;

  //V0
  double                  fV0PtMin; //
  double                 fV0CPAMin; //
  double            fV0RapidityMax; //
  double         fV0DecayLengthMin; //
  double         fV0DecayLengthMax; //
  double        fV0DCAToPrimVtxMax; //
  double fV0DcaBetweenDaughtersMax; //

  //V0 Daughter
  double           fDaughtersPtMax; //
  double          fDaughtersEtaMax; //
  double      fDaughtersTPCNclsMin; //
  double fDaughtersDCAToPrimVtxMin; //
  float      fV0PosProtonTPCNsigma; //
  float        fV0NegPionTPCNsigma; //
  float      fV0NegProtonTPCNsigma; //
  float        fV0PosPionTPCNsigma; //
  bool          IsV0DaughterUseTOF; //
  float      fV0PosProtonTOFNsigma; //
  float        fV0NegPionTOFNsigma; //
  float      fV0NegProtonTOFNsigma; //
  float        fV0PosPionTOFNsigma; //

  //Lambda Mass
  double            fLambdaMassCut; //


  // Global Variables Unchanged in an Evt
  double                 fMassMean; //
  const float              fEtaCut; // eta cut
  const float             fDedxCut; // dedx cut

  ///////////////////The following files are from the data//////////////////////////////////
  /////////////
  // Handles //
  /////////////
  AliAODEvent*                fAOD; // aod Event
  AliPIDResponse*     fPIDResponse; // PID Handler
  AliAnalysisUtils*         fUtils; // Event Selection Options

  ////////////////////////////////
  // Global Variables from data //
  ////////////////////////////////
  std::map<int,int>*    runNumList;
  double                fVertex[3];
  int                      fRunNum; // runnumber
  int                   fOldRunNum;
  int                   fRunNumBin; // runnumer bin; 10:139510...; 11:170387...; 15HIR:246994...
  int                       fVzBin; // vertex z bin
  int                     fCentBin; // centrality bin: 0-10
  double                     fCent; // value of centrality

  //PID Vector
  //p+ p-
  std::vector<int>           vecID[6]; 
  std::vector<float>         vecPt[6];
  std::vector<float>         vecPhi[6];

  //Lambda Vector
  std::vector<int>      vecDaughterPosID[6]; // Pos Daughter ID
  std::vector<int>      vecDaughterNegID[6]; // Neg Daughter ID


  //Buffer for Mix Event
  //p+ p- lambda lambdabar
  std::vector<int>           vecNumMixBuffer[6];
  std::vector<float>         vecPhiMixBuffer[6];


  ///////////////////The following files are read from external sources////////////////////
  ////////////////////////
  // Pile up Function
  ////////////////////////
  TF1*                       fSPDCutPU;
  TF1*                        fV0CutPU;
  TF1*                    fCenCutLowPU;
  TF1*                   fCenCutHighPU;
  TF1*                      fMultCutPU;

  ///////////////////The following files will be saved//////////////////////////////////
  ////////////
  //QA Plots//
  ////////////
  // Event-wise QA
  TList*             fOutputList;
  TH1D*                fEvtCount;
  TH1I*           fHistRunNumBin;
  TH1D*             fHistCent[2];
  TH1D*               fHistVz[2];
  TH2D*         fHist2DCentQA[8];
  TH2D*     fHist2DMultCentQA[2];
  TH2D*     fHist2DMultMultQA[6];
  // Track-wise QA
  TH1D*             fHistPt;
  TH1D*            fHistEta;
  TH1D*          fHistNhits;
  TH2D*        fHist2DPDedx;
  TH1D*          fHistDcaXY;
  TH1D*           fHistDcaZ;
  TH1D*         fHistPhi[2];
  TH2D*    fHist2DEtaPhi[2];
  
  //V0s QA
  TH1D*                 fHistV0Pt; // Raw V0s' pT
  TH1D*                fHistV0Eta; // Raw V0s' eta
  TH1D*    fHistV0DcatoPrimVertex; // Raw V0s' DcatoPV
  TH1D*                fHistV0CPA; // Raw V0s' CPA(cosine pointing angle)
  TH1D*        fHistV0DecayLength; // Raw V0s' DecayLength
  ///Lambda QA
  //[0]:Before the Mass Cut [1]:After the Mass Cut
  TH1D*                        fHistLambdaPt[2]; //
  TH1D*                       fHistLambdaEta[2]; //
  TH1D*                       fHistLambdaPhi[2]; //
  TH1D*           fHistLambdaDcaToPrimVertex[2]; //
  TH1D*                       fHistLambdaCPA[2]; //
  TH1D*               fHistLambdaDecayLength[2]; //
  TH1D*                      fHistLambdaMass[2]; //
  TH1D*                    fHistAntiLambdaPt[2]; //
  TH1D*                   fHistAntiLambdaEta[2]; //
  TH1D*                   fHistAntiLambdaPhi[2]; //
  TH1D*       fHistAntiLambdaDcaToPrimVertex[2]; //
  TH1D*                   fHistAntiLambdaCPA[2]; //
  TH1D*           fHistAntiLambdaDecayLength[2]; //
  TH1D*                  fHistAntiLambdaMass[2]; //
  TProfile*           fProfileLambdaMassVsPt[2]; //
  TProfile*       fProfileAntiLambdaMassVsPt[2]; //

  /////////////
  // Results //
  /////////////

  ///Lambda-X correlators
  TH1F*       fHistSig[6][6];
  TH1F*       fHistBkg[6][6];

  AliAnalysisTaskPIDAngCorr(const AliAnalysisTaskPIDAngCorr&);
  AliAnalysisTaskPIDAngCorr& operator=(const AliAnalysisTaskPIDAngCorr&);

  ClassDef(AliAnalysisTaskPIDAngCorr, 1);
};

#endif