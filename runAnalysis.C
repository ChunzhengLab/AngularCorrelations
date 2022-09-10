// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly
#include "AliAnalysisTaskPIDAngCorr.h"
#include "TInterpreter.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "TMacro.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliMultSelectionTask.h"
#include "TSystem.h"
#include "TChain.h"
#include "AliAnalysisAlien.h"
#include "AliPhysicsSelectionTask.h"


void runAnalysis()
{
    // Dataset
    TString dataset = "LHC18q";
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t local      = kFALSE;
    Bool_t gridTest   = kFALSE;
    Bool_t gridRun    = kFALSE;
    Bool_t gridMerge  = kFALSE;
    Bool_t localMerge = kTRUE;

    if (gridMerge || localMerge) gridRun = kFALSE;
    if (gridRun) {
      gridTest = kFALSE;
      local = kFALSE;
    }

    Int_t groupNumber = 1;
    
    // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif
     
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    TMacro PIDadd(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
    AliAnalysisTaskPIDResponse* PIDresponseTask = reinterpret_cast<AliAnalysisTaskPIDResponse*>(PIDadd.Exec());
   
    // if (dataset.EqualTo("LHC15o") || dataset.EqualTo("LHC18q") || dataset.EqualTo("LHC18r")) {
    //   // TMacro physicsSelection(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"));
    //   // AliPhysicsSelectionTask* physicsSelectionTask = reinterpret_cast<AliPhysicsSelectionTask*>(physicsSelection.Exec());
    //   TMacro multSelection(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
    //   AliMultSelectionTask* multSelectionTask = reinterpret_cast<AliMultSelectionTask*>(multSelection.Exec());
    // }
    
    // compile the class and load the add task macro
    // here we have to differentiate between using the just-in-time compiler
    // from root6, or the interpreter of root5
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AliAnalysisTaskPIDAngCorr.cxx++g");
    AliAnalysisTaskPIDAngCorr *task = nullptr;
    if(local)  task = reinterpret_cast<AliAnalysisTaskPIDAngCorr*>(gInterpreter->ExecuteMacro("AddTaskPIDAngCorr.C"));
    else       task = reinterpret_cast<AliAnalysisTaskPIDAngCorr*>(gInterpreter->ExecuteMacro("AddTaskPIDAngCorr.C"));
#else
    gROOT->LoadMacro("AliAnalysisTaskPIDAngCorr.cxx++g");
    gROOT->LoadMacro("AddTaskPIDAngCorr.C");
    AliAnalysisTaskPIDAngCorr *task = AddTaskPIDAngCorr();
#endif


    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        // add a few files to the chain (change this so that your local files are added)
        if(dataset.EqualTo("LHC10h")) chain->Add("/Users/wangchunzheng/alice/data/2010/LHC10h/000139510/ESDs/pass2/AOD160/0247/AliAOD.root");
        if(dataset.EqualTo("LHC15o")) chain->Add("/Users/wangchunzheng/alice/data/2015/LHC15o/000245151/pass2/AOD252/0008/AliAOD.root");
        if(dataset.EqualTo("LHC18q")) chain->Add("/Users/wangchunzheng/alice/data/2018/LHC18q/000295588/pass3/AOD252/AOD/001/AliAOD.root");
        if(dataset.EqualTo("LHC18r")) chain->Add("/Users/wangchunzheng/alice/data/2018/LHC18r/000296691/pass3/AOD252/0001/AliAOD.root");
        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskPIDAngCorr.cxx AliAnalysisTaskPIDAngCorr.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskPIDAngCorr.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20220316_ROOT6-1");
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");

        if (dataset.EqualTo("LHC10h")) {
          //10h
          // select the input data
          alienHandler->SetGridDataDir("/alice/data/2010/LHC10h");
          alienHandler->SetDataPattern("ESDs/pass2/AOD160/*/AliAOD.root");
          // MC has no prefix, data has prefix 000
          alienHandler->SetRunPrefix("000");
          // runnumber
          //1
          if(groupNumber == 1) {
            alienHandler->AddRunNumber(139510);
            alienHandler->AddRunNumber(139507);
            alienHandler->AddRunNumber(139505);
            // alienHandler->AddRunNumber(139503);
            // alienHandler->AddRunNumber(139465);
            // alienHandler->AddRunNumber(139438);
            // alienHandler->AddRunNumber(139437);
            // alienHandler->AddRunNumber(139360);
            // alienHandler->AddRunNumber(139329);
            // alienHandler->AddRunNumber(139328);
            // alienHandler->AddRunNumber(139314);
            // alienHandler->AddRunNumber(139310);
            }
          if(groupNumber == 2) {
            alienHandler->AddRunNumber(139309);
            alienHandler->AddRunNumber(139173);
            alienHandler->AddRunNumber(139107);
            alienHandler->AddRunNumber(139105);
            alienHandler->AddRunNumber(139038);
            alienHandler->AddRunNumber(139037);
            alienHandler->AddRunNumber(139036);
            alienHandler->AddRunNumber(139029);
            alienHandler->AddRunNumber(139028);
            alienHandler->AddRunNumber(138872);
            alienHandler->AddRunNumber(138871);
            alienHandler->AddRunNumber(138870);
          }
          if(groupNumber == 3) {
            alienHandler->AddRunNumber(138837);
            alienHandler->AddRunNumber(138732);
            alienHandler->AddRunNumber(138730);
            alienHandler->AddRunNumber(138666);
            alienHandler->AddRunNumber(138662);
            alienHandler->AddRunNumber(138653);
            alienHandler->AddRunNumber(138652);
            alienHandler->AddRunNumber(138638);
            alienHandler->AddRunNumber(138624);
            alienHandler->AddRunNumber(138621);
            alienHandler->AddRunNumber(138583);
            alienHandler->AddRunNumber(138582);
          }
          if(groupNumber == 4) {
            alienHandler->AddRunNumber(138578);
            alienHandler->AddRunNumber(138534);
            alienHandler->AddRunNumber(138469);
            alienHandler->AddRunNumber(138442);
            alienHandler->AddRunNumber(138439);
            alienHandler->AddRunNumber(138438);
            alienHandler->AddRunNumber(138396);
            alienHandler->AddRunNumber(138364);
            alienHandler->AddRunNumber(138275);
            alienHandler->AddRunNumber(138225);
            alienHandler->AddRunNumber(138201);
            alienHandler->AddRunNumber(138197);
          }
          if(groupNumber == 5) {
            alienHandler->AddRunNumber(138192);
            alienHandler->AddRunNumber(138190);
            alienHandler->AddRunNumber(137848);
            alienHandler->AddRunNumber(137844);
            alienHandler->AddRunNumber(137752);
            alienHandler->AddRunNumber(137751);
            alienHandler->AddRunNumber(137724);
            alienHandler->AddRunNumber(137722);
            alienHandler->AddRunNumber(137718);
            alienHandler->AddRunNumber(137704);
            alienHandler->AddRunNumber(137693);
            alienHandler->AddRunNumber(137692);
          }
          if(groupNumber == 6) {
            alienHandler->AddRunNumber(137691);
            alienHandler->AddRunNumber(137686);
            alienHandler->AddRunNumber(137685);
            alienHandler->AddRunNumber(137639);
            alienHandler->AddRunNumber(137638);
            alienHandler->AddRunNumber(137608);
            alienHandler->AddRunNumber(137595);
            alienHandler->AddRunNumber(137549);
            alienHandler->AddRunNumber(137546);
            alienHandler->AddRunNumber(137544);
            alienHandler->AddRunNumber(137541);
            alienHandler->AddRunNumber(137539);
          }
          if(groupNumber == 7) {
            alienHandler->AddRunNumber(137531);
            alienHandler->AddRunNumber(137530);
            alienHandler->AddRunNumber(137443);
            alienHandler->AddRunNumber(137441);
            alienHandler->AddRunNumber(137440);
            alienHandler->AddRunNumber(137439);
            alienHandler->AddRunNumber(137434);
            alienHandler->AddRunNumber(137432);
            alienHandler->AddRunNumber(137431);
            alienHandler->AddRunNumber(137430);
            alienHandler->AddRunNumber(137243);
            alienHandler->AddRunNumber(137236);
            alienHandler->AddRunNumber(137235);
          }
          if(groupNumber == 8) {
            alienHandler->AddRunNumber(137232);
            alienHandler->AddRunNumber(137231);
            alienHandler->AddRunNumber(137230);
            alienHandler->AddRunNumber(137162);
            alienHandler->AddRunNumber(137161);
            alienHandler->AddRunNumber(138732);
          }
        }
        if (dataset.EqualTo("LHC15o")) {
          // select the input data
          alienHandler->SetGridDataDir("/alice/data/2015/LHC15o");
          alienHandler->SetDataPattern("pass2/AOD252/AOD/*/AliAOD.root");
          // MC has no prefix, data has prefix 000
          alienHandler->SetRunPrefix("000");
          if(groupNumber == 1) {
            alienHandler->AddRunNumber(244917);
            alienHandler->AddRunNumber(244918);
            alienHandler->AddRunNumber(244975);
            alienHandler->AddRunNumber(244980);
            alienHandler->AddRunNumber(244982);
            alienHandler->AddRunNumber(244983);
            alienHandler->AddRunNumber(245061);
            alienHandler->AddRunNumber(245064);
          }
          if(groupNumber == 2) {
            alienHandler->AddRunNumber(245066);
            alienHandler->AddRunNumber(245068);
            alienHandler->AddRunNumber(245145);
            alienHandler->AddRunNumber(245146);
            alienHandler->AddRunNumber(245148);
            alienHandler->AddRunNumber(245151);
            alienHandler->AddRunNumber(245152);
            alienHandler->AddRunNumber(245231);
          }
          if(groupNumber == 3) {
            alienHandler->AddRunNumber(245232);
            alienHandler->AddRunNumber(245233);
            alienHandler->AddRunNumber(245259);
            alienHandler->AddRunNumber(245343);
            alienHandler->AddRunNumber(245345);
            alienHandler->AddRunNumber(245346);
            alienHandler->AddRunNumber(245347);
            alienHandler->AddRunNumber(245349);
          }
          if(groupNumber == 4) {
            alienHandler->AddRunNumber(245353);
            alienHandler->AddRunNumber(245396);
            alienHandler->AddRunNumber(245397);
            alienHandler->AddRunNumber(245401);
            alienHandler->AddRunNumber(245407);
            alienHandler->AddRunNumber(245409);
            alienHandler->AddRunNumber(245410);
            alienHandler->AddRunNumber(245411);
            alienHandler->AddRunNumber(245439);
          }
          if(groupNumber == 5) {
            alienHandler->AddRunNumber(245441);
            alienHandler->AddRunNumber(245446);
            alienHandler->AddRunNumber(245450);
            alienHandler->AddRunNumber(245452);
            alienHandler->AddRunNumber(245453);
            alienHandler->AddRunNumber(245454);
            alienHandler->AddRunNumber(245496);
            alienHandler->AddRunNumber(245497);
          }
          if(groupNumber == 6) {
            alienHandler->AddRunNumber(245501);
            alienHandler->AddRunNumber(245504);
            alienHandler->AddRunNumber(245505);
            alienHandler->AddRunNumber(245507);
            alienHandler->AddRunNumber(245535);
            alienHandler->AddRunNumber(245540);
            alienHandler->AddRunNumber(245542);
            alienHandler->AddRunNumber(245543);
          }
          if(groupNumber == 7) {
            alienHandler->AddRunNumber(245544);
            alienHandler->AddRunNumber(245545);
            alienHandler->AddRunNumber(245554);
            alienHandler->AddRunNumber(245683);
            alienHandler->AddRunNumber(245692);
            alienHandler->AddRunNumber(245700);
            alienHandler->AddRunNumber(245702);
            alienHandler->AddRunNumber(245705);
          }
          if(groupNumber == 8) {
            alienHandler->AddRunNumber(245729);
            alienHandler->AddRunNumber(245731);
            alienHandler->AddRunNumber(245738);
            alienHandler->AddRunNumber(245752);
            alienHandler->AddRunNumber(245759);
            alienHandler->AddRunNumber(245766);
            alienHandler->AddRunNumber(245775);
            alienHandler->AddRunNumber(245785);
          }
          if(groupNumber == 8) {
            alienHandler->AddRunNumber(245793);
            alienHandler->AddRunNumber(245829);
            alienHandler->AddRunNumber(245831);
            alienHandler->AddRunNumber(245833);
            alienHandler->AddRunNumber(245923);
            alienHandler->AddRunNumber(245949);
            alienHandler->AddRunNumber(245952);
            alienHandler->AddRunNumber(245954);
          }
          if(groupNumber == 9) {
            alienHandler->AddRunNumber(245963);
            alienHandler->AddRunNumber(246001);
            alienHandler->AddRunNumber(246003);
            alienHandler->AddRunNumber(246012);
            alienHandler->AddRunNumber(246036);
            alienHandler->AddRunNumber(246037);
            alienHandler->AddRunNumber(246042);
            alienHandler->AddRunNumber(246048);
          }
          if(groupNumber == 10) {
            alienHandler->AddRunNumber(246049);
            alienHandler->AddRunNumber(246052);
            alienHandler->AddRunNumber(246053);
            alienHandler->AddRunNumber(246087);
            alienHandler->AddRunNumber(246089);
            alienHandler->AddRunNumber(246113);
            alienHandler->AddRunNumber(246115);
            alienHandler->AddRunNumber(246148);
          }
          if(groupNumber == 11) {
            alienHandler->AddRunNumber(246151);
            alienHandler->AddRunNumber(246152);
            alienHandler->AddRunNumber(246153);
            alienHandler->AddRunNumber(246178);
            alienHandler->AddRunNumber(246180);
            alienHandler->AddRunNumber(246181);
            alienHandler->AddRunNumber(246182);
            alienHandler->AddRunNumber(246185);
          }
          if(groupNumber == 12) {
            alienHandler->AddRunNumber(246217);
            alienHandler->AddRunNumber(246222);
            alienHandler->AddRunNumber(246225);
            alienHandler->AddRunNumber(246271);
            alienHandler->AddRunNumber(246272);
            alienHandler->AddRunNumber(246275);
            alienHandler->AddRunNumber(246276);
            alienHandler->AddRunNumber(246390);
          }
          if(groupNumber == 13) {
            alienHandler->AddRunNumber(246391);
            alienHandler->AddRunNumber(246392);
            alienHandler->AddRunNumber(246424);
            alienHandler->AddRunNumber(246428);
            alienHandler->AddRunNumber(246431);
            alienHandler->AddRunNumber(246434);
            alienHandler->AddRunNumber(246487);
            alienHandler->AddRunNumber(246488);
          }
          if(groupNumber == 14) {
            alienHandler->AddRunNumber(246493);
            alienHandler->AddRunNumber(246495);
            alienHandler->AddRunNumber(246540);
            alienHandler->AddRunNumber(246543);
            alienHandler->AddRunNumber(246553);
            alienHandler->AddRunNumber(246567);
            alienHandler->AddRunNumber(246568);
            alienHandler->AddRunNumber(246575);
          }
          if(groupNumber == 15) {
            alienHandler->AddRunNumber(246583);
            alienHandler->AddRunNumber(246648);
            alienHandler->AddRunNumber(246671);
            alienHandler->AddRunNumber(246675);
            alienHandler->AddRunNumber(246676);
            alienHandler->AddRunNumber(246750);
            alienHandler->AddRunNumber(246751);
            alienHandler->AddRunNumber(246757);
            alienHandler->AddRunNumber(246758);
          }
          if(groupNumber == 16) {
            alienHandler->AddRunNumber(246759);
            alienHandler->AddRunNumber(246760);
            alienHandler->AddRunNumber(246763);
            alienHandler->AddRunNumber(246765);
            alienHandler->AddRunNumber(246766);
            alienHandler->AddRunNumber(246804);
            alienHandler->AddRunNumber(246805);
            alienHandler->AddRunNumber(246807);
          }
          if(groupNumber == 17) {
            alienHandler->AddRunNumber(246808);
            alienHandler->AddRunNumber(246809);
            alienHandler->AddRunNumber(246810);
            alienHandler->AddRunNumber(246844);
            alienHandler->AddRunNumber(246845);
            alienHandler->AddRunNumber(246846);
            alienHandler->AddRunNumber(246847);
            alienHandler->AddRunNumber(246851);
          }
          if(groupNumber == 18) {
            alienHandler->AddRunNumber(246858);
            alienHandler->AddRunNumber(246859);
            alienHandler->AddRunNumber(246864);
            alienHandler->AddRunNumber(246865);
            alienHandler->AddRunNumber(246867);
            alienHandler->AddRunNumber(246870);
            alienHandler->AddRunNumber(246871);
            alienHandler->AddRunNumber(246928);
            alienHandler->AddRunNumber(246945);
          }
          if(groupNumber == 19) {
            alienHandler->AddRunNumber(246948);
            alienHandler->AddRunNumber(246980);
            alienHandler->AddRunNumber(246982);
            alienHandler->AddRunNumber(246984);
            alienHandler->AddRunNumber(246989);
            alienHandler->AddRunNumber(246991);
            alienHandler->AddRunNumber(246994);
          }
        }
        if (dataset.EqualTo("LHC18q")) {
          // select the input data
          alienHandler->SetGridDataDir("/alice/data/2018/LHC18q");
          alienHandler->SetDataPattern("pass3/AOD252/AOD/*/AliAOD.root");
          // MC has no prefix, data has prefix 000
          alienHandler->SetRunPrefix("000");
          if(groupNumber == 1) {
            alienHandler->AddRunNumber(296269);
            alienHandler->AddRunNumber(296270); 
            alienHandler->AddRunNumber(296547);
            alienHandler->AddRunNumber(296472);
            alienHandler->AddRunNumber(295581);
            alienHandler->AddRunNumber(295584);
            alienHandler->AddRunNumber(295585);
            alienHandler->AddRunNumber(295586);
          }
          if(groupNumber == 2) {
            alienHandler->AddRunNumber(295587);
            alienHandler->AddRunNumber(295588);
            alienHandler->AddRunNumber(295589);
            alienHandler->AddRunNumber(295610);
            alienHandler->AddRunNumber(295611);
            alienHandler->AddRunNumber(295612);
            alienHandler->AddRunNumber(295615);
            alienHandler->AddRunNumber(295665);
          }
          if(groupNumber == 3) {
            alienHandler->AddRunNumber(295666);
            alienHandler->AddRunNumber(295667);
            alienHandler->AddRunNumber(295668);
            alienHandler->AddRunNumber(295671);
            alienHandler->AddRunNumber(295673);
            alienHandler->AddRunNumber(295675);
            alienHandler->AddRunNumber(295676);
            alienHandler->AddRunNumber(295677);
          }
          if(groupNumber == 4) {
            alienHandler->AddRunNumber(295712);
            alienHandler->AddRunNumber(295714);
            alienHandler->AddRunNumber(295717);
            alienHandler->AddRunNumber(295718);
            alienHandler->AddRunNumber(295719);
            alienHandler->AddRunNumber(295720);
            alienHandler->AddRunNumber(295721);
            alienHandler->AddRunNumber(295723);
          }
          if(groupNumber == 5) {
            alienHandler->AddRunNumber(295725);
            alienHandler->AddRunNumber(295753);
            alienHandler->AddRunNumber(295754);
            alienHandler->AddRunNumber(295755);
            alienHandler->AddRunNumber(295756);
            alienHandler->AddRunNumber(295758);
            alienHandler->AddRunNumber(295759);
            alienHandler->AddRunNumber(295762);
            alienHandler->AddRunNumber(295763);
          }
          if(groupNumber == 6) {
            alienHandler->AddRunNumber(295786);
            alienHandler->AddRunNumber(295788);
            alienHandler->AddRunNumber(295791);
            alienHandler->AddRunNumber(295816);
            alienHandler->AddRunNumber(295818);
            alienHandler->AddRunNumber(295819);
            alienHandler->AddRunNumber(295822);
            alienHandler->AddRunNumber(295825);
            alienHandler->AddRunNumber(295826);
          }
          if(groupNumber == 7) {
            alienHandler->AddRunNumber(295829);
            alienHandler->AddRunNumber(295831);
            alienHandler->AddRunNumber(295853);
            alienHandler->AddRunNumber(295854);
            alienHandler->AddRunNumber(295855);
            alienHandler->AddRunNumber(295856);
            alienHandler->AddRunNumber(295859);
            alienHandler->AddRunNumber(295860);
          }
          if(groupNumber == 8) {
            alienHandler->AddRunNumber(295861);
            alienHandler->AddRunNumber(295881);
            alienHandler->AddRunNumber(295908);
            alienHandler->AddRunNumber(295909);
            alienHandler->AddRunNumber(295910);
            alienHandler->AddRunNumber(295913);
            alienHandler->AddRunNumber(295936);
            alienHandler->AddRunNumber(295937);
            alienHandler->AddRunNumber(295941);
          }
          if(groupNumber == 9) {
            alienHandler->AddRunNumber(295942);
            alienHandler->AddRunNumber(295943);
            alienHandler->AddRunNumber(295945);
            alienHandler->AddRunNumber(295947);
            alienHandler->AddRunNumber(296016);
            alienHandler->AddRunNumber(296060);
            alienHandler->AddRunNumber(296061);
            alienHandler->AddRunNumber(296062);
            alienHandler->AddRunNumber(296063);
          }
          if(groupNumber == 10) {
            alienHandler->AddRunNumber(296065);
            alienHandler->AddRunNumber(296066);
            alienHandler->AddRunNumber(296068);
            alienHandler->AddRunNumber(296074);
            alienHandler->AddRunNumber(296123);
            alienHandler->AddRunNumber(296132);
            alienHandler->AddRunNumber(296133);
            alienHandler->AddRunNumber(296134);
            alienHandler->AddRunNumber(296135);
          }
          if(groupNumber == 11) {
            alienHandler->AddRunNumber(296142);
            alienHandler->AddRunNumber(296143);
            alienHandler->AddRunNumber(296191);
            alienHandler->AddRunNumber(296192);
            alienHandler->AddRunNumber(296195);
            alienHandler->AddRunNumber(296196);
            alienHandler->AddRunNumber(296197);
            alienHandler->AddRunNumber(296198);
          }
          if(groupNumber == 12) {
            alienHandler->AddRunNumber(296240);
            alienHandler->AddRunNumber(296241);
            alienHandler->AddRunNumber(296242);
            alienHandler->AddRunNumber(296243);
            alienHandler->AddRunNumber(296244);
            alienHandler->AddRunNumber(296246);
            alienHandler->AddRunNumber(296247);
            alienHandler->AddRunNumber(296269);
          }
          if(groupNumber == 13) {
            alienHandler->AddRunNumber(296270);
            alienHandler->AddRunNumber(296273);
            alienHandler->AddRunNumber(296279);
            alienHandler->AddRunNumber(296280);
            alienHandler->AddRunNumber(296303);
            alienHandler->AddRunNumber(296304);
            alienHandler->AddRunNumber(296307);
            alienHandler->AddRunNumber(296309);
          }
          if(groupNumber == 14) {
            alienHandler->AddRunNumber(296375);
            alienHandler->AddRunNumber(296376);
            alienHandler->AddRunNumber(296377);
            alienHandler->AddRunNumber(296378);
            alienHandler->AddRunNumber(296379);
            alienHandler->AddRunNumber(296380);
            alienHandler->AddRunNumber(296381);
            alienHandler->AddRunNumber(296383);
            alienHandler->AddRunNumber(296414);
          }
          if(groupNumber == 15) {
            alienHandler->AddRunNumber(296415);
            alienHandler->AddRunNumber(296419);
            alienHandler->AddRunNumber(296420);
            alienHandler->AddRunNumber(296423);
            alienHandler->AddRunNumber(296424);
            alienHandler->AddRunNumber(296433);
            alienHandler->AddRunNumber(296472);
            alienHandler->AddRunNumber(296509);
            alienHandler->AddRunNumber(296510);
          }
          if(groupNumber == 16) {
            alienHandler->AddRunNumber(296511);
            alienHandler->AddRunNumber(296512);
            alienHandler->AddRunNumber(296514);
            alienHandler->AddRunNumber(296516);
            alienHandler->AddRunNumber(296547);
            alienHandler->AddRunNumber(296548);
            alienHandler->AddRunNumber(296549);
            alienHandler->AddRunNumber(296550);
            alienHandler->AddRunNumber(296551);
          }
          if(groupNumber == 17) {
            alienHandler->AddRunNumber(296552);
            alienHandler->AddRunNumber(296553);
            alienHandler->AddRunNumber(296594);
            alienHandler->AddRunNumber(296615);
            alienHandler->AddRunNumber(296616);
            alienHandler->AddRunNumber(296618);
            alienHandler->AddRunNumber(296619);
            alienHandler->AddRunNumber(296621);
            alienHandler->AddRunNumber(296622);
            alienHandler->AddRunNumber(296623);
          }
        }
        if (dataset.EqualTo("LHC18r")) {
          // select the input data
          alienHandler->SetGridDataDir("/alice/data/2018/LHC18r");
          alienHandler->SetDataPattern("pass3/AOD252/AOD/*/AliAOD.root");
          // MC has no prefix, data has prefix 000
          alienHandler->SetRunPrefix("000");

          //groupNumber = 0;
          if(groupNumber == 1) {
            // alienHandler->AddRunNumber(296690);
            // alienHandler->AddRunNumber(296691);
            // alienHandler->AddRunNumber(296693);
            // alienHandler->AddRunNumber(296694);
            // alienHandler->AddRunNumber(296752);
            // alienHandler->AddRunNumber(296784);
            // alienHandler->AddRunNumber(296785);
            alienHandler->AddRunNumber(296786);
            alienHandler->AddRunNumber(296787);
          }
          if(groupNumber == 2) {
            alienHandler->AddRunNumber(296790);
            alienHandler->AddRunNumber(296793);
            alienHandler->AddRunNumber(296794);
            alienHandler->AddRunNumber(296799);
            alienHandler->AddRunNumber(296835);
            alienHandler->AddRunNumber(296836);
            alienHandler->AddRunNumber(296838);
            alienHandler->AddRunNumber(296839);
            alienHandler->AddRunNumber(296848);
          }
          if(groupNumber == 3) {
            alienHandler->AddRunNumber(296850);
            alienHandler->AddRunNumber(296851);
            alienHandler->AddRunNumber(296852);
            alienHandler->AddRunNumber(296894);
            alienHandler->AddRunNumber(296899);
            alienHandler->AddRunNumber(296900);
            alienHandler->AddRunNumber(296903);
            alienHandler->AddRunNumber(296930);
            alienHandler->AddRunNumber(296931);
          }
          if(groupNumber == 4) {
            alienHandler->AddRunNumber(296932);
            alienHandler->AddRunNumber(296934);
            alienHandler->AddRunNumber(296935);
            alienHandler->AddRunNumber(296938);
            alienHandler->AddRunNumber(296941);
            alienHandler->AddRunNumber(296966);
            alienHandler->AddRunNumber(297031);
            alienHandler->AddRunNumber(297035);
            alienHandler->AddRunNumber(297085);
            alienHandler->AddRunNumber(297117);
          }
          if(groupNumber == 5) {
            alienHandler->AddRunNumber(297118);
            alienHandler->AddRunNumber(297119);
            alienHandler->AddRunNumber(297123);
            alienHandler->AddRunNumber(297124);
            alienHandler->AddRunNumber(297128);
            alienHandler->AddRunNumber(297129);
            alienHandler->AddRunNumber(297132);
            alienHandler->AddRunNumber(297133);
            alienHandler->AddRunNumber(297193);
            alienHandler->AddRunNumber(297195);
          }
          if(groupNumber == 6) {
            alienHandler->AddRunNumber(297196);
            alienHandler->AddRunNumber(297218);
            alienHandler->AddRunNumber(297221);
            alienHandler->AddRunNumber(297222);
            alienHandler->AddRunNumber(297278);
            alienHandler->AddRunNumber(297310);
            alienHandler->AddRunNumber(297311);
            alienHandler->AddRunNumber(297312);
            alienHandler->AddRunNumber(297317);
            alienHandler->AddRunNumber(297332);
          }
          if(groupNumber == 7) {
            alienHandler->AddRunNumber(297333);
            alienHandler->AddRunNumber(297335);
            alienHandler->AddRunNumber(297336);
            alienHandler->AddRunNumber(297363);
            alienHandler->AddRunNumber(297366);
            alienHandler->AddRunNumber(297367);
            alienHandler->AddRunNumber(297372);
            alienHandler->AddRunNumber(297379);
            alienHandler->AddRunNumber(297380);
            alienHandler->AddRunNumber(297405);
          }
          if(groupNumber == 8) {
            alienHandler->AddRunNumber(297406);
            alienHandler->AddRunNumber(297413);
            alienHandler->AddRunNumber(297414);
            alienHandler->AddRunNumber(297415);
            alienHandler->AddRunNumber(297441);
            alienHandler->AddRunNumber(297442);
            alienHandler->AddRunNumber(297446);
            alienHandler->AddRunNumber(297450);
            alienHandler->AddRunNumber(297451);
          }
          if(groupNumber == 9) {
            alienHandler->AddRunNumber(297452);
            alienHandler->AddRunNumber(297479);
            // alienHandler->AddRunNumber(297483);
            // alienHandler->AddRunNumber(297512);
            // alienHandler->AddRunNumber(297537);
            // alienHandler->AddRunNumber(297540);
            // alienHandler->AddRunNumber(297541);
            // alienHandler->AddRunNumber(297542);
          }
          if(groupNumber == 10) {
            alienHandler->AddRunNumber(297544);
            alienHandler->AddRunNumber(297558);
            alienHandler->AddRunNumber(297588);
            alienHandler->AddRunNumber(297590);
            alienHandler->AddRunNumber(297595);
            alienHandler->AddRunNumber(297624);
          }
        }

        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(40);
        alienHandler->SetExecutable("myTask.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(10000);
        alienHandler->SetJDLName("myTask.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate") 
        // (see below) mode, set SetMergeViaJDL(kFALSE) 
        // to collect final results
        alienHandler->SetMaxMergeStages(1);
        if(gridRun || gridMerge) {
          alienHandler->SetMergeViaJDL(kTRUE);
        }
        if (localMerge) {
          alienHandler->SetMergeViaJDL(kFALSE);
        }

        // define the output folders
        if (dataset.EqualTo("LHC10h")) alienHandler->SetGridWorkingDir("PIDAngCorr_4050/LHC10h");
        if (dataset.EqualTo("LHC15o")) alienHandler->SetGridWorkingDir("PIDAngCorr_4050/LHC15o");
        if (dataset.EqualTo("LHC18q")) alienHandler->SetGridWorkingDir("PIDAngCorr_4050/LHC18q");
        if (dataset.EqualTo("LHC18r")) alienHandler->SetGridWorkingDir("PIDAngCorr_4050/LHC18r");
        alienHandler->SetGridOutputDir("output");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            if(gridMerge || localMerge) alienHandler->SetRunMode("terminate");
            if(gridRun)   alienHandler->SetRunMode("full");
            mgr->StartAnalysis("grid");
        }
    }
}

