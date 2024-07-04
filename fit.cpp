#include "TTree.h"
#include "TROOT.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include <RooFitHS3/RooJSONFactoryWSTool.h>
#include "RooStats/HistFactory/HistoToWorkspaceFactoryFast.h"
#include "RooStats/ModelConfig.h"
#include "TFile.h"
#include "utils.h"


using namespace RooStats;
using namespace HistFactory;

void histo()
{

    utils utils;

    std::string InputFile = "./histograms.root";
    RooStats::HistFactory::Measurement meas("meas", "meas");
    RooStats::HistFactory::Channel chan( "channel_4j1b_CR" );
    // chan.Set
    chan.SetData( "4j1b_pseudodata", InputFile);
    chan.SetStatErrorConfig( 0.01, "Gaussian" );

    meas.SetLumi( 1.0 );
    meas.SetLumiRelErr( 0.0 );

    Sample tt_bar = Sample("ttbar", "4j1b_ttbar", InputFile);
    tt_bar.AddOverallSys("Lumi", 0.97, 1.03);
    tt_bar.ActivateStatError();
    tt_bar.AddHistoSys("ME_variation", "4j1b_ttbar_ME_var_DOWN", InputFile, "", 
                            "4j1b_ttbar_ME_var", InputFile, "",
                            true, true); // ALL HIDDEN AUTOMATIZATION IS HERE :)
    tt_bar.AddHistoSys("PS_variation", "4j1b_ttbar_PS_var_DOWN", InputFile, "", 
                        "4j1b_ttbar_PS_var", InputFile, "",
                        true, true);
    tt_bar.AddHistoSys("tt_scale_var", "4j1b_ttbar_scaledown", InputFile, "", 
                        "4j1b_ttbar_scaleup", InputFile, "",
                        false, true);
    tt_bar.AddHistoSys("jet_energy_scale", "4j1b_ttbar_pt_scale_down", InputFile, "", 
                        "4j1b_ttbar_pt_scale_up", InputFile, "",
                        true, true);
    tt_bar.AddHistoSys("jet_energy_res", "4j1b_ttbar_pt_res_down", InputFile, "", 
                        "4j1b_ttbar_pt_res_up", InputFile, "",
                        true, true); 
    tt_bar.AddHistoSys("b_tag_NP_1", "4j1b_ttbar_btag_var_0_down", InputFile, "", 
                        "4j1b_ttbar_btag_var_0_up", InputFile, "",
                        false, true); 
    tt_bar.AddHistoSys("b_tag_NP_2", "4j1b_ttbar_btag_var_1_down", InputFile, "", 
                        "4j1b_ttbar_btag_var_1_up", InputFile, "",
                        false, true); 
    tt_bar.AddHistoSys("b_tag_NP_3", "4j1b_ttbar_btag_var_2_down", InputFile, "", 
                        "4j1b_ttbar_btag_var_2_up", InputFile, "",
                        false, true); 
    tt_bar.AddHistoSys("b_tag_NP_4", "4j1b_ttbar_btag_var_3_down", InputFile, "", 
                        "4j1b_ttbar_btag_var_3_up", InputFile, "",
                        false, true); 
    tt_bar.AddNormFactor("ttbar_norm", 1, 0, 10);
    chan.AddSample(tt_bar);

    Sample wjets = Sample("wjets", "4j1b_wjets", InputFile);
    wjets.AddOverallSys("Lumi", 0.97, 1.03);
    wjets.ActivateStatError();
    wjets.AddHistoSys("jet_energy_scale", "4j1b_wjets_pt_scale_down", InputFile, "", 
                        "4j1b_wjets_pt_scale_up", InputFile, "",
                        true, true);
    wjets.AddHistoSys("jet_energy_res", "4j1b_wjets_pt_res_down", InputFile, "", 
                        "4j1b_wjets_pt_res_up", InputFile, "",
                        true, true); 
    wjets.AddHistoSys("b_tag_NP_1", "4j1b_wjets_btag_var_0_down", InputFile, "", 
                        "4j1b_wjets_btag_var_0_up", InputFile, "",
                        false, true); 
    wjets.AddHistoSys("b_tag_NP_2", "4j1b_wjets_btag_var_1_down", InputFile, "", 
                        "4j1b_wjets_btag_var_1_up", InputFile, "",
                        false, true); 
    wjets.AddHistoSys("b_tag_NP_3", "4j1b_wjets_btag_var_2_down", InputFile, "", 
                        "4j1b_wjets_btag_var_2_up", InputFile, "",
                        false, true); 
    wjets.AddHistoSys("b_tag_NP_4", "4j1b_wjets_btag_var_3_down", InputFile, "", 
                        "4j1b_wjets_btag_var_3_up", InputFile, "",
                        false, true); 
    wjets.AddHistoSys("w_plus_jets_scale_var", "4j1b_wjets_scale_var_down", InputFile, "", 
                        "4j1b_wjets_scale_var_up", InputFile, "",
                        false, true);                
    chan.AddSample(wjets);

    Sample single_top_s_chan = Sample("single_top_s", "4j1b_single_top_s_chan", InputFile);
    single_top_s_chan.AddOverallSys("Lumi", 0.97, 1.03);
    single_top_s_chan.ActivateStatError();
    single_top_s_chan.AddHistoSys("jet_energy_scale", "4j1b_single_top_s_chan_pt_scale_down", InputFile, "", 
                        "4j1b_single_top_s_chan_pt_scale_up", InputFile, "",
                        true, true);
    single_top_s_chan.AddHistoSys("jet_energy_res", "4j1b_single_top_s_chan_pt_res_down", InputFile, "", 
                        "4j1b_single_top_s_chan_pt_res_up", InputFile, "",
                        true, true); 
    single_top_s_chan.AddHistoSys("b_tag_NP_1", "4j1b_single_top_s_chan_btag_var_0_down", InputFile, "", 
                        "4j1b_single_top_s_chan_btag_var_0_up", InputFile, "",
                        false, true); 
    single_top_s_chan.AddHistoSys("b_tag_NP_2", "4j1b_single_top_s_chan_btag_var_1_down", InputFile, "", 
                        "4j1b_single_top_s_chan_btag_var_1_up", InputFile, "",
                        false, true); 
    single_top_s_chan.AddHistoSys("b_tag_NP_3", "4j1b_single_top_s_chan_btag_var_2_down", InputFile, "", 
                        "4j1b_single_top_s_chan_btag_var_2_up", InputFile, "",
                        false, true); 
    single_top_s_chan.AddHistoSys("b_tag_NP_4", "4j1b_single_top_s_chan_btag_var_3_down", InputFile, "", 
                        "4j1b_single_top_s_chan_btag_var_3_up", InputFile, "",
                        false, true); 
    chan.AddSample(single_top_s_chan);

    Sample single_top_t_chan = Sample("single_top_t", "4j1b_single_top_t_chan", InputFile);
    single_top_t_chan.AddOverallSys("Lumi", 0.97, 1.03);
    single_top_t_chan.ActivateStatError();

    single_top_t_chan.AddHistoSys("jet_energy_scale", "4j1b_single_top_t_chan_pt_scale_down", InputFile, "", 
                        "4j1b_single_top_t_chan_pt_scale_up", InputFile, "",
                        true, true);
    single_top_t_chan.AddHistoSys("jet_energy_res", "4j1b_single_top_t_chan_pt_res_down", InputFile, "", 
                        "4j1b_single_top_t_chan_pt_res_up", InputFile, "",
                        true, true); 
    single_top_t_chan.AddHistoSys("b_tag_NP_1", "4j1b_single_top_t_chan_btag_var_0_down", InputFile, "", 
                        "4j1b_single_top_t_chan_btag_var_0_up", InputFile, "",
                        false, true); 
    single_top_t_chan.AddHistoSys("b_tag_NP_2", "4j1b_single_top_t_chan_btag_var_1_down", InputFile, "", 
                        "4j1b_single_top_t_chan_btag_var_1_up", InputFile, "",
                        false, true); 
    single_top_t_chan.AddHistoSys("b_tag_NP_3", "4j1b_single_top_t_chan_btag_var_2_down", InputFile, "", 
                        "4j1b_single_top_t_chan_btag_var_2_up", InputFile, "",
                        false, true); 
    single_top_t_chan.AddHistoSys("b_tag_NP_4", "4j1b_single_top_t_chan_btag_var_3_down", InputFile, "", 
                        "4j1b_single_top_t_chan_btag_var_3_up", InputFile, "",
                        false, true); 
    chan.AddSample(single_top_t_chan);

    Sample single_top_tW = Sample("single_top_tW", "4j1b_single_top_tW", InputFile);
    single_top_tW.AddOverallSys("Lumi", 0.97, 1.03);
    single_top_tW.ActivateStatError();
    single_top_tW.AddHistoSys("jet_energy_scale", "4j1b_single_top_tW_pt_scale_down", InputFile, "", 
                        "4j1b_single_top_tW_pt_scale_up", InputFile, "",
                        true, true);
    single_top_tW.AddHistoSys("jet_energy_res", "4j1b_single_top_tW_pt_res_down", InputFile, "", 
                        "4j1b_single_top_tW_pt_res_up", InputFile, "",
                        true, true); 
    single_top_tW.AddHistoSys("b_tag_NP_1", "4j1b_single_top_tW_btag_var_0_down", InputFile, "", 
                        "4j1b_single_top_tW_btag_var_0_up", InputFile, "",
                        false, true); 
    single_top_tW.AddHistoSys("b_tag_NP_2", "4j1b_single_top_tW_btag_var_1_down", InputFile, "", 
                        "4j1b_single_top_tW_btag_var_1_up", InputFile, "",
                        false, true); 
    single_top_tW.AddHistoSys("b_tag_NP_3", "4j1b_single_top_tW_btag_var_2_down", InputFile, "", 
                        "4j1b_single_top_tW_btag_var_2_up", InputFile, "",
                        false, true); 
    single_top_tW.AddHistoSys("b_tag_NP_4", "4j1b_single_top_tW_btag_var_3_down", InputFile, "", 
                        "4j1b_single_top_tW_btag_var_3_up", InputFile, "",
                        false, true); 
    chan.AddSample(single_top_tW);


    meas.AddChannel(chan);



    RooStats::HistFactory::Channel chan_2b( "channel_4j2b_SR" );
    chan_2b.SetData( "4j2b_pseudodata", InputFile);
    chan_2b.SetStatErrorConfig( 0.01, "Gaussian" );

    tt_bar = Sample("ttbar", "4j2b_ttbar", InputFile);
    tt_bar.ActivateStatError();
    tt_bar.AddOverallSys("Lumi", 0.97, 1.03);
    tt_bar.AddHistoSys("ME_variation", "4j2b_ttbar_ME_var_DOWN", InputFile, "", 
                            "4j2b_ttbar_ME_var", InputFile, "",
                            true, true); // ALL HIDDEN AUTOMATIZATION IS HERE :)
    tt_bar.AddHistoSys("PS_variation", "4j2b_ttbar_PS_var_DOWN", InputFile, "", 
                        "4j2b_ttbar_PS_var", InputFile, "",
                        true, true);
    tt_bar.AddHistoSys("tt_scale_var", "4j2b_ttbar_scaledown", InputFile, "", 
                        "4j2b_ttbar_scaleup", InputFile, "",
                        false, true);
    tt_bar.AddHistoSys("jet_energy_scale", "4j2b_ttbar_pt_scale_down", InputFile, "", 
                        "4j2b_ttbar_pt_scale_up", InputFile, "",
                        true, true);
    tt_bar.AddHistoSys("jet_energy_res", "4j2b_ttbar_pt_res_down", InputFile, "", 
                        "4j2b_ttbar_pt_res_up", InputFile, "",
                        true, true); 
    tt_bar.AddHistoSys("b_tag_NP_1", "4j2b_ttbar_btag_var_0_down", InputFile, "", 
                        "4j2b_ttbar_btag_var_0_up", InputFile, "",
                        false, true); 
    tt_bar.AddHistoSys("b_tag_NP_2", "4j2b_ttbar_btag_var_1_down", InputFile, "", 
                        "4j2b_ttbar_btag_var_1_up", InputFile, "",
                        false, true); 
    tt_bar.AddHistoSys("b_tag_NP_3", "4j2b_ttbar_btag_var_2_down", InputFile, "", 
                        "4j2b_ttbar_btag_var_2_up", InputFile, "",
                        false, true);
    tt_bar.AddHistoSys("b_tag_NP_4", "4j2b_ttbar_btag_var_3_down", InputFile, "", 
                        "4j2b_ttbar_btag_var_3_up", InputFile, "",
                        false, true); 
    tt_bar.AddNormFactor("ttbar_norm", 1, 0, 10);
    chan_2b.AddSample(tt_bar);

    wjets = Sample("wjets", "4j2b_wjets", InputFile);
    wjets.ActivateStatError();
    wjets.AddOverallSys("Lumi", 0.97, 1.03);
    wjets.AddHistoSys("jet_energy_scale", "4j2b_wjets_pt_scale_down", InputFile, "", 
                        "4j2b_wjets_pt_scale_up", InputFile, "",
                        true, true);
    wjets.AddHistoSys("jet_energy_res", "4j2b_wjets_pt_res_down", InputFile, "", 
                        "4j2b_wjets_pt_res_up", InputFile, "",
                        true, true); 
    wjets.AddHistoSys("b_tag_NP_1", "4j2b_wjets_btag_var_0_down", InputFile, "", 
                        "4j2b_wjets_btag_var_0_up", InputFile, "",
                        false, true); 
    wjets.AddHistoSys("b_tag_NP_2", "4j2b_wjets_btag_var_1_down", InputFile, "", 
                        "4j2b_wjets_btag_var_1_up", InputFile, "",
                        false, true); 
    wjets.AddHistoSys("b_tag_NP_3", "4j2b_wjets_btag_var_2_down", InputFile, "", 
                        "4j2b_wjets_btag_var_2_up", InputFile, "",
                        false, true); 
    wjets.AddHistoSys("b_tag_NP_4", "4j2b_wjets_btag_var_3_down", InputFile, "", 
                        "4j2b_wjets_btag_var_3_up", InputFile, "",
                        false, true);  
    wjets.AddHistoSys("w_plus_jets_scale_var", "4j2b_wjets_scale_var_down", InputFile, "", 
                        "4j2b_wjets_scale_var_up", InputFile, "",
                        false, true);               
    chan_2b.AddSample(wjets);

    single_top_s_chan = Sample("single_top_s", "4j2b_single_top_s_chan", InputFile);
    single_top_s_chan.ActivateStatError();
    single_top_s_chan.AddOverallSys("Lumi", 0.97, 1.03);
    single_top_s_chan.AddHistoSys("jet_energy_scale", "4j2b_single_top_s_chan_pt_scale_down", InputFile, "", 
                        "4j2b_single_top_s_chan_pt_scale_up", InputFile, "",
                        true, true);
    single_top_s_chan.AddHistoSys("jet_energy_res", "4j2b_single_top_s_chan_pt_res_down", InputFile, "", 
                        "4j2b_single_top_s_chan_pt_res_up", InputFile, "",
                        true, true); 
    single_top_s_chan.AddHistoSys("b_tag_NP_1", "4j2b_single_top_s_chan_btag_var_0_down", InputFile, "", 
                        "4j2b_single_top_s_chan_btag_var_0_up", InputFile, "",
                        false, true); 
    single_top_s_chan.AddHistoSys("b_tag_NP_2", "4j2b_single_top_s_chan_btag_var_1_down", InputFile, "", 
                        "4j2b_single_top_s_chan_btag_var_1_up", InputFile, "",
                        false, true); 
    single_top_s_chan.AddHistoSys("b_tag_NP_3", "4j2b_single_top_s_chan_btag_var_2_down", InputFile, "", 
                        "4j2b_single_top_s_chan_btag_var_2_up", InputFile, "",
                        false, true); 
    single_top_s_chan.AddHistoSys("b_tag_NP_4", "4j2b_single_top_s_chan_btag_var_3_down", InputFile, "", 
                        "4j2b_single_top_s_chan_btag_var_3_up", InputFile, "",
                        false, true); 
    chan_2b.AddSample(single_top_s_chan);

    single_top_t_chan = Sample("single_top_t", "4j2b_single_top_t_chan", InputFile);
    single_top_t_chan.ActivateStatError();
    single_top_t_chan.AddOverallSys("Lumi", 0.97, 1.03);
    single_top_t_chan.AddHistoSys("jet_energy_scale", "4j2b_single_top_t_chan_pt_scale_down", InputFile, "", 
                        "4j2b_single_top_t_chan_pt_scale_up", InputFile, "",
                        true, true);
    single_top_t_chan.AddHistoSys("jet_energy_res", "4j2b_single_top_t_chan_pt_res_down", InputFile, "", 
                        "4j2b_single_top_t_chan_pt_res_up", InputFile, "",
                        true, true); 
    single_top_t_chan.AddHistoSys("b_tag_NP_1", "4j2b_single_top_t_chan_btag_var_0_down", InputFile, "", 
                        "4j2b_single_top_t_chan_btag_var_0_up", InputFile, "",
                        false, true); 
    single_top_t_chan.AddHistoSys("b_tag_NP_2", "4j2b_single_top_t_chan_btag_var_1_down", InputFile, "", 
                        "4j2b_single_top_t_chan_btag_var_1_up", InputFile, "",
                        false, true); 
    single_top_t_chan.AddHistoSys("b_tag_NP_3", "4j2b_single_top_t_chan_btag_var_2_down", InputFile, "", 
                        "4j2b_single_top_t_chan_btag_var_2_up", InputFile, "",
                        false, true); 
    single_top_t_chan.AddHistoSys("b_tag_NP_4", "4j2b_single_top_t_chan_btag_var_3_down", InputFile, "", 
                        "4j2b_single_top_t_chan_btag_var_3_up", InputFile, "",
                        false, true);
    chan_2b.AddSample(single_top_t_chan);

    single_top_tW = Sample("single_top_tW", "4j2b_single_top_tW", InputFile);
    single_top_tW.ActivateStatError();
    single_top_tW.AddOverallSys("Lumi", 0.97, 1.03);
    single_top_tW.AddHistoSys("jet_energy_scale", "4j2b_single_top_tW_pt_scale_down", InputFile, "", 
                        "4j2b_single_top_tW_pt_scale_up", InputFile, "",
                        true, true);
    single_top_tW.AddHistoSys("jet_energy_res", "4j2b_single_top_tW_pt_res_down", InputFile, "", 
                        "4j2b_single_top_tW_pt_res_up", InputFile, "",
                        true, true); 
    single_top_tW.AddHistoSys("b_tag_NP_1", "4j2b_single_top_tW_btag_var_0_down", InputFile, "", 
                        "4j2b_single_top_tW_btag_var_0_up", InputFile, "",
                        false, true); 
    single_top_tW.AddHistoSys("b_tag_NP_2", "4j2b_single_top_tW_btag_var_1_down", InputFile, "", 
                        "4j2b_single_top_tW_btag_var_1_up", InputFile, "",
                        false, true); 
    single_top_tW.AddHistoSys("b_tag_NP_3", "4j2b_single_top_tW_btag_var_2_down", InputFile, "", 
                        "4j2b_single_top_tW_btag_var_2_up", InputFile, "",
                        false, true); 
    single_top_tW.AddHistoSys("b_tag_NP_4", "4j2b_single_top_tW_btag_var_3_down", InputFile, "", 
                        "4j2b_single_top_tW_btag_var_3_up", InputFile, "",
                        false, true); 
    chan_2b.AddSample(single_top_tW);


    meas.AddChannel(chan_2b);


    

    meas.SetPOI("ttbar_norm");
    
    // meas.SetRebinningConfig(2, 110, 550);
    // meas.ApplyRebinning();

    meas.CollectHistograms();
    meas.PrintXML( "xmlFromPy", meas.GetOutputFilePrefix() );
    meas.PrintJSON( "jsonFromPy", meas.GetOutputFilePrefix() );



    // RooStats::HistFactory::Sample signal( "Signal", "Nominal", InputFile, "SR/Signal" );
    // signal.AddOverallSys( "Lumi",  0.95, 1.05 );
    // signal.ActivateStatError();
    // signal.AddNormFactor( "Signal_norm", 1, 0, 3 );
    // chan.AddSample( signal );

    MakeModelAndMeasurementFast( meas );

    std::unique_ptr<RooWorkspace> ws{MakeModelAndMeasurementFast(meas)};

    RooJSONFactoryWSTool{*ws}.exportJSON("my_json.json");
    RooStats::ModelConfig *modelConfig = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));

    RooAbsPdf *pdf = modelConfig->GetPdf();
    RooArgSet globalObservables{*modelConfig->GetGlobalObservables()};

    using namespace RooFit;
    std::unique_ptr<RooFitResult> result{
        pdf->fitTo(*ws->data("obsData"), Save(), PrintLevel(-1), GlobalObservables(globalObservables))
    };

    // result->Na

    result->Print();


    utils.CreateAndSavePicture("picture.png", *ws);




// Now, do the measurement
//   std::unique_ptr<RooWorkspace> ws{MakeModelAndMeasurementFast(meas)};


//   RooStats::ModelConfig *modelConfig = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));

//   RooAbsPdf *pdf = modelConfig->GetPdf();
//   RooArgSet globalObservables{*modelConfig->GetGlobalObservables()};

//   using namespace RooFit;
//   std::unique_ptr<RooFitResult> result{
//       pdf->fitTo(*ws->data("obsData"), Save(), PrintLevel(-1), GlobalObservables(globalObservables))
//   };

//   result->Print();

//   std::cout << ws->var("SigXsecOverSM")->getVal() << std::endl;
//   std::cout << ws->var("SigXsecOverSM")->getError() << std::endl;

}

int fit()
{
    histo();
    return 0;
}