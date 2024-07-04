import ROOT
from root_utils import ROOT_AGC_Utils

root_utils = ROOT_AGC_Utils()

root_utils.set_xmin(110)
root_utils.set_rebinning(2)
root_utils.apply_rebinning("data/histograms.root", "data/temp_histos.root")

# root_utils.set_output_path("temp_histos.root")

meas = ROOT.RooStats.HistFactory.Measurement("meas", "meas")
meas.SetLumi(1.0)
meas.SetLumiRelErr(0.0)

input_file = "data/temp_histos.root"

channel = ROOT.RooStats.HistFactory.Channel("channel_4j1b_CR")
channel.SetData("4j1b_pseudodata", input_file)
channel.SetStatErrorConfig(0.001, "Gaussian")

ttbar = ROOT.RooStats.HistFactory.Sample("ttbar", "4j1b_ttbar", input_file)
ttbar.AddOverallSys("Lumi", 0.97, 1.03)
ttbar.ActivateStatError()
root_utils.AddNormPlusShapeHistoSys(ttbar, "ME_variation", "4j1b_ttbar_ME_var", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "PS_variation", "4j1b_ttbar_PS_var", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "tt_scale_var", "4j1b_ttbar_scaleup", input_file, "", 
                        "4j1b_ttbar_scaledown", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "jet_energy_scale", "4j1b_ttbar_pt_scale_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "jet_energy_res", "4j1b_ttbar_pt_res_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "b_tag_NP_1", "4j1b_ttbar_btag_var_0_up", input_file, "", 
                        "4j1b_ttbar_btag_var_0_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "b_tag_NP_2", "4j1b_ttbar_btag_var_1_up", input_file, "", 
                        "4j1b_ttbar_btag_var_1_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "b_tag_NP_3", "4j1b_ttbar_btag_var_2_up", input_file, "", 
                        "4j1b_ttbar_btag_var_2_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "b_tag_NP_4", "4j1b_ttbar_btag_var_3_up", input_file, "", 
                         "4j1b_ttbar_btag_var_3_down", input_file, "")
ttbar.AddNormFactor("ttbar_norm", 1, 0, 10)
channel.AddSample(ttbar)

wjets = ROOT.RooStats.HistFactory.Sample("wjets", "4j1b_wjets", input_file)
wjets.ActivateStatError()
wjets.AddOverallSys("Lumi", 0.97, 1.03)
root_utils.AddNormPlusShapeHistoSys(wjets, "jet_energy_scale", "4j1b_wjets_pt_scale_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(wjets, "jet_energy_res", "4j1b_wjets_pt_res_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(wjets, "b_tag_NP_1", "4j1b_wjets_btag_var_0_up", input_file, "", 
                        "4j1b_wjets_btag_var_0_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(wjets, "b_tag_NP_2", "4j1b_wjets_btag_var_1_up", input_file, "", 
                        "4j1b_wjets_btag_var_1_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(wjets, "b_tag_NP_3", "4j1b_wjets_btag_var_2_up", input_file, "", 
                        "4j1b_wjets_btag_var_2_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(wjets, "b_tag_NP_4", "4j1b_wjets_btag_var_3_up", input_file, "", 
                         "4j1b_wjets_btag_var_3_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(wjets, "w_plus_jets_scale_var", "4j1b_wjets_scale_var_up", input_file, "", 
                        "4j1b_wjets_scale_var_down", input_file, "")
channel.AddSample(wjets)

single_top_s_chan = ROOT.RooStats.HistFactory.Sample("single_top_s", "4j1b_single_top_s_chan", input_file)
single_top_s_chan.ActivateStatError()
single_top_s_chan.AddOverallSys("Lumi", 0.97, 1.03)
root_utils.AddNormPlusShapeHistoSys(single_top_s_chan, "jet_energy_scale", "4j1b_single_top_s_chan_pt_scale_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_s_chan, "jet_energy_res", "4j1b_single_top_s_chan_pt_res_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_s_chan, "b_tag_NP_1", "4j1b_single_top_s_chan_btag_var_0_up", input_file, "", 
                        "4j1b_single_top_s_chan_btag_var_0_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_s_chan, "b_tag_NP_2", "4j1b_single_top_s_chan_btag_var_1_up", input_file, "", 
                        "4j1b_single_top_s_chan_btag_var_1_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_s_chan, "b_tag_NP_3", "4j1b_single_top_s_chan_btag_var_2_up", input_file, "", 
                        "4j1b_single_top_s_chan_btag_var_2_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_s_chan, "b_tag_NP_4", "4j1b_single_top_s_chan_btag_var_3_up", input_file, "", 
                         "4j1b_single_top_s_chan_btag_var_3_down", input_file, "")

channel.AddSample(single_top_s_chan)

single_top_t_chan = ROOT.RooStats.HistFactory.Sample("single_top_t", "4j1b_single_top_t_chan", input_file)
single_top_t_chan.ActivateStatError()
single_top_t_chan.AddOverallSys("Lumi", 0.97, 1.03)
root_utils.AddNormPlusShapeHistoSys(single_top_t_chan, "jet_energy_scale", "4j1b_single_top_t_chan_pt_scale_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_t_chan, "jet_energy_res", "4j1b_single_top_t_chan_pt_res_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_t_chan, "b_tag_NP_1", "4j1b_single_top_t_chan_btag_var_0_up", input_file, "", 
                        "4j1b_single_top_t_chan_btag_var_0_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_t_chan, "b_tag_NP_2", "4j1b_single_top_t_chan_btag_var_1_up", input_file, "", 
                        "4j1b_single_top_t_chan_btag_var_1_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_t_chan, "b_tag_NP_3", "4j1b_single_top_t_chan_btag_var_2_up", input_file, "", 
                        "4j1b_single_top_t_chan_btag_var_2_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_t_chan, "b_tag_NP_4", "4j1b_single_top_t_chan_btag_var_3_up", input_file, "", 
                         "4j1b_single_top_t_chan_btag_var_3_down", input_file, "")

channel.AddSample(single_top_t_chan)

single_top_tW = ROOT.RooStats.HistFactory.Sample("single_top_tW", "4j1b_single_top_tW", input_file)
single_top_tW.ActivateStatError()
single_top_tW.AddOverallSys("Lumi", 0.97, 1.03)
root_utils.AddNormPlusShapeHistoSys(single_top_tW, "jet_energy_scale", "4j1b_single_top_tW_pt_scale_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_tW, "jet_energy_res", "4j1b_single_top_tW_pt_res_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_tW, "b_tag_NP_1", "4j1b_single_top_tW_btag_var_0_up", input_file, "", 
                        "4j1b_single_top_tW_btag_var_0_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_tW, "b_tag_NP_2", "4j1b_single_top_tW_btag_var_1_up", input_file, "", 
                        "4j1b_single_top_tW_btag_var_1_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_tW, "b_tag_NP_3", "4j1b_single_top_tW_btag_var_2_up", input_file, "", 
                        "4j1b_single_top_tW_btag_var_2_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_tW, "b_tag_NP_4", "4j1b_single_top_tW_btag_var_3_up", input_file, "", 
                         "4j1b_single_top_tW_btag_var_3_down", input_file, "")

channel.AddSample(single_top_tW)

meas.AddChannel(channel)

channel_2b = ROOT.RooStats.HistFactory.Channel("channel_4j2b_SR")
channel_2b.SetData("4j2b_pseudodata", input_file)
channel_2b.SetStatErrorConfig(0.001, "Gaussian")

ttbar = ROOT.RooStats.HistFactory.Sample("ttbar", "4j2b_ttbar", input_file)
ttbar.AddOverallSys("Lumi", 0.97, 1.03)
ttbar.ActivateStatError()
root_utils.AddNormPlusShapeHistoSys(ttbar, "ME_variation", "4j2b_ttbar_ME_var", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "PS_variation", "4j2b_ttbar_PS_var", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "tt_scale_var", "4j2b_ttbar_scaleup", input_file, "", 
                        "4j2b_ttbar_scaledown", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "jet_energy_scale", "4j2b_ttbar_pt_scale_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "jet_energy_res", "4j2b_ttbar_pt_res_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "b_tag_NP_1", "4j2b_ttbar_btag_var_0_up", input_file, "", 
                        "4j2b_ttbar_btag_var_0_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "b_tag_NP_2", "4j2b_ttbar_btag_var_1_up", input_file, "", 
                        "4j2b_ttbar_btag_var_1_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "b_tag_NP_3", "4j2b_ttbar_btag_var_2_up", input_file, "", 
                        "4j2b_ttbar_btag_var_2_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(ttbar, "b_tag_NP_4", "4j2b_ttbar_btag_var_3_up", input_file, "", 
                         "4j2b_ttbar_btag_var_3_down", input_file, "")
ttbar.AddNormFactor("ttbar_norm", 1, 0, 10)
channel_2b.AddSample(ttbar)

wjets = ROOT.RooStats.HistFactory.Sample("wjets", "4j2b_wjets", input_file)
wjets.ActivateStatError()
wjets.AddOverallSys("Lumi", 0.97, 1.03)
root_utils.AddNormPlusShapeHistoSys(wjets, "jet_energy_scale", "4j2b_wjets_pt_scale_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(wjets, "jet_energy_res", "4j2b_wjets_pt_res_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(wjets, "b_tag_NP_1", "4j2b_wjets_btag_var_0_up", input_file, "", 
                        "4j2b_wjets_btag_var_0_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(wjets, "b_tag_NP_2", "4j2b_wjets_btag_var_1_up", input_file, "", 
                        "4j2b_wjets_btag_var_1_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(wjets, "b_tag_NP_3", "4j2b_wjets_btag_var_2_up", input_file, "", 
                        "4j2b_wjets_btag_var_2_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(wjets, "b_tag_NP_4", "4j2b_wjets_btag_var_3_up", input_file, "", 
                         "4j2b_wjets_btag_var_3_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(wjets, "w_plus_jets_scale_var", "4j2b_wjets_scale_var_up", input_file, "", 
                        "4j2b_wjets_scale_var_down", input_file, "")
channel_2b.AddSample(wjets)

single_top_s_chan = ROOT.RooStats.HistFactory.Sample("single_top_s", "4j2b_single_top_s_chan", input_file)
single_top_s_chan.ActivateStatError()
single_top_s_chan.AddOverallSys("Lumi", 0.97, 1.03)
root_utils.AddNormPlusShapeHistoSys(single_top_s_chan, "jet_energy_scale", "4j2b_single_top_s_chan_pt_scale_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_s_chan, "jet_energy_res", "4j2b_single_top_s_chan_pt_res_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_s_chan, "b_tag_NP_1", "4j2b_single_top_s_chan_btag_var_0_up", input_file, "", 
                        "4j2b_single_top_s_chan_btag_var_0_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_s_chan, "b_tag_NP_2", "4j2b_single_top_s_chan_btag_var_1_up", input_file, "", 
                        "4j2b_single_top_s_chan_btag_var_1_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_s_chan, "b_tag_NP_3", "4j2b_single_top_s_chan_btag_var_2_up", input_file, "", 
                        "4j2b_single_top_s_chan_btag_var_2_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_s_chan, "b_tag_NP_4", "4j2b_single_top_s_chan_btag_var_3_up", input_file, "", 
                         "4j2b_single_top_s_chan_btag_var_3_down", input_file, "")
channel_2b.AddSample(single_top_s_chan)

single_top_t_chan = ROOT.RooStats.HistFactory.Sample("single_top_t", "4j2b_single_top_t_chan", input_file)
single_top_t_chan.ActivateStatError()
single_top_t_chan.AddOverallSys("Lumi", 0.97, 1.03)
root_utils.AddNormPlusShapeHistoSys(single_top_t_chan, "jet_energy_scale", "4j2b_single_top_t_chan_pt_scale_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_t_chan, "jet_energy_res", "4j2b_single_top_t_chan_pt_res_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_t_chan, "b_tag_NP_1", "4j2b_single_top_t_chan_btag_var_0_up", input_file, "", 
                        "4j2b_single_top_t_chan_btag_var_0_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_t_chan, "b_tag_NP_2", "4j2b_single_top_t_chan_btag_var_1_up", input_file, "", 
                        "4j2b_single_top_t_chan_btag_var_1_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_t_chan, "b_tag_NP_3", "4j2b_single_top_t_chan_btag_var_2_up", input_file, "", 
                        "4j2b_single_top_t_chan_btag_var_2_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_t_chan, "b_tag_NP_4", "4j2b_single_top_t_chan_btag_var_3_up", input_file, "", 
                         "4j2b_single_top_t_chan_btag_var_3_down", input_file, "")

channel_2b.AddSample(single_top_t_chan)

single_top_tW = ROOT.RooStats.HistFactory.Sample("single_top_tW", "4j2b_single_top_tW", input_file)
single_top_tW.ActivateStatError()
single_top_tW.AddOverallSys("Lumi", 0.97, 1.03)
root_utils.AddNormPlusShapeHistoSys(single_top_tW, "jet_energy_scale", "4j2b_single_top_tW_pt_scale_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_tW, "jet_energy_res", "4j2b_single_top_tW_pt_res_up", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_tW, "b_tag_NP_1", "4j2b_single_top_tW_btag_var_0_up", input_file, "", 
                        "4j2b_single_top_tW_btag_var_0_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_tW, "b_tag_NP_2", "4j2b_single_top_tW_btag_var_1_up", input_file, "", 
                        "4j2b_single_top_tW_btag_var_1_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_tW, "b_tag_NP_3", "4j2b_single_top_tW_btag_var_2_up", input_file, "", 
                        "4j2b_single_top_tW_btag_var_2_down", input_file, "")
root_utils.AddNormPlusShapeHistoSys(single_top_tW, "b_tag_NP_4", "4j2b_single_top_tW_btag_var_3_up", input_file, "", 
                         "4j2b_single_top_tW_btag_var_3_down", input_file, "")

channel_2b.AddSample(single_top_tW)
meas.AddChannel(channel_2b)

meas.SetPOI("ttbar_norm")
meas.CollectHistograms()

ws = ROOT.RooStats.HistFactory.MakeModelAndMeasurementFast(meas)

# Retrieve the ModelConfig
modelConfig = ws.obj("ModelConfig")

json_tool = ROOT.RooJSONFactoryWSTool(ws)
json_tool.exportJSON("my_json.json")

# Extract the PDF and global observables
pdf = modelConfig.GetPdf()
globalObservables = ROOT.RooArgSet(modelConfig.GetGlobalObservables())

meas.PrintXML("xmlFromPy", meas.GetOutputFilePrefix())

# Perform the fit
result = pdf.fitTo(ws.data("obsData"), ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(ROOT.RooFit.FATAL), ROOT.RooFit.GlobalObservables(globalObservables))

# Print the fit result
result.Print()