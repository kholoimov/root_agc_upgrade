import ROOT
import math

class ROOT_AGC_Utils:
    def __init__(self):
        self.fRebin = 1
        self.fHistoBinLow = None
        self.fHistoBinHigh = None
        self.output_path = 'ROOT_AGC_Utils/HistFactoryExtra.root'
        self.fRebinning = False
        file = ROOT.TFile(self.output_path, "RECREATE")
        file.Close()

    def set_output_path(self, output_path):
        self.output_path = output_path
        file = ROOT.TFile(self.output_path, "RECREATE")
        file.Close()

    def set_rebinning(self, rebinning):
        self.fRebin = rebinning
        self.fRebinning = True

    def set_xmin(self, xmin):
        self.fHistoBinLow = xmin
    
    def set_xmax(self, xmax):
        self.fHistoBinHigh = xmax

    def is_integer(self, value):
        return int(value) == value
    
    def list_histograms(self, directory, path=""):
    # Get the list of keys in the directory
        keys = directory.GetListOfKeys()

        # Loop over all keys
        for key in keys:
            # Get the object associated with the key
            obj = key.ReadObj()
            
            # Construct the full path of the object
            if path == "":
                obj_path = obj.GetName()
            else: 
                obj_path = path + "/" + obj.GetName()

            # Check if the object is a directory
            if obj.InheritsFrom("TDirectory"):
                # Recursively list histograms in subdirectories
                self.list_histograms(obj, obj_path)
            # Check if the object is a histogram
            elif obj.InheritsFrom("TH1"):
                self.all_histograms.append(obj_path)

    def mkdir_p(self, file, directory_path):
        file.cd()
        parts = directory_path.strip("/").split("/")
        current_dir = ROOT.gDirectory
        for part in parts:
            if not current_dir.GetDirectory(part):
                current_dir.mkdir(part)
            current_dir.cd(part)
        ROOT.gDirectory = current_dir

    def apply_rebinning(self, filename, output_file_path):
        file = ROOT.TFile(filename, "READ")

        output_file = ROOT.TFile(output_file_path, "RECREATE")
        output_file.Close()

        self.all_histograms = []
        self.list_histograms(file)

        for hist_path in self.all_histograms:
            hist = file.Get(hist_path)
            hist_rebinned = self.rebin_histogram(hist)
            output_file = ROOT.TFile(output_file_path, "UPDATE")

            if "/" in hist_path:
                dir_path = hist_path.rsplit("/", 1)[0]
                self.mkdir_p(output_file, dir_path)
                output_file.cd(dir_path)

            hist_rebinned.Write(hist_path.split("/")[-1])
            output_file.Close()

        file.Close()

    def rebin_histogram(self, original):
        if not self.fRebinning:
            return original

        left_edge = self.fHistoBinLow
        right_edge = self.fHistoBinHigh

        if left_edge is None:
            left_edge = original.GetXaxis().GetXmin()

        if right_edge is None:
            right_edge = original.GetXaxis().GetXmax()

        original_bin_width = original.GetXaxis().GetBinWidth(1)
        assert_check_left = (left_edge - original.GetXaxis().GetXmin()) / original_bin_width
        assert_check_right = (original.GetXaxis().GetXmax() - right_edge) / original_bin_width

        # print(f"fHistoBinLow: {left_edge}")
        # print(f"fHistoBinHigh: {right_edge}")
        # print(f"original_bin_width: {original_bin_width}")

        if not self.is_integer(assert_check_left) or not self.is_integer(assert_check_right):
            print("Error: The left_edge and right_edge are not multiples of the original bin width")
            return original

        number_of_remove_bins = ((left_edge - original.GetXaxis().GetXmin()) + 
                                 (original.GetXaxis().GetXmax() - right_edge)) / original_bin_width
        # print(f"number of remove bins: {number_of_remove_bins}")
        new_nbins = int(original.GetNbinsX() - number_of_remove_bins)
        # print(f"new_nbins: {new_nbins // self.fRebin}")
        
        original_new = ROOT.TH1F(original.GetName() + "temp_rebin_clone", original.GetTitle(), new_nbins, left_edge, right_edge)
        skipped_bins_left = int((left_edge - original.GetXaxis().GetXmin()) / original_bin_width)
        
        for i in range(1, new_nbins + 1):
            bin_idx = i + skipped_bins_left
            original_new.SetBinContent(i, original.GetBinContent(bin_idx))
            original_new.SetBinError(i, original.GetBinError(bin_idx))
        
        output = original_new.Rebin(self.fRebin)
        output.SetDirectory(ROOT.nullptr)

        return output
    
    def AddNormPlusShapeHistoSys(self, sample, name, histoname_up, histofile_up, histopath_up,
                                 histoname_down = None, histofile_down = None, histopath_down = None):
        if histoname_down is None:
            self.Symmetrize_AddNormPlusShapeHistoSys(name, sample, histoname_up, histofile_up, histopath_up)
        else:
            self.NonSymmetrize_AddNormPlusShapeHistoSys(sample, name, histoname_up, histofile_up, histopath_up,
                                 histoname_down, histofile_down, histopath_down)
    
    def Symmetrize_AddNormPlusShapeHistoSys(self, name, sample, histoname, histofile, histopath):

        channel_name = str(histoname.split("_")[0])

        file = ROOT.TFile(histofile, "READ")
        dir = file.GetDirectory(histopath)
        hist_top = (dir.Get(histoname))

        hist_nominal_file = ROOT.TFile(sample.GetInputFile(), "READ")
        hist_nominal_name = sample.GetHistoName()
        hist_nominal_directory = hist_nominal_file.GetDirectory(sample.GetHistoPath())
        hist_nominal = (hist_nominal_directory.Get(hist_nominal_name))

        norm_factor_up = hist_top.Integral() / hist_nominal.Integral()
        h_new = hist_top.Clone(channel_name + "_" + str(sample.GetName()) + "_" + name  + "_norm_plus_shape_up_clone")
        h_new.Scale(1/norm_factor_up)

        h_down = hist_nominal.Clone(channel_name + "_" + str(sample.GetName()) + "_" + name  + "_norm_plus_shape_down_clone")
        h_down.Scale(2)
        h_down.Add(h_new, -1)

        output_file = ROOT.TFile(self.output_path, "UPDATE")

        hist_up_name  = str(channel_name + "_" + str(sample.GetName()) + "_" + name  + "_norm_plus_shape_up")
        hist_down_name = str(channel_name + "_" + str(sample.GetName()) + "_" + name  + "_norm_plus_shape_down")

        print(h_down)
        print(h_new)

        h_new.Write(hist_up_name)
        h_down.Write(hist_down_name)

        output_file.Close()

        histSys = ROOT.RooStats.HistFactory.HistoSys()

        histSys.SetName(name)
        # histSys.SetHistoName(sample.GetSampleName())
        histSys.SetHistoPathHigh("")
        histSys.SetHistoPathLow("")

        histSys.SetHistoNameHigh(hist_up_name)
        histSys.SetHistoNameLow(hist_down_name)

        histSys.SetInputFileLow(self.output_path)
        histSys.SetInputFileHigh(self.output_path)

        overallSys = ROOT.RooStats.HistFactory.OverallSys()
        overallSys.SetName(name)
        overallSys.SetLow(2 - norm_factor_up)
        overallSys.SetHigh(norm_factor_up)

        sample.GetHistoSysList().push_back(histSys)
        sample.GetOverallSysList().push_back(overallSys)

    def NonSymmetrize_AddNormPlusShapeHistoSys(self, sample, name, histoname_up, histofile_up, histopath_up,
                                 histoname_down, histofile_down, histopath_down):
        channel_name = str(histoname_up.split("_")[0])

        file = ROOT.TFile(histofile_up, "READ")
        dir = file.GetDirectory(histopath_up)
        hist_top = (dir.Get(histoname_up))
        
        hist_nominal_file = ROOT.TFile(sample.GetInputFile(), "READ")
        hist_nominal_name = sample.GetHistoName()
        hist_nominal_directory = hist_nominal_file.GetDirectory(sample.GetHistoPath())
        hist_nominal = (hist_nominal_directory.Get(hist_nominal_name))

        norm_factor_up = hist_top.Integral() / hist_nominal.Integral()
        h_new_up = hist_top.Clone(channel_name + "_" + str(sample.GetName()) + "_" + name  + "_norm_plus_shape_up_clone")
        h_new_up.Scale(1/norm_factor_up)

        file_down = ROOT.TFile(histofile_down, "READ")
        dir_down = file_down.GetDirectory(histopath_down)
        hist_down = (dir_down.Get(histoname_down))

        norm_factor_down = hist_down.Integral() / hist_nominal.Integral()
        h_new_down = hist_down.Clone(channel_name + "_" + str(sample.GetName()) + "_" + name  + "_norm_plus_shape_down_clone")
        h_new_down.Scale(1/norm_factor_down)

        output_file = ROOT.TFile(self.output_path, "UPDATE")

        hist_up_name  = str(channel_name + "_" + str(sample.GetName()) + "_" + name  + "_norm_plus_shape_up")
        hist_down_name = str(channel_name + "_" + str(sample.GetName()) + "_" + name  + "_norm_plus_shape_down")

        print(h_new_down)
        print(h_new_up)

        h_new_up.Write(hist_up_name)
        h_new_down.Write(hist_down_name)

        output_file.Close()

        histSys = ROOT.RooStats.HistFactory.HistoSys()

        histSys.SetName(name)
        # histSys.SetHistoName(sample.GetSampleName())
        histSys.SetHistoPathHigh("")
        histSys.SetHistoPathLow("")

        histSys.SetHistoNameHigh(hist_up_name)
        histSys.SetHistoNameLow(hist_down_name)

        histSys.SetInputFileLow(self.output_path)
        histSys.SetInputFileHigh(self.output_path)

        overallSys = ROOT.RooStats.HistFactory.OverallSys()
        overallSys.SetName(name)
        overallSys.SetLow(norm_factor_down)
        overallSys.SetHigh(norm_factor_up)

        sample.GetHistoSysList().push_back(histSys)
        sample.GetOverallSysList().push_back(overallSys)


        # hist_nominal = sample.Get


    

#     cv_array = []
# hs_stacks = []
# sample_histograms = []
# bias_graphs = []
# bias_second_graphs = []
# normal_lines = []
# second_histos = []
# second_hs_stacks = []
# afterfit_histograms = []
# error_graphs = []
# error_pre_graphs = []
 
# predefined_colors = [
#     ROOT.TColor.GetColor("#F6C3C0"),  # Light Red
#     ROOT.TColor.GetColor("#A9D0F5"),  # Light Blue
#     ROOT.TColor.GetColor("#B6D7A8"),  # Light Green
#     ROOT.TColor.GetColor("#EAD1DC"),  # Light Magenta
#     ROOT.TColor.GetColor("#D0E9E8"),  # Light Cyan
#     ROOT.TColor.GetColor("#F6F9C6"),  # Light Yellow
#     ROOT.TColor.GetColor("#F6B8A3"),  # Light Orange
#     ROOT.TColor.GetColor("#F9C6E5"),  # Light Pink
#     ROOT.TColor.GetColor("#B6D9D7"),  # Light Teal
#     ROOT.TColor.GetColor("#C6E2B6"),  # Light Spring Green
#     ROOT.TColor.GetColor("#D0E1F9"),  # Light Azure
#     ROOT.TColor.GetColor("#E3D6F1"),  # Light Violet
#     ROOT.TColor.GetColor("#000000"),  # Black (for contrast)
#     ROOT.TColor.GetColor("#C0C0C0"),  # Gray
#     ROOT.TColor.GetColor("#F4A3A0"),  # Soft Red
#     ROOT.TColor.GetColor("#9AB9F5"),  # Soft Blue
#     ROOT.TColor.GetColor("#B5E5B0"),  # Soft Green
#     ROOT.TColor.GetColor("#F0A1B2"),  # Soft Magenta
#     ROOT.TColor.GetColor("#B0D6D5"),  # Soft Cyan
#     ROOT.TColor.GetColor("#F3F9A6")   # Soft Yellow
# ]

# def get_yields(pdf, result=None):
#     """Also get uncertainties if you pass a fit result.
#     """
#     yields = np.zeros(x.numBins())
#     yields_uncert = np.zeros(x.numBins())
#     sample_values = {}
#     for i_bin in range(x.numBins()):
#         x.setBin(i_bin)
#         bin_width = x.getBinWidth(i_bin)
#         for i_sample in range(pdf.funcList().size()):
#             sample_yield = ROOT.RooProduct("tmp", "tmp", [pdf.funcList()[i_sample], pdf.coefList()[i_sample]])
#             if i_sample in sample_values:
#                 sample_values[i_sample] += [bin_width * sample_yield.getVal()]
#             else:
#                 sample_values[i_sample] = [bin_width * sample_yield.getVal()]
#             yields[i_bin] += bin_width * sample_yield.getVal()
#             if result is not None:
#                 yields_uncert[i_bin] += bin_width * sample_yield.getPropagatedError(result, observables)
#     return yields, yields_uncert, sample_values
    
# for channel in meas.GetChannels():


#     channel_pdf = ws[str(channel.GetName()) + "_model"]
#     x = ws["obs_x_" + str(channel.GetName())]
#     observables = ROOT.RooArgSet(x)

    
#     divide_value = 0.3
#     # print(channel.GetName())
#     # channel_name = str(channel.GetName())
#     cv_array += [ROOT.TCanvas("canvas" +  str(channel.GetName()), "canvas" +  str(channel.GetName()), 1000, 600)]
#     # cv_array[-1].Divide(2, 1)
#     # cv_array[-1].cd(1)
#     pad1_upper = ROOT.TPad("pad1_upper" + str(channel.GetName()), "pad1_upper" + str(channel.GetName()), 0, divide_value, 0.5, 1)
#     pad1_upper.Draw()
#     pad1_bottom = ROOT.TPad("pad1_bottom" + str(channel.GetName()), "pad1_bottom" + str(channel.GetName()), 0, 0, 0.5, divide_value)
#     pad1_bottom.Draw()

#     pad2_upper = ROOT.TPad("pad2_upper" + str(channel.GetName()), "pad2_upper" + str(channel.GetName()), 0.5, divide_value, 1, 1)
#     pad2_upper.Draw()
#     pad2_bottom = ROOT.TPad("pad2_bottom" + str(channel.GetName()), "pad2_bottom" + str(channel.GetName()), 0.5, 0, 1, divide_value)    
#     pad2_bottom.Draw()
#     # pad2.cd()
#     # print(channel)
#     pad1_upper.cd()
#     data_histogram = channel.GetData().GetHisto()
#     data_histogram.SetStats(0)
#     data_histogram.SetMarkerStyle(8)
#     data_histogram.SetMarkerSize(0.5)


#     prefit_yields, prefit_unc, prefit_sample_values = get_yields(channel_pdf)

#     error_pre_graphs += [ROOT.TGraphErrors("prefit, prefit")]

#     bin_index = 1
#     for bin_value, bin_error in zip(prefit_yields, prefit_unc):
#         # afterfit_histograms[-1].SetBinContent(bin_index, bin_value)
#         # print(afterfit_histograms[-1].GetBinCenter(i), bin_value)
#         error_pre_graphs[-1].SetPoint(bin_index, data_histogram.GetBinCenter(bin_index), bin_value)
#         error_pre_graphs[-1].SetPointError(bin_index, 0, 0)
#         bin_index += 1


#     error_pre_graphs[-1].SetStats(0)
#     error_pre_graphs[-1].SetMarkerStyle(21)
#     error_pre_graphs[-1].SetMarkerSize(0.5)
#     error_pre_graphs[-1].SetLineStyle(0)

#     # pad1_upper.cd()


#     hs_stacks += [ROOT.THStack("hs" + str(channel.GetName()), "hs" + str(channel.GetName()))]
#     sample_histograms = []

#     for i, sample in enumerate(channel.GetSamples()):
#         # print(sample.GetName())
#         # sample_histograms += [sample.GetHisto()]
#         hist = sample.GetHisto()
#         hist.SetFillColor(predefined_colors[i])
#         hist.SetLineColor(ROOT.kBlack)
#         sample_histograms += [hist]
#         hs_stacks[-1].Add(sample_histograms[-1])

#     channel_name = "_".join(channel.GetName().split("_")[1:])
#     hs_stacks[-1].SetTitle(channel_name + " PREFIT")
#     hs_stacks[-1].Draw("hist")
#     error_pre_graphs[-1].Draw("same P")
#     # data_histogram.Draw("same")
#     # pad2.Draw()

#     pad1_bottom.cd()
#     number_of_bins = data_histogram.GetNbinsX()
#     bias_graphs += [ROOT.TGraph(number_of_bins)]
#     bias_graphs[-1].SetTitle("")
#     bias_graphs[-1].SetMarkerSize(0.4)
#     bias_graphs[-1].SetMarkerStyle(8)
#     bin_sums = [0] *  number_of_bins
#     for i in range(hs_stacks[-1].GetNhists()):
#         hist = hs_stacks[-1].GetHists().At(i)
#         for bin in range(1, number_of_bins + 1):
#             bin_sums[bin - 1] += hist.GetBinContent(bin)
    
#     for i in range(1, number_of_bins + 1):
#         bias_graphs[-1].SetPoint(i - 1, data_histogram.GetBinCenter(i), data_histogram.GetBinContent(i) / bin_sums[i - 1])

#     bias_graphs[-1].Draw("AP")
#     minimal_bin_value = data_histogram.GetBinLowEdge(1)
#     maximum_bin_value = data_histogram.GetBinLowEdge(number_of_bins) + data_histogram.GetBinWidth(number_of_bins)


#     bias_graphs[-1].GetYaxis().SetRangeUser(0.5, 1.5)
#     bias_graphs[-1].GetXaxis().SetRangeUser(minimal_bin_value, maximum_bin_value)

#     normal_lines += [ROOT.TLine(minimal_bin_value, 1, maximum_bin_value, 1)]
#     normal_lines[-1].SetLineStyle(2)
#     normal_lines[-1].SetLineWidth(1)
#     normal_lines[-1].Draw("same")
#     # dump_
#     # pad1.cd()
#     # pad1.Draw()
#     # cv_array[-1].Update()

#     pad2_upper.cd()

#     # second_histos += [data_histogram.Clone("second" + str(channel.GetName()))]
#     # second_histos[-1].Draw("same")

#     allSample = channel.GetSamples()

#     # print(ws.allPdfs())
#     # for pdf in ws.allPdfs():
#         # print(pdf.GetName())

#     temp_sample_histogram = []
#     # print(allSample)

#     # for ob in (ws.allVars()):
#         # print(ob.GetName())



#     # bookkeep prefit yields

#     postfit_yields, postfit_yields_uncert, postfit_sample_values = get_yields(channel_pdf, result)

#     afterfit_histograms += [ROOT.TH1F("afterfit", "afterfit", len(postfit_yields), minimal_bin_value, maximum_bin_value)]

#     error_graphs += [ROOT.TGraphErrors("postfit, postfit")]


#     bin_index = 1
#     for bin_value, bin_error in zip(postfit_yields, postfit_yields_uncert):
#         afterfit_histograms[-1].SetBinContent(bin_index, bin_value)
#         # print(afterfit_histograms[-1].GetBinCenter(i), bin_value)
#         error_graphs[-1].SetPoint(bin_index, afterfit_histograms[-1].GetBinCenter(bin_index), bin_value)
#         error_graphs[-1].SetPointError(bin_index, 0, bin_error)
#         # print(bin_value, bin_error)
#         bin_index += 1


#     error_graphs[-1].SetStats(0)
#     error_graphs[-1].SetMarkerStyle(21)
#     error_graphs[-1].SetMarkerSize(0.5)
#     error_graphs[-1].SetLineStyle(0)


#     second_hs_stacks += [ROOT.THStack("fitted stack", "fitted stack")]

#     # temp_histos
#     color_number = 0

#     for postfit in postfit_sample_values:
#         temp_histo = ROOT.TH1F("afterfit" + str(postfit), "afterfit" + str(postfit), len(postfit_yields), minimal_bin_value, maximum_bin_value)
#         temp_histo.SetFillColor(predefined_colors[color_number])
#         color_number += 1
#         bin_index = 1
#         for bin_value in postfit_sample_values[postfit]:
#             temp_histo.SetBinContent(bin_index, bin_value)
#             bin_index += 1
#         second_hs_stacks[-1].Add(temp_histo)

#     channel_name = "_".join(channel.GetName().split("_")[1:])
#     second_hs_stacks[-1].SetTitle(channel_name + " POSTFIT")
#     second_hs_stacks[-1].Draw("hist")
#     # afterfit_histograms[-1].Draw("same p")
#     error_graphs[-1].Draw("same P")

#     bias_second_graphs += [ROOT.TGraph(number_of_bins)]
#     bias_second_graphs[-1].SetTitle("")
#     bias_second_graphs[-1].SetMarkerSize(0.4)
#     bias_second_graphs[-1].SetMarkerStyle(8)

#     for i in range(1, number_of_bins + 1):
#         bias_second_graphs[-1].SetPoint(i - 1, data_histogram.GetBinCenter(i), data_histogram.GetBinContent(i) / postfit_yields[i - 1])

#     pad2_bottom.cd()
#     bias_second_graphs[-1].Draw("AP")
#     minimal_bin_value = data_histogram.GetBinLowEdge(1)
#     maximum_bin_value = data_histogram.GetBinLowEdge(number_of_bins) + data_histogram.GetBinWidth(number_of_bins)


#     bias_second_graphs[-1].GetYaxis().SetRangeUser(0.5, 1.5)
#     bias_second_graphs[-1].GetXaxis().SetRangeUser(minimal_bin_value, maximum_bin_value)

#     normal_lines += [ROOT.TLine(minimal_bin_value, 1, maximum_bin_value, 1)]
#     normal_lines[-1].SetLineStyle(2)
#     normal_lines[-1].SetLineWidth(1)
#     normal_lines[-1].Draw("same")

#     print(postfit_sample_values)
#     # print(prefit_yields)


#     cv_array[-1].SaveAs("histo.png")
#     cv_array[-1].Draw()
#     # break