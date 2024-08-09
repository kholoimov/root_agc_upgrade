import ROOT
import math
import numpy as np
class AGC_Sample(ROOT.RooStats.HistFactory.Sample):
    def __init__(self, name, histo_name, histo_file, histo_path = ""):
        ROOT.RooStats.HistFactory.Sample.__init__(self, name, histo_name, histo_file, histo_path)
        self.output_path = 'ROOT_AGC_Utils/HistFactoryExtra.root'
        self.fInputFile = ""

    def __del__(self):
        # Custom cleanup code here
        print(f"AGC_Sample instance {self.GetName()} is being destroyed")

    def SetSystematicsInputFile(self, file):
        self.fInputFile = file

    def AddHistoSys(self, name, histoname_up = None, histofile_up = None, histopath_up = "",
                    histoname_down = None, histofile_down = None, histopath_down = ""):
        
        histSys = ROOT.RooStats.HistFactory.HistoSys()

        if (histofile_up is None):
            histofile_down  =   self.fInputFile
            histofile_up    =   self.fInputFile

        histSys.SetName(name)
        histSys.SetHistoName(self.GetHistoName())
        histSys.SetHistoPathHigh(histopath_up)
        histSys.SetHistoPathLow(histopath_down)
        histSys.SetHistoNameHigh(histoname_up)
        histSys.SetHistoNameLow(histoname_down)
        histSys.SetInputFileHigh(histofile_up)
        histSys.SetInputFileLow(histofile_down)

        self.GetHistoSysList().push_back(histSys)

    def AddHistoFactor(self, name, histoname_up = None, histofile_up = None, histopath_up = "",
                    histoname_down = None, histofile_down = None, histopath_down = ""):
        
        histFactor = ROOT.RooStats.HistFactory.HistoFactor()

        if (histofile_up is None):
            histofile_down  =   self.fInputFile
            histofile_up    =   self.fInputFile

        histFactor.SetName(name)
        histFactor.SetHistoName(self.GetHistoName())
        histFactor.SetHistoPathHigh(histopath_up)
        histFactor.SetHistoPathLow(histopath_down)
        histFactor.SetHistoNameHigh(histoname_up)
        histFactor.SetHistoNameLow(histoname_down)
        histFactor.SetInputFileHigh(histofile_up)
        histFactor.SetInputFileLow(histofile_down)

        self.GetHistoFactorList().push_back(histFactor)

    def AddShapeSys(self, name, constraint_type = None, histoname = None, histofile = None, histopath = ""):
        
        shapeSys = ROOT.RooStats.HistFactory.ShapeSys()

        if histofile is None:
            histofile = self.fInputFile

        shapeSys.SetName(name)
        shapeSys.SetHistoName(self.GetHistoName())
        shapeSys.SetHistoPath(histopath)
        shapeSys.SetInputFile(histofile)

        shapeSys.SetConstraintType(constraint_type)

        self.GetShapeSysList().push_back(shapeSys)



    def AddNormPlusShapeHistoSys(self, name, histoname_up = None, histofile_up = None, histopath_up = "",
                                 histoname_down = None, histofile_down = None, histopath_down = ""):
        if histofile_up is None:
            histofile_up = self.fInputFile
            histofile_down = self.fInputFile
            assert histofile_up != "", "ERROR: You not specified input file for sample"

        if histoname_down is None:
            self.Symmetrize_AddNormPlusShapeHistoSys(name, histoname_up, histofile_up, histopath_up)
        else:
            self.NonSymmetrize_AddNormPlusShapeHistoSys(name, histoname_up, histofile_up, histopath_up,
                                 histoname_down, histofile_down, histopath_down)
            
    def Symmetrize_AddNormPlusShapeHistoSys(self, name, histoname, histofile, histopath):
        
        channel_name = str(histoname.split("_")[0])

        file = ROOT.TFile(histofile, "READ")
        dir = file.GetDirectory(histopath)
        hist_top = (dir.Get(histoname))

        hist_nominal_file = ROOT.TFile(self.GetInputFile(), "READ")
        hist_nominal_name = self.GetHistoName()
        hist_nominal_directory = hist_nominal_file.GetDirectory(self.GetHistoPath())
        hist_nominal = (hist_nominal_directory.Get(hist_nominal_name))

        norm_factor_up = hist_top.Integral() / hist_nominal.Integral()
        h_new = hist_top.Clone(channel_name + "_" + str(self.GetName()) + "_" + name  + "_norm_plus_shape_up_clone")
        h_new.Scale(1/norm_factor_up)

        h_down = hist_nominal.Clone(channel_name + "_" + str(self.GetName()) + "_" + name  + "_norm_plus_shape_down_clone")
        h_down.Scale(2)
        h_down.Add(h_new, -1)

        output_file = ROOT.TFile(self.output_path, "UPDATE")

        hist_up_name  = str(channel_name + "_" + str(self.GetName()) + "_" + name  + "_norm_plus_shape_up")
        hist_down_name = str(channel_name + "_" + str(self.GetName()) + "_" + name  + "_norm_plus_shape_down")

        # print(h_down)
        # print(h_new)

        h_new.Write(hist_up_name)
        h_down.Write(hist_down_name)

        output_file.Close()

        histSys = ROOT.RooStats.HistFactory.HistoSys()

        histSys.SetName(name)
        # histSys.SetHistoName(self.GetselfName())
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

        self.GetHistoSysList().push_back(histSys)
        self.GetOverallSysList().push_back(overallSys)


    def NonSymmetrize_AddNormPlusShapeHistoSys(self, name, histoname_up, histofile_up, histopath_up,
                                 histoname_down, histofile_down, histopath_down):
        channel_name = str(histoname_up.split("_")[0])

        file = ROOT.TFile(histofile_up, "READ")
        dir = file.GetDirectory(histopath_up)
        hist_top = (dir.Get(histoname_up))
        
        hist_nominal_file = ROOT.TFile(self.GetInputFile(), "READ")
        hist_nominal_name = self.GetHistoName()
        hist_nominal_directory = hist_nominal_file.GetDirectory(self.GetHistoPath())
        hist_nominal = (hist_nominal_directory.Get(hist_nominal_name))

        norm_factor_up = hist_top.Integral() / hist_nominal.Integral()
        h_new_up = hist_top.Clone(channel_name + "_" + str(self.GetName()) + "_" + name  + "_norm_plus_shape_up_clone")
        h_new_up.Scale(1/norm_factor_up)

        file_down = ROOT.TFile(histofile_down, "READ")
        dir_down = file_down.GetDirectory(histopath_down)
        hist_down = (dir_down.Get(histoname_down))

        norm_factor_down = hist_down.Integral() / hist_nominal.Integral()
        h_new_down = hist_down.Clone(channel_name + "_" + str(self.GetName()) + "_" + name  + "_norm_plus_shape_down_clone")
        h_new_down.Scale(1/norm_factor_down)

        output_file = ROOT.TFile(self.output_path, "UPDATE")

        hist_up_name  = str(channel_name + "_" + str(self.GetName()) + "_" + name  + "_norm_plus_shape_up")
        hist_down_name = str(channel_name + "_" + str(self.GetName()) + "_" + name  + "_norm_plus_shape_down")

        # print(h_new_down)
        # print(h_new_up)

        h_new_up.Write(hist_up_name)
        h_new_down.Write(hist_down_name)

        output_file.Close()

        histSys = ROOT.RooStats.HistFactory.HistoSys()

        histSys.SetName(name)
        # histSys.SetHistoName(self.GetselfName())
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

        self.GetHistoSysList().push_back(histSys)
        self.GetOverallSysList().push_back(overallSys)


        # hist_nominal = self.Get


class RebinningTool:

    def __init__(self):
        self.fHistoBinLow = None
        self.fHistoBinHigh = None
        self.fRebin = None
        self.fRebinning = False

    def set_xmin(self, xmin):
        self.fHistoBinLow = xmin

    def set_xmax(self, xmax):
        self.fHistoBinHigh = xmax

    def set_rebin(self, rebin):
        self.fRebinning = True
        self.fRebin = rebin

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

    def is_integer(self, value):
        return int(value) == value
    
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
       
    def apply_rebinning(self, input_file_path = None, output_file_path = None):
        if input_file_path is None:
            input_file_path = self.input_path
        if output_file_path is None:
            output_file_path = self.output_path
            
        file = ROOT.TFile(input_file_path, "READ")

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

    def set_input_path(self, input_path):
        self.input_path = input_path

    def set_output_path(self, output_path):
        self.output_path = output_path
        file = ROOT.TFile(self.output_path, "RECREATE")
        file.Close()

def MuteTool():
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Minimization);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Eval);


class Visualization:
    def CreateAndSavePicture(self, filename, fWorkspace):
        # std::vector<TString> names;
        # std::vector<double> values;
        # std::vector<double> errors;
        names = []
        values = []
        errors = []

        vars = fWorkspace.allVars();
        # for (auto var : vars)
        for var in vars:
            name = var.GetName();
            if ("alpha" in name):
                if (fWorkspace.var(name).getVal() == 0): continue
                names += [name];
                values += [fWorkspace.var(name).getVal()]
                errors += [fWorkspace.var(name).getError()]

        num = len(names)

        ROOT.gStyle.SetPalette(1);
        c1 =  ROOT.TCanvas("c1", "c1", 1200, 600)
        c1.SetLeftMargin(0.2);

        

        frame =  ROOT.TH2F("frame", "", 6, -3, 3, num, 0, num);


        frame.Draw("");


        box = ROOT.TBox(-2, 0, 2, num);
        box.SetFillColorAlpha(ROOT.kYellow - 9, 0.8);
        box.Draw("same A");

        box_internal = ROOT.TBox(-1, 0, 1, num);
        box_internal.SetFillColorAlpha(ROOT.kGreen - 9, 0.5);
        box_internal.Draw("same A");

        axis = ROOT.TGaxis(-3, num, 3, num,-3, 3, 510,"-");
        xaxis = ROOT.TGaxis(3, 0, 3, num, 0, num, num ,"+S");

        frame.GetYaxis().SetTickLength(0.);
        xaxis.SetTickLength(0.);

        # // frame.Draw("X+");

        graph = ROOT.TGraph();
        lines = []
        right_ticks = []
        left_ticks = []

        # for (int i = 0; i < num; ++i)
        for i in range(num):
            graph.SetPoint(i,  values[i], i + 0.5);
            frame.GetYaxis().SetBinLabel(i + 1, names[i]);

            lines += [ROOT.TLine(values[i] - errors[i], i + 0.5, values[i] + errors[i], i + 0.5)]
            lines[-1].SetLineColor(ROOT.kBlack);
            lines[-1].SetLineWidth(2);
            lines[-1].Draw();

            left_ticks += [ROOT.TLine(-3, i + 0.5, -2.95, i + 0.5)]
            left_ticks[-1].SetLineColor(ROOT.kBlack);
            left_ticks[-1].SetLineWidth(1);
            left_ticks[-1].Draw();

            right_ticks += [ROOT.TLine(2.95, i + 0.5, 3, i + 0.5)]
            right_ticks[-1].SetLineColor(ROOT.kBlack);
            right_ticks[-1].SetLineWidth(1);
            right_ticks[-1].Draw();


        tl = ROOT.TLine(0, 0, 0, num);
        tl.SetLineStyle(2);
        tl.Draw();


        graph.SetMarkerStyle(20);
        graph.SetMarkerSize(1);
        graph.SetMarkerColor(ROOT.kBlack);

        graph.Draw("P same");

        frame.SetStats(0);
        axis.SetLabelSize(0.0);
        axis.Draw();

        xaxis.SetLabelSize(0.0);
        xaxis.Draw();

        ROOT.gPad.RedrawAxis(); 

        c1.SaveAs(filename);
        c1.Draw();

        
class DrawModel:

    predefined_colors = [
        ROOT.TColor.GetColor("#F6C3C0"),  # Light Red
        ROOT.TColor.GetColor("#A9D0F5"),  # Light Blue
        ROOT.TColor.GetColor("#B6D7A8"),  # Light Green
        ROOT.TColor.GetColor("#EAD1DC"),  # Light Magenta
        ROOT.TColor.GetColor("#D0E9E8"),  # Light Cyan
        ROOT.TColor.GetColor("#F6F9C6"),  # Light Yellow
        ROOT.TColor.GetColor("#F6B8A3"),  # Light Orange
        ROOT.TColor.GetColor("#F9C6E5"),  # Light Pink
        ROOT.TColor.GetColor("#B6D9D7"),  # Light Teal
        ROOT.TColor.GetColor("#C6E2B6"),  # Light Spring Green
        ROOT.TColor.GetColor("#D0E1F9"),  # Light Azure
        ROOT.TColor.GetColor("#E3D6F1"),  # Light Violet
        ROOT.TColor.GetColor("#000000"),  # Black (for contrast)
        ROOT.TColor.GetColor("#C0C0C0"),  # Gray
        ROOT.TColor.GetColor("#F4A3A0"),  # Soft Red
        ROOT.TColor.GetColor("#9AB9F5"),  # Soft Blue
        ROOT.TColor.GetColor("#B5E5B0"),  # Soft Green
        ROOT.TColor.GetColor("#F0A1B2"),  # Soft Magenta
        ROOT.TColor.GetColor("#B0D6D5"),  # Soft Cyan
        ROOT.TColor.GetColor("#F3F9A6")   # Soft Yellow
    ]
    
    def __init__(self, meas, ws):
        self.meas = meas
        self.ws = ws
        self.cv_array = []
        self.hs_stacks = []
        self.sample_histograms = []
        self.bias_graphs = []
        self.bias_second_graphs = []
        self.normal_lines = []
        self.second_histos = []
        self.second_hs_stacks = []
        self.afterfit_histograms = []
        self.error_graphs = []
        self.error_pre_graphs = []

    def get_yields(self, variable, observables, pdf, result=None):
        """Also get uncertainties if you pass a fit result.
        """
        yields = np.zeros(variable.numBins())
        yields_uncert = np.zeros(variable.numBins())
        sample_values = {}
        for i_bin in range(variable.numBins()):
            variable.setBin(i_bin)
            bin_width = variable.getBinWidth(i_bin)
            for i_sample in range(pdf.funcList().size()):
                sample_yield = ROOT.RooProduct("tmp", "tmp", [pdf.funcList()[i_sample], pdf.coefList()[i_sample]])
                if i_sample in sample_values:
                    sample_values[i_sample] += [bin_width * sample_yield.getVal()]
                else:
                    sample_values[i_sample] = [bin_width * sample_yield.getVal()]
                yields[i_bin] += bin_width * sample_yield.getVal()
                if result is not None:
                    yields_uncert[i_bin] += bin_width * sample_yield.getPropagatedError(result, observables)
        return yields, yields_uncert, sample_values
    
    def get_my_propagated_error(self, roo_abs_real, fit_result, param_set):
        # Extract parameters from RooAbsReal and fit result
        all_params_in_abs_real = ROOT.RooArgSet()
        roo_abs_real.getParameters(param_set, all_params_in_abs_real)

        param_list = fit_result.floatParsInit()
        
        plus_var = []
        minus_var = []
        # plus_var.reserve(param_list.size())
        # minus_var.reserve(param_list.size())

        # V = fit_result.covarianceMatrix() if param_list.size() == fit_result.floatParsInit().size() else fit_result.reducedCovarianceMatrix(param_list)

        prefit_covariance_matrix = np.zeros((param_list.size(), param_list.size()))
        for i in range(param_list.size()):
            prefit_covariance_matrix[i][i] = param_list[i].getError()**2

        V = prefit_covariance_matrix



        for ivar in range(param_list.size()):
            rrv = param_list.at(ivar)

            cen_val = rrv.getVal()
            err_val = np.sqrt(V[ivar][ivar])

            # Make Plus variation
            rrv.setVal(cen_val + err_val)
            plus_var.append(all_params_in_abs_real[0].getVal())

            # Make Minus variation
            rrv.setVal(cen_val - err_val)
            minus_var.append(all_params_in_abs_real[0].getVal())

            rrv.setVal(cen_val)
        
        # Re-evaluate to ensure the state is consistent
        # roo_abs_real.getVal(param_set)

        # C = ROOT.TMatrixDSym(param_list.size())
        C = np.zeros((param_list.size(), param_list.size()))
        err_vec = np.zeros(param_list.size())

        for i in range(param_list.size()):
            err_vec[i] = np.sqrt(V[i][i])
            for j in range(i, param_list.size()):
                C[i][j] = V[i][j] / np.sqrt(V[i][i] * V[j][j])
                C[j][i] = C[i][j]
        
        F = np.zeros(len(plus_var))
        for j in range(len(plus_var)):
            F[j] = (plus_var[j] - minus_var[j]) * 0.5

        sum_error = np.dot(F, np.dot(V, F))

        print(sum_error)
        
        return np.sqrt(sum_error)

    # Example usage
    # Assuming `roo_abs_real`, `fit_result`, and `param_set` are properly defined RooFit objects
    # error = get_propagated_error(roo_abs_real, fit_result, param_set)
    # print("Propagated Error:", error)


    def get_prefit_yields(self, variable, observables, pdf, result=None):
        """Also get uncertainties if you pass a fit result.
        """
        yields = np.zeros(variable.numBins())
        yields_uncert = np.zeros(variable.numBins())
        sample_values = {}
        for i_bin in range(variable.numBins()):
            variable.setBin(i_bin)
            bin_width = variable.getBinWidth(i_bin)
            for i_sample in range(pdf.funcList().size()):
                sample_yield = ROOT.RooProduct("tmp", "tmp", [pdf.funcList()[i_sample], pdf.coefList()[i_sample]])
                if i_sample in sample_values:
                    sample_values[i_sample] += [bin_width * sample_yield.getVal()]
                else:
                    sample_values[i_sample] = [bin_width * sample_yield.getVal()]
                yields[i_bin] += bin_width * sample_yield.getVal()
                if result is not None:
                    yields_uncert[i_bin] += bin_width * self.get_my_propagated_error(sample_yield, result, observables)
        return yields, yields_uncert, sample_values

    
    def Draw(self, result):
        for channel in self.meas.GetChannels():
            channel_pdf = self.ws[str(channel.GetName()) + "_model"]
            obs_var = self.ws["obs_x_" + str(channel.GetName())]
            observables = ROOT.RooArgSet(obs_var)

            
            divide_value = 0.3
            # print(channel.GetName())
            # channel_name = str(channel.GetName())
            self.cv_array += [ROOT.TCanvas("canvas" +  str(channel.GetName()), "canvas" +  str(channel.GetName()), 1500, 600)]
            # self.cv_array[-1].Divide(2, 1)
            # self.cv_array[-1].cd(1)
            pad1_upper = ROOT.TPad("pad1_upper" + str(channel.GetName()), "pad1_upper" + str(channel.GetName()), 0, divide_value, 0.5, 1)
            pad1_upper.Draw()
            pad1_bottom = ROOT.TPad("pad1_bottom" + str(channel.GetName()), "pad1_bottom" + str(channel.GetName()), 0, 0, 0.5, divide_value)
            pad1_bottom.Draw()

            pad2_upper = ROOT.TPad("pad2_upper" + str(channel.GetName()), "pad2_upper" + str(channel.GetName()), 0.5, divide_value, 1, 1)
            pad2_upper.Draw()
            pad2_bottom = ROOT.TPad("pad2_bottom" + str(channel.GetName()), "pad2_bottom" + str(channel.GetName()), 0.5, 0, 1, divide_value)    
            pad2_bottom.Draw()
            # pad2.cd()
            # print(channel)
            pad1_upper.cd()
            data_histogram = channel.GetData().GetHisto()
            data_histogram.SetStats(0)
            data_histogram.SetMarkerStyle(8)
            data_histogram.SetMarkerSize(0.5)


            prefit_yields, prefit_unc, prefit_sample_values = self.get_prefit_yields(obs_var, observables, channel_pdf, result)

            self.error_pre_graphs += [ROOT.TGraphErrors("prefit, prefit")]

            bin_index = 1
            for bin_value, bin_error in zip(prefit_yields, prefit_unc):
                # afterfit_histograms[-1].SetBinContent(bin_index, bin_value)
                # print(afterfit_histograms[-1].GetBinCenter(i), bin_value)
                self.error_pre_graphs[-1].SetPoint(bin_index, data_histogram.GetBinCenter(bin_index), bin_value)
                self.error_pre_graphs[-1].SetPointError(bin_index, 0, bin_error)
                bin_index += 1


            self.error_pre_graphs[-1].SetStats(0)
            self.error_pre_graphs[-1].SetMarkerStyle(21)
            self.error_pre_graphs[-1].SetMarkerSize(0.5)
            self.error_pre_graphs[-1].SetLineStyle(0)

            # pad1_upper.cd()


            self.hs_stacks += [ROOT.THStack("hs" + str(channel.GetName()), "hs" + str(channel.GetName()))]
            sample_histograms = []

            for i, sample in enumerate(channel.GetSamples()):
                # print(sample.GetName())
                # sample_histograms += [sample.GetHisto()]
                hist = sample.GetHisto()
                hist.SetFillColor(self.predefined_colors[i])
                hist.SetLineColor(ROOT.kBlack)
                sample_histograms += [hist]
                self.hs_stacks[-1].Add(sample_histograms[-1])

            channel_name = "_".join(channel.GetName().split("_")[1:])
            self.hs_stacks[-1].SetTitle(channel_name + " PREFIT")
            self.hs_stacks[-1].Draw("hist")
            self.error_pre_graphs[-1].Draw("same P")
            # data_histogram.Draw("same")
            # pad2.Draw()

            pad1_bottom.cd()
            number_of_bins = data_histogram.GetNbinsX()
            self.bias_graphs += [ROOT.TGraph(number_of_bins)]
            self.bias_graphs[-1].SetTitle("")
            self.bias_graphs[-1].SetMarkerSize(0.4)
            self.bias_graphs[-1].SetMarkerStyle(8)
            bin_sums = [0] *  number_of_bins
            for i in range(self.hs_stacks[-1].GetNhists()):
                hist = self.hs_stacks[-1].GetHists().At(i)
                for bin in range(1, number_of_bins + 1):
                    bin_sums[bin - 1] += hist.GetBinContent(bin)
            
            for i in range(1, number_of_bins + 1):
                self.bias_graphs[-1].SetPoint(i - 1, data_histogram.GetBinCenter(i), data_histogram.GetBinContent(i) / bin_sums[i - 1])

            self.bias_graphs[-1].Draw("AP")
            minimal_bin_value = data_histogram.GetBinLowEdge(1)
            maximum_bin_value = data_histogram.GetBinLowEdge(number_of_bins) + data_histogram.GetBinWidth(number_of_bins)


            self.bias_graphs[-1].GetYaxis().SetRangeUser(0.5, 1.5)
            self.bias_graphs[-1].GetXaxis().SetRangeUser(minimal_bin_value, maximum_bin_value)

            self.normal_lines += [ROOT.TLine(minimal_bin_value, 1, maximum_bin_value, 1)]
            self.normal_lines[-1].SetLineStyle(2)
            self.normal_lines[-1].SetLineWidth(1)
            self.normal_lines[-1].Draw("same")
            # dump_
            # pad1.cd()
            # pad1.Draw()
            # self.cv_array[-1].Update()

            pad2_upper.cd()

            # second_histos += [data_histogram.Clone("second" + str(channel.GetName()))]
            # second_histos[-1].Draw("same")

            # allSample = channel.GetSamples()

            # print(self.ws.allPdfs())
            # for pdf in self.ws.allPdfs():
                # print(pdf.GetName())

            # temp_sample_histogram = []
            # print(allSample)

            # for ob in (self.ws.allVars()):
                # print(ob.GetName())



            # bookkeep prefit yields

            postfit_yields, postfit_yields_uncert, postfit_sample_values = self.get_yields(obs_var, observables, channel_pdf, result)

            self.afterfit_histograms += [ROOT.TH1F("afterfit", "afterfit", len(postfit_yields), minimal_bin_value, maximum_bin_value)]

            self.error_graphs += [ROOT.TGraphErrors("postfit, postfit")]


            bin_index = 1
            for bin_value, bin_error in zip(postfit_yields, postfit_yields_uncert):
                self.afterfit_histograms[-1].SetBinContent(bin_index, bin_value)
                # print(self.afterfit_histograms[-1].GetBinCenter(i), bin_value)
                self.error_graphs[-1].SetPoint(bin_index, self.afterfit_histograms[-1].GetBinCenter(bin_index), bin_value)
                self.error_graphs[-1].SetPointError(bin_index, 0, bin_error)
                # print(bin_value, bin_error)
                bin_index += 1


            self.error_graphs[-1].SetStats(0)
            self.error_graphs[-1].SetMarkerStyle(21)
            self.error_graphs[-1].SetMarkerSize(0.5)
            self.error_graphs[-1].SetLineStyle(0)


            self.second_hs_stacks += [ROOT.THStack("fitted stack", "fitted stack")]

            # temp_histos
            color_number = 0

            for postfit in postfit_sample_values:
                temp_histo = ROOT.TH1F("afterfit" + str(postfit), "afterfit" + str(postfit), len(postfit_yields), minimal_bin_value, maximum_bin_value)
                temp_histo.SetFillColor(self.predefined_colors[color_number])
                color_number += 1
                bin_index = 1
                for bin_value in postfit_sample_values[postfit]:
                    temp_histo.SetBinContent(bin_index, bin_value)
                    bin_index += 1
                self.second_hs_stacks[-1].Add(temp_histo)

            channel_name = "_".join(channel.GetName().split("_")[1:])
            self.second_hs_stacks[-1].SetTitle(channel_name + " POSTFIT")
            self.second_hs_stacks[-1].Draw("hist")
            # afterfit_histograms[-1].Draw("same p")
            self.error_graphs[-1].Draw("same P")

            self.bias_second_graphs += [ROOT.TGraph(number_of_bins)]
            self.bias_second_graphs[-1].SetTitle("")
            self.bias_second_graphs[-1].SetMarkerSize(0.4)
            self.bias_second_graphs[-1].SetMarkerStyle(8)

            for i in range(1, number_of_bins + 1):
                self.bias_second_graphs[-1].SetPoint(i - 1, data_histogram.GetBinCenter(i), data_histogram.GetBinContent(i) / postfit_yields[i - 1])

            pad2_bottom.cd()
            self.bias_second_graphs[-1].Draw("AP")
            minimal_bin_value = data_histogram.GetBinLowEdge(1)
            maximum_bin_value = data_histogram.GetBinLowEdge(number_of_bins) + data_histogram.GetBinWidth(number_of_bins)


            self.bias_second_graphs[-1].GetYaxis().SetRangeUser(0.5, 1.5)
            self.bias_second_graphs[-1].GetXaxis().SetRangeUser(minimal_bin_value, maximum_bin_value)

            self.normal_lines += [ROOT.TLine(minimal_bin_value, 1, maximum_bin_value, 1)]
            self.normal_lines[-1].SetLineStyle(2)
            self.normal_lines[-1].SetLineWidth(1)
            self.normal_lines[-1].Draw("same")

            # print(postfit_sample_values)
            # print(prefit_yields)


            self.cv_array[-1].SaveAs("histo.png")
            self.cv_array[-1].Draw()
