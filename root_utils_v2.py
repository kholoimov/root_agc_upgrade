import ROOT
import math
import numpy as np
class AGC_Sample(ROOT.RooStats.HistFactory.Sample):
    def __init__(self, name, histo_name, histo_file, histo_path = ""):
        ROOT.RooStats.HistFactory.Sample.__init__(self, name, histo_name, histo_file, histo_path)
        self.output_path = 'data/HistFactoryExtra.root'
        self.fInputFile = ""

    def __del__(self):
        # Custom cleanup code here
        pass
        # print(f"AGC_Sample instance {self.GetName()} is being destroyed")

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
    def CreateAndSaveResultsPicture(self, filename, fWorkspace):
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
        self.c1 =  ROOT.TCanvas("c1", "c1", 1200, 600)
        self.c1.SetLeftMargin(0.2);

        

        self.frame =  ROOT.TH2F("self.frame", "", 6, -3, 3, num, 0, num);


        self.frame.Draw("");


        self.box = ROOT.TBox(-2, 0, 2, num);
        self.box.SetFillColorAlpha(ROOT.kYellow - 9, 0.8);
        self.box.Draw("same A");

        self.box_internal = ROOT.TBox(-1, 0, 1, num);
        self.box_internal.SetFillColorAlpha(ROOT.kGreen, 0.5);
        self.box_internal.Draw("same A");

        self.axis = ROOT.TGaxis(-3, num, 3, num,-3, 3, 510,"-");
        self.xaxis = ROOT.TGaxis(3, 0, 3, num, 0, num, num ,"+S");

        self.frame.GetYaxis().SetTickLength(0.);
        self.xaxis.SetTickLength(0.);

        # // self.frame.Draw("X+");

        self.graph = ROOT.TGraph();
        self.lines = []
        self.right_ticks = []
        self.left_ticks = []

        # for (int i = 0; i < num; ++i)
        for i in range(num):
            self.graph.SetPoint(i,  values[i], i + 0.5);
            self.frame.GetYaxis().SetBinLabel(i + 1, names[i]);

            self.lines += [ROOT.TLine(values[i] - errors[i], i + 0.5, values[i] + errors[i], i + 0.5)]
            self.lines[-1].SetLineColor(ROOT.kBlack);
            self.lines[-1].SetLineWidth(2);
            self.lines[-1].Draw();

            self.left_ticks += [ROOT.TLine(-3, i + 0.5, -2.95, i + 0.5)]
            self.left_ticks[-1].SetLineColor(ROOT.kBlack);
            self.left_ticks[-1].SetLineWidth(1);
            self.left_ticks[-1].Draw();

            self.right_ticks += [ROOT.TLine(2.95, i + 0.5, 3, i + 0.5)]
            self.right_ticks[-1].SetLineColor(ROOT.kBlack);
            self.right_ticks[-1].SetLineWidth(1);
            self.right_ticks[-1].Draw();


        self.tl = ROOT.TLine(0, 0, 0, num);
        self.tl.SetLineStyle(2);
        self.tl.Draw();


        self.graph.SetMarkerStyle(20);
        self.graph.SetMarkerSize(1);
        self.graph.SetMarkerColor(ROOT.kBlack);

        self.graph.Draw("P same");

        self.frame.SetStats(0);
        self.axis.SetLabelSize(0.0);
        self.axis.Draw();

        self.xaxis.SetLabelSize(0.0);
        self.xaxis.Draw();

        ROOT.gPad.RedrawAxis(); 

        self.c1.SaveAs(filename);
        self.c1.Draw();

    def DrawCorrelationMatrix(self, filename, result):
        final_parameters = result.floatParsFinal()
        corr_matrix_before = result.correlationMatrix()

        number_of_inter_params = 0

        for i in range(len(final_parameters)):
            par = final_parameters.at(i)
            if "gamma" in par.GetName():
                continue

            number_of_inter_params += 1

        name = "CorrelationMatrix for fit results"

        n = corr_matrix_before.GetNcols()

        self.hh = ROOT.TH2D(name, name, number_of_inter_params, 0, number_of_inter_params, number_of_inter_params, 0, number_of_inter_params)

        internal_index = 0
        for i in range(n):

            par = final_parameters.at(i)
            if "gamma" in par.GetName():
                continue
            
            internal__internal_index = 0
            for j in range(n):
                par = final_parameters.at(j)
                if "gamma" in par.GetName():
                    continue
                self.hh.Fill(internal_index + 0.5, number_of_inter_params - internal__internal_index - 0.5, corr_matrix_before[i][j])
                internal__internal_index += 1    

            if par.GetName() == "tt_norm":
                print(i)

            self.hh.GetXaxis().SetBinLabel(internal_index + 1, final_parameters[i].GetName().split("alpha_")[-1])
            self.hh.GetYaxis().SetBinLabel(number_of_inter_params - internal_index, final_parameters[i].GetName().split("alpha_")[-1])
            internal_index += 1

        self.hh.SetMinimum(-1)
        self.hh.SetMaximum(+1)


        self.c = ROOT.TCanvas("self.c", "Canvas", number_of_inter_params * 100, number_of_inter_params * 60)
        self.hh.Draw("COLZ")
        self.hh.SetStats(0)

        ROOT.gStyle.SetPalette(87)
        palette = self.hh.GetListOfFunctions().FindObject("palette")
        if palette:
            palette.SetX1NDC(0.1)  # Adjust palette position
            palette.SetX2NDC(0.3)  # Adjust palette position

        # Show the canvas
        self.c.SaveAs(filename);
        self.c.Draw()


        
class DrawModel:

    predefined_colors = [
        ROOT.TColor.GetColor("#3F90DA"),  # Light Red
        ROOT.TColor.GetColor("#FFA90E"),  # Light Blue
        ROOT.TColor.GetColor("#BD1F01"),  # Light Green
        ROOT.TColor.GetColor("#94A4A2"),  # Light Magenta
        ROOT.TColor.GetColor("#832DB6"),  # Light Cyan
        ROOT.TColor.GetColor("#A96B59"),  # Light Yellow
        ROOT.TColor.GetColor("#E76300"),  # Light Orange
        ROOT.TColor.GetColor("#B9AC70"),  # Light Pink
        ROOT.TColor.GetColor("#717581"),  # Light Teal
        ROOT.TColor.GetColor("#92DADD"),  # Light Spring Green

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

    def get_yields(self, variable, observables, pdf, result, prefit = False):
        """Also get uncertainties if you pass a fit result.
        """
        yields = np.zeros(variable.numBins())
        yields_uncert = np.zeros(variable.numBins())
        sample_values = {}
        for i_bin in range(variable.numBins()):
            variable.setBin(i_bin)
            bin_width = variable.getBinWidth(i_bin)

            if prefit:
                fill_array = np.zeros((pdf.funcList().size(), len(pdf.getParameters(observables))))
                all_params = pdf.getParameters(observables)
                param_values = result.floatParsInit()

                for i_sample in range(pdf.funcList().size()):
                    sample_yield = ROOT.RooProduct("tmp", "tmp", [pdf.funcList()[i_sample], pdf.coefList()[i_sample]])

                    for j_parameter, par in enumerate(all_params):
                        
                        name = par.GetName()

                        original_ind = param_values.index(param_values.find(name))

                        postfit_cen_val = par.getVal()
                        cen_val = param_values[original_ind].getVal()
                        par_err = param_values[original_ind].getError()

                        par.setVal(cen_val + par_err)
                        par_upper_variation = sample_yield.getVal(observables)
                        par.setVal(cen_val - par_err)
                        par_bottom_variation = sample_yield.getVal(observables)

                        par.setVal(postfit_cen_val)

                        fill_array[i_sample, j_parameter] = (par_upper_variation - par_bottom_variation) / 2

                total_uncertanties_per_variation = np.sum(fill_array, axis=0)
                total_uncertainty = np.sum(np.power(total_uncertanties_per_variation, 2), axis=0)
                yields_uncert[i_bin] = np.sqrt(total_uncertainty) * bin_width        

            else:
                all_pdf_params = pdf.getParameters(observables)

                required_params = [i for i in all_pdf_params if i.getError() > 0.001]

                fill_array = np.zeros((pdf.funcList().size(), len(required_params)))

                number_of_parameters = len(required_params)
                corr_matrix = result.reducedCovarianceMatrix(required_params)

                numpy_corr_matrix = np.zeros((number_of_parameters, number_of_parameters))

                for i_index in range(number_of_parameters):
                    for j_index in range(i_index, number_of_parameters):
                        numpy_corr_matrix[i_index, j_index] = corr_matrix[i_index, j_index] / np.sqrt(corr_matrix[i_index, i_index] * corr_matrix[j_index, j_index])
                        numpy_corr_matrix[j_index, i_index] = numpy_corr_matrix[i_index, j_index]

                for i_sample in range(pdf.funcList().size()):
                    sample_yield = ROOT.RooProduct("tmp", "tmp", [pdf.funcList()[i_sample], pdf.coefList()[i_sample]])
                    yields[i_bin] += bin_width * sample_yield.getVal()
                    if i_sample in sample_values:
                        sample_values[i_sample] += [bin_width * sample_yield.getVal()]
                    else:
                        sample_values[i_sample] = [bin_width * sample_yield.getVal()]

                    for j_parameter, par in enumerate(required_params):
                        
                        cen_val = par.getVal()
                        par_err = par.getError()

                        par.setVal(cen_val + par_err)
                        par_upper_variation = sample_yield.getVal(observables)
                        par.setVal(cen_val - par_err)
                        par_bottom_variation = sample_yield.getVal(observables)


                        par.setVal(cen_val)
                        sample_yield.getVal(observables)

                        fill_array[i_sample, j_parameter] = (par_upper_variation - par_bottom_variation) / 2

                total_uncertanties_per_variation = np.sum(fill_array, axis=0)

                total_uncertainty = np.dot(total_uncertanties_per_variation, np.dot(numpy_corr_matrix, total_uncertanties_per_variation))

                yields_uncert[i_bin] = np.sqrt(total_uncertainty) * bin_width
                
        return yields, yields_uncert, sample_values
    
    def get_yields_no_fit(self, variable, observables, pdf, result):
        """Also get uncertainties if you pass a fit result.
        """
        yields = np.zeros(variable.numBins())
        yields_uncert = np.zeros(variable.numBins())
        sample_values = {}
        for i_bin in range(variable.numBins()):
            variable.setBin(i_bin)
            bin_width = variable.getBinWidth(i_bin)

            all_pdf_params = pdf.getParameters(observables)

            fit_result_params = result.floatParsFinal()

            required_params = all_pdf_params

            fill_array = np.zeros((pdf.funcList().size(), len(required_params)))

            number_of_parameters = len(required_params)

            parameters_indices_map = {}
            parameters_for_cov_matrix = []

            internal_index = 0
            for i_index, par in enumerate(required_params):
                par_index = fit_result_params.index(fit_result_params.find(par.GetName()))
                if (par_index == -1): continue
                parameters_indices_map[i_index] = internal_index
                internal_index += 1
                parameter_to_be_copied = fit_result_params[par_index]
                par.setVal(parameter_to_be_copied.getVal())
                par.setError(parameter_to_be_copied.getError())
                parameters_for_cov_matrix += [parameter_to_be_copied]

            corr_matrix = result.reducedCovarianceMatrix(parameters_for_cov_matrix)

            numpy_corr_matrix = np.zeros((number_of_parameters, number_of_parameters))

            for i_index in range(number_of_parameters):
                for j_index in range(i_index, number_of_parameters):
                    if i_index in parameters_indices_map and j_index in parameters_indices_map:
                        numpy_corr_matrix[i_index, j_index] = corr_matrix[parameters_indices_map[i_index], parameters_indices_map[j_index]] / \
                            np.sqrt(corr_matrix[parameters_indices_map[i_index], parameters_indices_map[i_index]] * corr_matrix[parameters_indices_map[j_index], parameters_indices_map[j_index]])
                        numpy_corr_matrix[j_index, i_index] = numpy_corr_matrix[i_index, j_index]
                    else:
                        numpy_corr_matrix[i_index, j_index] = int(i_index == j_index)

            for i_sample in range(pdf.funcList().size()):
                sample_yield = ROOT.RooProduct("tmp", "tmp", [pdf.funcList()[i_sample], pdf.coefList()[i_sample]])
                yields[i_bin] += bin_width * sample_yield.getVal()
                if i_sample in sample_values:
                    sample_values[i_sample] += [bin_width * sample_yield.getVal()]
                else:
                    sample_values[i_sample] = [bin_width * sample_yield.getVal()]

                for j_parameter, par in enumerate(required_params):
                    
                    cen_val = par.getVal()
                    par_err = par.getError()

                    par.setVal(cen_val + par_err)
                    par_upper_variation = sample_yield.getVal(observables)
                    par.setVal(cen_val - par_err)
                    par_bottom_variation = sample_yield.getVal(observables)


                    par.setVal(cen_val)
                    sample_yield.getVal(observables)

                    fill_array[i_sample, j_parameter] = (par_upper_variation - par_bottom_variation) / 2

            total_uncertanties_per_variation = np.sum(fill_array, axis=0)

            total_uncertainty = np.dot(total_uncertanties_per_variation, np.dot(numpy_corr_matrix, total_uncertanties_per_variation))

            yields_uncert[i_bin] = np.sqrt(total_uncertainty) * bin_width
                
        return yields, yields_uncert, sample_values
    
    def Draw(self, result, no_fit = False):
        self.boxes = []
        self.error_boxes = []
        self.error_boxes_prefit = []
        self.data_plots = []
        for channel in self.meas.GetChannels():
            channel_pdf = self.ws[str(channel.GetName()) + "_model"]
            obs_var = self.ws["obs_x_" + str(channel.GetName())]
            observables = ROOT.RooArgSet(obs_var)

            
            divide_value = 0.3
            self.cv_array += [ROOT.TCanvas("canvas" +  str(channel.GetName()), "canvas" +  str(channel.GetName()), 1500, 600)]
            pad1_upper = ROOT.TPad("pad1_upper" + str(channel.GetName()), "pad1_upper" + str(channel.GetName()), 0, divide_value, 0.5, 1)
            pad1_upper.Draw()
            pad1_bottom = ROOT.TPad("pad1_bottom" + str(channel.GetName()), "pad1_bottom" + str(channel.GetName()), 0, 0, 0.5, divide_value)
            pad1_bottom.Draw()

            pad2_upper = ROOT.TPad("pad2_upper" + str(channel.GetName()), "pad2_upper" + str(channel.GetName()), 0.5, divide_value, 1, 1)
            pad2_upper.Draw()
            pad2_bottom = ROOT.TPad("pad2_bottom" + str(channel.GetName()), "pad2_bottom" + str(channel.GetName()), 0.5, 0, 1, divide_value)    
            pad2_bottom.Draw()
            pad1_upper.cd()



            prefit_yields, prefit_unc, prefit_sample_values = self.get_yields(obs_var, observables, channel_pdf, result, True)

            self.hs_stacks += [ROOT.THStack("hs" + str(channel.GetName()), "hs" + str(channel.GetName()))]
            sample_histograms = []

            original_sample_bin_values = [0] * channel.GetData().GetHisto().GetNbinsX()

            for i, sample in enumerate(channel.GetSamples()):
                hist = sample.GetHisto()
                hist.SetFillColor(self.predefined_colors[i])
                hist.SetLineColor(ROOT.kBlack)
                sample_histograms += [hist]

                for i in range(hist.GetNbinsX()):
                    original_sample_bin_values[i] += hist.GetBinContent(i + 1)

                self.hs_stacks[-1].Add(sample_histograms[-1])

            channel_name = "_".join(channel.GetName().split("_")[1:])
            self.hs_stacks[-1].SetTitle(channel_name + " PREFIT")
            self.hs_stacks[-1].Draw("hist")

            maximum_y_val = ROOT.gPad.GetUymax()

            bin_index = 1

            for bin_index in range(1, sample_histograms[-1].GetNbinsX() + 1):
                leftEdge = sample_histograms[-1].GetBinLowEdge(bin_index)
                binWidth = sample_histograms[-1].GetBinWidth(bin_index)
                rightEdge = leftEdge + binWidth

                central_value = original_sample_bin_values[bin_index - 1]
                unc = prefit_unc[bin_index - 1]
                down_value = central_value - unc
                up_value = central_value + unc

                if up_value > maximum_y_val:
                    maximum_y_val = up_value

                self.boxes += [ROOT.TBox(leftEdge, down_value, rightEdge, up_value)]
                self.boxes[-1].SetFillStyle(3004)
                self.boxes[-1].SetFillColor(ROOT.kGray + 3)
                self.boxes[-1].Draw("same")

            self.hs_stacks[-1].SetMaximum(1.1 * maximum_y_val)

            self.data_histogram = channel.GetData().GetHisto()
            self.data_histogram.SetStats(0)
            self.data_histogram.SetMarkerStyle(3)
            self.data_histogram.SetMarkerSize(0.5)

            self.data_plots += [ROOT.TGraph()]
            for i in range(self.data_histogram.GetNbinsX()):
                self.data_plots[-1].SetPoint(i, self.data_histogram.GetBinCenter(i + 1), self.data_histogram.GetBinContent(i + 1))

            self.data_plots[-1].SetMarkerStyle(8)
            self.data_plots[-1].SetMarkerSize(0.5)
            self.data_plots[-1].Draw("same p")


            pad1_bottom.cd()
            number_of_bins = self.data_histogram.GetNbinsX()
            self.bias_graphs += [ROOT.TGraph(number_of_bins)]
            self.bias_graphs[-1].SetTitle("")
            self.bias_graphs[-1].SetMarkerSize(0.4)
            self.bias_graphs[-1].SetMarkerStyle(8)
            
            for i in range(1, number_of_bins + 1):
                original_value = original_sample_bin_values[i - 1]
                data_value = self.data_histogram.GetBinContent(i)
                unc = prefit_unc[i - 1]
                up_value = 1 + unc / (original_value)
                down_value = 1 - unc / (original_value)

                if down_value < 0.5:
                    down_value = 0.5
                if up_value > 1.5:
                    up_value = 1.5

                leftEdge = self.data_histogram.GetBinLowEdge(i)
                binWidth = self.data_histogram.GetBinWidth(i)
                rightEdge = leftEdge + binWidth

                self.bias_graphs[-1].SetPoint(i - 1, self.data_histogram.GetBinCenter(i),  data_value / original_value)
                self.error_boxes_prefit += [ROOT.TBox(leftEdge, down_value, rightEdge, up_value)]

                # print(leftEdge, down_value, rightEdge, up_value)

                self.error_boxes_prefit[-1].SetFillStyle(3004)
                self.error_boxes_prefit[-1].SetFillColor(ROOT.kGray + 3)
                self.error_boxes_prefit[-1].Draw("same")

            self.bias_graphs[-1].Draw("AP")

            # for box in self.error_boxes_prefit:
                # box.Draw("same")

            minimal_bin_value = self.data_histogram.GetBinLowEdge(1)
            maximum_bin_value = self.data_histogram.GetBinLowEdge(number_of_bins) + self.data_histogram.GetBinWidth(number_of_bins)


            self.bias_graphs[-1].GetYaxis().SetRangeUser(0.5, 1.5)
            self.bias_graphs[-1].GetXaxis().SetRangeUser(minimal_bin_value, maximum_bin_value)

            self.normal_lines += [ROOT.TLine(minimal_bin_value, 1, maximum_bin_value, 1)]
            self.normal_lines[-1].SetLineStyle(2)
            self.normal_lines[-1].SetLineWidth(1)
            self.normal_lines[-1].Draw("same")

            pad2_upper.cd()

            # bookkeep prefit yields

            if no_fit:
                postfit_yields, postfit_yields_uncert, postfit_sample_values = self.get_yields_no_fit(obs_var, observables, channel_pdf, result)
            else:
                postfit_yields, postfit_yields_uncert, postfit_sample_values = self.get_yields(obs_var, observables, channel_pdf, result, False)

            self.second_hs_stacks += [ROOT.THStack("fitted stack", "fitted stack")]

            # temp_histos
            color_number = 0

            for postfit in postfit_sample_values:
                temp_histo = ROOT.TH1F(channel_name + "afterfit" + str(postfit), channel_name + "afterfit" + str(postfit), len(postfit_yields), minimal_bin_value, maximum_bin_value)
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
            self.second_hs_stacks[-1].SetMaximum(1.1 * maximum_y_val)

            bin_index = 1
            for bin_index in range(1, temp_histo.GetNbinsX() + 1):
                leftEdge = temp_histo.GetBinLowEdge(bin_index)
                binWidth = temp_histo.GetBinWidth(bin_index)
                rightEdge = leftEdge + binWidth

                central_value = postfit_yields[bin_index - 1]
                unc = postfit_yields_uncert[bin_index - 1]
                down_value = central_value - unc
                up_value = central_value + unc

                self.boxes += [ROOT.TBox(leftEdge, down_value, rightEdge, up_value)]
                self.boxes[-1].SetFillStyle(3004)
                self.boxes[-1].SetFillColor(ROOT.kGray + 3)
                self.boxes[-1].Draw("same")


            self.data_plots[-1].Draw("same p")

            self.bias_second_graphs += [ROOT.TGraph(number_of_bins)]
            self.bias_second_graphs[-1].SetTitle("")
            self.bias_second_graphs[-1].SetMarkerSize(0.4)
            self.bias_second_graphs[-1].SetMarkerStyle(8)



            pad2_bottom.cd()

            for i in range(1, number_of_bins + 1):
                original_value = postfit_yields[i - 1]
                data_value = self.data_histogram.GetBinContent(i)
                unc = postfit_yields_uncert[i - 1]
                up_value = 1 + unc / (original_value)
                down_value = 1 - unc / (original_value)

                leftEdge = self.data_histogram.GetBinLowEdge(i)
                binWidth = self.data_histogram.GetBinWidth(i)
                rightEdge = leftEdge + binWidth

                self.bias_second_graphs[-1].SetPoint(i - 1, self.data_histogram.GetBinCenter(i),  data_value / original_value)
                self.error_boxes += [ROOT.TBox(leftEdge, down_value, rightEdge, up_value)]
                self.error_boxes[-1].SetFillStyle(3004)
                self.error_boxes[-1].SetFillColor(ROOT.kGray + 3)
                # self.error_boxes[-1].Draw("same")
            

            self.bias_second_graphs[-1].Draw("AP")

            for box in self.error_boxes:
                box.Draw("same")
            
            minimal_bin_value = self.data_histogram.GetBinLowEdge(1)
            maximum_bin_value = self.data_histogram.GetBinLowEdge(number_of_bins) + self.data_histogram.GetBinWidth(number_of_bins)


            self.bias_second_graphs[-1].GetYaxis().SetRangeUser(0.5, 1.5)
            self.bias_second_graphs[-1].GetXaxis().SetRangeUser(minimal_bin_value, maximum_bin_value)

            self.normal_lines += [ROOT.TLine(minimal_bin_value, 1, maximum_bin_value, 1)]
            self.normal_lines[-1].SetLineStyle(2)
            self.normal_lines[-1].SetLineWidth(1)
            self.normal_lines[-1].Draw("same")

            self.cv_array[-1].SaveAs(channel.GetName() + "_histo.png")
            self.cv_array[-1].Draw()
