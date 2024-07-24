import ROOT
import math

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

        c1.Draw();
        c1.SaveAs(filename);