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


    