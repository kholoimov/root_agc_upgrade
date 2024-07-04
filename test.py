import ROOT

# Open the ROOT file
file = ROOT.TFile("data/histograms.root", "READ")

# Get the list of keys in the file
keys = file.GetListOfKeys()

# Loop over all keys
for key in keys:
    # Get the object associated with the key
    obj = key.ReadObj()

    # Check if the object is a histogram
    if obj.InheritsFrom("TH1"):
        print(f"Found histogram: {obj.GetName()}")
        # Do something with the histogram
        # e.g., obj.Draw()
    else:
        print(f"Skipping non-histogram object: {obj.GetName()}")

# Close the file
file.Close()