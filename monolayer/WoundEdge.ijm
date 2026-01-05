// === Macro A: Wound area vs time ===
// 1. Open your wound stack (TIF/AVI/etc.)
// 2. Draw a rough rectangular ROI around the wound
// 3. Run this macro

run("Set Measurements...", "area mean min centroid redirect=None decimal=3");

// Get stack dimensions
getDimensions(width, height, channels, slices, frames);
// We will use "slices" directly as the number of time points

print("Total slices: " + slices);

// Duplicate only the ROI area to speed things up
roiManager("reset");
run("Duplicate...", "title=WoundROI");

selectWindow("WoundROI");

// Preview threshold on slice 1
Stack.setSlice(1);
run("Enhance Contrast...", "saturated=0.35");
run("Auto Threshold", "method=Otsu white");
waitForUser("Check threshold.\nIf the wound is not highlighted, modify 'Auto Threshold' settings.\nClick OK to process all frames.");

run("Clear Results");

for (i = 1; i <= slices; i++) {
    Stack.setSlice(i);

    // Work on temporary copy
    run("Duplicate...", "title=tmpSlice");
    run("8-bit");
    run("Auto Threshold", "method=Mean white");

    // Analyze particles — keep only the largest region = wound
    run("Analyze Particles...", "size=500-Infinity show=Nothing clear add");

    n = roiManager("count");

    if (n > 0) {
        maxA = -1;
        maxIndex = -1;

        for (r = 0; r < n; r++) {
            roiManager("select", r);
            run("Measure");
            aVal = getResult("Area", nResults - 1);

            if (aVal > maxA) {
                maxA = aVal;
                maxIndex = r;
            }
        }

        // Now measure only largest ROI as wound
        roiManager("select", maxIndex);
        run("Measure");
        row = nResults - 1;
        setResult("Slice", row, i);
    }

    roiManager("reset");
    close("tmpSlice");
}

// Save results
savePath = getDirectory("Choose a folder to save wound area CSV");
if (savePath != null) {
    saveAs("Results", savePath + "wound_area_vs_time.csv");
}