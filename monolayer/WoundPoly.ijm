// === Macro B: Extract wound edge polylines (universal Fiji version) ===
// 1. Open your wound stack.
// 2. Draw an ROI around the wound region (optional).
// 3. Run macro.

run("Set Measurements...", "area perimeter feret's redirect=None decimal=3");

// Get stack dimensions
getDimensions(width, height, channels, slices, frames);

// Use slices as number of time points
print("Total slices detected: " + slices);

roiManager("reset");

for (i = 1; i <= slices; i++) {

    setSlice(i);

    // Work on temporary copy
    run("Duplicate...", "title=tmpSlice");
    run("8-bit");

    // Change threshold method as needed (Triangle is often good for wounds)
    run("Auto Threshold", "method=Triangle white");
    run("Convert to Mask");

    // Find wound as largest particle
    run("Analyze Particles...", "size=500-Infinity show=Nothing clear add");
    count = roiManager("count");

    if (count > 0) {
        maxArea = -1;
        maxIdx = -1;
        for (r = 0; r < count; r++) {
            roiManager("select", r);
            run("Measure");
            aVal = getResult("Area", nResults - 1);
            if (aVal > maxArea) {
                maxArea = aVal;
                maxIdx = r;
            }
        }

        // Select only largest wound ROI
        roiManager("select", maxIdx);

        // Convert area ROI to boundary polyline
        run("Create Selection");

        // OPTIONAL: smooth boundary — replace this line with nothing if unwanted
        // run("Fit Spline Interpolator");

        // Add boundary to ROI manager
		roiManager("Add");
		// (Rename removed to avoid syntax / hidden-character issues)
    }

    roiManager("reset");
    close("tmpSlice");
}

print("Done. Edges stored in ROI Manager.");