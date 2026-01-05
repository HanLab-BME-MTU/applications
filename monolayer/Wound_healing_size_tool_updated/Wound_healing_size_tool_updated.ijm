//Variable declaration
var filter_radius = 20; 
var threshold = 100;
var sat_pix = 0.001;
var width = getWidth;
var height = getHeight;
var min_area = 100;
var ort = "No";
var scale = "No";
var diagonal= "No";

// ----------------------------------------------------------------------------------------------------------------------
// Auxiliary functions
// ----------------------------------------------------------------------------------------------------------------------
function edge_min_coordinates(arr_ROI) 
{
   arr_min = newArray();
   arr_min = Array.concat(arr_min, arr_ROI [0]);
   for (x = 1; x< arr_ROI.length - 2; x++)
   {
         if (arr_ROI [x-1] > arr_ROI [x] && arr_ROI [x] < arr_ROI [x+1])
         {
            arr_min = Array.concat(arr_min, arr_ROI [x]);
         }
   }
    
   return arr_min;   

}


// ----------------------------------------------------------------------------------------------------------------------
function edge_max_coordinates(arr_ROI) 
{
   arr_max = newArray();
   for (x = 1; x< arr_ROI.length - 2; x++)
   {
         if (arr_ROI [x-1] < arr_ROI [x] && arr_ROI [x] > arr_ROI [x+1])
         {
            arr_max = Array.concat(arr_max, arr_ROI [x]);
         }
   }
   
   arr_max = Array.concat(arr_max, arr_ROI [arr_ROI.length - 1]);
   
   return arr_max;   
}


// ----------------------------------------------------------------------------------------------------------------------
function diff_arrays(arr_max, arr_min)
{
    arr_diff = newArray();

    // Use the shorter of the two arrays to avoid out-of-range indexing
    len = arr_max.length;
    if (arr_min.length < len)
        len = arr_min.length;

    for (x = 0; x < len; x++)
    {
        arr_diff = Array.concat(arr_diff, arr_max[x] - arr_min[x]);
    }
    return arr_diff;
}




//------------------------------------------------------------------------------------------------------------------------
function scale_fixed(scale)
{
	getPixelSize(unit, pw, ph);
   if(scale=="Yes")
	{
   run("Set Scale...", "distance=1 known="+pw+" unit=" + unit+" global");
	}
	if(scale=="No")
	{
   run("Set Scale...", "distance=1 known="+pw+" unit=" + unit);
	}

}


//------------------------------------------------------------------------------------------------------------------------
function diagonal_fixed(diagonal, angle)
{
	
   if(diagonal=="Yes")
	{
	angle1=(90-(angle-90));
	beta = angle1 * (PI/180);
   correct_angle=sin(beta);
	}
	if(diagonal=="No")
	{
   correct_angle=1;
	}
	return correct_angle;


}



// ----------------------------------------------------------------------------------------------------------------------
// Body of the macro
// ----------------------------------------------------------------------------------------------------------------------
macro 'Wound healing Size Stacks Tool - C000Da8Da9Db9DcbDddDecDedDeeC9abD31C68aD23Dc0CcddDbfC566DbbCbbbDffC89bD12D25D4eD55D61DbdDd7CeeeD37D70D71D72D73D74D75D77D78D7aD80D81D83D84D86D87D88D89D8aD8bD8cD8dDf0Df1Df2Df3C112DaaDcaCaacD5bC79aD2eD4dD53D60D6cDe1CdddD4cD67D7cD93Df4C666Db8CbcdD33C9abD1cD22D26D29D3eD46D62D7bD92Db4Db6Dc2Dd4DdaDe2De3C111DccC9acD16D19D28D2bD3cD6dD95Db5Dc8C78aD14D1dD32CdddD13D85D90Da1C579D1aDd0CbcdD35D45Dc3DcfDd3C89bD17Da6C444D97CabcD42D4aD6eDd5C89bD02D10D1bD3bD43D68Da0Dc6De4CddeD2fD52D64D79D8fDa5Df6C777Df7CccdD0fD1eD21D27D39D44D4fD5cD5fD6aD6fD7dD82Db3Dd8CaaaD9fC000De9C9abD15D47D49D6bD91Da3Dd6De5C78aD06D07D34D36D3aD40D4bD51D69D7eDa4DceDe6CdddD1fD24D2aD3fD63D76D7fD8eC469D0aCbccD3dD5eDbcDc1C89bD48D57D59D66Db0Db1Db7C223Da7CabcD04D2dD54D5aD65C79bDd2C68aD0cD38Da2Dd9C111D98D99D9aD9bD9cD9dD9eC67aD03D0bDc4Dc7De0C8abD5dC555Df8Df9DfaDfbDfcDfdDfeCbbcDdfC89aDc9CaabDafC9abD11CcddD41D58Dc5Df5C358D0dCbccD01D56D94DbeC122De7C68aD09D20D50C579D05D18C444DdeC000DbaDdcDe8DebC579D00D0eD30C334DaeCabcDd1CaaaDefC557DdbCbbcD2cD96C000DeaC469Db2C333DabDacC578DcdC123DadC57aD08'

{
   run("Select None");
   Width = newArray();
   StandarDeviation = newArray();
   Area= newArray();
   AreaFraction = newArray();
   height_total=getHeight();
   width_total=getWidth();

   run("Wound healing size tool options");
   run("Options...", " black");
   roiManager("reset");
   
   getPixelSize(unit, pw, ph);
   getPixelSize(unit, pw, ph);
   title = getTitle();

   for (t = 1; t <= nSlices; t++) {
       setSlice(t);

       // Work on a temporary copy of this slice only
       run("Duplicate...", "title=tmp");
       selectWindow("tmp");

       setForegroundColor(0, 0, 0);
       setBackgroundColor(255, 255, 255);
       roiManager("reset");
       roiManager("Associate", "false"); // we only care about this slice

       run("8-bit");
       run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
       run("Enhance Contrast...", "saturated="+sat_pix+" normalize");
       run("Variance...", "radius="+filter_radius); // NOTE: no 'stack' here
       setThreshold(0, threshold);
       run("Convert to Mask", "black");
       run("Fill Holes");
       run("Select All");

       // Analyze only this slice
       run("Analyze Particles...",
           "size="+min_area+"-Infinity circularity=0.00-1.00 show=Nothing add");

       count = roiManager("count");
       if (count == 0) {
           // No wound found on this frame
           print("No wound ROI found on slice " + t + " – skipping.");
           close("tmp");
           selectWindow(title);

           // Record zeros or NaNs as you prefer:
           Area           = Array.concat(Area, 0);
           AreaFraction   = Array.concat(AreaFraction, 0);
           Width          = Array.concat(Width, 0);
           StandarDeviation = Array.concat(StandarDeviation, 0);
           continue;
       }

       // If multiple ROIs, pick the largest as wound
       if (count > 1) {
           area_large = newArray(count);
           for (i = 0; i < count; i++) {
               roiManager("select", i);
               getStatistics(aTmp, mean, min, max, std, histogram);
               area_large[i] = aTmp;
           }
           largest = 0;
           large   = 0;
           for (i = 0; i < count; i++) {
               if (area_large[i] > largest) {
                   largest = area_large[i];
                   large = i;
               }
           }
           roiManager("select", large);
       } else {
           roiManager("select", 0);
       }

       // --- EDGE SHAPE / WIDTH ---
       Roi.getContainedPoints(xpoints, ypoints);
       min_x   = edge_min_coordinates(xpoints);
       max_x   = edge_max_coordinates(xpoints);
       diff_vec = diff_arrays(max_x, min_x);
       Array.getStatistics(diff_vec, low, high, avg_widthp, std_distp);

       // --- AREA of that wound on this slice ---
       run("Set Measurements...", "area redirect=None decimal=3");
       getStatistics(area2, mean, min, max, std, histogram);

       // Clean up temp window and go back to original
       close("tmp");
       selectWindow(title);

       // Convert to physical units and store
       total_area = (height_total * pw) * (width_total * pw);
       Area           = Array.concat(Area, area2);
       AreaFraction   = Array.concat(AreaFraction, (area2 / total_area) * 100.0);
       Width          = Array.concat(Width, (avg_widthp * pw));
       StandarDeviation = Array.concat(StandarDeviation, (std_distp * pw));
   }
   roiManager("Show None");
   Array.show("Results (row numbers)", Area,AreaFraction,Width,StandarDeviation);
   scale_fixed(scale);
   setTool("rectangle");
   
   

  
   
}
// ----------------------------------------------------------------------------------------------------------------------
// Body of the macro 2
// ----------------------------------------------------------------------------------------------------------------------
macro 'Wound healing size Tool - C358Dc2C9acD1cD47D59D6cD72Db6Db7Dc4C89bD01D0aD10D4aD58D64D70CdddD5bD62DaaDe1C68aD0fCbccD28D6aD92Dd9C9abD2aCeeeD38D45D5dD73D74D80D81D82D83D84D85D87D88D89D8aD8bD90D91D93D94D95D97D98D99D9aD9bD9cD9dD9eD9fDa7Da9DabDacDadDaeDc9DfeC67aD02D04D07D0bD30DefCabcD8fC9abD17D1dD2fD31D34D46D50D67D79D7eDa2Da5Db3Db9Dc3DceDd2De5Df6CdddD4dD8eDa8DafC79aD1eD3dDb4DbdDeeCccdD41De6DecDf2C9abD11D24D27D3eD69Da1DbfDd8DdfDebC579Dc6DdaDe8CabcD12D6fD8cDbeDdeDedDf5C8abD29D43DbaDbcDe3Df9CdddD21D3aD66D78D96Db1Dc5C78aD7aDb8CccdDf7C9abD15DfbC68aD09D7fDfcCbbcD51DbbDc1Df4CdeeDd3C89bD0cD18D1aD20D2eD32D5fD7dDcaCcddD13D25Da4Db5DfdC579D06D08D19D60Df0CaacD3bC89bD23D2cD5aD65D68D6bD7cDa6Dd1Dd5DeaDf1DffC78aD05D14D22D35D36D42D44D57D63D76Dc8Dd4Dd7De2DfaCbcdD1fD3fD54De7C68aD0dD0eD5eDb0Dc0DdcDddCabcD26D3cD48D6eDd6DdbDe9CddeD2dD75D7bD86D8dDa0Da3C79bD16D4eD5cDb2CcddD77Dc7DcbC57aD49Dd0CccdD33D37D4fD52D55D6dDccDcfC469D03D40De0CbccD2bD53DcdDf3C68aD4bC579D1bD39D4cD61CabcD56D71De4Df8C469D00'

{
   run("Select None");
   snapshot();
   setupUndo(); 
   run("Wound healing size tool options");
   run("Options...", " black");
   run("Duplicate...", "duplicate");
   setForegroundColor(0, 0, 0);
   setBackgroundColor(255, 255, 255);
   roiManager("reset");
   roiManager("Associate", "true");
   run("8-bit");
   

   if (isOpen("ROI Manager")) 
   {
      selectWindow("ROI Manager");
      run("Select None");
      run("Close");
   }
   getPixelSize(unit, pw, ph);
   run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
  
   run("Enhance Contrast...", "saturated="+sat_pix+" normalize");
   run("Variance...", "radius="+filter_radius+" stack"); 
   setThreshold(0, threshold);
   run("Convert to Mask", "black");
   run("Fill Holes");
   run("Select All");
   
   run("Analyze Particles...", 
          "size="+min_area+"-Infinity circularity=0.00-1.00 show=Nothing add stack");
   run("Revert");
 
   roiManager("Show None");
   close();
  
   if (roiManager("count")>1)
   {
       area_large = newArray(roiManager("count"));
       for (i = 0; i<(roiManager("count")); i++)
       {
           roiManager("select", i)
           getStatistics(area_large[i], mean, min, max, std, histogram);
       }
       largest = 0;
       for (i = 0; i<(roiManager("count")); i++)
       {
           if (area_large[i]>largest)
           {
               largest = area_large[i];
               large = i;
           }
       }
       roiManager("select", large);
   }
   else
   {
       roiManager("select", 0);
   }
   reset();
   setupUndo();
   roiManager("Set Color", "cyan");
   
   Roi.getContainedPoints(xpoints, ypoints);
   min_x = edge_min_coordinates(xpoints);
   max_x = edge_max_coordinates(xpoints);
   diff_vec = diff_arrays(max_x , min_x);

   Array.getStatistics(diff_vec, low, high, avg_width, std_dist);
    
   List.setMeasurements;  
   angle = List.getValue("Angle");
   correct_angle = diagonal_fixed(diagonal, angle);
   
   avg_width = (avg_width*pw)*correct_angle;
   std_dist = (std_dist*pw);

   run("Set Scale...", "distance=1 known="+pw+" unit=" + unit);
   run("Set Measurements...", "redirect=None decimal=3");
   getStatistics(area, mean, min, max, std, histogram);
   height_total=getHeight()*pw;
   width_total=getWidth()*pw;
   total_area=(height_total*width_total);
   area_fraction=(area/total_area)*100;
  
   n1=getValue("results.count");
   title_image=getTitle();
   setResult("Label", n1, title_image);
   setResult("Area "+unit+"^2", n1, area);
   setResult("Area %", n1, area_fraction);
   setResult("Width "+unit, n1, avg_width);
   setResult("Standard deviation "+unit, n1, std_dist);
   scale_fixed(scale);
   setTool("rectangle");
  
  
}

// ----------------------------------------------------------------------------------------------------------------------
// Body of the macro 3 manual
// ----------------------------------------------------------------------------------------------------------------------

macro 'Wound healing size Manual Tool - C89bD19D5dD6eDd2C579D03D0cCbccDe2C123DdaC9abD04D2bD42D49D65D7cC78aD00D02D23D48D4dD5eDb0Df0DfdCdddD4cD5aC111Dc8De9DeaC9abD11D17D1dD1fD22D25D38D4fD62D6bD95Da3Dc5De4De7C68aD0bD60CccdD8eC469D0dD0eDb2DfcCabcD01D27D3aD50D6dD7bD7dC89aDd6CeeeD80D90C001D8bDddC8abD6aDd3C679Df1DfbCccdD47C333D98D99D9fC9acD15D29D35D66C79aD34D3bD4aD4bD7eD92Dd4CeeeD37D52D5cD72D73D74D77D7aD81D82D83D84D85D86D93D94Dc3C222D9aD9bD9cD9dDabDacDadDbbDbcDbdDcbDccDcdC78aDf4CdddD41D67D71Da5Dd1De3DefC555D97CabcD40C89bD2dD30D3cD3dD46Da2Db4Dc2De0C000D88D89D8aD9eDb6Db7Db8Db9DbaC89bD12D1bD1cD2eD53D68D8fDb3De1DeeC579D18Df9CbcdD10D39D6cD6fDe8C233Da9C78aD43D55Dc1Dc4CddeD13D2cD78C122DaaDceC68aD0fD14D36D57D5bD61D69Da4Df3DfaCcddD24D54D56Dd5C569DdeDf8C79bD32D51Da0DffC111Dc7Dc9C8abD16D7fDe5De6C67aD06D07D09D59Dc0Df2CccdD21D33D4eD5fD64D79Da1Db1C345DcfCabcD20D76C68aD3fC789Dd7CbbcD1eD44D63C579D05D1aDf6CbccD28D2aD91Df5C122Da7CdddD3eD45D58CdddDedC111D8cDcaDdbC89bD31C344Dc6C68aD2fDd0DfeC778D96CabcD26DdfC57aD08C233DebCdeeD70D75C579Df7C445Da6C899Dd9C568Dd8C333DafDbfC444D8dC888D87Db5C223Da8C778DecC469D0a'

{
	
   roiManager("Add");
   roiManager("select", 0);
   reset();
   setupUndo();
   roiManager("Set Color", "cyan");
   getPixelSize(unit, pw, ph);
   run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
   
   Roi.getContainedPoints(xpoints, ypoints);
   min_x = edge_min_coordinates(xpoints);
   max_x = edge_max_coordinates(xpoints);
   diff_vec = diff_arrays(max_x , min_x);

   Array.getStatistics(diff_vec, low, high, avg_width, std_dist);
    
   List.setMeasurements;  
   angle = List.getValue("Angle");
   correct_angle = diagonal_fixed(diagonal, angle);
   
   avg_width = (avg_width*pw)*correct_angle;
   std_dist = (std_dist*pw);

   run("Set Scale...", "distance=1 known="+pw+" unit=" + unit);
   run("Set Measurements...", "redirect=None decimal=3");
   getStatistics(area, mean, min, max, std, histogram);
   height_total=getHeight()*pw;
   width_total=getWidth()*pw;
   total_area=(height_total*width_total);
   area_fraction=(area/total_area)*100;
  
   n1=getValue("results.count");
   title_image=getTitle();
   setResult("Label", n1, title_image);
   setResult("Area "+unit+"^2", n1, area);
   setResult("Area %", n1, area_fraction);
   setResult("Width "+unit, n1, avg_width);
   setResult("Standard deviation "+unit, n1, std_dist);
   scale_fixed(scale);
   setTool("rectangle");
   roiManager("Delete");
  
}



//
// ----------------------------------------------------------------------------------------------------------------------
// Tooltips & Dialogs

macro 'Rotate 90° Tool - C000D14D15D16D17D23D32D42D52D62D72CaacD19D6fDa8DaaC79bD55D6cDffCdddD70D7cDbcDdaDddDf8C469DecCbccD83C9abD74CeeeDb0Dc0C112D33CabcD2cD3dDbfDd6C89bD1bD1eD1fD2dD2fD36D44D4cD75D87D8cDb1Db4DdcCeeeD01D02D09D0aD0bD0cD0dD0eD0fD1aD2aD39D3aD46D49D4aD59D5aD69D6aD79D7aD88D89D8aD98D99D9aDa9Db5Db8Db9Dc8Dc9DcaDd8Dd9Df9C777D41D51CcccD60C9abD37D81D8bD93De1DedCeeeD00D10D20D30D40D50D80D90Da0Dd0De0Df0C000D22CabcD3cD8fDa7De2DfbC89bD3eD4bDa2Dc1Dc5De5CddeD3fD68D94De9C68aD6dCccdD45D47D7eD86DbbDc6Dd2De6DefDfeC9abD2bD48D67D7bD92D9fDabDadDbdDccDd1DdbDdeDe7Df1Df2C334D82CbbcD38D4eD6bD84D95D9eDa3DaeDb7C8abD1cDd5Df4C79aD34D65D85Dc4DcbDd3CdddDa6CaabD21CabcD4fD91DbaDc3DceDcfDd7DdfDebDf6C89bD35D7dDa5DafDb6De4DeaC68aD4dD56D6eD76D77D7fD96Df7DfdCbcdD5eDfcC333D71C888D31CcddD66D78D9cDd4Df5DfaC9acD57D5dD9bDa1Db2DbeC222D05D06D07CdeeD29C78aD54D8eDa4Dc2De3DeeC666D61CdddD5cDe8Df3CabbD64C57aDacCbccD11D3bD58Dc7C233D26D27CccdD8dD97C000D13C455D25C444D04D63C78aD5fDcdC567D28CbbbD03C579D2eCbccD5bDb3C333D24C344D53C334D18C666D08C68aD1dD9dC111D73C555D43C667D12'

{
	snapshot();
    run("Rotate 90 Degrees Right");
    setTool("rectangle");
    
   
}

macro 'Wound healing size tool options'
{

    set_scale_global = newArray("Yes","No");
    is_diagonal = newArray("Yes","No");
    Dialog.create("Wound healing size options");
    Dialog.addNumber("Variance window radius", filter_radius);
    Dialog.addSlider("Threshold value", 0, 255,  threshold);
    Dialog.addNumber("Percentage of saturated pixels", sat_pix);
    Dialog.addChoice("Set Scale global?", set_scale_global);
    Dialog.addChoice("The scratch is diagonal??", is_diagonal);

    Dialog.show();
    
    filter_radius = Dialog.getNumber();
    threshold = Dialog.getNumber();
    sat_pix = Dialog.getNumber();

    scale = Dialog.getChoice();
    diagonal = Dialog.getChoice();
    scale_fixed(scale);
  
}


// ----------------------------------------------------------------------------------------------------------------------
macro 'User defined background'
{
    Dialog.create("Background specification")
    Dialog.addMessage("ROI belonging to the gap in the wound must be white, select " +  
                                         "Invert' to turn it to white, 'Keep' to keep it as it is");
    arr_res = newArray("Keep", "Invert");
    Dialog.addChoice("Invert the colour?", arr_res);
    Dialog.show();
    res = Dialog.getChoice();
    if (res=="Invert")
        {
          run("Invert", "stack");
        }   
} 