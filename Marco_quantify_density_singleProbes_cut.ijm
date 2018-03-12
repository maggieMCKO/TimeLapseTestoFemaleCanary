/* 
1. For all probes, make ROIs on 10x images 
2. save ROIs in 10x folder
6. then run the following 
	aim: keep only the selected pre-defined ROIs
		open images as sequences
*/


// set scale 10x
run("Set Scale...", "distance=1 known=0.5396815708705357 pixel=1 unit=um");

input_Project = "F:\\Gahr\\EXPPIX\\RNAscope\\TestoProject\\";
probe = "UTS2B";
print(probe);

input = input_Project+probe+"\\10x\\";
output = input_Project+probe+"\\10x_cut\\";
File.makeDirectory(output);    

image_list = getFileList(input);
temp_image = input+image_list[1] ;
print(temp_image);
run("Image Sequence...", "open="+temp_image+" convert_to_rgb sort"); 

// open ROIs 
roi_path = input+File.separator+"RoiSet.zip";
roiManager("Open", roi_path); //the same as: open(roi_path); 

// get no. of ROIs
n_roi = roiManager("count");
print("n_roi = " + n_roi);

// get no. of slices
Stack.getDimensions(width, height, channels, slices, frames);
print("n_slices = " + slices);

// get mainTitle for selecting the window
mainTitle = getTitle();
print(mainTitle);

// looping through slices
for ( j = 1 ; j <= slices; j++) {
    selectWindow(mainTitle);
    setSlice(j);
    
    // setting output name
    Slice_name = getInfo("slice.label");
    dirCropOutput=output+File.separator+Slice_name ;
    print("dirCropOutput = "+dirCropOutput);

    // current slice to a new window
    run("Duplicate...", "title=crop");

    // regular expression on ROI name, becasue ROI are not necessary in order
    // source of this for-loop: http://forum.imagej.net/t/selecting-roi-based-on-name/3809/2
        for (u = 0; u< n_roi ; u++) { 
        	print("u ="+u);
        	roiManager("Select", u); 
        	rName = Roi.getName(); 
        	print("rName: "+ rName);
        	target = "000*"+j+"-.*" ;
        	print("target: "+ target);
        	//print(target);
        	
        	if (matches(rName, target))  { 
        		ind = u; 	
        		print("ind = "+ind);	
        		print(rName);
        		break
        		}
        	}

        // select the correct ROI
        roiManager("Select", ind);
        run("Make Inverse"); // inverse selection
        run("Cut"); // cut away non-HVC
        saveAs("Tiff", dirCropOutput);
        close(); // close duplication window
        //selectWindow(mainTitle);
        //close();
        
   
}
run("Close");
//close("roiManager");
//roiManager.close();
//roiManager("delete");
selectWindow("ROI Manager"); 
run("Close"); 