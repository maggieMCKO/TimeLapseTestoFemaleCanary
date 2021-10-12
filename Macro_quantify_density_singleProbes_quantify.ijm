/* 
1. For all probes, make ROIs on 10x images 
2. save ROIs in 10x folder
6. then run the following 
	aim: keep only the selected pre-defined ROIs
		open images as sequences
*/

input_Project = "F:\\Gahr\\EXPPIX\\RNAscope\\TestoProject\\";
probe = "UTS2B";
print(probe);

input = input_Project+probe+"\\10x_cut\\";
output = input_Project+probe+"\\10x_cut\\";
//File.makeDirectory(output);    

image_list = getFileList(input);

//for (j = 0; j < 1; j++){
for (j = 0; j < image_list.length; j++){
	temp_image = input+image_list[j] ;
	print(temp_image);
	open(temp_image);
	// set scale 10x
	run("Set Scale...", "distance=1 known=0.5396815708705357 pixel=1 unit=um");

	// Color Thresholder 2.0.0-rc-65/1.51s
	// Autogenerated macro, single images only!
	min=newArray(3);
	max=newArray(3);
	filter=newArray(3);
	a=getTitle();
	run("RGB Stack");
	run("Convert Stack to Images");
	selectWindow("Red");
	rename("0");
	selectWindow("Green");
	rename("1");
	selectWindow("Blue");
	rename("2");
	min[0]=0;
	max[0]=47;
	filter[0]="pass";
	min[1]=0;
	max[1]=165;
	filter[1]="pass";
	min[2]=0;
	max[2]=160;
	filter[2]="pass";
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  setThreshold(min[i], max[i]);
	  run("Convert to Mask");
	  if (filter[i]=="stop")  run("Invert");
	}
	imageCalculator("AND create", "0","1");
	imageCalculator("AND create", "Result of 0","2");
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  close();
	}
	selectWindow("Result of 0");
	close();
	selectWindow("Result of Result of 0");
	rename(a);
	// Colour Thresholding-------------

	// analyze Particles
	if (j == 0){
		run("Analyze Particles...", "size=8-Infinity clear display add");
		}else{
			run("Analyze Particles...", "size=8-Infinity display add");}

	// get no. of ROIs
	n_roi = roiManager("count");
	print("n_roi = " + n_roi);

	//run("Close");
	//close("roiManager");
	//roiManager.close();
	//roiManager("delete");
	//selectWindow("ROI Manager"); 
	//run("Close"); 
}

/*
run("Close");
//close("roiManager");
//roiManager.close();
//roiManager("delete");
selectWindow("ROI Manager"); 
run("Close"); 
*/
