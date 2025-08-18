//define DAPI channel
Ch_DAPI = 1
size_cutoff = 10

//define I/O directories
dir = getDirectory("Choose Analysis Main Directory");

dir_layers_cropped = dir+"/S02_layers_cropped/";
dir_masks = dir+"/S03_DAPI_masks/";
File.makeDirectory(dir_masks); 

print("Main directory: ", dir);
print("Image directory: ", dir_layers_cropped);
print("Mask directory: ", dir_masks);

//get image file list
image_list = getFileList(dir_layers_cropped);

//open threshold tool
run("Threshold...");

//loop for processing of all images in input folder

setBatchMode("hide");

for(i = 0; i<image_list.length; i++) {
	open(dir_layers_cropped+image_list[i]);
	print("Creating DAPI mask "+image_list[i]);

	//Set DAPI channel
	run("Split Channels");
	selectWindow("C"+Ch_DAPI+"-"+image_list[i]);
	//background subtraction and median filter
	run("Subtract Background...", "rolling=20 slice");
	run("Median...", "radius=1 slice");
	//set threshold
	setAutoThreshold("Huang dark");
	getThreshold(lower, upper);
	setThreshold(lower, upper);
	setOption("BlackBackground", false);

	//convert to mask, watershed, identify nuclei (remove small particles (noise)
	run("Convert to Mask", "method=Default background=Light calculate only");
	run("Watershed", "slice");
	run("Erode");
	run("Erode");
	run("Dilate");
	run("Dilate");
	run("Analyze Particles...", "size="+size_cutoff+"-Infinity show=Masks clear include summarize add slice");

	saveAs("Tiff", dir_masks+"/mask_DAPI_"+image_list[i]);
	//close("C"+Ch_DAPI+"-"+image_list[i]);
	
	//re-load DAPI image, overlay mask with original DAPI channel
	//open(maxp_dir+image_list[i]);
	//run("Split Channels");
	//selectWindow("C"+Ch_DAPI+"-"+image_list[i]);
	run("Merge Channels...", "c2='mask_DAPI_"+image_list[i]+"'] c6=[C"+Ch_DAPI+"-"+image_list[i]+"] create keep");
	run("Stack to RGB");
	saveAs("Tiff", dir_masks+"/DAPI image mask overlay_"+image_list[i]);
	
	

 	run("Close All");

};

