//define DAPI channel
Ch_DAPI = 1

//define I/O directories
dir = getDirectory("Choose Main Analysis Directory");
dir_layers_cropped = dir+"/S02_layers_cropped/";
dir_masks = dir+"/S03_DAPI_masks/";

print("Main directory: ", dir);
print("Image directory: ", dir_layers_cropped);
print("Mask directory: ", dir_masks);

dir_results = dir+"/S04_FIJI_results_nucl_intens/";
File.makeDirectory(dir_results); 


setBatchMode("hide");

//get image file list
image_list = getFileList(dir_layers_cropped);

//open threshold tool
run("Threshold...");

//loop for processing of all images in input folder

for(i = 0; i<image_list.length; i++) {
	print(image_list[i]);
	//open DAPI mask and generate ROI
	open(dir_masks+"mask_DAPI_"+image_list[i]);

	run("Analyze Particles...", "size=0-Infinity clear include summarize add slice");

	//open original image, measure all slices (multi-measure), save in results directory
	open(dir_layers_cropped+image_list[i]);
	roiManager("multi-measure measure_all");
	saveAs("Results", dir_results+image_list[i]+"_results.csv");
	run("Close All");

};

