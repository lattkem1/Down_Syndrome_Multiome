// define directories
dir = getDirectory("Choose Main Directory");

dir_images_complete = dir+"/S01_greyscale_images_complete/";
dir_images_ROIs = dir+"/S02_images_with_layer_ROIs/";
dir_ROIs = dir+"/S02_layer_ROIs/";
dir_layers_cropped = dir+"/S02_layers_cropped/";

print("Main directory: ", dir);
print("Complete Images directory: ", dir_images_complete);
print("Complete Images with ROIs will be saved in: ", dir_images_ROIs);
print("ROIs will be saved in: ", dir_ROIs);
print("Cropped images will be saved in: ", dir_layers_cropped);

File.makeDirectory(dir_images_ROIs); 
File.makeDirectory(dir_ROIs); 
File.makeDirectory(dir_layers_cropped); 

//get image file list
image_list = getFileList(dir_images_complete);

n_images = image_list.length;


do {
        
        waitForUser("Open image","Open image, then click OK");

        title = getTitle;
        
        layer_set_idx = 0;
        
        //let user select 6 layer bins, save, repeat until user moves to next image
        
        do {
        	
        	layer_set_idx = layer_set_idx+1;
        	
        	makeRotatedRectangle(300, 300, 800, 800, 300);//makeRectangle(x1, y1, x2, y2, width)
			waitForUser("Set ROIs","Set 6 ROIs from layer 1 to 6, then click OK");
        	
        	
	        //-- Create and show a blocking dialog box
			Dialog.create("Choice");
			Dialog.addChoice("Do you want to save this set of ROIs?", newArray("Yes","No"));
			Dialog.show();
			//-- When OK is clicked, read out the status of the dialog
			save_ROIs=Dialog.getChoice();
			
			if (save_ROIs == "Yes") {
				
				for (j=1;j<=RoiManager.size;j++) {
					
					roiManager("Select", j-1);
					run("Duplicate...", "duplicate");
					ROI_title = replace(title,".tif","_set_"+layer_set_idx+"_layer_"+j);
					print("Save ", ROI_title);
	        		saveAs("tiff", dir_layers_cropped+ROI_title);
	        		roiManager("Save", dir_ROIs+ROI_title+".zip");
	        		close();
	        		selectImage(1);
	        		roiManager("Set Line Width", 3);
					roiManager("Draw");
					
				}
				}
		
		roiManager("Deselect");
		roiManager("Delete");
		
        //-- Create and show a blocking dialog box
		Dialog.create("Choice");
		Dialog.addChoice("Do you want to add another set of ROIs?", newArray("Yes","No"));
		Dialog.show();
		//-- When OK is clicked, read out the status of the dialog
		finished=Dialog.getChoice();
        		
        } while (finished=="Yes")
        
	saveAs("tiff", dir_images_ROIs+title);
	close();
	
	//-- Create and show a blocking dialog box
	Dialog.create("Choice");
	Dialog.addChoice("Do you want to analyse another image?", newArray("Yes","No"));
	Dialog.show();
	//-- When OK is clicked, read out the status of the dialog
	finished=Dialog.getChoice();
        		
    } while (finished=="Yes")









