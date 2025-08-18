
dir = getDirectory("Choose Main Directory (Images need to be opened)");
dir_images_complete = dir+"/S01_greyscale_images_complete/";
File.makeDirectory(dir_images_complete); 

// for all images, get initial file name and save as greyscale tif
n_images = nImages;
print("save N images: "+n_images);
for (i=1;i<=n_images;i++) {
        selectImage(1);
        Stack.setDisplayMode("grayscale");
        title = getTitle;
        title_split = split(title,".ndpis");
        title_new = title_split[0];
        //title=replace(title,"/","_"); //remove "/" if in title, else error saving
        print(title+" => save as "+title_new+".tif");
        saveAs("tiff", dir_images_complete+title_new);
        close();
}

