dir1=getDirectory("Select the folder that contains your TIF stacks to analyse"); 
print(dir1); 

// Creates a directory called /raw_values on the same level as the selected directory
// and sets it as the destination for output
dir2= dir1+ "/.." + "/binarised_raw_values/"; 
print(dir2); 
File.makeDirectory(dir2); 
selectWindow("Log");
run("Close")

//Generates a list of those files
list = getFileList(dir1);
Array.sort(list);
 
//loop the macro for all files in the folder
for(i=0; i<list.length; i++) {
filename = dir1 + list[i];

//makes sure the file in the folder it is accessing is a Tiff file
if(endsWith(filename, ".tiff")){

open(filename);

//Stores the initial image name into the nameStore variable
nameStore = getTitle();
width=(getWidth()/2);
run("Radial Profile Angle", "x_center=1024 y_center=1024 radius=" + width + " starting_angle=0 integration_angle=180 calculate_radial_profile_on_stack");
					run("Clear Results");
					for(j = 0; j != Ext.getStackSize; j++)
					{           
						for(k = 0; k != Ext.getBinSize; k++)
						{
							setResult("Slice", j * Ext.getBinSize + k, j + 1);
							setResult("X"    , j * Ext.getBinSize + k, Ext.getXValue(k));
							setResult("Y"    , j * Ext.getBinSize + k, Ext.getYValue(j, k));
						}
					}	
					// Saves results using the input save directory and closes all
saveAs("Results", dir2+nameStore+"_binary.csv");
run("Close All");
}
}
selectWindow("Results");
run("Close")
