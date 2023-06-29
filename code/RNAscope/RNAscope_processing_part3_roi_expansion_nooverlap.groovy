//#@ImagePlus imp (label="Input Image")
//#@Integer expansion (label="Expand By (pixels)", value=5)
#@RoiManager rm
expansion = 5


// Start with binary mask
// Get Voronoi
// Get Mask of voronoi
// Get label from Voronoi
// Get Distance map from voronoi mask
// set distance threshold base on expansion and binarize
// normalize binary image
// multiply binarized distance map with voronoi label
// remove initial binary mask -> cytoplasm
// do not remove mask -> cell
// multiply original binary mask with voronoi label -> nuclei



import ij.*
import ij.gui.*
import ij.process.*
import inra.ijpb.binary.*
import inra.ijpb.morphology.strel.*
import inra.ijpb.morphology.*
import ij.plugin.ImageCalculator
import ij.plugin.frame.RoiManager
import groovy.io.FileType

//import ij.process.ImageConverter

Dir = "workingDirectory"
inputDir = Dir + "/inputFolder/"
outputDir = Dir + "/outputFolder/"
ext = ".tif"



//loop over folder
def listFiles = []

def dir = new File(inputDir)
dir.eachFileRecurse (FileType.FILES) { file ->
  listFiles << file
}


for (int i = 0; i < listFiles.size; i++){
	
	def image = listFiles[i].getName() //obtain image name
	print(System.getProperty("line.separator") + "image: " + image + System.getProperty("line.separator"))
	fileWithoutExt = listFiles[i].getName().split("\\.")[0] //obtain image name without extension
	print("image no ext: " + fileWithoutExt + System.getProperty("line.separator"))
	if(image.endsWith(ext)){ //process onl files with ext defined up
		print("Processing image..." + System.getProperty("line.separator"))
		imp = IJ.openImage(inputDir + image)
		
		ImageConverter.setDoScaling(true) //convert image to 8bit. if already 8bit skip this
		IJ.run(imp, "8-bit", "")
			
		def rm = RoiManager.getInstance() == null ? new RoiManager() : RoiManager.getInstance()
		rm.reset()

		def rois = CellMaker.getCells(imp, expansion)

/// 	ROI expension script
		rois.each{rm.addRoi(it)}

		rm.runCommand("Sort")
		rm.runCommand("Save", outputDir + fileWithoutExt + "_expanded.zip")
		
		IJ.save(imp, "Tiff")
		imp.close()
	
	}
}
/// end loop
print("Finished!");


return
class CellMaker {
		
			static List<Roi> getCells(ImagePlus mask, int expansion) {
				// Make sure it is a mask
				if ( !mask.getProcessor().isBinary() ) return null
				def voronoi = mask.duplicate()
				// Get Voronoi
				IJ.run(voronoi, "Voronoi", "")
				//voronoi.show()
				voronoi.getProcessor().setThreshold(1, 255, ImageProcessor.NO_LUT_UPDATE)
		
				voronoi.setProcessor(voronoi.getProcessor().createMask())
				voronoi.getProcessor().invert()
		
				// Label Voronoi and store
				def voronoi_labels = BinaryImages.componentsLabeling(voronoi, 8, 16)
				//voronoi_labels.show()
		
			
				// Dilate by circular structuring element
				def image = mask.getProcessor()
		 		def se = DiskStrel.fromRadius(expansion)
				def cells = new ImagePlus("Cells", Morphology.dilation(image, se))
		
				// Assign same labels as voronoi with multiplication
				def cells_labels = assignLabels(cells, voronoi_labels)
				//cells_labels.show()
				// Assign labels to initial mask too
				def nuclear_labels = assignLabels(mask, voronoi_labels)
				//nuclear_labels.show()
				// Make rois
				def nuclear_rois = labelsToRois(nuclear_labels, "Nucleus")
				def cells_rois   = labelsToRois(cells_labels, "Cell")
		
				// Cytoplasm needs special Shape ROIs
				def padding = Math.ceil(Math.log10(nuclear_rois.size()))
				def cytoplasm_rois = nuclear_rois.collect{ id, nuclear_roi ->
					def cell_roi = cells_rois.get(id)
					
					def cyto_roi = new ShapeRoi(cell_roi)
					cyto_roi = cyto_roi.not(new ShapeRoi(nuclear_roi))
					cyto_roi.setProperty("ID", id as String)
					cyto_roi.setProperty("Class", "Cytoplasm")
					cyto_roi.setName("Cytoplasm #"+IJ.pad(id as int,padding as int))
		
					return cyto_roi
				}
				
				def rois = []
		
				rois.addAll(nuclear_rois.collect{ key, roi -> roi })
				rois.addAll(cells_rois.collect{ key, roi -> roi })
				rois.addAll(cytoplasm_rois)
				
				return rois		
			}
		
			static ImagePlus assignLabels(ImagePlus image, ImagePlus labels) {
				def copy = image.duplicate()
				copy.getProcessor().subtract(254)
				ImageCalculator ic = new ImageCalculator()
				return ic.run("Multiply create", labels, copy)
			}
			
		
			
			static Map<Integer, Roi> labelsToRois(ImagePlus labels, String prefix) {
			    def ip = labels.getProcessor()
			    def wand = new Wand(ip)
			    //def ov = new Overlay()
		
				// Go through the image and flood fill everything
				def width = labels.getWidth()
				def height = labels.getHeight()
		
				// Get all the labels
				def labelIds = LabelImages.findAllLabels(labels) as List
				def padding = Math.ceil(Math.log10(labelIds.size()))
				println(labelIds)
				def doneLabels = []
				def rois = [:]
				for (int i=0; i<width; i++) {
					for(int j=0; j<height; j++) {
						def val = ip.getf(i,j)
						if (val != 0 && !doneLabels.contains(val)) {
							doneLabels.add(val)
							def id = labelIds.findIndexOf{it == val}
							//ov.add(new PointRoi(i,j))
							wand.autoOutline(i, j)
							def roi = new PolygonRoi(wand.xpoints, wand.ypoints, wand.npoints, Roi.FREEROI);
							roi.setProperty("ID", id as String)
							roi.setProperty("Class", prefix)
							roi.setName(prefix+" #"+IJ.pad(id as int,padding as int))
							rois.put(id,roi)
						}
					}
				}
				//labels.setOverlay(ov)
				
				
				return rois
				
			}
			
		}
			
		//import ij.gui.Overlay
		//import ij.gui.PointRoi
		
