//#@ File    (label = "Input directory", style = "directory") srcFile
//#@ File    (label = "Output directory", style = "directory") dstFile
//#@ String  (label = "File extension", value=".tif") ext
//#@ String  (label = "File name contains", value = "") containString
//#@ Boolean (label = "Keep directory structure when saving", value = true) keepDirectories

#@ ImagePlus imp (label="Input Image")
//#@ Integer expansion (label="Expand By (pixels)", value=5)
#@ RoiManager rm


import ij.IJ

import ij.*
import ij.gui.*
import ij.process.*
import inra.ijpb.binary.*
import inra.ijpb.morphology.strel.*
import inra.ijpb.morphology.*
import ij.plugin.ImageCalculator
import ij.plugin.frame.RoiManager

//variables to set up
srcFile = "/Users/elena_brivio/ownCloud/PhD/PhD_MPI/Projects/Acute_Chronic/Experiments/Validations_2022/RNAscope/Analysis/Analysis_directory/Experiment4/Step_2_output/"
dstFile = "/Users/elena_brivio/ownCloud/PhD/PhD_MPI/Projects/Acute_Chronic/Experiments/Validations_2022/RNAscope/Analysis/Analysis_directory/Experiment4/Step_3_output/"
ext = ".tif"
keepDirectories = true
expansion = 5


def main() {
	srcFile.eachFileRecurse {
		name = it.getName()
		if (name.endsWith(ext)) {
			process(it, srcFile, dstFile, keepDirectories)
		}
	}
}

def process(file, src, dst, keep) {
	println "Processing $file"

	// Opening the image
	imp = IJ.openImage(file.getAbsolutePath())

	// Put your processing steps here
	def rm = RoiManager.getInstance() == null ? new RoiManager() : RoiManager.getInstance()
rm.reset()


	def rois = CellMaker.getCells(imp, expansion)

	rois.each{rm.addRoi(it)}

	rm.runCommand("Sort")

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

	// Saving the result
	relativePath = keep ?
			src.toPath().relativize(file.getParentFile().toPath()).toString()
			: "" // no relative path
	saveDir = new File(dst.toPath().toString(), relativePath)
	if (!saveDir.exists()) saveDir.mkdirs()
	saveFile = new File(saveDir, file.getName()) // customize name if needed
	IJ.saveAs(imp, "Tiff", saveFile.getAbsolutePath());

	// Clean up
	imp.close()
}

main()
