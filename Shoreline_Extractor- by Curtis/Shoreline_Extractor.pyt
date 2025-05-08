
import arcpy
import numpy as np
from arcpy.sa import *
from skimage.filters import threshold_otsu  # Importing the Otsu thresholding method from skimage
import os


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = "toolbox"

        # List of tool classes associated with this toolbox
        self.tools = [Tool]


class Tool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Shoreline Extractor"
        self.description = '''Extract shoreline using MNDWI and Otsu based on user-defined raster and study area.
        The tool will create a binary raster, convert it to polygons, and then to lines. It will also project the shoreline to the specified coordinate system.
        This tool is designed for use in ArcGIS Pro and requires the Spatial Analyst extension.
        The output will be a shapefile of the shoreline and MNDWI in tif format in the specified coordinate system.
        Developed by Curtis Amo Dwira, MSc in Civil Engineering, Louisiana State University.
        With support from Dr. Leiwang, PhD, Associate Professor, Louisiana State University.
        Dr. Ahmmed Abdalla, PhD, Assistant Professor, Louisiana State University.
        Prof. Nan Walker, PhD, Professor, Louisiana State University.
        Prof. Bernard Kumi-Boateng, PhD, Professor, University of Mines and Technology, Ghana.

        '''
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define the tool parameters."""
        params = [
            arcpy.Parameter(
                displayName="Composite Raster",
                name="composite_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input"),

            arcpy.Parameter(
                displayName="Green Band Name (e.g., B3, Band_3, Layer_3)",
                name="green_band_name",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),

            arcpy.Parameter(
                displayName="SWIR Band Name (e.g., B11, Band_11, Layer_11)",
                name="swir_band_name",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),

            arcpy.Parameter(
                displayName="Study Area",
                name="study_area",
                datatype="GPFeatureLayer",
                parameterType="Required",
                direction="Input"),

            arcpy.Parameter(
                displayName="Output Folder",
                name="output_folder",
                datatype="DEFolder",
                parameterType="Required",
                direction="Input"),

            arcpy.Parameter(
                displayName="Output File Prefix",
                name="output_prefix",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),

            arcpy.Parameter(
                displayName="Output Coordinate System",
                name="output_crs",
                datatype="GPCoordinateSystem",
                parameterType="Required",
                direction="Input")
        ]
        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        params = [
            arcpy.Parameter(
                displayName="Composite Raster",
                name="composite_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input"),

            arcpy.Parameter(
                displayName="Green Band Name (e.g., B3, Band_3, Layer_3)",
                name="green_band_name",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),

            arcpy.Parameter(
                displayName="SWIR Band Name (e.g., B11, Band_11, Layer_11)",
                name="swir_band_name",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),

            arcpy.Parameter(
                displayName="Study Area",
                name="study_area",
                datatype="GPFeatureLayer",
                parameterType="Required",
                direction="Input"),

            arcpy.Parameter(
                displayName="Output Folder",
                name="output_folder",
                datatype="DEFolder",
                parameterType="Required",
                direction="Input"),

            arcpy.Parameter(
                displayName="Output File Prefix",
                name="output_prefix",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),

            arcpy.Parameter(
                displayName="Geographic Coordinate System",
                name="output_crs",
                datatype="GPCoordinateSystem",
                parameterType="Required",
                direction="Input")
        ]
        return params

    def execute(self, parameters, messages):
        arcpy.CheckOutExtension("Spatial")

        # Input values
        composite = parameters[0].value
        green_band_name = parameters[1].valueAsText
        swir_band_name = parameters[2].valueAsText
        study_area = parameters[3].value
        output_folder = parameters[4].valueAsText
        prefix = parameters[5].valueAsText
        output_crs = parameters[6].value

        # File paths
        binary_raster_path = os.path.join(output_folder, f"{prefix}_binary.tif")
        mndwi_raster_path = os.path.join(output_folder, f"{prefix}_mndwi.tif")
        polygon_path = os.path.join(output_folder, f"{prefix}_water_polygon.shp")
        valid_polygons = os.path.join(output_folder, f"{prefix}_valid_polygons.shp")
        shoreline_raw_path = os.path.join(output_folder, f"{prefix}_shoreline_raw.shp")
        shoreline_proj_path = os.path.join(output_folder, f"{prefix}_shoreline_projected.shp")

        # Project study area if provided
        if study_area:
            projected_study_area = os.path.join(output_folder, f"{prefix}_projected_study_area.shp")
            arcpy.Project_management(study_area, projected_study_area, arcpy.Describe(composite).spatialReference)
            processed_raster = ExtractByMask(composite, projected_study_area)
            messages.addMessage("Study area reprojected and raster clipped.")
        else:
            processed_raster = arcpy.Raster(composite)
            messages.addMessage("No study area provided. Using full raster.")

        # Use Describe to get catalogPath
        raster_path = arcpy.Describe(composite).catalogPath
        green_raster = arcpy.Raster(f"{raster_path}/{green_band_name}")
        swir_raster = arcpy.Raster(f"{raster_path}/{swir_band_name}")

        green_array = arcpy.RasterToNumPyArray(green_raster)
        swir_array = arcpy.RasterToNumPyArray(swir_raster)

        # Compute MNDWI and threshold
        mndwi_array = (green_array.astype(float) - swir_array.astype(float)) / (green_array + swir_array + 1e-10)
        threshold = threshold_otsu(mndwi_array[~np.isnan(mndwi_array)])
        binary_array = np.where(mndwi_array > threshold, 1, 0).astype(np.uint8)

        # Save MNDWI raster (clipped to study area if provided)
        ll = processed_raster.extent.lowerLeft
        cs = processed_raster.meanCellWidth
        mndwi_raster = arcpy.NumPyArrayToRaster(mndwi_array, ll, cs, cs)
        if study_area:
            mndwi_raster = ExtractByMask(mndwi_raster, projected_study_area)
        mndwi_raster.save(mndwi_raster_path)
        messages.addMessage("MNDWI raster saved.")

        # Save binary raster (clipped to study area if provided)
        binary_raster = arcpy.NumPyArrayToRaster(binary_array, ll, cs, cs, value_to_nodata=255)
        if study_area:
            binary_raster = ExtractByMask(binary_raster, projected_study_area)
        binary_raster.save(binary_raster_path)
        messages.addMessage("Binary raster saved.")

        # Raster → Polygon
        arcpy.RasterToPolygon_conversion(binary_raster_path, polygon_path, "SIMPLIFY", "Value")
        messages.addMessage("Water polygons created (with simplification).")

        arcpy.RepairGeometry_management(polygon_path)
        messages.addMessage("Polygon geometry repaired.")

        # Filter only valid water polygons
        arcpy.MakeFeatureLayer_management(polygon_path, "poly_lyr")
        arcpy.SelectLayerByAttribute_management("poly_lyr", "NEW_SELECTION", '"gridcode" = 1')
        arcpy.CopyFeatures_management("poly_lyr", valid_polygons)
        messages.addMessage("Filtered valid water polygons for conversion.")

        try:
            arcpy.PolygonToLine_management(valid_polygons, shoreline_raw_path)
            messages.addMessage("Shoreline (unprojected) extracted.")
        except Exception as e:
            messages.addErrorMessage("Polygon to Line conversion failed: " + str(e))
            raise

        # Define spatial reference for the unprojected shoreline using the output CRS
        arcpy.DefineProjection_management(shoreline_raw_path, output_crs)
        messages.addMessage(f"Spatial reference defined for: {shoreline_raw_path} using output CRS.")

        # Project shoreline (now it's essentially projecting from the output CRS to the output CRS, which is fine)
        arcpy.Project_management(shoreline_raw_path, shoreline_proj_path, output_crs)
        messages.addMessage(f"Projected shoreline saved: {shoreline_proj_path}")

        # Add to map
        try:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            map_obj = aprx.activeMap
            map_obj.addDataFromPath(mndwi_raster_path)
            map_obj.addDataFromPath(shoreline_proj_path)
            messages.addMessage("MNDWI and shoreline layers added to current map.")
        except Exception as e:
            messages.addWarningMessage("Could not add to map: " + str(e))

        arcpy.CheckInExtension("Spatial")
        return

    def postExecute(self, parameters):
        """Just Believe."""
        return