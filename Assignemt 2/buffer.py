import arcpy
arcpy.analysis.Buffer(
    in_features="Assignemt 2/no_retail.shp",
    out_feature_class="D:\Curtis_assignment\geog4057\programming\geog4057_Curtis\Assignemt 2/retail_buffer_1.shp",
    buffer_distance_or_field="500 Meters",
    line_side="FULL",
    line_end_type="ROUND",
    dissolve_option="NONE",
    dissolve_field=None,
    method="PLANAR"
)