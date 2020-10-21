mkdir Water_Stress_Maps_jpg
mkdir Water_Stress_Maps_GTiff
mkdir Water_Stress_Maps_CSV
Rscript EDI_Values_Workflow_VLab.R
zip -r archivename.zip Water_Stress_Maps_jpg
zip -r archivename.zip Water_Stress_Maps_GTiff
zip -r archivename.zip Water_Stress_Maps_CSV
