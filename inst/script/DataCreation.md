# Data Preparation for Experiment Data

The `inst/extdata/data-2024-05-20` file was created using experimental data in FASTQ format within the current app.

## Steps to Create the Data File

1. **Ensure Origin Distmap Files:**
   - Verify that you have the origin distmap files.
   - Update the file locations in the `server.R` code to point to the correct paths.

2. **Upload Experiment Data:**
   - Open the app and navigate to the home page.
   - Select "Method 2: Upload samplesheet and fastq directory".
   - Follow the instructions to upload your samplesheet and FASTQ directory.

3. **Load Experiment Data:**
   - After uploading, your experiment data should be visible in the "Loxcodes Experiment" section of the app.

4. **Generate RDS Data File:**
   - Use the "Download current" button to generate the RDS data file containing your experiment data.