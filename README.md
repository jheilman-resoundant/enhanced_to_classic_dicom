# enhanced_to_classic_dicom
A simply python script to convert multi-frame (enhanced) dicom to single-frame (classic) files.

## Usage:

Requires python >-3.11

Install:

python -m pip install -r requirements.txt

Convert a single file:

python enhanced_to_classic.py <path_to_enhanced_dicom_file>

Convert a directory of files:

python enhanced_to_classic.py <path_to_directory_with_enhanced_dicom>


## Theory:

Creates a classic dicom file by making a copy of the enhanced input, then copying the contents of 
Shared Functional Group Sequence and Per-Frame Functional Group sequences into the dataset,
and then streaming in one frame of pixel data.


## Notes:

- Outputs are stored in <input_directory>_classic_dicoms
- Compressed input data is not well validated, recommend  non-compressed data only