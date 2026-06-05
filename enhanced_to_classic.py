'''
enh_to_classic.py

A script to convert multi-frame DICOM (aka "enhanced dicom") to single-frame DICOM (aka classic dicom)

'''

import os
import sys
import re
import argparse
import unicodedata
import pydicom, pydicom.uid
from copy import deepcopy
from pathlib import Path



MIN_PY_VERSION = (3, 11)
if sys.version_info < MIN_PY_VERSION:
    sys.stderr.write(
        f"Python {MIN_PY_VERSION[0]}.{MIN_PY_VERSION[1]}+ required, "
    )
    sys.exit(1)

def convert_dicom_enhanced_to_classic(enhanced_dcm:pydicom.Dataset) -> list[pydicom.Dataset]:
    if is_enhanced_dcm(enhanced_dcm):
        return _enhanced_to_classic(enhanced_dcm)
    return [enhanced_dcm] # wasn't enhanced, give it back

def _enhanced_to_classic(enh_dcm:pydicom.Dataset) -> list[pydicom.Dataset]:
    # Pull pixel data into an array of frames.
    # Check if it is already in VM>1 pixel data array or stored as one long string.
    #todo: use transfer syntax to determine if data is compressed, replacing pixeldata.VM
    num_frames = int(enh_dcm.get('NumberOfFrames',0))
    # compute frame sizes
    pixel_bytes = 2
    pixeldata_vr = 'OW'
    if enh_dcm.get('BitsAllocated',0) == 8:
        pixel_bytes = 1
        pixeldata_vr = "OB"
    frame_size = enh_dcm.Rows * enh_dcm.Columns * pixel_bytes 

    # extract image data.
    # if it is already multi-frame pixel data (typ for compressed), set the array directly.
    # if it is a long string, split it up
    if num_frames == enh_dcm['PixelData'].VM:
        image_frames = enh_dcm.PixelData
    elif frame_size > 0:
        image_frames=[enh_dcm.PixelData[k*frame_size:(k+1)*frame_size] for k in range(num_frames)]
    else:
        print('WARNING _enhanced_to_classic(): failed to separate image frames (frame size = 0, pixel data VM = 0)')
        return []


    # reference enh dicom sequences
    shared_seq = enh_dcm.get('SharedFunctionalGroupsSequence',[])
    per_frame_seq = enh_dcm.get('PerFrameFunctionalGroupsSequence',[])

    # Consistency checks:
    if len(image_frames) != num_frames:
        raise ValueError(f'WARNING convert_dcm_enhanced_to_classic(): num_frames={num_frames} but {len(image_frames)} image frames extracted from PixelData')
    if len(per_frame_seq) != num_frames:
        raise ValueError(f'NumberOfFrames ({num_frames}) != PerFrameFunctionalGroupsSequence.VM ({len(per_frame_seq)})')
    if frame_size*num_frames != len(enh_dcm.PixelData):
        raise ValueError(f'PixelData length ({len(enh_dcm.PixelData)}) != Rows*Columns ({frame_size})')

    # Preparing Base dataset.
    # We choose to start with a copy of the enhanced dicom, remove the shared and per-frame functional groups
    # We could alternatively start from scratch but it's likely that certain metadata 
    # isn't repeated in the shared group and would be missed.
    # Removing tags from dcm_base
    dcm_base = pydicom.Dataset()
    try:
        dcm_base=deepcopy(enh_dcm)  # read b/c deepcopy threw errors
    except:
        print(f"could not deepcopy enhanced input, attempting to re-read from source file")
        dcm_base=pydicom.dcmread(enh_dcm.filename)  # read b/c deepcopy threw errors
    if dcm_base.get('filename',None) == None:
        print(f"WARNING - could not start base dicom from input, using default Dataset")
    for tag in ['PerFrameFunctionalGroupsSequence','SharedFunctionalGroupsSequence']:
        try:
            del dcm_base[tag]
        except Exception:
            pass

    #crop pixel size down to one frame.
    dcm_base.PixelData = dcm_base.PixelData[:frame_size]
    dcm_base.NumberOfFrames = 1
    # todo: set pixel size VR (supposedly set by pydicom?  idk)

    # copy all shared sequence elements into base so we only have to do that once
    if shared_seq != None:
        for shared_elem in shared_seq[0]:
            update_from_sequence(dcm_base, shared_elem)

    #Creating classic DICOMs
    classic_datasets=[]
    for k, frame_entry in enumerate(per_frame_seq):
        dcm_classic = deepcopy(dcm_base)
        for frame_elem in frame_entry:
            update_from_sequence(dcm_classic,frame_elem) # MR Echo Sequence
        dcm_classic.PixelData = image_frames[k]
        dcm_classic.file_meta.MediaStorageSOPClassUID = dcm_classic.SOPClassUID
        slicelocation_from_imageposition(dcm_classic)
        acqution_date_and_time(dcm_classic)
        classic_datasets.append(dcm_classic)
    return classic_datasets

def acqution_date_and_time(dcm_classic):
    acquisition_datetime = getattr(dcm_classic, "AcquisitionDateTime", None)
    if acquisition_datetime:
        acquisition_datetime = str(acquisition_datetime)
        if len(acquisition_datetime) >= 8:
            dcm_classic.AcquisitionDate = acquisition_datetime[:8]
            acquisition_time = acquisition_datetime[8:]
            if acquisition_time:
                dcm_classic.AcquisitionTime = acquisition_time

def update_from_sequence(dcm_out:pydicom.Dataset, data_in):
    ''' this function is intended to copy contents of a per-frame or shared sequence:
    call it with the elements inside one of these sequences.
    If it is a sub-sequence, it will discard the sequence wrapper.
    '''
    # Case 1: top-level input is a DataElement sequence
    if isinstance(data_in, pydicom.DataElement) and data_in.VR == "SQ":
        for item in data_in.value:
            if isinstance(item, pydicom.DataElement):
                dcm_out.add(deepcopy(item))
            elif isinstance(item, pydicom.Dataset):
                for subelem in item:
                    dcm_out.add(deepcopy(subelem))
            else:
                raise TypeError(f"Unexpected SQ item type: {type(item)}")

    # Case 2: top-level input is already a Dataset
    elif isinstance(data_in, pydicom.Dataset):
        for elem in data_in:
            dcm_out.add(deepcopy(elem))

    # Case 3: top-level input is a single non-sequence DataElement
    elif isinstance(data_in, pydicom.DataElement):
        dcm_out.add(deepcopy(data_in))

    else:
        print(f"Unsupported input type not copied: {type(data_in)}")
        # raise TypeError(f"Unsupported input type: {type(data_in)}")
            


# def generate_new_series_uid(uid:str)->str:
#     """
#     generate_new_series_uid():
#     Returns uid.1 if there is space, otherwise returns a new random uid
#     """
#     if isinstance(uid,str) and (len(uid)<=62):
#         return f'{uid}.1'
#     return pydicom.uid.generate_uid()

# def generate_new_sop_uid(uid:str,instance_number:int)->str:
#     """
#     generate_new_sop_uid():
#     Returns uid.instance_number if there is space, otherwise returns a new random uid
#     """
#     if isinstance(uid,str) and (len(uid)<=(63-len(str(instance_number)))):
#         return f'{uid}.{instance_number}'
#     return pydicom.uid.generate_uid()


def is_enhanced_dcm(dcm:pydicom.Dataset)->bool:
    """
    is_enhanced_dcm():
    Determines if dcm is an enhanced DICOM
    """
    return dcm.get('PixelData') and int(dcm.get('NumberOfFrames',0))>1


def slicelocation_from_imageposition(ds:pydicom.Dataset):
    ''' Set SliceLocation of ds if image position and orienatation are standard planar'''
    if ds.get("SliceLocation") != None:
        return
    ipp = ds.get("ImagePositionPatient")
    iop = ds.get("ImageOrientationPatient")
    if (ipp) and (not iop):
        print("Missing image orientation, forcing axial view")  # HACK
        ds.SliceLocation = ipp[2]
        ds.ImageOrientationPatient = [1,0,0,0,1,0]
        return
    if (not ipp) or (not iop):
        return
    if   iop[2]==0 and iop[5]==0: # Rz = Cz = 0, axial
        ds.SliceLocation = ipp[2]
    elif iop[1]==0 and iop[4]==0: # Ry = Cy = 0, coronal
        ds.SliceLocation = ipp[1]
    elif iop[0]==0 and iop[3]==0: # Rx = Cx = 0, saggital
        ds.SliceLocation = ipp[0]
    return


def save_classic(classic_dcm:list[pydicom.Dataset], output_folder:str|Path, human_readable=False)->list[str]:
    """
    save_classic():
    Saves a list of dicom files (ostensibly classic)
    Returns list of save location file paths

    Saves each file by its SOPInstanceUID, unless base_filename is set, in which case it saves as 
        'f{base_filename}{InstanceNumber}.dcm'
    """
    # test create output directory
    output_folder = Path(output_folder)
    # try:
    #     output_folder.mkdir(parents=True, exist_ok=True)
    # except:
    #     raise RuntimeError(f"Cannot create output folder {output_folder}")
    
    #Saving classic DICOMs
    new_file_paths=[]
    for dcm in classic_dcm:
        output_filename = str(dcm.get('SOPInstanceUID','ERR_NO_SOP_UID'))
        # output_subfolder = str(dcm.get('SeriesInstanceUID','ERR_NO_SERIES_UID'))
        if human_readable:
            output_filename = f"IMG{int(dcm.get('InstanceNumber',0)):05d}_{output_filename}"
            # output_subfolder = f"{dcm.get('SeriesNumber','X')}_{dcm.get('SeriesDescription','NODESC')}"
        file_path_out = output_folder / Path(f'{output_filename}.dcm')
        # file_path_out = Path(sanitize_path(str(file_path_out)))
        if file_path_out.exists():
            print(f"overwriting {file_path_out}")
        try:
            file_path_out.parent.mkdir(parents=True, exist_ok=True)
        except:
            raise RuntimeError(f"Cannot create output folder {output_folder}")
        dcm.save_as(str(file_path_out))
        new_file_paths.append(file_path_out)
    return new_file_paths




_WINDOWS_RESERVED_NAMES = {
    "CON", "PRN", "AUX", "NUL",
    *(f"COM{i}" for i in range(1, 10)),
    *(f"LPT{i}" for i in range(1, 10)),
}

def sanitize_path_component(
    name: str,
    replacement: str = "_",
    for_windows: bool = True,
    for_posix: bool = True,
    max_length: int = 255,
) -> str:
    """
    Sanitize a single path component (file or directory name) so it is safe on
    both Windows and POSIX systems by default.

    - Replaces characters illegal on Windows and POSIX with `replacement`
    - Removes control characters
    - Trims trailing spaces and dots (Windows restriction)
    - Avoids Windows reserved names by appending an underscore
    - Collapses multiple replacement chars into one
    """
    if not name:
        return replacement or "_"

    # Normalize unicode to a common form
    name = unicodedata.normalize("NFKC", name)

    illegal_chars = set()
    if for_windows:
        illegal_chars.update('<>:"/\\|?*')
    if for_posix:
        illegal_chars.add("/")  # not allowed inside a single component

    # Always disallow null byte
    illegal_chars.add("\0")

    def _sanitize_char(ch: str) -> str:
        # Control chars (0-31 and 127) are not allowed on Windows
        if ord(ch) < 32 or ord(ch) == 127:
            return replacement
        if ch in illegal_chars:
            return replacement
        return ch

    sanitized = "".join(_sanitize_char(ch) for ch in name)

    # Collapse multiple replacement characters
    if replacement:
        sanitized = re.sub(re.escape(replacement) + r"+", replacement, sanitized)

    # Strip leading/trailing whitespace
    sanitized = sanitized.strip()

    # Windows: no trailing space or dot
    if for_windows:
        sanitized = sanitized.rstrip(" .")

    # Don't allow an empty name
    if not sanitized:
        sanitized = replacement or "_"

    # Avoid reserved names on Windows (case-insensitive, without extension)
    if for_windows:
        base, ext = os.path.splitext(sanitized)
        if base.upper() in _WINDOWS_RESERVED_NAMES:
            sanitized = f"{sanitized}_"

    # Enforce max length (typical filesystem limit per component)
    if max_length is not None and max_length > 0:
        sanitized = sanitized[:max_length]

    return sanitized


def sanitize_path(
    path: str|Path,
    replacement: str = "_",
    for_windows: bool = True,
    for_posix: bool = True,
    max_component_length: int = 255,
) -> str:
    """
    Sanitize an entire relative/absolute path, cleaning each component so that
    none contain characters illegal on Windows or POSIX (by default).

    - Keeps path structure but cleans each part between separators.
    - Rebuilds the path using the current platform's os.sep.
        (Windows accepts both '\\' and '/'; Linux/macOS use '/'.)
    """
    if not path:
        return sanitize_path_component(
            "",
            replacement,
            for_windows,
            for_posix,
            max_component_length,
        )

    drive = ""
    rest = path

    # Handle Windows drive letter, if present (e.g. C:\foo\bar)
    if len(rest) >= 2 and rest[1] == ":":
        drive, rest = rest[:2], rest[2:]

    # Split on both types of separators
    parts = re.split(r"[\\/]+", rest)

    sanitized_parts = [
        sanitize_path_component(
            part,
            replacement=replacement,
            for_windows=for_windows,
            for_posix=for_posix,
            max_length=max_component_length,
        )
        for part in parts
        if part != ""  # avoid empty segments from // or leading/trailing slashes
    ]

    # Reconstruct the path using os.sep for the current platform
    sep = os.sep
    sanitized = sep.join(sanitized_parts)

    # Reattach drive if any
    if drive:
        sanitized = drive + sep + sanitized

    # Preserve leading slash if original path was absolute (POSIX style or UNC-like)
    if not drive and (path.startswith("/") or path.startswith("\\")):
        sanitized = sep + sanitized

    return sanitized

def main(input_arg:str|Path, recursive:bool = False):
    if os.path.isfile(input_arg):
        print(f"Converting enhanced file to classic: {input_arg}")
        file_list = [Path(input_arg)]
        src_dir = Path(input_arg).parent
        output_dir = Path(f"{os.path.dirname(input_arg)}_classic_dicom")
    else:
        print(f"Creating classic DICOM from all enhanced files in directory {input_arg}")
        print(f"  Recursive scan is {recursive}")
        src_dir = Path(input_arg)
        output_dir = Path(f"{input_arg}_classic_dicom")
        paths = src_dir.rglob('*') if recursive else src_dir.iterdir()
        file_list = [p for p in paths if p.is_file()]

    print(f"Outputs will be saved to {output_dir}")
    print(f"Attempting to convert {len(file_list)} files...")
    output_dir.mkdir(exist_ok=True)
    for enh_file in file_list:
        try:
            enh_dcm = pydicom.dcmread(enh_file)
        except:
            continue # quietly continue if not a DICOM
        try:
            if 'PixelData' not in enh_dcm:
                continue
            classic_files = convert_dicom_enhanced_to_classic(enh_dcm)
            print(f'Saving {enh_file} as {len(classic_files)} classic DICOMs')
            ehn_file_output_dir = output_dir / enh_file.parent.relative_to(src_dir)
            print(f'   output directory {ehn_file_output_dir}')
            save_classic(classic_files, ehn_file_output_dir)
        except Exception as e: 
            print(f"ERROR converting and saving {enh_file}: {e}")
            continue

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert enhanced DICOM files to classic DICOM files."
    )
    parser.add_argument("input_arg", type=Path)
    parser.add_argument("-r", "--recursive", action="store_true")
    args = parser.parse_args()
    main(args.input_arg, args.recursive)


