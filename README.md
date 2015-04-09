stitch_fields

usage: stitch_images.py [-h] [-o [OUTPUT_FORMAT]] [-i [INPUT_FORMAT]]
                        [-f FIELD_STRING] [-d [DIRECTION]] [-r] [-n]
                        [path]

A small utility for field images exported from automated microscope platforms. Starts in the middle and stitches fields in a spiral pattern to create the well image.


positional arguments:
  path                  path to images  (default: current directory)

optional arguments:
  -h, --help            show this help message and exit
  -o [OUTPUT_FORMAT], --output-format [OUTPUT_FORMAT]
                        format for the stitched image (default: jpeg)
  -i [INPUT_FORMAT], --input-format [INPUT_FORMAT]
                        format for images to be stitched, can also be a list of formats (default: bmp)
  -f FIELD_STRING, --field-string FIELD_STRING
                        string immediately preceding the field number in the file name (default: _f)
  -d [DIRECTION], --direction [DIRECTION]
                        directions from the 1st field to the 2nd and 3rd, e.g. left_up =
                        3, 4, 5, 
                        2, 1, 6, 
                        9, 8, 7 (default: left_up)
  -r, --recursive       stitch images in subdirectories
  -n, --no-flip         do not flip the input images vertically
