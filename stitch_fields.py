#!/usr/bin/env python

'''
Small script for stitching together cellomics images. The user needs to define the path, 
input/output image extensions and the scan direction (fieldlookuptable for the microscope)

Joel 2015-04-03
'''
from __future__ import print_function
from __future__ import division
from skimage import exposure #only for rescaling now, maybe replace PIL completely in the future
from itertools import cycle
from PIL import Image # could be either pillow or PIL?
import numpy as np
import argparse
import logging
import shutil
import time
import sys
import re
import os

# Putting the main logic of the program into main() can actually make it run faster
# this is since the local scopes are implememnted as arrays, which is faster then
# the dictionary implementation of global scopes. More details:
# http://stackoverflow.com/questions/11241523/why-does-python-code-run-faster-in-a-function
def main():
    # Set up the command line arguments
    # the formatter_class adds the default value to optional arguments
    # argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # for positional arguments, one has to add the same string that the optional gets automatically
    parser = argparse.ArgumentParser(description='A small utility for field images ' \
        'exported from automated microscope platforms.\nStarts in the middle and ' \
        'stitches fields in a spiral pattern to create the well image.\n' \
        'Each well and channel need to be in a separate directory. Use `-c` to sort '\
        'images into\ndirectories automatically. Make sure to specify the correct field '\
        'and well string (-f, -w).\nExample usage when the images from all wells are in the '\
        'same directory:\n\npython stitch_fields.py -cr -f <field_prefix> -w <well_prefix>', 
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('path', default=os.getcwd(), nargs='?',
        help='path to images  (default: current directory)')
    parser.add_argument('-o', '--output-format', nargs='?', default='jpeg',
        help='format for the stitched image (default: %(default)s)')
    parser.add_argument('-i', '--input-format', nargs='?', default='bmp',
        help='format for images to be stitched, can also be a list of formats (default: %(default)s)')
    parser.add_argument('-f', '--field-prefix', default='f',
        help='string immediately preceding the field number in the file name (default: %(default)s)')
    parser.add_argument('-w', '--well-prefix', default='0001_',
        help='string immediately preceding the well id in the file name (default: %(default)s)')
    parser.add_argument('-c', '--channel-prefix', default='d',
        help='string immediately preceding the channel id in the file name (default: %(default)s)')
    parser.add_argument('-d', '--scan-direction', nargs='?', default='left_down',
        help='The directions from the 1st field to the 2nd and 3rd, e.g. left_down =\n' \
        '9, 8, 7, \n' \
        '2, 1, 6, \n' \
        '3, 4, 5 (default: %(default)s)')
    parser.add_argument('-e', '--sort-wells', action='store_true',
        help='if all images are in the same directory, subfolders MUST be created. ' \
        'Can be used to only sort files if -r is omitted')
    parser.add_argument('-a', '--sort-channels', action='store_true',
        help='Create subfolders for each channel based on the specified channel prefix.')
    parser.add_argument('-r', '--recursive', action='store_true',
        help='Stitch images in subdirectories')
    parser.add_argument('-s', '--rescale_intensity', action='store_true',
        help='Scales the bit range from the input image to 8 bits. Our cellomics has a 12-bit camera and ' \
	'values are usually not higher than ~1500 so even rescaling to 12 bit (4095) would render dark images. ' \
	'Find the ~max (99.999th percentile in the entire plate and scales everything accordingly.' \
	'Now we can export 16-bit tiff files instead of converting to 8 bit in cellomics, which would ' \
	'makes them appear patchy.')
    parser.add_argument('--flip', default='none', choices=['horizontal', 'vertical', 'both', 'none'])
    #Initialize some variables
    args = parser.parse_args()
    # PIL's image function takes 'jpeg' instead of 'jpg' as an argument. We want to be able to
    # specify the image format to this function while defining the image extensions as 'jpg'.
    if args.output_format.lower() == 'jpeg':
        output_format = 'jpg'
    else:
        output_format = args.output_format.lower()
   # timestamp = str(int(time.time()))[3:]
    input_format = set((args.input_format.lower(),)) #can add extra ext here is needed, remember to not have same as stiched
    logging.basicConfig(filename='well_stitch.log',level=logging.DEBUG, format='%(message)s')
    # Print out the runtime parameters and store them in a log file
    for key in vars(args): # Returns a dictionary instead of a Namespace object
        print(key +'\t', vars(args)[key])
        #logging.info(key +'\t', vars(args)[key])

    # Sort channels and create subfolders
    channel_names = [''] # if channel_sort is not specified, this helps 
    if args.sort_channels:
        print('\n Moving images to channel subfolders...')
        channel_names = sort_channels(args.path, args.channel_prefix, input_format)
    
    # Sort wells and create subfolders 
    if args.sort_wells:
        print('\n Moving images to well subfolders...')
        well_names = sort_wells(args.path, args.well_prefix, input_format, channel_names)
    
    # Main program
    if args.recursive:
        # Create a new directory. Append a number if it already exists.
	print('\nStitching wells...')
	stitched_dir = os.path.join(args.path, 'stitched_wells')
	dir_suffix = 1
	while os.path.exists(stitched_dir): 
	    dir_suffix += 1
	    stitched_dir = os.path.join(args.path, 'stitched_wells_' + str(dir_suffix))
	os.makedirs(stitched_dir)
	logging.info('Created directory ' + os.path.join(stitched_dir))
	# Loop through only the well directories, the current directory does not need to be
	# included as the files will already be sorted into subdirectories
	dirs = [name for name in os.listdir(args.path) if os.path.isdir(name)]
	well_dirs = [name for name in dirs if not name.startswith('stitched_wells')]
	num_dirs = len(well_dirs)
	for num, dir_name in enumerate(sorted(well_dirs, key=nat_key), start=1):
	    # Progress bar. The trailing space in the print function is needed to update that position.
	    # Otherwise that would be forzen when moving from a two digit to a one digit number.
	    progress = int(num / num_dirs * 100)
	    print('{0}% {1} '.format(progress, dir_name), end='\r')
	    sys.stdout.flush()
	    
	    # add for channel here
	    
	    imgs, zeroth_field, max_int = find_images(dir_name, input_format, args.flip, args.field_prefix)
	    # If there are images in the directory
	    if imgs:
		if args.rescale_intensity:
		    imgs = rescale_intensities(imgs, max_int)
		fields, arr_dim, moves, starting_point = spiral_structure(dir_name, input_format, args.scan_direction)
		img_layout = spiral_array(fields, arr_dim, moves, starting_point, zeroth_field)
		stitched_well = stitch_images(imgs, img_layout, dir_name, output_format, arr_dim, stitched_dir)
		stitched_well_name = os.path.join(stitched_dir, os.path.basename(dir_name) + '.' + output_format)
		stitched_well.save(stitched_well_name, format=args.output_format)
		logging.info('Stitched image saved to ' + stitched_well_name + '\n')
	    else:
		logging.info('No images found in this directory\n')
    else:
        imgs, zeroth_field = find_images(args.path, input_format, args.flip, args.field_prefix, args.rescale_intensity)
        stitched_dir = os.path.join(args.path, 'stitched_wells')
        #if there are images in the directory
        if imgs:
            fields, arr_dim, moves, starting_point = spiral_structure(args.path, input_format, args.scan_direction)
            img_layout = spiral_array(fields, arr_dim, moves, starting_point, zeroth_field)
            stitched_name = stitch_images(imgs, img_layout, args.path, args.output_format, arr_dim)
            logging.info('Stitched image saved to ' + stitched_name + '\n')
        else:
            logging.info('No images found\n')
    #import time
    #time.sleep(2)
    #os.rename('./well_stitch.log', os.path.join(stitched_dir, 'well_stitch.log')) 
    
    print('\n\nStitched well images can be found in ' + stitched_dir + '.\nPlease check the log file for which images and what field layout were used to create the stitched image.\nDone.')

def rescale_intensities(imgs, max_int):
    print(max_int)
    for fnum, img in imgs.iteritems():
	img_norm = np.array(img) - np.array(img).min()
	img_norm = img_norm / (max_int / 256)
	imgs[fnum] = Image.fromarray(np.uint8(img_norm))
	print(img_norm.ravel().min(), img_norm.ravel().max())
	
	#img_stretched = exposure.rescale_intensity(np.array(img), in_range=(0, max_int), out_range=('uint8'))
	#print(img_stretched.ravel().min(), img_stretched.ravel().max())
	#imgs[fnum] = Image.fromarray(np.uint8(img_stretched))
    return imgs



## Functions ##

# Define movement function for filling in the spiral array
def move_right(x,y):
    return x, y +1


def move_down(x,y):
    return x+1,y


def move_left(x,y):
    return x,y -1


def move_up(x,y):
    return x -1,y

    
def spiral_structure(dir_path, input_format, scan_direction):
    '''
    Define the movement scheme and starting point for the field layout
    '''
    #find the number of fields/images matching the specified extension(s)
    fields = len([fname for fname in os.listdir(dir_path) if fname[-3:].lower() in input_format])
    #size the array based on the field number, array will be squared
    arr_dim = int(np.ceil(np.sqrt(fields)))
    #define the movement schema and find the starting point (middle) of the array
    if scan_direction == 'down_left':
        moves = [move_down, move_left, move_up, move_right]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:        
            starting_point = (int(arr_dim/2)-1, int(arr_dim/2))
    elif scan_direction == 'left_down':
        moves = [move_left, move_down, move_right, move_up]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2)-1, int(arr_dim/2))
    elif scan_direction == 'down_right':
        moves = [move_down, move_right, move_up, move_left]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2)-1, int(arr_dim/2)-1)
    elif scan_direction == 'right_down':
        moves = [move_right, move_down, move_left, move_up]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2)-1, int(arr_dim/2)-1)
    elif scan_direction == 'up_left':
        moves = [move_up, move_left, move_down, move_right]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
    elif scan_direction == 'left_up':
        moves = [move_left, move_up, move_right, move_down]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
    elif scan_direction == 'up_right':
        moves = [move_up, move_right, move_down, move_left]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2), int(arr_dim/2)-1)
    elif scan_direction == 'right_up':
        moves = [move_right, move_up, move_left, move_down]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2), int(arr_dim/2)-1)

    return fields, arr_dim, moves, starting_point
        

def gen_points(end, moves, starting_point):
    '''
    Generate coordinates for each field number. From stackoverflow
    '''
    _moves = cycle(moves)
    n = 1
    pos = starting_point
    times_to_move = 1

    yield n,pos

    while True:
        for _ in range(2):
            move = next(_moves)
            for _ in range(times_to_move):
                if n >= end:
                    return
                pos = move(*pos)
                n+=1
                yield n,pos

        times_to_move+=1
              

def spiral_array(fields, arr_dim, moves, starting_point, zeroth_field):
    '''
    Fill the array in the given direction
    '''
    #create an array of zeros and then fill with a number not used for any field
    img_layout = np.zeros((arr_dim, arr_dim), dtype=int)
    img_layout[:] = -1 #TODO this means that the zeroth field will be put in multiple places... fixed?
    #create a different layout depending on the numbering of the first field
    if zeroth_field:
        for point, coord in list(gen_points(fields, moves, starting_point)):
            img_layout[coord] = point -1
        logging.info('\nField layout:')    
        for row in np.ma.masked_equal(img_layout, -1): #TODO fix so that lines are showing for unused field
            logging.info(' '.join(['{:>2}'.format(str(i)) for i in row]))

    else:
        for point, coord in list(gen_points(fields, moves, starting_point)):
            img_layout[coord] = point
        logging.info('Field layout:')    
        for row in np.ma.masked_equal(img_layout, -1):
            logging.info(' '.join(['{:>2}'.format(str(i)) for i in row]))

    return img_layout
              

def find_images(dir_path, input_format, flip, field_str):
    '''
    Create a dictionary with the field numbers as keys to the field images
    '''
    zeroth_field = False #changes if a zeroth field is found in 'find_images'
    imgs = {}
    max_ints = []
    logging.info('----------------------------------------------')
    logging.info(dir_path)
    #go through each directory
    for fname in os.listdir(dir_path):
        if fname[-3:].lower() in input_format:
            logging.info(fname)
            #find the index of the field identifier string in the file name
            field_ind = fname.index(field_str) + len(field_str)
            fnum = int(''.join([str(int(char)) for char in fname[field_ind:field_ind+2] if char.isdigit()]))
            # If field 0 is encountered, change start numbering of the array
            if fnum == 0:
                zeroth_field = True
            # The default is to flip horizontally since this is the most common case
            if flip == 'none':
                imgs[fnum] = Image.open(os.path.join(dir_path, fname))
            elif flip == 'horizontal':
                imgs[fnum] = Image.open(os.path.join(dir_path, fname)).transpose(Image.FLIP_LEFT_RIGHT)
            elif flip == 'vertical':
                # dunno why we need to flip...
                imgs[fnum] = Image.open(os.path.join(dir_path, fname)).transpose(Image.FLIP_TOP_BOTTOM)
            elif flip == 'both':
                # I don't think thei sould ever be the case, it could just be adjusted with another 
                # spiral rotation, but putting it here for completion
                imgs[fnum] = Image.open(os.path.join(dir_path, fname)).transpose(Image.FLIP_TOP_BOTTOM).transpose(Image.FLIP_LEFT_RIGHT)
	    # Collect max intensities here instead of looping through an extra time
	    max_ints.append(np.percentile(np.array(imgs[fnum]).ravel(), 99.999)) # make this a user variable
	 #   min_ints.append(np.percentile(np.array(imgs[fnum]).ravel(), 0.11)) # make this a user variable
#            if rescale_intensity:
#		img_stretched = exposure.rescale_intensity(np.array(imgs[fnum]), in_range='uint12', out_range=('uint8'))
#		imgs[fnum] = Image.fromarray(np.uint8(img_stretched))
    print(max_ints)
    if max_ints != []:
       # max_ints = max(max_ints) # the highest intensity in the entire plate
        max_ints = np.percentile(np.array(max_ints), 70) # the highest intensity in the entire plate
        print(max_ints)
    return imgs, zeroth_field, max_ints




#stitch the image row by row
def stitch_images(imgs, img_layout, dir_path, output_format, arr_dim, stiched_dir):
    '''
    Stitch images by going row and column wise in the img_layout and look up
    the number of the image to place at the (row, col) coordinate. So not filling
    in a spiral but using the spiral lookuptable instead.
    '''
    # Create the size of the well image to be filled in
    width, height = imgs[1].size
    num = 0
    stitched_well = Image.new('RGB', (width*arr_dim, height*arr_dim))
    for row in xrange(0, width*arr_dim, width):
        for col in xrange(0, height*arr_dim, height):
            #since the image is filled by row and col instead of sprial, this
            #error catching is needed for the empty places
            try:
                #'stitch' fields by pasting them at the appropriate place in the black background
                stitched_well.paste(imgs[img_layout.ravel()[num]], (col, row))
            except KeyError:
                pass
            num += 1
    #save image
#    stitched_name = os.path.join(dir_path, 'stitched_wells/stitched_' + timestamp + '.' + output_format)
#    stitched.save(stitched_name, format=output_format)
    # stitched.show()
    
    return stitched_well


def sort_channels(dir_path, channel_str, input_format):
    well_names = []
    num = 0
    for fname in os.listdir(dir_path):
        if fname[-3:].lower() in input_format:
	    well_ind = fname.index(channel_str) + len(channel_str)
            well_name = [fname[well_ind]]
            well_name = ''.join([char for sublist in well_name for char in sublist])
            well_names.append(well_name)
            #if the well name is not the same as the previous or if it is the first
            #file, make a new directory
            if well_name != well_names[num-1] or num == 0:
                num += 1
                if not os.path.exists(os.path.join(dir_path, well_name)):
                    os.makedirs(os.path.join(dir_path, well_name))
            logging.info('moving ./' + fname + ' to ./' + os.path.join(well_name, fname))
            shutil.move(os.path.join(dir_path, fname), os.path.join(os.path.join(dir_path, well_name), fname))       
    logging.info('created well directories ' + str(set(well_names)))

    return set(well_names)

def sort_wells(dir_path, well_str, input_format, channel_names):
    well_names = []
    num = 0
#    for channel_name in channel_names:
#        print('cahnnel name = ' + channel_name)
#	channel_path = os.path.join(dir_path, channel_name)
	#go through all the files to find the well names
    for fname in os.listdir(dir_path):
	print(fname)
	if fname[-3:].lower() in input_format:
	    #find the well id using the provided helper string and add it to the list
	    well_ind = fname.index(well_str) + len(well_str)
	    well_name = [fname[well_ind]]
	    well_name.append([str(int(char)) for char in fname[well_ind+1:well_ind+3] if char.isdigit()])
	    well_name = ''.join([char for sublist in well_name for char in sublist])
	    well_names.append(well_name)
	    #if the well name is not the same as the previous or if it is the first
	    #file, make a new directory
	    if well_name != well_names[num-1] or num == 0:
		num += 1
		if not os.path.exists(os.path.join(dir_path, well_name)):
		    os.makedirs(os.path.join(dir_path, well_name))
	    #move the current file to the new directory
	    logging.info('moving ./' + fname + ' to ./' + os.path.join(well_name, fname))
	    shutil.move(os.path.join(dir_path, fname), os.path.join(os.path.join(dir_path, well_name), fname))       
    logging.info('created well directories ' + str(set(well_names)))
	#Create a dir for each well and move all the matching files to it
    #    for well_name in set(well_names):
    #        if not os.path.exists(os.path.join(dir_path, well_name)):
    #            os.makedirs(os.path.join(dir_path, well_name))
    #            for fname in os.listdir(dir_path):
    #                if fname[-3:] in input_format:
    #                    well_ind = fname.index(well_str) + len(well_str)
    #                    if well_name == ''.join([str(int(char)) for char in fname[well_ind:well_ind+2] if char.isdigit()]):
    #                        shutil.move(os.path.join(dir_path, fname), os.path.join(os.path.join(dir_path, well_name), fname))
    #    
    return set(well_names)
        
#if I craete a dict with well names and the files, I can use the well name to create:
#the stitched file in the script level directory and use the well name to separate them instead
#come up with a new way to find the field name since the well string is now variable, either portion just 
#behind it or between it and the well string %w?       

def nat_key(key):
    '''
    A key to use with the `sorted()` function to sort naturally.
    '''
    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', key)]








# Runs main() if this script is called directly, but not if imported as a module from another script.
# Then the funcitons are just available in the namespace of tha tscript but nothing is executed.
if __name__ == '__main__':
    main()





#if args.recursive:
#    for dir_list in os.walk(dir_path):
#        print(dir_list[0])
#        for fname in os.listdir(dir_list[2]:
#            if fname[-3:] in input_format:
#                print(fname)
#                fnum = [int(char) for char in fname[33:35] if char.isdigit()][0]
#                #dunno why we need to flip...
#                imgs[fnum] = Image.open(os.path.join(dir_list[0], fname)).transpose(Image.FLIP_TOP_BOTTOM)
#else:
#    for fname in os.listdir(dir_path):
#        if fname[-3:] in input_format:
#            print(fname)
#            fnum = [int(char) for char in fname[33:35] if char.isdigit()][0]
#            #dunno why we need to flip...
#            imgs[fnum] = Image.open(os.path.join(dir_path, fname)).transpose(Image.FLIP_TOP_BOTTOM)
        


#could use this to create array of any size, but how to fill it in a spiral?
# array_dimension = int(np.ceil(np.sqrt(len(glob.glob(os.path.join(path, '*bmp'))))))
# img_layout = np.array([[3, 4, 5],
#                        [2, 1, 6],
#                        [9, 8, 7]]).ravel()


#set the layout
#img_layout = (13, 14, 15, 16, 17,
#              12,  3,  4,  5, 18,
#              11,  2,  1,  6, 19,
#              10,  9,  8,  7, 20,
#              25, 24, 23, 22, 21,)
#if args.layout_fields == 'right_down':
#    img_layout = img_layout[::-1]

