#!/usr/bin/env python

'''
Joel 2015-04-03

Small script for stitching together cellomics images. The user needs to define the path, 
input/output image extensions and the img_layout (fieldlookuptable for the microscope)
'''
from __future__ import division
import os
import Image
import argparse
import time
from itertools import cycle
import numpy as np

#Set up the command line arguments
#the formatter_class adds the default value to optional arguments
#argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#for positional arguments, one has to add the same string that the optional gets automatically
parser = argparse.ArgumentParser(description='A small utility for field images ' \
    'exported from automated microscope platforms.\nStarts in the middle and ' \
    'stitches fields in a spiral pattern to create the well image.\n' \
    'Each well and channel need to be in a separate directory', 
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('path', default=os.getcwd(), nargs='?',
    help='path to images  (default: current directory)')
parser.add_argument('-o', '--output-format', nargs='?', default='jpeg',
    help='format for the stitched image (default: %(default)s)')
parser.add_argument('-i', '--input-format', nargs='?', default='bmp',
    help='format for images to be stitched, can also be a list of formats (default: %(default)s)')
parser.add_argument('-f', '--field-string', default='_f',
    help='string immediately preceding the field number in the file name (default: %(default)s)')
parser.add_argument('-d', '--direction', nargs='?', default='left_up',
    help='directions from the 1st field to the 2nd and 3rd, e.g. left_up =\n' \
    '3, 4, 5, \n' \
    '2, 1, 6, \n' \
    '9, 8, 7 (default: %(default)s)')
parser.add_argument('-r', '--recursive', action='store_true',
    help='stitch images in subdirectories')
parser.add_argument('-n', '--no-flip', action='store_true',
    help='do not flip the input images vertically')
#Initialize some variables
args = parser.parse_args()
timestamp = str(int(time.time()))[3:]
img_types = set((args.input_format,)) #can add extra ext here is needed, remember to not have same as stiched
print('Stitching images...')

#Define movement function for filling in the sprial array
def move_right(x,y):
    return x, y +1


def move_down(x,y):
    return x+1,y


def move_left(x,y):
    return x,y -1


def move_up(x,y):
    return x -1,y

    
def spiral_structure(dir_path, img_types):
    '''
    Define the movement scheme and starting point for the field layout
    '''
    #find the number of fields/images matching the specified extension(s)
    fields = len([fname for fname in os.listdir(dir_path) if fname[-3:] in img_types])
    #size the array based on the field number, array will be squared
    arr_dim = int(np.ceil(np.sqrt(fields)))
    #define the movement schema and find the starting point (middle) of the array
    if args.direction == 'down_left':
        moves = [move_down, move_left, move_up, move_right]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:        
            starting_point = (int(arr_dim/2)-1, int(arr_dim/2))
    elif args.direction == 'left_down':
        moves = [move_left, move_down, move_right, move_up]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2)-1, int(arr_dim/2))
    elif args.direction == 'down_right':
        moves = [move_down, move_right, move_up, move_left]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2)-1, int(arr_dim/2)-1)
    elif args.direction == 'right_down':
        moves = [move_right, move_down, move_left, move_up]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2)-1, int(arr_dim/2)-1)
    elif args.direction == 'up_left':
        moves = [move_up, move_left, move_down, move_right]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
    elif args.direction == 'left_up':
        moves = [move_left, move_up, move_right, move_down]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
    elif args.direction == 'up_right':
        moves = [move_up, move_right, move_down, move_left]
        if arr_dim % 2 != 0:
            starting_point = (int(arr_dim/2), int(arr_dim/2))
        else:
            starting_point = (int(arr_dim/2), int(arr_dim/2)-1)
    elif args.direction == 'right_up':
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
        print('')
        print('Field layout:')    
        for row in np.ma.masked_equal(img_layout, -1): #TODO fix so that lines are showing for unused field
            print(' '.join(['{:>2}'.format(str(i)) for i in row]))

    else:
        for point, coord in list(gen_points(fields, moves, starting_point)):
            img_layout[coord] = point
        print('')
        print('Field layout:')    
        for row in np.ma.masked_equal(img_layout, -1):
            print(' '.join(['{:>2}'.format(str(i)) for i in row]))

    return img_layout
              

def find_images(dir_path, img_types, no_flip, field_str):
    '''
    Create a dictionary with the field numbers as keys to the field images
    '''
    zeroth_field = False #changes if a zeroth field is found in 'find_images'
    imgs = {}
    print('______________________________')
    print(dir_path)
    #go through each directory
    for fname in os.listdir(dir_path):
        if fname[-3:] in img_types:
            print(fname)
            #find the index of the field identifier string in the file name
            field_ind = fname.index(field_str) + len(field_str)
            fnum = int(''.join([str(int(char)) for char in fname[field_ind:field_ind+2] if char.isdigit()]))
            #If field 0 is encountered, change start numbering of the array
            if fnum == 0:
                zeroth_field = True
            #The default is to flip since this is the most common case
            if no_flip:
                imgs[fnum] = Image.open(os.path.join(dir_path, fname))
            else:
                #dunno why we need to flip...
                imgs[fnum] = Image.open(os.path.join(dir_path, fname)).transpose(Image.FLIP_TOP_BOTTOM) #transpose(Image.FLIP_LEFT_RIGHT)
                
    return imgs, zeroth_field

#stitch the image row by row
def stitch_images(imgs, img_layout, dir_path, output_format, arr_dim):
    '''
    Stitch images by going row and column wise in the img_layout and look up
    the number of the image to place at the (row, col) coordinate. So not filling
    in a spiral but using the spiral lookuptable instead.
    '''
    #create the size of the well image to be filled in
    width, height = imgs[1].size
    num = 0
    stitched = Image.new('RGB', (width*arr_dim, height*arr_dim))
    for row in xrange(0, width*arr_dim, width):
        for col in xrange(0, height*arr_dim, height):
            #since the image is filled by row and col instead of sprial, this
            #error catching is needed for the empty places
            try:
                #'stitch' fields by pasting them at the appropriate place in the black background
                stitched.paste(imgs[img_layout.ravel()[num]], (col, row))
            except KeyError:
                pass
            num += 1
    #save image
    stitched_name = os.path.join(dir_path, 'stitched_' + timestamp + '.' + output_format)
    stitched.save(stitched_name, format=output_format)
    # stitched.show()
    
    return stitched_name


#main program
if args.recursive:
    #walk through each subdirectory and the current direcotry
    for dir_list in os.walk(args.path):
        imgs, zeroth_field = find_images(dir_list[0], img_types, args.no_flip, args.field_string)
        #if there are images in the directory
        if imgs:
            fields, arr_dim, moves, starting_point = spiral_structure(dir_list[0], img_types)
            img_layout = spiral_array(fields, arr_dim, moves, starting_point, zeroth_field)
            stitched_name = stitch_images(imgs, img_layout, dir_list[0], args.output_format, arr_dim)
            print('Stitched image saved to ' + stitched_name + '\n')
        else:
            print('No images found')
else:
    imgs, zeroth_field = find_images(args.path, img_types, args.no_flip, args.field_string)
    #if there are images in the directory
    if imgs:
        fields, arr_dim, moves, starting_point = spiral_structure(args.path, img_types)
        img_layout = spiral_array(fields, arr_dim, moves, starting_point, zeroth_field)
        stitched_name = stitch_images(imgs, img_layout, args.path, args.output_format, arr_dim)
        print('Stitched image saved to ' + stitched_name + '\n')
    else:
        print('No images found')      

print('\nDone')

#if args.recursive:
#    for dir_list in os.walk(dir_path):
#        print(dir_list[0])
#        for fname in os.listdir(dir_list[2]:
#            if fname[-3:] in img_types:
#                print(fname)
#                fnum = [int(char) for char in fname[33:35] if char.isdigit()][0]
#                #dunno why we need to flip...
#                imgs[fnum] = Image.open(os.path.join(dir_list[0], fname)).transpose(Image.FLIP_TOP_BOTTOM)
#else:
#    for fname in os.listdir(dir_path):
#        if fname[-3:] in img_types:
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

