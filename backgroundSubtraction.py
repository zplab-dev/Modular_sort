# -*- coding: utf-8 -*-
"""
Created on Sun May 31 19:06:30 2015

@author: Willie
"""

import os
import freeimage
import numpy as np

from zplib.image import mask as zplib_image_mask

def subtract_one_time(background_model, focal_frame, ancillary_frames, lawn_mask):
	'''
	Does background subtraction for a single frame when given a background_model.
	'''
	focal_frame = focal_frame.copy()

	background_model = background_model.astype('int16')
	foreground_model = abs(focal_frame.astype('int16') - background_model.astype('int16')).astype('uint16')
	focal_mask = percentile_floor(foreground_model, threshold_proportion = 0.975)
	focal_mask[np.invert(lawn_mask)] = 0
	focal_mask = clean_dust_and_holes(focal_mask).astype('bool')

	ancillary_masks = []
	for ancillary_frame in ancillary_frames:
		foreground_frame = abs(ancillary_frame - background_model.astype('int16')).astype('uint16')
		ancillary_mask = percentile_floor(foreground_frame, threshold_proportion = 0.975)
		ancillary_mask = clean_dust_and_holes(ancillary_mask).astype('bool')
		ancillary_masks.append(ancillary_mask)
	
	return (focal_mask, ancillary_masks)

def background_frame(context_frames, focal_frame, ancillary_frames, lawn_mask, i):
	'''
	Does background subtraction for a single frame and prepares the context_frames for the next frame. This assumes that focal_frame is not yet included in context_frames.
	'''
	focal_frame = focal_frame.copy()
	context_frames.append(focal_frame)

	(foreground_file, background_file) = simple_running_median_subtraction(focal_frame, context_frames)
	background_file = background_file.astype('uint16')
	focal_mask = percentile_floor(foreground_file, threshold_proportion = 0.975)
	focal_mask[np.invert(lawn_mask)] = 0
	focal_mask= clean_dust_and_holes(focal_mask).astype('bool')

	if i > 15:
		focal_frame[focal_mask] = background_file[focal_mask]		

	ancillary_masks = []
	for ancillary_frame in ancillary_frames:
		foreground_frame = abs(ancillary_frame - background_file.astype('int16')).astype('uint16')
		ancillary_mask = percentile_floor(foreground_frame, threshold_proportion = 0.975)
		ancillary_mask = clean_dust_and_holes(ancillary_mask).astype('bool')
		ancillary_masks.append(ancillary_mask)
	
	context_frames = context_frames[1:]
	return (context_frames, background_file, focal_mask, ancillary_masks)

def overallBackgroundSubtract(data_dir, match_string, new_string, temporal_radius):
	'''
	Do background subtraction to find worms. This uses only past data, masking out the worms to create a background that won't disappear once the worm stops  moving.
	'''
	my_files = sorted(os.listdir(data_dir))
	my_files = [a_file for a_file in my_files if match_string in a_file]
	my_times = [a_file.split(os.path.sep)[-1].split(' ')[0] for a_file in my_files]
	ending_list = [' bf00.png', ' bf01.png', ' bf10.png', ' bf11.png']

	# Run the actual simple subtraction, saving out masked files.
	context_files = [freeimage.read(data_dir + os.path.sep + my_files[j]) for j in range(0, temporal_radius)]
	for i in range(temporal_radius, len(my_files)):
		real_raw_file = freeimage.read(data_dir + os.path.sep + my_files[i])
		raw_file = real_raw_file.copy()		
		context_files.append(raw_file)

		(foreground_file, background_file) = simple_running_median_subtraction(raw_file, context_files)
		thresholded_mask = percentile_floor(foreground_file, threshold_proportion = 0.975)
		final_mask = clean_dust_and_holes(thresholded_mask)

		raw_file[final_mask.astype('bool')] = background_file[final_mask.astype('bool')]		
		freeimage.write(final_mask, data_dir + os.path.sep + my_files[i].replace(' ', ' mask_'))
		freeimage.write(background_file, data_dir + os.path.sep + my_files[i].replace(match_string, 'background.png'))	

		for an_ending in ending_list:
			my_path = data_dir + os.path.sep + my_times[i] + an_ending
			save_path = data_dir + os.path.sep + my_times[i] + an_ending.replace(' ', ' mask_')
			if os.path.isfile(my_path):
				my_image = freeimage.read(my_path)
				foreground_file = abs(my_image.astype('int16') - background_file.astype('int16'))
				foreground_file = foreground_file.astype('uint16')
				thresholded_mask = percentile_floor(foreground_file, threshold_proportion = 0.975)
				final_mask = clean_dust_and_holes(thresholded_mask)
				freeimage.write(final_mask, save_path)
		
		context_files = context_files[1:]

	# Run another small chunk of background subtraction backwards to fill out the early range.
	context_files = [freeimage.read(data_dir + os.path.sep + my_files[j]) for j in reversed(range(temporal_radius, temporal_radius*2))]
	for i in reversed(range(0, temporal_radius)):
		real_raw_file = freeimage.read(data_dir + os.path.sep + my_files[i])
		raw_file = real_raw_file.copy()		
		context_files.append(raw_file)

		(foreground_file, background_file) = simple_running_median_subtraction(raw_file, context_files)
		thresholded_mask = percentile_floor(foreground_file, threshold_proportion = 0.975)
		final_mask = clean_dust_and_holes(thresholded_mask)

		raw_file[final_mask.astype('bool')] = background_file[final_mask.astype('bool')]		
		freeimage.write(final_mask, data_dir + os.path.sep + my_files[i].replace(' ', ' mask_'))
		freeimage.write(background_file, data_dir + os.path.sep + my_files[i].replace(match_string, 'background.png'))	
		
		for an_ending in ending_list:
			if os.path.isfile(data_dir + os.path.sep + my_times[i] + an_ending):				
				my_image = freeimage.read(my_path)
				foreground_file = abs(my_image.astype('int16') - background_file.astype('int16'))
				foreground_file = foreground_file.astype('uint16')
				thresholded_mask = percentile_floor(foreground_file, threshold_proportion = 0.975)
				final_mask = clean_dust_and_holes(thresholded_mask)
				freeimage.write(final_mask, save_path)
	return

def clean_dust_and_holes(dusty_pic):
	'''
	Picks out the largest object in dusty_mask and fills in its holes, returning cleaned_mask.
	'''
	my_dtype = dusty_pic.dtype
	dust_mask = np.invert(zplib_image_mask.get_largest_object(dusty_pic))		
	dusty_pic[dust_mask] = 0
	cleaned_mask = simple_floor(zplib_image_mask.fill_small_area_holes(dusty_pic, 90000).astype(my_dtype), 1)
	return cleaned_mask

def simple_floor(focal_image, threshold_value):
	'''	
	Takes a grayscale focal image (in the form of a numpy array), and sets to zero (black) all values below the threshold value.
	'''
	max_value = -1
	binary_image = focal_image.copy()
	binary_image[binary_image < threshold_value] = 0 
	binary_image[binary_image >= threshold_value] = max_value 
	return binary_image
	
def percentile_floor(focal_image, threshold_proportion):
	'''
	Takes a grayscale focal image (in the form of a numpy array), and sets to zero (black) all values below the percentile indicated by threshold_proportion.
	'''
	max_value = -1
	binary_image = focal_image.copy()
	threshold_value = int(np.percentile(binary_image, threshold_proportion*100))
	binary_image[binary_image < threshold_value] = 0 
	binary_image[binary_image >= threshold_value] = max_value 	
	return binary_image

def simple_running_median_subtraction(focal_image, background_images):
	'''	
	Takes a focal image and a list of background images (grayscale, in the form of numpy arrays), and returns the focal image with the background subtracted. This simply takes the median value of each pixel to construct a background.	
	'''
	median_image = median_image_from_list(background_images)
	foreground_only = abs(focal_image.astype('int16') - median_image.astype('int16')).astype('uint16')
	return (foreground_only, median_image)

def median_image_from_list(background_images):
	'''	
	Takes a list of background images (grayscale, in the form of numpy arrays), and returns an image constructed by taking the median value of each pixel.
	'''
	big_array = np.array(background_images)
	median_image = np.median(big_array, axis = 0).astype('uint16')
	return median_image


def main():
	return

if __name__ == "__main__":
	main()
