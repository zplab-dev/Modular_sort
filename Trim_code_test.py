#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created March 2018
@author: Matt Mosley, m.mosley@wustl.edu
Edited version of Tim's/Nick's Modular Sort code, hopefully more lightweight.
Classes to sort on GFP, autofluorescence, length.
"""

import iotool
from scope import scope_client
from pathlib import Path
from scipy import ndimage
import numpy
import time
import freeimage
import threading
import csv
import requests
import collections
from skimage.filters import sobel
from skimage.morphology import watershed
from skimage.measure import label

import backgroundSubtraction
from ris_widget import ris_widget
import zplib.image.write_movie
import zplib.image.colorize

#import scope_controls


IMAGE_SIZE = (1280, 1080)
BOILER_AREA = (slice(530,1100), slice(530,590))
DETECTION_AREA = (slice(690,1100), slice(545, 575))
POSITION_AREA = (slice(690,1000), slice(545, 575))  #Align left end of sewer corner to 690, 575
CLEARING_AREA = (slice(530,1000), slice(530, 580))
FLUORESCENT_AREA = (slice(500,1100), slice(530, 590))

#Textbelt key = 08a7e7d3335d92542dfec857461cfb14af5e0805HQINVtRpVqvDjo9O2wa2I6tTo

FLUOR_PIXEL_BRIGHT_VALUE = 500		#Used for finding brightness, newer/better way to do this?
DETECTION_THRES = 4
POSITION_THRES = 3  #Play with this? Or if it ain't broke don't fix it?
CLEARING_THRES = 3
LOST_CUTOFF = 1.5
DOUBLE_THRESH = 1.3

CYAN_EXPOSURE_TIME = 7  #10 for lin-4, 50 for mir-71?
GREEN_YELLOW_EXPOSURE_TIME = 50	#for autofluorescence

BRIGHT_FIELD_EXPOSURE_TIME = 2

PICTURE_DELAY = .01
SORTING_INTERVAL = .4

BACKGROUND_REFRESH_RATE = 100000	#Not actually currently getting called, could use this or reset by number of worms
PROGRESS_RATE = 100

MIN_GFP_THRESH = 400 #Will need to reset to account for brighter exposure
                     #TODO: Set this value based on non-fluorescent worms

#ioTool language
PUSH_CHANNEL_PRESSURE= 'sh D6'
PUSH_CHANNEL_STATIC = 'sl D6'
SEWER_CHANNEL_PRESSURE = 'sh D7'
SEWER_CHANNEL_SUCK = 'sl D7'
UP_CHANNEL_SUCK = 'sl D3'
UP_CHANNEL_PRESSURE = 'sh D3'
STRAIGHT_CHANNEL_SUCK = 'sl D4'
STRAIGHT_CHANNEL_PRESSURE = 'sh D4'
DOWN_CHANNEL_SUCK = 'sl D5'
DOWN_CHANNEL_PRESSURE = 'sh D5'

def boiler():		#TODO: Work out better way to detect worms than having set detection windows based on boolean arrays
    """
    Returns a boolean array of possible locations of a worm duing sorting
    """
    boiler = numpy.ones(IMAGE_SIZE, dtype=bool)
    boiler[BOILER_AREA] = False
    return boiler

class MicroDevice:
    """Superclass which contains basic functions for running the device, saving images,
    building the calibration histogram, etc. Sorting itself should be run from a subclass e.g.:
    test = Trim_code_test.GFP('save/data/here')
    test.run()
    """

    def __init__(self, exp_direct, af_filter=True):
        """Initalizes the scope and device
        """
        #Initialize scope:
        self.scope = scope_client.ScopeClient()
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
        self.scope.camera.readout_rate = '280 MHz'
        self.scope.camera.binning = '2x2'
        self.scope.tl.lamp.enabled = True
        self.scope.il.spectra.cyan.enabled = False
        self.scope.il.spectra.green_yellow.enabled = False

        #Initialize device:
        self.device = iotool.IOTool("/dev/ttyMicrofluidics")

        #Create folder for data:
        self.file_location = Path(exp_direct)
        self.file_location.mkdir(mode=0o777, parents=True, exist_ok=True)
        self.info = exp_direct.split('/')[-1]
        self.setup_dirs()

        #Make sure device is clear:
        self.device_clear_tubes()
        if af_filter:
            self.af_filter = True
        else:
            self.af_filter = False

    def write_csv_line(self,csv,data):
        csv.write(','.join(map(str, data)) + '\n')

    def lamp_off(self):
        self.scope.il.spectra.cyan.enabled = False
        self.scope.tl.lamp.enabled = False
        self.scope.il.spectra.green_yellow.enabled = False
        self.scope.il.spectra.red.enabled = False

    def bright(self):
        self.scope.tl.lamp.enabled = True
        self.scope.il.spectra.cyan.enabled = False
        self.scope.il.spectra.green_yellow.enabled = False
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
        time.sleep(PICTURE_DELAY)

    def cyan(self):
        self.scope.tl.lamp.enabled = False
        self.scope.il.spectra.cyan.enabled = True
        self.scope.il.spectra.green_yellow.enabled = False
        self.scope.camera.exposure_time = CYAN_EXPOSURE_TIME
        time.sleep(PICTURE_DELAY)

    def green_yellow(self):
        self.scope.tl.lamp.enabled = False
        self.scope.il.spectra.cyan.enabled = False
        self.scope.il.spectra.green_yellow.enabled = True
        self.scope.camera.exposure_time = GREEN_YELLOW_EXPOSURE_TIME
        time.sleep(PICTURE_DELAY)

    def capture_image(self, type_of_image):
        type_of_image()
        self.scope.camera.send_software_trigger()
        return self.scope.camera.next_image()

    def device_stop_run(self):
        self.scope.camera.end_image_sequence_acquisition()
        self.device.execute(PUSH_CHANNEL_PRESSURE, SEWER_CHANNEL_PRESSURE,
                            UP_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_PRESSURE,
                            DOWN_CHANNEL_PRESSURE)

    def device_start_load(self):
        self.device.execute(PUSH_CHANNEL_STATIC, SEWER_CHANNEL_SUCK,
                            UP_CHANNEL_PRESSURE,STRAIGHT_CHANNEL_PRESSURE,
                            DOWN_CHANNEL_PRESSURE)

    def device_stop_load(self):
        self.device.execute(PUSH_CHANNEL_PRESSURE)

    def device_clear_tubes(self):
        """Sets all input channels to push.
        """
        self.device.execute(SEWER_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_PRESSURE,
                            UP_CHANNEL_PRESSURE, DOWN_CHANNEL_PRESSURE,
                            PUSH_CHANNEL_PRESSURE)

    def save_image(self, image, name, worm_count, _type=None):
        """Saves image to directory provided when script is first loaded.
           If type is specified, saves to subdirectory.
           types: bf (brightfield), cyan (GFP), tritc (autofluor), dead, small,
           big, and bent.
           """
        if _type is not None:   #Type given -> save to subdirectory
            save_location = str(self.dir_dict[_type]) + '/' + name + '_' + str(worm_count) + '.png'
            freeimage.write(image, save_location, flags = freeimage.IO_FLAGS.PNG_Z_BEST_SPEED)
        else:   #Otherwise just save to sort folder
            save_location = str(self.file_location) + '/' + name + '_' + str(worm_count) + '.png'
            freeimage.write(image, save_location, flags=freeimage.IO_FLAGS.PNG_Z_BEST_SPEED)

    def setup_dirs(self):
        self.bf_dir = self.file_location.joinpath('brightfield_images')
        self.mask_dir = self.file_location.joinpath('mask_images')
        self.cyan_dir = self.file_location.joinpath('GFP_images')
        self.tritc_dir = self.file_location.joinpath('autofluor_images')
        self.dead_dir = self.file_location.joinpath('dead_imgaes')
        self.smol_dir = self.file_location.joinpath('small_worm_images')
        self.big_dir = self.file_location.joinpath('double_worm_images')
        self.bent_dir = self.file_location.joinpath('bent_worm_images')
        self.dir_dict = {'bf' : self.bf_dir,
                         'mask' : self.mask_dir,
                         'cyan' : self.cyan_dir,
                         'tritc' : self.tritc_dir,
                         'dead' : self.dead_dir,
                         'small' : self.smol_dir,
                         'big' : self.big_dir,
                         'bent': self.bent_dir}

        for key in self.dir_dict:
            self.dir_dict[key].mkdir(mode=0o777, parents=True, exist_ok=True)

    def set_background_areas(self, worm_count):
        self.set_bf_background(worm_count)
        self.lamp_off()
        self.set_cyan_background(worm_count)
        self.lamp_off()
        self.set_green_yellow_background(worm_count)
        self.lamp_off()
        self.scope.tl.lamp.enabled = True
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME

    def set_bf_background(self, worm_count):
        """Sets background for brightfield worm detection.
        """
        self.background = self.capture_image(self.bright)
        time.sleep(PICTURE_DELAY)
        self.background_2 = self.capture_image(self.bright)
        subtracted = abs(self.background.astype('int32') - self.background_2.astype('int32'))
        self.detect_background = numpy.sum(subtracted[DETECTION_AREA])
        self.positioned_background = numpy.sum(subtracted[POSITION_AREA])
        self.clear_background = numpy.sum(subtracted[CLEARING_AREA])
        self.save_image(self.background, 'brightfield_background', worm_count, _type = 'bf')

    def set_cyan_background(self, worm_count):
        self.cyan_background = self.capture_image(self.cyan)
        self.save_image(self.cyan_background, 'cyan_background', worm_count, _type = 'cyan')

    def set_green_yellow_background(self, worm_count):
        self.green_yellow_background = self.capture_image(self.green_yellow)
        self.save_image(self.green_yellow_background, 'green_yellow_background', worm_count, _type = 'tritc')

    def reset_background(self, worm_count):
        self.set_background_areas(worm_count)
        print('resetting backgrounds')
        time.sleep(1)

    def detect_worms(self, current_image, background):      #TODO: print out what the actual detect_background is?
        """
        Function that returns if a worm has been detected
        """
        #print('Detected Value:' + str(numpy.abs(numpy.sum(current_image[DETECTION_AREA].astype('int32') - background[DETECTION_AREA].astype('int32')))))
        #print('Required Value:' + str( 5 * self.detect_background))
        #TODO: Play with sensitivity using this?

        return ((numpy.sum(numpy.abs(current_image[DETECTION_AREA].astype('int32')
                - background[DETECTION_AREA].astype('int32')))
                -self.detect_background) > DETECTION_THRES  * self.detect_background)

    def positioned_worm(self, current_image, detected_image):
        """
        Function that determines if a worm has been positioned
        The worm is positioned because the change is small.
        """
        worm_movment = abs(current_image[POSITION_AREA].astype('int32') - detected_image[POSITION_AREA].astype('int32'))
        return  ((numpy.sum(worm_movment) - self.positioned_background) < POSITION_THRES * self.positioned_background)

    def lost_worm(self, current_image, background):     #Probably a better way to do this
        """
        Function that determines if a worm was lost
        Worm is deciced lost because image is close enough to background.
        """
        worm_visibility = abs(current_image[POSITION_AREA].astype('int32')
                              - background[POSITION_AREA].astype('int32'))
        return ((numpy.sum(worm_visibility) -
                 self.positioned_background) < LOST_CUTOFF * self.positioned_background)

    def size_of_worm(self, subtracted_image):
        """
        Function that saves an image of the mask of a worm and returns the size of mask (in number of pixels)
        Depending on how background.backgroundSubtraction.clean_dust_and_holes(image) works this function might neeed to be modified so that it returns the appropriate image.
        Currently the clean_dust_and_holes does not return the clean image and actually modifies the passed image.
        """
        mask = self.worm_mask(subtracted_image)
        a = (1*numpy.sum(mask.astype('bool')))
        return a

    def check_size_change(self, size_before, size_after):
        #Construct new worm_mask
        #Compare with previous size
        #Should be within ~10% of original size, otherwise new worm may have shown up (or worm may have been lost)
        #Will this mess me up with worms that get imaged before being fully in frame? Need to also implement edge-closeness rule

        percent_change = abs(size_after - size_before)/size_before
        if percent_change > 0.2 * size_before:
            print('Percent change = ' + str(percent_change))
            return True
        else:
            print('Percent change = ' + str(percent_change))
            pass

    def worm_mask(self, subtracted_image):      #Make this more consistent/make everything go through this
        floored_image = backgroundSubtraction.percentile_floor(subtracted_image, .99)
        floored_image[self.boiler] = 0
        backgroundSubtraction.clean_dust_and_holes(floored_image)

        #Cleans mask/removes jagged lines from sorter background
        for i in range(3):
            floored_image = ndimage.binary_erosion(floored_image)
        for i in range(3):
            floored_image = ndimage.binary_dilation(floored_image)

        return floored_image

    def find_dimensions(self, mask):
        """Takes boolean mask (cleaned) and returns length and width, treating the worm roughly as a box
        """
        xs, ys = ndimage.find_objects(mask)[0]
        cropped_image = mask[xs, ys]

        length, width = cropped_image.shape
        print('Length = ' + str(length))
        print('aspect ratio = ' + str(length/width))

        #How to make it wait longer if a worm hasn't finished entering the viewing area? Makes sense to gate
        #by how far the x axis extends? Need to return an extra variable to check for that?
        return length, width

    def analyze(self):
        """
        function that tells the device what sorting/analzying method to use:
        Is overwritten by a super class
        """
        raise NotImplementedError('No sorting method given')

    def device_sort(self, direction, background, worm_count):
        """
        Command that toggles the device to sort a worm in given direction
        Have yet to test if having teh pusher channel be set to static when sorting to help worms load faster.

        """
        if direction == 'up':
            self.device.execute(SEWER_CHANNEL_PRESSURE, UP_CHANNEL_SUCK ,PUSH_CHANNEL_STATIC)
            if self.check_cleared(background, worm_count):
                time.sleep(SORTING_INTERVAL)
                self.device.execute(UP_CHANNEL_PRESSURE)
                self.up_count += 1
        elif direction == 'down':
            self.device.execute(SEWER_CHANNEL_PRESSURE, DOWN_CHANNEL_SUCK, PUSH_CHANNEL_STATIC)
            if self.check_cleared(background, worm_count):
                time.sleep(SORTING_INTERVAL)
                self.device.execute(DOWN_CHANNEL_PRESSURE)
                self.down_count += 1
        elif direction == 'straight':
            self.device.execute(SEWER_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_SUCK, PUSH_CHANNEL_STATIC)
            if self.check_cleared(background, worm_count):
                time.sleep(SORTING_INTERVAL)
                self.device.execute(STRAIGHT_CHANNEL_PRESSURE)
                #TODO: I don't actually want to iterate this unless it's a "good" worm, currently going up for all worms sorted straight
                self.straight_count += 1
        else:
            raise Exception('Direction "' + str(direction) + '" not recognized')

    def check_cleared(self, background, worm_count):        #TODO: Fix force reset option, because flag never triggers
                                                            #how to make this better? Reject worm after too long?
        time_clear_start = time.time()
        while True:
            current_image = self.capture_image(self.bright)
            sorted_worm_difference = abs(current_image[CLEARING_AREA].astype('int32') - background[CLEARING_AREA].astype('int32'))
            if numpy.sum(sorted_worm_difference) < CLEARING_THRES * self.clear_background:
                print('Worm determined cleared')
                return True
            else:
                self.device.execute(PUSH_CHANNEL_PRESSURE)
                self.device.execute(SEWER_CHANNEL_PRESSURE)
                print('still pushing')
                time_stuck = time.time()
                if time_stuck - time_clear_start > 60:
                    self.send_text(message='Sort stuck trying to clear')
                    print('Sent text: Device stuck, reseting bg')
                    message_sent = True
                    self.device.execute(UP_CHANNEL_SUCK, STRAIGHT_CHANNEL_SUCK, DOWN_CHANNEL_SUCK)
                    time.sleep(1)
                    self.reset_background(worm_count)
                    return True

    def clear_double_worms(self):
        self.device.execute(SEWER_CHANNEL_PRESSURE,
                    UP_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_SUCK,
                    DOWN_CHANNEL_PRESSURE)
        time.sleep(0.5)


    def send_text(self, message):
        """
        TODO:replace this with discord web hook

        requests.post('https://textbelt.com/text',
                        {'phone': '6019538192',
                         'message': str(message),
                         'key':'08a7e7d3335d92542dfec857461cfb14af5e0805HQINVtRpVqvDjo9O2wa2I6tTo'})
        """
        pass

    #def main(self):    #Not entirely sure what this is
    #    self.device = iotool.IOTool("/dev/ttyMicrofluidics")
    #    return

    def acquire_worm_images(self):
        """Captures bright field and fluorescent images.
        """

        bf_image = self.capture_image(self.bright)
        time.sleep(PICTURE_DELAY)
        cyan_image = self.capture_image(self.cyan)
        time.sleep(PICTURE_DELAY)
        tritc_image = self.capture_image(self.green_yellow)
        #time.sleep(PICTURE_DELAY)
        self.bright()

        return bf_image, cyan_image, tritc_image

    def check_dead(self, cyan_subtracted, mask):
       """Compares 95th percentile with median (or mean) fluorescence, in order to determine if
       animal has full wave of death fluorescence.
       """
       fluor_95th = self.find_95th_fluor_amount(cyan_subtracted, mask)
       fluor_median = self.find_median_fluor_amount(cyan_subtracted, mask)

       if 2*fluor_median >= fluor_95th:
           return True
       else:
           pass

    def check_aspect_ratio(self, length, width):
        """Overwritten if sorting by worm length
        """
        pass

    def find_95th_fluor_amount(self, subtracted_image, worm_mask):
        fluor_worm = subtracted_image[worm_mask]
        fluor_95th = numpy.percentile(fluor_worm, 95)
        return fluor_95th
    
    def fluor_over_99_mean(self, subtracted_image, worm_mask):
        fluor_worm = subtracted_image[worm_mask]
        fluor_99 = numpy.percentile(fluor_worm, 99)
        pix_over_99 = [pixel for pixel in fluor_worm if pixel >= fluor_99]
        over_99_mean = numpy.mean(pix_over_99)
        return over_99_mean

    def find_mean_fluor_amount(self, subtracted_image, worm_mask):
        fluor_worm = subtracted_image[worm_mask]
        mean_fluor = numpy.mean(fluor_worm)
        return mean_fluor

    def find_median_fluor_amount(self, subtracted_image, worm_mask):
        fluor_worm = subtracted_image[worm_mask]
        fluor_median = numpy.median(fluor_worm)
        return fluor_median

    def build_hist(self, size = 100):
        """Function that builds initial histogram from which percentiles are taken to determine upper and lower thresholds
        for sorting. Also returns sorting parameter (fluorescence, length, etc) as list which can be updated.
        """

        #TODO: would it not make more sense to save each parameter (fluor, size, etc) as an array and then put them all
        #together in the summary file? Then I wouldn't need a separate list for hist values, would just pull from all
        #fluors taken

        self.sort(calibration = True, initial_hist_size = size)

        self.upper_threshold = numpy.percentile(self.hist_values, 90)
        self.lower_threshold = numpy.percentile(self.hist_values, 10)

        print('Value list = \n ' , str(self.hist_values))
        print('Avg =' + str(numpy.mean(self.hist_values)))
        print('90th percentile =' + str(self.upper_threshold))
        print('10th percentile =' + str(self.lower_threshold))

        return self.hist_values, self.upper_threshold, self.lower_threshold

    def update_hist(self, sort_param):
        if type(sort_param) != str:
            self.hist_values.append(sort_param)
            self.upper_threshold = numpy.percentile(self.hist_values, 90)
            self.lower_threshold = numpy.percentile(self.hist_values, 10)
            print(self.upper_threshold, self.lower_threshold)
        else:
            pass

    def check_metrics(self, current_image, worm_count, size1, size2, cyan_subtracted, worm_mask, cyan_image, af):
        """Checks all metrics that would cause a worm to be rejected in the sort,
        i.e. size, death fluorescence, aspect ratio (if sorting by length).
        How to include 2x size check to make sure worm didn't change too much
        (new worm showing up)? (Need to check images for new worms showing up to see if that's working)

        Returns note, which explains if worm was deemed good for sorting or explains
        reason why worm was rejected.
        """
        length, width = self.find_dimensions(worm_mask)
        
        #Autofluorescence (should be within an expected range for 6dph or older worms, if sorting 4dph or long-lived mutants, may need to turn off)
        if self.af_filter:
            if af <= 100:
                print('Autofluorescence suspiciously low')
                return 'low_af'

        #Worm size when imaged:
        if size1 > self.size_threshold:     #TODO: think on how to better decide size thresholds
            print('Detected double worm')
            self.save_image(current_image, 'doubled_worm_analyze', worm_count, _type = 'big')
            return 'doubled_worm'
        elif size1 < self.min_worm_size:     #TODO: Make sure I'm not rejecting slow, big worms
            print('Detected small worm')
            self.save_image(current_image, 'small_worm_analyze', worm_count, _type = 'small')
            return 'small'

        #Death fluorescence:    TODO: How often is this getting called and is it correct?
        elif self.check_dead(cyan_subtracted, worm_mask):
            print('Worm ' + str(worm_count) + ' determined dead')
            self.save_image(current_image, 'dead_worm', worm_count, _type = 'dead')
            self.save_image(cyan_image, 'dead_worm_cyan', worm_count, _type = 'dead')
            return 'dead'

        #Aspect ratio
        elif self.check_aspect_ratio(length, width) == 'bent':
            #Function should only be called if running from length class
            print('Worm bent over')
            self.save_image(current_image, 'bent_worm', worm_count, _type = 'bent')
            return 'bent'

        #Size change (New worm appeared or old worm dissapeared)
        elif self.check_size_change(size1, size2):
            print('Detected appreciable size change')
            print('****Feature was useful****')
            self.save_image(post_analysis_image, 'size_difference_analyze', worm_count, _type = 'big')
            #Filing these with worms above size threshold for now because i'm not sure how often this gets used
            return 'new_worm'

        #Passed all?
        else:
            return 'sort'
    
    def check_note(self, calibration, note, direction):
        
        if note != 'sort':
            direction = 'straight'
            self.bad_worm_count += 1
            #sort straight if worm is bad for some reason
        elif calibration and note == 'sort':
            direction = 'straight'
            note = 'calibration' #note worms used in initial histogram as calibration worms
        elif not calibration and note == 'sort':
            #if worm is deemed good and sort is running, note is sort
            pass
        
        return note, direction


    def sort(self, calibration, initial_hist_size=100):
        """
        Function that starts the device running with a given purpose
        #0 set background
        #1 load worm
        #2 detect worm
        #3 stop worms
        #4 position worms
        #5 picture worms
        #6 analze worms
        #7 sort worms
        #8 move worms
        #9 --> 1
        """

        self.up_count = 0
        self.down_count = 0
        self.straight_count = 0
        self.bad_worm_count = 0
        cycle_count = 0

        test_counter = 0

        #1 Loading Worms
        self.device_start_load()

        time_seen = self.time_start #Initial setting for time_seen before worm found
        message_sent = False

        self.film_reel = collections.deque(maxlen = 200)

        #rw = ris_widget.RisWidget() <was trying to play around with mask visualization

        #2 Detect Worms
        try:
            while True:
                cycle_count += 1

                if calibration:
                    if len(self.hist_values) >= initial_hist_size:
                        print('Finished initial histogram')
                        return self.hist_values
                        break

                time_without_worm = time.time()
                if time_without_worm - time_seen > 180 and message_sent == False:
                    self.send_text(message='No worm for 3 min')
                    print('Warning sent: No worm')
                    message_sent = True

                current_image = self.capture_image(self.bright)

                if self.detect_worms(current_image, self.background):
                    print(' ')
                    print(' ')
                    print('Worm has been detected')
                    time_seen = time.time()
                    time_between_worms = time_seen - self.time_start
                    print('Time since start: ' + str(round(time_between_worms / 60, 3))+ ' min')
                    detected_image = current_image
                    self.worm_count += 1
                    if message_sent == True:
                        self.send_text(message = 'Worms are back')
                        message_sent = False

                    #3 Stop worms
                    self.device_stop_load()

                    #4 Position worms
                    while True:
                        current_image = self.capture_image(self.bright)
                        if self.lost_worm(current_image, self.background):  #This doesnt seem to get called often
                            print('Worm was lost')
                            self.device.execute(SEWER_CHANNEL_PRESSURE)
                            time.sleep(.1)
                            note = 'lost'
                            self.worm_count -= 1
                            break

                        elif self.positioned_worm(current_image, detected_image):
                            #Maybe include some param that says if a worm is too close to the edge, let it get to the center

                            bf_image, cyan_image, tritc_image = self.acquire_worm_images()
                            bf_subtracted = (abs(bf_image.astype('int32') - self.background.astype('int32')))
                            worm_mask = self.worm_mask(bf_subtracted)
                            self.scale_image(worm_mask, 'mask')
                            size1 = self.size_of_worm(bf_subtracted)
                            #Used as worm size but also for comparison later
                            print('Size of worm before sorting: ' + str(size1))

                            if size1 == 0: #Dumb way of fixing empty array problem, since size thresholds being checked later
                                print('No worm')
                                note = 'lost'
                                self.worm_count -= 1
                                break

                            #rw.image = worm_mask <riswidget seems to just bug out, how to get mask to show up?

                            #numpy.clip gives limits for the output of an operation, so this avoids negative values and "bounce" from taking an absolute value.
                            #Having this standard across experiments is likely fairly important. For instance, 2018/4/26 autofluor sort was before this when absolute value was used.
                            #This was updated 2018/4/30
                            cyan_subtracted = numpy.clip(cyan_image.astype('int32') - self.cyan_background.astype('int32'), 0, 100000)
                            tritc_subtracted = numpy.clip(tritc_image.astype('int32') - self.green_yellow_background.astype('int32'), 0, 100000)

                            autofluorescence = self.find_95th_fluor_amount(tritc_subtracted, worm_mask)
                            print('Autofluorescence = ' + str(autofluorescence))

                            print('Worm positioned')

                            if calibration:
                                self.save_image(current_image, 'calibration_brightfield_worm', self.worm_count, _type = 'bf')
                                self.save_image(worm_mask.astype('uint8')*255, 'calibration_worm_mask', self.worm_count, _type = 'mask')
                                self.save_image(cyan_image, 'calibration_cyan_worm', self.worm_count, _type = 'cyan')
                                self.save_image(tritc_image, 'calibration_tritc_worm', self.worm_count, _type = 'tritc')
                                print('Images saved')

                                sort_param, direction = self.analyze(cyan_subtracted, tritc_subtracted, worm_mask, self.worm_count, calibration = True)
                                print('Calibration worm number: ' + str(self.worm_count))

                            elif not calibration:
                                #5 Save image
                                self.save_image(current_image, 'brightfield_worm', self.worm_count, _type = 'bf')    #TODO: better names?
                                self.save_image(worm_mask.astype('uint8')*255, 'worm_mask', self.worm_count, _type = 'mask')
                                self.save_image(cyan_image, 'cyan_worm', self.worm_count, _type = 'cyan')
                                self.save_image(tritc_image, 'tritc_worm', self.worm_count, _type = 'tritc')
                                print('Images saved')

                                sort_param, direction = self.analyze(cyan_subtracted, tritc_subtracted, worm_mask, self.worm_count, calibration = False)
                                print('Worm number: ' + str(self.worm_count))

                            #Second size change to make sure new worm hasn't shown up:
                            post_analysis_image = self.capture_image(self.bright)
                            difference_between_worm_background = (abs(post_analysis_image.astype('int32') - self.background.astype('int32')))
                            size2 = self.size_of_worm(difference_between_worm_background)

                            print('Size of worm after analysis: ' + str(size2))

                            note = self.check_metrics(current_image, self.worm_count, size1, size2, cyan_subtracted, worm_mask, cyan_image, autofluorescence)

                            note, direction = self.check_note(calibration, note, direction)

                            print('Note = ' + str(note))
                            
                            if note == 'calibration':
                            
                                self.hist_values.append(sort_param)  #building initial histogram
                                print('Good calibration worms: ' + str(len(self.hist_values)))
                                
                            if note == 'sort':
                                self.update_hist(sort_param)

                            #7 Sort worms
                            self.device_sort(direction, self.background, self.worm_count)
                            print('Worm ' + str(self.worm_count) + ' sorted ' + direction)

                            if not calibration:
                                print('Up: ' + str(self.up_count), ' Straight (sorted): ' + str(self.straight_count-self.bad_worm_count), ' Straight (bad): ' + str(self.bad_worm_count), ' Down: ' + str(self.down_count))
                            
                            break   #Neccessary to break while loop

                        else:
                            detected_image = current_image  #What is this doing? Pass back to front?

                    if note != 'lost':
                        #No need to save data for lost worms, and I don't think we were losing many anyway
                        worm_data = self.generate_data(self.worm_count, size1, autofluorescence, sort_param, time_between_worms, direction, note)
                        #TODO: Fix data generation across all classes (if sort_param == list etc)
                        self.write_csv_line(self.summary_csv, worm_data)
                        print('Worm ' + str(self.worm_count) + ' data saved')

                    print('_______________________________________________________________________________')

                    self.make_movie(self.worm_count, self.film_reel)

                    if self.worm_count % 100 == 0 and self.worm_count != 0:   #Resets background every 100 worms
                        self.reset_background(self.worm_count)
                        #Check background pics to make sure this isn't ever taking pics with worms

                    self.device_start_load()
                    #9 --> 1

                elif cycle_count % PROGRESS_RATE == 0:  #Gives update every 100 cycles to let you know it's still working
                    print(str(PROGRESS_RATE) + ' Cycles')

        except KeyboardInterrupt:
            pass

        finally:
            self.device_stop_run()
            print('finally')
            #Close csv file?

    def run(self, hist = []):

        self.setup_csv(self.file_location, self.info)

        self.scope.camera.start_image_sequence_acquisition(frame_count=None, trigger_mode='Software')

        self.boiler = boiler()

        self.worm_count = 0
        self.time_start = time.time()

        #0 Setting Background
        self.set_background_areas(self.worm_count)
        print('setting backgrounds')
        time.sleep(1)

        self.size_threshold = 5500   #hard coding sizes for now
        self.min_worm_size = 2300
        
        self.hist_values = hist
        self.build_hist(size=100)

        self.scope.camera.start_image_sequence_acquisition(frame_count=None, trigger_mode='Software')

        self.sort(calibration = False)

        self.summary_csv.close()

    def manual_set_up(self):
        self.hist_values = list(input('Histogram values: '))
        self.upper_threshold = numpy.percentile(self.hist_values, 90)
        self.lower_threshold = numpy.percentile(self.hist_values, 10)

    def scale_image(self, image, type_of_image):
        pass

    def make_movie(self, worm_count, movie_list):
        pass

class GFP_95(MicroDevice):
    """
    Sorts worms based on 95th percentile pixel intensity in the background-subtracted fluorescent image. Appropriate for lin-4p::GFP,
    mir-240-786p::GFP, and red autofluorescence.
    """

    def setup_csv(self, file_location, info):
        #TODO: add something to header line to let later reader know what kind of sort it was?
        #What about including a note line for max size and min size? exposure time? etc?
        self.summary_csv_location = file_location.joinpath('summary_csv' + info + '.csv')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'red_autofluor', 'fluorescence_95', 'time', 'direction', 'note']
        self.summary_csv.write(','.join(header) + '\n')

    def analyze(self, cyan_subtracted, tritc_subtracted, worm_mask, worm_count, calibration = False):
        """Class-specific method to determine measurement being sorted.
        Takes background image, bf image of worm, worm count, and returns measurement + sort direction
        Saves fluorescent (GFP) image
        """

        #Relevant param is fluor GFP (cyan)
        worm_fluor = self.find_95th_fluor_amount(cyan_subtracted, worm_mask)
        print('GFP value (95th percentile) = ' + str(worm_fluor))

        if calibration:
            direction = 'straight'
            return worm_fluor, direction

        elif not calibration:

            if worm_fluor >= self.upper_threshold:
                direction = 'up'
            elif worm_fluor <= self.lower_threshold:
                direction = 'down'
            else:
                direction = 'straight'

            return worm_fluor, direction

    def generate_data(self, worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note):
        """Function for returning the right csv lin config for a given class.
        Here, returns the number, size, fluorescence, time between worms, directions, and note (hist or sort)
        sort_param = GFP fluorescence (percentile 95)
        """

        worm_data = [worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note]
        return worm_data
    
class GFP_over_99_mean(MicroDevice):
    """Sorts by fluorescence value given by taking the mean of pixel intensities over the 99th percetnile. Appropriate for mir-47p::GFP.
    """

    def setup_csv(self, file_location, info):
        self.summary_csv_location = file_location.joinpath('summary_csv' + info + '.csv')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'red_autofluor', 'fluorescence_over_99_mean', 'time', 'direction', 'note']
        self.summary_csv.write(','.join(header) + '\n')

    def analyze(self, cyan_subtracted, tritc_subtracted, worm_mask, worm_count, calibration = False):
        """Class-specific method to determine measurement being sorted.
        Takes background image, bf image of worm, worm count, and returns measurement + sort direction
        Saves fluorescent (GFP) image
        """

        #Relevant param is fluor GFP (cyan)
        worm_fluor = self.fluor_over_99_mean(cyan_subtracted, worm_mask)
        print('GFP value (over 99th percentile mean) = ' + str(worm_fluor))

        if calibration:
            direction = 'straight'
            return worm_fluor, direction

        elif not calibration:

            if worm_fluor >= self.upper_threshold:
                direction = 'up'
            elif worm_fluor <= self.lower_threshold:
                direction = 'down'
            else:
                direction = 'straight'

            return worm_fluor, direction

    def generate_data(self, worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note):
        """Function for returning the right csv lin config for a given class.
        Here, returns the number, size, fluorescence, time between worms, directions, and note (hist or sort)
        sort_param = GFP fluorescence over 99 mean
        """

        worm_data = [worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note]
        return worm_data


class Autofluorescence(MicroDevice):

    def setup_csv(self, file_location, info):
        self.summary_csv_location = file_location.joinpath('summary_csv' + info + '.csv')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'autofluorescence', 'time', 'direction', 'note']
        self.summary_csv.write(','.join(header) + '\n')

    def analyze(self, cyan_subtracted, tritc_subtracted, worm_mask, worm_count, calibration = False):
        """
        """

        worm_fluor = self.find_95th_fluor_amount(tritc_subtracted, worm_mask)
        print('Autofluorescence value = ' + str(worm_fluor))

        if calibration:
            direction = 'straight'
            return worm_fluor, direction

        elif not calibration:

            if worm_fluor >= self.upper_threshold:
                direction = 'up'
            elif worm_fluor <= self.lower_threshold:
                direction = 'down'
            else:
                direction = 'straight'

        return worm_fluor, direction

    def generate_data(self, worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note):
        worm_data = [worm_count, worm_size, sort_param, time_between_worms, direction, note]
        return worm_data
    
class Filter(MicroDevice):
    """For generating a population of worms free of small, low af, or otherwise undesirable worms, and also simulating a sort in the process.
    Good worms go straight, other worms go up.
    """
    
    def setup_csv(self, file_location, info):
        self.summary_csv_location = file_location.joinpath('summary_csv' + info + '.csv')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'red_autofluor', 'fluorescence', 'time', 'direction', 'note']
        self.summary_csv.write(','.join(header) + '\n')
        
    def build_hist(self, size):
        pass
    
    def update_hist(self, sort_param):
        pass
    
        
    def analyze(self, cyan_subtracted, tritc_subtracted, worm_mask, worm_count, calibration = False):
        """Sends worm straight, and computes GFP value with 95th percentile brightest pixel, may be useful just for interest of a population.
        """

        worm_fluor = self.find_95th_fluor_amount(cyan_subtracted, worm_mask)
        print('GFP value (95th percentile) = ' + str(worm_fluor))
        direction = 'straight'

        return worm_fluor, direction
        
    def check_note(self, calibration, note, direction):
        """New version of this function specific to the Filter class, sends bad worms up and other worms straight.
        """
        
        if note != 'sort':
            direction = 'up'
        elif note == 'sort':
            #if worm is deemed good and sort is running, note is sort
            pass
        
        return note, direction
    
    def generate_data(self, worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note):
        """Function for returning the right csv lin config for a given class.
        Here, returns the number, size, fluorescence, time between worms, directions, and note (hist or sort)
        sort_param = GFP fluorescence
        """

        worm_data = [worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note]
        return worm_data
   """ 
    def run(self):

        self.setup_csv(self.file_location, self.info)

        self.scope.camera.start_image_sequence_acquisition(frame_count=None, trigger_mode='Software')

        self.boiler = boiler()

        self.worm_count = 0
        self.time_start = time.time()

        #0 Setting Background
        self.set_background_areas(self.worm_count)
        print('setting backgrounds')
        time.sleep(1)
        #input for build hist?

        self.size_threshold = 7600   #hard coding sizes for now
        self.min_worm_size = 2300

        self.build_hist(size=100)

        #self.scope.camera.start_image_sequence_acquisition(frame_count=None, trigger_mode='Software')

        self.sort(calibration = False)

        self.summary_csv.close
    """
    
class Simulate(MicroDevice):
    """Sends worms up, down, and straight repeatedly. Useful for simulating a sort in each direction. Unfortunately, undesirable worms (short, low af, etc) will still be sent straight, so those will need to be manually separated at the end before lifespan assay is carried out).
    """
    def setup_csv(self, file_location, info):
        self.summary_csv_location = file_location.joinpath('summary_csv' + info + '.csv')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'red_autofluor', 'fluorescence', 'time', 'direction', 'note']
        self.summary_csv.write(','.join(header) + '\n')
        
    def build_hist(self, size):
        pass
    
    def update_hist(self, size):
        pass
    
    def analyze(self, cyan_subtracted, tritc_subtracted, worm_mask, worm_count, calibration = False):
        
        worm_fluor = self.find_95th_fluor_amount(cyan_subtracted, worm_mask)
        
        counter = worm_count % 3
        if counter == 0:
            direction = 'straight'
        elif counter == 1:
            direction = 'up'
        elif counter == 2:
            direction = 'down'
            
        return worm_fluor, direction
    
    def generate_data(self, worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note):
        """Function for returning the right csv lin config for a given class.
        Here, returns the number, size, fluorescence, time between worms, directions, and note (hist or sort)
        sort_param = GFP fluorescence
        """

        worm_data = [worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note]
        return worm_data
        
    def run(self):

        self.setup_csv(self.file_location, self.info)

        self.scope.camera.start_image_sequence_acquisition(frame_count=None, trigger_mode='Software')

        self.boiler = boiler()

        self.worm_count = 0
        self.time_start = time.time()

        #0 Setting Background
        self.set_background_areas(self.worm_count)
        print('setting backgrounds')
        time.sleep(1)
        #input for build hist?

        self.size_threshold = 7600   #hard coding sizes for now
        self.min_worm_size = 2600

        self.build_hist(size=100)

        #self.scope.camera.start_image_sequence_acquisition(frame_count=None, trigger_mode='Software')

        self.sort(calibration = False)

        self.summary_csv.close()

class Background(MicroDevice):
    """For gathering data on nonfluorescent worms, background autofluorescence of the system
    """

    #Really just need to take fluor images of non fluor worms, then write code to compare different areas of each image
    #So things like 95th percentile & whole body mean for cyan and green_yellow in age-matched non-fluor worms
    #And for fluor measurements also looking at autofluorescence in the background of the device itself

    def setup_csv(self, file_location, info):
        self.summary_csv_location = file_location.joinpath('summary_csv' + info + '.csv')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'cyan_95th', 'cyan_mean', 'tritc_95th', 'tritc_mean', 'time', 'direction', 'note']
        self.summary_csv.write(','.join(header) + '\n')

    def analyze(self, cyan_subtracted, tritc_subtracted, worm_mask, worm_count, calibration = True):
        """Saves cyan and tritc images, calculates 95th percentile and mean for
        both. Automatically sorts worms straight after imaging.

        Returns sort param as a list.
        """

        cyan_95th = self.find_95th_fluor_amount(cyan_subtracted, worm_mask)
        cyan_mean = self.find_mean_fluor_amount(cyan_subtracted, worm_mask)

        tritc_95th = self.find_95th_fluor_amount(tritc_subtracted, worm_mask)
        tritc_mean = self.find_mean_fluor_amount(tritc_subtracted, worm_mask)

        direction = 'straight'

        sort_param = [cyan_95th, cyan_mean, tritc_95th, tritc_mean]

        print(str(sort_param), direction)

        return sort_param, direction

    def generate_data(worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note):
        """sort_param is a list here, necessitating specific function. Direction and note aren't really salient here.
        """
        worm_data = [worm_count, worm_size] + sort_param + [time_between_worms, direction, note]
        return worm_data
    
    def build_hist(self, size):
        pass
    

class Length(MicroDevice):

    def setup_csv(self, file_location, info):
        self.summary_csv_location = file_location.joinpath('summary_csv' + info + '.csv')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'red_autofluor(95th%)', 'length', 'time', 'direction', 'note']
        self.summary_csv.write(','.join(header) + '\n')

    def check_aspect_ratio(self, length, width):
        """Takes length and width measurements from find_dimensions and returns whether or not the worm is likely
        to have been folded over, meaning it's length measurement is not accurate. Straight worms typically have
        an "aspect ratio" (lengt/width) of >=6.
        """
        aspect_ratio = length/width

        if aspect_ratio >= 6:
            return 'straight'
        elif aspect_ratio < 6:
            print('Worm folded over')
            return 'bent'

    def analyze(self, cyan_subtracted, tritc_subtracted, worm_mask, worm_count, calibration = False):

        length, width = self.find_dimensions(worm_mask)

        print('Worm length = ' + str(length))

        if calibration:
            direction = 'straight'
            return length, direction

        elif not calibration:
            if length >= self.upper_threshold:
                direction = 'up'
            elif length <= self.lower_threshold:
                direction = 'down'
            else:
                direction = 'straight'

        return length, direction

    def generate_data(worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note):
        #Realizing this could be generalized as sort param is always a single var or a list
        #Just write as if type(sort_param) =! list
        worm_data = [worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note]
        return worm_data

class Isaac(MicroDevice):

    def setup_csv(self, file_location, info):
        #TODO: add something to header line to let later reader know what kind of sort it was?
        #What about including a note line for max size and min size? exposure time? etc?
        self.summary_csv_location = file_location.joinpath('summary_csv' + info + '.csv')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'red_autofluor', 'aggregate_count', 'time', 'direction', 'note']
        self.summary_csv.write(','.join(header) + '\n')

    def watershed_aggs(self, image):
        #Image is a GFP image as an ndarray
        high_thresh = numpy.percentile(image,99.99) #Alternatively, can use the 97th percentile mask pixel if the mask is available
        low_thresh = numpy.percentile(image,99)
        markers = numpy.zeros_like(image)
        markers[image > high_thresh] = 2
        markers[image < low_thresh] = 1
        watershed_im = watershed(sobel(image),markers)
        aggregates = numpy.max(label(watershed_im))
        return aggregates

    def analyze(self, cyan_subtracted, tritc_subtracted, worm_mask, worm_count, calibration = False):
        """Class-specific method to determine measurement being sorted.
        Takes background image, bf image of worm, worm count, and returns measurement + sort direction
        Saves fluorescent (GFP) image
        """

        #Relevant param is fluor GFP (cyan)
        aggregate_count = self.watershed_aggs(cyan_subtracted)
        print('Number of aggregates = ' + str(aggregate_count))

        if calibration:
            direction = 'straight'
            return aggregate_count, direction

        elif not calibration:

            if aggregate_count >= self.upper_threshold:
                direction = 'up'
            elif aggregate_count <= self.lower_threshold:
                direction = 'down'
            else:
                direction = 'straight'

            return aggregate_count, direction

    def generate_data(self, worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note):
        """Function for returning the right csv lin config for a given class.
        Here, returns the number, size, fluorescence, time between worms, directions, and note (hist or sort)
        sort_param = GFP fluorescence
        """

        worm_data = [worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note]
        return worm_data

    def go(self):

        self.setup_csv(self.file_location, self.info)

        self.scope.camera.start_image_sequence_acquisition(frame_count=None, trigger_mode='Software')

        self.boiler = boiler()

        worm_count = 0
        self.set_background_areas(worm_count)
        print('setting backgrounds')
        time.sleep(1)
        #input for build hist?

        self.size_threshold = 7500   #hard coding sizes for now
        self.min_worm_size = 2500

        self.build_hist()

        self.scope.camera.start_image_sequence_acquisition(frame_count=None, trigger_mode='Software')

        self.sort(calibration = False)

        self.summary_csv.close()

class Movie(MicroDevice):

    def setup_csv(self, file_location, info):
        self.summary_csv_location = file_location.joinpath('summary_csv' + info + '.csv')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'autofluorescence', 'GFP_fluorescence', 'time', 'direction', 'note']
        self.summary_csv.write(','.join(header) + '\n')

    def analyze(self, cyan_subtracted, tritc_subtracted, worm_mask, worm_count, calibration = False):
        """
        """

        worm_fluor = self.find_95th_fluor_amount(cyan_subtracted, worm_mask)

        if worm_count <= 5:
            direction = 'straight'

        else:
            counter = worm_count % 3
            if counter == 0:
                direction = 'straight'
            elif counter == 1:
                direction = 'up'
            elif counter == 2:
                direction = 'down'

        return worm_fluor, direction

    def generate_data(self, worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note):
        worm_data = [worm_count, worm_size, autofluorescence, sort_param, time_between_worms, direction, note]
        return worm_data

    def update_hist(self, sort_param):
        pass

    def capture_image(self, type_of_image):
        #overwriting the normal capture_image function for the movie maker class such that it scales and saves images.

        type_of_image()
        self.scope.camera.send_software_trigger()
        image = self.scope.camera.next_image()

        scaled_image = self.scale_image(image, type_of_image)
        self.film_reel.append(scaled_image.astype('uint8'))

        return image

    def scale_image(self, image, type_of_image):
        if type_of_image == self.bright:
            scaled_image = zplib.image.colorize.scale(image, min=0, max=53000)
        elif type_of_image == self.cyan:
            scaled_image = zplib.image.colorize.scale(image, min=0, max=20000)
        elif type_of_image == self.green_yellow:
            scaled_image = zplib.image.colorize.scale(image, min=0, max=2000)
        elif type_of_image == 'mask':
            scaled_image = zplib.image.colorize.scale(image, min=0, max=1)
            self.film_reel.append(scaled_image.astype('uint8'))

        return scaled_image


    def make_movie(self, worm_count, movie_list):
        zplib.image.write_movie.write_movie(movie_list, 'worm' + str(worm_count) + '.mp4')

    def go(self):

        self.setup_csv(self.file_location, self.info)

        self.scope.camera.start_image_sequence_acquisition(frame_count=None, trigger_mode='Software')

        self.boiler = boiler()

        self.film_reel = collections.deque(maxlen=200)

        worm_count = 0
        self.set_background_areas(worm_count)
        print('setting backgrounds')
        time.sleep(1)

        self.size_threshold = 8000   #hard coding sizes for now
        self.min_worm_size = 3000

        self.sort(calibration = False)

        self.summary_csv.close()
