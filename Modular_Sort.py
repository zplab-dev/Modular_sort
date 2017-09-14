#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 12:32:55 2017

@author: Tim

Code to sort Wild Type, Red fluorescent, and Green fluorescent worms

"""

import iotool
from scope import scope_client
from pathlib import Path
import scipy
import numpy
import time
import backgroundSubtraction
import freeimage
import threading
import csv



#Setting useful constants
#Areas for image analysis
IMAGE_SIZE = (1280, 1080)
BOILER_AREA = (slice(100,750), slice(530,590))
DETECTION_AREA = (slice(100,950),slice(535, 585))
POSITION_AREA = (slice(100,640), slice(535, 585))
CLEARING_AREA = POSITION_AREA #Modified clearing area from 620 to 520 to insure the worm has left the area before being determined to be sorted.
FLUORESCENT_AREA = POSITION_AREA
QUEUE_AREA = (slice(1030,1230),slice(520,590))

FLUOR_PIXEL_BRIGHT_VALUE = 500
DETECTION_THRES = 4
PUSH_THRESH = 2.5
QUEUE_THREH = .6
POSITION_THRES = 3
CLEARING_THRES = 3
LOST_CUTOFF = 1.5
MOVEMENT_THRES = 1.3
NULL_THRESH = 5
DOUBLE_THRES = 1.3

BUBBLE_THRES = 10
CLOG_THRESH = 10

FLUORESCENCE_PERCENTILE = 99
BACKGROUND_FRACTION = .99 #For Setting worm mask

CYAN_EXPOSURE_TIME = 4
YELLOW_EXPOSURE_TIME = 8
BRIGHT_FIELD_EXPOSURE_TIME = 2
MAX_PUSH_TIME = .2

LIGHT_DELAY = .05
PICTURE_DELAY = .01
SORTING_INTERVAL = .05
MAX_SORTING_TIME = 1.5

BACKGROUND_REFRESH_RATE = 100000
PROGRESS_RATE = 100

#Setting useful commands for device control

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
RELIEF_CHANNEL_SUCK = 'sl D2'
RELIEF_CHANNEL_PRESSURE = 'sh D2'

#Setting file locations for saving images

#File_location = '/mnt/iscopearray/Nonet_Tim/Test##_##_##'


def boiler():
    """
    Returns a boolean array of possible locations of a worm duing sorting
    """
    boiler = numpy.ones(IMAGE_SIZE, dtype=bool)
    boiler[BOILER_AREA] = False
    return boiler

class MicroDevice(threading.Thread):
    """
    Class for running a Microfluidic Device 
    """ 
    
    def __init__(self, exp_direct):
        """
        Initalizes the scope and device 
        """
        self.scope, scope_properties = scope_client.client_main()
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
        self.scope.camera.readout_rate = '280 MHz'
        self.scope.camera.binning = '2x2'
        self.scope.nosepiece.magnification = 5
        self.scope.tl.lamp.enabled = True
        
        self.device = iotool.IOTool("/dev/ttyMicrofluidics")
        
        self.file_location = Path(exp_direct)
        self.file_location.mkdir(mode=0o777, parents=True, exist_ok=True)
        self.summary_location = self.file_location.joinpath('summary.txt')
        self.summary_statistics = open(str(self.summary_location),'w')
        self.data_location = self.file_location.joinpath('wormdata.csv')
        
        self.device_stop_load()
        
        self.worm_count = 0
        self.up = 0
        self.down = 0
        self.straight = 0
        #worm_data = list() 
        
        #Pausing stuff
        self.running = False
        super().__init__(daemon=True)
        self.quitting = False
        self.cleared = False
        
    def write_csv_file(self, worm_data):
         with open(str(self.data_location), 'w', newline = '') as wormdata:
             wormwriter = csv.writer(wormdata, dialect = 'excel')
             wormwriter.writerow(['Worm Number', 'Worm Size', 'Worm Direction',
                                  'Detection Time', 'Analysis Time', 'Sort Time',
                                  'fluorMcherry', 'fluorGFP'])
             for worms in worm_data: 
                 wormwriter.writerow(worms)

        
    def set_scope(self):
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
        self.scope.camera.readout_rate = '280 MHz'
        self.scope.camera.binning = '2x2'
        self.scope.tl.lamp.enabled = True
        
    def pause(self):
        print('Paused')
        self.running = False
    
    def resume(self):
        print('unpausing')
        self.running = True
        self.pause_tell = False
        
    def quit(self):
        print('Quiting')
        self.quitting = True
        
    def clear(self):
        print('Cleared')
        self.cleared = True
            
    def device_stop_run(self):
        """
        Command that stops the device from sorting or loading worm.
        Device is set to a safe steady state
        """
        self.scope.camera.end_image_sequence_acquisition()
        self.device.execute(PUSH_CHANNEL_PRESSURE, SEWER_CHANNEL_PRESSURE,
                            UP_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_PRESSURE,
                            DOWN_CHANNEL_PRESSURE,RELIEF_CHANNEL_PRESSURE)
        
    def device_start_load(self):
        """
        Command that toggles the device to being loading worms into the device and 
        will continue to push worms into the device until given another command
        """
        print('starting loading')
        self.device.execute(PUSH_CHANNEL_STATIC, SEWER_CHANNEL_SUCK,
                            UP_CHANNEL_PRESSURE,STRAIGHT_CHANNEL_PRESSURE,
                            DOWN_CHANNEL_PRESSURE, RELIEF_CHANNEL_SUCK)
        
    def device_push_queue(self):
        print('pushing queue')
        self.time_queue_push_start = time.time()
        self.device.execute(RELIEF_CHANNEL_PRESSURE)
        
    def device_position_worm(self):
        print(' ')
        print(' ')
        print('Worm has been sucesfully pushed into the device')
        self.time_seen = time.time()
        self.device.execute(PUSH_CHANNEL_PRESSURE, RELIEF_CHANNEL_SUCK)
        
    def device_stop_load(self):
        """
        Command that toggles the device to stop loading new worms
        """
        self.device.execute(SEWER_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_PRESSURE, 
                            UP_CHANNEL_PRESSURE, DOWN_CHANNEL_PRESSURE,
                            PUSH_CHANNEL_PRESSURE, RELIEF_CHANNEL_PRESSURE)
        
    def device_sort(self, direction):
        """
        Command that toggles the device to sort a worm in given direction
        Have yet to test if having teh pusher channel be set to static when sorting to help worms load faster.
        
        """
        print('Sorting worm'+ direction)
        if direction == 'up':
            self.device.execute(SEWER_CHANNEL_PRESSURE,
                                UP_CHANNEL_SUCK,
                                RELIEF_CHANNEL_SUCK)
        elif direction == 'down':
            self.device.execute(SEWER_CHANNEL_PRESSURE,
                                DOWN_CHANNEL_SUCK,
                                RELIEF_CHANNEL_SUCK)
        else:
            self.device.execute(SEWER_CHANNEL_PRESSURE,
                                STRAIGHT_CHANNEL_SUCK,
                                RELIEF_CHANNEL_SUCK) 
        cleared_image = self.capture_image(self.bright)
        while not self.check_cleared(cleared_image):
            if self.cleared:
                break
            self.flutter_direction(direction)
            cleared_image = self.capture_image(self.bright)
        time.sleep(SORTING_INTERVAL)

    def flutter_direction(self, direction):
        print('fluttering')
        if direction == 'up':
            self.device.execute(UP_CHANNEL_PRESSURE)
            time.sleep(.05)
            self.device.execute(UP_CHANNEL_SUCK)
        elif direction == 'down':
            self.device.execute(DOWN_CHANNEL_PRESSURE)
            time.sleep(.05)
            self.device.execute(DOWN_CHANNEL_SUCK)  
        else:
            self.device.execute(STRAIGHT_CHANNEL_PRESSURE)
            time.sleep(.05)
            self.device.execute(STRAIGHT_CHANNEL_SUCK) 

            
    def device_clear_bubbles(self):
        """
        Command that toggles the device to alternate between blowing and sucking
        in an attempt to clear bubbles from the device
        """
        try:
            while True:
                self.device.execute(SEWER_CHANNEL_PRESSURE,
                                    STRAIGHT_CHANNEL_SUCK, 
                                    UP_CHANNEL_SUCK,
                                    DOWN_CHANNEL_SUCK,
                                    PUSH_CHANNEL_PRESSURE,
                                    RELIEF_CHANNEL_PRESSURE)
                time.sleep(1)
                self.device.execute(SEWER_CHANNEL_SUCK,
                                    PUSH_CHANNEL_STATIC,
                                    STRAIGHT_CHANNEL_PRESSURE,
                                    UP_CHANNEL_PRESSURE,
                                    DOWN_CHANNEL_PRESSURE,
                                    RELIEF_CHANNEL_SUCK)
        except KeyboardInterrupt:
            self.device.execute(PUSH_CHANNEL_PRESSURE,
                                SEWER_CHANNEL_PRESSURE,
                                UP_CHANNEL_PRESSURE, 
                                STRAIGHT_CHANNEL_PRESSURE,
                                DOWN_CHANNEL_PRESSURE,
                                RELIEF_CHANNEL_PRESSURE)
        
        
    def device_clear_tubes(self):
        """
        Command that toggles the deivce to set all tubes to push water and hopefully
        clear the tubes of debris
        """
        self.device.execute(SEWER_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_PRESSURE, 
                            UP_CHANNEL_PRESSURE, DOWN_CHANNEL_PRESSURE,
                            PUSH_CHANNEL_PRESSURE, RELIEF_CHANNEL_PRESSURE)       
    def lamp_off(self):
        """
        Turns the lamp off
        """
        self.scope.il.spectra.cyan.enabled = False
        self.scope.tl.lamp.enabled = False
        self.scope.il.spectra.green_yellow.enabled = False
        
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
        self.scope.camera.exposure_time = YELLOW_EXPOSURE_TIME
        time.sleep(PICTURE_DELAY)
  
    def capture_image(self, type_of_image):
        """
        Returns an int 32 image of with the features passed by the set up 
        function type_of_image
        """
        type_of_image()
        self.scope.camera.send_software_trigger()
        return self.scope.camera.next_image().astype('int32')
        
    def save_image(self, image, name):
        save_location = str(self.file_location) + '/' + name + '.png'
        freeimage.write(image.astype('uint16'), save_location,
                        flags=freeimage.IO_FLAGS.PNG_Z_BEST_SPEED)
                        
    def set_background_areas(self):
        """
        Function that sets the background values for areas of interest in an
        list of numpy arrays
        """
        raise NotImplementedError('No sorting method given') 

    def device_clear_lost_worm(self):
        print('Worm has been lost')
        self.device_sort('straight lost')
        self.worm_direction = 'straight'
        self.summary_statistics.write("\nWorm " + str(self.worm_count) + "was lost")
        
    def check_position(self, current_image, detected_image):
        """
        Function that determines if a worm has been positioned
        The worm is positioned because the change is small.
        """
        print('checking position')
        worm_movment = abs(current_image[POSITION_AREA]
         - detected_image[POSITION_AREA])
        return  ((numpy.sum(worm_movment) - self.positioned_background) 
            < POSITION_THRES * self.positioned_background)

    def check_queue(self, current_image):
        """
        print('Detected Value = ' + str((numpy.sum(numpy.abs(current_image[QUEUE_AREA].astype('int32')
                                     - background[QUEUE_AREA].astype('int32'))))))
        print('Required Value = ' + str(QUEUE_THREH  * self.detect_background))
        """
        #print('Checking Queue')
        return ((numpy.sum(numpy.abs(current_image[QUEUE_AREA].astype('int32')
                                     - self.background[QUEUE_AREA].astype('int32'))))
                > QUEUE_THREH  * self.detect_background)

    def check_lost(self, current_image):
        """
        Function that determines if a worm was lost
        Worm is deciced lost because image is close enough to background.
        """
        print('Checking lost')
        worm_visibility = abs(current_image[BOILER_AREA]
                              - self.background[BOILER_AREA])
        return ((numpy.sum(worm_visibility) - self.positioned_background) 
            < LOST_CUTOFF * self.positioned_background) 

    def check_cleared(self, current_image):
        print('Checking Clear')
        worm_visibility = abs(current_image[DETECTION_AREA]
                              - self.background[DETECTION_AREA])
        return ((numpy.sum(worm_visibility) - self.positioned_background) 
            < LOST_CUTOFF * self.positioned_background) 
    
    def check_worm(self, current_image):
        pass

    def check_pushed_forwards(self, current_image):
        """
        Function that returns if a worm has been detected
        """
        #print('Checking pushing forwards')
        #rint('Detected Value:' + str((numpy.sum(numpy.abs(current_image[DETECTION_AREA] - self.background[DETECTION_AREA]))- self.detect_background)))
        #print('Required Value:' + str( PUSH_THRESH * self.detect_background))
        if time.time() - self.time_queue_push_start > MAX_PUSH_TIME:
            return True
        return ((numpy.sum(numpy.abs(current_image[DETECTION_AREA] 
                                     - self.background[DETECTION_AREA]))-self.detect_background) 
        > PUSH_THRESH  * self.detect_background)
    
    def worm_mask(self, worm_iamge):
        """
        Function that saves an image of the mask of a worm and returns the size of mask (in number of pixels)
        Depending on how background.backgroundSubtraction.clean_dust_and_holes(image) works this function might neeed to be modified so that it returns the appropriate image.
        Currently the clean_dust_and_holes does not return the clean image and actually modifies the passed image.
        """
        subtracted_image = worm_image - self.bacgrkound
        floored_image = backgroundSubtraction.percentile_floor(subtracted_image, .99)
        floored_image[self.shade_area] = 0
        backgroundSubtraction.clean_dust_and_holes(floored_image)
        return floored_image.astype('bool')
        
    def mask_size(self, worm_mask):
        return numpy.count_nonzero(worm_mask)

    def analyze_worm(self):
        """
        function that tells the device what sorting/analzying method to use:
        Is overwritten by a super class
        """
        raise NotImplementedError('No sorting method given')

    def device_clear_and_reset(self):
        self.device_sort('straight')
        self.background = self.capture_image(self.bright)
        self.save_image(self.background, 'background' + str(self.worm_count))
        self.set_background_areas()
        self.cleared = False
        print('Reset background')

    def check_size_worm(self, current_image):
        worm_mask = self.worm_mask(current_image)
        worm_size = self.mask_size(worm_mask)
        print('Size of worm before sorting: ' + str(worm_size))
        if worm_size > self.size_threshold or worm_size < self.min_size_threshold:
            print('Detected Bad Worm')
            self.clear_worms()
        else:
            print('Properly sized mask')

    def cycle_background_reset(self):
        self.device_stop_load()
        self.device_sort('straight')
        self.background = self.capture_image(self.bright)
        self.save_image(self.background, 'background' + str(self.worm_count))
        self.set_background_areas()
        self.device_start_load()
      
    def initialize_sorting(self):
        """
        This function is run at the start of sorting 
        """
        self.set_background_areas()
        self.boiler = boiler()
        print('Backgrounds have been set.')

        
        #1 Loading Worms
        self.device_start_load()
        self.time_start = time.time()
        
        #2 Detect Worms
        self.time_between_worms = list()
        self.time_to_position_worms = list() 
        self.worm_count = 0

        self.background = self.capture_image(self.bright)
        self.save_image(self.background, 'background')   
                                    
    def run(self):
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
        self.scope.camera.start_image_sequence_acquisition(
            frame_count=None, trigger_mode='Software')
        cycle_count = 0
        self.initialize_sorting()
        #0 Setting Background
        try:
            print('entering loop')
            while not self.quitting:
                self.pause_tell = False
                while not self.running:
                    if self.pause_tell == False:
                        print('Paused')
                        self.pause_tell = True
                    paused_image = self.capture_image(self.bright)
                cycle_count += 1
                current_image = self.capture_image(self.bright)
                if self.check_queue(current_image):
                    self.device_push_queue()
                    while not self.quitting:
                        current_image = self.capture_image(self.bright)
                        if self.check_pushed_forwards(current_image):
                            #3 stop worms
                            self.device_position_worm()
                            while not self.quitting:
                                detected_image = current_image
                                current_image = self.capture_image(self.bright)
                                if self.check_lost(current_image):
                                    self.device_clear_lost_worm()
                                    break
                                elif self.check_position(current_image, detected_image):
                                    self.check_worm(current_image)
                                    self.analyze_worm(current_image)
                                    break
                                if self.cleared:
                                    self.device_clear_and_reset()
                                    break
                                else:
                                    detected_image = current_image
                            print('Breaking out of second loop')
                            self.device_start_load()
                            break
                    
                    
                if cycle_count % PROGRESS_RATE == 0:
                    print(str(PROGRESS_RATE) + ' Cycles')
                    
                elif cycle_count % BACKGROUND_REFRESH_RATE == 0:
                    print(str(cycle_count) + ' Cycles Reseting Background')
                    
        except KeyboardInterrupt:
            pass
        finally:
            self.summary_statistics.write('\n Average worm detection time :' 
                                          + str(numpy.mean(self.time_between_worms)) 
                                          + '\n Average worm positioning time :' 
                                          + str(numpy.mean(self.time_to_position_worms)))
            self.device_stop_run()
            print('fianlly went')
            self.summary_statistics.close()
                
    def main(self):
        self.device = iotool.IOTool("/dev/ttyMicrofluidics")
        return  

class NoSort(MicroDevice):
    """
    Sub class of MicroDevice that simply images the worms and sends them straight
    """
    
    def set_background_areas(self):
        """
        Function that sets the background values for areas of interest in an
        list of numpy arrays
        """
        image1 = self.capture_image(self.bright)
        image2 = self.capture_image(self.bright)
        subtracted = abs(image1 - image2)
        self.detect_background = numpy.sum(subtracted[DETECTION_AREA])
        self.positioned_background = numpy.sum(subtracted[POSITION_AREA])
        self.clear_background = numpy.sum(subtracted[CLEARING_AREA])

    def analyze_worm(self, current_image):
        """
        function that tells the device what sorting/analzying method to use:
        Is overwritten by a super class
        """
        self.device_sort('straight')
        self.worm_direction = 'straight'

class Alternate(MicroDevice):
    """
    """
    def set_background_areas(self):
        """
        Function that sets the background values for areas of interest in an
        list of numpy arrays
        """
        image1 = self.capture_image(self.bright)
        time.sleep(PICTURE_DELAY)
        image2 = self.capture_image(self.bright)
        subtracted = abs(image1.astype('int32') - image2.astype('int32'))
        self.detect_background = numpy.sum(subtracted[DETECTION_AREA])
        self.positioned_background = numpy.sum(subtracted[POSITION_AREA])
        self.clear_background = numpy.sum(subtracted[CLEARING_AREA])
        
    def analyze_worm(self, current_image):
        """
        function that tells the device what sorting/analzying method to use:
        Is overwritten by a super class
        """
        start = time.time()
        while time.time() - start < 1:
            current_image = self.capture_image(self.bright)
        direction = self.worm_count % 3
        dirlst = ['straight','up','down']
        self.device_sort(dirlst[direction])
        self.worm_direction = dirlst[direction]

class Mir71(MicroDevice):
    """
    """
    def set_background_areas(self):
        """
        Function that sets the background values for areas of interest in an
        list of numpy arrays
        """
        image1 = self.capture_image(self.bright)
        time.sleep(PICTURE_DELAY)
        image2 = self.capture_image(self.bright)
        subtracted = abs(image1.astype('int32') - image2.astype('int32'))
        self.detect_background = numpy.sum(subtracted[DETECTION_AREA])
        self.positioned_background = numpy.sum(subtracted[POSITION_AREA])
        self.clear_background = numpy.sum(subtracted[CLEARING_AREA])
        
        self.cyan_background = self.capture_image(self.cyan)
        self.save_image(self.cyan_background, 'cyan_background')

        self.lamp_off()
        self.scope.tl.lamp.enabled = True
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME

    
    def auto_set_up(self):
        top_threshold = input('What do you want the top X% of expressiont to be = ')
        self.upper_mir71_threshold = int(top_threshold)
        bottom_threshold = input('What do you want the bottom X% of expressiont to be = ')
        self.bottom_mir71_threshold = int(bottom_threshold)
        size = input('What do you want the Size Threshold = ')
        self.size_threshold = int(size)* DOUBLE_THRES
        min_size = int(input('What is the min size?'))
        self.min_size_threshold = min_size
            
    def update_thresholds(self, gfp_amount):
        self.run_fluorescence.append(gfp_amount)
        self.bottom_mir71_threshold = numpy.percentile(self.run_fluorescence, 10)
        self.upper_mir71_threshold = numpy.percentile(self.run_fluorescence, 90)
        
    def find_fluor_amount(self, subtracted_image, worm_mask):
        gfp_image = subtracted_image[worm_mask]
        gfp_count = numpy.percentile(gfp_image, 95)
        return gfp_count

class Mir71_Sort(Mir71):
    """
    """
    def __init__(self, exp_direct, min_size, max_size, 
        bottom_mir71_threshold, upper_mir71_threshold):
        super().__init__(exp_direct)
        self.max_size_threshold = max_size
        self.min_size_threshold = min_size
        self.upper_mir71_threshold = upper_mir71_threshold
        self.bottom_mir71_threshold = bottom_mir71_threshold
        
    def analyze_worm(self, worm_image):
        gfp_fluor_image = self.capture_image(self.cyan)
        self.save_image(gfp_fluor_image, 'fluor_gfp', True)
        gfp_subtracted = abs(gfp_fluor_image.astype('int32')
                             - self.cyan_background.astype('int32'))
        worm_mask = self.worm_mask(worm_image) 
        worm_fluor = self.find_fluor_amount(gfp_subtracted, worm_mask)
        
        self.summary_statistics.write("Gfp Fluorescence: " 
            + str(worm_fluor))
        
        print('GFP value = ' + str(worm_fluor))
        
        double_image = self.capture_image(self.bright)
        worm_size = self.mask_size(self.worm_mask(double_image))
        print("Size of worm after imaging :" + str(worm_size))
        
        if worm_size > self.max_size_threshold or worm_size < self.min_size_threshold:
            print('Detected Double Worm')
            self.save_image(double_image, 'doubled worm_analyze' + str(self.worm_count))
            self.summary_statistics.write( '\n doubled worm size of: '+ str(worm_size)) 
            self.device_sort('straight')
            self.worm_direction = 'straight'
            self.summary_statistics.write("Straight\n")
            print('Doubled worms sorted Straight')
        elif worm_fluor > self.upper_mir71_threshold:
            self.update_thresholds(worm_fluor)
            self.up_worms.append(worm_fluor)
            self.up += 1
            self.device_sort('up')
            self.worm_direction = 'up'
            self.summary_statistics.write("Up\n")
            print('Worm sorted Up    ' + str(self.up))
        elif worm_fluor < self.bottom_mir71_threshold:
            self.update_thresholds(worm_fluor)
            self.down += 1
            self.device_sort('down')
            self.worm_direction = 'down'
            self.summary_statistics.write("Down\n")
            print('Worm sorted Down   ' + str(self.down))
        else:
            self.update_thresholds(worm_fluor)
            self.straight += 1
            self.device_sort('straight')
            self.worm_direction = 'straight'
            self.summary_statistics.write("straight\n")
            print('Worm sorted straight    ' + str(self.straight))

class Mir71_SetUp(Mir71):
    """
    """
    
    def __init__(self, exp_direct):
        super().__init__(exp_direct)
        self.max_worm_size = int(input('Whats the initial size threshold?'))
        self.min_worm_size = int(input('What is the initial small size threshold?'))
        self.num_of_worms = int(input('How many worms to survey?'))
        self.size = list()
        self.fluorescence = list()

    def find_thresholds(self, num_of_worms):
        """
        Input desired number of worms to build histograms of fluorescence and size.
        """

        self.run()

        avg_size = numpy.mean(self.size)
        size_90 = numpy.percentile(self.size, 90)
        size_10 = numpy.percentile(self.size, 10)
        avg_gfp = numpy.mean(self.fluorescence)
        self.bottom_mir71_threshold = numpy.percentile(self.fluorescence, 10)
        self.upper_mir71_threshold = numpy.percentile(self.fluorescence, 90)
        self.size_threshold = size_90 * DOUBLE_THRES
        self.min_size_threshold = size_10 * .5

        print('Avg Size =' + str(avg_size))
        print('90_size =' + str(size_90))
        print('10_size =' + str(size_10))
        print('Avg Gfp =' + str(avg_gfp))
        print('10_gfp =' + str(self.bottom_mir71_threshold))
        print('90_gfp =' + str(self.upper_mir71_threshold))

    def analyze_worm(self, current_image):
        if self.worm_count > (num_of_worms - 1):
            self.quit()
        worm_mask = self.worm_mask(current_image)
        worm_size = self.mask_size(worm_mask)
        if worm_size > max_worm_size or worm_size < min_worm_size:
            self.worm_direction = 'straight'
        else:
            self.worm_count += 1
            print('Worm number ' + str(worm_count) + ' out of ' + str(num_of_worms))
            self.save_image(current_image, 'calibration_worm'+ str(worm_count))
            print('images_saved')
            current_image = self.capture_image(self.cyan)
            self.save_image(current_image, 'calibration_worm_fluor' + str(worm_count))
            gfp_image = abs(current_image.astype('int32')- self.cyan_background.astype('int32'))
            gfp_amount = self.find_fluor_amount(gfp_image, worm_mask)
            print('GFP amount = ' + str(gfp_amount))
            self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
            self.lamp_off()
            self.scope.tl.lamp.enabled = True
            self.size.append(worm_size)
            self.fluorescence.append(gfp_amount)
            self.device_sort('straight')
            self.worm_direction = 'straight'

class fluorRedGreen(MicroDevice):
    """
    """
    def __init__(self, exp_direct):
        super().__init__(exp_direct)
        gfp = input('What do you want as the GFP Threshold = ')
        self.gfp_threshold = int(gfp)
        mcherry = input('What do you want as the mcherry Threshold = ')
        self.mcherry_threshold = int(mcherry)
        size = input('What do you want the Size Threshold = ')
        self.size_threshold = int(size)
        min_size = input('What do you want the small size threshold = ')
        self.min_size_threshold = int(min_size)
        
    def find_fluor_amount(self, image):
        
        blurred = scipy.ndimage.gaussian_filter(image, sigma = 2)
        low_vales = blurred < FLUOR_PIXEL_BRIGHT_VALUE
        blurred[low_vales] = 0
        scipy.ndimage.morphology.binary_erosion(blurred, None, 4)
        scipy.ndimage.morphology.binary_dilation(blurred, None, 4)
        return numpy.sum(blurred[FLUORESCENT_AREA])
    
    def set_background_areas(self):
        """
        Function that sets the background values for areas of interest in an
        list of numpy arrays
        """
        image1 = self.capture_image(self.bright)
        time.sleep(PICTURE_DELAY)
        image2 = self.capture_image(self.bright)
        subtracted = abs(image1.astype('int32') - image2.astype('int32'))
        self.detect_background = numpy.sum(subtracted[DETECTION_AREA])
        self.positioned_background = numpy.sum(subtracted[POSITION_AREA])
        self.clear_background = numpy.sum(subtracted[CLEARING_AREA])
        
        self.cyan_background = self.capture_image(self.cyan)
        self.save_image(self.cyan_background, 'cyan_background')
        
        self.green_background = self.capture_image(self.green_yellow)
        self.save_image(self.green_background, 'green_background')

        self.lamp_off()
        self.scope.tl.lamp.enabled = True
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
        
    def analyze_worm(self, current_image):
        """
        function that tells the device what sorting/analzying method to use:
        Is overwritten by a super class
        """
        gfp_fluor_image = self.capture_image(self.cyan)
        gfp_subtracted = abs(gfp_fluor_image.astype('int32')
                             - self.cyan_background.astype('int32'))
        color_value_cyan = self.find_fluor_amount(gfp_subtracted)
        self.save_image(gfp_fluor_image, 'fluor_gfp' + str(self.worm_count))
        
        mcherry_fluor_image = self.capture_image(self.green_yellow)
        mcherry_subtracted = abs(mcherry_fluor_image.astype('int32')
                                 - self.green_background.astype('int32'))
        color_value_green = self.find_fluor_amount(mcherry_subtracted)
        self.save_image(mcherry_fluor_image, 'fluor_mcherry' + str(self.worm_count))

        print('GFP value = ' + str(color_value_cyan))
        print('mCherry value = ' + str(color_value_green))
        
        double_image = self.capture_image(self.bright)
        worm_size = self.mask_size(self.worm_mask(double_image))
        print("Size of worm after imaging :" + str(worm_size))
        
        self.lamp_off()
        self.scope.tl.lamp.enabled = True
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
        
        self.summary_statistics.write("Gfp Fluorescence: " 
            + str(color_value_cyan) 
            + "\nGfp Required: " 
            + str(self.gfp_threshold) 
            + "\nMcherry Fluorescence: " 
            + str(color_value_green) 
            + "\nMcherry Required :" 
            + str(self.mcherry_threshold) 
            + "\nSize of Worm after imaging: "
            + str(worm_size) 
            + "\n")

        if worm_size > self.size_threshold:
            print('Detected Double Worm')
            self.save_image(double_image, 'doubled worm_analyze' + str(self.worm_count))
            self.summary_statistics.write( '\n doubled worm size of: '+ str(worm_size)) 
            self.device_sort('straight')
            self.worm_direction = 'straight'
            self.summary_statistics.write("Straight\n")
            print('Worm sorted Straight')     

        elif ((color_value_cyan > self.gfp_threshold) 
        and (color_value_green < self.mcherry_threshold)):
            #Worm is determined to be green and not red.
            self.device_sort('up')
            self.worm_direction = 'up'
            self.summary_statistics.write("Up\n")
            print('Worm sorted Up')

        elif ((color_value_cyan < self.gfp_threshold) 
        and (color_value_green > self.mcherry_threshold)):
            #Worm is detremined to be red and not green.
            self.device_sort('down')
            self.worm_direction = 'down'
            self.summary_statistics.write("Down\n")
            print('Worm sorted Down')

        else:
            self.device_sort('straight')
            self.worm_direction = 'straight'
            self.summary_statistics.write("Straight\n")
            print('Worm sorted Straight')
            
            
class FluoRedGreen_Setup(fluorRedGreen):
    """
    """
    def __init__(self, exp_direct):
        super().__init__(exp_direct)
        self.max_worm_size = int(input('Whats the initial size threshold?'))
        self.min_worm_size = int(input('What is the initial small size threshold?'))
        self.num_of_worms = int(input('How many worms to survey?'))
        self.size = list()
        self.GFP_fluorescence = list()
        self.mCherry_fluorescence = list()

    def find_thresholds(self, num_of_worms):
        """
        Input desired number of worms to build histograms of fluorescence and size.
        """

        self.run()

        avg_size = numpy.mean(self.size)
        size_90 = numpy.percentile(self.size, 90)
        size_10 = numpy.percentile(self.size, 10)
        avg_gfp = numpy.mean(self.GFP_fluorescence)
        avg_mCherry = numpy.mean(self.mCherry_fluorescence)
       
        self.size_threshold = size_90 * DOUBLE_THRES
        self.min_size_threshold = size_10 * .5

        print('Avg Size =' + str(avg_size))
        print('90_size =' + str(size_90))
        print('10_size =' + str(size_10))
        print('Avg Gfp =' + str(avg_gfp))
        print('Avg mCherry = ' + str(avg_mCherry))
        print('10_gfp =' + str(self.bottom_mir71_threshold))
        print('90_gfp =' + str(self.upper_mir71_threshold))

    def analyze_worm(self, current_image):
        if self.worm_count > (num_of_worms - 1):
            self.quit()
        worm_mask = self.worm_mask(current_image)
        worm_size = self.mask_size(worm_mask)
        if worm_size > max_worm_size or worm_size < min_worm_size:
            print('Bad Worm')
            self.worm_direction = 'straight'
        else:
            self.worm_count += 1
            print('Worm number ' + str(worm_count) + ' out of ' + str(num_of_worms))
            self.save_image(current_image, 'calibration_wormGFP'+ str(worm_count))
            print('images_saved')
            
            current_image = self.capture_image(self.cyan)
            self.save_image(current_image, 'calibration_worm_GFPfluor' + str(worm_count))
            gfp_image = abs(current_image.astype('int32')- self.cyan_background.astype('int32'))
            gfp_amount = self.find_fluor_amount(gfp_image, worm_mask)
            print('GFP amount = ' + str(gfp_amount))
            
            current_image = self.capture_image(self.green_yellow)
            self.save_image(current_image, 'calibration_wormMcherry' + str(worm_count))
            mCherry_image = abs(current_image.astype('int32')-self.green_background.astype('int32'))
            mcherry_amount = self.find_fluor_amount(mcherry_image, worm_mask)
            print('mCherry amount = ' + str(mcherry_amount))
        
            self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
            self.lamp_off()
            self.scope.tl.lamp.enabled = True
            self.size.append(worm_size)
            self.GFP_fluorescence.append(gfp_amount)
            self.mCherry_fluorescence.append(mcherry_amount)
            self.device_sort('straight')
            self.worm_direction = 'straight'

