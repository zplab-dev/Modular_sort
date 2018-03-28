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
import scipy
import numpy
import time
import backgroundSubtraction
import freeimage
import threading
import csv
import requests

IMAGE_SIZE = (1280, 1080)
BOILER_AREA = (slice(530,1100), slice(530,590))
DETECTION_AREA = (slice(690,1100), slice(545, 575))
POSITION_AREA = (slice(690,1000), slice(545, 575))  #Align left end of sewer corner to 690, 575
CLEARING_AREA = (slice(530,1000), slice(530, 580)) 
FLUORESCENT_AREA = (slice(500,1100), slice(530, 590))

#Textbelt key = 08a7e7d3335d92542dfec857461cfb14af5e0805HQINVtRpVqvDjo9O2wa2I6tTo

FLUOR_PIXEL_BRIGHT_VALUE = 500		#Used for finding brightness, newer/better way to do this?
DETECTION_THRES = 4
POSITION_THRES = 3
CLEARING_THRES = 3
LOST_CUTOFF = 1.5
DOUBLE_THRESH = 1.3

CYAN_EXPOSURE_TIME = 7  #~20 for lin-4, 50 for mir-71?
GREEN_YELLOW_EXPOSURE_TIME = 50	#for mcherry

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

class MicroDevice(threading.Thread):
    """
    (Super?)class which contains basic functions for running the device, with specific methods defined below.
    """ 
    
    def __init__(self, exp_direct):
        """Initalizes the scope and device 
        """
        self.scope, scope_properties = scope_client.client_main()
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
        self.scope.camera.readout_rate = '280 MHz'
        self.scope.camera.binning = '2x2'
        self.scope.tl.lamp.enabled = True
        
        self.device = iotool.IOTool("/dev/ttyMicrofluidics")
        
        self.file_location = Path(exp_direct)
        self.file_location.mkdir(mode=0o777, parents=True, exist_ok=True)
        self.info = exp_direct.split('/')[-1]

        self.device_clear_tubes()
    
        #Create event flag, initially set to false
        self.running = threading.Event()
        super().__init__(daemon=True)

    def pause(self):
        """Sets 'running' flag to false
        """
        self.running.clear()

    def resume(self):
        """Sets 'running' flag to true
        """
        self.running.set()

    def quit(self):
        print('Quiting')
        self.quitting = True
        
    def clear(self):
        print('Cleared')
        self.cleared = True

    def reset(self):
        print('Reset')
        self.reset = True
		
    def write_csv_line(self,csv,data):
        csv.write(','.join(map(str, data)) + '\n')
        
    def set_scope(self):
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
        self.scope.camera.readout_rate = '280 MHz'
        self.scope.camera.binning = '2x2'
        self.scope.tl.lamp.enabled = True

    def lamp_off(self):
        """
        Turns the lamp off
        """
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
        """
        Command that stops the device from sorting or loading worm.
        Device is set to a safe steady state
        """
        self.scope.camera.end_image_sequence_acquisition()
        self.device.execute(PUSH_CHANNEL_PRESSURE, SEWER_CHANNEL_PRESSURE,
                            UP_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_PRESSURE,
                            DOWN_CHANNEL_PRESSURE)
        
    def device_start_load(self):
        """
        Command that toggles the device to being loading worms into the device and 
        will continue to push worms into the device until given another command
        """
        self.device.execute(PUSH_CHANNEL_STATIC, SEWER_CHANNEL_SUCK,
                            UP_CHANNEL_PRESSURE,STRAIGHT_CHANNEL_PRESSURE,
                            DOWN_CHANNEL_PRESSURE)

    def device_stop_load(self):
        """
        Command that toggles the device to stop loading new worms
        """
        self.device.execute(PUSH_CHANNEL_PRESSURE)

    def device_clear_tubes(self):
        """
        Command that toggles the deivce to set all tubes to push water and hopefully
        clear the tubes of debris
        """
        self.device.execute(SEWER_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_PRESSURE,
                            UP_CHANNEL_PRESSURE, DOWN_CHANNEL_PRESSURE,
                            PUSH_CHANNEL_PRESSURE)

    def save_image(self, image, name, worm_count):
        save_location = str(self.file_location) + '/' + name + '_' + str(worm_count) + '.png'
        freeimage.write(image, save_location, flags=freeimage.IO_FLAGS.PNG_Z_BEST_SPEED)

    def set_background_areas(self):
        """
        Function that sets the background values for areas of interest in an
        list of numpy arrays
        """
        raise NotImplementedError('No sorting method given')

    def set_bf_background(self, worm_count):        #Why does this take two images while fluor bgs only take one?
        """
        Function that sets the background values for areas of interest in an
        list of numpy arrays
        """
        self.background = self.capture_image(self.bright)
        time.sleep(PICTURE_DELAY)
        self.background_2 = self.capture_image(self.bright)
        subtracted = abs(self.background.astype('int32') - self.background_2.astype('int32'))
        self.detect_background = numpy.sum(subtracted[DETECTION_AREA])
        self.positioned_background = numpy.sum(subtracted[POSITION_AREA])
        self.clear_background = numpy.sum(subtracted[CLEARING_AREA])
        self.save_image(self.background, 'brightfield_background', worm_count)

    def set_cyan_background(self, worm_count):
        self.cyan_background = self.capture_image(self.cyan)
        self.save_image(self.cyan_background, 'cyan_background', worm_count)

    def set_green_yellow_background(self, worm_count):
        self.green_yellow_background = self.capture_image(self.green_yellow)
        self.save_image(self.green_yellow_background, 'green_yellow_background', worm_count)

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
        #self.save_image(floored_image.astype('uint16'), 'worm_mask', True)
        a = (-1*numpy.sum(mask))
        return (a)

    def worm_mask(self, subtracted_image):      #Make this more consistent/make everything go through this
        floored_image = backgroundSubtraction.percentile_floor(subtracted_image, .99)
        floored_image[self.boiler] = 0
        backgroundSubtraction.clean_dust_and_holes(floored_image)
        #return floored_image.astype('bool')
        return floored_image

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
                self.reset = False
                time.sleep(SORTING_INTERVAL)
                self.device.execute(UP_CHANNEL_PRESSURE)
                self.up_count += 1
        elif direction == 'down':
            self.device.execute(SEWER_CHANNEL_PRESSURE, DOWN_CHANNEL_SUCK, PUSH_CHANNEL_STATIC)
            if self.check_cleared(background, worm_count):
                self.reset = False
                time.sleep(SORTING_INTERVAL)
                self.device.execute(DOWN_CHANNEL_PRESSURE)
                self.down_count += 1
        elif direction == 'straight':
            self.device.execute(SEWER_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_SUCK, PUSH_CHANNEL_STATIC)
            if self.check_cleared(background, worm_count):
                self.reset = False
                time.sleep(SORTING_INTERVAL)
                self.device.execute(STRAIGHT_CHANNEL_PRESSURE)
                self.straight_count += 1

    def check_cleared(self, background, worm_count):        #TODO: Fix force reset option, because flag never triggers
                                                            #how to make this better? Reject worm after too long?
        time_clear_start = time.time()
        message_sent = False
        while not self.reset:
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
                if time_stuck - time_clear_start > 60 and message_sent == False:
                    self.send_text(message='Sort stuck trying to clear')
                    print('Sent text: Device stuck')
                    message_sent = True

        if self.reset == True: #Use if stuck in loop, or if background reset required
            self.reset_background(worm_count)
            return True

    def clear_double_worms(self):
        self.device.execute(SEWER_CHANNEL_PRESSURE,
                    UP_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_SUCK,
                    DOWN_CHANNEL_PRESSURE)
        time.sleep(0.5)

    def send_text(self, message):
        requests.post('https://textbelt.com/text',
                        {'phone': '6019538192',
                         'message': str(message),
                         'key':'08a7e7d3335d92542dfec857461cfb14af5e0805HQINVtRpVqvDjo9O2wa2I6tTo'})

    def main(self):    #Not entirely sure what this is
        self.device = iotool.IOTool("/dev/ttyMicrofluidics")
        return

    def check_size(self, worm_size, max_size, min_size):
        #Is there a way I could write this function here so it didn't clutter the code bellow?
        #Problem is that the bellow code is reliant on breaking the loop in the case that a worm is rejected--how to
        #return the same effect in a function?
        pass

    def find_95th_fluor_amount(self, subtracted_image, worm_mask):
        fluor_worm = subtracted_image[worm_mask]
        fluor_95th = numpy.percentile(fluor_worm, 95)
        return fluor_95th

    def find_mean_fluor_amount(self, subtracted_image, worm_mask):
        fluor_worm = subtracted_image[worm_mask]
        mean_fluor = numpy.mean(fluor_worm)
        return mean_fluor

    def build_hist(self):
        """Function that builds initial histogram from which percentiles are taken to determine upper and lower thresholds
        for sorting. Also returns sorting parameter (fluorescence, length, etc) as list which can be updated.
        """

        self.hist_values = []
        #TODO: would it not make more sense to save each parameter (fluor, size, etc) as an array and then put them all
        #together in the summary file? Then I wouldn't need a separate list for hist values, would just pull from all 
        #fluors taken

        self.sort(calibration = True)
        #Data should just be put into the csv file
        #How to determine initial min and max sizes, to reject progeny/doubled worms without a good sample?
        self.upper_threshold = numpy.percentile(self.hist_values, 90)
        self.lower_threshold = numpy.percentile(self.hist_values, 10)

        return self.hist_values, self.upper_threshold, self.lower_threshold


    def sort(self, calibration):      #TODO: break this into general sorting function, put other relevant details in run func.
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

        worm_count = 0
        self.up_count = 0
        self.down_count = 0
        self.straight_count = 0
        cycle_count = 0

        #1 Loading Worms
        self.device_start_load()
        time_start = time.time()
        time_seen = time_start #Initial setting for time_seen before worm found

        #2 Detect Worms
        try:
            while True:

                cycle_count += 1
                """
                if calibration:
                    if worm_count >= 100:   #Better way than to hard code 100 worms? passing too many arguments?
                        print('Finished initial histogram')
                        return self.hist_values
                        break
                """
                time_without_worm = time.time()
                if time_without_worm - time_seen > 120 and message_sent == False:
                    self.send_text(message='No worm for 2 min')
                    print('Warning sent: No worm')
                    message_sent = True

                current_image = self.capture_image(self.bright)

                if self.detect_worms(current_image, self.background):
                    print(' ')
                    print(' ')
                    print('Worm has been detected')
                    time_seen = time.time()
                    time_between_worms = time_seen - time_start
                    detected_image = current_image

                    #3 Stop worms
                    self.device_stop_load()

                    #4 Position worms
                    while True:
                        current_image = self.capture_image(self.bright)
                        if self.lost_worm(current_image, self.background):
                            print('Worm was lost')
                            self.device.execute(SEWER_CHANNEL_PRESSURE)
                            time.sleep(.1)
                            break

                        elif self.positioned_worm(current_image, detected_image):   #Does this do it's job? Sometimes rejects big worms
                            #Maybe include some param that says if a worm is too close to the edge, let it get to the center 
                            difference_between_worm_background = (abs(current_image.astype('int32') - self.background.astype('int32')))
                            worm_size = self.size_of_worm(difference_between_worm_background)
                            print('Size of worm before sorting: ' + str(worm_size))
                            
                            if worm_size > self.size_threshold:     #TODO: think on how to better decide size thresholds
                                print('Detected Double Worm')
                                self.save_image(current_image, 'doubled worm_analyze', worm_count)
                                self.device_sort('straight', self.background, worm_count)
                                print('Doubled worms sorted Straight')
                                time.sleep(0.5)
                                break

                            elif worm_size < self.min_worm_size:     #TODO: look back and decide if problem
                                                                     #Probably actually rejecting large (slow) worms here
                                print('Detected small worm')
                                self.save_image(current_image, 'small worm_analyze', worm_count)
                                self.device_sort('straight', self.background, worm_count)
                                time.sleep(0.5)
                                break

                            worm_count += 1
                            print('Worm positioned')

                            if calibration:
                                self.save_image(current_image, 'calibration_worm', worm_count)
                                print('Image saved')

                                sort_param, direction = self.analyze(self.background, current_image, worm_count, calibration = True)
                                print('Calibration worm number: ' + str(worm_count))

                                #self.hist_values.append(sort_param)  #building initial histogram
                                #How to make this generalizable if not using a histogram?
                            elif not calibration:
                                #5 Save image
                                self.save_image(current_image, 'positioned', worm_count)    #TODO: better names?
                                print('Image saved')

                                #6 Analyze worms
                                sort_param, direction = self.analyze(self.background, current_image, worm_count, calibration = False)
                                print('Worm number: ' + str(worm_count))

                            #Second size check to make sure new worms haven't shown up 
                            #TODO: can probably put this in function

                            post_analysis_image = self.capture_image(self.bright)
                            difference_between_worm_background = (abs(post_analysis_image.astype('int32') - self.background.astype('int32')))
                            worm_size_2 = self.size_of_worm(difference_between_worm_background)
                            print('Size of worm after analysis: ' + str(worm_size_2))
                            #Maybe instead of a second size threshold test, check if size has changed appreciably

                            if worm_size_2 > self.size_threshold:
                                print('Detected Double Worm')
                                self.save_image(current_image, 'doubled worm_analyze', worm_count)
                                self.device_sort('straight', self.background, worm_count)
                                print('Doubled worms sorted Straight')
                                time.sleep(0.5)
                                break

                            elif worm_size_2 < self.min_worm_size:
                                print('Worm was lost')
                                self.save_image(current_image, 'small worm_analyze', worm_count)
                                self.device_sort('straight', self.background, worm_count)
                                time.sleep(0.5)
                                break

                            #7 Sort worms
                            self.device_sort(direction, self.background, worm_count)
                            print('Worm ' + str(worm_count) + ' sorted ' + direction)

                            #8 Move worms??
                            #Old code here has an if statement to check if a reset is required, better place to put this?
                            #Not sure self.cleared was ever being set as true

         
                            #TODO: Add in 'reason'+measurements for worms that are rejected (e.g. for too small, doubled, etc.)
                            #Could just save worm data array, then write out in subsequent function
                            #This might be the easier way to save multiple files to different directories
                            #Currently, worm_data is written worm by worm, and discarded immediately afterward
                            if calibration:
                                note = 'hist'
                            elif not calibration:
                                note = 'sort'
                                print('Up: ' + str(self.up_count), ' Straight: ' + str(self.straight_count), ' Down: ' + str(self.down_count))
                                self.update_hist(sort_param)

                            worm_data = self.generate_data(worm_count, worm_size, sort_param, time_between_worms, direction, note)
                            self.write_csv_line(self.summary_csv, worm_data)

                            break

                        else:
                            detected_image = current_image  #What is this doing?
                        

                    if worm_count % 100 == 0 and worm_count != 0:   #Resets background every 100 worms
                        self.reset_background(worm_count)
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

    def run(self):
        
        print('Running from the superclass!')

        self.setup_csv(self.file_location, self.info)

        self.scope.camera.start_image_sequence_acquisition(frame_count=None, trigger_mode='Software')

        self.boiler = boiler()
        self.reset = False
        self.message_sent = False
        
        #0 Setting Background
        worm_count = 0
        self.set_background_areas(worm_count)
        print('setting backgrounds')
        time.sleep(1)
        #input for build hist?

        self.size_threshold = 7000   #hard coding sizes for now
        self.min_worm_size = 3000

        self.build_hist()

        self.sort(calibration = False)

        self.summary_csv.close()

    """
    Some pseudocode:
    def run(self):
        build_hist()    #class-specific, calls sort and analyze for worms 1-100
            #ask if build hist, check for variables, or input variables
            #Maybe even ask if hist is acceptable, putting up hist, images, etc (detection of extremes?)
            #Could even implement code to a) rebuild hist or b) add another 100 worms 
        sort()          #Still using class-specific analysis function
            #sort should only look for worms, then ask analyze what to do
        save_data()     #Should also be class-specific, as params saved change from case to case
            #Does this mean that if sort exits badly data can't be saved? Does data have to be saved continuously?
    """

class GFP(MicroDevice):

    def setup_csv(self, file_location, info):
        self.summary_csv_location = file_location.joinpath('summary_csv' + info + '.csv')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'fluorescence', 'time', 'direction', 'note']
        self.summary_csv.write(','.join(header) + '\n')

    def set_background_areas(self, worm_count):
        self.set_bf_background(worm_count)
        self.set_cyan_background(worm_count)
        self.lamp_off()
        self.scope.tl.lamp.enabled = True
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME

    def manual_set_up(self):
        self.upper_threshold = int(input('Upper GFP threshold = '))
        self.lower_threshold = int(input('Lower GFP threshold = '))
        self.size_threshold = int(input('Max size threshold = '))
        self.min_worm_size = int(input('Min size threshold = '))

    def analyze(self, background, worm_image, worm_count, calibration=False):
        """Class-specific method to determine measurement being sorted.
        Takes background image, bf image of worm, worm count, and returns measurement + sort direction
        Saves fluorescent (GFP) image
        """
        gfp_fluor_image = self.capture_image(self.cyan)
        gfp_subtracted = abs(gfp_fluor_image.astype('int32') - self.cyan_background.astype('int32'))
        worm_subtracted = abs(worm_image.astype('int32') - background.astype('int32'))
        worm_mask = self.worm_mask(worm_subtracted).astype('bool')      #Old function returned as bool automatically, see if bool works with size-finding

        worm_fluor = self.find_95th_fluor_amount(gfp_subtracted, worm_mask)
        print('GFP value = ' + str(worm_fluor))

        if calibration:
            self.save_image(worm_mask.astype('uint16'), 'calibration_worm_mask', worm_count)
            self.save_image(gfp_fluor_image, 'calibration_worm_fluor', worm_count)
            direction = 'straight'
            return worm_fluor, direction

        elif not calibration:
            self.save_image(worm_mask.astype('uint16'), 'worm_mask', worm_count)
            self.save_image(gfp_fluor_image, 'fluor_gfp', worm_count)

            if worm_fluor >= self.upper_threshold:
                direction = 'up'
            elif worm_fluor <= self.lower_threshold:
                direction = 'down'
            else:
                direction = 'straight'

        return worm_fluor, direction

    def update_hist(self, fluor_amount):
        self.hist_values.append(fluor_amount)
        self.upper_threshold = numpy.percentile(self.hist_values, 90)
        self.loewr_threshold = numpy.percentile(self.hist_values, 10)
        print(self.upper_threshold, self.lower_threshold)

    def generate_data(self, worm_count, size, sort_param, time, direction, note):
        """Function for returning the right csv lin config for a given class.
        Here, returns the number, size, fluorescence, time between worms, directions, and note (hist or sort)
        sort_param = GFP fluorescence
        """

        worm_data = [worm_count, size, sort_param, time, direction, note]
        return worm_data


class Autofluor(MicroDevice):

    def setup_csv(self, file_location, info):
        self.summary_csv_location = file_location.joinpath('summary_csv' + info + '.csv')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'fluorescence', 'time', 'direction', 'note']
        self.summary_csv.write(','.join(header) + '\n')

    def set_background_areas(self, worm_count):
        self.set_bf_background(worm_count)
        self.set_green_yellow_background(worm_count)
        self.lamp_off()
        self.scope.tl.lamp.enabled = True
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME

    def analyze(self, background, worm_image, worm_count, calibration=False):
        """
        """
        tritc_image = self.capture_image(self.green_yellow)
        tritc_subtracted = abs(tritc_image.astype('int32')- self.green_yellow_background.astype('int32'))
        worm_subtracted = abs(worm_image.astype('int32') - background.astype('int32'))
        worm_mask = self.worm_mask(worm_subtracted).astype('bool')      #Old function returned as bool automatically, see if bool works with size-finding

        worm_fluor = self.find_95th_fluor_amount(tritc_subtracted, worm_mask)
        print('Autofluorescence value = ' + str(worm_fluor))

        if calibration:
            self.save_image(worm_mask.astype('uint16'), 'calibration_worm_mask', worm_count)
            self.save_image(tritc_image, 'calibration_worm_fluor', worm_count)
            direction = 'straight'
            return worm_fluor, direction

        elif not calibration:
            self.save_image(worm_mask.astype('uint16'), 'worm_mask', worm_count)
            self.save_image(tritc_image, 'fluor_gfp', worm_count)

            if worm_fluor >= self.upper_threshold:
                direction = 'up'
            elif worm_fluor <= self.lower_threshold:
                direction = 'down'
            else:
                direction = 'straight'

        return worm_fluor, direction

    def generate_data(self, worm_count, size, sort_param, time, direction, note):
        worm_data = [worm_count, size, sort_param, time, direction, note]
        return worm_data

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

    def set_background_areas(self, worm_count):
        self.set_bf_background(worm_count)
        self.set_cyan_background(worm_count)
        self.set_green_yellow_background(worm_count)
        self.lamp_off()
        self.scope.tl.lamp.enabled = True
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME

    def analyze(self, background, worm_image, worm_count, calibration=True):
        """
        """
        worm_subtracted = abs(worm_image.astype('int32') - background.astype('int32'))
        worm_mask = self.worm_mask(worm_subtracted).astype('bool')
        self.save_image(worm_mask.astype('uint16'), 'worm_mask', worm_count)

        cyan_fluor_image = self.capture_image(self.cyan)
        cyan_subtracted = abs(cyan_fluor_image.astype('int32')- self.cyan_background.astype('int32'))
        self.save_image(cyan_fluor_image, 'cyan_worm', worm_count)

        cyan_95th = self.find_95th_fluor_amount(cyan_subtracted, worm_mask)
        cyan_mean = self.find_mean_fluor_amount(cyan_subtracted, worm_mask)

        tritc_fluor_image = self.capture_image(self.green_yellow)
        tritc_subtracted = abs(tritc_fluor_image.astype('int32') - self.green_yellow_background.astype('int32'))
        self.save_image(tritc_fluor_image, 'tritc_worm', worm_count)

        tritc_95th = self.find_95th_fluor_amount(tritc_subtracted, worm_mask)
        tritc_mean = self.find_mean_fluor_amount(tritc_subtracted, worm_mask)

        direction = 'straight'

        sort_param = [cyan_95th, cyan_mean, tritc_95th, tritc_mean]
        
        print(str(sort_param), direction)

        return sort_param, direction

    def generate_data(self, worm_count, size, sort_param, time, direction, note):
        """sort_param is a list here, necessitating specific function. Direction and note aren't really salient here.
        """
        worm_data = [worm_count, size] + sort_param + [time, direction, note]
        return worm_data

    #Need specific funciton for running worms throuhg here?
    """
    def nosort(self):
        self.setup_csv(self.file_location, self.info)
        
        self.size_threshold = 7500   #hardcoding sizes for now
        self.min_worm_size = 2500
        
        self.sort(calibration=True)
        self.summary_csv.close()
    """
    
    def nosort(self):

        self.setup_csv(self.file_location, self.info)

        self.scope.camera.start_image_sequence_acquisition(frame_count=None, trigger_mode='Software')

        self.boiler = boiler()
        self.reset = False
        self.message_sent = False
        
        #0 Setting Background
        worm_count = 0
        self.set_background_areas(worm_count)
        print('setting backgrounds')
        time.sleep(1)
        #input for build hist?

        self.size_threshold = 7500   #hard coding sizes for now
        self.min_worm_size = 2500

        self.sort(calibration = True)

        self.summary_csv.close()        






