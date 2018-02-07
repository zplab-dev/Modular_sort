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
import requests
#import IPython
#import prompt_toolkit


#Setting useful constants
#Areas for image anaylsis
#old
IMAGE_SIZE = (1280, 1080)
BOILER_AREA = (slice(530,1100), slice(530,590))
DETECTION_AREA = (slice(690,1100), slice(545, 575))
POSITION_AREA = (slice(690,1000), slice(545, 575))  #Align left end of sewer corner to 690, 575
CLEARING_AREA = (slice(530,1000), slice(530, 580))
FLUORESCENT_AREA = (slice(500,1100), slice(530, 590))
#DOUBLE_AREA = (slice(670,870), slice(570,580))
#TOP_PUSHER_CHANNEL_AREA = (slice(580,830),slice(55,300))
#STRAIGHT_PUSHER_CHANNEL_AREA = (slice(70,30),slice(350,500))
#BOTTOM_PUSHER_CHANNEL_AREA = (slice(580,830), slice(780,1030))
#BUBBLE_AREA = (slice(800,1000),slice(376,500))
#Textbetl key = 08a7e7d3335d92542dfec857461cfb14af5e0805HQINVtRpVqvDjo9O2wa2I6tTo

FLUOR_PIXEL_BRIGHT_VALUE = 500
DETECTION_THRES = 4
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

CYAN_EXPOSURE_TIME = 7  #7 for lin-4, 50 for mir-71?
YELLOW_EXPOSURE_TIME = 8
BRIGHT_FIELD_EXPOSURE_TIME = 2

LIGHT_DELAY = .05
PICTURE_DELAY = .01
SORTING_INTERVAL = .4
MAX_SORTING_TIME = 1.5

BACKGROUND_REFRESH_RATE = 100000
PROGRESS_RATE = 100

MIN_GFP_THRESH = 400 #Will need to reset to account for brighter exposure

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


#Setting file locations for saving images

#File_location = '/mnt/iscopearray/Nonet_Tim/Test##_##_##'

#def ip_input(message=''):
#    ip = IPython.get_ipython()
#    el = prompt_toolkit.shortcuts.create_eventloop(ip.inputhook)
#    return prompt_toolkit.prompt(message, eventloop=el)


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
        self.scope.tl.lamp.enabled = True

        self.device = iotool.IOTool("/dev/ttyMicrofluidics")

        self.file_location = Path(exp_direct)
        self.file_location.mkdir(mode=0o777, parents=True, exist_ok=True)
        self.date = exp_direct.split('/')[-1]
        
        self.summary_location = self.file_location.joinpath('summary_' + self.date + '.txt')
        self.summary_statistics = open(str(self.summary_location),'w')
        
        #summary_csv_location = self.file_location.joinpath('summary_csv' + self.date + '.csv')   #TODO: change to .csv, fix path/writing
        #self.summary_csv = open(str(summary_csv_location), 'w')
        #header = ['worm_number', 'size', 'fluorescence', 'time', 'direction'] #TODO: add reason and such
        #self.summary_csv.write(','.join(header) + '\n')

        #self.data_location = self.file_location.joinpath('wormdata.csv')   #old data saving, never implemented

        self.device_clear_tubes()

        """
        device_type = input('Which type of device is being used? [o] for Old, [n] for new?: ')
        if device_type == 'o':
            print('Using Old Set up')
        if device_type == 'n':
            print('Using New Set up')
            global IMAGE_SIZE, BOILER_AREA, DETECTION_AREA, POSITION_AREA
            global CLEARING_AREA, FLUORESCENT_AREA, DOUBLE_AREA
            IMAGE_SIZE = (1280, 1080)
            BOILER_AREA = (slice(480,1000), slice(540,585))
            DETECTION_AREA = (slice(500,1045), slice(550, 580))
            POSITION_AREA = (slice(510,930), slice(545, 575))
            CLEARING_AREA = (slice(375,750), slice(535, 590))
            FLUORESCENT_AREA = (slice(480,1000), slice(540,585))
            DOUBLE_AREA = (slice(480,1000), slice(540,585))
            print(DETECTION_AREA)
        """

        worm_count = 0
        self.up = 0
        self.down = 0
        self.straight = 0
        self.TIME_FOR_PUSHING = .8
        #worm_data = list()

        #Pausing stuff
        self.running = True
        super().__init__(daemon=True)
        self.quitting = False
        self.cleared = False

        self.resume() #(testing)

    def write_csv_line(self,csv,data):
        csv.write(','.join(map(str, data)) + '\n')

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
                            DOWN_CHANNEL_PRESSURE)

    def device_start_load(self):
        """
        Command that toggles the device to being loading worms into the device and
        will continue to push worms into the device until given another command
        """
        self.device.execute(PUSH_CHANNEL_STATIC, SEWER_CHANNEL_SUCK,
                            UP_CHANNEL_PRESSURE,STRAIGHT_CHANNEL_PRESSURE,
                            DOWN_CHANNEL_PRESSURE)


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
        elif direction == 'down':
            self.device.execute(SEWER_CHANNEL_PRESSURE, DOWN_CHANNEL_SUCK, PUSH_CHANNEL_STATIC)
            if self.check_cleared(background, worm_count):
                self.reset = False
                time.sleep(SORTING_INTERVAL)
                self.device.execute(DOWN_CHANNEL_PRESSURE)
        elif direction == 'straight':
            self.device.execute(SEWER_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_SUCK, PUSH_CHANNEL_STATIC)
            if self.check_cleared(background, worm_count):
                self.reset = False
                time.sleep(SORTING_INTERVAL)
                self.device.execute(STRAIGHT_CHANNEL_PRESSURE)


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
                                    PUSH_CHANNEL_PRESSURE)
                time.sleep(1)
                self.device.execute(SEWER_CHANNEL_SUCK,
                                    PUSH_CHANNEL_STATIC,
                                    STRAIGHT_CHANNEL_PRESSURE,
                                    UP_CHANNEL_PRESSURE,
                                    DOWN_CHANNEL_PRESSURE)
        except KeyboardInterrupt:
            self.device.execute(PUSH_CHANNEL_PRESSURE,
                                SEWER_CHANNEL_PRESSURE,
                                UP_CHANNEL_PRESSURE,
                                STRAIGHT_CHANNEL_PRESSURE,
                                DOWN_CHANNEL_PRESSURE)


    def device_clear_tubes(self):
        """
        Command that toggles the deivce to set all tubes to push water and hopefully
        clear the tubes of debris
        """
        self.device.execute(SEWER_CHANNEL_PRESSURE, STRAIGHT_CHANNEL_PRESSURE,
                            UP_CHANNEL_PRESSURE, DOWN_CHANNEL_PRESSURE,
                            PUSH_CHANNEL_PRESSURE)


    def device_stop_load(self):
        """
        Command that toggles the device to stop loading new worms
        """
        self.device.execute(PUSH_CHANNEL_PRESSURE)

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
        type_of_image()
        self.scope.camera.send_software_trigger()
        return self.scope.camera.next_image()

    def save_image(self, image, name, worm_count):
        save_location = str(self.file_location) + '/' + name + '_' + str(worm_count) + '.png'
        freeimage.write(image, save_location, flags=freeimage.IO_FLAGS.PNG_Z_BEST_SPEED)

    def set_background_areas(self):
        """
        Function that sets the background values for areas of interest in an
        list of numpy arrays
        """
        raise NotImplementedError('No sorting method given')

    def detect_worms(self, current_image, background):
        """
        Function that returns if a worm has been detected
        """
        #print('Detected Value:' + str(numpy.abs(numpy.sum(current_image[DETECTION_AREA].astype('int32') - background[DETECTION_AREA].astype('int32')))))
        #print('Required Value:' + str( 5 * self.detect_background))

        return ((numpy.sum(numpy.abs(current_image[DETECTION_AREA].astype('int32')
                                     - background[DETECTION_AREA].astype('int32')))
                 -self.detect_background) > DETECTION_THRES  * self.detect_background)

    def lost_worm(self, current_image, background):
        """
        Function that determines if a worm was lost
        Worm is deciced lost because image is close enough to background.
        """
        worm_visibility = abs(current_image[POSITION_AREA].astype('int32')
                              - background[POSITION_AREA].astype('int32'))
        return ((numpy.sum(worm_visibility) -
                 self.positioned_background) < LOST_CUTOFF * self.positioned_background)

    def positioned_worm(self, current_image, detected_image):
        """
        Function that determines if a worm has been positioned
        The worm is positioned because the change is small.
        """
        worm_movment = abs(current_image[POSITION_AREA].astype('int32') - detected_image[POSITION_AREA].astype('int32'))
        #print('detected value ' + str(numpy.sum(worm_movment)))
        #print('Needed value   ' + str(POSITION_THRES * self.positioned_background))
        return  ((numpy.sum(worm_movment) - self.positioned_background) < POSITION_THRES * self.positioned_background)

    def size_of_worm(self, subtracted_image):
        """
        Function that saves an image of the mask of a worm and returns the size of mask (in number of pixels)
        Depending on how background.backgroundSubtraction.clean_dust_and_holes(image) works this function might neeed to be modified so that it returns the appropriate image.
        Currently the clean_dust_and_holes does not return the clean image and actually modifies the passed image.
        """

        floored_image = backgroundSubtraction.percentile_floor(subtracted_image, .99)
        floored_image[self.boiler] = 0
        backgroundSubtraction.clean_dust_and_holes(floored_image)
        #self.save_image(floored_image.astype('uint16'), 'worm_mask', True)
        a = (-1*numpy.sum(floored_image))
        return (a)

    def worm_mask(self, subtracted_image):
        floored_image = backgroundSubtraction.percentile_floor(subtracted_image, .99)
        floored_image[self.boiler] = 0
        backgroundSubtraction.clean_dust_and_holes(floored_image)
        return floored_image.astype('bool')

    def analyze(self):
        """
        function that tells the device what sorting/analzying method to use:
        Is overwritten by a super class
        """
        raise NotImplementedError('No sorting method given')

    """                                                                         #Under construction
    def clear_worms(self, background, photos=False, calibration=False):
        ""
        Function that moves the worm from the sewere gratting and attempts
        to take a picutre of the worm

        Possible to change this so it only sorts until the worm has left the current frame.
        ""
        push_time = time.time()
        current_time = time.time()
        sucking = False
        pushing = False
        self.cleared = False
        try:
            while not self.cleared:
                current_time = time.time()
               #if self.worm_direction:
               #Have the direction the worm is being sorted in  flutter to prevent a worm from being stuck to the
               #directional channel when it is pulling a vacuum.
               #    directon =
               #    if self.worm_direction == 'up':
               #
               #    elif sort_direction
               #    self.device.execute(
                if (current_time - push_time > self.TIME_FOR_PUSHING
                    and not sucking):
                    sucking = True
                    print('Sewer Channel set to SUCK before sorted worm is cleared')
                    self.TIME_FOR_PUSHING += .01
                    self.device.execute(SEWER_CHANNEL_SUCK)
                if (current_time - push_time > MAX_SORTING_TIME and
                    not pushing):
                    print('Wormed has failed to clear well')
                    pushing = True
                    self.device.execute(SEWER_CHANNEL_PRESSURE)
                    self.device.execute(PUSH_CHANNEL_PRESSURE)
                if pushing:
                    print('still pushing')

                current_image = self.capture_image(self.bright)
                sorted_worm_difference = abs(current_image[CLEARING_AREA].astype('int32')
                                             - background[CLEARING_AREA].astype('int32'))
                if numpy.sum(sorted_worm_difference) < CLEARING_THRES * self.clear_background:
                    #Worm has been deteremiend cleared as the image is close eough to the background.
                    if not photos:
                        print('Worm sorted ' + str(self.worm_direction))
                    if photos:
                        print('Picutre of worm being sorted ' + str(self.worm_direction))
                        if calibration:
                            self.save_image(current_image, 'calibration_worm_sent', True)
                        else:
                            self.save_image(current_image, 'sent', True)
                    time.sleep(SORTING_INTERVAL)
                    break
        except KeyboardInterrupt:
            pass
    """

    def check_cleared(self, background, worm_count):
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
        time.sleep(1)

    def reset_background(self, worm_count):
        background = self.capture_image(self.bright)
        self.save_image(background, 'new_background_worm_', worm_count)
        self.set_background_areas(worm_count)
        print('resetting backgrounds')
        time.sleep(1)

    def send_text(self, message):
        requests.post('https://textbelt.com/text',
                        {'phone': '6019538192',
                        'message': str(message),
                        'key':'08a7e7d3335d92542dfec857461cfb14af5e0805HQINVtRpVqvDjo9O2wa2I6tTo'})

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
        self.summary_csv_location = self.file_location.joinpath('summary_csv' + self.date + '.txt')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'fluorescence', 'time', 'direction']
        self.summary_csv.write(','.join(header) + '\n')

        self.scope.camera.start_image_sequence_acquisition(
            frame_count=None, trigger_mode='Software')

        worm_count = 0
        cycle_count = 0
        self.boiler = boiler()

        #0 Setting Background
        background = self.capture_image(self.bright)
        self.save_image(background, 'background', worm_count)
        self.set_background_areas(worm_count)
        print('setting backgrounds')
        time.sleep(1)
        self.reset = False

        #1 Loading Worms
        self.device_start_load()
        time_start = time.time()
        time_seen = time_start #Initial setting for time_seen before worm found
        message_sent = False

        #2 Detect Worms
        time_between_worms = list()
        time_to_position_worms = list()
        try:
            while True:
                if self.quitting:
                    self.device_clear_tubes()
                    break
                while not self.running:
                    paused_image = self.capture_image(self.bright)
                    time.sleep(.05)

                cycle_count += 1
                time_without_worm = time.time()
                
                if time_without_worm - time_seen > 120 and message_sent == False:
                    self.send_text(message='No worm for 2 min')
                    print('Warning sent: No worm')
                    message_sent = True

                current_image = self.capture_image(self.bright)
                if self.detect_worms(current_image, background):
                    print(' ')
                    print(' ')
                    print('Worm has been detected')
                    message_sent = False
                    #Worm is detected because of a significant change from background.
                    time_seen = time.time()
                    detected_image = current_image

                    #3 stop worms
                    self.device_stop_load()

                    #4 position worms
                    while True:
                        if self.quitting or not self.running:
                            self.clear_double_worms()
                            break
                        current_image = self.capture_image(self.bright)
                        if self.lost_worm(current_image, background):
                            print('Worm has been lost')
                            self.device.execute(SEWER_CHANNEL_PRESSURE)
                            time.sleep(.1)
                            self.summary_statistics.write("\nWorm " + str(worm_count) + "was lost")
                            break
                        elif self.positioned_worm(current_image, detected_image):
                            difference_between_worm_background = (abs(current_image.astype('int32') - background.astype('int32')))
                            worm_size = self.size_of_worm(difference_between_worm_background)
                            print('Size of worm before sorting: ' + str(worm_size))
                            if worm_size > self.size_threshold:
                                print('Detected Double Worm')
                                self.save_image(current_image, 'doubled worm_analyze', worm_count)
                                self.summary_statistics.write( '\n doubled worm size of: '+ str(worm_size))
                                self.device_sort('straight', background, worm_count)
                                self.worm_direction = 'straight'
                                self.reason = 'double_worm'
                                self.summary_statistics.write("Straight\n")
                                print('Doubled worms sorted Straight')
                                time.sleep(0.5)
                                break
                            elif worm_size < self.min_worm_size:
                                print('Detected small worm')
                                self.save_image(current_image, 'small worm_analyze', worm_count)
                                self.summary_statistics.write( '\n small worm size of: ' + str(worm_size))
                                self.device_sort('straight', background, worm_count)
                                self.worm_direction = 'straight'
                                self.reason = 'small_worm'
                                self.summary_statistics.write('Straight\n')
                                time.sleep(0.5)
                                break
                            worm_count += 1
                            #self.worm_data = list()
                            #self.worm_data.append(worm_count)
                            #self.worm_data.append(worm_size)
                            time_positioned = time.time()
                            time_between_worms.append(time_seen - time_start)
                            time_to_position_worms.append(time_positioned - time_seen)
                            self.summary_statistics.write('\n' + str(worm_count)
                                                          + "\nTime Between Worms: "
                                                          + "{:10.5}".format(str(time_seen - time_start))
                                                          + "\n Worm size = :"
                                                          + str(worm_size)
                                                          + '\n')
                            #time_start = time_seen
                            print('Worm has been positioned')

                            #5 save worm image

                            self.save_image(current_image, 'position', worm_count)
                            #self.save_image(difference_between_worm_background.astype('uint16'), 'worm_subtracted', True)
                            print('Images saved')

                            #6/7 analyze worms/sort
                            gfp_amount, self.worm_direction = self.analyze(background, current_image, worm_count, calibration = False)
                            print('Worms analyzed and sorted')
                            print('Worm number: ' + str(worm_count))

                            #8 move worms
                            #self.clear_worms(background)                #Under construction
                            if self.cleared:                             # Not sure this is doing anything
                                self.device.execute(SEWER_CHANNEL_PRESSURE, UP_CHANNEL_PRESSURE,
                                                    STRAIGHT_CHANNEL_PRESSURE,DOWN_CHANNEL_PRESSURE)
                                time.sleep(2)
                                background = self.capture_image(self.bright)
                                self.save_image(background, 'background' + str(worm_count))
                                self.set_background_areas()
                                print('Reset background')

                            print('Sorted')
                            time_sorted = time.time()
                            self.summary_statistics.write("Time to Sorted: "
                                                          + "{:10.5}".format(str(time_sorted - time_positioned))
                                                          + "\n")
                            
                            time_between = time_seen - time_start

                            worm_data = [worm_count, worm_size, gfp_amount, time_between, self.worm_direction]     #TODO: Add in 'reason'+measurements for worms that are rejected (e.g. for too small, doubled, etc.)
                            self.write_csv_line(self.summary_csv, worm_data)
                            break
                        else:
                            detected_image = current_image
                        #9 --> 1

                    if worm_count % 100 == 0 and worm_count != 0:  #Periodically resets background
                        self.reset_background(worm_count)     #Need to make sure this is not happening while a worm is in view

                    self.device_start_load()

                elif cycle_count % PROGRESS_RATE == 0:
                    print(str(PROGRESS_RATE) + ' Cycles')

                elif cycle_count % BACKGROUND_REFRESH_RATE == 0:
                    print(str(cycle_count) + ' Cycles Reseting Background')
                    self.device_stop_load()
                    self.device.execute(SEWER_CHANNEL_PRESSURE,
                                        STRAIGHT_CHANNEL_PRESSURE,
                                        UP_CHANNEL_PRESSURE,
                                        DOWN_CHANNEL_PRESSURE,
                                        PUSH_CHANNEL_PRESSURE)
                    #Clears the device without required a human input
                    time.sleep(3)
                    background = self.capture_image(self.bright)
                    self.set_background_areas()
                    self.device_start_load()

        except KeyboardInterrupt:
            pass
        finally:
            self.summary_statistics.write('\n Average worm detection time :'
                                          + str(numpy.mean(time_between_worms))
                                          + '\n Average worm positioning time :'
                                          + str(numpy.mean(time_to_position_worms)))
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
        time.sleep(PICTURE_DELAY)
        image2 = self.capture_image(self.bright)
        subtracted = abs(image1.astype('int32') - image2.astype('int32'))
        self.detect_background = numpy.sum(subtracted[DETECTION_AREA])
        self.positioned_background = numpy.sum(subtracted[POSITION_AREA])
        self.clear_background = numpy.sum(subtracted[CLEARING_AREA])

    def analyze(self, background):
        """
        function that tells the device what sorting/analzying method to use:
        Is overwritten by a super class
        """
        self.device_sort('straight', background)
        self.worm_direction = 'straight'

class Alternate(MicroDevice):

    size_threshold = 6800

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

    def analyze(self,background, worm_image=False):
        """
        function that tells the device what sorting/analzying method to use:
        Is overwritten by a super class
        """
        direction = worm_count % 3
        dirlst = ['straight','up','down']
        self.device_sort(dirlst[direction])
        self.worm_direction = dirlst[direction]

class Mir71(MicroDevice):

    def set_background_areas(self, worm_count):
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
        self.save_image(self.cyan_background, 'cyan_background', worm_count)

        self.lamp_off()
        self.scope.tl.lamp.enabled = True
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME


    def auto_set_up(self):
        top_threshold = input('What do you want the top X% of expressiont to be = ')
        self.upper_mir71_threshold = int(top_threshold)
        bottom_threshold = input('What do you want the bottom X% of expressiont to be = ')
        self.bottom_mir71_threshold = int(bottom_threshold)
        self.size_threshold = int(input('What do you want the large size threshold to be = '))
        self.min_worm_size = int(input('What do you want the small size threshold to be = '))


    def find_thresholds(self, num_of_worms):
        self.scope.camera.start_image_sequence_acquisition(frame_count=None,
                                                           trigger_mode='Software')

        self.histogram_values_location = self.file_location.joinpath('histogram' + self.date + '.txt')
        self.histogram_values = open(str(self.histogram_values_location),'w')
        csv_location = self.file_location.joinpath('histogram_csv' + self.date + '.txt')
        histogram_csv = open(str(csv_location), 'w')
        header = ['worm_number', 'size', 'fluorescence', 'time', 'direction', 'reason']
        histogram_csv.write(','.join(header) + '\n')

        worm_count = 0
        background = self.capture_image(self.bright)
        self.boiler = boiler()
        self.save_image(background, 'calibration_background', worm_count)
        self.set_background_areas(worm_count)
        print('setting backgrounds')
        time.sleep(1)
        self.size = list()
        self.fluorescence = list()
        self.device_start_load()
        #initial_max_size = int(input('Initial size threshold: '))
        #initial_min_size = int(input('Initial min size threshold: '))
        initial_max_size = 8000
        initial_min_size = 3500

        cycle_count= 0
        time_start = time.time()
        self.reset = False
        try:
            while True:
                cycle_count += 1
                if not cycle_count % PROGRESS_RATE:
                    print(str(PROGRESS_RATE) + ' cycles')
                if worm_count >= num_of_worms:
                    print('Done')
                    break
                current_image = self.capture_image(self.bright)
                if self.detect_worms(current_image, background):
                    print('Worm has been detected')
                    #Worm is detected because of a significant change from background.
                    detected_image = current_image
                    time_seen = time.time() - time_start
                    #Stop Worms
                    self.device_stop_load()
                    print('stopped loading more worms')
                    #Position Worms:
                    while True:
                        current_image = self.capture_image(self.bright)
                        if self.lost_worm(current_image, background):
                            print('Worm has been lost')
                            break
                        elif self.positioned_worm(current_image, detected_image):
                            #worm_quality = input('Is this worm accetable: [y]/n')
                            #if worm_quality != 'y':
                            #    self.lamp_off()
                            #    self.scope.tl.lamp.enabled = True
                            #    self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
                            #    self.device_sort('straight', background)
                            #    time.sleep(2)
                            #    print('Not a worm')
                            #    break
                            worm_count += 1
                            difference_between_worm_background = (abs(current_image.astype('int32') - background.astype('int32')))
                            worm_size = self.size_of_worm(difference_between_worm_background)
                            print('worm size = ' + str(worm_size))
                            if worm_size > initial_max_size:
                                worm_count -= 1
                                print('Detected Double Worm')
                                self.worm_direction = 'straight'
                                reason = 'double'
                                self.device_sort('straight', background, worm_count)
                                #self.clear_worms(background)
                                break
                            if worm_size < initial_min_size:
                                worm_count -= 1
                                print('Small worm')
                                self.worm_direction = 'straight'
                                reason = 'small'
                                self.device_sort('straight', background, worm_count)
                                #self.clear_worms(background)
                                break

                            #worm_mask = self.worm_mask(difference_between_worm_background)
                            #self.save_image(worm_mask.astype('uint16'), 'calibration_worm_mask', True)
                            print('Worm number ' + str(worm_count) + ' out of ' + str(num_of_worms))

                            #self.save_image(worm_mask.astpye('int32'), 'calibration_worm_mask' + str(worm_count))
                            #self.save_image(current_image, 'calibration_worm'+ str(worm_count))
                            self.save_image(current_image, 'calibration_worm', worm_count)
                            print('images_saved')

                            #CYAN
                            #current_image = self.capture_image(self.cyan)
                            #self.save_image(current_image, 'calibration_worm_fluor' + str(worm_count))
                            #gfp_image = abs(current_image.astype('int32')- self.cyan_background.astype('int32'))
                            #gfp_amount = self.find_fluor_amount(gfp_image, worm_mask)
                            #print('GFP amount = ' + str(gfp_amount))

                            gfp_amount = self.analyze(background, current_image, worm_count, calibration = True)
                            if gfp_amount < MIN_GFP_THRESH:
                                worm_count -= 1
                                print('Worm bellow Mir-71 expression threshold')
                                self.worm_direction = 'straight'
                                reason = 'low_GFP'
                                self.device_sort('straight', background, worm_count)
                                #self.clear_worms(background)
                                break

                            print('Worm in appropriate size and fluorescence range')
                            self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
                            self.lamp_off()
                            self.scope.tl.lamp.enabled = True
                            self.size.append(worm_size)
                            self.fluorescence.append(gfp_amount)
                            self.device_sort('straight', background, worm_count)
                            self.worm_direction = 'straight'
                            reason = 'histogram'
                            #self.clear_worms(background,True, True)
                            worm_data = [worm_count, worm_size, gfp_amount, time_seen, self.worm_direction, reason]
                            self.write_csv_line(histogram_csv, worm_data)
                            break

                        else:
                            detected_image = current_image
                            #9 --> 1

                self.device_start_load()

        except KeyboardInterrupt:
            pass
        finally:
            print('finally')
            self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
            self.device_stop_run()
            avg_size = numpy.mean(self.size)
            size_90 = numpy.percentile(self.size, 90)
            size_10 = numpy.percentile(self.size, 10)
            avg_gfp = numpy.mean(self.fluorescence)
            self.bottom_mir71_threshold = numpy.percentile(self.fluorescence, 10)
            self.upper_mir71_threshold = numpy.percentile(self.fluorescence, 90)
            self.size_threshold = size_90 * DOUBLE_THRES
            self.min_worm_size = size_10 * 0.7

            print('Fluorescence list = \n ' , str(self.fluorescence))
            print('Avg Gfp =' + str(avg_gfp))
            print('90_gfp =' + str(self.upper_mir71_threshold))
            print('10_gfp =' + str(self.bottom_mir71_threshold))

            print('Size list = \n' , str(self.size))

            print('Avg Size =' + str(avg_size))
            print('90_size =' + str(size_90))
            print('10_size =' + str(size_10))

            self.histogram_values.write('Fluorescence list = '+ str(self.fluorescence)
                                        +' based off of '+ str(len(self.fluorescence))+ ' worms'
                                        +'\n Avg GFP = '+ str(avg_gfp)
                                        +'\n 90th percentile GFP = '+ str(self.upper_mir71_threshold)
                                        +'\n 10th percentile GFP = '+ str(self.bottom_mir71_threshold)
                                        +'\n'
                                        +'\n Size list ='+ str(self.size)
                                        +'\n Avg size = '+ str(avg_size)
                                        +'\n 90th percentile size = '+ str(size_90)
                                        +'\n 10th percentile size = '+ str(size_10)
                                        )

            self.histogram_values.close()
            histogram_csv.close()


    def analyze(self, background, worm_image, worm_count, calibration=False):     #TODO: add in autofluorescence            measurements & repeated fluor imaging to test for inter-image variance...

        gfp_fluor_image = self.capture_image(self.cyan)
        gfp_subtracted = abs(gfp_fluor_image.astype('int32')- self.cyan_background.astype('int32'))
        worm_subtracted = abs(worm_image.astype('int32') - background.astype('int32'))
        worm_mask = self.worm_mask(worm_subtracted)

        worm_fluor = self.find_fluor_amount(gfp_subtracted, worm_mask)
        print('GFP value = ' + str(worm_fluor))

        if calibration:
            self.save_image(worm_mask.astype('uint16'), 'calibration_worm_mask', worm_count)
            self.save_image(gfp_fluor_image, 'calibration_worm_fluor', worm_count)
            return worm_fluor

        if not calibration:
            self.save_image(worm_mask.astype('uint16'), 'worm_mask', worm_count)
            self.save_image(gfp_fluor_image, 'fluor_gfp', worm_count)
            self.summary_statistics.write("Gfp Fluorescence: " + str(worm_fluor))

            after_image = self.capture_image(self.bright)
            difference_between_worm_background = (abs(after_image.astype('int32') - background.astype('int32')))
            worm_size_2 = self.size_of_worm(difference_between_worm_background)
            print("Size of worm after imaging :" + str(worm_size_2))

            if worm_size_2 > self.size_threshold:
                worm_count -= 1
                print('Detected Double Worm')
                self.save_image(after_image, 'doubled worm_analyze', worm_count)
                self.summary_statistics.write( '\n doubled worm size of: '+ str(worm_size_2))
                self.device_sort('straight', background, worm_count)
                self.worm_direction = 'straight'
                reason = 'New_worm'
                self.summary_statistics.write("Straight\n")
                print('Doubled worms sorted Straight')
            elif worm_size_2 < self.min_worm_size:
                worm_count -= 1
                print('Detected small worm')
                self.save_image(after_image, 'small worm_analyze', worm_count)
                self.summary_statistics.write( '\n small worm size of: ' + str(worm_size_2))
                self.device_sort('straight', background, worm_count)
                self.worm_direction = 'straight'
                reason = 'lost?'
                self.summary_statistics.write('Straight\n')

            elif worm_fluor > self.upper_mir71_threshold:
                self.up += 1
                self.device_sort('up', background, worm_count)
                self.worm_direction = 'up'
                self.summary_statistics.write("Up\n")
                print('Worm sorted Up    ' + str(self.up))
            elif worm_fluor < self.bottom_mir71_threshold and worm_fluor > 300:
                self.down += 1
                self.device_sort('down', background, worm_count)
                self.worm_direction = 'down'
                self.summary_statistics.write("Down\n")
                print('Worm sorted Down   ' + str(self.down))
            else:
                self.straight += 1
                self.device_sort('straight', background, worm_count)
                self.worm_direction = 'straight'
                self.summary_statistics.write("straight\n")
                print('Worm sorted straight    ' + str(self.straight))

            return worm_fluor, self.worm_direction

    def find_fluor_amount(self, subtracted_image, worm_mask):
        gfp_image = subtracted_image[worm_mask]
        gfp_count = numpy.percentile(gfp_image, 95)
        return gfp_count


class FluorRedGreen(MicroDevice):

    def auto_set_up(self):
        gfp = input('What do you want as the GFP Threshold = ')
        self.gfp_threshold = int(gfp)
        mcherry = input('What do you want as the mcherry Threshold = ')
        self.mcherry_threshold = int(mcherry)
        size = input('What do you want the Size Threshold = ')
        self.size_threshold = int(size)

    def device_set_up(self, num_of_worms):
        """
        Function that runs a given amount of worms through the device and using a humans
        discretion to determine if they are fluorsecent will set threshold values
        for the automation of sorting fluorescent worms
        """
        self.scope.camera.start_image_sequence_acquisition(frame_count=None,
                                                           trigger_mode='Software')
        background = self.capture_image(self.bright)
        self.save_image(background, 'calibration_background')
        self.set_background_areas()
        self.device_start_load()
        gfp_thresh = list()
        mcherry_thresh = list()
        null_gfp_thresh = list()
        null_mcherry_thresh = list()
        size = list()
        i = 0
        cycle_count = 0
        try:
            while True:
                cycle_count += 1
                if not cycle_count % PROGRESS_RATE:
                    print(str(PROGRESS_RATE) + ' cycles')
                if i >= num_of_worms:
                    break
                current_image = self.capture_image(self.bright)
                if self.detect_worms(current_image, background):
                    print('Worm has been detected')
                    #Worm is detected because of a significant change from background.
                    detected_image = current_image
                    #Stop Worms
                    self.device_stop_load()
                    #Position Worms:
                    while True:
                        current_image = self.capture_image(self.bright)
                        if self.lost_worm(current_image, background):
                            print('Worm has been lost')
                            break
                        elif self.positioned_worm(current_image, detected_image):
                            i += 1
                            print('Worm number ' + str(i) + ' out of ' + str(num_of_worms))
                            difference_between_worm_background = (abs(current_image.astype('int32') - background.astype('int32')))
                            self.save_image(current_image, 'calibration_worm'+ str(i))

                            #CYAN
                            current_image = self.capture_image(self.cyan)
                            gfp_image = abs(current_image.astype('int32')- self.cyan_background.astype('int32'))
                            GFP_amount = self.find_fluor_amount(gfp_image)
                            print('GFP amount = ' + str(GFP_amount))
                            is_worm_gfp = input('Is this worm fluorescent in GFP? [y] for yes, [c] for not a worm: ')
                            if is_worm_gfp == 'y':
                                self.save_image(gfp_image.astype('uint16'), 'calibration_gfp_image'+ str(i))
                                gfp_thresh.append(GFP_amount)
                            if is_worm_gfp == 'c':
                                i -= 1
                                self.lamp_off()
                                self.scope.tl.lamp.enabled = True
                                self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
                                self.device_sort('straight', background)
                                time.sleep(2)
                                print('Not a worm')
                                break
                            elif is_worm_gfp != 'y':
                                self.save_image(current_image, 'calibration_null_gfp_worm'+ str(i))
                                null_gfp_image = abs(current_image.astype('int32')- self.cyan_background.astype('int32'))
                                null_gfp_thresh.append(numpy.percentile(null_gfp_image[FLUORESCENT_AREA].astype('int32'),
                                                                        FLUORESCENCE_PERCENTILE))
                            #Green_yellow
                            current_image = self.capture_image(self.green_yellow)
                            mcherry_image = abs(current_image.astype('int32')- self.green_background.astype('int32'))
                            MCHERRY_amount = self.find_fluor_amount(mcherry_image)
                            print('Amount of Mcherry = ' + str(MCHERRY_amount))
                            is_worm_mcherry = input('Is this worm fluorescent in mcherry? [y] for yes ')
                            if is_worm_mcherry == 'y':
                                self.save_image(mcherry_image.astype('uint16'), 'calibration_mcherry_image'+ str(i))
                                mcherry_thresh.append(MCHERRY_amount)
                            elif is_worm_mcherry != 'y':
                                self.save_image(current_image, 'calibration_null_mcherry_worm' + str(i))
                                null_mcherry_image = abs(current_image.astype('int32')- self.green_background.astype('int32'))
                                null_mcherry_thresh.append(numpy.percentile(null_mcherry_image[FLUORESCENT_AREA].astype('int32'),
                                                                            FLUORESCENCE_PERCENTILE))
                            self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
                            self.lamp_off()
                            self.scope.tl.lamp.enabled = True
                            size.append(self.size_of_worm(difference_between_worm_background))
                            self.device_sort('straight', background)
                            self.worm_direction = 'straight'
                            self.clear_worms(background,True, True)
                            break
                        else:
                            detected_image = current_image
                            #9 --> 1

                self.device_start_load()
        except KeyboardInterrupt:
            pass
        finally:
            print('GFP Data: ' + str(gfp_thresh))
            print('Mcherry Data: ' + str(mcherry_thresh))
            print('GFP Null Data: ' + str(null_gfp_thresh))
            print('Mcherry Null Data: ' + str(null_mcherry_thresh))
            print('Worm Size Data: ' + str(size))
            self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME
            self.device_stop_run()
            self.size_threshold = DOUBLE_THRES * numpy.max(size)
            if null_gfp_thresh:
                self.null_gfp_threshold = numpy.mean(null_gfp_thresh)
            if not null_gfp_thresh:
                self.null_gfp_threshold = 0
            if null_mcherry_thresh:
                self.null_mcherry_threshold = numpy.mean(null_mcherry_thresh)
            if not null_mcherry_thresh:
                self.null_mcherry_threshold = 0
            if gfp_thresh:
                self.gfp_threshold = numpy.min(gfp_thresh) - self.null_gfp_threshold
            if not gfp_thresh:
                self.gfp_threshold = NULL_THRESH * self.null_gfp_threshold
            if mcherry_thresh:
                self.mcherry_threshold = numpy.min(mcherry_thresh) -  self.null_mcherry_threshold
            if not mcherry_thresh:
                self.mcherry_threshold = NULL_THRESH * self.null_mcherry_threshold
            self.summary_statistics.write('\nGFP Threshold = '
                                          + str(self.gfp_threshold)
                                          + ' based off of '
                                        + str(len(gfp_thresh))
                                        + ' worms'
                                        +'\n Mcherry Threshold = '
                                        + str(self.mcherry_threshold)
                                        +' based off of '
                                        + str(len(mcherry_thresh))
                                        + ' worms'
                                        +'\n Max worm size seen: ='
                                        + str(self.size_threshold))

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

    def analyze(self, background, current_image):
        """
        function that tells the device what sorting/analzying method to use:
        Is overwritten by a super class
        """
        gfp_fluor_image = self.capture_image(self.cyan)
        gfp_subtracted = abs(gfp_fluor_image.astype('int32')
                             - self.cyan_background.astype('int32'))
        color_value_cyan = self.find_fluor_amount(gfp_subtracted)
        self.save_image(gfp_fluor_image, 'fluor_gfp', True)

        mcherry_fluor_image = self.capture_image(self.green_yellow)
        mcherry_subtracted = abs(mcherry_fluor_image.astype('int32')
                                 - self.green_background.astype('int32'))
        color_value_green = self.find_fluor_amount(mcherry_subtracted)
        self.save_image(mcherry_fluor_image, 'fluor_mcherry', True)

        after_image = self.capture_image(self.bright)
        difference_between_worm_background = (abs(after_image.astype('int32') - background.astype('int32')))
        worm_size = self.size_of_worm(difference_between_worm_background)
        print("Size of worm after imaging :" + str(worm_size))
        Double = False


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
            self.save_image(after_image, 'doubled worm_analyze', True)
            self.summary_statistics.write( '\n doubled worm size of: '+ str(worm_size))
            Double = True

        if ((color_value_cyan > self.gfp_threshold)
        and (color_value_green < self.mcherry_threshold)
        and not Double):
            #Worm is determined to be green and not red.
            self.device_sort('up', background)
            self.worm_direction = 'up'
            self.summary_statistics.write("Up\n")
            print('Worm sorted Up')
        elif ((color_value_cyan < self.gfp_threshold)
        and (color_value_green > self.mcherry_threshold)
        and not Double):
            #Worm is detremined to be red and not green.
            self.device_sort('down', background)
            self.worm_direction = 'down'
            self.summary_statistics.write("Down\n")
            print('Worm sorted Down')
        else:
            self.device_sort('straight', background)
            self.worm_direction = 'straight'
            self.summary_statistics.write("Straight\n")
            print('Worm sorted Straight')
