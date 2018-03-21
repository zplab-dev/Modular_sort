
#Merging the code, trying to slim it down
#TODO: Split into more modular files

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

CYAN_EXPOSURE_TIME = 50  #7 for lin-4, 50 for mir-71?
YELLOW_EXPOSURE_TIME = 8	#for mcherry
RED_EXPOSURE_TIME = 50 #For autofluorescence
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
        self.scope.camera.exposure_time = YELLOW_EXPOSURE_TIME
        time.sleep(PICTURE_DELAY)

    def red(self):
        self.scope.tl.lamp.enabled = False
        self.scope.il.spectra.cyan.enabled = False
        self.scope.il.spectra.green_yellow.enabled = False
        self.scope.il.spectra.red.enabled = True
        self.scope.camera.exposure_time = RED_EXPOSURE_TIME
        time.sleep(PICTURE_DELAY)
            
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

    def set_bf_background(self):        #Why does this take two images while fluor bgs only take one?
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

    def set_green_background(self, worm_count):
        self.green_yellow_background = self.capture_image(self.green_yellow)
        self.save_image(self.green_yellow_background, 'green_yellow_background', worm_count)

    def set_red_background(self, worm_count):
        self.red_background = self.capture_image(self.red)
        self.save_image(self.red_background, 'red_background', worm_count)

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
        mask = worm_mask(subtracted_image)
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

    def check_cleared(self, background, worm_count):        #TODO: Fix force reset option, because flag never triggers
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

    def run(self):      #TODO: break this into general sorting function, put other relevant details in run func.
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

        self.setup_csv(self.file_location, self.info)

        self.scope.camera.start_image_sequence_acquisition(frame_count=None, trigger_mode='Software')

        worm_count = 0
        self.up = 0
        self.down = 0
        self.straight = 0
        cycle_count = 0
        self.boiler = boiler()
        self.reset = False
        self.message_sent = False
        
        #0 Setting Background
        self.set_background_areas(worm_count)
        print('setting backgrounds')
        time.sleep(1)

        #1 Loading Worms
        self.device_start_load()
        time_start = time.time()
        time_seen = time_start #Initial setting for time_seen before worm found

        #2 Detect Worms
        try:
            while True:

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
                    time_seen = time.time()
                    detected_image = current_image

                    #3 Stop worms
                    self.device_stop_load()

                    #4 Position worms
                    while True:
                        current_image = self.capture_image(self.bright)
                        if self.lost_worm(current_image, background):       #Is this actually working
                            print('Worm was lost')
                            self.device.execute(SEWER_CHANNEL_PRESSURE)
                            time.sleep(.1)
                            break
                        elif self.positioned_worm(current_image, detected_image):   #Does this do it's job? Sometimes rejects big worms
                            #Maybe include some param that says if a worm is too close to the edge, let it get to the center 

    """
    Some pseudocode:
    def run(self):
        build_hist()    #class-specific, calls sort and analyze for worms 1-100
            #ask if build hist, check for variables, or input variables
        sort()          #Still using class-specific analysis function
            #sort should only look for worms, then ask analyze what to do
        save_data()     #Should also be class-specific, as params saved change from case to case
            #Does this mean that if sort exits badly data can't be saved? Does data have to be saved continuously?
    """




class GFP(MicroDevice):

    def setup_csv(self, file_location, info):
        self.summary_csv_location = file_location.joinpath('summary_csv' + info + '.csv')
        self.summary_csv = open(str(self.summary_csv_location), 'w')
        header = ['worm_number', 'size', 'fluorescence', 'time', 'direction']
        self.summary_csv.write(','.join(header) + '\n')

    def set_background_areas(self, worm_count):
        self.set_bf_background(worm_count)
        self.set_cyan_background(worm_count)
        self.lamp_off()
        self.scope.tl.lamp.enabled = True
        self.scope.camera.exposure_time = BRIGHT_FIELD_EXPOSURE_TIME

    def manual_set_up(self):
        self.upper_mir71_threshold = int(input('Upper threshold = '))
        self.bottom_mir71_threshold = int(input('Lower threshold = '))
        self.size_threshold = int(input('Max size threshold = '))
        self.min_worm_size = int(input('Small size threshold = '))

