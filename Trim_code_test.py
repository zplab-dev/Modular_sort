
#Merging the code, trying to slim it down

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

#Textbetl key = 08a7e7d3335d92542dfec857461cfb14af5e0805HQINVtRpVqvDjo9O2wa2I6tTo

FLUOR_PIXEL_BRIGHT_VALUE = 500		#Used for finding brightness, newer/better way to do this?
DETECTION_THRES = 4
POSITION_THRES = 3
CLEARING_THRES = 3
LOST_CUTOFF = 1.5
DOUBLE_THRESH = 1.3

CYAN_EXPOSURE_TIME = 50  #7 for lin-4, 50 for mir-71?
YELLOW_EXPOSURE_TIME = 8	#Yellow_green for red fluor?
BRIGHT_FIELD_EXPOSURE_TIME = 3

PICTURE_DELAY = .01
SORTING_INTERVAL = .4

BACKGROUND_REFRESH_RATE = 100000	#Not actually currently getting called, could use this or reset by number of worms
PROGRESS_RATE = 100

MIN_GFP_THRESH = 400 #Will need to reset to account for brighter exposure


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
        self.date = exp_direct.split('/')[-1]
        self.summary_location = self.file_location.joinpath('summary_' + self.date + '.txt')
        self.summary_statistics = open(str(self.summary_location),'w')
        summary_csv_location = self.file_location.joinpath('summary_csv.txt')   #TODO: change to .csv, fix path/writing
        summary_csv = open(str(summary_csv_location), 'w')
        header = ['worm_number', 'size', 'fluorescence', 'time', 'direction', 'reason']
        summary_csv.write(','.join(header) + '\n')
        
        #self.data_location = self.file_location.joinpath('wormdata.csv')   #old data saving, never implemented

        self.device_clear_tubes()

        worm_count = 0
        self.up = 0
        self.down = 0
        self.straight = 0
        self.TIME_FOR_PUSHING = .8
        #worm_data = list()
        
        """					TODO: Fix threading/pausing to make more useful
        #Pausing stuff
        self.running = True
        super().__init__(daemon=True)
        self.quitting = False
        self.cleared = False
        
        self.resume() #(testing)
		"""
		
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
        
        

