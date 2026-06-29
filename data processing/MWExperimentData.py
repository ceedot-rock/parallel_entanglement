#external imports
import numpy as np #type: ignore
import os
import matplotlib.pyplot as plt #type: ignore
import matplotlib as mpl #type: ignore

#Internal imports
from savedata import SaveData

class MWExperimentData(SaveData):

    def __init__(self, filename, check_line_ID=False, normalize_counts=False):

        #==================================================================================================
        # (LEGACY) DAQ Counter Assumptions
        # All counters used should be added to the columnNames of the data object as 'ctrN' where N=0,1,2,3
        # 'ctr0' : open counter.
        # 'ctr1' : DAQ clock, can not be used as counter currently.
        # 'ctr2' : APD counts containing spin dependent information.
        # 'ctr3' : If used, Line ID. if not used, no issues, but use is highly suggested.
        # 
        # This is legacy before DAQ counters were configured and saved in the sequence generation process.
        # Other allocations of counters may be different from above if set in the seq gen config file.
        #==================================================================================================

        if ".mat" not in filename:
            filename += ".mat"

        full_filename = filename

        file_exists = os.path.isfile(full_filename)
        if not file_exists:
            raise ValueError(f"file '{full_filename}' not found.")

        #Initialize parent SaveData class to process .mat file
        super().__init__(full_filename)

        #Create usual data sets based on matlab file
        self.parse_data() 
        self.reshape_data()
        self.rescale_counts(normalize_counts)
        if check_line_ID:
            self.check_line_ID()

        #since rescaling the data can change what we mean by "norm counts", cast the main data sets by generic names 
        self.x_data = self.independent_var
        self.y_data = self.norm_counts
        self.y_data_err = self.norm_err

        return
    
#================================
#Functions for handeling raw data
#================================

    def parse_data(self):
        # parse out data saved that will be used in the reshape process

        #Cast as numpy array to use np.where
        self.sDataHeader = np.array(self.sDataHeader)

        #Get nVars
        nVars_idx = np.where(self.sDataHeader == 'nVars')[0]
        if len(nVars_idx) == 0:
            raise ValueError("Could not locate 'nVars' in sDataHeader")
        elif len(nVars_idx) > 1:
            raise ValueError("Multiple entries of 'nVars' in sDataHeader")
        else:
            if self.nSweeps > 1:
                self.nVars = self.sData[0][nVars_idx[0]] #index sData to first entry, all are same
            else:
                self.nVars = self.sData[nVars_idx[0]]

        #Get total repeats aka samples per sweep (?)
        total_repeats_idx = np.where(self.sDataHeader == 'samples per sweep')
        if len(total_repeats_idx) == 0:
            raise ValueError("Could not locate 'samples per sweep' in sDataHeader")
        elif len(total_repeats_idx) > 1:
            raise ValueError("Multiple entries of 'samples per sweep' in sDataHeader")
        else:
            if self.nSweeps > 1:
                self.total_repeats = self.sData[-1][total_repeats_idx[0]] #take last of sData since it has the true total repeats
            else:
                self.total_repeats = self.sData[total_repeats_idx[0]]

        #Get which counters were used
        if len(self.columnNames) == 0:
            raise ValueError("No column names were found. These are expected to be the counters used on the DAQ")
        
        valid_counters = ['ctr0','ctr1','ctr2','ctr3']
        self.num_counters = len(valid_counters) #avoiding hard coding the number 4

        for name in self.columnNames:
            if name not in valid_counters:
                raise ValueError("Invalid name in column names, expecting 'ctr0','ctr1','ctr2' or 'ctr3'")
            
        self.ctr_idxs = [] #doesn't assume columns are ordered but creates an ordered list of indices for ctr0,ctr1,...
        for ctr in valid_counters:
            where_ctr = np.where(np.array(self.columnNames) == ctr)[0]
            if len(where_ctr) == 0:
                self.ctr_idxs.append(None)
            elif len(where_ctr) > 1:
                raise ValueError("more than 1 column named "+ctr+". Each ctr should appear at most once")
            else:
                self.ctr_idxs.append(where_ctr[0])

        #if self.ctr_idxs[2] == None:
        #    raise ValueError("ctr2 needs to be taking data as it is assumed this is spin dependent PL")
        
        self.num_counters_used = 0
        for ctr_idx in self.ctr_idxs:
            if ctr_idx != None:
                self.num_counters_used += 1

        #Get independent variable information
        if 'var' in self.extraInfo.keys():
            if 'values' in self.extraInfo['var'].keys():
                self.independent_var = self.extraInfo['var']['values']
                if isinstance(self.independent_var,int):
                    self.independent_var = np.array([self.independent_var])
            else:
                raise ValueError("Could not locate 'values' in extraInfo['vals'] keys")
            if 'name' in self.extraInfo['var'].keys():
                self.independent_var_name = self.extraInfo['var']['name']
            else:
                raise ValueError("Could not locate 'name' in extraInfo['vals'] keys")
        else:
            raise ValueError("Could not locate 'var' in extraInfo keys")
        
        #Get experiment name
        if 'experiment_name' in self.extraInfo.keys():
            self.experiment_name = str(self.extraInfo['experiment_name'])
        else:
            self.experiment_name = ''
            print("WARNING: 'experiment_name' could not be located in extraInfo, setting to ''")
    
        #Get readout calibration flag and change count data accordingly
        if 'readout_calibration' in self.extraInfo.keys():
            self.readout_calib = int(self.extraInfo['readout_calibration'])
        else:
            self.readout_calib = False
            print("WARNING: 'readout_calibration' could not be located in extraInfo, setting to False")

        #Get type of readout
        if 'readout' in self.extraInfo.keys():
            self.readout = self.extraInfo['readout']
        else:
            self.readout = None
            print("WARNING: 'readout' could not be located in extraInfo, setting to None")

        #Get SRS center frequency used
        if 'center_frequency' in self.extraInfo.keys():
            self.center_freq = self.extraInfo['center_frequency']
        else:
            self.center_freq = None
            print("WARNING SRS center frequency was not saved in dataObj")

        #Get MW attenuation used as set by the py_seq_gen_acquire script
        if 'attenuation' in self.extraInfo.keys():
            self.attenuation = self.extraInfo['attenuation']
        else:
            self.attenuation = None
            print("WARNING attenuation was not saved in dataObj")
    
        #Get data partition based on nuclear readout lines interelaeaved
        if 'data_partition' in self.extraInfo.keys():
            self.data_partition = self.extraInfo['data_partition']
        else:
            self.data_partition = None

        #Get magnetic field information used in sequence generation
        if 'B_mag_tuple' in self.extraInfo.keys() and 'B_theta_tuple' in self.extraInfo.keys():
            self.B_mag_tuple = self.extraInfo['B_mag_tuple'] #ordering is [val, upper error, lower error]
            self.B_theta_tuple = self.extraInfo['B_theta_tuple'] #ordering is [val, upper error, lower error]
        else:
            self.B_mag_tuple = None
            self.B_theta_tuple = None

        #Get DAQ counter numbers connection from config
        if 'APD_DAQ_ctr' in self.extraInfo.keys():
            self.APD_DAQ_ctr = self.extraInfo['APD_DAQ_ctr']
        else:
            self.APD_DAQ_ctr = 2 #for legacy compatibility
        if 'line_ID_ctr' in self.extraInfo.keys():
            self.line_ID_ctr = self.extraInfo['line_ID_ctr']
        else:
            # self.line_ID_ctr = 3 #for legacy compatibility
            self.line_ID_ctr = None
        if 'ref_ctr' in self.extraInfo.keys():
            self.ref_ctr = self.extraInfo['ref_ctr']
        else:
            self.ref_ctr = None # haven't consistently used this in the past

        #compare against counters found in the data header
        #if self.ctr_idxs[self.APD_DAQ_ctr] is None:
        #    raise ValueError((f"The APD DAQ ctr was set to {self.APD_DAQ_ctr} in the config (python end), but this counter wasn't used to take data (matlab end).\n"
        #                       "Either update the config for the counters used, or add new counters to the DAQ instance on matlab."))
        #if self.ctr_idxs[self.line_ID_ctr] is None:
        #    raise ValueError((f"The line ID DAQ ctr was set to {self.line_ID_ctr} in the config (python end), but this counter wasn't used to take data (matlab end).\n"
        #                       "Either update the config for the counters used, or add new counters to the DAQ instance on matlab."))

        return
    
    def reshape_data(self):

        #Process photon count data
        repeats_per_line = self.total_repeats / self.nVars
        if isinstance(self.data, np.ndarray):
            raw_counts = np.squeeze(self.data)
        else:
            raw_counts = np.array(self.data)
        
        # Calculate the difference in counts

        # Single sweep experiment
        if self.nSweeps == 1:

            diff_counts = np.diff(raw_counts,axis = 0)
            
            # 1 counter used case isn't technically 2 dimensional, so reshape
            diff_counts = diff_counts.reshape(-1,int(self.num_counters_used))
            # append 0 to account for taking diff
            diff_counts = np.append(diff_counts, np.zeros((1,int(self.num_counters_used))),axis = 0)
        
        # Multi-sweep experiment
        elif self.nSweeps > 1:
            
            #NOTE IMPORTANT between sweeps, the diff is no longer defined. Only within one sweep does the relative counts on the DAQ carry meaningful information.
            
            #single counter
            if len(np.shape(raw_counts)) == 2:

                if self.num_counters_used != 1:
                    raise ValueError("Expecting only a single countered based on shape of the data.")
                
                filtered_rows = [row[~np.isnan(row)] for row in raw_counts]

                diff_counts = []
                for row in filtered_rows:
                    row_diff = np.diff(row)
                    row_diff = np.append(row_diff,0)
                    diff_counts.extend(list(row_diff))
                diff_counts = np.array(diff_counts)
                diff_counts = diff_counts.reshape(-1,1)

            #multi counter
            elif len(np.shape(raw_counts)) == 3: 

                if self.num_counters_used <= 1:
                    raise ValueError("Expecting multiple counters based on the shape of the data.")

                diff_count_list = []
                for i_raw_counts in raw_counts:

                    i_diff = np.diff(i_raw_counts,axis=0)

                    #each sweep can be a different "true" length, remove the place holding nan values generated by saveData
                    filtered_diff = i_diff[~np.isnan(i_diff)] #returns a 1 dimensional array
                    filtered_diff = np.reshape(filtered_diff,(len(filtered_diff)//self.num_counters_used,self.num_counters_used))
                    # append 0 to account for taking diff
                    filtered_diff = np.append(filtered_diff, np.zeros((1,int(self.num_counters_used))),axis = 0)
                    diff_count_list.append(filtered_diff)

                diff_counts = np.concatenate(diff_count_list,axis=0)

            else:
                raise ValueError("Unexpected data shape.")

        else:
            raise ValueError("'nSweeps' should be a positive integer. Check where this is being written to saveData obj.")
        
        #initialize shaped empty arrays to place data into
        self.norm_counts = np.empty((self.num_counters,int(self.nVars)))
        self.norm_err = np.empty_like(self.norm_counts)
        self.raw_counts_per_var = np.empty((self.num_counters,int(self.nVars),int(repeats_per_line)))

        for iCol,ctr_idx in enumerate(self.ctr_idxs):
            if ctr_idx != None:
                diff_count_col = diff_counts[:,ctr_idx]
                reshaped_counts_matrix = np.reshape(diff_count_col,(int(repeats_per_line),int(self.nVars)))
                self.raw_counts_per_var[iCol] = np.transpose(reshaped_counts_matrix)
                norm_counts_row = np.sum(reshaped_counts_matrix, axis = 0) / repeats_per_line
                self.norm_counts[iCol] = norm_counts_row
                norm_err_row = np.std(reshaped_counts_matrix, axis = 0) / np.sqrt(repeats_per_line)
                self.norm_err[iCol] = norm_err_row
            else: #this counter was not used, set to 0 counts
                self.norm_counts[iCol] = np.zeros(int(self.nVars)).astype(int)
                self.norm_err[iCol] = np.zeros(int(self.nVars)).astype(int)

        #If ctr 3 was used, it was for Line ID. Seperate out that data before truncating / rescaling when calibrating
        # if self.ctr_idxs[self.line_ID_ctr] != None:
        #     self.line_ID_data = self.norm_counts[self.line_ID_ctr]
    
        return
    
    def check_line_ID(self,plot=False):
        # If ctr3 wasn't used, Line ID can't be checked or plotted
        if self.ctr_idxs[self.line_ID_ctr] == None:
            if plot:
                print(f"Line ID cannot be plotted because the counter corresponding to (counter {self.line_ID_ctr}) it was not found.")
                return
            else:
                print(f"Line ID cannot be check because the counter corresponding to (counter {self.line_ID_ctr}) it was not found.")
                return
        
        # Check Line ID

        zero_diff_count = 0 #slope zero points can happen, but only a certain number of times (see below for case)
        for diff in (np.diff(self.line_ID_data)).tolist():
            if round(diff,2) == 0:
                zero_diff_count += 1
            if round(diff,1) != 1 and round(diff,1) != 0: #this needs to be rounded due to the off by 1 error on the very last point of the last shot
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(self.line_ID_data,"-o",label="line ID data")
                ax.set_xlabel("Independent variable")
                ax.set_ylabel("Line ID number")
                ax.set_title("This should be a slope 1 line")
                ax.legend()
                raise ValueError("Line ID check failed.\nThis could be due to a mismatch between AWG lines and expected data or due to a data reshaping issue.\nIf plot above looks okay, it could be a round issue due to the off by one error on the very last data point.")

        #readout calib lines always get 1 count on Line ID
        #we assume there are only 2 readout lines per experiment, so should be at most 2 slope zero points in Line ID
        if self.ctr_idxs[self.line_ID_ctr] != None:
            if self.readout_calib:
                if zero_diff_count != 2:
                    raise ValueError("Readout calibration and ctr3 were used, but "+str(zero_diff_count)+" slope zero points were detected in Line ID.")
            else:
                if zero_diff_count > 0:
                    raise ValueError("Readout calibration was not used, but ctr3 was. There should only be slope 1 data in Line ID, but detected "+str(zero_diff_count)+" slope zero points.")

        print("Line ID check passed")

        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.line_ID_data,"-o",label="ctr3 / Line ID Data")
            ax.set_xlabel("Line Number")
            ax.set_ylabel("Average Counts")
            ax.legend()

        return
    
    def rescale_counts(self, normalize_counts):

        # NOTE IMPORTANT!!!
        #if readout is calibrated, the first 2 data points *should* be the high and low PL levels
        #to be sure, check your sequence. First MW line should be empty and the second MW line should be an AP IQ pulse / pi pulse
    
        if self.readout_calib and not normalize_counts:

            self.e0_PL_level = self.norm_counts[self.APD_DAQ_ctr][0] #|0>_e
            self.e1_PL_level = self.norm_counts[self.APD_DAQ_ctr][1] #|1>_e

            calibrated_norm_counts = np.empty((self.num_counters,int(self.nVars - 2)))
            calibrated_norm_err = np.empty_like(calibrated_norm_counts)

            for i,ctr_idx in enumerate(self.ctr_idxs):
                if ctr_idx == None or i == self.line_ID_ctr:
                    #if ctr wasn't used, or it it's the Line ID counter, don't calibrate (WILL fuck up Line ID)
                    calibrated_norm_counts[i] = self.norm_counts[i][2:] 
                    calibrated_norm_err[i] = self.norm_err[i][2:]
                else:
                    if self.readout in ["PL", "SCC"]:
                        calibrated_norm_counts[i] = (self.norm_counts[i][2:] - self.e1_PL_level) / (self.e0_PL_level - self.e1_PL_level)
                        calibrated_norm_err[i] = self.norm_err[i][2:] / np.abs(self.e0_PL_level - self.e1_PL_level)
                    else:
                        raise ValueError("unexpected readout method. See options above. If 'readout'=None, this is an error on the sequence side.")
        
            self.norm_counts = calibrated_norm_counts
            self.norm_err = calibrated_norm_err

            # create y-label for plotting functions
            self.readout_label = r"P($|0\rangle_e$)"

        elif self.readout_calib and normalize_counts:
            raise ValueError("Readout was calibrated, there's no reason to normalize counts also.")
        
        elif not self.readout_calib and normalize_counts:

            for i,ctr_idx in enumerate(self.ctr_idxs):
                max_count = np.max(self.norm_counts[i])
                if max_count != 0: #don't update unused counters
                    self.norm_counts[i] *= (1/max_count)
                    self.norm_err[i] *= (1/max_count)
            self.readout_label = "Normalized counts"

        else:
            # create y-label for plotting functions
            self.readout_label = "Average counts"

        return
    
#===========================
#Functions for plotting data
#===========================

    def plot_data(self,ax,ctrs_to_plot=None,**kwargs):

        #check which counters to plot
        if ctrs_to_plot is None:
            ctrs_to_plot = [self.APD_DAQ_ctr]

        #kwargs for this function call
        x_scale = kwargs.pop("x_scale",1)
        x_shift = kwargs.pop("x_shift",0)
        y_scale = kwargs.pop("y_scale",1)
        exp_val = kwargs.pop("exp_val",0)
        self.abs_freq = kwargs.pop("abs_freq",False)
        colors = kwargs.pop("colors",[])

        # default plot call kwargs if any of these aren't passed specifically
        label = kwargs.setdefault("label",self.experiment_name+" data")
        color = kwargs.setdefault("color","dimgrey")
        kwargs["color"] = color
        kwargs.setdefault("fmt",".")
        kwargs.setdefault("linestyle","None")
        kwargs.setdefault("capsize", mpl.rcParams.get("errorbar.capsize", 0))
        kwargs.setdefault("markersize", mpl.rcParams.get("lines.markersize", 4))

        # X axis manipulation
        x_label = self.independent_var_name
        if self.abs_freq:
            self.x_data = self.center_freq*np.ones(len(self.x_data)) + self.x_data
            x_label = "MW Frequency [MHz]"
        self.x_data = self.x_data * x_scale
        self.x_data = self.x_data - np.ones_like(self.x_data)*x_shift
        
        # Y axis manipulation
        y_label = self.readout_label
        if exp_val:
            # don't transform the Line ID data
            self.y_data[:3] = (2 * self.y_data[:3] - 1)
            self.y_data_err[:3] = 2 * self.y_data_err[:3]
            y_label = r"$\langle Z \rangle$"

        if colors == []:
            for ctr_i in ctrs_to_plot:
                ebar = ax.errorbar(self.x_data,y_scale*self.y_data[ctr_i],self.y_data_err[ctr_i],**kwargs)
        else:
            if len(colors) != len(self.x_data):
                raise ValueError("'colors' list should be the length of the data, in this case length: "+str(len(self.x_data)))
            else:
                for ctr_i in ctrs_to_plot:
                    for i,color in enumerate(colors):
                        kwargs["color"] = self.color_parser(color)
                        kwargs["label"] = ""
                        ebar = ax.errorbar(self.x_data[i],y_scale*self.y_data[ctr_i][i],self.y_data_err[ctr_i][i],**kwargs)


        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        if label != "":
            ax.legend(loc="best")
        return ebar   

#=============
#Magic Methods
#=============
    
    def __iadd__(self,other):
        self.independent_var = np.append(self.independent_var,other.independent_var)
        self.norm_counts = np.append(self.norm_counts,other.norm_counts,axis=1)
        self.norm_err = np.append(self.norm_err,other.norm_err,axis=1)
        self.x_data = np.append(self.x_data,other.x_data)
        self.y_data = np.append(self.y_data,other.y_data,axis=1)
        self.y_data_err = np.append(self.y_data_err,other.y_data_err,axis=1)
        return self