# -*- coding: utf-8 -*-
"""
Created on Dec 12 2022
@author:  BIONICS_LAB
@company: DZNE
"""

import sys
sys.path.append('../')
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter
import h5py
import pandas as pd
import LFP_Denoising as LFP_denoising
import Spikes_Denoising as Spikes_Denoising
import numpy as np
import os

class DENOISING:
    """
    The main class for denoising, including LFPs and Spikes, and saves the denoising results.
    """
    def __init__(self, srcfilepath):
        self.srcfilepath = srcfilepath  # main path

    def filter_data(self, data, SamplingRate, low, high, order=2):
        '''
        Filter the data by the bandpass filter. The bandpass filter is used to filter the LFPs/Spikes data.
         Parameters
            ----------
            :param data: The original data.
            :param SamplingRate: The sampling rate of the data.
            :param low: The low frequency of the bandpass filter.
            :param high: The high frequency of the bandpass filter.
            :param order: The order of the bandpass filter.
        Returns
        ----------------
            :return: The filtered signal.
        '''
        # Determine Nyquist frequency
        nyq = SamplingRate / 2
        # Set bands
        low = low / nyq
        high = high / nyq
        # Calculate coefficients
        b, a = butter(order, [low, high], btype='band')
        # Filter signal
        filtered_data = lfilter(b, a, data)
        return filtered_data

    def getRawdata(self, start, stop, expFile, channel, thr=2000):
        '''
        Get the raw data from the brw file. The raw data is the original data without any processing.
         Parameters
            ----------
            :param start: float
                    The start time of the raw data.
            :param stop: float
                    The stop time of the raw data.
            :param expFile: string
                    The file name of the .brw file.
            :param channel: int
                    The channel of the raw data.
            :param thr: float
                    The threshold of the raw data.
         Returns
            ----------------
            :return: raw: numpy array
                The waveform generated from the .brw file.
        '''
        filehdf5 = h5py.File(expFile, 'r')
        # CONSTANT
        MaxValue = 4096.
        MaxVolt = np.asarray(filehdf5["3BRecInfo"]["3BRecVars"]["MaxVolt"])[0]
        MinVolt = np.asarray(filehdf5["3BRecInfo"]["3BRecVars"]["MinVolt"])[0]
        NRecFrames = np.asarray(filehdf5["3BRecInfo"]["3BRecVars"]["NRecFrames"])[0]
        SignInversion = np.asarray(filehdf5["3BRecInfo"]["3BRecVars"]["SignalInversion"])[0]
        stepVolt = (MaxVolt - MinVolt) / MaxValue
        version = int(filehdf5["3BData"].attrs['Version'])
        rawData = filehdf5["3BData"]["Raw"]
        if start < 0:
            start = 0
        if stop >= NRecFrames:
            stop = NRecFrames - 1

        if isinstance(channel, int) or isinstance(channel, float):
            # get one Single channel
            if version == 100:
                raw = ((rawData[int(start):int(stop), channel] - (4096.0 / 2)) * stepVolt * SignInversion)
            else:
                raw = rawData[int(start) * 4096:int(stop) * 4096]
                raw = raw.reshape((len(raw) // 4096, 4096))
                raw = (raw[:, channel] - (4096.0 / 2)) * stepVolt * SignInversion

        elif isinstance(channel, str):
            # Get all channels
            if version == 100:
                raw = ((rawData[int(start):int(stop), :] - (4096.0 / 2))) * SignInversion
            else:
                raw = rawData[int(start) * 4096:int(stop) * 4096]
                raw = raw.reshape((len(raw) // 4096, 4096))
                raw = (raw - (4096.0 / 2)) * SignInversion
            #           Put to 0 saturation sample
            index = np.where(raw > thr)
            raw[index[0], index[1]] = 0
            raw = raw * float(stepVolt)
        elif isinstance(channel, np.ndarray):
            # Get an array of channels
            if version == 100:
                raw = ((rawData[int(start):int(stop), channel] - (
                        4096.0 / 2))) * stepVolt * SignInversion
            else:
                raw = rawData[int(start) * 4096:int(stop) * 4096]
                raw = raw.reshape((len(raw) // 4096, 4096))
                raw = (raw[:, channel] - (4096.0 / 2)) * stepVolt * SignInversion
        else:
            print("Error: channel must be an int, a float, a str or a numpy array")

        return raw


    def get_filename_path(self,filepath, filetype):
        '''
        Get the filename and the path of the file, which is the same type.
        Parameters
            ----------
        :param filepath: string
                The path of the file.
        :param filetype: string
                The type of the file.
         Returns
        ----------------
            :return: filename: list
                    The list of the filename.
                    Root: list
                    The list of the path of the file.
        '''
        filename = []
        Root = []
        for root, dirs, files in os.walk(filepath):
            for i in files:
                if filetype in i:
                    filename.append(i)
                    Root.append(root)
        return filename, Root

    def DENOISING(self, Analysis_Item='Spikes'):
        '''
        The main function for denoising, including LFPs and Spikes, and save the denoising results.
        Parameters
        ----------
            :param Analysis_Item:
                    'Spikes' or 'LFPs'
        Returns
        ----------------
        '''
        filetype_bxr = '.bxr'
        filename_bxr, Root = self.get_filename_path(self.srcfilepath, filetype_bxr)

        for expFile in filename_bxr:
            if expFile[0] != '.':
                print(expFile)
                if Analysis_Item == 'LFPs':
                    Analysis = LFP_denoising.LFP_DENOISING_Function(self.srcfilepath)
                    lfpChId, lfpTimes, LfpForms = Analysis.AnalyzeExp(expFile=expFile)
                else:
                    Analysis = Spikes_Denoising.Spikes_DENOISING_Function(self.srcfilepath)
                    SpikesChId, SpikesTimes = Analysis.AnalyzeExp(expFile=expFile)
    def Raster_plot(self,Analysis_Item='Spikes',denoising=False):
        '''
        The main function for raster plot, including LFPs and Spikes, and save the raster plot results.
        :param Analysis_Item: String
                'Spikes' or 'LFPs'
        :param denoising: Boolean
                True: the raster plot for DENOISING; False: the raster plot for Raw.
        :return:
        '''
        filetype_bxr = '.bxr'
        filename_bxr, Root = self.get_filename_path(self.srcfilepath, filetype_bxr)

        for expFile in filename_bxr:
            if expFile[0] != '.':
                print(expFile)
                # Read the LFPs/Spikes data
                filehdf5_bxr = h5py.File(self.srcfilepath + expFile, 'r')  # read LFPs bxr files
                NRecFrames = np.asarray(filehdf5_bxr["3BRecInfo"]["3BRecVars"]["NRecFrames"])[0]
                samplingRate = np.asarray(filehdf5_bxr["3BRecInfo"]["3BRecVars"]["SamplingRate"])[0]
                recordingLength = NRecFrames / samplingRate  # recording duraton in [s]

                if Analysis_Item == "LFPs":
                    LfpForms = np.asarray(filehdf5_bxr["3BResults"]["3BChEvents"]["LfpForms"])
                    lfpChId = np.asarray(filehdf5_bxr["3BResults"]["3BChEvents"]["LfpChIDs"])
                    lfpTimes = np.asarray(filehdf5_bxr["3BResults"]["3BChEvents"]["LfpTimes"])/ samplingRate
                    binsDistr = np.arange(0, np.ceil(max(lfpTimes)), 0.1)  # fixed bin size
                    if denoising:
                        lfpTimes = np.load(self.srcfilepath + expFile[:-4] + '_denoised_LfpTimes' + '.npy')
                        lfpChId = np.load(self.srcfilepath + expFile[:-4] + '_denoised_LfpChIDs' + '.npy')

                else:
                    lfpTimes = np.asarray(filehdf5_bxr["3BResults"]["3BChEvents"]["SpikeTimes"])/ samplingRate
                    lfpChId = np.asarray(filehdf5_bxr["3BResults"]["3BChEvents"]["SpikeChIDs"])
                    binsDistr = np.arange(0, np.ceil(max(lfpTimes)), 0.1)  # fixed bin size
                    if denoising:
                        lfpTimes = np.load(self.srcfilepath + expFile[:-4] + '_denoised_SpikesTimes' + '.npy')
                        lfpChId = np.load(self.srcfilepath + expFile[:-4] +'_denoised_SpikesChIDs' + '.npy')


                lfpTimes_last_time = lfpTimes
                averageEventRateHist, averageEventRateBinsEdge = np.histogram(lfpTimes_last_time, bins=binsDistr,normed=False, weights=None, density=None)
                averageEventRateXScale = (averageEventRateBinsEdge[1:] - averageEventRateBinsEdge[0:-1]) / 2.0 + averageEventRateBinsEdge[0:-1]
                # Plot the raster plot

                fig, ax = plt.subplots(2, 1, figsize=(15, 10))

                ax[0].plot(lfpTimes_last_time, lfpChId, linestyle='', marker='|', c='black', markersize=2,
                           markeredgewidth=1)

                ax[0].set_ylabel('')

                ax[0].set_xlim((averageEventRateXScale[0], averageEventRateXScale[-1]))

                for tick in ax[0].yaxis.get_major_ticks():
                    tick.label.set_fontsize(5)

                ax[1].plot(averageEventRateXScale, averageEventRateHist, color='black')

                ax[1].set_ylabel('# events')
                ax[1].set_xlabel('Time[Sec]')

                ax[1].set_xlim((averageEventRateXScale[0], averageEventRateXScale[-1]))
                fig.tight_layout()
                if denoising:
                    colorMapTitle = expFile[:-4] + '_Raster_plot_DENOISING_' + Analysis_Item
                else:
                    colorMapTitle = expFile[:-4] + '_Raster_plot_Raw_' + Analysis_Item
                fig.savefig(self.srcfilepath + colorMapTitle + ".png", format='png', dpi=600)
                plt.close()

if __name__ == '__main__':
    srcfilepath = r'H:/DENOISING/Spike/'
    analysis = DENOISING(srcfilepath)
    analysis.DENOISING(Analysis_Item='Spikes')
    analysis.Raster_plot(Analysis_Item='Spikes',denoising=True)
