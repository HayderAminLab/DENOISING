# **DENOISING: Dynamic Enhancement and Noise Overcoming in Multimodal Neural Observations via High-density CMOS-based Biosensors**

<p align="center">
  <img width="1090" height="614" src="https://github.com/HayderAminLab/DENOISING/assets/158823360/1c5ad1e1-5746-44b4-8a2c-92ceecde429b">
</p>

 A computational framework designed to enhance the fidelity and clarity of large-scale bioelectrical signals by dynamically mitigating noise. It addresses the challenge of 
 pervasive noise in neural recordings from high-density microelectrode arrays, offering an adaptable solution for various noise types and signal dimensions. By improving the signal-to-noise ratio 
 and facilitating the separation of signal and noise components, DENOISING aids in the accurate analysis and interpretation 
 of complex neural connectivity and dynamics without relying on specific data transformations.

## **Contents**

  - README - current file
  - LICENSE - usage and redistrubution conditions
  - DENOISING.py - main top-level function used to denoise multimodal datasets
  - LFP_Denoising.py - top-level function used for denoising detected-LFP events
  - Spikes_Denoising.py - top-level function used for denoising detected-spikes

## **Getting Started**

  - Download data from Zenodo and separate into subfolders based on datatype i.e. LFP and Spike.
  - If using personal data, format folder structure as instructed for the provided data. 
  - Clone the respository locally.
       > git clone https://github.com/HayderAminLab/DENOISING.git
  - Go to DENOISING.py and work through it as follows:
    1. Change srcfilepath to folder with dataset
    2. Indicate data type (Analysis_Item='LFPs') OR (Analysis_Item='Spikes')
    3. Change parameters based on dataset using LFP_Denoising.py OR Spikes_Denoising.py as needed.
    4. Use Raster_plot function to compare raw (denoising=False) and denoised (denoising=True) raster data.
    5. Refine denoised datsets through repeption of steps 3-4 if needed. 

## **Data**

The following datasets have been provided for using and test running the DENOISING script.
  - Sample LFP-detected Dataset with Noise (XXX)
  - Sample Spike-detected Dataset with Noise (XXX)

## **Requirements**

Python >= 3.7; all analyses and testing were performed using Python 3.7 within PyCharm V.2023.2
